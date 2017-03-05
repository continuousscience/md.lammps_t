#include <stdio.h>
#include <sil_ext.h>
#include <proto_prim.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

// Required parts of LAMMPS library.h C-API.
void  lammps_open_no_mpi(int, char **, void **);
char *lammps_command(void *, char *);
void lammps_file(void *, char *);

#define DAT_BUFLEN 4096
typedef struct LmpDatum LmpDatum;
struct LmpDatum {
    LmpDatum *next;
    int fd;
    char name[18]; // "/tmp/lmp.XXXXXXXX"
                   // name is 17 chars + 1 zero-byte
};

// LmpDatum are used for loading coords and exporting vectors,
// as well as for serializing the LAMMPS state
// via a call to write_restart.

typedef struct {
    void *lmp;
    char *err; // NULL unless "An error occurred."
    LmpDatum *dat; // linked list
    int initialized;
} LAMMPS;

static const unsigned char lammps_hash[HASH_SIZE+1] = /*!hash!*/;

// This takes over responsibility for freeing lmp
static inline int sil_pushlammps(sil_State *S, LAMMPS *lmp) {
    return sil_newuserdata(S, lammps_hash, lmp);
}

// File Management API impl.
static int new_datum(LmpDatum **datp, const char *buf, size_t len) {
    LmpDatum *dat;
    char name[] = "/tmp/lmp.XXXXXXXX";
    int fd = mkstemp(name);
    int n = 0;

    if(fd < 0) {
        fprintf(stderr, "lammps_t: Error creating temp file.\n");
        return -1;
    }

    for(; *datp != NULL; datp = &(*datp)->next) n++;

    dat = (LmpDatum *)malloc(sizeof(LmpDatum));
    *datp = dat;
    dat->next = NULL;
    dat->fd = fd;
    memcpy(dat->name, name, 18);

    if(buf != NULL) {
        write(fd, buf, len);
    }
    return n;
}

static void free_datums(LmpDatum *dat) {
    LmpDatum *next;
    for(; dat != NULL; dat = next) {
        next = dat->next;
        close(dat->fd);
        unlink(dat->name);
        free(dat);
    }
}

static LmpDatum **get_datum(LmpDatum **dat, int n) {
    for(; n > 0; n--) {
        if(*dat == NULL)
            return dat;
        dat = &(*dat)->next;
    }
    return dat;
}

static size_t datum_size(LmpDatum *dat, int n) {
    LmpDatum **datp = get_datum(&dat, n);
    if(!datp) return 0;
    struct stat st;
    dat = *datp;

    if(fstat(dat->fd, &st) < 0) {
        return 0;
    }
    return st.st_size;
}

static char *read_datum(struct LmpDatum *dat, int n, size_t *lenp) {
    size_t len = datum_size(dat, n);
    char *buf = (char *)malloc(len);
    *lenp = len;
    if(len == 0) return buf;

    dat = *get_datum(&dat, n);
    lseek(dat->fd, 0, SEEK_SET);
    read(dat->fd, buf, len); // should be len...
    return buf;
}

// Files can be cleared, but not deleted (until dtor of course).
static int put_datum(LmpDatum *dat, int n, const char *buf, size_t len) {
    dat = *get_datum(&dat, n);
    if(!dat) return 1;
    if(dat->fd >= 0) {
        ftruncate(dat->fd, 0);
    }
    if(buf != NULL && len > 0)
        return write(dat->fd, buf, len) != len;
    return 0;
}

static LmpDatum *copy_datums(LmpDatum *dat) {
    char buf[DAT_BUFLEN];
    size_t len;
    LmpDatum *out = NULL;
    LmpDatum *end;
    for(; dat != NULL; dat = dat->next) {
        int n = new_datum(&out, NULL, 0);
        if(n < 0) {
            goto err;
        }
        end = *get_datum(&out, n);
        lseek(dat->fd, 0, SEEK_SET);
        while( (len = read(dat->fd, buf, DAT_BUFLEN)) > 0) {
            if(write(end->fd, buf, len) != len) {
                goto err;
            }
        }
    }
    return out;
err:
    free_datums(out);
    return NULL;
}

/*** Initial LAMMPS Startup ***/
// returns 0 on success or -1 on failure
// This routine is unable to detect lammps issues!
//
// Note: The caller is responsible for creating datum '0'
// which holds restart information.
// The returned LAMMPS structure is set to 'initialized = 0'.
static LAMMPS *open_lammps() {
    LAMMPS *lmp = (LAMMPS *)calloc(1, sizeof(LAMMPS));
    lammps_open_no_mpi(0, NULL, &lmp->lmp);

    return lmp;
}

// lmp->dat == NULL -- create a new 'restart' file and ignore init flag.
//   init == 0 : run commands in file directly
//   init != 0 : read 'restart' file
static void startup_lammps(LAMMPS *lmp, int init) {
    if(lmp->dat == NULL) {
        new_datum(&lmp->dat, NULL, 0);
        lmp->initialized = 0;
        return;
    }
    if(init == 0) {
        lammps_file(lmp->lmp, lmp->dat->name);
        lmp->initialized = 0;
    } else {
        char cmd[] = "read_restart /tmp/lmp.XXXXXXXX";
        memcpy(cmd+13, lmp->dat->name, 17);
        lammps_command(lmp->lmp, cmd);
        lmp->initialized = 1;
    }
}

