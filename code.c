#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>

// Extra API bits required for this file.
void  lammps_close(void *);
int   lammps_get_natoms(void *);

// FIXME! storing lmp->fd is not thread-safe for (rare)
// simultaneous read/write/copy operations on the same LAMMPS object.
// We should add a mutex to the internal LAMMPS struct.

static void scrub_restart(LAMMPS *lmp) {
    if(lmp->fd >= 0) {
        close(lmp->fd);
        lmp->fd = -1;
        unlink(lmp->name);
    }
}

// name is 17 chars + 1 nil
static size_t fd_size(LAMMPS *lmp) {
    struct stat st;
    if(lmp->fd < 0) return 0;

    if(fstat(lmp->fd, &st) < 0) {
        scrub_restart(lmp);
        return 0;
    }
    return st.st_size;
}

// Create a restart file for the lammps object if one does not exist.
// returns 0 on success or -1 otherwise.
//
// TODO: create fake restart file if LAMMPS is not initialized...
static int mk_restart(LAMMPS *lmp) {
    if(lmp->fd >= 0) return 0;

    char cmd[] = "write_restart /tmp/lmp->XXXXXXXX";
    int fd = mkstemp(cmd+14);

    if(fd < 0) {
        fprintf(stderr, "lammps_t: Error creating temp file.\n");
        return -1;
    }

    lammps_command(lmp->lmp, cmd);
    lmp->fd = fd;
    memcpy(lmp->name, cmd+14, 18);
    return 0;
}

char *show(const LAMMPS *lmp) {
    char *out;
    asprintf(&out, "LAMMPS object (%d atoms)", lammps_get_natoms(lmp->lmp));
    return out;
}

size_t size(LAMMPS *lmp) {
    if(mk_restart(lmp)) {
        return 0;
    }
    return fd_size(lmp);
}

#define BUFLEN 4096
void serialize(SWriter *s, LAMMPS *lmp) {
    char buf[BUFLEN];
    int r;
    if(mk_restart(lmp)) {
        return;
    }
    while( (r = read(lmp->fd, buf, BUFLEN)) > 0) {
        s->write(s->stream, buf, r);
    }
    scrub_restart(lmp);
}

void parse(sil_State *S, const void *buf, size_t len) {
    LAMMPS *lmp = open_lammps(buf, len);
    if(lmp == NULL) {
        sil_err(S, "Unable to read LAMMPS restart file.");
        return;
    }
    sil_pushlammps(S, lmp);
    return;
}

void handler(LAMMPS *lmp) {
    printf("Closing LAMMPS %p.\n", lmp);
    lammps_close(lmp->lmp);
}

void copy(sil_State *S) {
    size_t len; // Always assume len is wrong!
    void *buf = NULL;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    LAMMPS *lmp_copy = NULL;

    if(lmp == NULL) {
        sil_err(S, "Can't copy - no LAMMPS present.");
        return;
    }
    len = size(lmp);
    if(lmp->fd < 0) {
        sil_err(S, "Can't copy - LAMMPS write error.");
        return;
    }
    if( (buf = malloc(len)) == NULL) {
        sil_err(S, "Can't copy - memory error.");
        return;
    }
    if(read(lmp->fd, buf, len) != len) {
        scrub_restart(lmp);
        sil_err(S, "Can't copy - read error.");
        goto err;
    }
    scrub_restart(lmp);

    if( (lmp_copy = open_lammps(buf, len)) == NULL) {
        sil_err(S, "Can't copy - LAMMPS restart error.");
        goto err;
    }

    sil_setST(S, lmp_copy, sizeof(LAMMPS));
err:
    if(buf) free(buf);
}

