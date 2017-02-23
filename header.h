#include <stdio.h>
#include <sil_ext.h>
#include <proto_prim.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct {
    void *lmp;
    char name[20]; // these two are for serializing the LAMMPS structure
    int fd; // via a call to write_restart
} LAMMPS;

// basic required API for this file
char *lammps_command(void *, char *);
void  lammps_open_no_mpi(int, char **, void **);

static const unsigned char lammps_hash[HASH_SIZE+1] = /*!hash!*/;

// This takes over responsibility for freeing lmp
static inline int sil_pushlammps(sil_State *S, LAMMPS *lmp) {
    return sil_newuserdata(S, lammps_hash, lmp);
}

// returns 0 on success or -1 on failure
// This routine is unable to detect lammps issues!
static LAMMPS *open_lammps(const void *buf, size_t len) {
    LAMMPS *lmp = NULL;
    if(buf == NULL) {
        lmp = (LAMMPS *)malloc(sizeof(LAMMPS));
        lmp->fd = -1;
        lammps_open_no_mpi(0, NULL, &lmp->lmp);
        return lmp;
    }
    char cmd[] = "read_restart /tmp/lmp_r.XXXXXXXX";
    int fd = mkstemp(cmd+14);

    if(fd < 0) {
        fprintf(stderr, "lammps_t: Error creating temp file.\n");
        return NULL;
    }

    if(write(fd, buf, len) != len) {
        goto err;
    }

    lmp = (LAMMPS *)malloc(sizeof(LAMMPS));
    lmp->fd = -1;
    lammps_open_no_mpi(0, NULL, &lmp->lmp);

    lammps_command(lmp->lmp, cmd);

err:
    close(fd);
    unlink(cmd+14);
    return lmp;
}

