#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

// Basic required API for this file.
void  lammps_close(void *);
int   lammps_get_natoms(void *);

// Returns 0 on success or -1 otherwise.
static int write_dump(LAMMPS *lmp) {
    char cmd[] = "write_dump all custom ${datum1} id type x y z vx vy vz";

    if(!lmp->initialized) return 0; // not needed.

    if(lmp->dat == NULL || lmp->dat->next == NULL) {
        fprintf(stderr, "lammps_t: Error creating temp file.\n");
        return -1;
    }

    put_datum(lmp->dat, 1, NULL, 0);
    memcpy(cmd+22, lmp->dat->next->name, 17);
    lammps_command(lmp->lmp, cmd);
    return 0;
}


char *show(const LAMMPS *lmp) {
    char *out;
    int n = 0;
    for(LmpDatum *dat = lmp->dat; dat != NULL; dat = dat->next) n++;
    asprintf(&out, "LAMMPS object (%d atoms, %d dats)",
                    lammps_get_natoms(lmp->lmp), n);
    return out;
}

// [steps] [init?] [uint64 size] [size bytes] [uint64 size] [size bytes] ...
// file 0 is a list of commands if !init
// else file 0 is write_restart data.
size_t size(LAMMPS *lmp) {
    size_t len = size_uint64(lmp->steps)
               + size_uint32(lmp->initialized);
    if(write_dump(lmp)) {
        return 0;
    }
    
    int n = 0;
    for(LmpDatum *dat = lmp->dat; dat != NULL; dat = dat->next) {
        size_t tot = datum_size(lmp->dat, n++);
        len += size_uint64(tot) + tot;
    }
    return len;
}

// FIXME: If files change betwen size() and serialize(),
// then this will overflow/underflow the output message.
// This shouldn't happen in practice unless lammps (or the system)
// does something asynchronously.
void serialize(SWriter *s, LAMMPS *lmp) {
    char buf[DAT_BUFLEN];

    write_uint64(s, lmp->steps);
    write_uint32(s, lmp->initialized != 0);

    // Serialize all files (size is always called first).
    int n = 0;
    for(LmpDatum *dat = lmp->dat; dat != NULL; dat = dat->next) {
        size_t tot = datum_size(lmp->dat, n++);
        size_t len = 0;
        write_uint64(s, tot);
        lseek(dat->fd, 0, SEEK_SET);
        while(tot > 0) {
            if( (len = read(dat->fd, buf, DAT_BUFLEN)) == 0) { // cut short?
                len = tot < DAT_BUFLEN ? tot : DAT_BUFLEN;
                memset(buf, 0, len);
            }
            if(len > tot) len = tot;
            s->write(s->stream, buf, len);
            tot -= len;
        }
    }
    put_datum(lmp->dat, 1, NULL, 0);
}

void handler(LAMMPS *lmp) {
    LmpDatum *next;

    printf("Closing LAMMPS %p.\n", lmp);
    free_datums(lmp->dat);
    lammps_close(lmp->lmp);
}

void parse(sil_State *S, const uint8_t *buf, size_t len) {
    LAMMPS *lmp;
    unsigned k;
    uint64_t steps;
    uint32_t init;

    if(!len) {
        sil_err(S, "Unable to read LAMMPS restart file.");
        return;
    }
    k = read_uint64(&steps, buf, len);
    len -= k; buf += k;
    k = read_uint32(&init, buf, len);
    len -= k; buf += k;
    if(!len) {
        sil_err(S, "Unable to read LAMMPS restart file.");
        return;
    }

    if( (lmp = open_lammps()) == NULL) {
        sil_err(S, "Unable to open LAMMPS instance.");
        return;
    }

    while(len) {
        uint64_t tot;
        k = read_uint64(&tot, buf, len);
        len -= k; buf += k;
        tot = tot < len ? tot : len;
        if(new_datum(&lmp->dat, (char *)buf, tot) < 0) {
            goto err;
        }
        
        len -= tot; buf += tot;
    }
    startup_lammps(lmp, init, steps);

    sil_pushlammps(S, lmp);
    return;
err:
    handler(lmp);
    sil_err(S, "Unable to open LAMMPS instance.");
    return;
}

void copy(sil_State *S) {
    size_t len; // Always assume len is wrong!
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    LAMMPS *lmp_copy = NULL;
    LmpDatum *dat = NULL;

    if(lmp == NULL) {
        sil_err(S, "Can't copy - no LAMMPS present.");
        goto err;
    }

    if(write_dump(lmp)) {
        sil_err(S, "Can't copy - lammps dump error.");
        goto err;
    }
    if( (dat = copy_datums(lmp->dat)) == NULL) {
        sil_err(S, "Can't copy - datum copy error.");
        goto err;
    }
    if( (lmp_copy = open_lammps()) == NULL) {
        sil_err(S, "Can't copy - LAMMPS restart error.");
        free_datums(dat);
        goto err;
    }
    lmp_copy->dat = dat;
    startup_lammps(lmp_copy, lmp->initialized, lmp->steps);
    put_datum(lmp->dat, 1, NULL, 0);

err:
    sil_setST(S, lmp_copy, sizeof(LAMMPS));
}

