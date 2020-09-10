#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
struct io_header
{
    int npart[3];
    double mass[3];
    double time;
    int npartTotal[3];
    int flag_cooling;
    int flag_entropy_instead_u;
    char fill[184]; /* fills to 256 Bytes */
} header;

int NumPart, Ngas;

struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;

    float U, Rho;
} * P;


int* Id;

double Time; 

/* this routine allocates the memory for the
 * particle data.
 */
int allocate_memory(void)
{
    printf("allocating memory...\n");

    if (!(P = malloc(NumPart * sizeof(struct particle_data))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }

    P--; /* start with offset 1 */


    if (!(Id = malloc(NumPart * sizeof(int))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }

    Id--; /* start with offset 1 */

    printf("allocating memory...done\n");
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char* fname, int files)
{
    FILE* fd;
    char buf[200];
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

    for (i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
        if (files > 1)
            sprintf(buf, "%s.%d", fname, i);
        else
            sprintf(buf, "%s", fname);

        if (!(fd = fopen(buf, "r")))
        {
            printf("can't open file `%s`\n", buf);
            exit(0);
        }

        printf("reading `%s' ...\n", buf);
        fflush(stdout);

        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&header, sizeof(header), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);

        if (files == 1)
        {
            for (k = 0, NumPart = 0, ntot_withmasses = 0; k < 3; k++)
                NumPart += header.npart[k];
            Ngas = header.npart[0];
        }
        else
        {
            for (k = 0, NumPart = 0, ntot_withmasses = 0; k < 3; k++)
                NumPart += header.npartTotal[k];
            Ngas = header.npartTotal[0];
        }

        for (k = 0, ntot_withmasses = 0; k < 3; k++)
        {
            if (header.mass[k] == 0)
                ntot_withmasses += header.npart[k];
        }

        if (i == 0)
            allocate_memory();

        SKIP;
        for (k = 0, pc_new = pc; k < 3; k++)
        {
            for (n = 0; n < header.npart[k]; n++)
            {
                fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for (k = 0, pc_new = pc; k < 3; k++)
        {
            for (n = 0; n < header.npart[k]; n++)
            {
                fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;


        SKIP;
        for (k = 0, pc_new = pc; k < 3; k++)
        {
            for (n = 0; n < header.npart[k]; n++)
            {
                fread(&Id[pc_new], sizeof(int), 1, fd);
                pc_new++;
            }
        }
        SKIP;


        if (ntot_withmasses > 0)
            SKIP;
        for (k = 0, pc_new = pc; k < 3; k++)
        {
            for (n = 0; n < header.npart[k]; n++)
            {
                P[pc_new].Type = k;   
                fread(&P[pc_new].Mass, sizeof(float), 1, fd);              
                pc_new++;
            }
        }
        if (ntot_withmasses > 0)
            SKIP;


        if (header.npart[0] > 0)
        {
            SKIP;
            for (n = 0, pc_sph = pc; n < header.npart[0]; n++)
            {
                fread(&P[pc_sph].U, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;
        }

        SKIP;
        for (k = 0, pc_new = pc; k < 3; k++)
        {
            for (n = 0; n < header.npart[k]; n++)
            {
                fread(&P[pc_new].Rho, sizeof(float), 1, fd);
                pc_new++;
            }
        }
        SKIP;       
        fclose(fd);
    }

    Time = header.time;
}


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char** argv)
{
    char inname[200], outname[200];
    inname[0] = '\0';
    strcat(inname, argv[1]);
    /* number of files per snapshot */
    load_snapshot(inname, 1);   
    FILE* out;
    outname[0] = '\0';
    strcat(outname, "./dat/");
    strcat(outname, argv[1]);
    strcat(outname, ".dat");
    printf("%s", outname);
    out = fopen(outname, "w");
    int i;
    for (i = 1; i <= NumPart; i++)
    {
        fprintf(out, "%f %f %e %f %f %i %i\n", P[i].Pos[0], P[i].Vel[0], P[i].Rho, P[i].U,
            P[i].U * 0.4 * P[i].Rho, P[i].Type, Id[i]);
    }
    fclose(out);
}

