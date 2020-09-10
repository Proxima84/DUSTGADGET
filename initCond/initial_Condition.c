#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define AE 1.4959787e13
#define YR 3.155815e7

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

int NumPart;

struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;

    float U, Rho;
} * P;

int* Id;
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

    if (!(Id = malloc(NumPart * sizeof(int))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }

    printf("allocating memory...done\n");
}

int create_init(void)
{
    FILE* out;
    out = fopen("outbeg", "wb");
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;
    pc = 0;
    int Number = 0;
#define SKIP fwrite(&dummy, sizeof(dummy), 1, out);
    dummy = sizeof(header);
    printf("dummy=%i ", dummy);
    SKIP;
    fwrite(&header, sizeof(header), 1, out);
    SKIP;
    for (k = 0; k < 3; k++)
    {
        Number += header.npart[k];
    }
    dummy = sizeof(float) * Number * 3;
    SKIP;
    for (k = 0, pc_new = pc; k < 3; k++)
    {
        for (n = 0; n < header.npart[k]; n++)
        {
            fwrite(&P[pc_new].Pos[0], sizeof(float), 3, out);
            pc_new++;
        }
    }
    SKIP;
    SKIP;
    for (k = 0, pc_new = pc; k < 3; k++)
    {
        for (n = 0; n < header.npart[k]; n++)
        {
            fwrite(&P[pc_new].Vel[0], sizeof(float), 3, out);
            pc_new++;
        }
    }
    SKIP;
    dummy = sizeof(int) * Number;
    SKIP;
    for (k = 0, pc_new = pc; k < 3; k++)
    {
        for (n = 0; n < header.npart[k]; n++)
        {
            fwrite(&Id[pc_new], sizeof(int), 1, out);
            pc_new++;
        }
    }
    SKIP;
    dummy = sizeof(float) * Number;
    SKIP;
    for (k = 0, pc_new = pc; k < 3; k++)
    {
        for (n = 0; n < header.npart[k]; n++)
        {
            fwrite(&P[pc_new].Mass, sizeof(float), 1, out);
            pc_new++;
        }
    }
    SKIP;

    if (header.npart[0] > 0)
    {
        dummy = sizeof(float) * header.npart[0];
        SKIP;
        for (n = 0, pc_sph = pc; n < header.npart[0]; n++)
        {
            fwrite(&P[pc_sph].U, sizeof(float), 1, out);
            pc_sph++;
        }
        SKIP;
        SKIP;
        for (n = 0, pc_sph = pc; n < header.npart[0]; n++)
        {
            fwrite(&P[pc_sph].Rho, sizeof(float), 1, out);
            pc_sph++;
        }
        SKIP;
    }
}

/**************************************************************************/
int main(void)
{
    FILE *inGas, *inDust, *inStars;
    char gas[45], dust[45], star[45];
    printf("Enter file name with gas particles properties "
           "(x,y,z,vx,vy,vz,mass,U,id), if no gas particles enter any key:");
    scanf("%s", &gas);
    printf("Enter file name with dust particles properties "
           "(x,y,z,vx,vy,vz,mass,id), if no dust particles enter any key:");
    scanf("%s", &dust);
    printf("Enter file name with stars properties (x,y,z,vx,vy,vz,mass,id), if "
           "no stars particles enter any key:");
    scanf("%s", &star);
    printf("%s %s %s\n", gas, dust, star);
    int i, Ngas, Ndust, Nstars, n, id;
    Ngas = Ndust = Nstars = 0;
    double x, y, z, vx, vy, vz, U, m;
    inGas = fopen(gas, "r");
    if (inGas != NULL)
    {
        i = 0;
        while (!feof(inGas))
        {
            fscanf(inGas, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz,
                &m, &U, &id);
            i++;
        }
        Ngas = i;
        rewind(inGas);
    }

    inDust = fopen(dust, "r");
    if (inDust != NULL)
    {
        i = 0;
        while (!feof(inDust))
        {
            fscanf(inDust, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz, &m,
                &id);
            i++;
        }
        Ndust = i;
        rewind(inDust);
    }
    inStars = fopen(star, "r");
    if (inStars != NULL)
    {
        i = 0;
        while (!feof(inStars))
        {
            fscanf(inStars, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz,
                &m, &id);
            i++;
        }
        Nstars = i;
        rewind(inStars);
    }
    NumPart = Ngas + Ndust + Nstars;
    printf("Ngas=%i Ndust=%i Nstars=%i NumPart=%i\n", Ngas, Ndust, Nstars, NumPart);
    for (i = 0; i < 3; i++)
    {
        header.npart[i] = 0;
        header.mass[i] = 0.0;
        header.npartTotal[i] = 0;
    }
    header.npart[0] = Ngas;
    header.npart[1] = Ndust;
    header.npart[2] = Nstars;
    header.npartTotal[0] = Ngas;
    header.npartTotal[1] = Ndust;
    header.npartTotal[2] = Nstars;
    header.time = 0;
    header.flag_cooling = 0;
    header.flag_entropy_instead_u = 0;
    allocate_memory();
    n = 0;
    for (i = 0; i < Ngas; i++)
    {
        fscanf(inGas, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz, &m,
            &U, &id);
        P[n].Pos[0] = x;
        P[n].Pos[1] = y;
        P[n].Pos[2] = z;
        P[n].Vel[0] = vx;
        P[n].Vel[1] = vy;
        P[n].Vel[2] = vz;
        P[n].Mass = m;
        P[n].U = U;
        P[n].Type = 0;
        Id[n] = id;
        n++;
    }
    for (i = 0; i < Ndust; i++)
    {
        fscanf(
            inDust, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz, &m, &id);
        P[n].Pos[0] = x;
        P[n].Pos[1] = y;
        P[n].Pos[2] = z;
        P[n].Vel[0] = vx;
        P[n].Vel[1] = vy;
        P[n].Vel[2] = vz;
        P[n].Mass = m;
        P[n].U = 0;
        P[n].Type = 1;
        Id[n] = id;
        n++;
    }
    for (i = 0; i < Nstars; i++)
    {
        fscanf(
            inStars, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", &x, &y, &z, &vx, &vy, &vz, &m, &id);
        P[n].Pos[0] = x;
        P[n].Pos[1] = y;
        P[n].Pos[2] = z;
        P[n].Vel[0] = vx;
        P[n].Vel[1] = vy;
        P[n].Vel[2] = vz;
        P[n].Mass = m;
        P[n].U = 0;
        P[n].Type = 2;
        Id[n] = id;
        n++;
    }
    create_init();
}

