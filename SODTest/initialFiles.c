#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define A 0.0001
#define eps 0.01
#define GasDensity1 1
#define GasDensity2 0.125
#define limit 0.5
#define N 2400

int NumPart, Ngas;
int* Id;
struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;
    float U;
} * P;

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

void fill_x(double* x, int N1, int N2, double Ri, double Re)
{
    // step = length / number
    double step_left = fabs(Ri) / N1;
    double step_right = Re / N2;
    int i;
    // write left "real" into both arrays
    x[0] = Ri + (step_left / 2);
    for (i = 1; i < N1; i++)
    {
        x[i] = x[i - 1] + step_left;
    }
    // write right real into both arrays
    x[N1] = (step_right / 2);
    for (i = 1; i < N2; i++)
    {
        x[N1 + i] = x[N1 + i - 1] + step_right;
    }
}


/**************************************************************************/
int main(void)
{

    NumPart = 4 * N;
    allocate_memory();
    double Ri, Re;
    Ri = -limit - 0.1;
    Re = limit + 0.1;
    double mass, step;
    step = 1. / N;
    mass = 1. / N * Re;
    double R, *x;
    int i, n, frac;
    frac = GasDensity1 / GasDensity2;
    x = malloc(NumPart * sizeof(double));
    fill_x(x, N, N, Ri, Re);
    n = 0;
    for (i = 0; i < N; i++)
    {
        P[n].Pos[0] = x[i];
        P[n].Pos[1] = 0.;
        P[n].Pos[2] = 0.;
        P[n].Vel[0] = 0.;
        P[n].Vel[1] = 0.;
        P[n].Vel[2] = 0.;
        P[n].U = 2.5;
        P[n].Type = 0;
        P[n].Mass = mass;
        if (x[i] < -limit)
            Id[n] = 0;
        else
            Id[n] = n;
        n++;
    }
    for (i = 0; i < N; i++)
    {
        P[n].Pos[0] = x[N + i];
        P[n].Pos[1] = 0.;
        P[n].Pos[2] = 0.;
        P[n].Vel[0] = 0.;
        P[n].Vel[1] = 0.;
        P[n].Vel[2] = 0.;
        P[n].U = 2.;
        P[n].Type = 0;
        P[n].Mass = mass / (1. * frac);
        if (x[N + i] > limit)
            Id[n] = 0;
        else
            Id[n] = n;
        n++;
    }
    for (i = 0; i < N; i++)
    {
        P[n].Pos[0] = x[i];
        P[n].Pos[1] = 0.;
        P[n].Pos[2] = 0.;
        P[n].Vel[0] = 0.;
        P[n].Vel[1] = 0.;
        P[n].Vel[2] = 0.;
        P[n].U = 0;
        P[n].Type = 1;
        P[n].Mass = mass;
        if (x[i] < -limit)
            Id[n] = 0;
        else
            Id[n] = n;
        n++;
    }
    for (i = 0; i < N; i++)
    {
        P[n].Pos[0] = x[N + i];
        P[n].Pos[1] = 0.;
        P[n].Pos[2] = 0.;
        P[n].Vel[0] = 0.;
        P[n].Vel[1] = 0.;
        P[n].Vel[2] = 0.;
        P[n].U = 0;
        P[n].Type = 1;
        P[n].Mass = mass / (1. * frac);
        if (x[N + i] > limit)
            Id[n] = 0;
        else
            Id[n] = n;
        n++;
    }
    free(x);
    FILE *outG, *outD;
    outG = fopen("gas.in", "wt");
    outD = fopen("dust.in", "wt");
    for (i = 0; i < NumPart; i++)
    {
        if (P[i].Type == 0)
            fprintf(outG, "%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t%i\n", P[i].Pos[0], P[i].Pos[1],
                P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass, P[i].U, Id[i]);
        else if (P[i].Type == 1)
            fprintf(outD, "%f\t%f\t%f\t%f\t%f\t%f\t%e\t%i\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass, Id[i]);
    }
    fclose(outG);
    fclose(outD);
}

