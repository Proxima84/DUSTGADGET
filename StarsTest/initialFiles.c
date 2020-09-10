#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define Msun 1.98892e33
#define AE 1.4959787e13
#define YR 3.155815e7
#define SOLAR_RADIUS 0.004652
#define flaredFac 0.33437
#define hydroMassFrac 0.23 /*mass fraction of hydrogen*/
#define heightFac                                                                                  \
    0.305e-2 /*sqrt(kappa*AU^3/(G*Msun*mp))/AU kappa-Boltzman constant,          \                 \
mp-proton mass */
#define rhoFac 2.5066 /*sqrt(2PI)*/
#define GAMMA 1
#define dtg 0.01 /* dust to gas mass ratio*/
#define NG 65 /*initial gas number*/
#define ND 65 /*initial dust number*/
#define Rin 0.1 /*inner border of the disk*/
#define Rout 5. /*outer border of the disk*/
#define diskMass 0.01 /*the disk mass as a fraction of the stars mass sum*/
#define U_koef 0.367e-3 /*kappa*YR^2/mp/AE^2 kappa-Boltzman constant, mp-proton mass*/

int NumPart, Ngas;

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
    printf("allocating memory...done\n");
}

/**************************************************************************/
int main(void)
{
    int k, i, l, n, Rs, hs, N0, Ni;
    double t_unit, m_unit, v_unit, r_unit;
    double R, Rstar, hr, hz, z, v, Sigma;
    float period, mass, mass1, mass2, e, ang, Ts, molecWigth;
    double r_m, r_m_2, dr, dz, rr, fi, sinfi, cosfi, a;
    double Mp, Md, Tmid;
    r_m = RAND_MAX;
    r_m_2 = (RAND_MAX - 1) / 2 + 1;
    FILE *out, *outS, *outG, *outD;
    out = fopen("description.in", "wt");
    printf("mass of the primary M1 (Msun):");
    scanf("%f", &mass1);
    fprintf(out, "M1=%f Msun\n", mass1);
    printf("mass of the secondary M2 (Msun):");
    scanf("%f", &mass2);
    fprintf(out, "M2=%f Msun\n", mass2);
    printf("eccentricity e:");
    scanf("%f", &e);
    fprintf(out, "e=%1.3lf\n", e);
    printf("the orbit inclination i:");
    scanf("%f", &ang);
    fprintf(out, "i=%1.3lf\n", ang);
    ang *= M_PI / 180.;
    printf("the period P (year):");
    scanf("%f", &period);
    fprintf(out, "P=%1.3lf\n", period);
    mass = mass1 + mass2;
    a = pow(mass * period * period, 0.33333333);
    fprintf(out, "a=%1.3lf\n", a);
    printf("Central Temperature (K):");
    scanf("%f", &Ts);
    fprintf(out, "Ts=%lf\n", Ts);
    molecWigth = 4.0 / (1 + 3 * hydroMassFrac);
    t_unit = period;
    r_unit = a;
    m_unit = mass / (4 * M_PI * M_PI);
    v_unit = a / period;
    ;
    printf("\n                            UNITS\n\n");
    printf("%f UNIT of MASS     = %.2lf MASS of SUN\n", mass, m_unit);
    printf("UNIT of DISTANCE = %.16lf\n", r_unit);
    printf("UNIT of VELOCITY = %.16lf\n", v_unit);
    printf("UNIT of TIME     = %.3lf years\n\n", t_unit);
    printf("MASS of BINARY   = %.2lf MASS of SUN\n", 4 * M_PI * M_PI * m_unit);
    printf("PERIOD           = %.3lf years\n\n", t_unit);
    Rstar = SOLAR_RADIUS / r_unit;
    dr = 0.05;
    dz = 0.01;
    Rs = floor((Rout - Rin) / dr);
    N0 = NG;
    n = 0;
    float S0, S;
    S0 = 2 * Rin + dr;
    R = Rin;
    for (i = 0; i < Rs; i++)
    {
        Tmid = flaredFac * sqrt(Rstar / R) * Ts;
        hz = heightFac
            * sqrt(GAMMA * R * r_unit * R * r_unit * R * r_unit * Tmid / molecWigth / mass)
            / r_unit;
        hr = 3 * hz;
        hs = 2 * floor(hr / dz);
        z = -(hs - 1) * dz / 2.;
        S = 2 * R + dr;
        Sigma = floor(N0 * S / S0 * Rin / R);
        for (k = 0; k < hs; k++)
        {
            Ni = floor(Sigma / rhoFac / hz * exp(-z * z / (2. * hz * hz)));
            n += Ni;
            z += dz;
        }
        R += dr;
    }
    Ngas = n;
    int Ndust;
    N0 = ND;
    n = 0;
    R = Rin;
    for (i = 0; i < Rs; i++)
    {
        Tmid = flaredFac * sqrt(Rstar / R) * Ts;
        hz = heightFac
            * sqrt(GAMMA * R * r_unit * R * r_unit * R * r_unit * Tmid / molecWigth / mass)
            / r_unit;
        hr = 3 * hz;
        hs = 2 * floor(hr / dz);
        z = -(hs - 1) * dz / 2.;
        S = 2 * R + dr;
        Sigma = floor(N0 * S / S0 * Rin / R);
        for (k = 0; k < hs; k++)
        {
            Ni = floor(Sigma / rhoFac / hz * exp(-z * z / (2. * hz * hz)));
            n += Ni;
            z += dz;
        }
        R += dr;
    }
    Ndust = n;
    NumPart = Ngas + Ndust + 2;
    allocate_memory();
    Mp = diskMass * mass / Ngas;
    printf("Ngas=%i Ndust=%i NumPart=%i\n", Ngas, Ndust, NumPart);
    fprintf(out, "Mp=%1.3e\n", Mp * m_unit * Msun);
    fprintf(out, "Ngas=%i\n", Ngas);
    fprintf(out, "Md=%1.3e\n", Mp * m_unit * Msun * dtg);
    fprintf(out, "Ndust=%i\n", Ndust);
    fprintf(out, "UnitVelocity_in_cm_per_s %e\n", v_unit * AE / YR);
    fprintf(out, "UnitMass_in_g %e\n", m_unit * Msun);
    fprintf(out, "UnitLength_in_cm %e\n", a * AE);
    fclose(out);
    srand(time(0));
    n = 0;
    float pos1, pos2, v1, v2, VA;
    pos1 = -mass2 / mass * (1 + e);
    pos2 = mass1 / mass * (1 + e);
    VA = 2 * M_PI * sqrt((1 - e) / (1 + e));
    v1 = mass2 / mass * VA;
    v2 = -mass1 / mass * VA;
    int Nstar;
    Nstar = Ngas + Ndust;
    P[Nstar].Pos[0] = pos1 * cos(ang);
    P[Nstar + 1].Pos[0] = pos2 * cos(ang);
    P[Nstar].Pos[1] = 0.;
    P[Nstar + 1].Pos[1] = 0.;
    P[Nstar].Pos[2] = pos1 * sin(ang);
    P[Nstar + 1].Pos[2] = pos2 * sin(ang);
    P[Nstar].Vel[0] = 0;
    P[Nstar + 1].Vel[0] = 0;
    P[Nstar].Vel[1] = v1;
    P[Nstar + 1].Vel[1] = v2;
    P[Nstar].Vel[2] = 0;
    P[Nstar + 1].Vel[2] = 0;
    P[Nstar].Mass = mass1 / mass * 4 * M_PI * M_PI;
    P[Nstar + 1].Mass = mass2 / mass * 4 * M_PI * M_PI;
    P[Nstar].Type = 2;
    P[Nstar + 1].Type = 2;
    P[Nstar].U = 0;
    P[Nstar + 1].U = 0;
    double TmR, TaR, RR;
    N0 = NG;
    R = Rin;
    for (i = 0; i < Rs; i++)
    {
        Tmid = flaredFac * sqrt(Rstar / R) * Ts;
        hz = heightFac
            * sqrt(GAMMA * R * r_unit * R * r_unit * R * r_unit * Tmid / molecWigth / mass)
            / r_unit;
        hr = 3 * hz;
        hs = 2 * floor(hr / dz);
        z = -(hs - 1) * dz / 2.;
        S = 2 * R + dr;
        Sigma = floor(N0 * S / S0 * Rin / R);
        for (k = 0; k < hs; k++)
        {
            Ni = floor(Sigma / rhoFac / hz * exp(-z * z / (2. * hz * hz)));
            for (l = 0; l < Ni; l++)
            {
                rr = R + dr * (rand() / r_m);
                fi = 2 * M_PI * l / Ni + 2 * M_PI / Ni * (rand() / r_m);
                sinfi = sin(fi);
                cosfi = cos(fi);
                P[n].Pos[0] = rr * cosfi;
                P[n].Pos[1] = rr * sinfi;
                P[n].Pos[2] = z - dz / 2. + dz * (rand() / r_m);
                RR = sqrt(rr * rr + P[n].Pos[2] * P[n].Pos[2]);
                v = sqrt(4 * M_PI * M_PI / RR);
                P[n].Vel[0] = v * sinfi;
                P[n].Vel[1] = -v * cosfi;
                P[n].Vel[2] = 0;
                P[n].Mass = Mp / m_unit;
                TmR = flaredFac * sqrt(Rstar / rr) * Ts;
                P[n].U = U_koef * TmR / molecWigth / (v_unit * v_unit);
                P[n].Type = 0;
                n++;
            }
            z += dz;
        }
        R += dr;
    }
    N0 = ND;
    R = Rin;
    for (i = 0; i < Rs; i++)
    {
        Tmid = flaredFac * sqrt(Rstar / R) * Ts;
        hz = heightFac
            * sqrt(GAMMA * R * r_unit * R * r_unit * R * r_unit * Tmid / molecWigth / mass)
            / r_unit;
        hr = 3 * hz;
        hs = 2 * floor(hr / dz);
        z = -(hs - 1) * dz / 2.;
        S = 2 * R + dr;
        Sigma = floor(N0 * S / S0 * Rin / R);
        for (k = 0; k < hs; k++)
        {
            Ni = floor(Sigma / rhoFac / hz * exp(-z * z / (2. * hz * hz)));
            for (l = 0; l < Ni; l++)
            {
                rr = R + dr * (rand() / r_m);
                fi = 2 * M_PI * l / Ni + 2 * M_PI / Ni * (rand() / r_m);
                sinfi = sin(fi);
                cosfi = cos(fi);
                P[n].Pos[0] = rr * cosfi;
                P[n].Pos[1] = rr * sinfi;
                P[n].Pos[2] = z - dz / 2. + dz * (rand() / r_m);
                RR = sqrt(rr * rr + P[n].Pos[2] * P[n].Pos[2]);
                v = sqrt(4 * M_PI * M_PI / RR);
                P[n].Vel[0] = v * sinfi;
                P[n].Vel[1] = -v * cosfi;
                P[n].Vel[2] = 0.;
                P[n].Mass = Mp / m_unit * dtg;
                P[n].U = 0;
                P[n].Type = 1;
                n++;
            }
            z += dz;
        }
        R += dr;
    }
    outS = fopen("stars.in", "wt");
    outG = fopen("gas.in", "wt");
    outD = fopen("dust.in", "wt");
    for (i = 0; i < NumPart; i++)
    {
        if (P[i].Type == 0)
            fprintf(outG, "%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t%i\n", P[i].Pos[0], P[i].Pos[1],
                P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass, P[i].U, i);
        else if (P[i].Type == 1)
            fprintf(outD, "%f\t%f\t%f\t%f\t%f\t%f\t%e\t%i\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass, i);
        else
            fprintf(outS, "%f\t%f\t%f\t%f\t%f\t%f\t%e\t%i\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass, i);
    }
    fclose(out);
    fclose(outG);
    fclose(outD);
    fclose(outS);
}

