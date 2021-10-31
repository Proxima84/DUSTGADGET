#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define Xpart 100
#define Rgas 1.0
#define Vgas 1.0
#define E0 18.6
#define GasMass 1000.0	
#define eps 0.8
#define H 0.1


int     NumPart, Ngas;

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  double U;
} *P;


int* Id;

/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
 printf("allocating memory...\n");
 if(!(P=malloc(NumPart*sizeof(struct particle_data))))
 {
  fprintf(stderr,"failed to allocate memory.\n");
  exit(0);
 }
 if (!(Id = malloc(NumPart * sizeof(int))))
 {
   fprintf(stderr, "failed to allocate memory.\n");
   exit(0);
 }

 printf("allocating memory...done\n");
}


/**************************************************************************/
int main(void)                          
{
 int i,j,k,n,l,Ngas,Ndust; 
 n=0;
 double step,R,mass,Rball;
 double x,y,z;
 Rball=Rgas+2.*H;
 step=2.*Rball/Xpart; 
 l=0;
 for(i = 0; i < Xpart; i++)
  for(j = 0; j < Xpart; j++)
   for(k = 0; k < Xpart; k++)    
   {
    x = -Rball+(i+0.5)*step;
    y = -Rball+(j+0.5)*step;
    z = -Rball+(k+0.5)*step;
    R=sqrt(x*x+y*y+z*z);   
    if(R<Rball)
     n++;   
    if(R<Rgas)
     l++;
   }
  printf("n=%i l=%i %f \n",n,l,Rball);
  Ngas=n;
  Ndust=n;
  NumPart=2*n;
  allocate_memory();
  mass=step*step*step;
  n=0; 
  for(i = 0; i < Xpart; i++)
  for(j = 0; j < Xpart; j++)
   for(k = 0; k < Xpart; k++)    
   {
    P[n].Pos[0] = -Rball+(i+0.5)*step;
    P[n].Pos[1] = -Rball+(j+0.5)*step;
    P[n].Pos[2] = -Rball+(k+0.5)*step;
    R=sqrt(P[n].Pos[0]*P[n].Pos[0]+P[n].Pos[1]*P[n].Pos[1]+P[n].Pos[2]*P[n].Pos[2]);    
    P[n].Vel[0] = Vgas*P[n].Pos[0]/Rgas;    
    P[n].Vel[1] = Vgas*P[n].Pos[1]/Rgas;
    P[n].Vel[2] = Vgas*P[n].Pos[2]/Rgas;
    P[n].Mass=GasMass*2./NumPart;
    P[n].Type=0;  
    P[n].U=E0*(1.0 - (R/Rgas)*(R/Rgas)); 
    if(R<Rgas)  
    {
     Id[n] = n+1;       
     n++;
    }    
    else if(R<Rball)
    {
     Id[n]=0;
     n++;
    }     
   }
 for(i = 0; i < Xpart; i++)
  for(j = 0; j < Xpart; j++)
   for(k = 0; k < Xpart; k++)    
   {
    P[n].Pos[0] = -Rball+(i+0.5)*step+1.0e-4*(rand()/RAND_MAX-0.5);
    P[n].Pos[1] = -Rball+(j+0.5)*step+1.0e-4*(rand()/RAND_MAX-0.5);
    P[n].Pos[2] = -Rball+(k+0.5)*step+1.0e-4*(rand()/RAND_MAX-0.5);
    R=sqrt(P[n].Pos[0]*P[n].Pos[0]+P[n].Pos[1]*P[n].Pos[1]+P[n].Pos[2]*P[n].Pos[2]);     
    P[n].Vel[0] = Vgas*P[n].Pos[0]/Rgas;
    P[n].Vel[1] = Vgas*P[n].Pos[1]/Rgas;
    P[n].Vel[2] = Vgas*P[n].Pos[2]/Rgas;
    P[n].Mass=eps*GasMass*2./NumPart;
    P[n].Type=1;  
    P[n].U=0.;       
   if(R<Rgas)  
    {
     Id[n] = n+1;       
     n++;
    }    
    else if(R<Rball)
    {
     Id[n]=0;
     n++;
    }      
   }
 FILE *outG,*outD;
 outG=fopen("gas.in","wt");
 outD=fopen("dust.in","wt");
 for(i = 0; i < NumPart; i++)
 {
  if(P[i].Type==0)
   fprintf(outG,"%f\t%f\t%f\t%f\t%f\t%f\t%e\t%e\t%i\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],P[i].Vel[0],P[i].Vel[1],P[i].Vel[2],P[i].Mass,P[i].U,Id[i]);
  else if(P[i].Type==1)
   fprintf(outD,"%f\t%f\t%f\t%f\t%f\t%f\t%e\t%i\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],P[i].Vel[0],P[i].Vel[1],P[i].Vel[2],P[i].Mass,Id[i]);  
 }
 fclose(outG);
 fclose(outD);
}


