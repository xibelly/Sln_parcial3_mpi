/* Xibelly Eliseth Mosquera Escobar - Lunes 21 de Nov de 2016

   En este codigo resolvemos el problema gravitacional de N-cuerpos
   
   Esto corresponde a calcular la fuerza gravitacional que sienten las N particulas, 
   haciendo uso de la siente formula

   F = -G sum_{i!=j}^{N} [-G mi * mj (ri - rj) ] / |ri - rj|^3

 */

/*
Diseño
1) Leemos el conjunto de  datos y almacenamos en una estructura 

2) se calculan las fuerzas
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<time.h>

#define G 1

//variables globales
clock_t tini, tend, tacum;
double cpu_time_used;

//Estructuras

struct DATA
{
  double *x;           /*coordenadas particulas*/
  double *y;
  double *z;
  double *m;           /*masa particulas*/
  double *Fx;
  double *Fy;
  double *Fz;
  
};
struct DATA data;


//Sub-rutinas

#include "distancia.c"
#include "fuerza.c"

/*----------------------------------------------------------------PROGRAMA PRINCIPAL*/

int main(int argc, char **argv)
{
  int i,j,k;
  int nread;
  int Npart;
  int Nprocess = 1;

  double xi, yi, zi, mi;
  double xj, yj, zj, mj;
  double dist;
  
  char buff[200];
  char *file;

  FILE *read = NULL;  
  FILE *write = NULL;
  FILE *write2 = NULL;
  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("mpiexec -n 1 %s FILE NUMBER_OF_PARTCLES \n",argv[0]);
      exit(0);  
    }

  //Carga de parametros

  file = argv[1];
  Npart = atof(argv[2]);
     
  if(Npart < 0)
    printf("ERROR: NO PARTICLES\n");

  data.x = (double *) malloc(Npart *sizeof(double));  /*Global Particles*/
  data.y = (double *) malloc(Npart *sizeof(double));
  data.z = (double *) malloc(Npart *sizeof(double));
  data.m = (double *) malloc(Npart *sizeof(double));
  data.Fx = (double *) malloc(Npart *sizeof(double));
  data.Fy = (double *) malloc(Npart *sizeof(double));
  data.Fz = (double *) malloc(Npart *sizeof(double));

   
  printf("NUMBER OF TOTAL PARTICLES IS=%d\n",  Npart);
  printf("----------------------------------------------------------------------------\n");

  
  printf("READING FILE... \n");
  
  read = fopen(file,"r");
  
  for(i=0; i<Npart; i++)
    {
      nread = fscanf(read,"%lf %lf %lf %lf",&data.x[i],&data.y[i],&data.z[i],&data.m[i]);
      
    }
  
  sprintf(buff,"forces_seire_Nparts%d.dat",Npart);
  write = fopen(buff,"w");

  tini = clock();
  
  for(i=0; i<Npart; i++)
    {
      
      for(j=0; j<Npart; j++)
	{
	  if(i != j)
	    {
	      xi = data.x[i];
	      yi = data.y[i];
	      zi = data.z[i];
	      
	      xj = data.x[j];
	      yj = data.y[j];
	      zj = data.z[j];
	      
	      mi = data.m[i];
	      mj = data.m[j];
	      
	      dist = distance(xi,yi,zi,xj,yj,zj);
	      
	      data.Fx[i] += fuerza(G, mi,mj, xi, xj,dist);
		  
	      data.Fy[i] += fuerza(G, mi,mj, yi, yj,dist);
	      
	      data.Fz[i] += fuerza(G, mi,mj, zi, zj,dist);			  
	      
	      fprintf(write,"%lf %e %e %e %e %e %e \n",dist,xi,yi,zi,data.Fx[i],data.Fy[i],data.Fz[i]);
	    }
	}
    }
  
  
  printf("DATA HAS BEEN READING\n");
  printf("THE TOTAL FORCES HAS BEEN CALCULATED\n");
  
  tend = clock();
  
  cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
  printf("CPU TIME USED TO CALCULATE THE FORCE IS: %g\n",cpu_time_used);
  
  sprintf(buff,"execution_time_serie_Nparts%d",Npart);
  write2 = fopen(buff,"w");

  fprintf(write2,"%lf %d",cpu_time_used,Nprocess);

  return 0;
}


  
