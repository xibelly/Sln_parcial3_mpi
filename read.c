/* Xibelly Eliseth Mosquera Escobar -  Nov de 2016

   Este programa lee el archivo origianl con el conjunto total de particulas
   -Ntotal, lee las coordenadas x,y,z y la masa asociada a cada una de ellas, 
   luego selecciona de manera aleatoria un subconjunto de Nsub particulas, 
   este ultimo parametro se recive por linea de comandos.

   

 */



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_errno.h>

//variables globales
long int Ntotal;
int Nsub;
int *value;

//sub-rutinas

int random_taus()
{
  
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

   
  rng = gsl_rng_alloc (gsl_rng_taus);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<Nsub; i++)
    {
      value[i] = gsl_rng_uniform_int(rng,Ntotal); //genera numeros entre 0 y Ntotal-1
    }
 
  gsl_rng_free (rng);

  printf("STATE OF RANDOM NUMBER GENERATION BY 'GFSR4': SUCCESS\n");
 
  return 0;
  
}


int main(int argc, char **argv)
{
  int i,j,nread,step;

  double *x, *y, *z, *M;
  double posx, posy, posz, masa;
  
  FILE *read = NULL;
  FILE *write = NULL;

  
  char *file;
  char buff[200];
    
  if(argc != 4)
    {
      printf("ERROR--> use as:\n");
      printf("%s FILE N_part N_sub_part \n",argv[0]);
      exit(0);  
    }

  //Carga de parametros

  file = argv[1];
  Ntotal = atoi(argv[2]);
  Nsub = atoi(argv[3]);

  value = (int *)malloc(Nsub *sizeof(int));

  M = (double *)malloc(Ntotal *sizeof(double));

  x = (double *) malloc(Ntotal *sizeof(double));

  y = (double *) malloc(Ntotal *sizeof(double));

  z = (double *) malloc(Ntotal *sizeof(double));
  
  //lectura del archivo de datos

  read = fopen(file,"r");
  
  for(i=0; i<Ntotal; i++)
    {
      nread = fscanf(read,"%lf %lf %lf %lf", &x[i], &y[i], &z[i], &M[i]);
      
    }

  printf("THE READING HAS BEEN SUCCESS\n");
  
  //seleccion del subconjunto de datos

  random_taus();

  sprintf(buff,"data_parts_%d.dat",Nsub);

  write = fopen(buff,"w");
  
  for(i=0; i<Nsub; i++)
    {
      step = value[i];
               	    
      posx = x[step];
      
      posy = y[step];
      
      posz = z[step];
      
      masa = M[step];
      
      fprintf(write,"%lf %lf %lf %lf\n",posx, posy, posz, masa);
      //printf("%lf %lf %lf %lf %d\n",posx, posy, posz, masa, step);
    }
      
  printf("THE SUB-SAMPLE HAS BEEN CREATED WITH SUCCESS\n");
    
  return 0;
}

 
