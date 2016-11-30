/* Xibelly Eliseth Mosquera Escobar - Lunes 21 de Nov de 2016

   En este codigo resolvemos el problema gravitacional de N-cuerpos
   
   Esto corresponde a calcular la fuerza gravitacional que sienten las N particulas, 
   haciendo uso de la siente formula

   F = -G sum_{i!=j}^{N} [-G mi * mj (ri - rj) ] / |ri - rj|^3

 */

/*
Diseño

1) Con el programa read.c generamos un sub-muestreo del conjunto total de particulas

2) Leemos dicho sub-conjunto de datos y almacenamos en una estructura 

3) Definimos el dominio computacional para cada proceso, como:

 Domain_size = ceil (Nparts / Number_of_process)

Nota: En el caso en que Nparts no sea divisible por Number_of_process, el ultimo proceso
      sera sobre cargado con las particulas que no son asignadas.

El incio y el fin del dominio se define como:

  istar = j * Domain_size;
  iend = (j+1) * Domain_size ;  con j =0,..,Number_of_process

4) En el interior de cada dominio se hacen dos copias de las particulas una local y otra 
 que se rotara entre los procesadores para el calculo de la fuerza. Con la copia local
 se procede a calcular el campo de fuerza entre las particulas del mismo grupo.

5) Para el calculo del campo entre particulas de diferentes grupos, se hara uso de la comunicacion
 en anillo en forma dextrogira, rotanto entre los  procesadores la copia rotante de las 
 coordenadas de las particulas.

Nota: Todo el trabajos e hace bajo el proceso raiz el cual se encargara de recibir los calculos
que efectuan los demas procesos. Este a su vez imprimira en disco el campo de fuerza total entre 
las particulas.

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

#define G 1

//Variables globales

int Number_of_process;
int Npart;
int Parts_domain;


FILE *out1 = NULL;
FILE *out2 = NULL;

//Estructuras

struct particles   
{
  double *x;           /*coordenadas particulas*/
  double *y;
  double *z;
  double *m;           /*masa particulas*/

  double *Fx, *Fy, *Fz; /*Conponentes de la fuerza*/

};  
struct particles Part_rot, *Part_local;

struct DATA
{
  double *x;           /*coordenadas particulas*/
  double *y;
  double *z;
  double *m;           /*masa particulas*/
  int *id;
};
struct DATA data;


//Sub-rutinas

//funcion que calcula la distancia entre las particulas
double distance(double xi, double yi, double zi, double xj, double yj, double zj)
{

  double dist, dx, dy, dz;
  
  dx = (xi-xj)*(xi-xj);
  dy = (yi-yj)*(yi-yj);
  dz = (zi-zj)*(zi-zj);
  
  dist = sqrt(dx + dy + dz); 
  
  return dist;

}


/*----------------------------------------------------------------PROGRAMA PRINCIPAL*/

int main(int argc, char **argv)
{
  int i,j,k,l;
  int ii, jj;
  int nread;
  int istar, iend;
  int err, dest, remit;
  int error, eclass, len;
  int task;

  int Residuo;
  int Domain_size;
  double recvbuf=0.0;
  double xi, yi, zi, mi;
  double xj, yj, zj, mj;
  double dist;
  

  char estring[MPI_MAX_ERROR_STRING];
  char *file;
  char buff[200];
  
  MPI_Status status;

  FILE *read = NULL;  
  
  
  //Se inicializa MPI
  
  err = MPI_Init(&argc, &argv); 
  
  MPI_Comm_size(MPI_COMM_WORLD, &Number_of_process);
  MPI_Comm_rank(MPI_COMM_WORLD, &task); 
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("%s FILE NUMBER_OF_PARTCLES \n",argv[0]);
      exit(0);  
    }

  //Carga de parametros

  file = argv[1];
  Npart = atof(argv[2]);
  
  //Tamaño del Dominio -particulas en el dominio-

  Parts_domain = ceil(Npart / Number_of_process);
  
  if(Npart < 0)
    printf("ERROR: NO PARTICLES\n");

  data.x = (double *) malloc(Npart *sizeof(double));  /*Global Particles*/
  data.y = (double *) malloc(Npart *sizeof(double));
  data.z = (double *) malloc(Npart *sizeof(double));
  data.m = (double *) malloc(Npart *sizeof(double));
  data.id = (int *) malloc(Npart *sizeof(int));
  
  Part_rot.x = (double *) malloc( Parts_domain *sizeof(double));  /*Rotate Particles*/
  Part_rot.y = (double *) malloc( Parts_domain *sizeof(double));
  Part_rot.z = (double *) malloc( Parts_domain *sizeof(double));
  Part_rot.m = (double *) malloc( Parts_domain *sizeof(double));


  Part_local = (struct particles *) malloc( Number_of_process *sizeof(struct particles));  /*Local Particles*/
  
  for(i=0; i<Number_of_process; i++)
    {
      Part_local[i].x = (double *) malloc( Parts_domain *sizeof(double)); 
      Part_local[i].y = (double *) malloc( Parts_domain *sizeof(double));
      Part_local[i].z = (double *) malloc( Parts_domain *sizeof(double));
      Part_local[i].m = (double *) malloc( Parts_domain *sizeof(double));

      Part_local[i].Fx = (double *) malloc( Parts_domain *sizeof(double));

      Part_local[i].Fy = (double *) malloc( Parts_domain *sizeof(double));

      Part_local[i].Fz = (double *) malloc( Parts_domain *sizeof(double)); 
    }
   
 
 
  //-------------------------------------------------------------------------Calculo campo de fuerza por comunicacion en anillo

  
   
    
  if(task == 0)
    {
      Residuo = (Npart % Number_of_process);//indica el # de particulas que no son asignadas cuando Npart no es divisible sobre  Number_of_process.
      printf("THE RESIDUE IS: %d\n",Residuo);
      
      printf("PARTICLES BY TASK: %d\n",Parts_domain);
      //-------------------------------------------------------------------------Calculo campo de fuerza CASO 1

      if(Residuo == 0)//Todos los task se cargan equitativamente
	{            
	  read = fopen(file,"r");
	  
	  for(k=0; k<Npart; k++)//Se leen las particulas
	    {
	      nread = fscanf(read,"%lf %lf %lf %lf", &data.x[k], &data.y[k], &data.z[k], &data.m[k]);
	      data.id[k] = k;
	    }
	  
	  for(l=0; l<Number_of_process; l++)//Se establece cada Dominio
	    {
	      
	      istar = l * (int) ceil(Parts_domain);     //DOMAIN
	      iend = (l + 1) * (int) ceil( Parts_domain);
	      
	      for(i=istar; i<iend; i++)//Se distribuyen las parts en cada Dominio
		{
		  Part_local[l].x[i] = data.x[i];
		  Part_local[l].y[i] = data.y[i];
		  Part_local[l].z[i] = data.z[i];
		  Part_local[l].m[i] = data.m[i];
		  
		}

	      //Calculo campo de fuerza entre parts de un mismo dominio
	      
	      sprintf(buff,"field_#process%d_task%d.dat",Number_of_process,l);
	      out1 = fopen(buff,"w");
	      
	      for(ii=istar; ii< iend; ii++)
		{
		  for(jj=istar; jj< iend; jj++)
		    {
		      if(jj != ii)
			{
			  xi = Part_local[l].x[ii];
			  yi = Part_local[l].y[ii];
			  zi = Part_local[l].z[ii];
			  
			  xj = Part_local[l].x[jj];
			  yj = Part_local[l].y[jj];
			  zj = Part_local[l].z[jj];

			  mi = Part_local[l].m[ii];
			  mj = Part_local[l].m[jj];

			  dist = distance(xi,yi,zi,xj,yj,zj);
			  
			  Part_local[l].Fx[ii] += -(G * mi *mj* (xi-xj) )/(dist*dist*dist);

			  Part_local[l].Fy[ii] += -(G * mi *mj* (yi-yj) )/(dist*dist*dist);

			  Part_local[l].Fz[ii] += -(G * mi *mj* (zi-zj) )/(dist*dist*dist);			  

			  fprintf(out1,"%lf %lf %lf %lf\n",dist,Part_local[l].Fx[ii],Part_local[l].Fy[ii],Part_local[l].Fz[ii]);
			  
			}
		      
		    }
		  
		}	    
	      
	    }//close cada dominio

	  printf("CASE 1: SUCESS\n");

	  fclose(out1);
	}
      

      //-------------------------------------------------------------------------Calculo campo de fuerza CASO 2

      if(Residuo != 0)//Se sobre carga el ultimo proceso con las particulas restante
	{
	  read = fopen(file,"r");
	  
	  for(k=0; k<Npart; k++)//Se leen las particulas
	    {
	      nread = fscanf(read,"%lf %lf %lf %lf", &data.x[k], &data.y[k], &data.z[k], &data.m[k]);
	      data.id[k] = k;
	    }
	  
	  for(l=0; l<Number_of_process; l++)//Se establece cada Dominio
	    {
	      
	      istar = l * (int) ceil(Parts_domain);     //DOMAIN
	      iend = (l + 1) * (int) ceil( Parts_domain);
	      
	      for(i=istar; i<iend; i++)//Se distribuyen las parts en cada Dominio
		{
		  Part_local[l].x[i] = data.x[i];
		  Part_local[l].y[i] = data.y[i];
		  Part_local[l].z[i] = data.z[i];
		  Part_local[l].m[i] = data.m[i];
		  
		}

	      //Calculo campo de fuerza entre parts de un mismo dominio
	      

	      sprintf(buff,"field_#process%d_task%d.dat",Number_of_process,l);
	      out1 = fopen(buff,"w");
	      
	      for(ii=istar; ii< iend; ii++)
		{
		  for(jj=istar; jj< iend; jj++)
		    {
		      if(jj != ii)
			{
			  xi = Part_local[l].x[ii];
			  yi = Part_local[l].y[ii];
			  zi = Part_local[l].z[ii];
			  
			  xj = Part_local[l].x[jj];
			  yj = Part_local[l].y[jj];
			  zj = Part_local[l].z[jj];
			  
			  mi = Part_local[l].m[ii];
			  mj = Part_local[l].m[jj];

			  dist = distance(xi,yi,zi,xj,yj,zj);
			  
			  Part_local[l].Fx[ii] += -(G * mi *mj* (xi-xj) )/(dist*dist*dist);

			  Part_local[l].Fy[ii] += -(G * mi *mj* (yi-yj) )/(dist*dist*dist);

			  Part_local[l].Fz[ii] += -(G * mi *mj* (zi-zj) )/(dist*dist*dist);			  

			  fprintf(out1,"%lf %lf %lf %lf\n",dist,Part_local[l].Fx[ii],Part_local[l].Fy[ii],Part_local[l].Fz[ii]);
			}
		      
		    }
		 
		}

	      //printf("istar:%d iend:%d\n",istar,iend);
	      if(l == (Number_of_process - 1))//Caso en que Nparts no es divisble entre Number_of_process
		{
		  sprintf(buff,"field_#process%d_task%d.dat",Number_of_process,l);

		  out2 = fopen(buff,"w");
		  
		  iend = Npart;

		  //printf("istar:%d iend:%d\n",istar,iend);
		  
		  Domain_size = Residuo + Parts_domain ; 

		  printf("THE TASK %d HAS %d PARTICLES\n",l,Domain_size);

		  //Se alocan las particulas que contendra este proceso

		  Part_local[l].x = (double *) malloc( Domain_size *sizeof(double)); 
		  Part_local[l].y = (double *) malloc( Domain_size *sizeof(double));
		  Part_local[l].z = (double *) malloc( Domain_size *sizeof(double));
		  Part_local[l].m = (double *) malloc( Domain_size *sizeof(double));
		  
		  for(j=istar; j<iend; j++)//Se distribuyen las parts en el Dominio
		    {
		      Part_local[l].x[j] = data.x[j];
		      Part_local[l].y[j] = data.y[j];
		      Part_local[l].z[j] = data.z[j];
		      Part_local[l].m[j] = data.m[j];
		  
		    }
		  
		  //Calculo campo de fuerza entre parts de un mismo dominio
		  for(ii=istar; ii< iend; ii++)
		    {
		      for(jj=istar; jj< iend; jj++)
			{
			  if(jj != ii)
			    {
			      xi = Part_local[l].x[ii];
			      yi = Part_local[l].y[ii];
			      zi = Part_local[l].z[ii];
			      
			      xj = Part_local[l].x[jj];
			      yj = Part_local[l].y[jj];
			      zj = Part_local[l].z[jj];
			      
			      mi = Part_local[l].m[ii];
			      mj = Part_local[l].m[jj];
			      
			      dist = distance(xi,yi,zi,xj,yj,zj);
			  
			      Part_local[l].Fx[ii] += -(G * mi *mj* (xi-xj) )/(dist*dist*dist);

			      Part_local[l].Fy[ii] += -(G * mi *mj* (yi-yj) )/(dist*dist*dist);
			      
			      Part_local[l].Fz[ii] += -(G * mi *mj* (zi-zj) )/(dist*dist*dist);			  
			      
			      fprintf(out2,"%lf %lf %lf %lf\n",dist,Part_local[l].Fx[ii],Part_local[l].Fy[ii],Part_local[l].Fz[ii]);
			    }
			  
			}
		      
		    }
		}//close sobre carga 
	      
	    }//close task      

	  printf("CASE 2: SUCESS\n");

	  fclose(out1);
	  fclose(out2);

	}//close reidui!=0

      

    }//close root process

  
  MPI_Barrier(MPI_COMM_WORLD);
  err = MPI_Finalize();
}

	  


  


