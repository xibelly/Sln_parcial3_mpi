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
#include<time.h>

#define G 1


//Variables globales

int Number_of_process;
int Npart;
int Parts_domain;
int local_parts;
clock_t tini, tend, tacum;
double cpu_time_used;

FILE *out1 = NULL;
FILE *out2 = NULL;
FILE *write = NULL;

//Estructuras

struct particles   
{
  double x;           /*coordenadas particulas*/
  double y;
  double z;
  double m;  /*masa particulas*/
  double Fx,Fy,Fz;

};  
struct particles *Part_temp, *Part_recv,*Part_local0,*Part_local;

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

#include "distancia.c"
#include "fuerza.c"

/*----------------------------------------------------------------PROGRAMA PRINCIPAL*/

int main(int argc, char **argv)
{
  int i,j,k,l,p;
  int ii, jj;
  int itask=0;
  int vecino, send, recv;
  int nread;
  int istar, iend;
  int err, dest, remit;
  int error, eclass, len;
  int task;
  int tag1 = 1;
  int tag2 = 2;
    
  int Residuo;
  int Domain_size;
  double recvbuf=0.0;
  double xi, yi, zi, mi;
  double xj, yj, zj, mj;
  double dist;
  

  char estring[MPI_MAX_ERROR_STRING];
  char *file;
  char buff[200];

  MPI_Datatype matrix_type;
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
      printf("mpiexec -n Number_process %s FILE NUMBER_OF_PARTCLES \n",argv[0]);
      exit(0);  
    }

  //Carga de parametros

  //Tamaño del Dominio -particulas en el dominio-
  file = argv[1];
  Npart = atof(argv[2]);
  Parts_domain = (int) ceil((1.0*Npart) /Number_of_process);
  Residuo = Npart % Number_of_process;

    
  if(Npart < 0)
    printf("ERROR: NO PARTICLES\n");

  data.x = (double *) malloc(Npart *sizeof(double));  /*Global Particles*/
  data.y = (double *) malloc(Npart *sizeof(double));
  data.z = (double *) malloc(Npart *sizeof(double));
  data.m = (double *) malloc(Npart *sizeof(double));
  data.id = (int *) malloc(Npart *sizeof(int));

   
  if(task == 0)
    {
      printf("NUMBER OF TOTAL PARTICLES IS=%d\n",  Npart);
      printf("PARTICLES BY DOMAIN: %d\n",Parts_domain);
      printf("TOTAL DOMAINS=%d\n", Number_of_process);
      printf("EXTRA PARTICLES THAT NEED TO BE ALLOCATED : %d\n",Residuo);
      
      printf("----------------------------------------------------------------------------\n");
    }

  MPI_Barrier(MPI_COMM_WORLD);


  if(task == 0)
    {
      printf("READING FILE... \n");
      
      read = fopen(file,"r");

      for(i=0; i<Npart; i++)
	{
	  nread = fscanf(read,"%lf %lf %lf %lf",&data.x[i],&data.y[i],&data.z[i],&data.m[i]);
	  //printf("%lf %lf %lf %lf\n",data.x[i],data.y[i],data.z[i],data.m[i]);
	}
      fclose(read);

      Part_local0 = (struct particles *) malloc(Parts_domain *sizeof(struct particles));
          
      for(i=0; i<Parts_domain; i++)
	{
	  Part_local0[i].x =data.x[i];
	  Part_local0[i].y =data.y[i];
	  Part_local0[i].z =data.z[i];
	  Part_local0[i].m =data.m[i];

	  //printf("%lf %lf %lf %lf\n",Part_local[task].x[i],Part_local[task].y[i],Part_local[task].z[i],Part_local[task].m [i]);
	}

      sprintf(buff,"forces_task%d_Nparts%d.dat",task,Npart);
      out2 = fopen(buff,"w");

      tini = clock();
      
      for(i=0; i<Parts_domain; i++)
	{
	  //printf("%lf %lf %lf %lf\n",Part_local0[i].x,Part_local0[i].y,Part_local0[i].z,Part_local0[i].m);
	  for(j=0; j<Parts_domain; j++)
	    {
	      if(i != j)
		{
		  xi = Part_local0[i].x;
		  yi = Part_local0[i].y;
		  zi = Part_local0[i].z;
		  
		  xj = Part_local0[j].x;
		  yj = Part_local0[j].y;
		  zj = Part_local0[j].z;
		  
		  mi = Part_local0[i].m;
		  mj = Part_local0[j].m;
		  
		  dist = distance(xi,yi,zi,xj,yj,zj);
	      
		  Part_local0[i].Fx += fuerza(G, mi,mj, xi, xj,dist);//-(G * mi *mj* (xi-xj) )/(dist*dist*dist);
		  
		  Part_local0[i].Fy += fuerza(G, mi,mj, yi, yj,dist);//-(G * mi *mj* (yi-yj) )/(dist*dist*dist);
		  
		  Part_local0[i].Fz += fuerza(G, mi,mj, zi, zj,dist);//-(G * mi *mj* (zi-zj) )/(dist*dist*dist);			  
		  
		  fprintf(out2,"%lf %e %e %e %e %e %e \n",dist,xi,yi,zi,Part_local0[i].Fx,Part_local0[i].Fy,Part_local0[i].Fz);
		}
	    }
	}
      tend = clock();
      
      printf("DATA HAS BEEN READING\n");
      printf("task %d HAS BEEN CALCULATED ITS INTERNAL FORCES\n",task);

      cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
      printf("CPU TIME USED TO CALCULATE THE FORCE BY taks %d IS: %g\n",task,cpu_time_used);

      sprintf(buff,"execution_time_task%d_Nparts%d",task,Parts_domain);
      write = fopen(buff,"w");

      fprintf(write,"%lf %d %d",cpu_time_used,task,Number_of_process);
      
    }//close task0

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0)
    {
      printf("INITIALIZING THE DOMAIN DESCOMPOSITION \n");
      
      printf("----------------------------------------------------------------------------\n");
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0)
    {
            
      //------------------------------------------------------------------------Descomposicion de dominio 
      int size_dom;
      char buffer[200];
      int pos=0;
      
      
      for(itask=1; itask<Number_of_process; itask++)
	{
	  if(Residuo == 0)
	    {
	      size_dom = Parts_domain;
	      Part_temp = (struct particles *) malloc(size_dom *sizeof(struct particles)); 
	      
	     
	      istar = itask*size_dom ;
	      iend = (itask+1)* size_dom;

	      i=0;
	      for(j=istar; j<iend; j++)
		{
		  Part_temp[i].x = data.x[j];
		  Part_temp[i].y = data.y[j];
		  Part_temp[i].z = data.z[j];
		  Part_temp[i].m = data.m[j];
		  
		  //printf("%lf %lf %lf %lf %d\n",Part_temp.x[j],Part_temp.y[j],Part_temp.z[j],Part_temp.m[j],itask);
		  i+=1;
		}
	      //Se envian las particulas a cada dominio
	      
	      printf("Sending to task %d\n", itask);fflush(stdout);
	      printf("----------------------------------------------------------------------------\n");
	      MPI_Send(&size_dom, 1, MPI_INT, itask, tag1, MPI_COMM_WORLD);
	      
	      MPI_Send(Part_temp, size_dom *sizeof(struct particles), MPI_BYTE, itask, tag2, MPI_COMM_WORLD);
	      	      
	    }//close if1
	  	  
	  else if(itask < Residuo)
	    {
	      size_dom = Parts_domain;
	      Part_temp = (struct particles *) malloc(size_dom *sizeof(struct particles)); 
	      
	      istar = itask*size_dom ;
	      iend = (itask+1)* size_dom;
	      
	      i=0;
	      for(j=istar; j<iend; j++)
		{
		  Part_temp[i].x = data.x[j];
		  Part_temp[i].y = data.y[j];
		  Part_temp[i].z = data.z[j];
		  Part_temp[i].m = data.m[j];
		  
		  //printf("%lf %lf %lf %lf %d\n",Part_temp.x[j],Part_temp.y[j],Part_temp.z[j],Part_temp.m[j],itask);
		  i+=1;
		}
	      //Se envian las particulas a cada dominio
	      
	      printf("Sending to task %d\n", itask);fflush(stdout);
	      printf("----------------------------------------------------------------------------\n");
	      MPI_Send(&size_dom, 1, MPI_INT, itask, tag1, MPI_COMM_WORLD);
	      
	      MPI_Send(Part_temp, size_dom *sizeof(struct particles), MPI_BYTE, itask, tag2, MPI_COMM_WORLD);
	    }//close if2
      
	  else if(itask >= Residuo)
	    {
	      size_dom = Parts_domain - 1;

	      istar = itask*size_dom ;
	      iend = (itask+1)* size_dom;
	      
	      i=0;
	      for(j=istar; j<iend; j++)
		{
		  Part_temp[i].x = data.x[j];
		  Part_temp[i].y = data.y[j];
		  Part_temp[i].z = data.z[j];
		  Part_temp[i].m = data.m[j];
		  
		  //printf("%lf %lf %lf %lf %d\n",Part_temp.x[j],Part_temp.y[j],Part_temp.z[j],Part_temp.m[j],itask);
		  i+=1;
		}
	      //Se envian las particulas a cada dominio
	      
	      printf("Sending to task %d\n", itask);fflush(stdout);
	      printf("----------------------------------------------------------------------------\n");
	      MPI_Send(&size_dom, 1, MPI_INT, itask, tag1, MPI_COMM_WORLD);
	      
	      MPI_Send(Part_temp, size_dom *sizeof(struct particles), MPI_BYTE, itask, tag2, MPI_COMM_WORLD);

	    }//close if3

	}//close domains

    }//total close 


    
   //Se reciben las Particulas de cada dominio

  if(task !=0)
    {
      MPI_Recv(&local_parts, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);
      printf("Sent value from task %d to task %d is %d\n", 0, task, local_parts);
      printf("Size domain in bits is: %lu\n",sizeof(local_parts)*sizeof(struct particles));
      fflush(stdout);
 
      Part_recv = (struct particles *) malloc(local_parts *sizeof(struct particles)); 

      error = MPI_Recv(Part_recv,local_parts *sizeof(struct particles), MPI_BYTE, 0, tag2, MPI_COMM_WORLD, &status);
      
      MPI_Error_class(error, &eclass);
      MPI_Error_string(error, estring, &len);
      printf("Error %d:%s\n",eclass, estring);


      out1 = fopen("received_data.dat","a");
      for(i=1; i<Number_of_process; i++)
	{
	  
	  if(task == i)
	    {	      
	      for(j=0; j<local_parts; j++)
		{
		  fprintf(out1,"%lf %lf %lf %lf %d\n",Part_recv[j].x,Part_recv[j].y,Part_recv[j].z,Part_recv[j].m,task);
		  //printf("%lf %lf %lf %lf %d\n",Part_recv[j].x,Part_recv[j].y,Part_recv[j].z,Part_recv[j].m,task);
		}
	    } 
	  
	}
      fclose(out1);

      //Alocacion las particulas locales   
      
     
      Part_local = (struct particles *) malloc(local_parts *sizeof(struct particles));

      //Copia de local de las particulas en cada dominio

      for(i=1; i<Number_of_process; i++)
	{
	  
	  if(task == i)
	    {
	      //l=local_parts;
	      for(j=0; j<local_parts; j++)
		{
		  Part_local[j].x = Part_recv[j].x;
		  Part_local[j].y = Part_recv[j].y;
		  Part_local[j].z = Part_recv[j].z;
		  Part_local[j].m = Part_recv[j].m;

		  //printf("%lf %lf %lf %lf %d\n",Part_local[j].x,Part_local[j].y,Part_local[j].z,Part_local[j].m,task);
		  //l +=1; 
		}
	    } 
	  
	}
      tini = clock();
      for(i=1; i<Number_of_process; i++)
	{
	  
	  if(task == i)
	    {
	      sprintf(buff,"forces_task%d_Nparts%d.dat",task,Npart);
	      out2 = fopen(buff,"w");
		      
	      for(ii=0; ii<local_parts; ii++)
		{
		  
		  for(jj=0; jj<local_parts; jj++)
		    {
		      if(ii != jj)
			{
			  xi = Part_local[ii].x;
			  yi = Part_local[ii].y;
			  zi = Part_local[ii].z;
			  
			  xj = Part_local[jj].x;
			  yj = Part_local[jj].y;
			  zj = Part_local[jj].z;
			  
			  mi = Part_local[ii].m;
			  mj = Part_local[jj].m;
			  			  
			  dist =distance(xi,yi,zi,xj,yj,zj);
			  
			  Part_local[ii].Fx += fuerza(G, mi,mj, xi, xj,dist);//-(G * mi *mj* (xi-xj) )/(dist*dist*dist);
			  
			  Part_local[ii].Fy += fuerza(G, mi,mj, yi, yj,dist);//-(G * mi *mj* (yi-yj) )/(dist*dist*dist);
			  
			  Part_local[ii].Fz += fuerza(G, mi,mj, zi, zj,dist);//-(G * mi *mj* (zi-zj) )/(dist*dist*dist);			  
		  
			  fprintf(out2,"%lf %e %e %e %e %e %e\n",dist,xi,yi,zi,Part_local[ii].Fx,Part_local[ii].Fy,Part_local[ii].Fz);
			  //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",dist,xi,yi,zi,xj,yj,zj,mi,mj,task);
			}
		    }
		  
		}
	      
	    }
	}

      printf("task %d HAS BEEN CALCULATED ITS INTERNAL FORCES\n",task);

      tend = clock();
      
      cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
  
      printf("CPU TIME USED TO CALCULATE THE FORCE BY task %d IS: %g\n",task,cpu_time_used);

      sprintf(buff,"execution_time_task%d_Nparts%d",task,local_parts);
      write = fopen(buff,"w");

      fprintf(write,"%lf %d %d",cpu_time_used,task,Number_of_process);

    }//close task!=0

  
  //------------------------------------------------------------------------CALCULO DE FUERZAS EN CADA GRUPO

  if(task == (Number_of_process- 1))
    {
      printf("Data is going to rotate between the tasks\n");
      printf("----------------------------------------------------------------------------\n");
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  int left = (task - 1 + Number_of_process) % Number_of_process;
  int right = (task + 1) % Number_of_process;


  MPI_Barrier(MPI_COMM_WORLD);
  err = MPI_Finalize();
}


  
