/* 

   DGSOL solves the distance geometry problem of finding positions
   x_1,...,x_n for atoms that satisfy the distance contraints

        l_{i,j} <= || x_i - x_j || <= u_{i,j}

   See the INFO file in the docs directory for additional information.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int main (argc, argv)

  int argc;
  char *argv[];

  {

   FILE *input, *output, *summary;
   char namein[200], nameout[200], namesum[200];

   int numAtoms = 0, numSols = 0, numDs = 0, errInfo = 0; 
   int lintWkSp, ldblWkSp, *iInd, *jInd, *intWkSp;
   int i, k, l, m, n, c, seed, mod, sol, myid, nprocs, msgtag;

   double u, v, *lbd, *ubd, *xOpt, *yOpt, *fOpt, *dErr, *dblWkSp;
   MPI_Status stat;


   /*
   MPI initialization.
   */

   MPI_Init (&argc, &argv);

   MPI_Comm_rank (MPI_COMM_WORLD, &myid);
   MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

   if (myid == 0) 
   
   {

   /*
   Process command line arguments.
   */

   if (argc > 7) {
      fprintf (stderr, "Too many arguments!\n");
      errInfo = 1; goto Error_Checking;
      }

   argc --; argv ++;
   while (argc > 0 && (*argv)[0] == '-') {
     c = (*argv)[1];
     switch (c) {
        case 'h':
	  fprintf (stderr,
          "Usage: dgsol [-a] [-d] [-s] [input output summary]\n");
	  errInfo = 2; goto Error_Checking;
        case 'a':
          sscanf (&(*argv)[2], "%d", &numAtoms);
          break;
        case 'd':
          sscanf (&(*argv)[2], "%d", &numDs);
          break;
        case 's':
          sscanf (&(*argv)[2], "%d", &numSols);
          break;
        default:
          argc = -1;
        }
     argc --; argv ++;
     }

   sprintf (namein, "%s", " ");
   if (argc > 0) {
      sprintf (namein, "%s", argv[0]);
      argc --; argv ++;
      }
  
   sprintf (nameout, "%s", " ");
   if (argc > 0) {
      sprintf (nameout, "%s", argv[0]);
      argc --; argv ++;
      }
  
   sprintf (namesum, "%s", " ");
   if (argc > 0) {
      sprintf (namesum, "%s", argv[0]);
      argc --; argv ++;
      }

   if (argc != 0) {
      fprintf (stderr, "Usage: dgsol [-a] [-d] [-s] [input output summary]\n");
      errInfo = 3; goto Error_Checking;
      }

   /*
   Process distance. Default in and out files are dg.data and dg.sol.
   */

   if (namein[0] == ' ') {
      sprintf (namein, "%s", "dg.data");
      }

   if (nameout[0] == ' ') {
      sprintf (nameout, "%s", "dg.sol");
      }

   if (namesum[0] == ' ') {
      sprintf (namesum, "%s", "dg.sum");
      } 

   if (!(input = fopen (namein, "r"))) {
      fprintf (stderr, "Cannot open input file! \n");
      errInfo = 4; goto Error_Checking;
      }
  
   if (!(output = fopen (nameout, "w"))) {
      fprintf (stderr, "Cannot open output file! \n");
      errInfo = 5; goto Error_Checking;
      }
 
   if (!(summary = fopen (namesum, "w"))) {
      fprintf (stderr, "Cannot open summary file! \n");
      errInfo = 6; goto Error_Checking;
      }

   fprintf (summary,
   "\n===================  PERFORMANCE SUMMARY  =====================\n");
   fprintf (summary,
   "\n n_atoms   n_dist     f_err     derr_min    derr_avg    derr_max\n\n");

   /*
   Process the distance data.
   */ 

   printf ("\nProcessing distance data\n");

   if (!numDs) {
      for (numDs = 0; EOF != fscanf 
          (input, "%d %d %lf %lf \n", &k, &l, &u, &v); numDs++); 
      rewind (input);
      }

   }

   Error_Checking:	;

   MPI_Bcast (&errInfo, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if (errInfo > 0) goto End_of_Program;
      
   MPI_Bcast (&numDs, 1, MPI_INT, 0, MPI_COMM_WORLD);

   iInd = (int *) malloc (numDs*sizeof(int)); 
   jInd = (int *) malloc (numDs*sizeof(int)); 

   lbd = (double *) malloc (numDs*sizeof(double)); 
   ubd = (double *) malloc (numDs*sizeof(double)); 

   if (myid == 0) 
   
   {
  
   for (k = 0; k < numDs; k++)
       fscanf (input,"%d %d %lf %lf \n",&iInd[k],&jInd[k],&lbd[k],&ubd[k]);
   
   if (!numAtoms) 
      for (k = 0; k < numDs; k ++) {
          if (numAtoms < iInd[k]) numAtoms = iInd[k];
          if (numAtoms < jInd[k]) numAtoms = jInd[k];
          }

   if (!numSols) numSols = nprocs;

   }

   MPI_Bcast (&numAtoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&numSols, 1, MPI_INT, 0, MPI_COMM_WORLD);
   
   MPI_Bcast (iInd, numDs, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (jInd, numDs, MPI_INT, 0, MPI_COMM_WORLD);

   MPI_Bcast (lbd, numDs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast (ubd, numDs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   for (k = 0; k < numDs; k ++) {
       lbd [k] = lbd [k] * lbd [k];
       ubd [k] = ubd [k] * ubd [k];
       }

   /*
   Allocate memory.
   */

   n = 3*numAtoms;
   lintWkSp = numDs;
   ldblWkSp = 30*n;

   fOpt = (double *) malloc (sizeof(double));
   dErr = (double *) malloc (3*sizeof(double));
   xOpt = (double *) malloc (n*sizeof(double));
   yOpt = (double *) malloc (n*sizeof(double));

   intWkSp = (int *) malloc (lintWkSp*sizeof(int)); 
   dblWkSp = (double *) malloc (ldblWkSp*sizeof(double));

   /*
   Call the DG optimizer to obtain a set of solutions.
   */

   if (myid == 0) {
      seed = 1;
      printf ("\nStarting optimization with DGSOL\n");
      }

   mod = numSols % nprocs;
   sol = numSols / nprocs;

   for (k = 0; k < sol || (k == sol && myid < mod); k++) {
     
     m = sol ? k / sol : 1;

     /*
     Call the DG optimizer.
     */

     msgtag = 100;

     if (myid > 0) 

        MPI_Recv (xOpt,n,MPI_DOUBLE,0,msgtag+myid,MPI_COMM_WORLD,&stat);

     else {

        dginitial (&numAtoms,xOpt,&numDs,iInd,jInd,lbd,ubd,&seed,intWkSp);

        for (l = 1; l < m * (mod - nprocs) + nprocs; l ++) {
            dginitial (&numAtoms,yOpt,&numDs,iInd,jInd,lbd,ubd,&seed,intWkSp);
            MPI_Send (yOpt,n,MPI_DOUBLE,l,msgtag+l,MPI_COMM_WORLD);
            }

        }

     dgopt (&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,intWkSp,dblWkSp);
     dgerr (&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,dErr);

     /*
     Collect solutions.
     */

     msgtag = 200;

     if (myid > 0) {

        MPI_Send (xOpt,n,MPI_DOUBLE,0,1*msgtag+myid,MPI_COMM_WORLD);
        MPI_Send (fOpt,1,MPI_DOUBLE,0,2*msgtag+myid,MPI_COMM_WORLD);
        MPI_Send (dErr,3,MPI_DOUBLE,0,3*msgtag+myid,MPI_COMM_WORLD);

        }

     else {

        for (l = 0; l < m * (mod - nprocs) + nprocs; l ++) {

            if (l > 0) {

               MPI_Recv (xOpt,n,MPI_DOUBLE,l,1*msgtag+l,MPI_COMM_WORLD,&stat);
               MPI_Recv (fOpt,1,MPI_DOUBLE,l,2*msgtag+l,MPI_COMM_WORLD,&stat);
               MPI_Recv (dErr,3,MPI_DOUBLE,l,3*msgtag+l,MPI_COMM_WORLD,&stat);

               }

            for (i = 0; i < numAtoms; i ++) {
                fprintf (output, 
                "%16.8e %16.8e %16.8e \n",xOpt[3*i],xOpt[3*i+1],xOpt[3*i+2]);
                }

            fprintf (output, "\n %20.12e \n\n", *fOpt);
            
            fprintf (summary,
            "%6d  %8d  %10.2e  %10.2e  %10.2e  %10.2e \n",
            numAtoms, numDs, *fOpt, dErr[0], dErr[1], dErr[2]);


            printf ("\nOptimization %d - Largest distance error is %12.3e\n",
                   k*nprocs+l+1, dErr[2]);

            }

        }

     }

   if (myid == 0) { 
      printf ("\nEnd of optimization with DGSOL\n");
      }

   /*
   Clean up.
   */


   free (iInd); free (jInd);
   free (lbd); free (ubd); 
   free (fOpt); free(dErr); free (xOpt); 
   free (intWkSp); free (dblWkSp); 

   End_of_Program:	;

   if (myid == 0) { 
      fclose (input); fclose (output); fclose (summary);
      }

   MPI_Finalize ();
  
   }

