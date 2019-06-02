/* 

   DGSOL solves the distance geometry problem of finding positions
   x_1,...,x_n for atoms that satisfy the distance contraints

        l_{i,j} <= || x_i - x_j || <= u_{i,j}

   See the INFO file in the docs directory for additional information.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>



int main (argc, argv)

  int argc;
  char *argv[];

  {

   FILE *input, *output, *summary;
   char namein[200], nameout[200], namesum[200];

   int numAtoms = 0, numSols = 0, numDs = 0; 
   int lintWkSp, ldblWkSp, *iInd, *jInd, *intWkSp;
   int n, k, l, c, seed;

   double u, v, *lbd, *ubd, *xOpt, *fOpt, *dErr, *dblWkSp;

   /*
   Process command line arguments.
   */

   if (argc > 7) {
      fprintf (stderr, "Too many arguments!\n");
      exit (-1);
      }

   argc --; argv ++;
   while (argc > 0 && (*argv)[0] == '-') {
     c = (*argv)[1];
     switch (c) {
        case 'h':
	  fprintf (stderr, 
	  "Usage: dgsol [-a] [-d] [-s] [input output summary]\n");
	  exit (-1);
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
      exit (-1);
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
      exit (-1);
      }
  
   if (!(output = fopen (nameout, "w"))) {
      fprintf (stderr, "Cannot open output file! \n");
      exit (-1);
      }
 
   if (!(summary = fopen (namesum, "w"))) {
      fprintf (stderr, "Cannot open summary file! \n");
      exit (-1);
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

   iInd = (int *) malloc (numDs*sizeof(int)); 
   jInd = (int *) malloc (numDs*sizeof(int)); 

   lbd = (double *) malloc (numDs*sizeof(double)); 
   ubd = (double *) malloc (numDs*sizeof(double)); 

   for (k = 0; k < numDs; k++)
   fscanf (input,"%d %d %lf %lf \n",&iInd[k],&jInd[k],&lbd[k],&ubd[k]);
   
   if (!numAtoms) 
      for (k = 0; k < numDs; k ++) {
          if (numAtoms < iInd[k]) numAtoms = iInd[k];
          if (numAtoms < jInd[k]) numAtoms = jInd[k];
          }

   if (!numSols) numSols = 1;

   for (k = 0; k < numDs; k ++) {
       lbd [k] = lbd [k] * lbd [k];
       ubd [k] = ubd [k] * ubd [k];
       }

   /*
   Allocate memory.
   */
   n = 3*numAtoms;
   lintWkSp = 2*numAtoms > numDs ? 2*numAtoms : numDs;

   /*
   In general ldblWkSp is 2*(m+1)*n + 2*m + 1, where m is the number
   of vectors used by VMLM.
   */
 
   ldblWkSp = 22*n + 21;

   fOpt = (double *) malloc (sizeof(double));
   dErr = (double *) malloc (3*sizeof(double));
   xOpt = (double *) malloc (n*sizeof(double));

   intWkSp = (int *) malloc (lintWkSp*sizeof(int)); 
   dblWkSp = (double *) malloc (ldblWkSp*sizeof(double));

   /*
   Call the DG optimizer to obtain a set of solutions.
   */

   seed = 1;

   printf ("\nStarting optimization with DGSOL\n");

   for (k = 0; k < numSols; k++) {

     /*
     Call the DG optimizer.
     */

     dginitial(&numAtoms,xOpt,&numDs,iInd,jInd,lbd,ubd,&seed,intWkSp);
     dgopt(&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,intWkSp,dblWkSp);
     dgerr(&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,dErr);

     printf 
     ("\nOptimization %d - Largest distance error is %12.3e\n", k+1, dErr[2]);

     /*
     Collect solutions.
     */

     for (l = 0; l < numAtoms; l ++) 
         fprintf (output, 
         "%16.8e %16.8e %16.8e \n", xOpt[3*l],xOpt[3*l+1],xOpt[3*l+2]); 
     fprintf (output, "\n %20.12e \n\n", *fOpt);
 
     fprintf (summary,
     "%6d  %8d  %10.2e  %10.2e  %10.2e  %10.2e \n",
     numAtoms, numDs, *fOpt, dErr[0], dErr[1], dErr[2]);

     }

   printf ("\nEnd of optimization with DGSOL\n");

   /*
   Clean up.
   */

   free (iInd); 
   free (jInd);
 
   free (lbd);  
   free (ubd);

   free (fOpt);  
   free (xOpt);
   free (dErr);

   free (intWkSp);  
   free (dblWkSp);

   fclose (input);
   fclose (output);
   fclose (summary);

   }

