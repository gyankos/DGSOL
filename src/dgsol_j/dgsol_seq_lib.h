//
// Created by giacomo on 01/06/19.
//

#include "dgsol.h"

#ifndef DGSOL_DGSOL_SEQ_LIB_H
#define DGSOL_DGSOL_SEQ_LIB_H

/**
 *
 * @param numDs        Number of all the known distances
 * @param numAtoms        Number of all the involved nodes
 * @param iInd              Array of all the source nodes
 * @param jInd              Array of all the target nodes
 * @param lbd               Array of all the lower bounds
 * @param ubd               Array of all the upper bounds
 * @param mpi               Choose whether the computation has to be parallelized using openmpi or not.
 *
 * @return            The three dimensional solution, where a_1,a_2,a_3,b_1,b_2,b_3,....
 */
double* dgsol_j(int numDs, int numAtoms, int*iInd, int* jInd, double* lbd, double* ubd,int* seed, struct dgsol_summary* summary);

#endif //DGSOL_DGSOL_SEQ_LIB_H
