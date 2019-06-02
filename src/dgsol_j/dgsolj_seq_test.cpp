//
// Created by giacomo on 01/06/19.
//

extern "C" {
#include "dgsol_seq_lib.h"
#include <stdio.h>
}

int main(int argc, char* argv[]) {
    int src[3] = {1,2,3};
    int dst[3] = {2,3,1};
    double min[3] = {1.0, 99999, 2.0};
    double max[3] = {1.0, 99999, 2.0};
    int seed;
    struct dgsol_summary summary1;
    dgsol_j(3, 3, (int*)src, (int*)dst, (double*)min, (double*)max, &seed, &summary1);
    double* result =  dgsol_j(3, 3, src, dst, min, min, &seed, &summary1);
    for (int i = 0; i<3; i++) {
        printf(
                "%16.8e %16.8e %16.8e \n",result[3*i], result[3*i+1],result[3*i+2]);
    }
}