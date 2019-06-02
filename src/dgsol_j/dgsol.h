//
// Created by giacomo on 01/06/19.
//

#ifndef DGSOL_DGSOL_H
#define DGSO
/**
 * This provides the summary for the current solution.
 */
struct dgsol_summary {
    double f_err;       // Value of the merit function
    double derr_min;    // Smallest error in the distances
    double derr_avg;    // Average error in the distances
    double derr_max;    // Largest error in the distances
};

#endif //DGSOL_DGSOL_H
