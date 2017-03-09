//! \file Parameters_ABCD.h
//! \brief Constants of the ABCD class
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Constants used in the ABCDs class.
//!

#ifndef PARAMETERS_ABCD_H
#define PARAMETERS_ABCD_H

#include "defaults.h"

namespace parab {
    // JOB
    const int JOB_INIT(-1), JOB_END(-2), JOB_ANALYSIS(1), JOB_FACTO(2),
    JOB_ANAL_FACTO(4), JOB_SOLVE(3), JOB_FACTO_SOL(5), JOB_ALL(6);
    // SYM
    const int SYM_UNSYM(false), SYM_SYM(true);

    // Number of partitions
    const int NBPARTS(Controls::nbparts);
    const int NBPARTS_DEFAULT(8);

    // Partitioning strategy
    const int PART_TYPE(Controls::part_type);
    const int PART_TYPE_MANUAL(1);
    const int PART_TYPE_AUTO_UNIF(2);
    const int PART_TYPE_AUTO_HYPER(3);

    // Guess the number of partitions
    const int PART_GUESS(Controls::part_guess);
    const int PART_GUESS_NO(0);
    const int PART_GUESS_YES(1);

    // Scaling
    const int SCALING(Controls::scaling);
    const int SCALE_NO(0);
    const int SCALE_YES(1);

    // Max number of iterations
    const int ITMAX(Controls::itmax);
    const int ITMAX_DEFAULT(1000);

    // Block-CG block size
    const int BLOCK_SIZE(Controls::block_size);
    const int BLOCK_SIZE_DEFAULT(1);    // If >1, Stabilized Block-CG

    // Verbose level
    const int VERBOSE(Controls::verbose_level);
    const int VERBOSE_DEFAULT(6);
    const int VERBOSE_NO(0);

    // Augmentation type
    const int AUG(Controls::aug_type);
    const int AUG_NO(0);
    const int AUG_Cij(1);   // Should be used with scaling
    const int AUG_Aij(0);

    // Blocking factor
    const int AUG_BLOCK(Controls::aug_blocking);
    const int AUG_BLOCK_DEFAULT(128);
}

#endif /* PARAMETERS_ABCD_H */

