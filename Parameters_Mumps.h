/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Parameters_Mumps.h
 * Author: filou
 *
 * Created on 13 novembre 2016, 18:34
 */

#ifndef PARAMETERS_MUMPS_H
#define PARAMETERS_MUMPS_H

namespace parm {
    // JOB
    const int JOB_INIT(-1), JOB_END(-2), JOB_ANALYSIS(1), JOB_FACTO(2),
    JOB_ANAL_FACTO(4), JOB_SOLVE(3), JOB_FACTO_SOL(5), JOB_ALL(6);
    // SYM
    const int SYM_UNSYM(0), SYM_GENERAL(1), SYM_DEFPOS(2);
    // PAR
    const int WORKING_HOST(1), NOT_WORKING_HOST(0);
    const int HOST_ID(0);
    
    //Output stream
    const int OUT_ERROR(1);
    const int OUT_DIAGNOSTIC(2);
    const int OUT_GINFO(3);
    const int OUT_LEVEL(4);

    //Input Matrix
    //Format
    const int A_FORMAT(5);
    const int A_ASSEMBLED_FORMAT(0);
    const int A_ELEMENTAL_FORMAT(1);
    //Distribution
    const int A_DISTRIBUTION(18);
    const int A_CENTRALIZED(0);
    const int A_DISTR_FACTO_MAPPING(1);
    const int A_DISTR_FACTO(2);
    const int A_DISTR_ANALYSIS(3);
    
    //RHS Matrix
    //Format
    const int RHS_FORMAT(20);
    const int RHS_DENSE_FORMAT(0);
    const int RHS_YES_SPARSITY(1);
    const int RHS_NO_SPARSITY(2);
    const int RHS_AUTO_SPARSITY(3);

    //Permutation
    const int PERMUTATION(6);
    const int PERM_NO(0);
    const int PERM_MAX_CARD_DIAG(1);
    const int PERM_MAX_MIN_DIAG(2);
    const int PERM_VAR_MAX_MIN_DIAG(3); //Needs A at Analysis
    const int PERM_MAX_SUM_DIAG(4);
    const int PERM_MAX_MULT_DIAG(5);
    const int PERM_VAR_MAX_MULT_DIAG(6); //Needs A at Analysis
    const int PERM_AUTO_PERM(7);

    //Scaling
    const int SCALING(8);
    const int SCALE_ANALYSIS(-2);
    const int SCALE_USER(-1); //Provided in COLSCA and ROWSCA
    const int SCALE_NO(0);
    const int SCALE_DIAG(1); //during facto
    const int SCALE_COL(3); //during facto
    const int SCALE_ROW_COL(4); //based on infinite row/col norm, during facto
    const int SCALE_SIM_ROW_COL(7); //during facto
    const int SCALE_VAR_SIM_ROW_COL(8); //during facto, more rigorous & expensive
    const int SCALE_AUTO(77);

    //Symmetric permutation
    //is analysis sequential or parallel
    const int SEQPAR_ANALYSIS(28);
    const int ANAL_AUTO_SEQPAR(0);
    const int ANAL_SEQ(1);
    const int ANAL_PAR(2);
    //sequential permutation (ICNTL(SEQPAR_ANALYSIS)=PAR_ANAL)
    const int SYMPERM_SEQ(7);
    const int SYMPERM_AMD(0);
    const int SYMPERM_USER(1);
    const int SYMPERM_AMF(2);
    const int SYMPERM_SCOTCH(3);
    const int SYMPERM_PORD(4);
    const int SYMPERM_METIS(5);
    const int SYMPERM_QAMD(6);
    const int SYMPERM_SEQAUTO(7);
    //Parallel permutation (ICNTL(SEQPAR_ANALYSIS)=SEQ_ANAL)
    const int SYMPERM_PAR(29);
    const int SYMPERM_PARAUTO(0);
    const int SYMPERM_PTSCOTCH(1);
    const int SYMPERM_PARMETIS(2);

    //Iterative refinement
    const int ITER_REFINEMENT = 10;
    const int ITERREF_NUM = 10;
    const int ITERREF_NO = 0;
    const int ITERREF_MAX_NUM = -3;
    const int ITERREF_CRITER_MAX_NUM = 10;
    const int ITERREF_CNTL = 2; // CNTL(2) is alpha

    //Error Analysis
    const int ERR_ANALYSIS = 11;
    const int ERRANAL_NO = 0;
    //om1=RINFOG(7), om2=RINFOG(8), ||x||inf=RINFOG(5), scaled res=RINFOG(6)
    const int ERRANAL_PART = 2;
    //+ cond1=RINFOG(10), cond2=RINFOG(11), forward error upper bond=RINFOG(9)
    const int ERRANAL_FULL = 1;

    //Out of Core
    const int INOUT_CORE = 22;
    const int IN_CORE = 0;
    const int OUT_CORE = 1;

    //Null pivot row detection
    const int NULL_PIVOT = 24;
    const int NULL_PIVOT_NO = 0;
    const int NULL_PIVOT_YES = 1;

    //Discard matrix factor
    const int DISCARD = 31;
    const int DISC_NO = 0;
    const int DISC_NO_SOLVE = 1;
    const int DISC_SOLVE = 2;
    
    //Determinant
    const int DETERMINANT = 33;
    const int DET_NO = 0;
    const int DET_YES = 1;
    
    //Forward elimination during facto
    const int FACTO_ELIM = 32;
    const int FELIM_NO = 0;
    const int FELIM_YES = 1;
    
    //Entries of A-1
    const int INV_ENTRIES = 30;
    const int INV_NO = 0;
    const int INV_YES = 1;
    
    //SCHUR TODO
    //WORKSPACE
    const int MEMORY_PERCENT_INC = 14;
    const int MEMORY_DEFAULT_PERCENT_INC = 80;
    const int MEMORY_SIZE = 23;
    const int MEMORY_DEFAULT_SIZE = 1000;
}

#endif /* PARAMETERS_MUMPS_H */

