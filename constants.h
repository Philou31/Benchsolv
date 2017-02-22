//! \file constants.h
//! \brief Constants of the whole project
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Constants valid through whole software.
//!

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace cst
{
    const std::string ALPHANUM =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const std::string EMPTY_ID = "01234567";
    const std::string RUN_COUNT = "run.count";
 
    //MPI
    const int USE_COMM_WORLD = -987654;
    const int HOST_ID = 0;
    const int TAG_GO = 0;
    const int TAG_SIZE = 1;
    const int TAG_ARRAY = 2;
    const int TAG_IRN = 3;
    const int TAG_JCN = 4;
    
    //Defaults arguments
    const std::string EMPTY_STRING_OPT_KEY = "";
    const int EMPTY_INT_OPT_KEY = -42;
    const int EMPTY_INT_OPT_VALUE = -42;
    //Solver
    const std::string SOLVER = "solver";
    const char SOLVER_C = 's';
    const std::string MUMPS = "mumps";
    const std::string MUMPS_C = "m";
    const std::string QR_MUMPS = "qr_mumps";
    const std::string QR_MUMPS_C = "qrm";
    const int NIL_MEM_RELAX = 0;
    const int NIL_MEM_FACTOR = 1;
    // Phases
    const std::string ANALYSIS_PHASE = "analysis";
    const std::string FACTORIZATION_PHASE = "factorization";
    const std::string SOLVE_PHASE = "solve";
    const std::string METRICS_PHASE = "metrics";
    //Bool option
    const std::string TRUE = "true";
    const std::string FALSE = "false";
    //Multiple bench
    const std::string MULTIPLE_BENCH = "multiple_bench";
    const char MULTIPLE_BENCH_C = 'q';
    const std::string OPTION = "option";
    const std::string SINGLE = "single";
    const std::string MULTIPLE = "multiple";
    //Block Distribution Option
    const int DISTR_ROW_UNEVEN = -1;
    const int DISTR_ROW_BLOCK = 0;
    const int DISTR_COL_BLOCK = 1;
    const int DISTR_ARROW_BLOCK = 2;
    
    //FILES
    const std::string EMPTY_FILE = "";
    const char OPTS_DELIMITER = ' ';
    const int LOC_CHUNK_NUMBER = 128;
    //SHEMAT-Suite/mphase/small/
    const std::string DATA_FOLDER = "/home/filou/Bureau/CERFACS/Data/";
    const std::string A_WATER_MPHASE_SMALL_FILE = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/small/jac.mm";
    const std::string A_WATER_MPHASE_SMALL_FILE_LOC0 = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/small/jac.mm.loc_0";
    const std::string A_WATER_MPHASE_SMALL_FILE_LOC1 = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/small/jac.mm.loc_1";
    const std::string A_WATER_MPHASE_SMALL_FILE_LOC2 = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/small/jac.mm.loc_2";
    const std::string RHS_WATER_MPHASE_SMALL_FILE = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/small/F.mm";
    //SHEMAT-Suite/mphase/big/
    const std::string A_WATER_MPHASE_BIG_FILE = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/big/jac.mm";
    const std::string RHS_WATER_MPHASE_BIG_FILE = cst::DATA_FOLDER + "SHEMAT-Suite/mphase/big/F.mm";
    //Fares/
    const std::string A_FARES_FILE_PREFIX = cst::DATA_FOLDER + "Fares/PROCESS_00";
    const std::string RHS_FARES_FILE = cst::DATA_FOLDER + "Fares/RHS";
    //MHD
    //10
    const std::string A_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/A.txt";
    const std::string RHS_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/RHS.txt";
    const std::string Sol_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/Sol.txt";
    //50
    const std::string A_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/A.txt";
    const std::string RHS_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/RHS.txt";
    const std::string Sol_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/Sol.txt";
    
    //MATHS
    const double ABS_EPSILON = 1e-12;
    const double REL_EPSILON = 1e-8;
    const double TIME_RATIO = 1e-6;
}

#endif /* CONSTANTS_H */