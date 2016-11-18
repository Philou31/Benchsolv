//! \file constants.h
//! \brief Constants for the project
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Constants for the project:
//!     - file names
//!     - physical constants
//!     - parameter values
//!     - mathematical constants
//!
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace cst
{
    //MPI
    const int USE_COMM_WORLD = -987654;
    const int HOST_ID = 0;
    
    //FILES
    const std::string EMPTY_FILE = "";
    const char OPTS_DELIMITER = ' ';
    //data.eocoe.eu/water/SHEMAT-Suite/mphase/small/
    const std::string DATA_FOLDER = "/home/filou/Bureau/CERFACS/Data/";
    const std::string A_WATER_MPHASE_SMALL_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mphase/small/jac.mm";
    const std::string RHS_WATER_MPHASE_SMALL_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mphase/small/F.mm";
    //data.eocoe.eu/water/SHEMAT-Suite/mphase/big/
    const std::string A_WATER_MPHASE_BIG_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mphase/big/jac.mm";
    const std::string RHS_WATER_MPHASE_BIG_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mphase/big/F.mm";
    //data.eocoe.eu/water/SHEMAT-Suite/mpmc/small/
    const std::string A_WATER_MPMC_SMALL_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mpmc/small/jac.mm";
    const std::string RHS_WATER_MPMC_SMALL_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mpmc/small/F.mm";
    //data.eocoe.eu/water/SHEMAT-Suite/mpmc/big/
    const std::string A_WATER_MPMC_BIG_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mpmc/big/jac.mm";
    const std::string RHS_WATER_MPMC_BIG_FILE = cst::DATA_FOLDER + "data.eocoe.eu/water/SHEMAT-Suite/mpmc/big/F.mm";
    //Fares/
    const std::string A_FARES_FILE_PREFIX = cst::DATA_FOLDER + "Fares/PROCESS_00";
    const std::string RHS_FARES_FILE = cst::DATA_FOLDER + "Fares/RHS";
    //MHD
    //10
    const std::string A_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/A.txt";
    const std::string RHS_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/RHS.txt";
    const std::string Sol_MHD_10_FILE = cst::DATA_FOLDER + "MHD/10/Sol.txt";
    //100
    const std::string A_MHD_100_FILE = cst::DATA_FOLDER + "MHD/100/A.txt";
    const std::string RHS_MHD_100_FILE = cst::DATA_FOLDER + "MHD/100/RHS.txt";
    const std::string Sol_MHD_100_FILE = cst::DATA_FOLDER + "MHD/100/Sol.txt";
    //150
    const std::string A_MHD_150_FILE = cst::DATA_FOLDER + "MHD/150/A.txt";
    const std::string RHS_MHD_150_FILE = cst::DATA_FOLDER + "MHD/150/RHS.txt";
    const std::string Sol_MHD_150_FILE = cst::DATA_FOLDER + "MHD/150/Sol.txt";
    //50
    const std::string A_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/A.txt";
    const std::string RHS_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/RHS.txt";
    const std::string Sol_MHD_50_FILE = cst::DATA_FOLDER + "MHD/50/Sol.txt";
    //Commented Lines in MatrixMarket format
    const int DEFAULT_NB_COMMENTED_LINES = 0;
    
    //MATHS
    const double ABS_EPSILON = 1e-12;
    const double REL_EPSILON = 1e-8;
}

#endif /* CONSTANTS_H */

