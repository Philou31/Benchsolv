//! \file Parameters_QR_Mumps.h
//! \brief Constants of the QR_Mumps class
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Constants used in the QR_Mumps class.
//!

#ifndef PARAMETER_QR_MUMPS_H
#define PARAMETER_QR_MUMPS_H

namespace parqrm {
    // Output
    const std::string OUNIT = "qrm_ounit";
    const std::string EUNIT = "qrm_eunit";
    const std::string DUNIT = "qrm_dunit";

    // Ordering
    const std::string ORDERING = "qrm_ordering";
//    ords ORD_AUTO = qrm_auto;
//    ords ORD_NATURAL = qrm_natural_;
//    ords ORD_GIVEN = qrm_given_;
//    ords ORD_CAMD = qrm_colamd_;
//    ords ORD_SCOTCH = qrm_scotch_;
//    ords ORD_METIS = qrm_metis_;
    
    // Keeph
    const std::string KEEPH = "qrm_keeph";
    const std::string MB = "qrm_mb";
    const std::string NB = "qrm_nb";
    const std::string IB = "qrm_ib";
    const std::string BH = "qrm_bh";
    const std::string RHSNB = "qrm_rhsnb";
    const std::string MEM_RELAX = "qrm_mem_relax";
    
    //Constants
//    const std::string YES = "qrm_yes_";
//    const std::string NO = "qrm_no_";
    
    //Infos
    //Global
    const std::string MAX_MEM = "qrm_max_mem";
    const std::string TOT_MEM = "qrm_tot_mem";
    
    //Local
//    const std::string E_R_NZ = "qrm_e_nnz_r";
//    const std::string R_NZ = "qrm_nnz_r";
//    const std::string E_H_NZ = "qrm_e_nnz_h";
//    const std::string H_NZ = "qrm_e_nnz_h";
//    const std::string E_FLOPS = "qrm_e_facto_flops";
//    const std::string E_MEM_PEAK = "qrm_e_facto_mem_peak";
}

#endif /* PARAMETER_QR_MUMPS_H */

