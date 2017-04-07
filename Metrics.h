//! \file Metrics.h
//! \brief Metrics class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class Metrics computing metrics on the problem and the found solution. This
//! class is meant to be fair in the results gotten from the solvers: metrics 
//! are computed in the exact same way.
//!

#ifndef METRICS_H
#define METRICS_H

#include "QR_Mumps.h"

class Metrics {
private:
    dqrm_spmat_type_c _qrm_mat; // QR_Mumps data structure
    double *_r, *_x;    // residual (Ax-b) and solution arrays
public:
    //!
    //! \brief Constructor of the Metrics class
    //!
    //! \return the Metrics instance
    //!
    Metrics();
    
    //!
    //! \brief Destructor of the Metrics class
    //!
    ~Metrics();
    
    //!
    //! \fn void init_metrics(int m, int n, int nz, double *values, int *irn,
    //! int*jcn, double *x, double *r)
    //! \brief Initialize the data (QR_Mumps structure, residual and solution)
    //!
    //! \param m, n, nz: number of rows, columns and non-zero values
    //! \param values, irn/jcn: arrays of non-zero values, row/column indices
    //! \param x: solution array
    //! \param r: residual array (Ax-b)
    //!
    void init_metrics(int m, int n, int nz, double *values, int *irn,
        int*jcn, double *x, double *r);
    
    //!
    //! \fn A_norm(double &anrm)
    //! \brief Compute the norm of the matrix A
    //!
    //! \param anrm: norm of A
    //!
    void A_norm(double &anrm);
    
    //!
    //! \fn b_norm(double &bnrm)
    //! \brief Compute the norm of the right hand side b
    //!
    //! \param bnrm: norm of the right hand side b
    //!
    void b_norm(double &bnrm);
    
    //!
    //! \fn sol_norm(double &xnrm)
    //! \brief Compute the norm of the solution x
    //!
    //! \param xnrm: norm of the solution x
    //!
    void sol_norm(double &xnrm);
    
    //!
    //! \fn r_norm(double &rnrm)
    //! \brief Compute the norm of the residual r=Ax-b
    //!
    //! \param rnrm: norm of the residual
    //!
    void residual_norm(double &rnrm);
    
    //!
    //! \fn o_norm(double &onrm)
    //! \brief Compute the norm of the orthogonal residual onrm=||A'r||/||r||
    //!
    //! \param onrm: norm of the orthogonal residual
    //!
    void residual_orth(double &onrm);
    
    //!
    //! \fn compute_metrics(double rnrm, double onrm, double anrm, double xnrm, 
    //!     double bnrm)
    //! \brief Compute the different metrics
    //!
    //! This method will compute all the different metrics:
    //!     - norm of A
    //!     - norm of b
    //!     - norm of the solution
    //!     - norm of the residual r=Ax-b
    //!     - norm of the orthogonal residual: ||A'r||/||r||
    //!
    //! \param onrm: norm of the orthogonal residual
    //!
    void compute_metrics(double &rnrm, double &onrm, double &anrm, double &xnrm, 
        double &bnrm);
};

#endif /* METRICS_H */

