/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.0)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2010 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: fastlss.hpp,v 1.10 2010/12/21 09:53:47 cxsc Exp $ */

/*
**  FastPLSS: A library of (parallel) verified linear (interval) system 
**  solvers using C-XSC (V 0.2)
**
**  Author: Michael Zimmer
**
**  This software is based on:
**    - Module LinSys of the C-XSC-Toolbox
**      Authors: Rolf Hammer, Matthias Hocks, Dietmar Ratz
**    - Self-verifying solver for a dense system of linear equations
**      Authors: Carlos Holbig, Walter Kraemer, Paulo Sergio Morandi Junior,
**               Bernardo Frederes Kramer Alcalde, 
*/

 
#ifndef _CXSC_FASTLSS_HPP
#define _CXSC_FASTLSS_HPP

#include <rmatrix.hpp>
#include <imatrix.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp> 
//Computation of approximate inverse
#include <matinv_aprx.hpp>   


#if !defined(_CXSC_LSS_HPP) && !defined(_CXSC_ILSS_HPP) && !defined(_CXSC_CLSS_HPP) && !defined(_CXSC_CILSS_HPP)
#define _CXSC_LSS_HPP
#define _CXSC_ILSS_HPP
#define _CXSC_CLSS_HPP
#define _CXSC_CILSS_HPP
#define _CXSC_LSS_UNDEFINED
#endif

namespace cxsc {

#ifndef _CXSC_LSS_CONSTANTS_DEFINED
#define _CXSC_LSS_CONSTANTS_DEFINED
//Control constants
static const int
  LSS_ONLY_PART_ONE = 0,
  LSS_ONLY_PART_TWO = 1,
  LSS_BOTH_PARTS    = 2;
#endif

//! Translates the error codes of the solver into corresponding error messages
std::string LinSolveErrMsg(int);

#ifdef _CXSC_LSS_UNDEFINED
//! Entry function of real linear sytems solver
static inline void  lss(cxsc::rmatrix&, cxsc::rmatrix&, cxsc::imatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(cxsc::rmatrix&, cxsc::rvector&, cxsc::ivector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(cxsc::rmatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
#endif

#ifdef _CXSC_LSS_UNDEFINED
//! Entry function of interval linear sytems solver
static inline void  ilss(cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of interval linear sytems solver
static inline void  ilss(cxsc::imatrix&, cxsc::ivector&, cxsc::ivector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
#endif

#ifdef _CXSC_LSS_UNDEFINED
//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cmatrix&, cxsc::cimatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cvector&, cxsc::civector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
#endif

#ifdef _CXSC_LSS_UNDEFINED
//! Entry function of complex interval linear sytems solver
static inline void  cilss(cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex interval linear sytems solver
static inline void  cilss(cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);
#endif


#ifdef CXSC_USE_BLAS
static inline void IminusRA(const rmatrix& A, const rmatrix& B, imatrix& C) {
  int rnd = getround();

   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubC_j - lbC_j + 1;
   
   double *DA = new double[m*n];
   double *DB = new double[n*o];
   double *DC = new double[m*o];
   
   //Copy A and B into double array   
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(i-1)*n+(j-1)] = _double(-A[lbA_i+i-1][lbA_j+j-1]);
      }
   }        

   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DB[(i-1)*o+(j-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);         
      }
   }        


   int LDA=n,LDB=o,LDC=m;
   double alpha = 1.0, beta = 0.0;
    
   
   //Compute Infimum
   setround(-1);
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DA, 
          LDA, DB, LDB, beta, DC, LDC);
   
   //Copy infimum into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
           C[lbC_i+j-1][lbC_j+i-1] = DC[(i-1)*m+(j-1)];
      }
      C[lbC_i+i-1][lbC_j+i-1] = DC[(i-1)*m+(i-1)] + 1.0;
   }        
   
   //Compute supremum
   setround(1);
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DA, 
          LDA, DB, LDB, beta, DC, LDC);


   //copy supremum into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         SetSup(C[lbC_i+j-1][lbC_j+i-1], DC[(i-1)*m+(j-1)]);
      }
      SetSup(C[lbC_i+i-1][lbC_j+i-1], DC[(i-1)*m+(i-1)] + 1.0);
   }        
  
   delete[] DA;
   delete[] DB;
   delete[] DC; 

  setround(rnd);
}

static inline void IminusRA(const rmatrix& R, const imatrix& A, imatrix& C) {
  int rnd = getround();

  C = -R*A;

  for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++)
    C[i][i] += 1.0;

  setround(rnd);
}

static inline void IminusRA(const cmatrix& A, const cmatrix& B, cimatrix& C) {
   int rnd = getround();

   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubC_j - lbC_j + 1;
   
   double *DAR  = new double[m*n]; //A.re
   double *DAI  = new double[m*n]; //A.Im
   double *DAIm = new double[m*n]; //-A.Im
   double *DBR = new double[n*o];  //B.re
   double *DBI = new double[n*o];  //B.Im
   double *DCR = new double[m*o];  //C.re
   double *DCI = new double[m*o];  //C.Im
   
   //Copy A and B into corresponding double arrays
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((i-1)*n+(j-1));
         DAR[ind] = _double(-Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(-Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] =  -DAI[ind];         
      }
   }        

   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));         
         DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));                  
      }
   }        

   int LDA=n,LDB=o,LDC=m;
   double alpha = 1.0, beta = 0.0;
    
   
   //Compute lower bound of result    
   setround(-1);
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAR, 
          LDA, DBR, LDB, beta, DCR, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAIm, 
          LDA, DBI, LDB, beta, DCR, LDC);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAR, 
          LDA, DBI, LDB, beta, DCI, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAI, 
          LDA, DBR, LDB, beta, DCI, LDC);

   //Copy lower bound into C  
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         complex c = complex(DCR[ind],DCI[ind]);
         if(i==j) c += 1.0;
         C[lbC_i+j-1][lbC_j+i-1] = c;
      }
   }        
 
   //Compute upper bound of result
   setround(1);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAR, 
          LDA, DBR, LDB, beta, DCR, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAIm, 
          LDA, DBI, LDB, beta, DCR, LDC);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAR, 
          LDA, DBI, LDB, beta, DCI, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m, o, n, alpha, DAI, 
          LDA, DBR, LDB, beta, DCI, LDC);

   //Copy upper bound into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         complex c = complex(DCR[ind],DCI[ind]);
         if(i==j) c += 1.0;         
         SetSup(C[lbC_i+j-1][lbC_j+i-1], c);
      }
   }       

   setround(rnd);
   
   delete[] DAR;
   delete[] DAI;
   delete[] DAIm;
   delete[] DBR;
   delete[] DBI;   
   delete[] DCR;    
   delete[] DCI; 
}

static inline void IminusRA(const cmatrix& R, const cimatrix& A, cimatrix& C) {
  int rnd = getround();

  C = -R*A;

  for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++)
    C[i][i] += 1.0;

  setround(rnd);
}

#endif

} //namespace cxsc

#include <fastlss.inl>

namespace cxsc {

#ifdef _CXSC_LSS_HPP
//! Start function for solver
/*
  This function is called by the user to start the solver for a real system. 

    \param A Matrix of the system
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param K Precision to use for DotK-Algorithm (K-fold double precision).
             If K<2 the C-XSC accumulator will be used (=maximum precision, 
             but slower)
    \param np Number of thread to be used for OpenMP.
    \param msg If true, status messages will be put out to cout during 
               computation
    \param lsscfg Determines, if only part 1, only part 2 or both parts
                  of the solver are to be executed
*/
static inline void lss( rmatrix& A, rmatrix& b, imatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void lss( rmatrix& A, rvector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  rmatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  lss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

static inline void lss( rmatrix& A, imatrix& b, imatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<rmatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void lss( rmatrix& A, ivector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  lss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}
#endif

//----------------------------------------------------------------------------                      

#ifdef _CXSC_ILSS_HPP
//! Start function for solver
/*
  This function is called by the user to start the solver for a real interval system. 

    \param A Matrix of the system
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param K Precision to use for DotK-Algorithm (K-fold double precision).
             If K<2 the C-XSC accumulator will be used (=maximum precision, 
             but slower)
    \param np Number of thread to be used for OpenMP.
    \param msg If true, status messages will be put out to cout during 
               computation
    \param lsscfg Determines, if only part 1, only part 2 or both parts
                  of the solver are to be executed
*/
static inline void ilss( imatrix& A, imatrix& b, imatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void ilss( imatrix& A, ivector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  ilss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}
#endif

//----------------------------------------------------------------------------                      

#ifdef _CXSC_CLSS_HPP
//! Start function for solver
/*
  This function is called by the user to start the solver for a complex system. 

    \param A Matrix of the system
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param K Precision to use for DotK-Algorithm (K-fold double precision).
             If K<2 the C-XSC accumulator will be used (=maximum precision, 
             but slower)
    \param np Number of thread to be used for OpenMP.
    \param msg If true, status messages will be put out to cout during 
               computation
    \param lsscfg Determines, if only part 1, only part 2 or both parts
                  of the solver are to be executed
*/
static inline void clss( cmatrix& A, cmatrix& b, cimatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void clss( cmatrix& A, cvector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cmatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  clss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

static inline void clss( cmatrix& A, cimatrix& b, cimatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<cmatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void clss( cmatrix& A, civector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  clss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}
#endif

//----------------------------------------------------------------------------                      

#ifdef _CXSC_CILSS_HPP
//! Start function for solver
/*
  This function is called by the user to start the solver for a complex interval system. 

    \param A Matrix of the system
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param K Precision to use for DotK-Algorithm (K-fold double precision).
             If K<2 the C-XSC accumulator will be used (=maximum precision, 
             but slower)
    \param np Number of thread to be used for OpenMP.
    \param msg If true, status messages will be put out to cout during 
               computation
    \param lsscfg Determines, if only part 1, only part 2 or both parts
                  of the solver are to be executed
*/
static inline void cilss( cimatrix& A, cimatrix& b, cimatrix& xx, int& Err, int K, bool msg, int lsscfg, int np) {
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, b, xx, Err, K, np, msg, lsscfg);
}

static inline void cilss( cimatrix& A, civector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  X[Col(1)] = x;
  cilss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}
#endif

} //namespace cxsc

#endif 


