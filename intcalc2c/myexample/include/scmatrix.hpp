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

/* CVS $Id: scmatrix.hpp,v 1.12 2010/12/21 09:54:01 cxsc Exp $ */

#ifndef _CXSC_SCMATRIX_HPP_INCLUDED
#define _CXSC_SCMATRIX_HPP_INCLUDED

#include <complex.hpp>
#include <cmatrix.hpp>
#include <scvector.hpp>
#include <cidot.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sparsecdot.hpp>
#include <sparsematrix.hpp>
#include <srmatrix.hpp>

namespace cxsc {

//definiert in srmatrix.hpp
//enum STORAGE_TYPE{triplet,compressed_row,compressed_column};

class scmatrix_slice;
class scmatrix_subv;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;


inline bool comp_pair_c(std::pair<int,complex> p1, std::pair<int,complex> p2) {
  return p1.first < p2.first;
}

class scmatrix {

  private:
    std::vector<int> p;
    std::vector<int> ind;
    std::vector<complex> x;
    int m;
    int n;
    int lb1,ub1,lb2,ub2;

  public:

    std::vector<int>& column_pointers() {
      return p;
    }

    std::vector<int>& row_indices() {
      return ind;
    }

    std::vector<complex>& values() {
      return x;
    }

    const std::vector<int>& column_pointers() const {
      return p;
    }

    const std::vector<int>& row_indices() const {
      return ind;
    }

    const std::vector<complex>& values() const {
      return x;
    }

    scmatrix() {
      p.push_back(0);
      m = n = 0;
      lb1 = lb2 = ub1 = ub2 = 0;
    }

    scmatrix(const int r, const int c) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(2*(m+n));
      x.reserve(2*(m+n));

      p[0] = 0;
    }

    scmatrix(const int r, const int c, const int e) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(e);
      x.reserve(e);

      p[0] = 0;
    }

    scmatrix(const int m, const int n, const int nnz, const intvector& rows, const intvector& cols, const cvector& values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<complex> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<complex>(rows[Lb(rows)+k],cols[Lb(cols)+k],values[Lb(values)+k]));
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
         
      } else if(t == compressed_row) {

         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[Lb(rows)+i];

         std::vector<triplet_store<complex> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<complex>(j,cols[Lb(cols)+k],values[Lb(values)+k]));
           }
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
    
      } else if(t == compressed_column) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[Lb(rows)+i];

         std::vector<std::pair<int,complex> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[Lb(cols)+k],values[Lb(values)+k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_c);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    scmatrix(const int m, const int n, const int nnz, const int* rows, const int* cols, const complex* values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<complex> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<complex>(rows[k],cols[k],values[k]));
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
         
      } else if(t == compressed_row) {

         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[i];

         std::vector<triplet_store<complex> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<complex>(j,cols[k],values[k]));
           }
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
    
      } else if(t == compressed_column) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[i];

         std::vector<std::pair<int,complex> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[k],values[k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_c);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }


    scmatrix(const srmatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(complex(A.x[i]));
    }


    scmatrix(const rmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(complex(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scmatrix(const cmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(complex(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scmatrix(const srmatrix_slice&);
    scmatrix(const scmatrix_slice&);

    void full(cmatrix& A) const {
       A = cmatrix(lb1,ub1,lb2,ub2);
       A = 0.0;
       for(int j=0 ; j<n ; j++) {
          for(int k=p[j] ; k<p[j+1] ; k++) {
             A[ind[k]+lb1][j+lb2] = x[k];
          }
       }
    }

    void dropzeros() {
      std::vector<int> pnew(n+1,0);
      std::vector<int> indnew;
      std::vector<complex> xnew;
      int nnznew = 0;

      for(int j=0 ; j<n ; j++) {
        for(int k=p[j] ; k<p[j+1] ; k++) {
          if(x[k] != 0.0) {
            xnew.push_back(x[k]);
            indnew.push_back(ind[k]);
            nnznew++;
          }
        }
        pnew[j+1] = nnznew;
      }

      p = pnew;
      ind = indnew;
      x = xnew;
    }


    scmatrix& operator=(const real& A) {
      return sp_ms_assign<scmatrix,real,complex>(*this,A);
    }

    scmatrix& operator=(const complex& A) {
      return sp_ms_assign<scmatrix,complex,complex>(*this,A);
    }

    scmatrix& operator=(const rmatrix& A) {
      return spf_mm_assign<scmatrix,rmatrix,complex>(*this,A);
    }

    scmatrix& operator=(const cmatrix& A) {
      return spf_mm_assign<scmatrix,cmatrix,complex>(*this,A);
    }

    scmatrix& operator=(const rmatrix_slice& A) {
      return spf_mm_assign<scmatrix,rmatrix_slice,complex>(*this,A);
    }

    scmatrix& operator=(const cmatrix_slice& A) {
      return spf_mm_assign<scmatrix,cmatrix_slice,complex>(*this,A);
    }

    scmatrix& operator=(const srmatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(complex(A.x[i]));
      return *this;
    }

    /* scmatrix& operator=(const scmatrix& A) {
      p = A.p;
      ind = A.ind;
      x = A.x;
      return *this;
    } */

    scmatrix& operator=(const srmatrix_slice&);
    scmatrix& operator=(const scmatrix_slice&);

    const complex operator()(int i, int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator()(int, int)"));
#endif
      complex r(0.0);
      for(int k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  r = x[k];
      }
      return r;
    }

    complex& element(int i, int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::element()(int, int)"));
#endif
      int k;
      for(k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  return x[k];
      }

      //Nicht gefunden, Element muss angelegt werden, da Schreibzugriff moeglich
      std::vector<int>::iterator ind_it = ind.begin() + k;
      std::vector<complex>::iterator x_it  = x.begin() + k;
      ind.insert(ind_it, i-lb1);
      x_it = x.insert(x_it, complex(0.0));
      for(k=j-lb2+1 ; k<(int)p.size() ; k++)
        p[k]++;

      return *x_it;
    }

    scmatrix_subv operator[](const cxscmatrix_column&);
    scmatrix_subv operator[](const int);
    const scmatrix_subv operator[](const cxscmatrix_column&) const;
    const scmatrix_subv operator[](const int) const;

    scmatrix_slice operator()(const int, const int , const int, const int);

    scmatrix operator()(const intvector& pervec, const intvector& q) {
      scmatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      int nnz=0;
      for(int k=0 ; k<n ; k++) {
        A.p[k] = nnz;

        std::map<int,complex> work;
        for(int j=p[q[Lb(q)+k]] ; j<p[q[Lb(q)+k]+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,complex>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }

        nnz += work.size();
 
      }

      A.p[n] = nnz;

      return A;
    }

    scmatrix operator()(const intvector& pervec) {
      scmatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      for(int k=0 ; k<n ; k++) {
        A.p[k] = p[k];

        std::map<int,complex> work;
        for(int j=p[k] ; j<p[k+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,complex>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }
 
      }

      A.p[n] = p[n];

      return A;
    }

    scmatrix operator()(const intmatrix& P, const intmatrix& Q) {
      intvector p = permvec(P);
      intvector q = perminv(permvec(Q));
      return (*this)(p,q);
    }

    scmatrix operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    real density() const {
      return p[n]/((double)m*n);
    }

    int get_nnz() const {
      return p[n];
    }

    scmatrix& operator+=(const rmatrix& B) {
      return spf_mm_addassign<scmatrix,rmatrix,cmatrix>(*this,B);
    }

    scmatrix& operator+=(const cmatrix& B) {
      return spf_mm_addassign<scmatrix,cmatrix,cmatrix>(*this,B);
    }

    scmatrix& operator+=(const rmatrix_slice& B) {
      return spf_mm_addassign<scmatrix,rmatrix_slice,cmatrix>(*this,B);
    }

    scmatrix& operator+=(const cmatrix_slice& B) {
      return spf_mm_addassign<scmatrix,cmatrix_slice,cmatrix>(*this,B);
    }

    scmatrix& operator+=(const srmatrix& B) {
      return spsp_mm_addassign<scmatrix,srmatrix,complex>(*this,B);
    }

    scmatrix& operator+=(const scmatrix& B) {
      return spsp_mm_addassign<scmatrix,scmatrix,complex>(*this,B);
    }

    scmatrix& operator-=(const rmatrix& B) {
      return spf_mm_subassign<scmatrix,rmatrix,cmatrix>(*this,B);
    }

    scmatrix& operator-=(const cmatrix& B) {
      return spf_mm_subassign<scmatrix,cmatrix,cmatrix>(*this,B);
    }

    scmatrix& operator-=(const rmatrix_slice& B) {
      return spf_mm_subassign<scmatrix,rmatrix_slice,cmatrix>(*this,B);
    }

    scmatrix& operator-=(const cmatrix_slice& B) {
      return spf_mm_subassign<scmatrix,cmatrix_slice,cmatrix>(*this,B);
    }

    scmatrix& operator-=(const srmatrix& B) {
      return spsp_mm_subassign<scmatrix,srmatrix,complex>(*this,B);
    }

    scmatrix& operator-=(const scmatrix& B) {
      return spsp_mm_subassign<scmatrix,scmatrix,complex>(*this,B);
    }

    scmatrix& operator*=(const cmatrix& B) {
      return spf_mm_multassign<scmatrix,cmatrix,sparse_cdot,cmatrix>(*this,B);
    }

    scmatrix& operator*=(const rmatrix& B) {
      return spf_mm_multassign<scmatrix,rmatrix,sparse_cdot,cmatrix>(*this,B);
    }

    scmatrix& operator*=(const rmatrix_slice& B) {
      return spf_mm_multassign<scmatrix,rmatrix_slice,sparse_cdot,cmatrix>(*this,B);
    }

    scmatrix& operator*=(const cmatrix_slice& B) {
      return spf_mm_multassign<scmatrix,cmatrix_slice,sparse_cdot,cmatrix>(*this,B);
    }

    scmatrix& operator*=(const srmatrix& B) {
      return spsp_mm_multassign<scmatrix,srmatrix,sparse_cdot,complex>(*this,B);
    }

    scmatrix& operator*=(const scmatrix& B) {
      return spsp_mm_multassign<scmatrix,scmatrix,sparse_cdot,complex>(*this,B);
    }

    scmatrix& operator*=(const real& r) {
      return sp_ms_multassign(*this,r);
    }

    scmatrix& operator*=(const complex& r) {
      return sp_ms_multassign(*this,r);
    }

    scmatrix& operator/=(const real& r) {
      return sp_ms_divassign(*this,r);
    }

    scmatrix& operator/=(const complex& r) {
      return sp_ms_divassign(*this,r);
    }

    friend void SetLb(scmatrix&, const int, const int);
    friend void SetUb(scmatrix&, const int, const int);    
    friend int Lb(const scmatrix&, int);
    friend int Ub(const scmatrix&, int);
    friend int RowLen(const scmatrix&);
    friend int ColLen(const scmatrix&);
    friend srmatrix Re(const scmatrix&);
    friend srmatrix Im(const scmatrix&);
    friend scmatrix Inf(const scimatrix&);
    friend scmatrix Sup(const scimatrix&);
    friend scmatrix mid(const scimatrix&);
    friend scmatrix diam(const scimatrix&);

    friend scmatrix transp(const scmatrix&);
    friend scmatrix Id(const scmatrix&);

    friend std::istream& operator>>(std::istream&, scmatrix_slice&);
    friend std::istream& operator>>(std::istream&, scmatrix_subv&);

    friend class srmatrix_slice;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix_slice;
    friend class scmatrix_subv;
    friend class scvector;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class cmatrix;
    friend class cimatrix;


#include "matrix_friend_declarations.inl"
};

inline cmatrix::cmatrix(const srmatrix& A) {
  dat = new complex[A.m*A.n];
  lb1 = A.lb1; lb2 = A.lb2; ub1 = A.ub1; ub2 = A.ub2;
  xsize = A.n;
  ysize = A.m;
  *this = 0.0;
  for(int j=0 ; j<A.n ; j++) {
     for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        dat[A.ind[k]*A.n+j] = A.x[k];
     }
  }
}

inline cmatrix::cmatrix(const scmatrix& A) {
  dat = new complex[A.m*A.n];
  lb1 = A.lb1; lb2 = A.lb2; ub1 = A.ub1; ub2 = A.ub2;
  xsize = A.n;
  ysize = A.m;
  *this = 0.0;
  for(int j=0 ; j<A.n ; j++) {
     for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        dat[A.ind[k]*A.n+j] = A.x[k];
     }
  }
}

inline scmatrix Id(const scmatrix& A) {
  scmatrix I(A.m, A.n, (A.m>A.n) ? A.m : A.n);
  I.lb1 = A.lb1; I.lb2 = A.lb2;
  I.ub1 = A.ub1; I.ub2 = A.ub2;

  if(A.m < A.n) {
    for(int i=0 ; i<A.m ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(complex(1.0));
    }
  } else {
    for(int i=0 ; i<A.n ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(complex(1.0));
    }
  }

  return I;
}

inline scmatrix transp(const scmatrix& A) {
  scmatrix B(A.n, A.m, A.get_nnz());
    
  //Nichtnullen pro Zeile bestimmen
  std::vector<int> w(A.m,0);
  for(unsigned int i=0 ; i<A.ind.size() ; i++) 
    w[A.ind[i]]++;

  //Spalten"pointer" setzen
  B.p.resize(A.m+1);
  B.p[0] = 0;
  for(unsigned int i=1 ; i<B.p.size() ; i++)
    B.p[i] = w[i-1] + B.p[i-1];

  //w vorbereiten
  w.insert(w.begin(), 0); 
  for(unsigned int i=1 ; i<w.size() ; i++) {
    w[i] += w[i-1];
  }

  //neuer zeilenindex und wert wird gesetzt
  int q;
  B.ind.resize(A.get_nnz());
  B.x.resize(A.get_nnz());
  for(int j=0 ; j<A.n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      q = w[A.ind[k]]++;
      B.ind[q] = j;
      B.x[q] = A.x[k];
    }
  }

  return B;
}

inline void SetLb(scmatrix& A, const int i, const int j) {
  if(i==1) {
    A.lb1 = j;
    A.ub1 = j + A.m - 1;
  } else if(i==2) {
    A.lb2 = j;
    A.ub2 = j + A.n - 1;
  }
}

inline void SetUb(scmatrix& A, const int i, const int j) {
  if(i==1) {
    A.ub1 = j;
    A.lb1 = j - A.m + 1;
  } else if(i==2) {
    A.ub2 = j;
    A.lb2 = j - A.n + 1;
  }
}


inline int Lb(const scmatrix& A, int i) {
  if(i==1) 
    return A.lb1;
  else if(i==2)
    return A.lb2;
  else
    return 1;
}

inline int Ub(const scmatrix& A, int i) {
  if(i==1) 
    return A.ub1;
  else if(i==2)
    return A.ub2;
  else
    return 1;
}

inline int RowLen(const scmatrix& A) {
  return A.n;
}

inline int ColLen(const scmatrix& A) {
  return A.m;
}

inline void Resize(scmatrix& A) {
  sp_m_resize(A);
}

inline void Resize(scmatrix& A, const int m, const int n) {
  sp_m_resize(A,m,n);
}

inline void Resize(scmatrix& A, const int l1, const int u1, const int l2, const int u2) {
  sp_m_resize(A,l1,u1,l2,u2);
}

inline srmatrix Re(const scmatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Re(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix Im(const scmatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Im(A.x[i]));

  res.dropzeros();

  return res; 
}


inline cmatrix operator*(const cmatrix& A, const srmatrix& B) {
  return fsp_mm_mult<cmatrix,srmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const rmatrix& A, const scmatrix& B) {
  return fsp_mm_mult<rmatrix,scmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const cmatrix& A, const scmatrix& B) {
  return fsp_mm_mult<cmatrix,scmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const scmatrix& A, const rmatrix& B) {
  return spf_mm_mult<scmatrix,rmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const srmatrix& A, const cmatrix& B) {
  return spf_mm_mult<srmatrix,cmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const scmatrix& A, const cmatrix& B) {
  return spf_mm_mult<scmatrix,cmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const cmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_mult<cmatrix_slice,srmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const rmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_mult<rmatrix_slice,scmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const cmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_mult<cmatrix_slice,scmatrix,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const scmatrix& A, const rmatrix_slice& B) {
  return spf_mm_mult<scmatrix,rmatrix_slice,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const srmatrix& A, const cmatrix_slice& B) {
  return spf_mm_mult<srmatrix,cmatrix_slice,cmatrix,sparse_cdot>(A,B);
}

inline cmatrix operator*(const scmatrix& A, const cmatrix_slice& B) {
  return spf_mm_mult<scmatrix,cmatrix_slice,cmatrix,sparse_cdot>(A,B);
}

inline scmatrix operator*(const scmatrix& A, const srmatrix& B) {
  return spsp_mm_mult<scmatrix,srmatrix,scmatrix,sparse_cdot,complex>(A,B);
}

inline scmatrix operator*(const srmatrix& A, const scmatrix& B) {
  return spsp_mm_mult<srmatrix,scmatrix,scmatrix,sparse_cdot,complex>(A,B);
}

inline scmatrix operator*(const scmatrix& A, const scmatrix& B) {
  return spsp_mm_mult<scmatrix,scmatrix,scmatrix,sparse_cdot,complex>(A,B);
}

inline scmatrix operator/(const scmatrix& A, const real& r) {
  return sp_ms_div<scmatrix,real,scmatrix>(A,r);
}

inline scmatrix operator/(const scmatrix& A, const complex& r) {
  return sp_ms_div<scmatrix,complex,scmatrix>(A,r);
}

inline scmatrix operator/(const srmatrix& A, const complex& r) {
  return sp_ms_div<srmatrix,complex,scmatrix>(A,r);
}

inline scmatrix operator*(const scmatrix& A, const real& r) {
  return sp_ms_mult<scmatrix,real,scmatrix>(A,r);
}

inline scmatrix operator*(const scmatrix& A, const complex& r) {
  return sp_ms_mult<scmatrix,complex,scmatrix>(A,r);
}

inline scmatrix operator*(const srmatrix& A, const complex& r) {
  return sp_ms_mult<srmatrix,complex,scmatrix>(A,r);
}

inline scmatrix operator*(const real& r, const scmatrix& A) {
  return sp_sm_mult<real,scmatrix,scmatrix>(r,A);
}

inline scmatrix operator*(const complex& r, const scmatrix& A) {
  return sp_sm_mult<complex,scmatrix,scmatrix>(r,A);
}

inline scmatrix operator*(const complex& r, const srmatrix& A) {
  return sp_sm_mult<complex,srmatrix,scmatrix>(r,A);
}

inline cvector operator*(const scmatrix& A, const rvector& v) {
  return spf_mv_mult<scmatrix,rvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const srmatrix& A, const cvector& v) {
  return spf_mv_mult<srmatrix,cvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const scmatrix& A, const cvector& v) {
  return spf_mv_mult<scmatrix,cvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const scmatrix& A, const rvector_slice& v) {
  return spf_mv_mult<scmatrix,rvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const srmatrix& A, const cvector_slice& v) {
  return spf_mv_mult<srmatrix,cvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const scmatrix& A, const cvector_slice& v) {
  return spf_mv_mult<scmatrix,cvector_slice,cvector,sparse_cdot>(A,v);
}

inline scvector operator*(const scmatrix& A, const srvector& v) {
  return spsp_mv_mult<scmatrix,srvector,scvector,sparse_cdot,complex>(A,v);
}

inline scvector operator*(const srmatrix& A, const scvector& v) {
  return spsp_mv_mult<srmatrix,scvector,scvector,sparse_cdot,complex>(A,v);
}

inline scvector operator*(const scmatrix& A, const scvector& v) {
  return spsp_mv_mult<scmatrix,scvector,scvector,sparse_cdot,complex>(A,v);
}

inline scvector operator*(const scmatrix& A, const srvector_slice& v) {
  return spsl_mv_mult<scmatrix,srvector_slice,scvector,sparse_cdot,complex>(A,v);
}

inline scvector operator*(const srmatrix& A, const scvector_slice& v) {
  return spsl_mv_mult<srmatrix,scvector_slice,scvector,sparse_cdot,complex>(A,v);
}

inline scvector operator*(const scmatrix& A, const scvector_slice& v) {
  return spsl_mv_mult<scmatrix,scvector_slice,scvector,sparse_cdot,complex>(A,v);
}

inline cvector operator*(const cmatrix& A, const srvector& v) {
  return fsp_mv_mult<cmatrix,srvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const rmatrix& A, const scvector& v) {
  return fsp_mv_mult<rmatrix,scvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix& A, const scvector& v) {
  return fsp_mv_mult<cmatrix,scvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix_slice& A, const srvector& v) {
  return fsp_mv_mult<cmatrix_slice,srvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const rmatrix_slice& A, const scvector& v) {
  return fsp_mv_mult<rmatrix_slice,scvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix_slice& A, const scvector& v) {
  return fsp_mv_mult<cmatrix_slice,scvector,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix& A, const srvector_slice& v) {
  return fsl_mv_mult<cmatrix,srvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const rmatrix& A, const scvector_slice& v) {
  return fsl_mv_mult<rmatrix,scvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix& A, const scvector_slice& v) {
  return fsl_mv_mult<cmatrix,scvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix_slice& A, const srvector_slice& v) {
  return fsl_mv_mult<cmatrix_slice,srvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const rmatrix_slice& A, const scvector_slice& v) {
  return fsl_mv_mult<rmatrix_slice,scvector_slice,cvector,sparse_cdot>(A,v);
}

inline cvector operator*(const cmatrix_slice& A, const scvector_slice& v) {
  return fsl_mv_mult<cmatrix_slice,scvector_slice,cvector,sparse_cdot>(A,v);
}

inline cmatrix operator+(const cmatrix& A, const srmatrix& B) {
  return fsp_mm_add<cmatrix,srmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const rmatrix& A, const scmatrix& B) {
  return fsp_mm_add<rmatrix,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const cmatrix& A, const scmatrix& B) {
  return fsp_mm_add<cmatrix,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const scmatrix& A, const rmatrix& B) {
  return spf_mm_add<scmatrix,rmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const srmatrix& A, const cmatrix& B) {
  return spf_mm_add<srmatrix,cmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const scmatrix& A, const cmatrix& B) {
  return spf_mm_add<scmatrix,cmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const cmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_add<cmatrix_slice,srmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const rmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_add<rmatrix_slice,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const cmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_add<cmatrix_slice,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator+(const scmatrix& A, const rmatrix_slice& B) {
  return spf_mm_add<scmatrix,rmatrix_slice,cmatrix>(A,B);
}

inline cmatrix operator+(const srmatrix& A, const cmatrix_slice& B) {
  return spf_mm_add<srmatrix,cmatrix_slice,cmatrix>(A,B);
}

inline cmatrix operator+(const scmatrix& A, const cmatrix_slice& B) {
  return spf_mm_add<scmatrix,cmatrix_slice,cmatrix>(A,B);
}

inline scmatrix operator+(const scmatrix& A, const srmatrix& B) {
  return spsp_mm_add<scmatrix,srmatrix,scmatrix,complex>(A,B);
}

inline scmatrix operator+(const srmatrix& A, const scmatrix& B) {
  return spsp_mm_add<srmatrix,scmatrix,scmatrix,complex>(A,B);
}

inline scmatrix operator+(const scmatrix& A, const scmatrix& B) {
  return spsp_mm_add<scmatrix,scmatrix,scmatrix,complex>(A,B);
}

inline cmatrix operator-(const cmatrix& A, const srmatrix& B) {
  return fsp_mm_sub<cmatrix,srmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const rmatrix& A, const scmatrix& B) {
  return fsp_mm_sub<rmatrix,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const cmatrix& A, const scmatrix& B) {
  return fsp_mm_sub<cmatrix,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const scmatrix& A, const rmatrix& B) {
  return spf_mm_sub<scmatrix,rmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const srmatrix& A, const cmatrix& B) {
  return spf_mm_sub<srmatrix,cmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const scmatrix& A, const cmatrix& B) {
  return spf_mm_sub<scmatrix,cmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const cmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_sub<cmatrix_slice,srmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const rmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_sub<rmatrix_slice,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const cmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_sub<cmatrix_slice,scmatrix,cmatrix>(A,B);
}

inline cmatrix operator-(const scmatrix& A, const rmatrix_slice& B) {
  return spf_mm_sub<scmatrix,rmatrix_slice,cmatrix>(A,B);
}

inline cmatrix operator-(const srmatrix& A, const cmatrix_slice& B) {
  return spf_mm_sub<srmatrix,cmatrix_slice,cmatrix>(A,B);
}

inline cmatrix operator-(const scmatrix& A, const cmatrix_slice& B) {
  return spf_mm_sub<scmatrix,cmatrix_slice,cmatrix>(A,B);
}

inline scmatrix operator-(const scmatrix& A, const srmatrix& B) {
  return spsp_mm_sub<scmatrix,srmatrix,scmatrix,complex>(A,B);
}

inline scmatrix operator-(const srmatrix& A, const scmatrix& B) {
  return spsp_mm_sub<srmatrix,scmatrix,scmatrix,complex>(A,B);
}

inline scmatrix operator-(const scmatrix& A, const scmatrix& B) {
  return spsp_mm_sub<scmatrix,scmatrix,scmatrix,complex>(A,B);
}

inline scmatrix operator-(const scmatrix& M) {
  return sp_m_negative<scmatrix,scmatrix>(M);
}

inline scmatrix& operator+(scmatrix& A) {
  return A;
}

inline cmatrix& cmatrix::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline cmatrix& cmatrix::operator=(const scmatrix& B) {
  *this = cmatrix(B);
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator=(const scmatrix& B) {
  *this = cmatrix(B);
  return *this;
}

inline cmatrix& cmatrix::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cmatrix& cmatrix::operator+=(const scmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator+=(const scmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cmatrix& cmatrix::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cmatrix& cmatrix::operator-=(const scmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator-=(const scmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cmatrix& cmatrix::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<cmatrix,srmatrix,sparse_cdot,cmatrix>(*this,B);
}

inline cmatrix& cmatrix::operator*=(const scmatrix& B) {
  return fsp_mm_multassign<cmatrix,scmatrix,sparse_cdot,cmatrix>(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<cmatrix_slice,srmatrix,sparse_cdot,cmatrix>(*this,B);
}

inline cmatrix_slice& cmatrix_slice::operator*=(const scmatrix& B) {
  return fsp_mm_multassign<cmatrix_slice,scmatrix,sparse_cdot,cmatrix>(*this,B);
}

inline bool operator==(const scmatrix& A, const srmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const scmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const scmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const rmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const cmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const cmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const cmatrix& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cmatrix& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const rmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const cmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const cmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const srmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const scmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const scmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const rmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const cmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const cmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const cmatrix& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cmatrix& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cmatrix_slice& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix_slice& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cmatrix_slice& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const rmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const cmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const cmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!(const scmatrix& A) {
  return sp_m_not(A);
}

inline std::ostream& operator<<(std::ostream& os, const scmatrix& A) {
  return sp_m_output<scmatrix,complex>(os,A);
}

inline std::istream& operator>>(std::istream& is, scmatrix& A) {
  return sp_m_input<scmatrix,complex>(is,A);
}

class scmatrix_slice {
  public:
    scmatrix  A;
    scmatrix* M; //Originalmatrix

  private:
    scmatrix_slice(scmatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
        A.lb1 = sl1l;
        A.lb2 = sl2l;
        A.ub1 = sl1u;
        A.ub2 = sl2u;
        A.m   = sl1u-sl1l+1;
        A.n   = sl2u-sl2l+1;
 
        //Kopieren der Werte aus A
        A.p = std::vector<int>(A.n+1, 0);
        A.ind.reserve(A.m + A.n);
        A.x.reserve(A.m + A.n);

        for(int i=0 ; i<A.n ; i++) {
           A.p[i+1] = A.p[i];
           for(int j=Mat.p[sl2l-Mat.lb2+i] ; j<Mat.p[sl2l-Mat.lb2+i+1] ; j++) {
              if(Mat.ind[j] >= sl1l-Mat.lb1  &&  Mat.ind[j] <= sl1u-Mat.lb1) {
                A.ind.push_back(Mat.ind[j]-(sl1l-Mat.lb1));
                A.x.push_back(Mat.x[j]);
                A.p[i+1]++;
             }
           }
        }

        //Zeiger auf A fuer Datenmanipulationen
        M = &Mat;
    }

    scmatrix_slice(const scmatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
        A.lb1 = sl1l;
        A.lb2 = sl2l;
        A.ub1 = sl1u;
        A.ub2 = sl2u;
        A.m   = sl1u-sl1l+1;
        A.n   = sl2u-sl2l+1;
 
        //Kopieren der Werte aus A
        A.p = std::vector<int>(A.n+1, 0);
        A.ind.reserve(A.m + A.n);
        A.x.reserve(A.m + A.n);

        for(int i=0 ; i<A.n ; i++) {
           A.p[i+1] = A.p[i];
           for(int j=Mat.p[sl2l-Mat.lb2+i] ; j<Mat.p[sl2l-Mat.lb2+i+1] ; j++) {
              if(Mat.ind[j] >= sl1l-Mat.lb1  &&  Mat.ind[j] <= sl1u-Mat.lb1) {
                A.ind.push_back(Mat.ind[j]-(sl1l-Mat.lb1));
                A.x.push_back(Mat.x[j]);
                A.p[i+1]++;
             }
           }
        }

        //Zeiger auf A fuer Datenmanipulationen
        M = const_cast<scmatrix*>(&Mat);
    }

  public:
    scmatrix_slice& operator=(const real& C) {
      return sl_ms_assign<scmatrix_slice, real, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const complex& C) {
      return sl_ms_assign<scmatrix_slice, complex, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const srmatrix& C) {
      return slsp_mm_assign<scmatrix_slice, srmatrix, std::vector<complex>::iterator>(*this,C);
    }

    scmatrix_slice& operator=(const scmatrix& C) {
      return slsp_mm_assign<scmatrix_slice, scmatrix, std::vector<complex>::iterator>(*this,C);
    }

    scmatrix_slice& operator=(const rmatrix& C) {
      return slf_mm_assign<scmatrix_slice, rmatrix, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const cmatrix& C) {
      return slf_mm_assign<scmatrix_slice, cmatrix, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const rmatrix_slice& C) {
      return slf_mm_assign<scmatrix_slice, rmatrix_slice, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const cmatrix_slice& C) {
      return slf_mm_assign<scmatrix_slice, cmatrix_slice, std::vector<complex>::iterator, complex>(*this,C);
    }

    scmatrix_slice& operator=(const srmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scmatrix_slice& operator=(const scmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scmatrix_slice& operator*=(const srmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scmatrix_slice& operator*=(const scmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scmatrix_slice& operator*=(const srmatrix& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const scmatrix& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const rmatrix& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const cmatrix& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const rmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const cmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scmatrix_slice& operator*=(const real& r) {
      *this = A*r;
      return *this;
    }

    scmatrix_slice& operator*=(const complex& r) {
      *this = A*r;
      return *this;
    }

    scmatrix_slice& operator/=(const real& r) {
      *this = A/r;
      return *this;
    }

    scmatrix_slice& operator/=(const complex& r) {
      *this = A/r;
      return *this;
    }

    scmatrix_slice& operator+=(const srmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scmatrix_slice& operator+=(const scmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scmatrix_slice& operator+=(const srmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator+=(const scmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator+=(const rmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator+=(const cmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator+=(const rmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator+=(const cmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scmatrix_slice& operator-=(const srmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scmatrix_slice& operator-=(const scmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scmatrix_slice& operator-=(const srmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scmatrix_slice& operator-=(const scmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scmatrix_slice& operator-=(const rmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scmatrix_slice& operator-=(const cmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scmatrix_slice& operator-=(const rmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    scmatrix_slice& operator-=(const cmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    const complex operator()(const int i, const int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_slice::operator()(int, int)"));
#endif
      complex r = A(i,j);
      return r;
    }

    complex& element(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_slice::element(int, int)"));
#endif
      return M->element(i,j);
    }

    scmatrix_subv operator[](const int);
    scmatrix_subv operator[](const cxscmatrix_column&);
    const scmatrix_subv operator[](const int) const;
    const scmatrix_subv operator[](const cxscmatrix_column&) const;

    friend std::ostream& operator<<(std::ostream&, const scmatrix_slice&);
    friend std::istream& operator>>(std::istream&, const scmatrix_subv&);

    friend int Lb(const scmatrix_slice&, const int);
    friend int Ub(const scmatrix_slice&, const int);
    friend srmatrix Re(const scmatrix_slice&);
    friend srmatrix Im(const scmatrix_slice&);
    friend int RowLen(const scmatrix_slice&);
    friend int ColLen(const scmatrix_slice&);

    friend class srmatrix;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix;
    friend class scmatrix_subv;
    friend class scvector;
    friend class scimatrix;
    friend class scimatrix_subv;
    friend class scimatrix_slice;
    friend class scivector;
    friend class cmatrix;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"    
};

inline cmatrix::cmatrix(const srmatrix_slice& A) {
  dat = new complex[A.A.m*A.A.n];
  lb1 = A.A.lb1; lb2 = A.A.lb2; ub1 = A.A.ub1; ub2 = A.A.ub2;
  xsize = A.A.n;
  ysize = A.A.m;
  *this = 0.0;
  for(int j=0 ; j<A.A.n ; j++) {
     for(int k=A.A.p[j] ; k<A.A.p[j+1] ; k++) {
        dat[A.A.ind[k]*A.A.n+j] = A.A.x[k];
     }
  }
}

inline cmatrix::cmatrix(const scmatrix_slice& A) {
  dat = new complex[A.A.m*A.A.n];
  lb1 = A.A.lb1; lb2 = A.A.lb2; ub1 = A.A.ub1; ub2 = A.A.ub2;
  xsize = A.A.n;
  ysize = A.A.m;
  *this = 0.0;
  for(int j=0 ; j<A.A.n ; j++) {
     for(int k=A.A.p[j] ; k<A.A.p[j+1] ; k++) {
        dat[A.A.ind[k]*A.A.n+j] = A.A.x[k];
     }
  }
}

inline int RowLen(const scmatrix_slice& S) {
  return RowLen(S.A);
}

inline int ColLen(const scmatrix_slice& S) {
  return ColLen(S.A);
}

inline scmatrix_slice scmatrix::operator()(const int i, const int j, const int k, const int l) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator()(int, int)"));
#endif
  return scmatrix_slice(*this, i, j, k, l);
}

inline scmatrix& scmatrix::operator=(const srmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline scmatrix& scmatrix::operator=(const scmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline int Lb(const scmatrix_slice& S, const int i) {
  return Lb(S.A, i);
}

inline int Ub(const scmatrix_slice& S, const int i) {
  return Ub(S.A, i);
}

inline srmatrix Re(const scmatrix_slice& S) {
  return Re(S.A);
}

inline srmatrix Im(const scmatrix_slice& S) {
  return Im(S.A);
}

inline scmatrix::scmatrix(const srmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scmatrix::scmatrix(const scmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scmatrix operator-(const scmatrix_slice& M) {
  return sp_m_negative<scmatrix,scmatrix>(M.A);
}

inline scmatrix operator+(const scmatrix_slice& M) {
  return M.A;
}

inline scmatrix operator*(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,srmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2.A);
}

inline scmatrix operator*(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2.A);
}

inline scmatrix operator*(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2.A);
}

inline scmatrix operator*(const scmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_mult<scmatrix,srmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2);
}

inline scmatrix operator*(const srmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_mult<srmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2);
}

inline scmatrix operator*(const scmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_mult<scmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1.A,M2);
}

inline scmatrix operator*(const scmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,srmatrix,scmatrix,sparse_cdot,complex>(M1,M2.A);
}

inline scmatrix operator*(const srmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1,M2.A);
}

inline scmatrix operator*(const scmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,scmatrix,scmatrix,sparse_cdot,complex>(M1,M2.A);
}

inline cmatrix operator*(const scmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_mult<scmatrix,rmatrix,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const srmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_mult<srmatrix,cmatrix,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const scmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_mult<scmatrix,cmatrix,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const cmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,srmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline cmatrix operator*(const rmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,scmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline cmatrix operator*(const cmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,scmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline cmatrix operator*(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_mult<scmatrix,rmatrix_slice,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_mult<srmatrix,cmatrix_slice,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_mult<scmatrix,cmatrix_slice,cmatrix,sparse_cdot>(M1.A,M2);
}

inline cmatrix operator*(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,srmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline cmatrix operator*(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,scmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline cmatrix operator*(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,scmatrix,cmatrix,sparse_cdot>(M1,M2.A);
}

inline scvector operator*(const scmatrix_slice& M, const srvector& v) {
  return spsp_mv_mult<scmatrix,srvector,scvector,sparse_cdot,complex>(M.A,v);
}

inline scvector operator*(const srmatrix_slice& M, const scvector& v) {
  return spsp_mv_mult<srmatrix,scvector,scvector,sparse_cdot,complex>(M.A,v);
}

inline scvector operator*(const scmatrix_slice& M, const scvector& v) {
  return spsp_mv_mult<scmatrix,scvector,scvector,sparse_cdot,complex>(M.A,v);
}

inline scvector operator*(const scmatrix_slice& M, const srvector_slice& v) {
  return spsl_mv_mult<scmatrix,srvector_slice,scvector,sparse_cdot,complex>(M.A,v);
}

inline scvector operator*(const srmatrix_slice& M, const scvector_slice& v) {
  return spsl_mv_mult<srmatrix,scvector_slice,scvector,sparse_cdot,complex>(M.A,v);
}

inline scvector operator*(const scmatrix_slice& M, const scvector_slice& v) {
  return spsl_mv_mult<scmatrix,scvector_slice,scvector,sparse_cdot,complex>(M.A,v);
}

inline cvector operator*(const scmatrix_slice& M, const rvector& v) {
  return spf_mv_mult<scmatrix,rvector,cvector,sparse_cdot>(M.A,v);
}

inline cvector operator*(const srmatrix_slice& M, const cvector& v) {
  return spf_mv_mult<srmatrix,cvector,cvector,sparse_cdot>(M.A,v);
}

inline cvector operator*(const scmatrix_slice& M, const cvector& v) {
  return spf_mv_mult<scmatrix,cvector,cvector,sparse_cdot>(M.A,v);
}

inline cvector operator*(const scmatrix_slice& M, const rvector_slice& v) {
  return spf_mv_mult<scmatrix,rvector_slice,cvector,sparse_cdot>(M.A,v);
}

inline cvector operator*(const srmatrix_slice& M, const cvector_slice& v) {
  return spf_mv_mult<srmatrix,cvector_slice,cvector,sparse_cdot>(M.A,v);
}

inline cvector operator*(const scmatrix_slice& M, const cvector_slice& v) {
  return spf_mv_mult<scmatrix,cvector_slice,cvector,sparse_cdot>(M.A,v);
}

inline scmatrix operator*(const scmatrix_slice& M, const real& r) {
  return sp_ms_mult<scmatrix,real,scmatrix>(M.A,r);
}

inline scmatrix operator*(const scmatrix_slice& M, const complex& r) {
  return sp_ms_mult<scmatrix,complex,scmatrix>(M.A,r);
}

inline scmatrix operator*(const srmatrix_slice& M, const complex& r) {
  return sp_ms_mult<srmatrix,complex,scmatrix>(M.A,r);
}

inline scmatrix operator/(const scmatrix_slice& M, const real& r) {
  return sp_ms_div<scmatrix,real,scmatrix>(M.A,r);
}

inline scmatrix operator/(const scmatrix_slice& M, const complex& r) {
  return sp_ms_div<scmatrix,complex,scmatrix>(M.A,r);
}

inline scmatrix operator/(const srmatrix_slice& M, const complex& r) {
  return sp_ms_div<srmatrix,complex,scmatrix>(M.A,r);
}

inline scmatrix operator*(const real& r, const scmatrix_slice& M) {
  return sp_sm_mult<real,scmatrix,scmatrix>(r,M.A);
}

inline scmatrix operator*(const complex& r, const srmatrix_slice& M) {
  return sp_sm_mult<complex,srmatrix,scmatrix>(r,M.A);
}

inline scmatrix operator*(const complex& r, const scmatrix_slice& M) {
  return sp_sm_mult<complex,scmatrix,scmatrix>(r,M.A);
}

inline scmatrix operator+(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<scmatrix,srmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator+(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<srmatrix,scmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator+(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<scmatrix,scmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator+(const scmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_add<scmatrix,srmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator+(const srmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_add<srmatrix,scmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator+(const scmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_add<scmatrix,scmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator+(const scmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<scmatrix,srmatrix,scmatrix,complex>(M1,M2.A);
}

inline scmatrix operator+(const srmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<srmatrix,scmatrix,scmatrix,complex>(M1,M2.A);
}

inline scmatrix operator+(const scmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<scmatrix,scmatrix,scmatrix,complex>(M1,M2.A);
}

inline cmatrix operator+(const scmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_add<scmatrix,rmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const srmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_add<srmatrix,cmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const scmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_add<scmatrix,cmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const cmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<cmatrix,srmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator+(const rmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<rmatrix,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator+(const cmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<cmatrix,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator+(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_add<scmatrix,rmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_add<srmatrix,cmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_add<scmatrix,cmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator+(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<cmatrix_slice,srmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator+(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<rmatrix_slice,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator+(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<cmatrix_slice,scmatrix,cmatrix>(M1,M2.A);
}

inline scmatrix operator-(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,srmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator-(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,scmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator-(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,scmatrix,scmatrix,complex>(M1.A,M2.A);
}

inline scmatrix operator-(const scmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_sub<scmatrix,srmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator-(const srmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_sub<srmatrix,scmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator-(const scmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_sub<scmatrix,scmatrix,scmatrix,complex>(M1.A,M2);
}

inline scmatrix operator-(const scmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,srmatrix,scmatrix,complex>(M1,M2.A);
}

inline scmatrix operator-(const srmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,scmatrix,scmatrix,complex>(M1,M2.A);
}

inline scmatrix operator-(const scmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,scmatrix,scmatrix,complex>(M1,M2.A);
}

inline cmatrix operator-(const scmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_sub<scmatrix,rmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const srmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_sub<srmatrix,cmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const scmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_sub<scmatrix,cmatrix,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const cmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<cmatrix,srmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator-(const rmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<rmatrix,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator-(const cmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<cmatrix,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator-(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_sub<scmatrix,rmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_sub<srmatrix,cmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_sub<scmatrix,cmatrix_slice,cmatrix>(M1.A,M2);
}

inline cmatrix operator-(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<cmatrix_slice,srmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator-(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<rmatrix_slice,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix operator-(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<cmatrix_slice,scmatrix,cmatrix>(M1,M2.A);
}

inline cmatrix& cmatrix::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline cmatrix& cmatrix::operator=(const scmatrix_slice& M) {
  *this = cmatrix(M);
  return *this;
}

inline cmatrix& cmatrix::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cmatrix& cmatrix::operator+=(const scmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator+=(const scmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cmatrix& cmatrix::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cmatrix& cmatrix::operator-=(const scmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator-=(const scmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cmatrix& cmatrix::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cmatrix& cmatrix::operator*=(const scmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cmatrix_slice& cmatrix_slice::operator*=(const scmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline bool operator==(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const srmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const cmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scmatrix_slice& M1, const srmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const scmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const scmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const srmatrix& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scmatrix& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scmatrix_slice& M1, const rmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const cmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const cmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const cmatrix& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cmatrix& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!(const scmatrix_slice& M) {
  return sp_m_not(M.A);
}

inline std::ostream& operator<<(std::ostream& os, const scmatrix_slice& M) {
  return sp_m_output<scmatrix,complex>(os, M.A);
}

inline std::istream& operator>>(std::istream& is, scmatrix_slice& M) {
  scmatrix tmp(M.A.m, M.A.n);
  is >> tmp;
  M = tmp;
  return is;
}


class scmatrix_subv {
  private:
    scmatrix_slice dat;
    bool row;
    int index;

    scmatrix_subv(scmatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

    scmatrix_subv(const scmatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

  
  public:
    complex& operator[](const int i) {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_subv::operator[](int)"));
#endif
        return dat.element(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_subv::operator[](int)"));
#endif
        return dat.element(i,index);
      }
    }

    const complex operator[](const int i) const {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_subv::operator[](int)"));
#endif
        return dat(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("scmatrix_subv::operator[](int)"));
#endif
        return dat(i,index);
      }
    }

    scmatrix_subv& operator=(const real& v) {
      return sv_vs_assign(*this,v);
    }

    scmatrix_subv& operator=(const complex& v) {
      return sv_vs_assign(*this,v);
    }

    scmatrix_subv& operator=(const srvector& v) {
      return svsp_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const scvector& v) {
      return svsp_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const srvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const scvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const rvector& v) {
      return svf_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const cvector& v) {
      return svf_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const rvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const cvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scmatrix_subv& operator=(const srmatrix_subv& v) {
      return svsp_vv_assign(*this,srvector(v));
    }

    scmatrix_subv& operator=(const scmatrix_subv& v) {
      return svsp_vv_assign(*this,scvector(v));
    }


    scmatrix_subv& operator*=(const real&);
    scmatrix_subv& operator*=(const complex&);
    scmatrix_subv& operator/=(const real&);
    scmatrix_subv& operator/=(const complex&);
    scmatrix_subv& operator+=(const srvector&);
    scmatrix_subv& operator+=(const srvector_slice&);
    scmatrix_subv& operator+=(const rvector&);
    scmatrix_subv& operator+=(const rvector_slice&);
    scmatrix_subv& operator-=(const srvector&);
    scmatrix_subv& operator-=(const srvector_slice&);
    scmatrix_subv& operator-=(const rvector&);
    scmatrix_subv& operator-=(const rvector_slice&);
    scmatrix_subv& operator+=(const scvector&);
    scmatrix_subv& operator+=(const scvector_slice&);
    scmatrix_subv& operator+=(const cvector&);
    scmatrix_subv& operator+=(const cvector_slice&);
    scmatrix_subv& operator-=(const scvector&);
    scmatrix_subv& operator-=(const scvector_slice&);
    scmatrix_subv& operator-=(const cvector&);
    scmatrix_subv& operator-=(const cvector_slice&);

    friend scvector operator-(const scmatrix_subv&);
    friend std::istream& operator>>(std::istream&, scmatrix_subv&);

    friend int Lb(const scmatrix_subv&);
    friend int Ub(const scmatrix_subv&);
    friend int VecLen(const scmatrix_subv&);
    friend srvector Re(const scmatrix_subv&);
    friend srvector Im(const scmatrix_subv&);

    friend class srvector;
    friend class srmatrix;
    friend class srmatrix_slice;
    friend class scvector;
    friend class scmatrix;
    friend class scmatrix_slice;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;

#include "vector_friend_declarations.inl"
};

inline int Lb(const scmatrix_subv& S) {
  if(S.row)
    return Lb(S.dat, 2);
  else
    return Lb(S.dat, 1);
}

inline int Ub(const scmatrix_subv& S) {
  if(S.row)
    return Ub(S.dat, 2);
  else
    return Ub(S.dat, 1);
}

inline int VecLen(const scmatrix_subv& S) {
  return Ub(S)-Lb(S)+1;
}

inline srvector Re(const scmatrix_subv& S) {
  return Re(scvector(S));
}

inline srvector Im(const scmatrix_subv& S) {
  return Im(scvector(S));
}

inline std::ostream& operator<<(std::ostream& os, const scmatrix_subv& v) {
  os << scvector(v);
  return os;
}

inline std::istream& operator>>(std::istream& is, scmatrix_subv& v) {
  int n=0;
  if(v.row) n=v.dat.A.n; else n=v.dat.A.m;
  scvector tmp(n);
  is >> tmp;
  v = tmp;
  return is;
}

inline scmatrix_subv scmatrix::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator[](const cxscmatrix_column&)"));
#endif
  return scmatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline scmatrix_subv scmatrix::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator[](const int)"));
#endif
  return scmatrix_subv(*this, true, i, i, lb2, ub2);
}

inline const scmatrix_subv scmatrix::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator[](const cxscmatrix_column&)"));
#endif
  return scmatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline const scmatrix_subv scmatrix::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix::operator[](const int)"));
#endif
  return scmatrix_subv(*this, true, i, i, lb2, ub2);
}

inline scmatrix_subv scmatrix_slice::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix_slice::operator[](const int)"));
#endif
  return scmatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline scmatrix_subv scmatrix_slice::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return scmatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline const scmatrix_subv scmatrix_slice::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix_slice::operator[](const int)"));
#endif
  return scmatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline const scmatrix_subv scmatrix_slice::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scmatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return scmatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline scvector::scvector(const scmatrix_subv& A) {
  int nnz = A.dat.A.get_nnz();
  p.reserve(nnz);
  x.reserve(nnz);

  if(A.row) {
    lb = A.dat.A.lb2;
    ub = A.dat.A.ub2;
    n = ub-lb+1; 

    for(int j=0 ; j<n ; j++) {
      for(int k=A.dat.A.p[j] ; k<A.dat.A.p[j+1] ; k++) {
        p.push_back(j);
        x.push_back(A.dat.A.x[k]);
      }
    }

  } else {
    lb = A.dat.A.lb1;
    ub = A.dat.A.ub1;
    n = ub-lb+1; 

    for(unsigned int k=0 ; k<A.dat.A.ind.size() ; k++) {
        p.push_back(A.dat.A.ind[k]);
        x.push_back(A.dat.A.x[k]);
    }
  }
}

inline scvector operator-(const scmatrix_subv& v) {
 scvector s(v);
 return -s;
}

inline scvector operator/(const scmatrix_subv& v1, const real& v2) {
  return scvector(v1) / v2;
}

inline scvector operator/(const scmatrix_subv& v1, const complex& v2) {
  return scvector(v1) / v2;
}

inline scvector operator/(const srmatrix_subv& v1, const complex& v2) {
  return srvector(v1) / v2;
}

inline scvector operator*(const scmatrix_subv& v1, const real& v2) {
  return scvector(v1) * v2;
}

inline scvector operator*(const scmatrix_subv& v1, const complex& v2) {
  return scvector(v1) * v2;
}

inline scvector operator*(const srmatrix_subv& v1, const complex& v2) {
  return srvector(v1) * v2;
}

inline scvector operator*(const real& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline scvector operator*(const complex& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline scvector operator*(const complex& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline complex operator*(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) * v2;
}

inline complex operator*(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) * v2;
}

inline complex operator*(const scvector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline complex operator*(const srvector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const scvector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline complex operator*(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const cvector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline complex operator*(const rvector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const cvector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline complex operator*(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline complex operator*(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline scvector operator+(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) + v2;
}

inline scvector operator+(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) + v2;
}

inline scvector operator+(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) + v2;
}

inline scvector operator+(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) + v2;
}

inline scvector operator+(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) + v2;
}

inline scvector operator+(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) + v2;
}

inline cvector operator+(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) + v2;
}

inline cvector operator+(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) + v2;
}

inline cvector operator+(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) + v2;
}

inline cvector operator+(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) + v2;
}

inline cvector operator+(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) + v2;
}

inline cvector operator+(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) + v2;
}

inline scvector operator+(const scvector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline scvector operator+(const srvector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scvector operator+(const scvector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scvector operator+(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline scvector operator+(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scvector operator+(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline cvector operator+(const cvector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline cvector operator+(const rvector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline cvector operator+(const cvector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline cvector operator+(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline cvector operator+(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline cvector operator+(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scvector operator-(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) - v2;
}

inline scvector operator-(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) - v2;
}

inline scvector operator-(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) - v2;
}

inline scvector operator-(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) - v2;
}

inline scvector operator-(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) - v2;
}

inline scvector operator-(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) - v2;
}

inline cvector operator-(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) - v2;
}

inline cvector operator-(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) - v2;
}

inline cvector operator-(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) - v2;
}

inline cvector operator-(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) - v2;
}

inline cvector operator-(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) - v2;
}

inline cvector operator-(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) - v2;
}

inline scvector operator-(const scvector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline scvector operator-(const srvector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scvector operator-(const scvector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scvector operator-(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline scvector operator-(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scvector operator-(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline cvector operator-(const cvector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline cvector operator-(const rvector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline cvector operator-(const cvector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline cvector operator-(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline cvector operator-(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline cvector operator-(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scmatrix_subv& scmatrix_subv::operator*=(const real& v) {
  *this = *this * v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator*=(const complex& v) {
  *this = *this * v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator/=(const real& v) {
  *this = *this / v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator/=(const complex& v) {
  *this = *this / v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const srvector& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const srvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const rvector& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const rvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const srvector& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const srvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const rvector& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const rvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const scvector& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const scvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const cvector& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator+=(const cvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const scvector& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const scvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const cvector& v) {
  *this = *this - v;
  return *this;
}

inline scmatrix_subv& scmatrix_subv::operator-=(const cvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const srmatrix_subv& v) {
  *this += rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const scmatrix_subv& v) {
  *this += cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const srvector& v) {
  *this += rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const scvector& v) {
  *this += cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const srvector_slice& v) {
  *this += rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator+=(const scvector_slice& v) {
  *this += cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator-=(const srmatrix_subv& v) {
  *this -= rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator-=(const scmatrix_subv& v) {
  *this -= cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const scvector& v) {
  *this = cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const scvector_slice& v) {
  *this = cvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const srmatrix_subv& v) {
  *this = rvector(v);
  return *this;
}

inline cmatrix_subv& cmatrix_subv::operator=(const scmatrix_subv& v) {
  *this = cvector(v);
  return *this;
}

inline bool operator==(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const scvector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const srvector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const scvector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const cvector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const cvector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator!=(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const scvector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const srvector& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const scvector& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const cvector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const cvector& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!(const scmatrix_subv& x) {
  return sv_v_not(x);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const scmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,scvector,scvector,sparse_cdot>(dot, scvector(v1), scvector(v2));
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,scvector,srvector,sparse_cdot>(dot, scvector(v1), srvector(v2));
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const scmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,srvector,scvector,sparse_cdot>(dot, srvector(v1), scvector(v2));
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const scvector& v2) {
  spsp_vv_accu<cdotprecision,scvector,scvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const srvector& v2) {
  spsp_vv_accu<cdotprecision,scvector,srvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const scvector& v2) {
  spsp_vv_accu<cdotprecision,srvector,scvector,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const scvector_slice& v2) {
  spsl_vv_accu<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accu<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const scvector_slice& v2) {
  spsl_vv_accu<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const cvector& v2) {
  spf_vv_accu<cdotprecision,scvector,cvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const rvector& v2) {
  spf_vv_accu<cdotprecision,scvector,rvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const cvector& v2) {
  spf_vv_accu<cdotprecision,srvector,cvector,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const cvector_slice& v2) {
  spf_vv_accu<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scmatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accu<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const cvector_slice& v2) {
  spf_vv_accu<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate(cdotprecision& dot, const scvector& v1, const scmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,scvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const scvector& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,scvector,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate(cdotprecision& dot, const srvector& v1, const scmatrix_subv& v2) {
  spsp_vv_accu<cdotprecision,srvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const scvector_slice& v1, const scmatrix_subv& v2) {
  slsp_vv_accu<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const scvector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accu<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate(cdotprecision& dot, const srvector_slice& v1, const scmatrix_subv& v2) {
  slsp_vv_accu<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const cvector& v1, const scmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,cvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const cvector& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,cvector,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate(cdotprecision& dot, const rvector& v1, const scmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,rvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const cvector_slice& v1, const scmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cdotprecision& dot, const cvector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate(cdotprecision& dot, const rvector_slice& v1, const scmatrix_subv& v2) {
  fsp_vv_accu<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const scmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,scvector,sparse_cdot>(dot, scvector(v1), scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,srvector,sparse_cdot>(dot, scvector(v1), srvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const scmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,srvector,scvector,sparse_cdot>(dot, srvector(v1), scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const scvector& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,scvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const srvector& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,srvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const scvector& v2) {
  spsp_vv_accuapprox<cdotprecision,srvector,scvector,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const scvector_slice& v2) {
  spsl_vv_accuapprox<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accuapprox<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const scvector_slice& v2) {
  spsl_vv_accuapprox<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const cvector& v2) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const rvector& v2) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const cvector& v2) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const cvector_slice& v2) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scmatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot, scvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const cvector_slice& v2) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& v1, const scmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scvector& v1, const srmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,scvector,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const srvector& v1, const scmatrix_subv& v2) {
  spsp_vv_accuapprox<cdotprecision,srvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& v1, const scmatrix_subv& v2) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& v1, const scmatrix_subv& v2) {
  slsp_vv_accuapprox<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const cvector& v1, const scmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,cvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const cvector& v1, const srmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,cvector,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const rvector& v1, const scmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,rvector,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& v1, const scmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(cdotprecision& dot, const rvector_slice& v1, const scmatrix_subv& v2) {
  fsp_vv_accuapprox<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot, v1, scvector(v2));
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const srmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const scvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const srvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const scvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const scvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const srvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const scvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const cvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const rvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const cvector& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const cvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scmatrix_subv& v1, const rvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,scvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const cvector_slice& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const srmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const srmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const srmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const srmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const scmatrix_subv& v2) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,scvector(v2));
  dot += tmp;
}

}  //namespace cxsc;

#include "sparsematrix.inl"

#endif 
 
