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

/* CVS $Id: simatrix.hpp,v 1.12 2010/12/21 09:54:02 cxsc Exp $ */

#ifndef _CXSC_SIMATRIX_HPP_INCLUDED
#define _CXSC_SIMATRIX_HPP_INCLUDED

#include <interval.hpp>
#include <imatrix.hpp>
#include <ivector.hpp>
#include <sivector.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cidot.hpp>
#include <sparseidot.hpp>
#include <sparsematrix.hpp>
#include <srmatrix.hpp>

namespace cxsc {

//definiert in srmatrix.hpp
//enum STORAGE_TYPE{triplet,compressed_row,compressed_column};

class simatrix_slice;
class simatrix_subv;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;

inline bool comp_pair_i(std::pair<int,interval> p1, std::pair<int,interval> p2) {
  return p1.first < p2.first;
}

class simatrix {

  private:
    std::vector<int> p;
    std::vector<int> ind;
    std::vector<interval> x;
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

    std::vector<interval>& values() {
      return x;
    }

    const std::vector<int>& column_pointers() const {
      return p;
    }

    const std::vector<int>& row_indices() const {
      return ind;
    }

    const std::vector<interval>& values() const {
      return x;
    }

    simatrix() {
      p.push_back(0);
      m = n = 0;
      lb1 = lb2 = ub1 = ub2 = 0;
    }

    simatrix(const int r, const int c) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(2*(m+n));
      x.reserve(2*(m+n));

      p[0] = 0;
    }

    simatrix(const int r, const int c, const int e) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(e);
      x.reserve(e);

      p[0] = 0;
    }

    simatrix(const int m, const int n, const int nnz, const intvector& rows, const intvector& cols, const ivector& values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<interval>(rows[Lb(rows)+k],cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<interval>(j,cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<std::pair<int,interval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[Lb(cols)+k],values[Lb(values)+k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_i);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    simatrix(const int m, const int n, const int nnz, const int* rows, const int* cols, const interval* values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<interval>(rows[k],cols[k],values[k]));
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

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<interval>(j,cols[k],values[k]));
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

         std::vector<std::pair<int,interval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[k],values[k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_i);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }


    simatrix(const srmatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(interval(A.x[i]));
    }


    simatrix(const rmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(interval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    simatrix(const imatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(interval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    simatrix(const srmatrix_slice&);
    simatrix(const simatrix_slice&);

    void full(imatrix& A) const {
       A = imatrix(lb1,ub1,lb2,ub2);
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
      std::vector<interval> xnew;
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


    simatrix& operator=(const real& A) {
      return sp_ms_assign<simatrix,real,interval>(*this,A);
    }

    simatrix& operator=(const interval& A) {
      return sp_ms_assign<simatrix,interval,interval>(*this,A);
    }

    simatrix& operator=(const rmatrix& A) {
      return spf_mm_assign<simatrix,rmatrix,interval>(*this,A);
    }

    simatrix& operator=(const imatrix& A) {
      return spf_mm_assign<simatrix,imatrix,interval>(*this,A);
    }

    simatrix& operator=(const rmatrix_slice& A) {
      return spf_mm_assign<simatrix,rmatrix_slice,interval>(*this,A);
    }

    simatrix& operator=(const imatrix_slice& A) {
      return spf_mm_assign<simatrix,imatrix_slice,interval>(*this,A);
    }

    simatrix& operator=(const srmatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(interval(A.x[i]));
      return *this;
    }

    /* simatrix& operator=(const simatrix& A) {
      p = A.p;
      ind = A.ind;
      x = A.x;
      return *this;
    } */

    simatrix& operator=(const srmatrix_slice&);
    simatrix& operator=(const simatrix_slice&);

    const interval operator()(int i, int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int)"));
#endif
      interval r(0.0);
      for(int k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  r = x[k];
      }
      return r;
    }

    interval& element(int i, int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int)"));
#endif
      int k;
      for(k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  return x[k];
      }

      //Nicht gefunden, Element muss angelegt werden, da Schreibzugriff moeglich
      std::vector<int>::iterator ind_it = ind.begin() + k;
      std::vector<interval>::iterator x_it  = x.begin() + k;
      ind.insert(ind_it, i-lb1);
      x_it = x.insert(x_it, interval(0.0));
      for(k=j-lb2+1 ; k<(int)p.size() ; k++)
        p[k]++;

      return *x_it;
    }

    simatrix_subv operator[](const cxscmatrix_column&);
    simatrix_subv operator[](const int);
    const simatrix_subv operator[](const cxscmatrix_column&) const;
    const simatrix_subv operator[](const int) const;

    simatrix_slice operator()(const int, const int , const int, const int);

    simatrix operator()(const intvector& pervec, const intvector& q) {
      simatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      int nnz=0;
      for(int k=0 ; k<n ; k++) {
        A.p[k] = nnz;

        std::map<int,interval> work;
        for(int j=p[q[Lb(q)+k]] ; j<p[q[Lb(q)+k]+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,interval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }

        nnz += work.size();
 
      }

      A.p[n] = nnz;

      return A;
    }

    simatrix operator()(const intvector& pervec) {
      simatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      for(int k=0 ; k<n ; k++) {
        A.p[k] = p[k];

        std::map<int,interval> work;
        for(int j=p[k] ; j<p[k+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,interval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }
 
      }

      A.p[n] = p[n];

      return A;
    }

    simatrix operator()(const intmatrix& P, const intmatrix& Q) {
      intvector p = permvec(P);
      intvector q = perminv(permvec(Q));
      return (*this)(p,q);
    }

    simatrix operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    real density() const {
      return p[n]/((double)m*n);
    }

    int get_nnz() const {
      return p[n];
    }

    simatrix& operator+=(const rmatrix& B) {
      return spf_mm_addassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    simatrix& operator+=(const imatrix& B) {
      return spf_mm_addassign<simatrix,imatrix,imatrix>(*this,B);
    }

    simatrix& operator+=(const rmatrix_slice& B) {
      return spf_mm_addassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator+=(const imatrix_slice& B) {
      return spf_mm_addassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator+=(const srmatrix& B) {
      return spsp_mm_addassign<simatrix,srmatrix,interval>(*this,B);
    }

    simatrix& operator+=(const simatrix& B) {
      return spsp_mm_addassign<simatrix,simatrix,interval>(*this,B);
    }

    simatrix& operator-=(const rmatrix& B) {
      return spf_mm_subassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    simatrix& operator-=(const imatrix& B) {
      return spf_mm_subassign<simatrix,imatrix,imatrix>(*this,B);
    }

    simatrix& operator-=(const rmatrix_slice& B) {
      return spf_mm_subassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator-=(const imatrix_slice& B) {
      return spf_mm_subassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator-=(const srmatrix& B) {
      return spsp_mm_subassign<simatrix,srmatrix,interval>(*this,B);
    }

    simatrix& operator-=(const simatrix& B) {
      return spsp_mm_subassign<simatrix,simatrix,interval>(*this,B);
    }

    simatrix& operator*=(const imatrix& B) {
      return spf_mm_multassign<simatrix,imatrix,sparse_idot,imatrix>(*this,B);
    }

    simatrix& operator*=(const rmatrix& B) {
      return spf_mm_multassign<simatrix,rmatrix,sparse_idot,imatrix>(*this,B);
    }

    simatrix& operator*=(const rmatrix_slice& B) {
      return spf_mm_multassign<simatrix,rmatrix_slice,sparse_idot,imatrix>(*this,B);
    }

    simatrix& operator*=(const imatrix_slice& B) {
      return spf_mm_multassign<simatrix,imatrix_slice,sparse_idot,imatrix>(*this,B);
    }

    simatrix& operator*=(const srmatrix& B) {
      return spsp_mm_multassign<simatrix,srmatrix,sparse_idot,interval>(*this,B);
    }

    simatrix& operator*=(const simatrix& B) {
      return spsp_mm_multassign<simatrix,simatrix,sparse_idot,interval>(*this,B);
    }

    simatrix& operator*=(const real& r) {
      return sp_ms_multassign(*this,r);
    }

    simatrix& operator*=(const interval& r) {
      return sp_ms_multassign(*this,r);
    }

    simatrix& operator/=(const real& r) {
      return sp_ms_divassign(*this,r);
    }

    simatrix& operator/=(const interval& r) {
      return sp_ms_divassign(*this,r);
    }

    simatrix& operator|=(const rmatrix& B) {
      return spf_mm_hullassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    simatrix& operator|=(const imatrix& B) {
      return spf_mm_hullassign<simatrix,imatrix,imatrix>(*this,B);
    }

    simatrix& operator|=(const rmatrix_slice& B) {
      return spf_mm_hullassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator|=(const imatrix_slice& B) {
      return spf_mm_hullassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator|=(const srmatrix& B) {
      return spsp_mm_hullassign<simatrix,srmatrix,interval>(*this,B);
    }

    simatrix& operator|=(const simatrix& B) {
      return spsp_mm_hullassign<simatrix,simatrix,interval>(*this,B);
    }

    simatrix& operator&=(const imatrix& B) {
      return spf_mm_intersectassign<simatrix,imatrix,imatrix>(*this,B);
    }

    simatrix& operator&=(const imatrix_slice& B) {
      return spf_mm_intersectassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    simatrix& operator&=(const simatrix& B) {
      return spsp_mm_intersectassign<simatrix,simatrix,interval>(*this,B);
    }

    friend void SetLb(simatrix&, const int, const int);
    friend void SetUb(simatrix&, const int, const int);    
    friend int Lb(const simatrix&, int);
    friend int Ub(const simatrix&, int);
    friend int RowLen(const simatrix&);
    friend int ColLen(const simatrix&);
    friend srmatrix Inf(const simatrix&);
    friend srmatrix Sup(const simatrix&);
    friend simatrix Re(const scimatrix&);
    friend simatrix Im(const scimatrix&);
    friend simatrix abs(const simatrix&);
    friend srmatrix mid(const simatrix&);
    friend srmatrix diam(const simatrix&);
    friend simatrix abs(const scimatrix&);
    friend srmatrix absmin(const simatrix&);
    friend srmatrix absmax(const simatrix&);

    friend simatrix transp(const simatrix&);
    friend simatrix Id(const simatrix&);
    friend srmatrix CompMat(const simatrix&);

    friend std::istream& operator>>(std::istream&, simatrix_slice&);
    friend std::istream& operator>>(std::istream&, simatrix_subv&);

    friend class srmatrix_slice;
    friend class srmatrix_subv;
    friend class srvector;
    friend class simatrix_slice;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class scivector;
    friend class rmatrix;
    friend class imatrix;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"
};

inline imatrix::imatrix(const srmatrix& A) {
  dat = new interval[A.m*A.n];
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

inline imatrix::imatrix(const simatrix& A) {
  dat = new interval[A.m*A.n];
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

inline simatrix Id(const simatrix& A) {
  simatrix I(A.m, A.n, (A.m>A.n) ? A.m : A.n);
  I.lb1 = A.lb1; I.lb2 = A.lb2;
  I.ub1 = A.ub1; I.ub2 = A.ub2;

  if(A.m < A.n) {
    for(int i=0 ; i<A.m ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(interval(1.0));
    }
  } else {
    for(int i=0 ; i<A.n ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(interval(1.0));
    }
  }

  return I;
}

inline simatrix transp(const simatrix& A) {
  simatrix B(A.n, A.m, A.get_nnz());
     
  //NIchtnullen pro Zeile bestimmen
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

inline void SetLb(simatrix& A, const int i, const int j) {
  if(i==1) {
    A.lb1 = j;
    A.ub1 = j + A.m - 1;
  } else if(i==2) {
    A.lb2 = j;
    A.ub2 = j + A.n - 1;
  }
}

inline void SetUb(simatrix& A, const int i, const int j) {
  if(i==1) {
    A.ub1 = j;
    A.lb1 = j - A.m + 1;
  } else if(i==2) {
    A.ub2 = j;
    A.lb2 = j - A.n + 1;
  }
}

inline int Lb(const simatrix& A, int i) {
  if(i==1) 
    return A.lb1;
  else if(i==2)
    return A.lb2;
  else
    return 1;
}

inline int Ub(const simatrix& A, int i) {
  if(i==1) 
    return A.ub1;
  else if(i==2)
    return A.ub2;
  else
    return 1;
}

inline int RowLen(const simatrix& A) {
  return A.n;
}

inline int ColLen(const simatrix& A) {
  return A.m;
}

inline void Resize(simatrix& A) {
  sp_m_resize(A);
}

inline void Resize(simatrix& A, const int m, const int n) {
  sp_m_resize(A,m,n);
}

inline void Resize(simatrix& A, const int l1, const int u1, const int l2, const int u2) {
  sp_m_resize(A,l1,u1,l2,u2);
}

inline srmatrix Inf(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Inf(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix Sup(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Sup(A.x[i]));

  res.dropzeros();

  return res; 
}

inline simatrix abs(const simatrix& A) {
  simatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(abs(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix absmin(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(AbsMin(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix absmax(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(AbsMax(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix CompMat(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;
  res.p[A.n] = A.p[A.n];

  for(int j=0 ; j<res.n ; j++) {
    res.p[j] = A.p[j];
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      res.ind.push_back(A.ind[k]);
      if(A.ind[k] == j)
        res.x.push_back(AbsMin(A.x[k]));
      else
        res.x.push_back(-AbsMax(A.x[k]));
    }
  }

  res.dropzeros();

  return res; 
}

inline srmatrix mid(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(mid(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix diam(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(diam(A.x[i]));

  res.dropzeros();

  return res; 
}


inline imatrix operator*(const imatrix& A, const srmatrix& B) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const rmatrix& A, const simatrix& B) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const imatrix& A, const simatrix& B) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const simatrix& A, const rmatrix& B) {
  return spf_mm_mult<simatrix,rmatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const srmatrix& A, const imatrix& B) {
  return spf_mm_mult<srmatrix,imatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const simatrix& A, const imatrix& B) {
  return spf_mm_mult<simatrix,imatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_mult<imatrix_slice,srmatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<rmatrix_slice,simatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<imatrix_slice,simatrix,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_mult<simatrix,rmatrix_slice,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<srmatrix,imatrix_slice,imatrix,sparse_idot>(A,B);
}

inline imatrix operator*(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<simatrix,imatrix_slice,imatrix,sparse_idot>(A,B);
}

inline simatrix operator*(const simatrix& A, const srmatrix& B) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(A,B);
}

inline simatrix operator*(const srmatrix& A, const simatrix& B) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(A,B);
}

inline simatrix operator*(const simatrix& A, const simatrix& B) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(A,B);
}

inline simatrix operator/(const simatrix& A, const real& r) {
  return sp_ms_div<simatrix,real,simatrix>(A,r);
}

inline simatrix operator/(const simatrix& A, const interval& r) {
  return sp_ms_div<simatrix,interval,simatrix>(A,r);
}

inline simatrix operator/(const srmatrix& A, const interval& r) {
  return sp_ms_div<srmatrix,interval,simatrix>(A,r);
}

inline simatrix operator*(const simatrix& A, const real& r) {
  return sp_ms_mult<simatrix,real,simatrix>(A,r);
}

inline simatrix operator*(const simatrix& A, const interval& r) {
  return sp_ms_mult<simatrix,interval,simatrix>(A,r);
}

inline simatrix operator*(const srmatrix& A, const interval& r) {
  return sp_ms_mult<srmatrix,interval,simatrix>(A,r);
}

inline simatrix operator*(const real& r, const simatrix& A) {
  return sp_sm_mult<real,simatrix,simatrix>(r,A);
}

inline simatrix operator*(const interval& r, const simatrix& A) {
  return sp_sm_mult<interval,simatrix,simatrix>(r,A);
}

inline simatrix operator*(const interval& r, const srmatrix& A) {
  return sp_sm_mult<interval,srmatrix,simatrix>(r,A);
}

inline ivector operator*(const simatrix& A, const rvector& v) {
  return spf_mv_mult<simatrix,rvector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const srmatrix& A, const ivector& v) {
  return spf_mv_mult<srmatrix,ivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const simatrix& A, const ivector& v) {
  return spf_mv_mult<simatrix,ivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const simatrix& A, const rvector_slice& v) {
  return spf_mv_mult<simatrix,rvector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const srmatrix& A, const ivector_slice& v) {
  return spf_mv_mult<srmatrix,ivector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const simatrix& A, const ivector_slice& v) {
  return spf_mv_mult<simatrix,ivector_slice,ivector,sparse_idot>(A,v);
}

inline sivector operator*(const simatrix& A, const srvector& v) {
  return spsp_mv_mult<simatrix,srvector,sivector,sparse_idot,interval>(A,v);
}

inline sivector operator*(const srmatrix& A, const sivector& v) {
  return spsp_mv_mult<srmatrix,sivector,sivector,sparse_idot,interval>(A,v);
}

inline sivector operator*(const simatrix& A, const sivector& v) {
  return spsp_mv_mult<simatrix,sivector,sivector,sparse_idot,interval>(A,v);
}

inline sivector operator*(const simatrix& A, const srvector_slice& v) {
  return spsl_mv_mult<simatrix,srvector_slice,sivector,sparse_idot,interval>(A,v);
}

inline sivector operator*(const srmatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<srmatrix,sivector_slice,sivector,sparse_idot,interval>(A,v);
}

inline sivector operator*(const simatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<simatrix,sivector_slice,sivector,sparse_idot,interval>(A,v);
}

inline ivector operator*(const imatrix& A, const srvector& v) {
  return fsp_mv_mult<imatrix,srvector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const rmatrix& A, const sivector& v) {
  return fsp_mv_mult<rmatrix,sivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix& A, const sivector& v) {
  return fsp_mv_mult<imatrix,sivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix_slice& A, const srvector& v) {
  return fsp_mv_mult<imatrix_slice,srvector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const rmatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<rmatrix_slice,sivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<imatrix_slice,sivector,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix& A, const srvector_slice& v) {
  return fsl_mv_mult<imatrix,srvector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const rmatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<rmatrix,sivector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<imatrix,sivector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix_slice& A, const srvector_slice& v) {
  return fsl_mv_mult<imatrix_slice,srvector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const rmatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<rmatrix_slice,sivector_slice,ivector,sparse_idot>(A,v);
}

inline ivector operator*(const imatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<imatrix_slice,sivector_slice,ivector,sparse_idot>(A,v);
}

inline imatrix operator+(const imatrix& A, const srmatrix& B) {
  return fsp_mm_add<imatrix,srmatrix,imatrix>(A,B);
}

inline imatrix operator+(const rmatrix& A, const simatrix& B) {
  return fsp_mm_add<rmatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator+(const imatrix& A, const simatrix& B) {
  return fsp_mm_add<imatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator+(const simatrix& A, const rmatrix& B) {
  return spf_mm_add<simatrix,rmatrix,imatrix>(A,B);
}

inline imatrix operator+(const srmatrix& A, const imatrix& B) {
  return spf_mm_add<srmatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator+(const simatrix& A, const imatrix& B) {
  return spf_mm_add<simatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator+(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_add<imatrix_slice,srmatrix,imatrix>(A,B);
}

inline imatrix operator+(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<rmatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator+(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<imatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator+(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_add<simatrix,rmatrix_slice,imatrix>(A,B);
}

inline imatrix operator+(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_add<srmatrix,imatrix_slice,imatrix>(A,B);
}

inline imatrix operator+(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_add<simatrix,imatrix_slice,imatrix>(A,B);
}

inline simatrix operator+(const simatrix& A, const srmatrix& B) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(A,B);
}

inline simatrix operator+(const srmatrix& A, const simatrix& B) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(A,B);
}

inline simatrix operator+(const simatrix& A, const simatrix& B) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(A,B);
}

inline imatrix operator-(const imatrix& A, const srmatrix& B) {
  return fsp_mm_sub<imatrix,srmatrix,imatrix>(A,B);
}

inline imatrix operator-(const rmatrix& A, const simatrix& B) {
  return fsp_mm_sub<rmatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator-(const imatrix& A, const simatrix& B) {
  return fsp_mm_sub<imatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator-(const simatrix& A, const rmatrix& B) {
  return spf_mm_sub<simatrix,rmatrix,imatrix>(A,B);
}

inline imatrix operator-(const srmatrix& A, const imatrix& B) {
  return spf_mm_sub<srmatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator-(const simatrix& A, const imatrix& B) {
  return spf_mm_sub<simatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator-(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_sub<imatrix_slice,srmatrix,imatrix>(A,B);
}

inline imatrix operator-(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<rmatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator-(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<imatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator-(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_sub<simatrix,rmatrix_slice,imatrix>(A,B);
}

inline imatrix operator-(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<srmatrix,imatrix_slice,imatrix>(A,B);
}

inline imatrix operator-(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<simatrix,imatrix_slice,imatrix>(A,B);
}

inline simatrix operator-(const simatrix& A, const srmatrix& B) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(A,B);
}

inline simatrix operator-(const srmatrix& A, const simatrix& B) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(A,B);
}

inline simatrix operator-(const simatrix& A, const simatrix& B) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(A,B);
}

inline imatrix operator|(const imatrix& A, const srmatrix& B) {
  return fsp_mm_hull<imatrix,srmatrix,imatrix>(A,B);
}

inline imatrix operator|(const rmatrix& A, const simatrix& B) {
  return fsp_mm_hull<rmatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator|(const imatrix& A, const simatrix& B) {
  return fsp_mm_hull<imatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator|(const simatrix& A, const rmatrix& B) {
  return spf_mm_hull<simatrix,rmatrix,imatrix>(A,B);
}

inline imatrix operator|(const srmatrix& A, const imatrix& B) {
  return spf_mm_hull<srmatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator|(const simatrix& A, const imatrix& B) {
  return spf_mm_hull<simatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator|(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<imatrix_slice,srmatrix,imatrix>(A,B);
}

inline imatrix operator|(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<rmatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator|(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<imatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator|(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<simatrix,rmatrix_slice,imatrix>(A,B);
}

inline imatrix operator|(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<srmatrix,imatrix_slice,imatrix>(A,B);
}

inline imatrix operator|(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<simatrix,imatrix_slice,imatrix>(A,B);
}

inline simatrix operator|(const simatrix& A, const srmatrix& B) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(A,B);
}

inline simatrix operator|(const srmatrix& A, const simatrix& B) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(A,B);
}

inline simatrix operator|(const simatrix& A, const simatrix& B) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(A,B);
}

inline imatrix operator|(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_hull<rmatrix,srmatrix,imatrix>(A,B);
}

inline imatrix operator|(const srmatrix& A, const rmatrix& B) {
  return spf_mm_hull<srmatrix,rmatrix,imatrix>(A,B);
}

inline imatrix operator|(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<rmatrix_slice,srmatrix,imatrix>(A,B);
}

inline imatrix operator|(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<srmatrix,rmatrix_slice,imatrix>(A,B);
}

inline simatrix operator|(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(A,B);
}

inline imatrix operator&(const imatrix& A, const simatrix& B) {
  return fsp_mm_intersect<imatrix,simatrix,imatrix>(A,B);
}

inline imatrix operator&(const simatrix& A, const imatrix& B) {
  return spf_mm_intersect<simatrix,imatrix,imatrix>(A,B);
}

inline imatrix operator&(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_intersect<imatrix_slice,simatrix,imatrix>(A,B);
}

inline imatrix operator&(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_intersect<simatrix,imatrix_slice,imatrix>(A,B);
}

inline simatrix operator&(const simatrix& A, const simatrix& B) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(A,B);
}

inline simatrix operator-(const simatrix& M) {
  return sp_m_negative<simatrix,simatrix>(M);
}

inline simatrix& operator+(simatrix& A) {
  return A;
}

inline imatrix& imatrix::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline imatrix& imatrix::operator=(const simatrix& B) {
  *this = imatrix(B);
  return *this;
}

inline imatrix& imatrix::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix& imatrix::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix& imatrix::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix& imatrix::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix& imatrix::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<imatrix,srmatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix& imatrix::operator*=(const simatrix& B) {
  return fsp_mm_multassign<imatrix,simatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix_slice& imatrix_slice::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<imatrix_slice,srmatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix_slice& imatrix_slice::operator*=(const simatrix& B) {
  return fsp_mm_multassign<imatrix_slice,simatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix& imatrix::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix& imatrix::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix& imatrix::operator&=(const srmatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix& imatrix::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator&=(const srmatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline bool operator==(const simatrix& A, const srmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const simatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const simatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const rmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const imatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const imatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const imatrix& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const imatrix& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const srmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const simatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const simatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const rmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const imatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const imatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const imatrix& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const imatrix& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const imatrix_slice& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix_slice& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const imatrix_slice& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const rmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const imatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const imatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator<(const srmatrix& A, const simatrix& B) {
  return spsp_mm_less<srmatrix,simatrix,interval>(A,B);
}

inline bool operator<(const simatrix& A, const simatrix& B) {
  return spsp_mm_less<simatrix,simatrix,interval>(A,B);
}

inline bool operator<(const srmatrix& A, const imatrix& B) {
  return spf_mm_less<srmatrix,imatrix,interval>(A,B);
}

inline bool operator<(const simatrix& A, const imatrix& B) {
  return spf_mm_less<simatrix,imatrix,interval>(A,B);
}

inline bool operator<(const rmatrix& A, const simatrix& B) {
  return fsp_mm_less<rmatrix,simatrix,interval>(A,B);
}

inline bool operator<(const imatrix& A, const simatrix& B) {
  return fsp_mm_less<imatrix,simatrix,interval>(A,B);
}

inline bool operator<(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_less<rmatrix_slice,simatrix,interval>(A,B);
}

inline bool operator<(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_less<imatrix_slice,simatrix,interval>(A,B);
}

inline bool operator<(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_less<srmatrix,imatrix_slice,interval>(A,B);
}

inline bool operator<(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_less<simatrix,imatrix_slice,interval>(A,B);
}

inline bool operator<=(const srmatrix& A, const simatrix& B) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(A,B);
}

inline bool operator<=(const simatrix& A, const simatrix& B) {
  return spsp_mm_leq<simatrix,simatrix,interval>(A,B);
}

inline bool operator<=(const srmatrix& A, const imatrix& B) {
  return spf_mm_leq<srmatrix,imatrix,interval>(A,B);
}

inline bool operator<=(const simatrix& A, const imatrix& B) {
  return spf_mm_leq<simatrix,imatrix,interval>(A,B);
}

inline bool operator<=(const rmatrix& A, const simatrix& B) {
  return fsp_mm_leq<rmatrix,simatrix,interval>(A,B);
}

inline bool operator<=(const imatrix& A, const simatrix& B) {
  return fsp_mm_leq<imatrix,simatrix,interval>(A,B);
}

inline bool operator<=(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_leq<rmatrix_slice,simatrix,interval>(A,B);
}

inline bool operator<=(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_leq<imatrix_slice,simatrix,interval>(A,B);
}

inline bool operator<=(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_leq<srmatrix,imatrix_slice,interval>(A,B);
}

inline bool operator<=(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_leq<simatrix,imatrix_slice,interval>(A,B);
}

inline bool operator>(const simatrix& A, const srmatrix& B) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(A,B);
}

inline bool operator>(const simatrix& A, const simatrix& B) {
  return spsp_mm_greater<simatrix,simatrix,interval>(A,B);
}

inline bool operator>(const simatrix& A, const rmatrix& B) {
  return spf_mm_greater<simatrix,rmatrix,interval>(A,B);
}

inline bool operator>(const simatrix& A, const imatrix& B) {
  return spf_mm_greater<simatrix,imatrix,interval>(A,B);
}

inline bool operator>(const imatrix& A, const srmatrix& B) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(A,B);
}

inline bool operator>(const imatrix& A, const simatrix& B) {
  return fsp_mm_greater<imatrix,simatrix,interval>(A,B);
}

inline bool operator>(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_greater<imatrix_slice,srmatrix,interval>(A,B);
}

inline bool operator>(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_greater<imatrix_slice,simatrix,interval>(A,B);
}

inline bool operator>(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_greater<simatrix,rmatrix_slice,interval>(A,B);
}

inline bool operator>(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_greater<simatrix,imatrix_slice,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const srmatrix& B) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const simatrix& B) {
  return spsp_mm_geq<simatrix,simatrix,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const rmatrix& B) {
  return spf_mm_geq<simatrix,rmatrix,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const imatrix& B) {
  return spf_mm_geq<simatrix,imatrix,interval>(A,B);
}

inline bool operator>=(const imatrix& A, const srmatrix& B) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(A,B);
}

inline bool operator>=(const imatrix& A, const simatrix& B) {
  return fsp_mm_geq<imatrix,simatrix,interval>(A,B);
}

inline bool operator>=(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_geq<imatrix_slice,srmatrix,interval>(A,B);
}

inline bool operator>=(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_geq<imatrix_slice,simatrix,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_geq<simatrix,rmatrix_slice,interval>(A,B);
}

inline bool operator>=(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_geq<simatrix,imatrix_slice,interval>(A,B);
}

inline bool operator!(const simatrix& A) {
  return sp_m_not(A);
}

inline std::ostream& operator<<(std::ostream& os, const simatrix& A) {
  return sp_m_output<simatrix,interval>(os,A);
}

inline std::istream& operator>>(std::istream& is, simatrix& A) {
  return sp_m_input<simatrix,interval>(is,A);
}

class simatrix_slice {
  public:
    simatrix  A;
    simatrix* M; //Originalmatrix

  private:
    simatrix_slice(simatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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

    simatrix_slice(const simatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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
        M = const_cast<simatrix*>(&Mat); //Vorgehen noetig um Schreibweise A[i][j] bei auslesen von const A zu ermoeglichen
    }


  public:
    simatrix_slice& operator=(const real& C) {
      return sl_ms_assign<simatrix_slice, real, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const interval& C) {
      return sl_ms_assign<simatrix_slice, interval, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const srmatrix& C) {
      return slsp_mm_assign<simatrix_slice, srmatrix, std::vector<interval>::iterator>(*this,C);
    }

    simatrix_slice& operator=(const simatrix& C) {
      return slsp_mm_assign<simatrix_slice, simatrix, std::vector<interval>::iterator>(*this,C);
    }

    simatrix_slice& operator=(const rmatrix& C) {
      return slf_mm_assign<simatrix_slice, rmatrix, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const imatrix& C) {
      return slf_mm_assign<simatrix_slice, imatrix, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const rmatrix_slice& C) {
      return slf_mm_assign<simatrix_slice, rmatrix_slice, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const imatrix_slice& C) {
      return slf_mm_assign<simatrix_slice, imatrix_slice, std::vector<interval>::iterator, interval>(*this,C);
    }

    simatrix_slice& operator=(const srmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    simatrix_slice& operator=(const simatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    simatrix_slice& operator*=(const srmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    simatrix_slice& operator*=(const simatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    simatrix_slice& operator*=(const srmatrix& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const simatrix& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const rmatrix& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const imatrix& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const rmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const imatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    simatrix_slice& operator*=(const real& r) {
      *this = A*r;
      return *this;
    }

    simatrix_slice& operator*=(const interval& r) {
      *this = A*r;
      return *this;
    }

    simatrix_slice& operator/=(const real& r) {
      *this = A/r;
      return *this;
    }

    simatrix_slice& operator/=(const interval& r) {
      *this = A/r;
      return *this;
    }

    simatrix_slice& operator+=(const srmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    simatrix_slice& operator+=(const simatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    simatrix_slice& operator+=(const srmatrix& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator+=(const simatrix& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator+=(const rmatrix& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator+=(const imatrix& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator+=(const rmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator+=(const imatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    simatrix_slice& operator-=(const srmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    simatrix_slice& operator-=(const simatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    simatrix_slice& operator-=(const srmatrix& M) {
      *this = A-M;
      return *this;
    } 

    simatrix_slice& operator-=(const simatrix& M) {
      *this = A-M;
      return *this;
    } 

    simatrix_slice& operator-=(const rmatrix& M) {
      *this = A-M;
      return *this;
    } 

    simatrix_slice& operator-=(const imatrix& M) {
      *this = A-M;
      return *this;
    } 

    simatrix_slice& operator-=(const rmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    simatrix_slice& operator-=(const imatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    simatrix_slice& operator|=(const srmatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    simatrix_slice& operator|=(const simatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    simatrix_slice& operator|=(const srmatrix& M) {
      *this = A|M;
      return *this;
    } 

    simatrix_slice& operator|=(const simatrix& M) {
      *this = A|M;
      return *this;
    } 

    simatrix_slice& operator|=(const rmatrix& M) {
      *this = A|M;
      return *this;
    } 

    simatrix_slice& operator|=(const imatrix& M) {
      *this = A|M;
      return *this;
    } 

    simatrix_slice& operator|=(const rmatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    simatrix_slice& operator|=(const imatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    const interval operator()(const int i, const int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_slice::operator()(int, int)"));
#endif
      interval r = A(i,j);
      return r;
    }

    interval& element(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_slice::element(int, int)"));
#endif
      return M->element(i,j);
    }

    simatrix_subv operator[](const int);
    simatrix_subv operator[](const cxscmatrix_column&);
    const simatrix_subv operator[](const int) const;
    const simatrix_subv operator[](const cxscmatrix_column&) const;

    friend std::ostream& operator<<(std::ostream&, const simatrix_slice&);

    friend int Lb(const simatrix_slice&, const int);
    friend int Ub(const simatrix_slice&, const int);
    friend srmatrix Inf(const simatrix_slice&);
    friend srmatrix Sup(const simatrix_slice&);
    friend simatrix abs(const simatrix_slice&);
    friend srmatrix mid(const simatrix_slice&);
    friend srmatrix diam(const simatrix_slice&);
    friend int RowLen(const simatrix_slice&);
    friend int ColLen(const simatrix_slice&);

    friend class srmatrix;
    friend class srmatrix_subv;
    friend class srvector;
    friend class simatrix;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_subv;
    friend class scimatrix_slice;
    friend class scivector;
    friend class imatrix;
    friend class cimatrix;


#include "matrix_friend_declarations.inl"    
};

inline imatrix::imatrix(const srmatrix_slice& A) {
  dat = new interval[A.A.m*A.A.n];
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

inline imatrix::imatrix(const simatrix_slice& A) {
  dat = new interval[A.A.m*A.A.n];
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

inline int RowLen(const simatrix_slice& S) {
  return RowLen(S.A);
}

inline int ColLen(const simatrix_slice& S) {
  return ColLen(S.A);
}

inline simatrix_slice simatrix::operator()(const int i, const int j, const int k, const int l) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int)"));
#endif
  return simatrix_slice(*this, i, j, k, l);
}

inline int Lb(const simatrix_slice& S, const int i) {
  return Lb(S.A, i);
}

inline int Ub(const simatrix_slice& S, const int i) {
  return Ub(S.A, i);
}

inline srmatrix Inf(const simatrix_slice& S) {
  return Inf(S.A);
}

inline srmatrix Sup(const simatrix_slice& S) {
  return Sup(S.A);
}

inline simatrix abs(const simatrix_slice& S) {
  return abs(S.A);
}

inline srmatrix mid(const simatrix_slice& S) {
  return mid(S.A);
}

inline srmatrix diam(const simatrix_slice& S) {
  return diam(S.A);
}

inline simatrix::simatrix(const srmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline simatrix::simatrix(const simatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline simatrix& simatrix::operator=(const srmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline simatrix& simatrix::operator=(const simatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline imatrix& imatrix::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline imatrix& imatrix::operator=(const simatrix_slice& M) {
  *this = imatrix(M);
  return *this;
}

inline imatrix& imatrix::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix& imatrix::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix& imatrix::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix& imatrix::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix& imatrix::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix& imatrix::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix& imatrix::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix& imatrix::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix& imatrix::operator&=(const srmatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix& imatrix::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator&=(const srmatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline simatrix operator-(const simatrix_slice& M) {
  return sp_m_negative<simatrix,simatrix>(M.A);
}

inline simatrix operator+(const simatrix_slice& M) {
  return M.A;
}

inline simatrix operator*(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

inline simatrix operator*(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

inline simatrix operator*(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

inline simatrix operator*(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

inline simatrix operator*(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

inline simatrix operator*(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

inline simatrix operator*(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

inline simatrix operator*(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

inline simatrix operator*(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

inline imatrix operator*(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_mult<simatrix,rmatrix,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<srmatrix,imatrix,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<simatrix,imatrix,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline imatrix operator*(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline imatrix operator*(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline imatrix operator*(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_mult<simatrix,rmatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<srmatrix,imatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<simatrix,imatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

inline imatrix operator*(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline imatrix operator*(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline imatrix operator*(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

inline sivector operator*(const simatrix_slice& M, const srvector& v) {
  return spsp_mv_mult<simatrix,srvector,sivector,sparse_idot,interval>(M.A,v);
}

inline sivector operator*(const srmatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<srmatrix,sivector,sivector,sparse_idot,interval>(M.A,v);
}

inline sivector operator*(const simatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<simatrix,sivector,sivector,sparse_idot,interval>(M.A,v);
}

inline sivector operator*(const simatrix_slice& M, const srvector_slice& v) {
  return spsl_mv_mult<simatrix,srvector_slice,sivector,sparse_idot,interval>(M.A,v);
}

inline sivector operator*(const srmatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<srmatrix,sivector_slice,sivector,sparse_idot,interval>(M.A,v);
}

inline sivector operator*(const simatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<simatrix,sivector_slice,sivector,sparse_idot,interval>(M.A,v);
}

inline ivector operator*(const simatrix_slice& M, const rvector& v) {
  return spf_mv_mult<simatrix,rvector,ivector,sparse_idot>(M.A,v);
}

inline ivector operator*(const srmatrix_slice& M, const ivector& v) {
  return spf_mv_mult<srmatrix,ivector,ivector,sparse_idot>(M.A,v);
}

inline ivector operator*(const simatrix_slice& M, const ivector& v) {
  return spf_mv_mult<simatrix,ivector,ivector,sparse_idot>(M.A,v);
}

inline ivector operator*(const simatrix_slice& M, const rvector_slice& v) {
  return spf_mv_mult<simatrix,rvector_slice,ivector,sparse_idot>(M.A,v);
}

inline ivector operator*(const srmatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<srmatrix,ivector_slice,ivector,sparse_idot>(M.A,v);
}

inline ivector operator*(const simatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<simatrix,ivector_slice,ivector,sparse_idot>(M.A,v);
}

inline simatrix operator/(const simatrix_slice& M, const real& r) {
  return sp_ms_div<simatrix,real,simatrix>(M.A,r);
}

inline simatrix operator/(const simatrix_slice& M, const interval& r) {
  return sp_ms_div<simatrix,interval,simatrix>(M.A,r);
}

inline simatrix operator/(const srmatrix_slice& M, const interval& r) {
  return sp_ms_div<srmatrix,interval,simatrix>(M.A,r);
}

inline simatrix operator*(const simatrix_slice& M, const real& r) {
  return sp_ms_mult<simatrix,real,simatrix>(M.A,r);
}

inline simatrix operator*(const simatrix_slice& M, const interval& r) {
  return sp_ms_mult<simatrix,interval,simatrix>(M.A,r);
}

inline simatrix operator*(const srmatrix_slice& M, const interval& r) {
  return sp_ms_mult<srmatrix,interval,simatrix>(M.A,r);
}

inline simatrix operator*(const real& r, const simatrix_slice& M) {
  return sp_sm_mult<real,simatrix,simatrix>(r,M.A);
}

inline simatrix operator*(const interval& r, const srmatrix_slice& M) {
  return sp_sm_mult<interval,srmatrix,simatrix>(r,M.A);
}

inline simatrix operator*(const interval& r, const simatrix_slice& M) {
  return sp_sm_mult<interval,simatrix,simatrix>(r,M.A);
}

inline simatrix operator+(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator+(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator+(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator+(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator+(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator+(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator+(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator+(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator+(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline imatrix operator+(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_add<simatrix,rmatrix,imatrix>(M1.A,M2);
}

inline imatrix operator+(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<srmatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator+(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<simatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator+(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<imatrix,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator+(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<rmatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator+(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<imatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator+(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_add<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator+(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator+(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator+(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator+(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator+(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline simatrix operator-(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator-(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator-(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator-(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator-(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator-(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator-(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator-(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator-(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline imatrix operator-(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_sub<simatrix,rmatrix,imatrix>(M1.A,M2);
}

inline imatrix operator-(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<srmatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator-(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<simatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator-(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<imatrix,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator-(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<rmatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator-(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<imatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator-(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_sub<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator-(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator-(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator-(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator-(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator-(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline simatrix operator|(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator|(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator|(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator|(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator|(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator|(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator|(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator|(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline simatrix operator|(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline imatrix operator|(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<simatrix,rmatrix,imatrix>(M1.A,M2);
}

inline imatrix operator|(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<srmatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator|(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<simatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator|(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<imatrix,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<imatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator|(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator|(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator|(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline simatrix operator|(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator|(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator|(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

inline imatrix operator|(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<srmatrix,rmatrix,imatrix>(M1.A,M2);
}

inline imatrix operator|(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,srmatrix,imatrix>(M1,M2.A);
}

inline imatrix operator|(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<srmatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator|(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

inline simatrix operator&(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

inline simatrix operator&(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

inline simatrix operator&(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

inline imatrix operator&(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_intersect<simatrix,imatrix,imatrix>(M1.A,M2);
}

inline imatrix operator&(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<imatrix,simatrix,imatrix>(M1,M2.A);
}

inline imatrix operator&(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_intersect<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

inline imatrix operator&(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const simatrix_slice& M1, const srmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const simatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const simatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const srmatrix& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const simatrix& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const simatrix_slice& M1, const rmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const imatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const imatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const imatrix& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const imatrix& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator<(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator<(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator<(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_less<srmatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_less<simatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator<(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_less<rmatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_less<imatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_less<rmatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_less<imatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_less<srmatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_less<simatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator<=(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator<=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator<=(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_leq<srmatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_leq<simatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator<=(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<rmatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<imatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<rmatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<imatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_leq<srmatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_leq<simatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1.A,M2.A);
}

inline bool operator>(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator>(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator>(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator>(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_greater<simatrix,rmatrix,interval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_greater<simatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator>(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<imatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator>(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<imatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator>(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_greater<simatrix,rmatrix_slice,interval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_greater<simatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1.A,M2.A);
}

inline bool operator>=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1.A,M2.A);
}

inline bool operator>=(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1.A,M2);
}

inline bool operator>=(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>=(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator>=(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_geq<simatrix,rmatrix,interval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_geq<simatrix,imatrix,interval>(M1.A,M2);
}

inline bool operator>=(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>=(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<imatrix,simatrix,interval>(M1,M2.A);
}

inline bool operator>=(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(M1,M2.A);
}

inline bool operator>=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<imatrix_slice,simatrix,interval>(M1,M2.A);
}

inline bool operator>=(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_geq<simatrix,rmatrix_slice,interval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_geq<simatrix,imatrix_slice,interval>(M1.A,M2);
}

inline bool operator!(const simatrix_slice& M) {
  return sp_m_not(M.A);
}

inline std::ostream& operator<<(std::ostream& os, const simatrix_slice& M) {
  return sp_m_output<simatrix,interval>(os, M.A);
}

inline std::istream& operator>>(std::istream& is, simatrix_slice& M) {
  simatrix tmp(M.A.m, M.A.n);
  sp_m_input<simatrix,interval>(is, tmp);
  M = tmp;
  return is;
}


class simatrix_subv {
  private:
    simatrix_slice dat;
    bool row;
    int index;

    simatrix_subv(simatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

    simatrix_subv(const simatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

  public:
    interval& operator[](const int i) {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat.element(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat.element(i,index);
      }
    }

    const interval operator[](const int i) const{
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat(i,index);
      }
    }

    simatrix_subv& operator=(const real& v) {
      return sv_vs_assign(*this,v);
    }

    simatrix_subv& operator=(const interval& v) {
      return sv_vs_assign(*this,v);
    }

    simatrix_subv& operator=(const srvector& v) {
      return svsp_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const sivector& v) {
      return svsp_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const srvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const sivector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const rvector& v) {
      return svf_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const ivector& v) {
      return svf_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const rvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const ivector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    simatrix_subv& operator=(const srmatrix_subv& v) {
      return svsp_vv_assign(*this,srvector(v));
    }

    simatrix_subv& operator=(const simatrix_subv& v) {
      return svsp_vv_assign(*this,sivector(v));
    }

    simatrix_subv& operator*=(const real&);
    simatrix_subv& operator*=(const interval&);
    simatrix_subv& operator/=(const real&);
    simatrix_subv& operator/=(const interval&);
    simatrix_subv& operator+=(const srvector&);
    simatrix_subv& operator+=(const srvector_slice&);
    simatrix_subv& operator+=(const rvector&);
    simatrix_subv& operator+=(const rvector_slice&);
    simatrix_subv& operator-=(const srvector&);
    simatrix_subv& operator-=(const srvector_slice&);
    simatrix_subv& operator-=(const rvector&);
    simatrix_subv& operator-=(const rvector_slice&);
    simatrix_subv& operator+=(const sivector&);
    simatrix_subv& operator+=(const sivector_slice&);
    simatrix_subv& operator+=(const ivector&);
    simatrix_subv& operator+=(const ivector_slice&);
    simatrix_subv& operator-=(const sivector&);
    simatrix_subv& operator-=(const sivector_slice&);
    simatrix_subv& operator-=(const ivector&);
    simatrix_subv& operator-=(const ivector_slice&);
    simatrix_subv& operator|=(const srvector&);
    simatrix_subv& operator|=(const srvector_slice&);
    simatrix_subv& operator|=(const rvector&);
    simatrix_subv& operator|=(const rvector_slice&);
    simatrix_subv& operator|=(const sivector&);
    simatrix_subv& operator|=(const sivector_slice&);
    simatrix_subv& operator|=(const ivector&);
    simatrix_subv& operator|=(const ivector_slice&);

    friend sivector operator-(const simatrix_subv&);

    friend std::istream& operator>>(std::istream&, simatrix_subv&);

    friend int Lb(const simatrix_subv&);
    friend int Ub(const simatrix_subv&);
    friend int VecLen(const simatrix_subv&);
    friend srvector Inf(const simatrix_subv&);
    friend srvector Sup(const simatrix_subv&);

    friend class srvector;
    friend class srmatrix;
    friend class srmatrix_slice;
    friend class sivector;
    friend class simatrix;
    friend class simatrix_slice;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;

#include "vector_friend_declarations.inl"
};

inline int Lb(const simatrix_subv& S) {
  if(S.row)
    return Lb(S.dat, 2);
  else
    return Lb(S.dat, 1);
}

inline int Ub(const simatrix_subv& S) {
  if(S.row)
    return Ub(S.dat, 2);
  else
    return Ub(S.dat, 1);
}

inline int VecLen(const simatrix_subv& S) {
  return Ub(S)-Lb(S)+1;
}

inline srvector Inf(const simatrix_subv& S) {
  return Inf(sivector(S));
}

inline srvector Sup(const simatrix_subv& S) {
  return Sup(sivector(S));
}

inline srvector mid(const simatrix_subv& S) {
  return mid(sivector(S));
}

inline srvector diam(const simatrix_subv& S) {
  return diam(sivector(S));
}

inline sivector abs(const simatrix_subv& S) {
  return abs(sivector(S));
}

inline std::ostream& operator<<(std::ostream& os, const simatrix_subv& v) {
  os << sivector(v);
  return os;
}

inline std::istream& operator>>(std::istream& is, simatrix_subv& v) {
  int n = 0;
  if(v.row) n=v.dat.A.n; else n=v.dat.A.m;
  sivector tmp(n);
  is >> tmp;
  v = tmp;
  return is;
}

inline simatrix_subv simatrix::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline simatrix_subv simatrix::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int)"));
#endif
  return simatrix_subv(*this, true, i, i, lb2, ub2);
}

inline const simatrix_subv simatrix::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline const simatrix_subv simatrix::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int)"));
#endif
  return simatrix_subv(*this, true, i, i, lb2, ub2);
}

inline simatrix_subv simatrix_slice::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int"));
#endif
  return simatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline simatrix_subv simatrix_slice::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline const simatrix_subv simatrix_slice::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int"));
#endif
  return simatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline const simatrix_subv simatrix_slice::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline sivector::sivector(const simatrix_subv& A) {
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

inline sivector operator-(const simatrix_subv& v) {
 sivector s(v);
 return -s;
}

inline sivector operator/(const simatrix_subv& v1, const real& v2) {
  return sivector(v1) / v2;
}

inline sivector operator/(const simatrix_subv& v1, const interval& v2) {
  return sivector(v1) / v2;
}

inline sivector operator/(const srmatrix_subv& v1, const interval& v2) {
  return srvector(v1) / v2;
}

inline sivector operator*(const simatrix_subv& v1, const real& v2) {
  return sivector(v1) * v2;
}

inline sivector operator*(const simatrix_subv& v1, const interval& v2) {
  return sivector(v1) * v2;
}

inline sivector operator*(const srmatrix_subv& v1, const interval& v2) {
  return srvector(v1) * v2;
}

inline sivector operator*(const real& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline sivector operator*(const interval& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline sivector operator*(const interval& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline interval operator*(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) * v2;
}

inline interval operator*(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) * v2;
}

inline interval operator*(const sivector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline interval operator*(const srvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const sivector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline interval operator*(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const ivector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline interval operator*(const rvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const ivector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline interval operator*(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline interval operator*(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline sivector operator+(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) + v2;
}

inline sivector operator+(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) + v2;
}

inline sivector operator+(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) + v2;
}

inline sivector operator+(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) + v2;
}

inline sivector operator+(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) + v2;
}

inline sivector operator+(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) + v2;
}

inline ivector operator+(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) + v2;
}

inline ivector operator+(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) + v2;
}

inline ivector operator+(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) + v2;
}

inline ivector operator+(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) + v2;
}

inline ivector operator+(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) + v2;
}

inline ivector operator+(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) + v2;
}

inline sivector operator+(const sivector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline sivector operator+(const srvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline sivector operator+(const sivector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline sivector operator+(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline sivector operator+(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline sivector operator+(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline ivector operator+(const ivector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline ivector operator+(const rvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline ivector operator+(const ivector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline ivector operator+(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline ivector operator+(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline ivector operator+(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline sivector operator-(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) - v2;
}

inline sivector operator-(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) - v2;
}

inline sivector operator-(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) - v2;
}

inline sivector operator-(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) - v2;
}

inline sivector operator-(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) - v2;
}

inline sivector operator-(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) - v2;
}

inline ivector operator-(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) - v2;
}

inline ivector operator-(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) - v2;
}

inline ivector operator-(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) - v2;
}

inline ivector operator-(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) - v2;
}

inline ivector operator-(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) - v2;
}

inline ivector operator-(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) - v2;
}

inline sivector operator-(const sivector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline sivector operator-(const srvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline sivector operator-(const sivector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline sivector operator-(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline sivector operator-(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline sivector operator-(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline ivector operator-(const ivector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline ivector operator-(const rvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline ivector operator-(const ivector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline ivector operator-(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline ivector operator-(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline ivector operator-(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline sivector operator|(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) | v2;
}

inline sivector operator|(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) | v2;
}

inline sivector operator|(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) | v2;
}

inline sivector operator|(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) | v2;
}

inline sivector operator|(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) | v2;
}

inline sivector operator|(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) | v2;
}

inline ivector operator|(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) | v2;
}

inline ivector operator|(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) | v2;
}

inline ivector operator|(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) | v2;
}

inline ivector operator|(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) | v2;
}

inline ivector operator|(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) | v2;
}

inline ivector operator|(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) | v2;
}

inline sivector operator|(const sivector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline sivector operator|(const srvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline sivector operator|(const sivector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline sivector operator|(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline sivector operator|(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline sivector operator|(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline ivector operator|(const ivector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline ivector operator|(const rvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline ivector operator|(const ivector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline ivector operator|(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline ivector operator|(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline ivector operator|(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline sivector operator|(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) | v2;
}

inline sivector operator|(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) | v2;
}

inline ivector operator|(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) | v2;
}

inline ivector operator|(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) | v2;
}

inline sivector operator|(const srvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline sivector operator|(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline ivector operator|(const rvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline ivector operator|(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline simatrix_subv& simatrix_subv::operator*=(const real& v) {
  *this = *this * v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator*=(const interval& v) {
  *this = *this * v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator/=(const real& v) {
  *this = *this / v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator/=(const interval& v) {
  *this = *this / v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const srvector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const srvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const rvector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const rvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const srvector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const srvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const rvector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const rvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const sivector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const sivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const ivector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const ivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const sivector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const sivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const ivector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const ivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const srvector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const srvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const rvector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const rvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const sivector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const sivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const ivector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const ivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srmatrix_subv& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const simatrix_subv& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srvector& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const sivector& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srvector_slice& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const sivector_slice& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srmatrix_subv& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const simatrix_subv& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srvector& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const sivector& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srvector_slice& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const sivector_slice& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srmatrix_subv& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const simatrix_subv& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srvector& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const sivector& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srvector_slice& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const sivector_slice& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srmatrix_subv& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const simatrix_subv& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srvector& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const sivector& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srvector_slice& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const sivector_slice& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srmatrix_subv& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const simatrix_subv& v) {
  *this = ivector(v);
  return *this;
}

inline bool operator==(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const sivector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const srvector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const sivector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const ivector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const ivector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator!=(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const sivector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const srvector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const sivector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const ivector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const ivector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator<(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const srvector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const sivector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const rvector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const ivector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<=(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const srvector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const sivector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const rvector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const ivector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator>(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const sivector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const sivector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const ivector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const ivector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>=(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const sivector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const sivector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const ivector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const ivector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator!(const simatrix_subv& x) {
  return sv_v_not(x);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, sivector(v1), sivector(v2));
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, sivector(v1), srvector(v2));
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, srvector(v1), sivector(v2));
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const sivector& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srvector& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const sivector& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, srvector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const sivector_slice& v2) {
  spsl_vv_accu<idotprecision,sivector,sivector_slice,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accu<idotprecision,sivector,srvector_slice,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const sivector_slice& v2) {
  spsl_vv_accu<idotprecision,srvector,sivector_slice,sparse_idot>(dot, srvector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const ivector& v2) {
  spf_vv_accu<idotprecision,sivector,ivector,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const rvector& v2) {
  spf_vv_accu<idotprecision,sivector,rvector,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const ivector& v2) {
  spf_vv_accu<idotprecision,srvector,ivector,sparse_idot>(dot, srvector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const ivector_slice& v2) {
  spf_vv_accu<idotprecision,sivector,ivector_slice,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accu<idotprecision,sivector,rvector_slice,sparse_idot>(dot, sivector(v1), v2);
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const ivector_slice& v2) {
  spf_vv_accu<idotprecision,srvector,ivector_slice,sparse_idot>(dot, srvector(v1), v2);
}

inline void accumulate(idotprecision& dot, const sivector& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const sivector& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, v1, srvector(v2));
}

inline void accumulate(idotprecision& dot, const srvector& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const sivector_slice& v1, const simatrix_subv& v2) {
  slsp_vv_accu<idotprecision,sivector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const sivector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accu<idotprecision,sivector_slice,srvector,sparse_idot>(dot, v1, srvector(v2));
}

inline void accumulate(idotprecision& dot, const srvector_slice& v1, const simatrix_subv& v2) {
  slsp_vv_accu<idotprecision,srvector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const ivector& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const ivector& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector,srvector,sparse_idot>(dot, v1, srvector(v2));
}

inline void accumulate(idotprecision& dot, const rvector& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,rvector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const ivector_slice& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(idotprecision& dot, const ivector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector_slice,srvector,sparse_idot>(dot, v1, srvector(v2));
}

inline void accumulate(idotprecision& dot, const rvector_slice& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,rvector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const sivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srvector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const sivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const sivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srvector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const sivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const ivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const rvector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const ivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const ivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const rvector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const ivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

}  //namespace cxsc;

#include "sparsematrix.inl"

#endif 
 
