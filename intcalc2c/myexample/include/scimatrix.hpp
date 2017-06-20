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

/* CVS $Id: scimatrix.hpp,v 1.12 2010/12/21 09:54:01 cxsc Exp $ */

#ifndef _CXSC_SCIMATRIX_HPP_INCLUDED
#define _CXSC_SCIMATRIX_HPP_INCLUDED

#include <cinterval.hpp>
#include <cimatrix.hpp>
#include <civector.hpp>
#include <cidot.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sparsecidot.hpp>
#include <sparsematrix.hpp>
#include <srmatrix.hpp>
#include <scmatrix.hpp>
#include <simatrix.hpp>
#include <srvector.hpp>
#include <sivector.hpp>
#include <scvector.hpp>
#include <scivector.hpp>

namespace cxsc {

//definiert in srmatrix.hpp
//enum STORAGE_TYPE{triplet,compressed_row,compressed_column};

class scimatrix_slice;
class scimatrix_subv;

inline bool comp_pair_ci(std::pair<int,cinterval> p1, std::pair<int,cinterval> p2) {
  return p1.first < p2.first;
}

class scimatrix {

  private:
    std::vector<int> p;
    std::vector<int> ind;
    std::vector<cinterval> x;
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

    std::vector<cinterval>& values() {
      return x;
    }

    const std::vector<int>& column_pointers() const {
      return p;
    }

    const std::vector<int>& row_indices() const {
      return ind;
    }

    const std::vector<cinterval>& values() const {
      return x;
    }

    scimatrix() {
      p.push_back(0);
      m = n = 0;
      lb1 = lb2 = ub1 = ub2 = 0;
    }

    scimatrix(const int r, const int c) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(2*(m+n));
      x.reserve(2*(m+n));

      p[0] = 0;
    }

    scimatrix(const int r, const int c, const int e) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(e);
      x.reserve(e);

      p[0] = 0;
    }

    scimatrix(const int m, const int n, const int nnz, const intvector& rows, const intvector& cols, const civector& values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<cinterval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<cinterval>(rows[Lb(rows)+k],cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<triplet_store<cinterval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<cinterval>(j,cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<std::pair<int,cinterval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[Lb(cols)+k],values[Lb(values)+k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_ci);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    scimatrix(const int m, const int n, const int nnz, const int* rows, const int* cols, const cinterval* values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<cinterval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<cinterval>(rows[k],cols[k],values[k]));
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

         std::vector<triplet_store<cinterval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<cinterval>(j,cols[k],values[k]));
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

         std::vector<std::pair<int,cinterval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[k],values[k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_ci);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    scimatrix(const srmatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
    }

    scimatrix(const scmatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
    }

    scimatrix(const simatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
    }

    scimatrix(const rmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(cinterval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scimatrix(const cmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(cinterval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scimatrix(const imatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(cinterval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scimatrix(const cimatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(A[i+lb1][j+lb2]);
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    scimatrix(const srmatrix_slice&);
    scimatrix(const scmatrix_slice&);
    scimatrix(const simatrix_slice&);
    scimatrix(const scimatrix_slice&);

    void full(cimatrix& A) const {
       A = cimatrix(lb1,ub1,lb2,ub2);
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
      std::vector<cinterval> xnew;
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


    scimatrix& operator=(const real& A) {
      return sp_ms_assign<scimatrix,real,cinterval>(*this,A);
    }

    scimatrix& operator=(const interval& A) {
      return sp_ms_assign<scimatrix,interval,cinterval>(*this,A);
    }

    scimatrix& operator=(const complex& A) {
      return sp_ms_assign<scimatrix,complex,cinterval>(*this,A);
    }

    scimatrix& operator=(const cinterval& A) {
      return sp_ms_assign<scimatrix,cinterval,cinterval>(*this,A);
    }

    scimatrix& operator=(const rmatrix& A) {
      return spf_mm_assign<scimatrix,rmatrix,cinterval>(*this,A);
    }

    scimatrix& operator=(const cmatrix& A) {
      return spf_mm_assign<scimatrix,cmatrix,cinterval>(*this,A);
    }

    scimatrix& operator=(const imatrix& A) {
      return spf_mm_assign<scimatrix,imatrix,cinterval>(*this,A);
    }

    scimatrix& operator=(const cimatrix& A) {
      return spf_mm_assign<scimatrix,cimatrix,cinterval>(*this,A);
    }

    scimatrix& operator=(const rmatrix_slice& A) {
      return spf_mm_assign<scimatrix,rmatrix_slice,cinterval>(*this,A);
    }

    scimatrix& operator=(const cmatrix_slice& A) {
      return spf_mm_assign<scimatrix,cmatrix_slice,cinterval>(*this,A);
    }

    scimatrix& operator=(const imatrix_slice& A) {
      return spf_mm_assign<scimatrix,imatrix_slice,cinterval>(*this,A);
    }

    scimatrix& operator=(const cimatrix_slice& A) {
      return spf_mm_assign<scimatrix,cimatrix_slice,cinterval>(*this,A);
    }

    scimatrix& operator=(const srmatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
      return *this;
    }

    scimatrix& operator=(const scmatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
      return *this;
    }

    scimatrix& operator=(const simatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(cinterval(A.x[i]));
      return *this;
    }

    /* scimatrix& operator=(const scimatrix& A) {
      p = A.p;
      ind = A.ind;
      x = A.x;
      return *this;
    } */

    scimatrix& operator=(const srmatrix_slice&);
    scimatrix& operator=(const scmatrix_slice&);
    scimatrix& operator=(const simatrix_slice&);
    scimatrix& operator=(const scimatrix_slice&);

    const cinterval operator()(int i, int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix::operator()(int, int)"));
#endif
      cinterval r(0.0);
      for(int k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  r = x[k];
      }
      return r;
    }

    cinterval& element(int i, int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix::element(int, int)"));
#endif
      int k;
      for(k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  return x[k];
      }

      //Nicht gefunden, Element muss angelegt werden, da Schreibzugriff moeglich
      std::vector<int>::iterator ind_it = ind.begin() + k;
      std::vector<cinterval>::iterator x_it  = x.begin() + k;
      ind.insert(ind_it, i-lb1);
      x_it = x.insert(x_it, cinterval(0.0));
      for(k=j-lb2+1 ; k<(int)p.size() ; k++)
        p[k]++;

      return *x_it;
    }

    scimatrix_subv operator[](const cxscmatrix_column&);
    scimatrix_subv operator[](const int);
    const scimatrix_subv operator[](const cxscmatrix_column&) const;
    const scimatrix_subv operator[](const int) const;

    scimatrix_slice operator()(const int, const int , const int, const int);

    scimatrix operator()(const intvector& pervec, const intvector& q) {
      scimatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      int nnz=0;
      for(int k=0 ; k<n ; k++) {
        A.p[k] = nnz;

        std::map<int,cinterval> work;
        for(int j=p[q[Lb(q)+k]] ; j<p[q[Lb(q)+k]+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,cinterval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }

        nnz += work.size();
 
      }

      A.p[n] = nnz;

      return A;
    }

    scimatrix operator()(const intvector& pervec) {
      scimatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      for(int k=0 ; k<n ; k++) {
        A.p[k] = p[k];

        std::map<int,cinterval> work;
        for(int j=p[k] ; j<p[k+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,cinterval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }
 
      }

      A.p[n] = p[n];

      return A;
    }

    scimatrix operator()(const intmatrix& P, const intmatrix& Q) {
      intvector p = permvec(P);
      intvector q = perminv(permvec(Q));
      return (*this)(p,q);
    }

    scimatrix operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    real density() const {
      return p[n]/((double)m*n);
    }

    int get_nnz() const {
      return p[n];
    }

    scimatrix& operator+=(const rmatrix& B) {
      return spf_mm_addassign<scimatrix,rmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const cmatrix& B) {
      return spf_mm_addassign<scimatrix,cmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const imatrix& B) {
      return spf_mm_addassign<scimatrix,imatrix,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const cimatrix& B) {
      return spf_mm_addassign<scimatrix,cimatrix,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const rmatrix_slice& B) {
      return spf_mm_addassign<scimatrix,rmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const cmatrix_slice& B) {
      return spf_mm_addassign<scimatrix,cmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const imatrix_slice& B) {
      return spf_mm_addassign<scimatrix,imatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const cimatrix_slice& B) {
      return spf_mm_addassign<scimatrix,cimatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator+=(const srmatrix& B) {
      return spsp_mm_addassign<scimatrix,srmatrix,cinterval>(*this,B);
    }

    scimatrix& operator+=(const scmatrix& B) {
      return spsp_mm_addassign<scimatrix,scmatrix,cinterval>(*this,B);
    }

    scimatrix& operator+=(const simatrix& B) {
      return spsp_mm_addassign<scimatrix,simatrix,cinterval>(*this,B);
    }

    scimatrix& operator+=(const scimatrix& B) {
      return spsp_mm_addassign<scimatrix,scimatrix,cinterval>(*this,B);
    }

    scimatrix& operator-=(const rmatrix& B) {
      return spf_mm_subassign<scimatrix,rmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const cmatrix& B) {
      return spf_mm_subassign<scimatrix,cmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const imatrix& B) {
      return spf_mm_subassign<scimatrix,imatrix,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const cimatrix& B) {
      return spf_mm_subassign<scimatrix,cimatrix,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const rmatrix_slice& B) {
      return spf_mm_subassign<scimatrix,rmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const cmatrix_slice& B) {
      return spf_mm_subassign<scimatrix,cmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const imatrix_slice& B) {
      return spf_mm_subassign<scimatrix,imatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const cimatrix_slice& B) {
      return spf_mm_subassign<scimatrix,cimatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator-=(const srmatrix& B) {
      return spsp_mm_subassign<scimatrix,srmatrix,cinterval>(*this,B);
    }

    scimatrix& operator-=(const scmatrix& B) {
      return spsp_mm_subassign<scimatrix,scmatrix,cinterval>(*this,B);
    }

    scimatrix& operator-=(const simatrix& B) {
      return spsp_mm_subassign<scimatrix,simatrix,cinterval>(*this,B);
    }

    scimatrix& operator-=(const scimatrix& B) {
      return spsp_mm_subassign<scimatrix,scimatrix,cinterval>(*this,B);
    }

    scimatrix& operator|=(const rmatrix& B) {
      return spf_mm_hullassign<scimatrix,rmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const cmatrix& B) {
      return spf_mm_hullassign<scimatrix,cmatrix,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const imatrix& B) {
      return spf_mm_hullassign<scimatrix,imatrix,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const cimatrix& B) {
      return spf_mm_hullassign<scimatrix,cimatrix,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const rmatrix_slice& B) {
      return spf_mm_hullassign<scimatrix,rmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const cmatrix_slice& B) {
      return spf_mm_hullassign<scimatrix,cmatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const imatrix_slice& B) {
      return spf_mm_hullassign<scimatrix,imatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const cimatrix_slice& B) {
      return spf_mm_hullassign<scimatrix,cimatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator|=(const srmatrix& B) {
      return spsp_mm_hullassign<scimatrix,srmatrix,cinterval>(*this,B);
    }

    scimatrix& operator|=(const scmatrix& B) {
      return spsp_mm_hullassign<scimatrix,scmatrix,cinterval>(*this,B);
    }

    scimatrix& operator|=(const simatrix& B) {
      return spsp_mm_hullassign<scimatrix,simatrix,cinterval>(*this,B);
    }

    scimatrix& operator|=(const scimatrix& B) {
      return spsp_mm_hullassign<scimatrix,scimatrix,cinterval>(*this,B);
    }

    scimatrix& operator&=(const imatrix& B) {
      return spf_mm_intersectassign<scimatrix,imatrix,cimatrix>(*this,B);
    }

    scimatrix& operator&=(const cimatrix& B) {
      return spf_mm_intersectassign<scimatrix,cimatrix,cimatrix>(*this,B);
    }

    scimatrix& operator&=(const imatrix_slice& B) {
      return spf_mm_intersectassign<scimatrix,imatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator&=(const cimatrix_slice& B) {
      return spf_mm_intersectassign<scimatrix,cimatrix_slice,cimatrix>(*this,B);
    }

    scimatrix& operator&=(const simatrix& B) {
      return spsp_mm_intersectassign<scimatrix,simatrix,cinterval>(*this,B);
    }

    scimatrix& operator&=(const scimatrix& B) {
      return spsp_mm_intersectassign<scimatrix,scimatrix,cinterval>(*this,B);
    }

    scimatrix& operator*=(const cmatrix& B) {
      return spf_mm_multassign<scimatrix,cmatrix,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const rmatrix& B) {
      return spf_mm_multassign<scimatrix,rmatrix,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const imatrix& B) {
      return spf_mm_multassign<scimatrix,imatrix,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const cimatrix& B) {
      return spf_mm_multassign<scimatrix,cimatrix,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const rmatrix_slice& B) {
      return spf_mm_multassign<scimatrix,rmatrix_slice,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const cmatrix_slice& B) {
      return spf_mm_multassign<scimatrix,cmatrix_slice,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const imatrix_slice& B) {
      return spf_mm_multassign<scimatrix,imatrix_slice,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const cimatrix_slice& B) {
      return spf_mm_multassign<scimatrix,cimatrix_slice,sparse_cidot,cimatrix>(*this,B);
    }

    scimatrix& operator*=(const srmatrix& B) {
      return spsp_mm_multassign<scimatrix,srmatrix,sparse_cidot,cinterval>(*this,B);
    }

    scimatrix& operator*=(const scmatrix& B) {
      return spsp_mm_multassign<scimatrix,scmatrix,sparse_cidot,cinterval>(*this,B);
    }

    scimatrix& operator*=(const simatrix& B) {
      return spsp_mm_multassign<scimatrix,simatrix,sparse_cidot,cinterval>(*this,B);
    }

    scimatrix& operator*=(const scimatrix& B) {
      return spsp_mm_multassign<scimatrix,scimatrix,sparse_cidot,cinterval>(*this,B);
    }

    scimatrix& operator*=(const real& r) {
      return sp_ms_multassign(*this,r);
    }

    scimatrix& operator*=(const complex& r) {
      return sp_ms_multassign(*this,r);
    }

    scimatrix& operator*=(const interval& r) {
      return sp_ms_multassign(*this,r);
    }

    scimatrix& operator*=(const cinterval& r) {
      return sp_ms_multassign(*this,r);
    }

    scimatrix& operator/=(const real& r) {
      return sp_ms_divassign(*this,r);
    }

    scimatrix& operator/=(const complex& r) {
      return sp_ms_divassign(*this,r);
    }

    scimatrix& operator/=(const interval& r) {
      return sp_ms_divassign(*this,r);
    }

    scimatrix& operator/=(const cinterval& r) {
      return sp_ms_divassign(*this,r);
    }

    friend void SetLb(scimatrix&, const int, const int);
    friend void SetUb(scimatrix&, const int, const int);
    friend int Lb(const scimatrix&, int);
    friend int Ub(const scimatrix&, int);
    friend int RowLen(const scimatrix&);
    friend int ColLen(const scimatrix&);
    friend simatrix Re(const scimatrix&);
    friend simatrix Im(const scimatrix&);
    friend scmatrix Sup(const scimatrix&);
    friend scmatrix Inf(const scimatrix&);
    friend srmatrix InfRe(const scimatrix&);
    friend srmatrix InfIm(const scimatrix&);
    friend srmatrix SupRe(const scimatrix&);
    friend srmatrix SupIm(const scimatrix&);
    friend scimatrix conj(const scimatrix&);
    friend simatrix abs(const scimatrix&);
    friend scmatrix mid(const scimatrix&);
    friend scmatrix diam(const scimatrix&);

    friend scimatrix transp(const scimatrix&);
    friend scimatrix Id(const scimatrix&);

    friend std::istream& operator>>(std::istream&, scimatrix_slice&);
    friend std::istream& operator>>(std::istream&, scimatrix_subv&);

    friend class srmatrix_slice;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix_slice;
    friend class scmatrix_subv;
    friend class scvector;
    friend class simatrix_slice;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class scivector;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"
};

inline cimatrix::cimatrix(const srmatrix& A) {
  dat = new cinterval[A.m*A.n];
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

inline cimatrix::cimatrix(const simatrix& A) {
  dat = new cinterval[A.m*A.n];
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

inline cimatrix::cimatrix(const scmatrix& A) {
  dat = new cinterval[A.m*A.n];
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

inline cimatrix::cimatrix(const scimatrix& A) {
  dat = new cinterval[A.m*A.n];
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

inline scimatrix Id(const scimatrix& A) {
  scimatrix I(A.m, A.n, (A.m>A.n) ? A.m : A.n);
  I.lb1 = A.lb1; I.lb2 = A.lb2;
  I.ub1 = A.ub1; I.ub2 = A.ub2;

  if(A.m < A.n) {
    for(int i=0 ; i<A.m ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(cinterval(1.0));
    }
  } else {
    for(int i=0 ; i<A.n ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(cinterval(1.0));
    }
  }

  return I;
}

inline scimatrix transp(const scimatrix& A) {
  scimatrix B(A.n, A.m, A.get_nnz());
     
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

inline void SetLb(scimatrix& A, const int i, const int j) {
  if(i==1) {
    A.lb1 = j;
    A.ub1 = j + A.m - 1;
  } else if(i==2) {
    A.lb2 = j;
    A.ub2 = j + A.n - 1;
  }
}

inline void SetUb(scimatrix& A, const int i, const int j) {
  if(i==1) {
    A.ub1 = j;
    A.lb1 = j - A.m + 1;
  } else if(i==2) {
    A.ub2 = j;
    A.lb2 = j - A.n + 1;
  }
}


inline int Lb(const scimatrix& A, int i) {
  if(i==1) 
    return A.lb1;
  else if(i==2)
    return A.lb2;
  else
    return 1;
}

inline int Ub(const scimatrix& A, int i) {
  if(i==1) 
    return A.ub1;
  else if(i==2)
    return A.ub2;
  else
    return 1;
}

inline int RowLen(const scimatrix& A) {
  return A.n;
}

inline int ColLen(const scimatrix& A) {
  return A.m;
}

inline void Resize(scimatrix& A) {
  sp_m_resize(A);
}

inline void Resize(scimatrix& A, const int m, const int n) {
  sp_m_resize(A,m,n);
}

inline void Resize(scimatrix& A, const int l1, const int u1, const int l2, const int u2) {
  sp_m_resize(A,l1,u1,l2,u2);
}

inline simatrix Re(const scimatrix& A) {
  simatrix res(A.m,A.n,A.get_nnz());
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

inline simatrix Im(const scimatrix& A) {
  simatrix res(A.m,A.n,A.get_nnz());
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

inline scmatrix Inf(const scimatrix& A) {
  scmatrix res(A.m,A.n,A.get_nnz());
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

inline scmatrix Sup(const scimatrix& A) {
  scmatrix res(A.m,A.n,A.get_nnz());
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

inline scimatrix conj(const scimatrix& A) {
  scimatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(conj(A.x[i]));

  res.dropzeros();

  return res; 
}

inline simatrix abs(const scimatrix& A) {
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

inline scmatrix mid(const scimatrix& A) {
  scmatrix res(A.m,A.n,A.get_nnz());
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

inline scmatrix diam(const scimatrix& A) {
  scmatrix res(A.m,A.n,A.get_nnz());
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

inline srmatrix InfRe(const scimatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(InfRe(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix InfIm(const scimatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(InfIm(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix SupRe(const scimatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(SupRe(A.x[i]));

  res.dropzeros();

  return res; 
}

inline srmatrix SupIm(const scimatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(SupIm(A.x[i]));

  res.dropzeros();

  return res; 
}

inline cimatrix operator*(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_mult<cimatrix,srmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_mult<cimatrix,scmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix& A, const simatrix& B) {
  return fsp_mm_mult<cimatrix,simatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_mult<cimatrix,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_mult<rmatrix,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_mult<cmatrix,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const imatrix& A, const scimatrix& B) {
  return fsp_mm_mult<imatrix,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cmatrix& A, const simatrix& B) {
  return fsp_mm_mult<cmatrix,simatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const imatrix& A, const scmatrix& B) {
  return fsp_mm_mult<imatrix,scmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const rmatrix& B) {
  return spf_mm_mult<scimatrix,rmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const cmatrix& B) {
  return spf_mm_mult<scimatrix,cmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const imatrix& B) {
  return spf_mm_mult<scimatrix,imatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const cimatrix& B) {
  return spf_mm_mult<scimatrix,cimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const srmatrix& A, const cimatrix& B) {
  return spf_mm_mult<srmatrix,cimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scmatrix& A, const cimatrix& B) {
  return spf_mm_mult<scmatrix,cimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const simatrix& A, const cimatrix& B) {
  return spf_mm_mult<simatrix,cimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scmatrix& A, const imatrix& B) {
  return spf_mm_mult<scmatrix,imatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const simatrix& A, const cmatrix& B) {
  return spf_mm_mult<simatrix,cmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_mult<cimatrix_slice,srmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_mult<cimatrix_slice,scmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<cimatrix_slice,simatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_mult<cimatrix_slice,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_mult<rmatrix_slice,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_mult<imatrix_slice,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_mult<cmatrix_slice,scimatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const cmatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<cmatrix_slice,simatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const imatrix_slice& A, const scmatrix& B) {
  return fsp_mm_mult<imatrix_slice,scmatrix,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_mult<scimatrix,rmatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_mult<scimatrix,cmatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<scimatrix,imatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_mult<scimatrix,cimatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_mult<srmatrix,cimatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_mult<scmatrix,cimatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_mult<simatrix,cimatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const scmatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<scmatrix,imatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline cimatrix operator*(const simatrix& A, const cmatrix_slice& B) {
  return spf_mm_mult<simatrix,cmatrix_slice,cimatrix,sparse_cidot>(A,B);
}

inline scimatrix operator*(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_mult<scimatrix,srmatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_mult<scimatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const scimatrix& A, const simatrix& B) {
  return spsp_mm_mult<scimatrix,simatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_mult<scimatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_mult<srmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_mult<scmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const simatrix& A, const scimatrix& B) {
  return spsp_mm_mult<simatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const scmatrix& A, const simatrix& B) {
  return spsp_mm_mult<scmatrix,simatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator*(const simatrix& A, const scmatrix& B) {
  return spsp_mm_mult<simatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(A,B);
}

inline scimatrix operator/(const scimatrix& A, const real& r) {
  return sp_ms_div<scimatrix,real,scimatrix>(A,r);
}

inline scimatrix operator/(const scimatrix& A, const complex& r) {
  return sp_ms_div<scimatrix,complex,scimatrix>(A,r);
}

inline scimatrix operator/(const scimatrix& A, const interval& r) {
  return sp_ms_div<scimatrix,interval,scimatrix>(A,r);
}

inline scimatrix operator/(const scimatrix& A, const cinterval& r) {
  return sp_ms_div<scimatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator/(const srmatrix& A, const cinterval& r) {
  return sp_ms_div<srmatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator/(const simatrix& A, const cinterval& r) {
  return sp_ms_div<simatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator/(const scmatrix& A, const cinterval& r) {
  return sp_ms_div<scmatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator/(const scmatrix& A, const interval& r) {
  return sp_ms_div<scmatrix,interval,scimatrix>(A,r);
}

inline scimatrix operator/(const simatrix& A, const complex& r) {
  return sp_ms_div<simatrix,complex,scimatrix>(A,r);
}


inline scimatrix operator*(const scimatrix& A, const real& r) {
  return sp_ms_mult<scimatrix,real,scimatrix>(A,r);
}

inline scimatrix operator*(const scimatrix& A, const complex& r) {
  return sp_ms_mult<scimatrix,complex,scimatrix>(A,r);
}

inline scimatrix operator*(const scimatrix& A, const interval& r) {
  return sp_ms_mult<scimatrix,interval,scimatrix>(A,r);
}

inline scimatrix operator*(const scimatrix& A, const cinterval& r) {
  return sp_ms_mult<scimatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator*(const srmatrix& A, const cinterval& r) {
  return sp_ms_mult<srmatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator*(const simatrix& A, const cinterval& r) {
  return sp_ms_mult<simatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator*(const scmatrix& A, const cinterval& r) {
  return sp_ms_mult<scmatrix,cinterval,scimatrix>(A,r);
}

inline scimatrix operator*(const scmatrix& A, const interval& r) {
  return sp_ms_mult<scmatrix,interval,scimatrix>(A,r);
}

inline scimatrix operator*(const simatrix& A, const complex& r) {
  return sp_ms_mult<simatrix,complex,scimatrix>(A,r);
}

inline scimatrix operator*(const real& r, const scimatrix& A) {
  return sp_sm_mult<real,scimatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const complex& r, const scimatrix& A) {
  return sp_sm_mult<complex,scimatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const interval& r, const scimatrix& A) {
  return sp_sm_mult<interval,scimatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const cinterval& r, const scimatrix& A) {
  return sp_sm_mult<cinterval,scimatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const cinterval& r, const srmatrix& A) {
  return sp_sm_mult<cinterval,srmatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const cinterval& r, const simatrix& A) {
  return sp_sm_mult<cinterval,simatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const cinterval& r, const scmatrix& A) {
  return sp_sm_mult<cinterval,scmatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const complex& r, const simatrix& A) {
  return sp_sm_mult<complex,simatrix,scimatrix>(r,A);
}

inline scimatrix operator*(const interval& r, const scmatrix& A) {
  return sp_sm_mult<interval,scmatrix,scimatrix>(r,A);
}

inline civector operator*(const scimatrix& A, const rvector& v) {
  return spf_mv_mult<scimatrix,rvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const cvector& v) {
  return spf_mv_mult<scimatrix,cvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const ivector& v) {
  return spf_mv_mult<scimatrix,ivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const civector& v) {
  return spf_mv_mult<scimatrix,civector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const srmatrix& A, const civector& v) {
  return spf_mv_mult<srmatrix,civector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scmatrix& A, const civector& v) {
  return spf_mv_mult<scmatrix,civector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const simatrix& A, const civector& v) {
  return spf_mv_mult<simatrix,civector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scmatrix& A, const ivector& v) {
  return spf_mv_mult<scmatrix,ivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const simatrix& A, const cvector& v) {
  return spf_mv_mult<simatrix,cvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const rvector_slice& v) {
  return spf_mv_mult<scimatrix,rvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const ivector_slice& v) {
  return spf_mv_mult<scimatrix,ivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const cvector_slice& v) {
  return spf_mv_mult<scimatrix,cvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scimatrix& A, const civector_slice& v) {
  return spf_mv_mult<scimatrix,civector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const srmatrix& A, const civector_slice& v) {
  return spf_mv_mult<srmatrix,civector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scmatrix& A, const civector_slice& v) {
  return spf_mv_mult<scmatrix,civector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const simatrix& A, const civector_slice& v) {
  return spf_mv_mult<simatrix,civector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const simatrix& A, const cvector_slice& v) {
  return spf_mv_mult<simatrix,cvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const scmatrix& A, const ivector_slice& v) {
  return spf_mv_mult<scmatrix,ivector_slice,civector,sparse_cidot>(A,v);
}

inline scivector operator*(const scimatrix& A, const srvector& v) {
  return spsp_mv_mult<scimatrix,srvector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const sivector& v) {
  return spsp_mv_mult<scimatrix,sivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const scvector& v) {
  return spsp_mv_mult<scimatrix,scvector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const scivector& v) {
  return spsp_mv_mult<scimatrix,scivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const srmatrix& A, const scivector& v) {
  return spsp_mv_mult<srmatrix,scivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scmatrix& A, const scivector& v) {
  return spsp_mv_mult<scmatrix,scivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const simatrix& A, const scivector& v) {
  return spsp_mv_mult<simatrix,scivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scmatrix& A, const sivector& v) {
  return spsp_mv_mult<scmatrix,sivector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const simatrix& A, const scvector& v) {
  return spsp_mv_mult<simatrix,scvector,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const srvector_slice& v) {
  return spsl_mv_mult<scimatrix,srvector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const scvector_slice& v) {
  return spsl_mv_mult<scimatrix,scvector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<scimatrix,sivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scimatrix& A, const scivector_slice& v) {
  return spsl_mv_mult<scimatrix,scivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const srmatrix& A, const scivector_slice& v) {
  return spsl_mv_mult<srmatrix,scivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scmatrix& A, const scivector_slice& v) {
  return spsl_mv_mult<scmatrix,scivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const simatrix& A, const scivector_slice& v) {
  return spsl_mv_mult<simatrix,scivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const simatrix& A, const scvector_slice& v) {
  return spsl_mv_mult<simatrix,scvector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline scivector operator*(const scmatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<scmatrix,sivector_slice,scivector,sparse_cidot,cinterval>(A,v);
}

inline civector operator*(const cimatrix& A, const srvector& v) {
  return fsp_mv_mult<cimatrix,srvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const sivector& v) {
  return fsp_mv_mult<cimatrix,sivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const scvector& v) {
  return fsp_mv_mult<cimatrix,scvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const scivector& v) {
  return fsp_mv_mult<cimatrix,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const rmatrix& A, const scivector& v) {
  return fsp_mv_mult<rmatrix,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix& A, const scivector& v) {
  return fsp_mv_mult<cmatrix,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix& A, const scivector& v) {
  return fsp_mv_mult<imatrix,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix& A, const sivector& v) {
  return fsp_mv_mult<cmatrix,sivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix& A, const scvector& v) {
  return fsp_mv_mult<imatrix,scvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const srvector& v) {
  return fsp_mv_mult<cimatrix_slice,srvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const scvector& v) {
  return fsp_mv_mult<cimatrix_slice,scvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<cimatrix_slice,sivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const scivector& v) {
  return fsp_mv_mult<cimatrix_slice,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const rmatrix_slice& A, const scivector& v) {
  return fsp_mv_mult<rmatrix_slice,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix_slice& A, const scivector& v) {
  return fsp_mv_mult<cmatrix_slice,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix_slice& A, const scivector& v) {
  return fsp_mv_mult<imatrix_slice,scivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<cmatrix_slice,sivector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix_slice& A, const scvector& v) {
  return fsp_mv_mult<imatrix_slice,scvector,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const srvector_slice& v) {
  return fsl_mv_mult<cimatrix,srvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const scvector_slice& v) {
  return fsl_mv_mult<cimatrix,scvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<cimatrix,sivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix& A, const scivector_slice& v) {
  return fsl_mv_mult<cimatrix,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const rmatrix& A, const scivector_slice& v) {
  return fsl_mv_mult<rmatrix,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix& A, const scivector_slice& v) {
  return fsl_mv_mult<cmatrix,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix& A, const scivector_slice& v) {
  return fsl_mv_mult<imatrix,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<cmatrix,sivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix& A, const scvector_slice& v) {
  return fsl_mv_mult<imatrix,scvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const srvector_slice& v) {
  return fsl_mv_mult<cimatrix_slice,srvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const scvector_slice& v) {
  return fsl_mv_mult<cimatrix_slice,scvector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<cimatrix_slice,sivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cimatrix_slice& A, const scivector_slice& v) {
  return fsl_mv_mult<cimatrix_slice,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const rmatrix_slice& A, const scivector_slice& v) {
  return fsl_mv_mult<rmatrix_slice,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix_slice& A, const scivector_slice& v) {
  return fsl_mv_mult<imatrix_slice,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix_slice& A, const scivector_slice& v) {
  return fsl_mv_mult<cmatrix_slice,scivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const cmatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<cmatrix_slice,sivector_slice,civector,sparse_cidot>(A,v);
}

inline civector operator*(const imatrix_slice& A, const scvector_slice& v) {
  return fsl_mv_mult<imatrix_slice,scvector_slice,civector,sparse_cidot>(A,v);
}

inline cimatrix operator+(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_add<cimatrix,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_add<cimatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix& A, const simatrix& B) {
  return fsp_mm_add<cimatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_add<cimatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_add<rmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_add<cmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const imatrix& A, const scimatrix& B) {
  return fsp_mm_add<imatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const imatrix& A, const scmatrix& B) {
  return fsp_mm_add<imatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cmatrix& A, const simatrix& B) {
  return fsp_mm_add<cmatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const rmatrix& B) {
  return spf_mm_add<scimatrix,rmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const cmatrix& B) {
  return spf_mm_add<scimatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const imatrix& B) {
  return spf_mm_add<scimatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const cimatrix& B) {
  return spf_mm_add<scimatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const srmatrix& A, const cimatrix& B) {
  return spf_mm_add<srmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scmatrix& A, const cimatrix& B) {
  return spf_mm_add<scmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const simatrix& A, const cimatrix& B) {
  return spf_mm_add<simatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const simatrix& A, const cmatrix& B) {
  return spf_mm_add<simatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scmatrix& A, const imatrix& B) {
  return spf_mm_add<scmatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_add<cimatrix_slice,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<cimatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_add<cimatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_add<cimatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_add<rmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_add<imatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_add<cmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const cmatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<cmatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const imatrix_slice& A, const scmatrix& B) {
  return fsp_mm_add<imatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_add<scimatrix,rmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_add<scimatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_add<scimatrix,imatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_add<scimatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_add<srmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_add<scmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_add<simatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const simatrix& A, const cmatrix_slice& B) {
  return spf_mm_add<simatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator+(const scmatrix& A, const imatrix_slice& B) {
  return spf_mm_add<scmatrix,imatrix_slice,cimatrix>(A,B);
}

inline scimatrix operator+(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_add<scimatrix,srmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_add<scimatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const scimatrix& A, const simatrix& B) {
  return spsp_mm_add<scimatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_add<scimatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_add<srmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_add<scmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const simatrix& A, const scimatrix& B) {
  return spsp_mm_add<simatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const simatrix& A, const scmatrix& B) {
  return spsp_mm_add<simatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator+(const scmatrix& A, const simatrix& B) {
  return spsp_mm_add<scmatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline cimatrix operator-(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_sub<cimatrix,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_sub<cimatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix& A, const simatrix& B) {
  return fsp_mm_sub<cimatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_sub<cimatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_sub<rmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_sub<cmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const imatrix& A, const scimatrix& B) {
  return fsp_mm_sub<imatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const imatrix& A, const scmatrix& B) {
  return fsp_mm_sub<imatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cmatrix& A, const simatrix& B) {
  return fsp_mm_sub<cmatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const rmatrix& B) {
  return spf_mm_sub<scimatrix,rmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const cmatrix& B) {
  return spf_mm_sub<scimatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const imatrix& B) {
  return spf_mm_sub<scimatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const cimatrix& B) {
  return spf_mm_sub<scimatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const srmatrix& A, const cimatrix& B) {
  return spf_mm_sub<srmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scmatrix& A, const cimatrix& B) {
  return spf_mm_sub<scmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const simatrix& A, const cimatrix& B) {
  return spf_mm_sub<simatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const simatrix& A, const cmatrix& B) {
  return spf_mm_sub<simatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scmatrix& A, const imatrix& B) {
  return spf_mm_sub<scmatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_sub<cimatrix_slice,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<cimatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_sub<cimatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_sub<cimatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_sub<rmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_sub<imatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_sub<cmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const cmatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<cmatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const imatrix_slice& A, const scmatrix& B) {
  return fsp_mm_sub<imatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_sub<scimatrix,rmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_sub<scimatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<scimatrix,imatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_sub<scimatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_sub<srmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_sub<scmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_sub<simatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const simatrix& A, const cmatrix_slice& B) {
  return spf_mm_sub<simatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator-(const scmatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<scmatrix,imatrix_slice,cimatrix>(A,B);
}

inline scimatrix operator-(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_sub<scimatrix,srmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_sub<scimatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scimatrix& A, const simatrix& B) {
  return spsp_mm_sub<scimatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_sub<scimatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_sub<srmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_sub<scmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const simatrix& A, const scimatrix& B) {
  return spsp_mm_sub<simatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const simatrix& A, const scmatrix& B) {
  return spsp_mm_sub<simatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scmatrix& A, const simatrix& B) {
  return spsp_mm_sub<scmatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline cimatrix operator|(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_hull<cimatrix,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_hull<cimatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix& A, const simatrix& B) {
  return fsp_mm_hull<cimatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_hull<cimatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_hull<rmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_hull<cmatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const imatrix& A, const scimatrix& B) {
  return fsp_mm_hull<imatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const imatrix& A, const scmatrix& B) {
  return fsp_mm_hull<imatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix& A, const simatrix& B) {
  return fsp_mm_hull<cmatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const rmatrix& B) {
  return spf_mm_hull<scimatrix,rmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const cmatrix& B) {
  return spf_mm_hull<scimatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const imatrix& B) {
  return spf_mm_hull<scimatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const cimatrix& B) {
  return spf_mm_hull<scimatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const srmatrix& A, const cimatrix& B) {
  return spf_mm_hull<srmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const cimatrix& B) {
  return spf_mm_hull<scmatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const simatrix& A, const cimatrix& B) {
  return spf_mm_hull<simatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const simatrix& A, const cmatrix& B) {
  return spf_mm_hull<simatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const imatrix& B) {
  return spf_mm_hull<scmatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<cimatrix_slice,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<cimatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_hull<cimatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_hull<cimatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_hull<rmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_hull<imatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_hull<cmatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<cmatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const imatrix_slice& A, const scmatrix& B) {
  return fsp_mm_hull<imatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<scimatrix,rmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_hull<scimatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<scimatrix,imatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_hull<scimatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_hull<srmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_hull<scmatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_hull<simatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const simatrix& A, const cmatrix_slice& B) {
  return spf_mm_hull<simatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<scmatrix,imatrix_slice,cimatrix>(A,B);
}

inline scimatrix operator|(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_hull<scimatrix,srmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_hull<scimatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scimatrix& A, const simatrix& B) {
  return spsp_mm_hull<scimatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_hull<scimatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_hull<srmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_hull<scmatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const simatrix& A, const scimatrix& B) {
  return spsp_mm_hull<simatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const simatrix& A, const scmatrix& B) {
  return spsp_mm_hull<simatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scmatrix& A, const simatrix& B) {
  return spsp_mm_hull<scmatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline cimatrix operator|(const cmatrix& A, const srmatrix& B) {
  return fsp_mm_hull<cmatrix,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const rmatrix& A, const scmatrix& B) {
  return fsp_mm_hull<rmatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix& A, const scmatrix& B) {
  return fsp_mm_hull<cmatrix,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const rmatrix& B) {
  return spf_mm_hull<scmatrix,rmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const srmatrix& A, const cmatrix& B) {
  return spf_mm_hull<srmatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const cmatrix& B) {
  return spf_mm_hull<scmatrix,cmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<cmatrix_slice,srmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const rmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_hull<rmatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const cmatrix_slice& A, const scmatrix& B) {
  return fsp_mm_hull<cmatrix_slice,scmatrix,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<scmatrix,rmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const srmatrix& A, const cmatrix_slice& B) {
  return spf_mm_hull<srmatrix,cmatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator|(const scmatrix& A, const cmatrix_slice& B) {
  return spf_mm_hull<scmatrix,cmatrix_slice,cimatrix>(A,B);
}

inline scimatrix operator|(const scmatrix& A, const srmatrix& B) {
  return spsp_mm_hull<scmatrix,srmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const srmatrix& A, const scmatrix& B) {
  return spsp_mm_hull<srmatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator|(const scmatrix& A, const scmatrix& B) {
  return spsp_mm_hull<scmatrix,scmatrix,scimatrix,cinterval>(A,B);
}

inline cimatrix operator&(const cimatrix& A, const simatrix& B) {
  return fsp_mm_intersect<cimatrix,simatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_intersect<cimatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const imatrix& A, const scimatrix& B) {
  return fsp_mm_intersect<imatrix,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const scimatrix& A, const imatrix& B) {
  return spf_mm_intersect<scimatrix,imatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const scimatrix& A, const cimatrix& B) {
  return spf_mm_intersect<scimatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const simatrix& A, const cimatrix& B) {
  return spf_mm_intersect<simatrix,cimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_intersect<cimatrix_slice,simatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_intersect<cimatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_intersect<imatrix_slice,scimatrix,cimatrix>(A,B);
}

inline cimatrix operator&(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_intersect<scimatrix,imatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator&(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_intersect<scimatrix,cimatrix_slice,cimatrix>(A,B);
}

inline cimatrix operator&(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_intersect<simatrix,cimatrix_slice,cimatrix>(A,B);
}

inline scimatrix operator&(const scimatrix& A, const simatrix& B) {
  return spsp_mm_intersect<scimatrix,simatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator&(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_intersect<scimatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator&(const simatrix& A, const scimatrix& B) {
  return spsp_mm_intersect<simatrix,scimatrix,scimatrix,cinterval>(A,B);
}

inline scimatrix operator-(const scimatrix& M) {
  return sp_m_negative<scimatrix,scimatrix>(M);
}

inline scimatrix& operator+(scimatrix& A) {
  return A;
}

inline cimatrix& cimatrix::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline cimatrix& cimatrix::operator=(const scmatrix& B) {
  *this = cmatrix(B);
  return *this;
}

inline cimatrix& cimatrix::operator=(const simatrix& B) {
  *this = imatrix(B);
  return *this;
}

inline cimatrix& cimatrix::operator=(const scimatrix& B) {
  *this = cimatrix(B);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const scmatrix& B) {
  *this = cmatrix(B);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const simatrix& B) {
  *this = imatrix(B);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const scimatrix& B) {
  *this = cimatrix(B);
  return *this;
}

inline cimatrix& cimatrix::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix& cimatrix::operator+=(const scmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix& cimatrix::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix& cimatrix::operator+=(const scimatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator+=(const scmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator+=(const scimatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline cimatrix& cimatrix::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix& cimatrix::operator-=(const scmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix& cimatrix::operator-=(const scimatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix& cimatrix::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator-=(const scmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator-=(const scimatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline cimatrix& cimatrix::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix& cimatrix::operator|=(const scmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix& cimatrix::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix& cimatrix::operator|=(const scimatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator|=(const scmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator|=(const scimatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline cimatrix& cimatrix::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline cimatrix& cimatrix::operator&=(const scimatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator&=(const scimatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline cimatrix& cimatrix::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<cimatrix,srmatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix& cimatrix::operator*=(const scmatrix& B) {
  return fsp_mm_multassign<cimatrix,scmatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix& cimatrix::operator*=(const simatrix& B) {
  return fsp_mm_multassign<cimatrix,simatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix& cimatrix::operator*=(const scimatrix& B) {
  return fsp_mm_multassign<cimatrix,scimatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<cimatrix_slice,srmatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator*=(const scmatrix& B) {
  return fsp_mm_multassign<cimatrix_slice,scmatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator*=(const simatrix& B) {
  return fsp_mm_multassign<cimatrix_slice,simatrix,sparse_cidot,cimatrix>(*this,B);
}

inline cimatrix_slice& cimatrix_slice::operator*=(const scimatrix& B) {
  return fsp_mm_multassign<cimatrix_slice,scimatrix,sparse_cidot,cimatrix>(*this,B);
}

inline bool operator==(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const simatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const scimatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const rmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const imatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const cmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const cimatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const cimatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const cimatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const cimatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const imatrix& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const srmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const scmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const simatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const scimatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const scimatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const scimatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const scimatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const rmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const imatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const cmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const cimatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const cimatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const cimatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const cimatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const cimatrix& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cmatrix& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const imatrix& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix_slice& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix_slice& A, const scmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix_slice& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cimatrix_slice& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix_slice& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const cmatrix_slice& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const imatrix_slice& A, const scimatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const rmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const imatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const cmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scimatrix& A, const cimatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const cimatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const scmatrix& A, const cimatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const simatrix& A, const cimatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator<(const scimatrix& A, const simatrix& B) {
  return spsp_mm_less<scimatrix,simatrix,cinterval>(A,B);
}

inline bool operator<(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_less<scimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_less<srmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_less<scmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const simatrix& A, const scimatrix& B) {
  return spsp_mm_less<simatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const scimatrix& A, const imatrix& B) {
  return spf_mm_less<scimatrix,imatrix,cinterval>(A,B);
}

inline bool operator<(const scimatrix& A, const cimatrix& B) {
  return spf_mm_less<scimatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<(const srmatrix& A, const cimatrix& B) {
  return spf_mm_less<srmatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<(const scmatrix& A, const cimatrix& B) {
  return spf_mm_less<scmatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<(const simatrix& A, const cimatrix& B) {
  return spf_mm_less<simatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<(const cimatrix& A, const simatrix& B) {
  return fsp_mm_less<cimatrix,simatrix,cinterval>(A,B);
}

inline bool operator<(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_less<cimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_less<rmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_less<cmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const imatrix& A, const scimatrix& B) {
  return fsp_mm_less<imatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_less<cimatrix_slice,simatrix,cinterval>(A,B);
}

inline bool operator<(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_less<cimatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_less<rmatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_less<cmatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_less<imatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_less<scimatrix,imatrix_slice,cinterval>(A,B);
}

inline bool operator<(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_less<scimatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_less<srmatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_less<scmatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_less<simatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const simatrix& B) {
  return spsp_mm_leq<scimatrix,simatrix,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_leq<scimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const srmatrix& A, const scimatrix& B) {
  return spsp_mm_leq<srmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const scmatrix& A, const scimatrix& B) {
  return spsp_mm_leq<scmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const simatrix& A, const scimatrix& B) {
  return spsp_mm_leq<simatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const imatrix& B) {
  return spf_mm_leq<scimatrix,imatrix,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const cimatrix& B) {
  return spf_mm_leq<scimatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<=(const srmatrix& A, const cimatrix& B) {
  return spf_mm_leq<srmatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<=(const scmatrix& A, const cimatrix& B) {
  return spf_mm_leq<scmatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<=(const simatrix& A, const cimatrix& B) {
  return spf_mm_leq<simatrix,cimatrix,cinterval>(A,B);
}

inline bool operator<=(const cimatrix& A, const simatrix& B) {
  return fsp_mm_leq<cimatrix,simatrix,cinterval>(A,B);
}

inline bool operator<=(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_leq<cimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const rmatrix& A, const scimatrix& B) {
  return fsp_mm_leq<rmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const cmatrix& A, const scimatrix& B) {
  return fsp_mm_leq<cmatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const imatrix& A, const scimatrix& B) {
  return fsp_mm_leq<imatrix,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_leq<cimatrix_slice,simatrix,cinterval>(A,B);
}

inline bool operator<=(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_leq<cimatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const rmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_leq<rmatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const cmatrix_slice& A, const scimatrix& B) {
  return fsp_mm_leq<cmatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_leq<imatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_leq<scimatrix,imatrix_slice,cinterval>(A,B);
}

inline bool operator<=(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_leq<scimatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<=(const srmatrix& A, const cimatrix_slice& B) {
  return spf_mm_leq<srmatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<=(const scmatrix& A, const cimatrix_slice& B) {
  return spf_mm_leq<scmatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator<=(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_leq<simatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_greater<scimatrix,srmatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_greater<scimatrix,scmatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const simatrix& B) {
  return spsp_mm_greater<scimatrix,simatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_greater<scimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>(const simatrix& A, const scimatrix& B) {
  return spsp_mm_greater<simatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const rmatrix& B) {
  return spf_mm_greater<scimatrix,rmatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const imatrix& B) {
  return spf_mm_greater<scimatrix,imatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const cmatrix& B) {
  return spf_mm_greater<scimatrix,cmatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const cimatrix& B) {
  return spf_mm_greater<scimatrix,cimatrix,cinterval>(A,B);
}

inline bool operator>(const simatrix& A, const cimatrix& B) {
  return spf_mm_greater<simatrix,cimatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_greater<cimatrix,srmatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix& A, const simatrix& B) {
  return fsp_mm_greater<cimatrix,simatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_greater<cimatrix,scmatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_greater<cimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>(const imatrix& A, const scimatrix& B) {
  return fsp_mm_greater<imatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_greater<cimatrix_slice,srmatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_greater<cimatrix_slice,scmatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_greater<cimatrix_slice,simatrix,cinterval>(A,B);
}

inline bool operator>(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_greater<cimatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator>(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_greater<imatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_greater<scimatrix,rmatrix_slice,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_greater<scimatrix,imatrix_slice,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_greater<scimatrix,cmatrix_slice,cinterval>(A,B);
}

inline bool operator>(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_greater<scimatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator>(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_greater<simatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const srmatrix& B) {
  return spsp_mm_geq<scimatrix,srmatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const scmatrix& B) {
  return spsp_mm_geq<scimatrix,scmatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const simatrix& B) {
  return spsp_mm_geq<scimatrix,simatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const scimatrix& B) {
  return spsp_mm_geq<scimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const simatrix& A, const scimatrix& B) {
  return spsp_mm_geq<simatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const rmatrix& B) {
  return spf_mm_geq<scimatrix,rmatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const imatrix& B) {
  return spf_mm_geq<scimatrix,imatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const cmatrix& B) {
  return spf_mm_geq<scimatrix,cmatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const cimatrix& B) {
  return spf_mm_geq<scimatrix,cimatrix,cinterval>(A,B);
}

inline bool operator>=(const simatrix& A, const cimatrix& B) {
  return spf_mm_geq<simatrix,cimatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix& A, const srmatrix& B) {
  return fsp_mm_geq<cimatrix,srmatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix& A, const simatrix& B) {
  return fsp_mm_geq<cimatrix,simatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix& A, const scmatrix& B) {
  return fsp_mm_geq<cimatrix,scmatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix& A, const scimatrix& B) {
  return fsp_mm_geq<cimatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const imatrix& A, const scimatrix& B) {
  return fsp_mm_geq<imatrix,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix_slice& A, const srmatrix& B) {
  return fsp_mm_geq<cimatrix_slice,srmatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix_slice& A, const scmatrix& B) {
  return fsp_mm_geq<cimatrix_slice,scmatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix_slice& A, const simatrix& B) {
  return fsp_mm_geq<cimatrix_slice,simatrix,cinterval>(A,B);
}

inline bool operator>=(const cimatrix_slice& A, const scimatrix& B) {
  return fsp_mm_geq<cimatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const imatrix_slice& A, const scimatrix& B) {
  return fsp_mm_geq<imatrix_slice,scimatrix,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const rmatrix_slice& B) {
  return spf_mm_geq<scimatrix,rmatrix_slice,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const imatrix_slice& B) {
  return spf_mm_geq<scimatrix,imatrix_slice,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const cmatrix_slice& B) {
  return spf_mm_geq<scimatrix,cmatrix_slice,cinterval>(A,B);
}

inline bool operator>=(const scimatrix& A, const cimatrix_slice& B) {
  return spf_mm_geq<scimatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator>=(const simatrix& A, const cimatrix_slice& B) {
  return spf_mm_geq<simatrix,cimatrix_slice,cinterval>(A,B);
}

inline bool operator!(const scimatrix& A) {
  return sp_m_not(A);
}

inline std::ostream& operator<<(std::ostream& os, const scimatrix& A) {
  return sp_m_output<scimatrix,cinterval>(os,A);
}

inline std::istream& operator>>(std::istream& is, scimatrix& A) {
  return sp_m_input<scimatrix,cinterval>(is,A);
}

class scimatrix_slice {
  public:
    scimatrix  A;
    scimatrix* M; //Originalmatrix

  private:
    scimatrix_slice(scimatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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

    scimatrix_slice(const scimatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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
        M = const_cast<scimatrix*>(&Mat);
    }


  public:
    scimatrix_slice& operator=(const real& C) {
      return sl_ms_assign<scimatrix_slice, real, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const interval& C) {
      return sl_ms_assign<scimatrix_slice, interval, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const complex& C) {
      return sl_ms_assign<scimatrix_slice, complex, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const cinterval& C) {
      return sl_ms_assign<scimatrix_slice, cinterval, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const srmatrix& C) {
      return slsp_mm_assign<scimatrix_slice, srmatrix, std::vector<cinterval>::iterator>(*this,C);
    }

    scimatrix_slice& operator=(const scmatrix& C) {
      return slsp_mm_assign<scimatrix_slice, scmatrix, std::vector<cinterval>::iterator>(*this,C);
    }

    scimatrix_slice& operator=(const simatrix& C) {
      return slsp_mm_assign<scimatrix_slice, simatrix, std::vector<cinterval>::iterator>(*this,C);
    }

    scimatrix_slice& operator=(const scimatrix& C) {
      return slsp_mm_assign<scimatrix_slice, scimatrix, std::vector<cinterval>::iterator>(*this,C);
    }

    scimatrix_slice& operator=(const rmatrix& C) {
      return slf_mm_assign<scimatrix_slice, rmatrix, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const cmatrix& C) {
      return slf_mm_assign<scimatrix_slice, cmatrix, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const imatrix& C) {
      return slf_mm_assign<scimatrix_slice, imatrix, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const cimatrix& C) {
      return slf_mm_assign<scimatrix_slice, cimatrix, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const rmatrix_slice& C) {
      return slf_mm_assign<scimatrix_slice, rmatrix_slice, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const cmatrix_slice& C) {
      return slf_mm_assign<scimatrix_slice, cmatrix_slice, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const imatrix_slice& C) {
      return slf_mm_assign<scimatrix_slice, imatrix_slice, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const cimatrix_slice& C) {
      return slf_mm_assign<scimatrix_slice, cimatrix_slice, std::vector<cinterval>::iterator, cinterval>(*this,C);
    }

    scimatrix_slice& operator=(const srmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scimatrix_slice& operator=(const scmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scimatrix_slice& operator=(const simatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scimatrix_slice& operator=(const scimatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    scimatrix_slice& operator*=(const srmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scimatrix_slice& operator*=(const scmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scimatrix_slice& operator*=(const simatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scimatrix_slice& operator*=(const scimatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    scimatrix_slice& operator*=(const srmatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const scmatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const simatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const scimatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const rmatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const cmatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const imatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const cimatrix& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const rmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const cmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const imatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const cimatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    scimatrix_slice& operator*=(const real& r) {
      *this = A*r;
      return *this;
    }

    scimatrix_slice& operator*=(const complex& r) {
      *this = A*r;
      return *this;
    }

    scimatrix_slice& operator*=(const interval& r) {
      *this = A*r;
      return *this;
    }

    scimatrix_slice& operator*=(const cinterval& r) {
      *this = A*r;
      return *this;
    }

    scimatrix_slice& operator/=(const real& r) {
      *this = A/r;
      return *this;
    }

    scimatrix_slice& operator/=(const complex& r) {
      *this = A/r;
      return *this;
    }

    scimatrix_slice& operator/=(const interval& r) {
      *this = A/r;
      return *this;
    }

    scimatrix_slice& operator/=(const cinterval& r) {
      *this = A/r;
      return *this;
    }

    scimatrix_slice& operator+=(const srmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scimatrix_slice& operator+=(const scmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scimatrix_slice& operator+=(const simatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scimatrix_slice& operator+=(const scimatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    scimatrix_slice& operator+=(const srmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const scmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const simatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const scimatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const rmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const cmatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const imatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const cimatrix& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const rmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const cmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const imatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator+=(const cimatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    scimatrix_slice& operator-=(const srmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scimatrix_slice& operator-=(const scmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scimatrix_slice& operator-=(const simatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scimatrix_slice& operator-=(const scimatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    scimatrix_slice& operator-=(const srmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const scmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const simatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const scimatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const rmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const cmatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const imatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const cimatrix& M) {
      *this = A-M;
      return *this;
    } 

    scimatrix_slice& operator-=(const rmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    scimatrix_slice& operator-=(const cmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    scimatrix_slice& operator-=(const imatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    scimatrix_slice& operator-=(const cimatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    scimatrix_slice& operator|=(const srmatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    scimatrix_slice& operator|=(const scmatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    scimatrix_slice& operator|=(const simatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    scimatrix_slice& operator|=(const scimatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    scimatrix_slice& operator|=(const srmatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const scmatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const simatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const scimatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const rmatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const cmatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const imatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const cimatrix& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const rmatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const cmatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const imatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    scimatrix_slice& operator|=(const cimatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    const cinterval operator()(const int i, const int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_slice::operator()(int, int)"));
#endif
      cinterval r = A(i,j);
      return r;
    }

    cinterval& element(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_slice::element(int, int)"));
#endif
      return M->element(i,j);
    }

    scimatrix_subv operator[](const int);
    scimatrix_subv operator[](const cxscmatrix_column&);
    const scimatrix_subv operator[](const int) const;
    const scimatrix_subv operator[](const cxscmatrix_column&) const;

    friend std::ostream& operator<<(std::ostream&, const scimatrix_slice&);

    friend int Lb(const scimatrix_slice&, const int);
    friend int Ub(const scimatrix_slice&, const int);
    friend simatrix Re(const scimatrix_slice&);
    friend simatrix Im(const scimatrix_slice&);
    friend scmatrix Inf(const scimatrix_slice&);
    friend scmatrix Sup(const scimatrix_slice&);
    friend srmatrix InfRe(const scimatrix_slice&);
    friend srmatrix InfIm(const scimatrix_slice&);
    friend srmatrix SupRe(const scimatrix_slice&);
    friend srmatrix SupIm(const scimatrix_slice&);
    friend int RowLen(const scimatrix_slice&);
    friend int ColLen(const scimatrix_slice&);

    friend class srmatrix;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix;
    friend class scmatrix_subv;
    friend class scvector;
    friend class simatrix;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_subv;
    friend class scivector;

#include "matrix_friend_declarations.inl"    
};

inline cimatrix::cimatrix(const srmatrix_slice& A) {
  dat = new cinterval[A.A.m*A.A.n];
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

inline cimatrix::cimatrix(const scmatrix_slice& A) {
  dat = new cinterval[A.A.m*A.A.n];
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

inline cimatrix::cimatrix(const simatrix_slice& A) {
  dat = new cinterval[A.A.m*A.A.n];
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

inline cimatrix::cimatrix(const scimatrix_slice& A) {
  dat = new cinterval[A.A.m*A.A.n];
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

inline int RowLen(const scimatrix_slice& S) {
  return RowLen(S.A);
}

inline int ColLen(const scimatrix_slice& S) {
  return ColLen(S.A);
}

inline scimatrix_slice scimatrix::operator()(const int i, const int j, const int k, const int l) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix::operator()(int, int)"));
#endif
  return scimatrix_slice(*this, i, j, k, l);
}

inline scimatrix& scimatrix::operator=(const srmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline scimatrix& scimatrix::operator=(const scmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline scimatrix& scimatrix::operator=(const simatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline scimatrix& scimatrix::operator=(const scimatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline int Lb(const scimatrix_slice& S, const int i) {
  return Lb(S.A, i);
}

inline int Ub(const scimatrix_slice& S, const int i) {
  return Ub(S.A, i);
}

inline simatrix Re(const scimatrix_slice& S) {
  return Re(S.A);
}

inline simatrix Im(const scimatrix_slice& S) {
  return Im(S.A);
}

inline scimatrix conj(const scimatrix_slice& S) {
  return conj(S.A);
}

inline simatrix abs(const scimatrix_slice& S) {
  return abs(S.A);
}

inline scmatrix mid(const scimatrix_slice& S) {
  return mid(S.A);
}

inline scmatrix diam(const scimatrix_slice& S) {
  return diam(S.A);
}

inline scmatrix Inf(const scimatrix_slice& S) {
  return Inf(S.A);
}

inline scmatrix Sup(const scimatrix_slice& S) {
  return Sup(S.A);
}

inline srmatrix InfRe(const scimatrix_slice& S) {
  return InfRe(S.A);
}

inline srmatrix InfIm(const scimatrix_slice& S) {
  return InfIm(S.A);
}

inline srmatrix SupRe(const scimatrix_slice& S) {
  return SupRe(S.A);
}

inline srmatrix SupIm(const scimatrix_slice& S) {
  return SupIm(S.A);
}

inline scimatrix::scimatrix(const srmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scimatrix::scimatrix(const scmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scimatrix::scimatrix(const simatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scimatrix::scimatrix(const scimatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline scimatrix operator-(const scimatrix_slice& M) {
  return sp_m_negative<scimatrix,scimatrix>(M.A);
}

inline scimatrix operator+(const scimatrix_slice& M) {
  return M.A;
}

inline scimatrix operator*(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,srmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<simatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const simatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2.A);
}

inline scimatrix operator*(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_mult<scimatrix,srmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_mult<scimatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<scimatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_mult<scimatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_mult<srmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_mult<simatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_mult<scmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const simatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_mult<simatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<scmatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1.A,M2);
}

inline scimatrix operator*(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,srmatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<scimatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_mult<simatrix,scimatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const simatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,scmatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline scimatrix operator*(const scmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<scmatrix,simatrix,scimatrix,sparse_cidot,cinterval>(M1,M2.A);
}

inline cimatrix operator*(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_mult<scimatrix,rmatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<scimatrix,imatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_mult<scimatrix,cmatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_mult<scimatrix,cimatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_mult<srmatrix,cimatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_mult<simatrix,cimatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_mult<scmatrix,cimatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const simatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_mult<simatrix,cmatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<scmatrix,imatrix,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,srmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,scmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,simatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<imatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const imatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,scmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,simatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_mult<scimatrix,rmatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_mult<scimatrix,cmatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<scimatrix,imatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_mult<scimatrix,cimatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_mult<srmatrix,cimatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_mult<simatrix,cimatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_mult<scmatrix,cimatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const simatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_mult<simatrix,cmatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const scmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<scmatrix,imatrix_slice,cimatrix,sparse_cidot>(M1.A,M2);
}

inline cimatrix operator*(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,srmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,scmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,simatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<cimatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<imatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,scimatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const imatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,scmatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline cimatrix operator*(const cmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<cmatrix,simatrix,cimatrix,sparse_cidot>(M1,M2.A);
}

inline scivector operator*(const scimatrix_slice& M, const srvector& v) {
  return spsp_mv_mult<scimatrix,srvector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<scimatrix,sivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const scvector& v) {
  return spsp_mv_mult<scimatrix,scvector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const scivector& v) {
  return spsp_mv_mult<scimatrix,scivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const srmatrix_slice& M, const scivector& v) {
  return spsp_mv_mult<srmatrix,scivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const simatrix_slice& M, const scivector& v) {
  return spsp_mv_mult<simatrix,scivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scmatrix_slice& M, const scivector& v) {
  return spsp_mv_mult<scmatrix,scivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const simatrix_slice& M, const scvector& v) {
  return spsp_mv_mult<simatrix,scvector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scmatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<scmatrix,sivector,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const srvector_slice& v) {
  return spsl_mv_mult<scimatrix,srvector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<scimatrix,sivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const scvector_slice& v) {
  return spsl_mv_mult<scimatrix,scvector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scimatrix_slice& M, const scivector_slice& v) {
  return spsl_mv_mult<scimatrix,scivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const srmatrix_slice& M, const scivector_slice& v) {
  return spsl_mv_mult<srmatrix,scivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scmatrix_slice& M, const scivector_slice& v) {
  return spsl_mv_mult<scmatrix,scivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const simatrix_slice& M, const scivector_slice& v) {
  return spsl_mv_mult<simatrix,scivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const simatrix_slice& M, const scvector_slice& v) {
  return spsl_mv_mult<simatrix,scvector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline scivector operator*(const scmatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<scmatrix,sivector_slice,scivector,sparse_cidot,cinterval>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const rvector& v) {
  return spf_mv_mult<scimatrix,rvector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const ivector& v) {
  return spf_mv_mult<scimatrix,ivector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const cvector& v) {
  return spf_mv_mult<scimatrix,cvector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const civector& v) {
  return spf_mv_mult<scimatrix,civector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const srmatrix_slice& M, const civector& v) {
  return spf_mv_mult<srmatrix,civector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const simatrix_slice& M, const civector& v) {
  return spf_mv_mult<simatrix,civector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scmatrix_slice& M, const civector& v) {
  return spf_mv_mult<scmatrix,civector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const simatrix_slice& M, const cvector& v) {
  return spf_mv_mult<simatrix,cvector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scmatrix_slice& M, const ivector& v) {
  return spf_mv_mult<scmatrix,ivector,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const rvector_slice& v) {
  return spf_mv_mult<scimatrix,rvector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<scimatrix,ivector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const cvector_slice& v) {
  return spf_mv_mult<scimatrix,cvector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scimatrix_slice& M, const civector_slice& v) {
  return spf_mv_mult<scimatrix,civector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const srmatrix_slice& M, const civector_slice& v) {
  return spf_mv_mult<srmatrix,civector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scmatrix_slice& M, const civector_slice& v) {
  return spf_mv_mult<scmatrix,civector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const simatrix_slice& M, const civector_slice& v) {
  return spf_mv_mult<simatrix,civector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const simatrix_slice& M, const cvector_slice& v) {
  return spf_mv_mult<simatrix,cvector_slice,civector,sparse_cidot>(M.A,v);
}

inline civector operator*(const scmatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<scmatrix,ivector_slice,civector,sparse_cidot>(M.A,v);
}

inline scimatrix operator/(const scimatrix_slice& M, const real& r) {
  return sp_ms_div<scimatrix,real,scimatrix>(M.A,r);
}

inline scimatrix operator/(const scimatrix_slice& M, const complex& r) {
  return sp_ms_div<scimatrix,complex,scimatrix>(M.A,r);
}

inline scimatrix operator/(const scimatrix_slice& M, const interval& r) {
  return sp_ms_div<scimatrix,interval,scimatrix>(M.A,r);
}

inline scimatrix operator/(const scimatrix_slice& M, const cinterval& r) {
  return sp_ms_div<scimatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator/(const srmatrix_slice& M, const cinterval& r) {
  return sp_ms_div<srmatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator/(const simatrix_slice& M, const cinterval& r) {
  return sp_ms_div<simatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator/(const scmatrix_slice& M, const cinterval& r) {
  return sp_ms_div<scmatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator/(const simatrix_slice& M, const complex& r) {
  return sp_ms_div<simatrix,complex,scimatrix>(M.A,r);
}

inline scimatrix operator/(const scmatrix_slice& M, const interval& r) {
  return sp_ms_div<scmatrix,interval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scimatrix_slice& M, const real& r) {
  return sp_ms_mult<scimatrix,real,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scimatrix_slice& M, const complex& r) {
  return sp_ms_mult<scimatrix,complex,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scimatrix_slice& M, const interval& r) {
  return sp_ms_mult<scimatrix,interval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scimatrix_slice& M, const cinterval& r) {
  return sp_ms_mult<scimatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const srmatrix_slice& M, const cinterval& r) {
  return sp_ms_mult<srmatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const simatrix_slice& M, const cinterval& r) {
  return sp_ms_mult<simatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scmatrix_slice& M, const cinterval& r) {
  return sp_ms_mult<scmatrix,cinterval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const simatrix_slice& M, const complex& r) {
  return sp_ms_mult<simatrix,complex,scimatrix>(M.A,r);
}

inline scimatrix operator*(const scmatrix_slice& M, const interval& r) {
  return sp_ms_mult<scmatrix,interval,scimatrix>(M.A,r);
}

inline scimatrix operator*(const real& r, const scimatrix_slice& M) {
  return sp_sm_mult<real,scimatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const complex& r, const scimatrix_slice& M) {
  return sp_sm_mult<complex,scimatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const interval& r, const scimatrix_slice& M) {
  return sp_sm_mult<interval,scimatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const cinterval& r, const scimatrix_slice& M) {
  return sp_sm_mult<cinterval,scimatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const cinterval& r, const srmatrix_slice& M) {
  return sp_sm_mult<cinterval,srmatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const cinterval& r, const simatrix_slice& M) {
  return sp_sm_mult<cinterval,simatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const cinterval& r, const scmatrix_slice& M) {
  return sp_sm_mult<cinterval,scmatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const complex& r, const simatrix_slice& M) {
  return sp_sm_mult<complex,simatrix,scimatrix>(r,M.A);
}

inline scimatrix operator*(const interval& r, const scmatrix_slice& M) {
  return sp_sm_mult<interval,scmatrix,scimatrix>(r,M.A);
}

inline scimatrix operator+(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const simatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator+(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_add<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_add<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_add<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_add<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_add<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_add<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const simatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_add<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator+(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<scimatrix,srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<scimatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<scimatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<scimatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<srmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<simatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_add<scmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const simatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_add<simatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator+(const scmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<scmatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline cimatrix operator+(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_add<scimatrix,rmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<scimatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_add<scimatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_add<scimatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_add<srmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_add<simatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_add<scmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const simatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_add<simatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<scmatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<cimatrix,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<cimatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<cimatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<cimatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<rmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<imatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<cmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const imatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<imatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<cmatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_add<scimatrix,rmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_add<scimatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<scimatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_add<scimatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_add<srmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_add<simatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_add<scmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const simatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_add<simatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const scmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<scmatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator+(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<cimatrix_slice,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<cimatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<cimatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<cimatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<rmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_add<cmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const imatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator+(const cmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<cmatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline scimatrix operator-(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const simatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator-(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_sub<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_sub<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_sub<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_sub<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_sub<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_sub<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const simatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_sub<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator-(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<scimatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<simatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const simatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator-(const scmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<scmatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline cimatrix operator-(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_sub<scimatrix,rmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<scimatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_sub<scimatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_sub<scimatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_sub<srmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_sub<simatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_sub<scmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const simatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_sub<simatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<scmatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<cimatrix,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<cimatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<cimatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<cimatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<rmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<imatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<cmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const imatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<imatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<cmatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_sub<scimatrix,rmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_sub<scimatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<scimatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_sub<scimatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_sub<srmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_sub<simatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_sub<scmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const simatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_sub<simatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const scmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<scmatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator-(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<cimatrix_slice,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<cimatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<cimatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<cimatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<rmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_sub<cmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const imatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator-(const cmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<cmatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline scimatrix operator|(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const simatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<scimatrix,srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_hull<scimatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_hull<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_hull<srmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_hull<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_hull<scmatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const simatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_hull<simatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<scmatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<scimatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<simatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const simatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline cimatrix operator|(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<scimatrix,rmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<scimatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_hull<scimatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_hull<scimatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_hull<srmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_hull<simatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_hull<scmatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const simatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_hull<simatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<scmatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<cimatrix,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<cimatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<cimatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<cimatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<imatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<cmatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const imatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<imatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<cmatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<scimatrix,rmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_hull<scimatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<scimatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_hull<scimatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_hull<srmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_hull<simatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_hull<scmatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const simatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_hull<simatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<scmatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<cimatrix_slice,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<cimatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<cimatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<cimatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_hull<cmatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const imatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<cmatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline scimatrix operator|(const scmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const srmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scmatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator|(const scmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<scmatrix,srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const srmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_hull<srmatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scmatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_hull<scmatrix,scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator|(const scmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const srmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator|(const scmatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_hull<scmatrix,scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline cimatrix operator|(const scmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<scmatrix,rmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const srmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_hull<srmatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_hull<scmatrix,cmatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const cmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<cmatrix,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const rmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<cmatrix,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const scmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<scmatrix,rmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const srmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_hull<srmatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const scmatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_hull<scmatrix,cmatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator|(const cmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<cmatrix_slice,srmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const rmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator|(const cmatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_hull<cmatrix_slice,scmatrix,cimatrix>(M1,M2.A);
}

inline scimatrix operator&(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator&(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_intersect<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator&(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline scimatrix operator&(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_intersect<scimatrix,simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator&(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_intersect<scimatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator&(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_intersect<simatrix,scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline scimatrix operator&(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<scimatrix,simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator&(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_intersect<scimatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline scimatrix operator&(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline cimatrix operator&(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_intersect<scimatrix,imatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_intersect<scimatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_intersect<simatrix,cimatrix,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<cimatrix,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator&(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_intersect<cimatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator&(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_intersect<imatrix,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator&(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_intersect<scimatrix,imatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_intersect<scimatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_intersect<simatrix,cimatrix_slice,cimatrix>(M1.A,M2);
}

inline cimatrix operator&(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<cimatrix_slice,simatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator&(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_intersect<cimatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix operator&(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_intersect<imatrix_slice,scimatrix,cimatrix>(M1,M2.A);
}

inline cimatrix& cimatrix::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline cimatrix& cimatrix::operator=(const simatrix_slice& M) {
  *this = imatrix(M);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const simatrix_slice& M) {
  *this = imatrix(M);
  return *this;
}

inline cimatrix& cimatrix::operator=(const scmatrix_slice& M) {
  *this = cmatrix(M);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const scmatrix_slice& M) {
  *this = cmatrix(M);
  return *this;
}

inline cimatrix& cimatrix::operator=(const scimatrix_slice& M) {
  *this = cimatrix(M);
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator=(const scimatrix_slice& M) {
  *this = cimatrix(M);
  return *this;
}

inline cimatrix& cimatrix::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix& cimatrix::operator+=(const scmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator+=(const scmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix& cimatrix::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix& cimatrix::operator+=(const scimatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator+=(const scimatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline cimatrix& cimatrix::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator-=(const scmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator-=(const scmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator-=(const scimatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator-=(const scimatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator*=(const scmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator*=(const scmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator*=(const scimatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator*=(const scimatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator|=(const scmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator|=(const scmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator|=(const scimatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator|=(const scimatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline cimatrix& cimatrix::operator&=(const scimatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline cimatrix_slice& cimatrix_slice::operator&=(const scimatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline bool operator==(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const srmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const scmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const simatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const scimatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const scimatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const scimatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const scimatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scimatrix& M1, const scmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scimatrix& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scimatrix& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const srmatrix& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scmatrix& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const simatrix& M1, const scimatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const rmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const cmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const imatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const cimatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const cimatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const cimatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const cimatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const cimatrix& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const imatrix& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cmatrix& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator<(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_less<scimatrix,simatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_less<scimatrix,simatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_less<scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_less<srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_less<scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_less<simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_less<scimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_less<simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_less<scimatrix,imatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_less<scimatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_less<srmatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_less<scmatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_less<simatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_less<cimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<cimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<rmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<imatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<cmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_less<cimatrix_slice,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<cimatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<rmatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<imatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_less<cmatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_less<scimatrix,imatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_less<scimatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_less<srmatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_less<simatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_less<scmatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<=(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<scimatrix,simatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<=(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<=(const scmatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<scmatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<=(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator<=(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_leq<scimatrix,simatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_leq<scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const srmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_leq<srmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const scmatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_leq<scmatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_leq<simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<scimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const srmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const scmatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<scmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_leq<simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_leq<scimatrix,imatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_leq<scimatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const srmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_leq<srmatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const scmatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_leq<scmatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_leq<simatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator<=(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<cimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<cimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const rmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<rmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<imatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const cmatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<cmatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<cimatrix_slice,simatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<cimatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const rmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<rmatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<imatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const cmatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_leq<cmatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator<=(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_leq<scimatrix,imatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<=(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_leq<scimatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<=(const srmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_leq<srmatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<=(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_leq<simatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator<=(const scmatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_leq<scmatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,srmatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,scmatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,simatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_greater<simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_greater<scimatrix,srmatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_greater<scimatrix,scmatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_greater<scimatrix,simatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_greater<scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_greater<simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_greater<scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_greater<simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_greater<scimatrix,rmatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_greater<scimatrix,cmatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_greater<scimatrix,imatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_greater<scimatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_greater<simatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator>(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<cimatrix,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<cimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_greater<cimatrix,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_greater<cimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_greater<imatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<cimatrix_slice,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<cimatrix_slice,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_greater<cimatrix_slice,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_greater<cimatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_greater<imatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_greater<scimatrix,rmatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_greater<scimatrix,cmatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_greater<scimatrix,imatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_greater<scimatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_greater<simatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,srmatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const scmatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,scmatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,simatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>=(const simatrix_slice& M1, const scimatrix_slice& M2) {
  return spsp_mm_geq<simatrix,scimatrix,cinterval>(M1.A,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_geq<scimatrix,srmatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const scmatrix& M2) {
  return spsp_mm_geq<scimatrix,scmatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_geq<scimatrix,simatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_geq<scimatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const scimatrix& M2) {
  return spsp_mm_geq<simatrix,scimatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const scimatrix& M1, const scmatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const scimatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const scimatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_geq<scimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const simatrix& M1, const scimatrix_slice& M2) {
  return spsp_mm_geq<simatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_geq<scimatrix,rmatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const cmatrix& M2) {
  return spf_mm_geq<scimatrix,cmatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const imatrix& M2) {
  return spf_mm_geq<scimatrix,imatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_geq<scimatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const cimatrix& M2) {
  return spf_mm_geq<simatrix,cimatrix,cinterval>(M1.A,M2);
}

inline bool operator>=(const cimatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<cimatrix,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<cimatrix,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix& M1, const scmatrix_slice& M2) {
  return fsp_mm_geq<cimatrix,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_geq<cimatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const imatrix& M1, const scimatrix_slice& M2) {
  return fsp_mm_geq<imatrix,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<cimatrix_slice,srmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<cimatrix_slice,simatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix_slice& M1, const scmatrix_slice& M2) {
  return fsp_mm_geq<cimatrix_slice,scmatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const cimatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_geq<cimatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const imatrix_slice& M1, const scimatrix_slice& M2) {
  return fsp_mm_geq<imatrix_slice,scimatrix,cinterval>(M1,M2.A);
}

inline bool operator>=(const scimatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_geq<scimatrix,rmatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const cmatrix_slice& M2) {
  return spf_mm_geq<scimatrix,cmatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_geq<scimatrix,imatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>=(const scimatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_geq<scimatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator>=(const simatrix_slice& M1, const cimatrix_slice& M2) {
  return spf_mm_geq<simatrix,cimatrix_slice,cinterval>(M1.A,M2);
}

inline bool operator!(const scimatrix_slice& M) {
  return sp_m_not(M.A);
}

inline std::ostream& operator<<(std::ostream& os, const scimatrix_slice& M) {
  return sp_m_output<scimatrix,cinterval>(os, M.A);
}

inline std::istream& operator>>(std::istream& is, scimatrix_slice& M) {
  scimatrix tmp(M.A.m,M.A.n);
  is >> tmp;
  M = tmp;
  return is;
}



class scimatrix_subv {
  private:
    scimatrix_slice dat;
    bool row;
    int index;

    scimatrix_subv(scimatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

    scimatrix_subv(const scimatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

  public:
    cinterval& operator[](const int i) {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_subv::operator[](int)"));
#endif
        return dat.element(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_subv::operator[](int)"));
#endif
        return dat.element(i,index);
      }
    }

    const cinterval operator[](const int i) const {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_subv::operator[](int)"));
#endif
        return dat(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("scimatrix_subv::operator[](int)"));
#endif
        return dat(i,index);
      }
    }

    scimatrix_subv& operator=(const real& v) {
      return sv_vs_assign(*this,v);
    }

    scimatrix_subv& operator=(const complex& v) {
      return sv_vs_assign(*this,v);
    }

    scimatrix_subv& operator=(const interval& v) {
      return sv_vs_assign(*this,v);
    }

    scimatrix_subv& operator=(const cinterval& v) {
      return sv_vs_assign(*this,v);
    }

    scimatrix_subv& operator=(const srvector& v) {
      return svsp_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const scvector& v) {
      return svsp_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const sivector& v) {
      return svsp_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const scivector& v) {
      return svsp_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const srvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const scvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const sivector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const scivector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const rvector& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const cvector& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const ivector& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const civector& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const rvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const cvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const ivector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const civector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    scimatrix_subv& operator=(const srmatrix_subv& v) {
      return svsp_vv_assign(*this,srvector(v));
    }

    scimatrix_subv& operator=(const simatrix_subv& v) {
      return svsp_vv_assign(*this,sivector(v));
    }

    scimatrix_subv& operator=(const scmatrix_subv& v) {
      return svsp_vv_assign(*this,scvector(v));
    }

    scimatrix_subv& operator=(const scimatrix_subv& v) {
      return svsp_vv_assign(*this,scivector(v));
    }

    scimatrix_subv& operator*=(const real&);
    scimatrix_subv& operator*=(const complex&);
    scimatrix_subv& operator*=(const interval&);
    scimatrix_subv& operator*=(const cinterval&);
    scimatrix_subv& operator/=(const real&);
    scimatrix_subv& operator/=(const complex&);
    scimatrix_subv& operator/=(const interval&);
    scimatrix_subv& operator/=(const cinterval&);
    scimatrix_subv& operator+=(const srvector&);
    scimatrix_subv& operator+=(const srvector_slice&);
    scimatrix_subv& operator+=(const rvector&);
    scimatrix_subv& operator+=(const rvector_slice&);
    scimatrix_subv& operator-=(const srvector&);
    scimatrix_subv& operator-=(const srvector_slice&);
    scimatrix_subv& operator-=(const rvector&);
    scimatrix_subv& operator-=(const rvector_slice&);
    scimatrix_subv& operator+=(const scvector&);
    scimatrix_subv& operator+=(const scvector_slice&);
    scimatrix_subv& operator+=(const cvector&);
    scimatrix_subv& operator+=(const cvector_slice&);
    scimatrix_subv& operator-=(const scvector&);
    scimatrix_subv& operator-=(const scvector_slice&);
    scimatrix_subv& operator-=(const cvector&);
    scimatrix_subv& operator-=(const cvector_slice&);
    scimatrix_subv& operator+=(const sivector&);
    scimatrix_subv& operator+=(const sivector_slice&);
    scimatrix_subv& operator+=(const ivector&);
    scimatrix_subv& operator+=(const ivector_slice&);
    scimatrix_subv& operator-=(const sivector&);
    scimatrix_subv& operator-=(const sivector_slice&);
    scimatrix_subv& operator-=(const ivector&);
    scimatrix_subv& operator-=(const ivector_slice&);
    scimatrix_subv& operator+=(const scivector&);
    scimatrix_subv& operator+=(const scivector_slice&);
    scimatrix_subv& operator+=(const civector&);
    scimatrix_subv& operator+=(const civector_slice&);
    scimatrix_subv& operator-=(const scivector&);
    scimatrix_subv& operator-=(const scivector_slice&);
    scimatrix_subv& operator-=(const civector&);
    scimatrix_subv& operator-=(const civector_slice&);
    scimatrix_subv& operator|=(const srvector&);
    scimatrix_subv& operator|=(const srvector_slice&);
    scimatrix_subv& operator|=(const rvector&);
    scimatrix_subv& operator|=(const rvector_slice&);
    scimatrix_subv& operator|=(const scvector&);
    scimatrix_subv& operator|=(const scvector_slice&);
    scimatrix_subv& operator|=(const cvector&);
    scimatrix_subv& operator|=(const cvector_slice&);
    scimatrix_subv& operator|=(const sivector&);
    scimatrix_subv& operator|=(const sivector_slice&);
    scimatrix_subv& operator|=(const ivector&);
    scimatrix_subv& operator|=(const ivector_slice&);
    scimatrix_subv& operator|=(const scivector&);
    scimatrix_subv& operator|=(const scivector_slice&);
    scimatrix_subv& operator|=(const civector&);
    scimatrix_subv& operator|=(const civector_slice&);

    friend scivector operator-(const scimatrix_subv&);

    friend std::istream& operator>>(std::istream&, scimatrix_subv&);

    friend int Lb(const scimatrix_subv&);
    friend int Ub(const scimatrix_subv&);
    friend int VecLen(const scimatrix_subv&);
    friend sivector Re(const scimatrix_subv&);
    friend sivector Im(const scimatrix_subv&);
    friend scvector Inf(const scimatrix_subv&);
    friend scvector Sup(const scimatrix_subv&);
    friend srvector InfRe(const scimatrix_subv&);
    friend srvector InfIm(const scimatrix_subv&);
    friend srvector SupRe(const scimatrix_subv&);
    friend srvector SupIm(const scimatrix_subv&);

    friend class srvector;
    friend class srmatrix;
    friend class srmatrix_slice;
    friend class scvector;
    friend class scmatrix;
    friend class scmatrix_slice;
    friend class sivector;
    friend class simatrix;
    friend class simatrix_slice;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;

#include "vector_friend_declarations.inl"
};

inline int Lb(const scimatrix_subv& S) {
  if(S.row)
    return Lb(S.dat, 2);
  else
    return Lb(S.dat, 1);
}

inline int Ub(const scimatrix_subv& S) {
  if(S.row)
    return Ub(S.dat, 2);
  else
    return Ub(S.dat, 1);
}

inline int VecLen(const scimatrix_subv& S) {
  return Ub(S)-Lb(S)+1;
}

inline sivector Re(const scimatrix_subv& S) {
  return Re(scivector(S));
}

inline sivector Im(const scimatrix_subv& S) {
  return Im(scivector(S));
}

inline scivector conj(const scimatrix_subv& S) {
  return conj(scivector(S));
}

inline sivector abs(const scimatrix_subv& S) {
  return abs(scivector(S));
}

inline scvector mid(const scimatrix_subv& S) {
  return mid(scivector(S));
}

inline scvector diam(const scimatrix_subv& S) {
  return diam(scivector(S));
}

inline scvector Inf(const scimatrix_subv& S) {
  return Inf(scivector(S));
}

inline scvector Sup(const scimatrix_subv& S) {
  return Sup(scivector(S));
}

inline srvector InfRe(const scimatrix_subv& S) {
  return InfRe(scivector(S));
}

inline srvector InfIm(const scimatrix_subv& S) {
  return InfIm(scivector(S));
}

inline srvector SupRe(const scimatrix_subv& S) {
  return SupRe(scivector(S));
}

inline srvector SupIm(const scimatrix_subv& S) {
  return SupIm(scivector(S));
}

inline std::ostream& operator<<(std::ostream& os, const scimatrix_subv& v) {
  os << scivector(v);
  return os;
}

inline std::istream& operator>>(std::istream& is, scimatrix_subv& v) {
  int n=0;
  if(v.row) n=v.dat.A.n; else n=v.dat.A.m;
  scivector tmp(n);
  is >> tmp;
  v = tmp;
  return is;
}

inline scimatrix_subv scimatrix::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix::operator[](const cxscmatrix_column&)"));
#endif
  return scimatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline scimatrix_subv scimatrix::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix::operator[](const int)"));
#endif
  return scimatrix_subv(*this, true, i, i, lb2, ub2);
}

inline const scimatrix_subv scimatrix::operator[](const cxscmatrix_column& c) const{
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix::operator[](const cxscmatrix_column&)"));
#endif
  return scimatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline const scimatrix_subv scimatrix::operator[](const int i) const{
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix::operator[](const int)"));
#endif
  return scimatrix_subv(*this, true, i, i, lb2, ub2);
}

inline scimatrix_subv scimatrix_slice::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix_slice::operator[](const int)"));
#endif
  return scimatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline scimatrix_subv scimatrix_slice::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return scimatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline const scimatrix_subv scimatrix_slice::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix_slice::operator[](const int)"));
#endif
  return scimatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline const scimatrix_subv scimatrix_slice::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("scimatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return scimatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline scivector::scivector(const scimatrix_subv& A) {
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

inline scivector::scivector(const srmatrix_subv& A) {
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
        x.push_back(cinterval(A.dat.A.x[k]));
      }
    }

  } else {
    lb = A.dat.A.lb1;
    ub = A.dat.A.ub1;
    n = ub-lb+1; 

    for(unsigned int k=0 ; k<A.dat.A.ind.size() ; k++) {
        p.push_back(A.dat.A.ind[k]);
        x.push_back(cinterval(A.dat.A.x[k]));
    }
  }
}

inline scivector::scivector(const simatrix_subv& A) {
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
        x.push_back(cinterval(A.dat.A.x[k]));
      }
    }

  } else {
    lb = A.dat.A.lb1;
    ub = A.dat.A.ub1;
    n = ub-lb+1; 

    for(unsigned int k=0 ; k<A.dat.A.ind.size() ; k++) {
        p.push_back(A.dat.A.ind[k]);
        x.push_back(cinterval(A.dat.A.x[k]));
    }
  }
}

inline scivector::scivector(const scmatrix_subv& A) {
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
        x.push_back(cinterval(A.dat.A.x[k]));
      }
    }

  } else {
    lb = A.dat.A.lb1;
    ub = A.dat.A.ub1;
    n = ub-lb+1; 

    for(unsigned int k=0 ; k<A.dat.A.ind.size() ; k++) {
        p.push_back(A.dat.A.ind[k]);
        x.push_back(cinterval(A.dat.A.x[k]));
    }
  }
}


inline scivector operator-(const scimatrix_subv& v) {
 scivector s(v);
 return -s;
}


inline scivector operator/(const scimatrix_subv& v1, const real& v2) {
  return scivector(v1) / v2;
}

inline scivector operator/(const scimatrix_subv& v1, const complex& v2) {
  return scivector(v1) / v2;
}

inline scivector operator/(const scimatrix_subv& v1, const interval& v2) {
  return scivector(v1) / v2;
}

inline scivector operator/(const scimatrix_subv& v1, const cinterval& v2) {
  return scivector(v1) / v2;
}

inline scivector operator/(const srmatrix_subv& v1, const cinterval& v2) {
  return srvector(v1) / v2;
}

inline scivector operator/(const scmatrix_subv& v1, const cinterval& v2) {
  return scvector(v1) / v2;
}

inline scivector operator/(const simatrix_subv& v1, const cinterval& v2) {
  return sivector(v1) / v2;
}

inline scivector operator/(const simatrix_subv& v1, const complex& v2) {
  return sivector(v1) / v2;
}

inline scivector operator/(const scmatrix_subv& v1, const interval& v2) {
  return scvector(v1) / v2;
}

inline scivector operator*(const scimatrix_subv& v1, const real& v2) {
  return scivector(v1) * v2;
}

inline scivector operator*(const scimatrix_subv& v1, const complex& v2) {
  return scivector(v1) * v2;
}

inline scivector operator*(const scimatrix_subv& v1, const interval& v2) {
  return scivector(v1) * v2;
}

inline scivector operator*(const scimatrix_subv& v1, const cinterval& v2) {
  return scivector(v1) * v2;
}

inline scivector operator*(const srmatrix_subv& v1, const cinterval& v2) {
  return srvector(v1) * v2;
}

inline scivector operator*(const scmatrix_subv& v1, const cinterval& v2) {
  return scvector(v1) * v2;
}

inline scivector operator*(const simatrix_subv& v1, const cinterval& v2) {
  return sivector(v1) * v2;
}

inline scivector operator*(const simatrix_subv& v1, const complex& v2) {
  return sivector(v1) * v2;
}

inline scivector operator*(const scmatrix_subv& v1, const interval& v2) {
  return scvector(v1) * v2;
}

inline scivector operator*(const real& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline scivector operator*(const complex& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline scivector operator*(const interval& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline scivector operator*(const cinterval& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline scivector operator*(const cinterval& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline scivector operator*(const cinterval& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline scivector operator*(const cinterval& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline scivector operator*(const complex& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline scivector operator*(const interval& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const sivector& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const scvector& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const sivector_slice& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const scvector_slice& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const civector& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const ivector& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const cvector& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) * v2;
}

inline cinterval operator*(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const civector_slice& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const civector_slice& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scmatrix_subv& v1, const ivector_slice& v2) {
  return scvector(v1) * v2;
}

inline cinterval operator*(const simatrix_subv& v1, const cvector_slice& v2) {
  return sivector(v1) * v2;
}

inline cinterval operator*(const scivector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline cinterval operator*(const scivector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const scivector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const scivector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const srvector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const scvector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const sivector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const scvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const sivector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline cinterval operator*(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const scvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const sivector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const civector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline cinterval operator*(const civector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const civector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const civector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const rvector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const cvector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const ivector& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const ivector& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const cvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline cinterval operator*(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline cinterval operator*(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 * scivector(v2);
}

inline cinterval operator*(const ivector_slice& v1, const scmatrix_subv& v2) {
  return v1 * scvector(v2);
}

inline cinterval operator*(const cvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

inline scivector operator+(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) + v2;
}

inline scivector operator+(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) + v2;
}

inline scivector operator+(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) + v2;
}

inline scivector operator+(const scmatrix_subv& v1, const sivector& v2) {
  return scvector(v1) + v2;
}

inline scivector operator+(const simatrix_subv& v1, const scvector& v2) {
  return sivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) + v2;
}

inline scivector operator+(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) + v2;
}

inline scivector operator+(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) + v2;
}

inline scivector operator+(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) + v2;
}

inline scivector operator+(const simatrix_subv& v1, const scvector_slice& v2) {
  return sivector(v1) + v2;
}

inline scivector operator+(const scmatrix_subv& v1, const sivector_slice& v2) {
  return scvector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) + v2;
}

inline civector operator+(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) + v2;
}

inline civector operator+(const scmatrix_subv& v1, const civector& v2) {
  return scvector(v1) + v2;
}

inline civector operator+(const scmatrix_subv& v1, const ivector& v2) {
  return scvector(v1) + v2;
}

inline civector operator+(const simatrix_subv& v1, const cvector& v2) {
  return sivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) + v2;
}

inline civector operator+(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) + v2;
}

inline civector operator+(const scmatrix_subv& v1, const civector_slice& v2) {
  return scvector(v1) + v2;
}

inline civector operator+(const simatrix_subv& v1, const civector_slice& v2) {
  return sivector(v1) + v2;
}

inline civector operator+(const simatrix_subv& v1, const cvector_slice& v2) {
  return sivector(v1) + v2;
}

inline civector operator+(const scmatrix_subv& v1, const ivector_slice& v2) {
  return scvector(v1) + v2;
}

inline scivector operator+(const scivector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline scivector operator+(const scivector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scivector operator+(const scivector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline scivector operator+(const scivector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const srvector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const scvector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const sivector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const sivector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scivector operator+(const scvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline scivector operator+(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline scivector operator+(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scivector operator+(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline scivector operator+(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline scivector operator+(const sivector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline scivector operator+(const scvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline civector operator+(const civector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline civector operator+(const civector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline civector operator+(const civector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline civector operator+(const civector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const rvector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const cvector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const ivector& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const ivector& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline civector operator+(const cvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline civector operator+(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline civector operator+(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline civector operator+(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline civector operator+(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 + scivector(v2);
}

inline civector operator+(const ivector_slice& v1, const scmatrix_subv& v2) {
  return v1 + scvector(v2);
}

inline civector operator+(const cvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

inline scivector operator-(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) - v2;
}

inline scivector operator-(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) - v2;
}

inline scivector operator-(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) - v2;
}

inline scivector operator-(const scmatrix_subv& v1, const sivector& v2) {
  return scvector(v1) - v2;
}

inline scivector operator-(const simatrix_subv& v1, const scvector& v2) {
  return sivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) - v2;
}

inline scivector operator-(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) - v2;
}

inline scivector operator-(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) - v2;
}

inline scivector operator-(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) - v2;
}

inline scivector operator-(const simatrix_subv& v1, const scvector_slice& v2) {
  return sivector(v1) - v2;
}

inline scivector operator-(const scmatrix_subv& v1, const sivector_slice& v2) {
  return scvector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) - v2;
}

inline civector operator-(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) - v2;
}

inline civector operator-(const scmatrix_subv& v1, const civector& v2) {
  return scvector(v1) - v2;
}

inline civector operator-(const scmatrix_subv& v1, const ivector& v2) {
  return scvector(v1) - v2;
}

inline civector operator-(const simatrix_subv& v1, const cvector& v2) {
  return sivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) - v2;
}

inline civector operator-(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) - v2;
}

inline civector operator-(const scmatrix_subv& v1, const civector_slice& v2) {
  return scvector(v1) - v2;
}

inline civector operator-(const simatrix_subv& v1, const civector_slice& v2) {
  return sivector(v1) - v2;
}

inline civector operator-(const simatrix_subv& v1, const cvector_slice& v2) {
  return sivector(v1) - v2;
}

inline civector operator-(const scmatrix_subv& v1, const ivector_slice& v2) {
  return scvector(v1) - v2;
}

inline scivector operator-(const scivector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline scivector operator-(const scivector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scivector operator-(const scivector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline scivector operator-(const scivector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const srvector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const scvector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const sivector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const sivector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scivector operator-(const scvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline scivector operator-(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline scivector operator-(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scivector operator-(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline scivector operator-(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline scivector operator-(const sivector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline scivector operator-(const scvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline civector operator-(const civector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline civector operator-(const civector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline civector operator-(const civector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline civector operator-(const civector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const rvector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const cvector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const ivector& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const ivector& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline civector operator-(const cvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline civector operator-(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline civector operator-(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline civector operator-(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline civector operator-(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 - scivector(v2);
}

inline civector operator-(const ivector_slice& v1, const scmatrix_subv& v2) {
  return v1 - scvector(v2);
}

inline civector operator-(const cvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

inline scivector operator|(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const sivector& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const simatrix_subv& v1, const scvector& v2) {
  return sivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) | v2;
}

inline scivector operator|(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) | v2;
}

inline scivector operator|(const simatrix_subv& v1, const scvector_slice& v2) {
  return sivector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const sivector_slice& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) | v2;
}

inline civector operator|(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const civector& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const ivector& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const simatrix_subv& v1, const cvector& v2) {
  return sivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) | v2;
}

inline civector operator|(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const civector_slice& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const simatrix_subv& v1, const civector_slice& v2) {
  return sivector(v1) | v2;
}

inline civector operator|(const simatrix_subv& v1, const cvector_slice& v2) {
  return sivector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const ivector_slice& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const scivector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline scivector operator|(const scivector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scivector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline scivector operator|(const scivector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const srvector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const scvector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const sivector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const sivector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline scivector operator|(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline scivector operator|(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline scivector operator|(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline scivector operator|(const sivector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline civector operator|(const civector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline civector operator|(const civector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const civector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline civector operator|(const civector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const rvector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const cvector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const ivector& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const ivector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline civector operator|(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline civector operator|(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline civector operator|(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 | scivector(v2);
}

inline civector operator|(const ivector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

inline scivector operator|(const scmatrix_subv& v1, const srvector& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const srmatrix_subv& v1, const scvector& v2) {
  return srvector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const scvector& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const srvector_slice& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const srmatrix_subv& v1, const scvector_slice& v2) {
  return srvector(v1) | v2;
}

inline scivector operator|(const scmatrix_subv& v1, const scvector_slice& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const rvector& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const srmatrix_subv& v1, const cvector& v2) {
  return srvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const cvector& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const rvector_slice& v2) {
  return scvector(v1) | v2;
}

inline civector operator|(const srmatrix_subv& v1, const cvector_slice& v2) {
  return srvector(v1) | v2;
}

inline civector operator|(const scmatrix_subv& v1, const cvector_slice& v2) {
  return scvector(v1) | v2;
}

inline scivector operator|(const scvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline scivector operator|(const srvector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scvector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline scivector operator|(const srvector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scivector operator|(const scvector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline civector operator|(const rvector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline civector operator|(const rvector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline civector operator|(const cvector_slice& v1, const scmatrix_subv& v2) {
  return v1 | scvector(v2);
}

inline scimatrix_subv& scimatrix_subv::operator*=(const real& v) {
  *this = *this * v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator*=(const complex& v) {
  *this = *this * v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator*=(const interval& v) {
  *this = *this * v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator*=(const cinterval& v) {
  *this = *this * v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator/=(const real& v) {
  *this = *this / v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator/=(const complex& v) {
  *this = *this / v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator/=(const interval& v) {
  *this = *this / v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator/=(const cinterval& v) {
  *this = *this / v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const srvector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const srvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const rvector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const rvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const srvector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const srvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const rvector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const rvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const scvector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const scvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const cvector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const cvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const scvector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const scvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const cvector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const cvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const sivector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const sivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const ivector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const ivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const sivector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const sivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const ivector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const ivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const scivector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const scivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const civector& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator+=(const civector_slice& v) {
  *this = *this + v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const scivector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const scivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const civector& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator-=(const civector_slice& v) {
  *this = *this - v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const srvector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const srvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const rvector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const rvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const sivector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const sivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const ivector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const ivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const scvector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const scvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const cvector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const cvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const scivector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const scivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const civector& v) {
  *this = *this | v;
  return *this;
}

inline scimatrix_subv& scimatrix_subv::operator|=(const civector_slice& v) {
  *this = *this | v;
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const srmatrix_subv& v) {
  *this += rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scmatrix_subv& v) {
  *this += cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const simatrix_subv& v) {
  *this += ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scimatrix_subv& v) {
  *this += civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const srvector& v) {
  *this += rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scvector& v) {
  *this += cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const sivector& v) {
  *this += ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scivector& v) {
  *this += civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const srvector_slice& v) {
  *this += rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scvector_slice& v) {
  *this += cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const sivector_slice& v) {
  *this += ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator+=(const scivector_slice& v) {
  *this += civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const srmatrix_subv& v) {
  *this -= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scmatrix_subv& v) {
  *this -= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const simatrix_subv& v) {
  *this -= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scimatrix_subv& v) {
  *this -= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const srvector& v) {
  *this -= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scvector& v) {
  *this -= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const sivector& v) {
  *this -= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scivector& v) {
  *this -= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const srvector_slice& v) {
  *this -= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scvector_slice& v) {
  *this -= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const sivector_slice& v) {
  *this -= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator-=(const scivector_slice& v) {
  *this -= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const srmatrix_subv& v) {
  *this |= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scmatrix_subv& v) {
  *this |= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const simatrix_subv& v) {
  *this |= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scimatrix_subv& v) {
  *this |= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const srvector& v) {
  *this |= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scvector& v) {
  *this |= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const sivector& v) {
  *this |= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scivector& v) {
  *this |= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const srvector_slice& v) {
  *this |= rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scvector_slice& v) {
  *this |= cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const sivector_slice& v) {
  *this |= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator|=(const scivector_slice& v) {
  *this |= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scvector& v) {
  *this = cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scivector& v) {
  *this = civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scvector_slice& v) {
  *this = cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scivector_slice& v) {
  *this = civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const srmatrix_subv& v) {
  *this = rvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const simatrix_subv& v) {
  *this = ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scmatrix_subv& v) {
  *this = cvector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator=(const scimatrix_subv& v) {
  *this = civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const simatrix_subv& v) {
  *this &= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const scimatrix_subv& v) {
  *this &= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const sivector& v) {
  *this &= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const scivector& v) {
  *this &= civector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const sivector_slice& v) {
  *this &= ivector(v);
  return *this;
}

inline cimatrix_subv& cimatrix_subv::operator&=(const scivector_slice& v) {
  *this &= civector(v);
  return *this;
}

inline bool operator==(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const civector& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const scmatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) == v2;
}

inline bool operator==(const scivector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const scivector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const scivector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const scivector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const srvector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const scvector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const sivector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const civector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const civector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const civector& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const civector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const rvector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const cvector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const ivector& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

inline bool operator==(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 == scvector(v2);
}

inline bool operator==(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator==(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 == scivector(v2);
}

inline bool operator!=(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const civector& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const scmatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) != v2;
}

inline bool operator!=(const scivector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const scivector& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const scivector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const scivector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const srvector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const scvector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const sivector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const civector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const civector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const civector& v1, const scmatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const civector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const rvector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const cvector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const ivector& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

inline bool operator!=(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 != scvector(v2);
}

inline bool operator!=(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator!=(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 != scivector(v2);
}

inline bool operator<(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) < v2;
}

inline bool operator<(const scmatrix_subv& v1, const civector& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const scmatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) < v2;
}

inline bool operator<(const scivector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const scivector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const srvector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const scvector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const sivector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const civector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const civector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const rvector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const cvector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const ivector& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

inline bool operator<(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 < scivector(v2);
}

inline bool operator<=(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const scivector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const scmatrix_subv& v1, const scivector& v2) {
  return scvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const scivector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const scmatrix_subv& v1, const scivector_slice& v2) {
  return scvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const civector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) <= v2;
}

inline bool operator<=(const scmatrix_subv& v1, const civector& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const civector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const scmatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) <= v2;
}

inline bool operator<=(const scivector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const scivector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const srvector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const scvector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const sivector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const srvector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const scvector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const civector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const civector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const rvector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const cvector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const ivector& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

inline bool operator<=(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const rvector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const cvector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator<=(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 <= scivector(v2);
}

inline bool operator>(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) > v2;
}

inline bool operator>(const scivector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const scivector& v1, const scmatrix_subv& v2) {
  return v1 > scvector(v2);
}

inline bool operator>(const scivector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const scivector& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const sivector& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 > scvector(v2);
}

inline bool operator>(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const civector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const civector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const civector& v1, const scmatrix_subv& v2) {
  return v1 > scvector(v2);
}

inline bool operator>(const civector& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const ivector& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

inline bool operator>(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 > scvector(v2);
}

inline bool operator>(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 > scivector(v2);
}

inline bool operator>=(const scimatrix_subv& v1, const srvector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const scvector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const sivector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const scivector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const scivector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const srvector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const sivector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const scvector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const scivector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const scivector_slice& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const rvector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const cvector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const ivector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const civector& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const civector& v2) {
  return sivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const rvector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const cvector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const ivector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scimatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const simatrix_subv& v1, const civector_slice& v2) {
  return scivector(v1) >= v2;
}

inline bool operator>=(const scivector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const scivector& v1, const scmatrix_subv& v2) {
  return v1 >= scvector(v2);
}

inline bool operator>=(const scivector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const scivector& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const sivector& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const scivector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const scivector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const scivector_slice& v1, const scmatrix_subv& v2) {
  return v1 >= scvector(v2);
}

inline bool operator>=(const scivector_slice& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const sivector_slice& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const civector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const civector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const civector& v1, const scmatrix_subv& v2) {
  return v1 >= scvector(v2);
}

inline bool operator>=(const civector& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const ivector& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const civector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const civector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

inline bool operator>=(const civector_slice& v1, const scmatrix_subv& v2) {
  return v1 >= scvector(v2);
}

inline bool operator>=(const civector_slice& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator>=(const ivector_slice& v1, const scimatrix_subv& v2) {
  return v1 >= scivector(v2);
}

inline bool operator!(const scimatrix_subv& x) {
  return sv_v_not(x);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scimatrix_subv& v2) {
  accumulate(dot,scivector(v1),scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scmatrix_subv& v2) {
  accumulate(dot,scivector(v1),scvector(v2));
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const simatrix_subv& v2) {
  accumulate(dot,scivector(v1),sivector(v2));
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const srmatrix_subv& v2) {
  accumulate(dot,scivector(v1),srvector(v2));
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const scimatrix_subv& v2) {
  accumulate(dot,scvector(v1),scivector(v2));
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const scimatrix_subv& v2) {
  accumulate(dot,sivector(v1),scivector(v2));
}

inline void accumulate(cidotprecision& dot,const srmatrix_subv& v1, const scimatrix_subv& v2) {
  accumulate(dot,srvector(v1),scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const srvector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scvector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const sivector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scivector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const srmatrix_subv& v1, const scivector& v2) {
  accumulate(dot,srvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const scivector& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const scivector& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const sivector& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const scvector& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const srvector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scvector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const sivector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const scivector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const srmatrix_subv& v1, const scivector_slice& v2) {
  accumulate(dot,srvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const scivector_slice& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const scivector_slice& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const sivector_slice& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const scvector_slice& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const rvector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const ivector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const cvector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const civector& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const srmatrix_subv& v1, const civector& v2) {
  accumulate(dot,srvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const civector& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const civector& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const ivector& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const cvector& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const rvector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const ivector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const cvector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scimatrix_subv& v1, const civector_slice& v2) {
  accumulate(dot,scivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const srmatrix_subv& v1, const civector_slice& v2) {
  accumulate(dot,srvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const civector_slice& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const civector_slice& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scmatrix_subv& v1, const ivector_slice& v2) {
  accumulate(dot,scvector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const simatrix_subv& v1, const cvector_slice& v2) {
  accumulate(dot,sivector(v1),v2);
}

inline void accumulate(cidotprecision& dot,const scivector& v1, const srmatrix_subv& v2) {
  accumulate(dot,v1,srvector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const srvector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scvector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const sivector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scvector& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const sivector& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector_slice& v1, const srmatrix_subv& v2) {
  accumulate(dot,v1,srvector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector_slice& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector_slice& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const scivector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const srvector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const sivector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scvector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const scvector_slice& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const sivector_slice& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const civector& v1, const srmatrix_subv& v2) {
  accumulate(dot,v1,srvector(v2));
}

inline void accumulate(cidotprecision& dot,const civector& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const civector& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const civector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const rvector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const cvector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const ivector& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const ivector& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const cvector& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const civector_slice& v1, const srmatrix_subv& v2) {
  accumulate(dot,v1,srvector(v2));
}

inline void accumulate(cidotprecision& dot,const civector_slice& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const civector_slice& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}

inline void accumulate(cidotprecision& dot,const civector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const rvector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const cvector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const ivector_slice& v1, const scimatrix_subv& v2) {
  accumulate(dot,v1,scivector(v2));
}

inline void accumulate(cidotprecision& dot,const ivector_slice& v1, const scmatrix_subv& v2) {
  accumulate(dot,v1,scvector(v2));
}

inline void accumulate(cidotprecision& dot,const cvector_slice& v1, const simatrix_subv& v2) {
  accumulate(dot,v1,sivector(v2));
}


}  //namespace cxsc;

#include "sparsematrix.inl"

#endif 
 
