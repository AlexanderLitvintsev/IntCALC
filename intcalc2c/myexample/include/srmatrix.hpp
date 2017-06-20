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

/* CVS $Id: srmatrix.hpp,v 1.13 2010/12/21 09:54:02 cxsc Exp $ */

#ifndef _CXSC_SRMATRIX_HPP_INCLUDED
#define _CXSC_SRMATRIX_HPP_INCLUDED

#include <real.hpp>
#include <rmatrix.hpp>
#include <srvector.hpp>
#include <cidot.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sparsedot.hpp>
#include <sparsematrix.hpp>

namespace cxsc {

enum STORAGE_TYPE{triplet,compressed_row,compressed_column};

class srmatrix_slice;
class srmatrix_subv;
class scmatrix;
class scmatrix_slice;
class scmatrix_subv;
class simatrix;
class simatrix_slice;
class simatrix_subv;
class sivector_slice;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;
class scivector_slice;


inline bool comp_pair(std::pair<int,real> p1, std::pair<int,real> p2) {
  return p1.first < p2.first;
}


class srmatrix {

  private:
    std::vector<int> p;
    std::vector<int> ind;
    std::vector<real> x;
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

    std::vector<real>& values() {
      return x;
    }

    const std::vector<int>& column_pointers() const {
      return p;
    }

    const std::vector<int>& row_indices() const {
      return ind;
    }

    const std::vector<real>& values() const {
      return x;
    }

    srmatrix() {
      p.push_back(0);
      m = n = 0;
      lb1 = lb2 = ub1 = ub2 = 0;
    }

    srmatrix(const int r, const int c) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(2*(m+n));
      x.reserve(2*(m+n));

      p[0] = 0;
    }

    srmatrix(const int r, const int c, const int e) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(e);
      x.reserve(e);

      p[0] = 0;
    }

    srmatrix(const int m, const int n, const int nnz, const intvector& rows, const intvector& cols, const rvector& values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<real> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<real>(rows[Lb(rows)+k],cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<triplet_store<real> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<real>(j,cols[Lb(cols)+k],values[Lb(values)+k]));
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

         std::vector<std::pair<int,real> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[Lb(cols)+k],values[Lb(values)+k]));
           }

           std::sort(work.begin(),work.end(),comp_pair);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    srmatrix(const int m, const int n, const int nnz, const int* rows, const int* cols, const real* values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<real> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<real>(rows[k],cols[k],values[k]));
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

         std::vector<triplet_store<real> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<real>(j,cols[k],values[k]));
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

         std::vector<std::pair<int,real> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[k],values[k]));
           }

           std::sort(work.begin(),work.end(),comp_pair);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    srmatrix(const rmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
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

    srmatrix(const srmatrix_slice&);

    void full(rmatrix& A) const {
       A = rmatrix(lb1,ub1,lb2,ub2);
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
      std::vector<real> xnew;
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


    srmatrix& operator=(const real& A) {
      return sp_ms_assign<srmatrix,real,real>(*this,A);
    }

    srmatrix& operator=(const rmatrix& A) {
      return spf_mm_assign<srmatrix,rmatrix,real>(*this,A);
    }

    srmatrix& operator=(const rmatrix_slice& A) {
      return spf_mm_assign<srmatrix,rmatrix,real>(*this,A);
    }

/*    srmatrix& operator=(const srmatrix& A) {
      p = A.p;
      ind = A.ind;
      x = A.x;
      return *this;
    }*/

    srmatrix& operator=(const srmatrix_slice&);

    const real operator()(int i, int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator()(int, int)"));
#endif
      real r = 0.0;
      for(int k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  r = x[k];
      }
      return r;
    }

    real& element(int i, int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::element(int, int)"));
#endif
      int k;
      for(k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  return x[k];
      }

      //Nicht gefunden, Element muss angelegt werden, da Schreibzugriff moeglich
      std::vector<int>::iterator ind_it = ind.begin() + k;
      std::vector<real>::iterator x_it  = x.begin() + k;
      ind.insert(ind_it, i-lb1);
      x_it = x.insert(x_it, 0.0);
      for(k=j-lb2+1 ; k<(int)p.size() ; k++)
        p[k]++;

      return *x_it;
    }

    srmatrix_subv operator[](const cxscmatrix_column&);
    srmatrix_subv operator[](const int);
    const srmatrix_subv operator[](const cxscmatrix_column&) const;
    const srmatrix_subv operator[](const int) const;

    srmatrix_slice operator()(const int, const int , const int, const int);

    srmatrix operator()(const intvector& pervec, const intvector& q) {
      srmatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      int nnz=0;
      for(int k=0 ; k<n ; k++) {
        A.p[k] = nnz;

        std::map<int,real> work;
        for(int j=p[q[Lb(q)+k]] ; j<p[q[Lb(q)+k]+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,real>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }

        nnz += work.size();
 
      }

      A.p[n] = nnz;

      return A;
    }

    srmatrix operator()(const intvector& pervec) {
      srmatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      for(int k=0 ; k<n ; k++) {
        A.p[k] = p[k];

        std::map<int,real> work;
        for(int j=p[k] ; j<p[k+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,real>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }
 
      }

      A.p[n] = p[n];

      return A;
    }

    srmatrix operator()(const intmatrix& P, const intmatrix& Q) {
      intvector p = permvec(P);
      intvector q = perminv(permvec(Q));
      return (*this)(p,q);
    }

    srmatrix operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    real density() const {
      return p[n]/((double)m*n);
    }

    int get_nnz() const {
      return p[n];
    }

    srmatrix& operator+=(const rmatrix& B) {
      return spf_mm_addassign<srmatrix,rmatrix,rmatrix>(*this,B);
    }

    srmatrix& operator+=(const rmatrix_slice& B) {
      return spf_mm_addassign<srmatrix,rmatrix_slice,rmatrix>(*this,B);
    }

    srmatrix& operator+=(const srmatrix& B) {
      return spsp_mm_addassign<srmatrix,srmatrix,real>(*this,B);
    }

    srmatrix& operator-=(const rmatrix& B) {
      return spf_mm_subassign<srmatrix,rmatrix,rmatrix>(*this,B);
    }

    srmatrix& operator-=(const rmatrix_slice& B) {
      return spf_mm_subassign<srmatrix,rmatrix_slice,rmatrix>(*this,B);
    }

    srmatrix& operator-=(const srmatrix& B) {
      return spsp_mm_subassign<srmatrix,srmatrix,real>(*this,B);
    }

    srmatrix& operator*=(const rmatrix& B) {
      return spf_mm_multassign<srmatrix,rmatrix,sparse_dot,rmatrix>(*this,B);
    }

    srmatrix& operator*=(const rmatrix_slice& B) {
      return spf_mm_multassign<srmatrix,rmatrix,sparse_dot,rmatrix>(*this,B);
    }

    srmatrix& operator*=(const srmatrix& B) {
      return spsp_mm_multassign<srmatrix,srmatrix,sparse_dot,real>(*this,B);
    }

    srmatrix& operator*=(const real& r) {
      return sp_ms_multassign(*this,r);
    }

    srmatrix& operator/=(const real& r) {
      return sp_ms_divassign(*this,r);
    }

    friend int Lb(const srmatrix&, int);
    friend int Ub(const srmatrix&, int);
    friend void SetLb(srmatrix&, const int, const int);
    friend void SetUb(srmatrix&, const int, const int);
    friend int RowLen(const srmatrix&);
    friend int ColLen(const srmatrix&);
    friend srmatrix Re(const scmatrix&);
    friend srmatrix Im(const scmatrix&);
    friend srmatrix Inf(const simatrix&);
    friend srmatrix Sup(const simatrix&);
    friend srmatrix InfRe(const scimatrix&);
    friend srmatrix SupRe(const scimatrix&);
    friend srmatrix InfIm(const scimatrix&);
    friend srmatrix SupIm(const scimatrix&);
    friend srmatrix mid(const simatrix&);
    friend srmatrix diam(const simatrix&);
    friend srmatrix absmin(const simatrix&);
    friend srmatrix absmax(const simatrix&);
    friend srmatrix abs(const srmatrix&);

    friend srmatrix CompMat(const simatrix&);
    friend srmatrix transp(const srmatrix&);
    friend srmatrix Id(const srmatrix&);

    friend std::istream& operator>>(std::istream&, srmatrix_slice&);
    friend std::istream& operator>>(std::istream&, srmatrix_subv&);

    friend class srmatrix_slice;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix;
    friend class scmatrix_slice;
    friend class scmatrix_subv;
    friend class simatrix;
    friend class simatrix_slice;
    friend class simatrix_subv;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class rmatrix;
    friend class imatrix;
    friend class cmatrix;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"
};

inline rmatrix::rmatrix(const srmatrix& A) {
  dat = new real[A.m*A.n];
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

inline srmatrix Id(const srmatrix& A) {
  srmatrix I(A.m, A.n, (A.m>A.n) ? A.m : A.n);
  I.lb1 = A.lb1; I.lb2 = A.lb2;
  I.ub1 = A.ub1; I.ub2 = A.ub2;

  if(A.m < A.n) {
    for(int i=0 ; i<A.m ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(1.0);
    }
  } else {
    for(int i=0 ; i<A.n ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(1.0);
    }
  }

  return I;
}

inline srmatrix transp(const srmatrix& A) {
  srmatrix B(A.n, A.m, A.get_nnz());

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


inline srmatrix abs(const srmatrix& A) {
  srmatrix ret(A);
  for(unsigned int i=0 ; i<ret.x.size() ; i++) 
    ret.x[i] = abs(ret.x[i]);
  return ret;
}


inline void SetLb(srmatrix& A, const int i, const int j) {
  if(i==1) {
    A.lb1 = j;
    A.ub1 = j + A.m - 1;
  } else if(i==2) {
    A.lb2 = j;
    A.ub2 = j + A.n - 1;
  }
}

inline void SetUb(srmatrix& A, const int i, const int j) {
  if(i==1) {
    A.ub1 = j;
    A.lb1 = j - A.m + 1;
  } else if(i==2) {
    A.ub2 = j;
    A.lb2 = j - A.n + 1;
  }
}

inline int Lb(const srmatrix& A, int i) {
  if(i==1) 
    return A.lb1;
  else if(i==2)
    return A.lb2;
  else
    return 1;
}

inline int Ub(const srmatrix& A, int i) {
  if(i==1) 
    return A.ub1;
  else if(i==2)
    return A.ub2;
  else
    return 1;
}

inline int RowLen(const srmatrix& A) {
  return A.n;
}

inline int ColLen(const srmatrix& A) {
  return A.m;
}

inline void Resize(srmatrix& A) {
  sp_m_resize(A);
}

inline void Resize(srmatrix& A, const int m, const int n) {
  sp_m_resize(A,m,n);
}

inline void Resize(srmatrix& A, const int l1, const int u1, const int l2, const int u2) {
  sp_m_resize(A,l1,u1,l2,u2);
}

inline rmatrix operator*(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_mult<rmatrix,srmatrix,rmatrix,sparse_dot>(A,B);
}

inline rmatrix operator*(const srmatrix& A, const rmatrix& B) {
  return spf_mm_mult<srmatrix,rmatrix,rmatrix,sparse_dot>(A,B);
}

inline rmatrix operator*(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_mult<rmatrix_slice,srmatrix,rmatrix,sparse_dot>(A,B);
}

inline rmatrix operator*(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_mult<srmatrix,rmatrix_slice,rmatrix,sparse_dot>(A,B);
}

inline srmatrix operator*(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_mult<srmatrix,srmatrix,srmatrix,sparse_dot,real>(A,B);
}

inline srmatrix operator*(const srmatrix& A, const real& r) {
  return sp_ms_mult<srmatrix,real,srmatrix>(A,r);
}

inline srmatrix operator/(const srmatrix& A, const real& r) {
  return sp_ms_div<srmatrix,real,srmatrix>(A,r);
}

inline srmatrix operator*(const real& r, const srmatrix& A) {
  return sp_sm_mult<real,srmatrix,srmatrix>(r,A);
}

inline rvector operator*(const srmatrix& A, const rvector& v) {
  return spf_mv_mult<srmatrix,rvector,rvector,sparse_dot>(A,v);
}

inline rvector operator*(const srmatrix& A, const rvector_slice& v) {
  return spf_mv_mult<srmatrix,rvector_slice,rvector,sparse_dot>(A,v);
}

inline srvector operator*(const srmatrix& A, const srvector& v) {
  return spsp_mv_mult<srmatrix,srvector,srvector,sparse_dot,real>(A,v);
}

inline srvector operator*(const srmatrix& A, const srvector_slice& v) {
  return spsl_mv_mult<srmatrix,srvector_slice,srvector,sparse_dot,real>(A,v);
}

inline rvector operator*(const rmatrix& A, const srvector& v) {
  return fsp_mv_mult<rmatrix,srvector,rvector,sparse_dot>(A,v);
}

inline rvector operator*(const rmatrix_slice& A, const srvector& v) {
  return fsp_mv_mult<rmatrix_slice,srvector,rvector,sparse_dot>(A,v);
}

inline rvector operator*(const rmatrix& A, const srvector_slice& v) {
  return fsl_mv_mult<rmatrix,srvector_slice,rvector,sparse_dot>(A,v);
}

inline rvector operator*(const rmatrix_slice& A, const srvector_slice& v) {
  return fsl_mv_mult<rmatrix_slice,srvector_slice,rvector,sparse_dot>(A,v);
}

inline rmatrix operator+(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_add<rmatrix,srmatrix,rmatrix>(A,B);
}

inline rmatrix operator+(const srmatrix& A, const rmatrix& B) {
  return spf_mm_add<srmatrix,rmatrix,rmatrix>(A,B);
}

inline rmatrix operator+(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_add<rmatrix_slice,srmatrix,rmatrix>(A,B);
}

inline rmatrix operator+(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_add<srmatrix,rmatrix_slice,rmatrix>(A,B);
}

inline srmatrix operator+(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_add<srmatrix,srmatrix,srmatrix,real>(A,B);
}

inline rmatrix operator-(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_sub<rmatrix,srmatrix,rmatrix>(A,B);
}

inline rmatrix operator-(const srmatrix& A, const rmatrix& B) {
  return spf_mm_sub<srmatrix,rmatrix,rmatrix>(A,B);
}

inline rmatrix operator-(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_sub<rmatrix_slice,srmatrix,rmatrix>(A,B);
}

inline rmatrix operator-(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_sub<srmatrix,rmatrix_slice,rmatrix>(A,B);
}

inline srmatrix operator-(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_sub<srmatrix,srmatrix,srmatrix,real>(A,B);
}

inline srmatrix operator-(const srmatrix& M) {
  return sp_m_negative<srmatrix,srmatrix>(M);
}

inline srmatrix& operator+(srmatrix& A) {
  return A;
}

inline rmatrix& rmatrix::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline rmatrix_slice& rmatrix_slice::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline rmatrix& rmatrix::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline rmatrix_slice& rmatrix_slice::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline rmatrix& rmatrix::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline rmatrix_slice& rmatrix_slice::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline rmatrix& rmatrix::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<rmatrix,srmatrix,sparse_dot,rmatrix>(*this,B);
}

inline rmatrix_slice& rmatrix_slice::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<rmatrix_slice,srmatrix,sparse_dot,rmatrix>(*this,B);
}

inline bool operator==(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const rmatrix& B) {
  return spf_mm_comp(A,B);
}

inline bool operator==(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

inline bool operator==(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const srmatrix& B) {
  return !spsp_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const rmatrix& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator!=(const rmatrix& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const rmatrix_slice& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

inline bool operator!=(const srmatrix& A, const rmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

inline bool operator<(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_less<srmatrix,srmatrix,real>(A,B);
}

inline bool operator<(const srmatrix& A, const rmatrix& B) {
  return spf_mm_less<srmatrix,rmatrix,real>(A,B);
}

inline bool operator<(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_less<rmatrix,srmatrix,real>(A,B);
}

inline bool operator<(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_less<rmatrix_slice,srmatrix,real>(A,B);
}

inline bool operator<(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_less<srmatrix,rmatrix_slice,real>(A,B);
}

inline bool operator<=(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_leq<srmatrix,srmatrix,real>(A,B);
}

inline bool operator<=(const srmatrix& A, const rmatrix& B) {
  return spf_mm_leq<srmatrix,rmatrix,real>(A,B);
}

inline bool operator<=(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_leq<rmatrix,srmatrix,real>(A,B);
}

inline bool operator<=(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_leq<rmatrix_slice,srmatrix,real>(A,B);
}

inline bool operator<=(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_leq<srmatrix,rmatrix_slice,real>(A,B);
}

inline bool operator>(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_greater<srmatrix,srmatrix,real>(A,B);
}

inline bool operator>(const srmatrix& A, const rmatrix& B) {
  return spf_mm_greater<srmatrix,rmatrix,real>(A,B);
}

inline bool operator>(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_greater<rmatrix,srmatrix,real>(A,B);
}

inline bool operator>(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_greater<rmatrix_slice,srmatrix,real>(A,B);
}

inline bool operator>(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_greater<srmatrix,rmatrix_slice,real>(A,B);
}

inline bool operator>=(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_geq<srmatrix,srmatrix,real>(A,B);
}

inline bool operator>=(const srmatrix& A, const rmatrix& B) {
  return spf_mm_geq<srmatrix,rmatrix,real>(A,B);
}

inline bool operator>=(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_geq<rmatrix,srmatrix,real>(A,B);
}

inline bool operator>=(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_geq<rmatrix_slice,srmatrix,real>(A,B);
}

inline bool operator>=(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_geq<srmatrix,rmatrix_slice,real>(A,B);
}

inline bool operator!(const srmatrix& A) {
  return sp_m_not(A);
}

inline std::ostream& operator<<(std::ostream& os, const srmatrix& A) {
  return sp_m_output<srmatrix,real>(os,A);
}

inline std::istream& operator>>(std::istream& is, srmatrix& A) {
  return sp_m_input<srmatrix,real>(is,A);
}

class srmatrix_slice {
  public:
    srmatrix  A;
    srmatrix* M; //Originalmatrix

  private:
    srmatrix_slice(srmatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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

    srmatrix_slice(const srmatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
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
        M = const_cast<srmatrix*>(&Mat);
    }

  public:
    srmatrix_slice& operator=(const real& C) {
      return sl_ms_assign<srmatrix_slice, real, std::vector<real>::iterator, real>(*this,C);
    }

    srmatrix_slice& operator=(const srmatrix& C) {
      return slsp_mm_assign<srmatrix_slice, srmatrix, std::vector<real>::iterator>(*this,C);
    }

    srmatrix_slice& operator=(const rmatrix& C) {
      return slf_mm_assign<srmatrix_slice, rmatrix, std::vector<real>::iterator, real>(*this,C);
    }

    srmatrix_slice& operator=(const rmatrix_slice& C) {
      return slf_mm_assign<srmatrix_slice, rmatrix_slice, std::vector<real>::iterator, real>(*this,C);
    }

    srmatrix_slice& operator=(const srmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    srmatrix_slice& operator*=(const srmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    srmatrix_slice& operator*=(const srmatrix& M) {
      *this = A*M;
      return *this;
    }

    srmatrix_slice& operator*=(const rmatrix& M) {
      *this = A*M;
      return *this;
    }

    srmatrix_slice& operator*=(const rmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    srmatrix_slice& operator*=(const real& r) {
      *this = A*r;
      return *this;
    }

    srmatrix_slice& operator/=(const real& r) {
      *this = A/r;
      return *this;
    }

    srmatrix_slice& operator+=(const srmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    srmatrix_slice& operator+=(const srmatrix& M) {
      *this = A+M;
      return *this;
    } 

    srmatrix_slice& operator+=(const rmatrix& M) {
      *this = A+M;
      return *this;
    } 

    srmatrix_slice& operator+=(const rmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    srmatrix_slice& operator-=(const srmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    srmatrix_slice& operator-=(const srmatrix& M) {
      *this = A-M;
      return *this;
    } 

    srmatrix_slice& operator-=(const rmatrix& M) {
      *this = A-M;
      return *this;
    } 

    srmatrix_slice& operator-=(const rmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    const real operator()(const int i, const int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix_slice::operator()(int, int)"));
#endif
      real r = A(i,j);
      return r;
    }

    real& element(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix::element(int, int)"));
#endif
      return M->element(i,j);
    }

    srmatrix_subv operator[](const int);
    srmatrix_subv operator[](const cxscmatrix_column&);
    const srmatrix_subv operator[](const int) const;
    const srmatrix_subv operator[](const cxscmatrix_column&) const;

    friend int Lb(const srmatrix_slice&, const int);
    friend int Ub(const srmatrix_slice&, const int);
    friend int RowLen(const srmatrix_slice&);
    friend int ColLen(const srmatrix_slice&);

    friend class srmatrix;
    friend class srmatrix_subv;
    friend class srvector;
    friend class scmatrix;
    friend class scmatrix_slice;
    friend class scmatrix_subv;
    friend class scvector;
    friend class simatrix;
    friend class simatrix_slice;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class scivector;
    friend class rmatrix;
    friend class imatrix;
    friend class cmatrix;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"    
};

inline rmatrix::rmatrix(const srmatrix_slice& A) {
  dat = new real[A.A.m*A.A.n];
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

inline int Lb(const srmatrix_slice& S, const int i) {
  return Lb(S.A, i);
}

inline int Ub(const srmatrix_slice& S, const int i) {
  return Ub(S.A, i);
}

inline int RowLen(const srmatrix_slice& S) {
  return RowLen(S.A);
}

inline int ColLen(const srmatrix_slice& S) {
  return ColLen(S.A);
}

inline srmatrix_slice srmatrix::operator()(const int i, const int j, const int k, const int l) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator()(int, int)"));
#endif
  return srmatrix_slice(*this, i, j, k, l);
}

inline srmatrix& srmatrix::operator=(const srmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline srmatrix::srmatrix(const srmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline srmatrix operator-(const srmatrix_slice& M) {
  return sp_m_negative<srmatrix,srmatrix>(M.A);
}

inline srmatrix operator+(const srmatrix_slice& M) {
  return M.A;
}

inline srmatrix operator*(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,srmatrix,srmatrix,sparse_dot,real>(M1.A,M2.A);
}

inline srmatrix operator*(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_mult<srmatrix,srmatrix,srmatrix,sparse_dot,real>(M1.A,M2);
}

inline srmatrix operator*(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,srmatrix,srmatrix,sparse_dot,real>(M1,M2.A);
}

inline rmatrix operator*(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_mult<srmatrix,rmatrix,rmatrix,sparse_dot>(M1.A,M2);
}

inline rmatrix operator*(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,srmatrix,rmatrix,sparse_dot>(M1,M2.A);
}

inline rmatrix operator*(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_mult<srmatrix,rmatrix_slice,rmatrix,sparse_dot>(M1.A,M2);
}

inline rmatrix operator*(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,srmatrix,rmatrix,sparse_dot>(M1,M2.A);
}

inline srvector operator*(const srmatrix_slice& M, const srvector& v) {
  return spsp_mv_mult<srmatrix,srvector,srvector,sparse_dot,real>(M.A,v);
}

inline srvector operator*(const srmatrix_slice& M, const srvector_slice& v) {
  return spsl_mv_mult<srmatrix,srvector_slice,srvector,sparse_dot,real>(M.A,v);
}

inline rvector operator*(const srmatrix_slice& M, const rvector& v) {
  return spf_mv_mult<srmatrix,rvector,rvector,sparse_dot>(M.A,v);
}

inline rvector operator*(const srmatrix_slice& M, const rvector_slice& v) {
  return spf_mv_mult<srmatrix,rvector_slice,rvector,sparse_dot>(M.A,v);
}

inline srmatrix operator*(const srmatrix_slice& M, const real& r) {
  return sp_ms_mult<srmatrix,real,srmatrix>(M.A,r);
}

inline srmatrix operator/(const srmatrix_slice& M, const real& r) {
  return sp_ms_div<srmatrix,real,srmatrix>(M.A,r);
}

inline srmatrix operator*(const real& r, const srmatrix_slice& M) {
  return sp_sm_mult<real,srmatrix,srmatrix>(r,M.A);
}

inline srmatrix operator+(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<srmatrix,srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline srmatrix operator+(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_add<srmatrix,srmatrix,srmatrix,real>(M1.A,M2);
}

inline srmatrix operator+(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<srmatrix,srmatrix,srmatrix,real>(M1,M2.A);
}

inline rmatrix operator+(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_add<srmatrix,rmatrix,rmatrix>(M1.A,M2);
}

inline rmatrix operator+(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<rmatrix,srmatrix,rmatrix>(M1,M2.A);
}

inline rmatrix operator+(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_add<srmatrix,rmatrix_slice,rmatrix>(M1.A,M2);
}

inline rmatrix operator+(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<rmatrix_slice,srmatrix,rmatrix>(M1,M2.A);
}

inline srmatrix operator-(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline srmatrix operator-(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_sub<srmatrix,srmatrix,srmatrix,real>(M1.A,M2);
}

inline srmatrix operator-(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,srmatrix,srmatrix,real>(M1,M2.A);
}

inline rmatrix operator-(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_sub<srmatrix,rmatrix,rmatrix>(M1.A,M2);
}

inline rmatrix operator-(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<rmatrix,srmatrix,rmatrix>(M1,M2.A);
}

inline rmatrix operator-(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_sub<srmatrix,rmatrix_slice,rmatrix>(M1.A,M2);
}

inline rmatrix operator-(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<rmatrix_slice,srmatrix,rmatrix>(M1,M2.A);
}

inline rmatrix& rmatrix::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline rmatrix& rmatrix::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline rmatrix_slice& rmatrix_slice::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline rmatrix& rmatrix::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline rmatrix_slice& rmatrix_slice::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline rmatrix& rmatrix::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline rmatrix_slice& rmatrix_slice::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline bool operator==(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

inline bool operator==(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator==(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

inline bool operator==(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const srmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

inline bool operator!=(const srmatrix& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const rmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator!=(const rmatrix& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

inline bool operator!=(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

inline bool operator<(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_less<srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_less<srmatrix,srmatrix,real>(M1.A,M2);
}

inline bool operator<(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_less<srmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_less<srmatrix,rmatrix,real>(M1.A,M2);
}

inline bool operator<(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_less<rmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator<(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_less<rmatrix_slice,srmatrix,real>(M1,M2.A);
}

inline bool operator<(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_less<srmatrix,rmatrix_slice,real>(M1.A,M2);
}

inline bool operator<=(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_leq<srmatrix,srmatrix,real>(M1.A,M2);
}

inline bool operator<=(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_leq<srmatrix,rmatrix,real>(M1.A,M2);
}

inline bool operator<=(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_leq<rmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator<=(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_leq<rmatrix_slice,srmatrix,real>(M1,M2.A);
}

inline bool operator<=(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_leq<srmatrix,rmatrix_slice,real>(M1.A,M2);
}

inline bool operator>(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline bool operator>(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_greater<srmatrix,srmatrix,real>(M1.A,M2);
}

inline bool operator>(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<srmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator>(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_greater<srmatrix,rmatrix,real>(M1.A,M2);
}

inline bool operator>(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<rmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator>(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<rmatrix_slice,srmatrix,real>(M1,M2.A);
}

inline bool operator>(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_greater<srmatrix,rmatrix_slice,real>(M1.A,M2);
}

inline bool operator>=(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<srmatrix,srmatrix,real>(M1.A,M2.A);
}

inline bool operator>=(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_geq<srmatrix,srmatrix,real>(M1.A,M2);
}

inline bool operator>=(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<srmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator>=(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_geq<srmatrix,rmatrix,real>(M1.A,M2);
}

inline bool operator>=(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<rmatrix,srmatrix,real>(M1,M2.A);
}

inline bool operator>=(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<rmatrix_slice,srmatrix,real>(M1,M2.A);
}

inline bool operator>=(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_geq<srmatrix,rmatrix_slice,real>(M1.A,M2);
}

inline bool operator!(const srmatrix_slice& M) {
  return sp_m_not(M.A);
}

inline std::ostream& operator<<(std::ostream& os, const srmatrix_slice& M) {
  return sp_m_output<srmatrix,real>(os, M.A);
}

inline std::istream& operator>>(std::istream& is, srmatrix_slice& M) {
  srmatrix tmp(M.A.m,M.A.n);
  sp_m_input<srmatrix,real>(is, tmp);
  M = tmp;
  return is;
}


class srmatrix_subv {
  private:
    srmatrix_slice dat;
    bool row;
    int index;

    srmatrix_subv(srmatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

    srmatrix_subv(const srmatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r){
       if(row) index=i; else index=k;
    }

  public:
    real& operator[](const int i) {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix_subv::operator[](int)"));
#endif
        return dat.element(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix_subv::operator[](int)"));
#endif
        return dat.element(i,index);
      }
    }

    const real operator[](const int i) const {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix_subv::operator[](int)"));
#endif
        return dat(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("srmatrix_subv::operator[](int)"));
#endif
        return dat(i,index);
      }
    }

    srmatrix_subv& operator=(const srmatrix_subv& v) {
      return svsp_vv_assign(*this,srvector(v));
    }

    srmatrix_subv& operator=(const real& v) {
      return sv_vs_assign(*this,v);
    }

    srmatrix_subv& operator=(const srvector& v) {
      return svsp_vv_assign(*this,v);
    }

    srmatrix_subv& operator=(const srvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    srmatrix_subv& operator=(const rvector& v) {
      return svf_vv_assign(*this,v);
    }

    srmatrix_subv& operator=(const rvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    srmatrix_subv& operator*=(const real&);
    srmatrix_subv& operator/=(const real&);
    srmatrix_subv& operator+=(const srvector&);
    srmatrix_subv& operator+=(const srvector_slice&);
    srmatrix_subv& operator+=(const rvector&);
    srmatrix_subv& operator+=(const rvector_slice&);
    srmatrix_subv& operator-=(const srvector&);
    srmatrix_subv& operator-=(const srvector_slice&);
    srmatrix_subv& operator-=(const rvector&);
    srmatrix_subv& operator-=(const rvector_slice&);

    friend srvector operator-(const srmatrix_subv&);

    friend std::istream& operator>>(std::istream&, srmatrix_subv&);

    friend int Lb(const srmatrix_subv&);
    friend int Ub(const srmatrix_subv&);
    friend int VecLen(const srmatrix_subv&);

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

inline int Lb(const srmatrix_subv& S) {
  if(S.row)
    return Lb(S.dat, 2);
  else
    return Lb(S.dat, 1);
}

inline int Ub(const srmatrix_subv& S) {
  if(S.row)
    return Ub(S.dat, 2);
  else
    return Ub(S.dat, 1);
}

inline int VecLen(const srmatrix_subv& S) {
  return Ub(S)-Lb(S)+1;
}

inline std::ostream& operator<<(std::ostream& os, const srmatrix_subv& v) {
  os << srvector(v);
  return os;
}

inline std::istream& operator>>(std::istream& is, srmatrix_subv& v) {
  int n = 0;
  if(v.row) n=v.dat.A.n; else n=v.dat.A.m;
  srvector tmp(n);
  is >> tmp;
  v = tmp;
  return is;
}

inline srmatrix_subv srmatrix::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator[](const cxscmatrix_column&)"));
#endif
  return srmatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline srmatrix_subv srmatrix::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator[](const int)"));
#endif
  return srmatrix_subv(*this, true, i, i, lb2, ub2);
}

inline const srmatrix_subv srmatrix::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator[](const cxscmatrix_column&)"));
#endif
  return srmatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline const srmatrix_subv srmatrix::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix::operator[](const int)"));
#endif
  return srmatrix_subv(*this, true, i, i, lb2, ub2);
}

inline srmatrix_subv srmatrix_slice::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix_slice::operator[](const int)"));
#endif
  return srmatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline srmatrix_subv srmatrix_slice::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return srmatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline const srmatrix_subv srmatrix_slice::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix_slice::operator[](const int)"));
#endif
  return srmatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline const srmatrix_subv srmatrix_slice::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("srmatrix_slice::operator[](const cxscmatrix_column&)"));
#endif
  return srmatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline srvector::srvector(const srmatrix_subv& A) {
  p.reserve(A.dat.A.get_nnz());
  x.reserve(A.dat.A.get_nnz());

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

inline srvector operator-(const srmatrix_subv& v) {
 srvector s(v);
 return -s;
}

inline srvector operator*(const srmatrix_subv& v1, const real& v2) {
  return srvector(v1) * v2;
}

inline srvector operator/(const srmatrix_subv& v1, const real& v2) {
  return srvector(v1) / v2;
}

inline srvector operator*(const real& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline real operator*(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) * v2;
}

inline real operator*(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) * v2;
}

inline real operator*(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) * v2;
}

inline real operator*(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) * v2;
}

inline real operator*(const srvector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline real operator*(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline real operator*(const rvector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline real operator*(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

inline srvector operator+(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) + v2;
}

inline srvector operator+(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) + v2;
}

inline rvector operator+(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) + v2;
}

inline rvector operator+(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) + v2;
}

inline srvector operator+(const srvector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline srvector operator+(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline rvector operator+(const rvector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline rvector operator+(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

inline srvector operator-(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) - v2;
}

inline srvector operator-(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) - v2;
}

inline rvector operator-(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) - v2;
}

inline rvector operator-(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) - v2;
}

inline srvector operator-(const srvector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline srvector operator-(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline rvector operator-(const rvector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline rvector operator-(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

inline srmatrix_subv& srmatrix_subv::operator*=(const real& v) {
  *this = *this * v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator/=(const real& v) {
  *this = *this * v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator+=(const srvector& v) {
  *this = *this + v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator+=(const srvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator+=(const rvector& v) {
  *this = *this + v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator+=(const rvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator-=(const srvector& v) {
  *this = *this - v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator-=(const srvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator-=(const rvector& v) {
  *this = *this - v;
  return *this;
}

inline srmatrix_subv& srmatrix_subv::operator-=(const rvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator+=(const srmatrix_subv& v) {
  *this += rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator+=(const srvector& v) {
  *this += rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator+=(const srvector_slice& v) {
  *this += rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator-=(const srmatrix_subv& v) {
  *this -= rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline rmatrix_subv& rmatrix_subv::operator=(const srmatrix_subv& v) {
  *this = rvector(v);
  return *this;
}

inline bool operator==(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) == v2;
}

inline bool operator==(const srvector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator==(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

inline bool operator!=(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) != v2;
}

inline bool operator!=(const srvector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator!=(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

inline bool operator<(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) < v2;
}

inline bool operator<(const srvector& v1, const srmatrix_subv& v2) {
  return v1 < srvector(v2);
}

inline bool operator<(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 < srvector(v2);
}

inline bool operator<(const rvector& v1, const srmatrix_subv& v2) {
  return v1 < srvector(v2);
}

inline bool operator<(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 < srvector(v2);
}

inline bool operator<=(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) <= v2;
}

inline bool operator<=(const srvector& v1, const srmatrix_subv& v2) {
  return v1 <= srvector(v2);
}

inline bool operator<=(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 <= srvector(v2);
}

inline bool operator<=(const rvector& v1, const srmatrix_subv& v2) {
  return v1 <= srvector(v2);
}

inline bool operator<=(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 <= srvector(v2);
}

inline bool operator>(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) > v2;
}

inline bool operator>(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) > v2;
}

inline bool operator>(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) > v2;
}

inline bool operator>(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) > v2;
}

inline bool operator>(const srvector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const rvector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

inline bool operator>=(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) >= v2;
}

inline bool operator>=(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) >= v2;
}

inline bool operator>=(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) >= v2;
}

inline bool operator>=(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) >= v2;
}

inline bool operator>=(const srvector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const rvector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator>=(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

inline bool operator!(const srmatrix_subv& v) {
  return sv_v_not(v);
}

inline void accumulate(dotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<dotprecision,srvector,srvector,sparse_dot>(dot, srvector(v1), srvector(v2));
}

inline void accumulate(dotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  spsp_vv_accu<dotprecision,srvector,srvector,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate(dotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accu<dotprecision,srvector,srvector_slice,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate(dotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  spf_vv_accu<dotprecision,srvector,rvector,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate(dotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accu<dotprecision,srvector,rvector_slice,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate(dotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<dotprecision,srvector,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate(dotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accu<dotprecision,srvector_slice,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate(dotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<dotprecision,rvector,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate(dotprecision& dot, const rvector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<dotprecision,rvector_slice,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(dotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accuapprox<dotprecision,srvector,srvector,sparse_dot>(dot, srvector(v1), srvector(v2));
}

inline void accumulate_approx(dotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  spsp_vv_accuapprox<dotprecision,srvector,srvector,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(dotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accuapprox<dotprecision,srvector,srvector_slice,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(dotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  spf_vv_accuapprox<dotprecision,srvector,rvector,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(dotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accuapprox<dotprecision,srvector,rvector_slice,sparse_dot>(dot, srvector(v1), v2);
}

inline void accumulate_approx(dotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  spsp_vv_accuapprox<dotprecision,srvector,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(dotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accuapprox<dotprecision,srvector_slice,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(dotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  fsp_vv_accuapprox<dotprecision,rvector,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate_approx(dotprecision& dot, const rvector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accuapprox<dotprecision,rvector_slice,srvector,sparse_dot>(dot, v1, srvector(v2));
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),srvector(v2));
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),srvector(v2));
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,srvector(v1),srvector(v2));
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const srvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const srvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const rvector& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const rvector_slice& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const srmatrix_subv& v2) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  dot += tmp;
}

}  //namespace cxsc;

#include "sparsematrix.inl"

#endif 
