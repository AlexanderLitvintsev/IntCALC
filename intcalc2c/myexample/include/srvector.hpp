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

/* CVS $Id: srvector.hpp,v 1.7 2010/12/21 09:54:02 cxsc Exp $ */

#ifndef _CXSC_SRVECTOR_HPP_INCLUDED
#define _CXSC_SRVECTOR_HPP_INCLUDED

#include <real.hpp>
#include <intvector.hpp>
#include <rvector.hpp>
#include <intmatrix.hpp>
#include <vector>
#include <map>
#include <iostream>
#include <except.hpp>
#include <cidot.hpp>
#include <sparsedot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;

class scvector;
class sivector;
class sivector_slice;
class scivector;

class srvector {
  private:
    std::vector<int> p;
    std::vector<real> x;
    int lb;
    int ub;
    int n; 

  public:
    srvector() : lb(0), ub(-1) , n(0) {
    }

    explicit srvector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    srvector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    srvector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    srvector(const int n, const int nnz, const intvector& index, const rvector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[Lb(values)+i] != 0.0) {
          p.push_back(index[Lb(index)+i]);
          x.push_back(values[Lb(values)+i]);
        }
      }
      
    }

    srvector(const srvector_slice&);
    srvector(const srmatrix_subv& A);

    std::vector<int>& row_indices() {
      return p;
    }

    std::vector<real>& values() {
      return x;
    }

    const std::vector<int>& row_indices() const {
      return p;
    }

    const std::vector<real>& values() const {
      return x;
    }

    int get_nnz() const {
      return x.size();
    }

    real density() const {
      return (double)x.size()/n;
    }

    void dropzeros() {
      for(int i=0 ; i<get_nnz() ; i++) {
        if(x[i] == 0.0) {
           x.erase(x.begin()+i);
           p.erase(p.begin()+i);
        }
      }
    }

    srvector& operator=(const real& v) {
      return sp_vs_assign<srvector,real,real>(*this,v);
    }

    srvector& operator=(const rvector& v) {
      return spf_vv_assign<srvector,rvector,real>(*this,v);
    }

    srvector& operator=(const rvector_slice& v) {
      return spf_vv_assign<srvector,rvector_slice,real>(*this,v);
    }

    srvector& operator=(const srvector_slice&);

    real& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, 0.0);

      return x[k];
    }

    real operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    real operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator()(const int)"));
#endif
      real r = 0.0;

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    srvector_slice operator()(const int, const int);

    srvector operator()(const intvector& per) {
      srvector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,real> work;
      for(int i=0 ; i<get_nnz() ; i++)
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
 
      for(std::map<int,real>::iterator it=work.begin() ; it!=work.end() ; it++) {
         v.p.push_back(it->first);
         v.x.push_back(it->second);
      }

      return v;
    }

    srvector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    srvector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    srvector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    srvector& operator+=(const rvector& v) 
    {
      return spf_vv_addassign(*this,v);
    }

    srvector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    srvector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    srvector& operator+=(const srvector_slice&);

    srvector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    srvector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    srvector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    srvector& operator-=(const srvector_slice&);

    friend int Lb(const srvector&);
    friend int Ub(const srvector&);
    friend void SetLb(srvector&, const int);
    friend void SetUb(srvector&, const int);

    friend int VecLen(const srvector&);
    friend srvector Re(const scvector&);
    friend srvector Im(const scvector&);
    friend srvector Inf(const sivector&);
    friend srvector Sup(const sivector&);
    friend srvector InfRe(const scivector&);
    friend srvector SupRe(const scivector&);
    friend srvector InfIm(const scivector&);
    friend srvector SupIm(const scivector&);
    friend srvector mid(const sivector&);
    friend srvector diam(const sivector&);
    friend srvector absmin(const sivector&);
    friend srvector absmax(const sivector&);
    friend srvector mid(const sivector_slice&);
    friend srvector diam(const sivector_slice&);


    friend class srvector_slice;
    friend class scvector_slice;
    friend class scvector;
    friend class sivector_slice;
    friend class sivector;
    friend class scivector_slice;
    friend class scivector;
    friend class srmatrix_subv;
    friend class rvector;
    friend class rvector_slice;
    friend class ivector;
    friend class ivector_slice;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;


#include "vector_friend_declarations.inl"
};

inline rvector::rvector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new real[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline rvector& rvector::operator=(const srvector& v) {
  return fsp_vv_assign<rvector,srvector,real>(*this,v);
}

inline rvector& rvector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<rvector,srvector_slice,real>(*this,v);
}


inline void SetLb(srvector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

inline void SetUb(srvector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

inline int Lb(const srvector& v) {
  return v.lb;
}

inline int Ub(const srvector& v) {
  return v.ub;
}

inline int VecLen(const srvector& v) {
  return v.n;
}

inline void Resize(srvector& v) {
  sp_v_resize(v);
}

inline void Resize(srvector& v, const int n) {
  sp_v_resize(v,n);
}

inline void Resize(srvector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

inline srvector operator-(const srvector& v) {
  return sp_v_negative(v);
}

inline real operator*(const srvector& v1, const rvector& v2) {
  return spf_vv_mult<srvector,rvector,real,sparse_dot>(v1,v2);
}

inline real operator*(const rvector& v1, const srvector& v2) {
  return fsp_vv_mult<rvector,srvector,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_mult<srvector,rvector_slice,real,sparse_dot>(v1,v2);
}

inline real operator*(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<rvector_slice,srvector,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector& v1, const srvector& v2) {
  return spsp_vv_mult<srvector,srvector,real,sparse_dot>(v1,v2);
}

inline srvector operator*(const srvector& v, const real& s) {
  return sp_vs_mult<srvector,real,srvector>(v,s);
}

inline srvector operator/(const srvector& v, const real& s) {
  return sp_vs_div<srvector,real,srvector>(v,s);
}

inline srvector operator*(const real& s, const srvector& v) {
  return sp_sv_mult<real,srvector,srvector>(s,v);
}

inline rvector operator+(const rvector& v1, const srvector& v2) {
  return fsp_vv_add<rvector,srvector,rvector>(v1,v2);
}

inline rvector operator+(const srvector& v1, const rvector& v2) {
  return spf_vv_add<srvector,rvector,rvector>(v1,v2);
}

inline rvector operator+(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_add<rvector_slice,srvector,rvector>(v1,v2);
}

inline rvector operator+(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_add<srvector,rvector_slice,rvector>(v1,v2);
}

inline srvector operator+(const srvector& v1, const srvector& v2) {
  return spsp_vv_add<srvector,srvector,srvector,real>(v1,v2);
}

inline rvector operator-(const rvector& v1, const srvector& v2) {
  return fsp_vv_sub<rvector,srvector,rvector>(v1,v2);
}

inline rvector operator-(const srvector& v1, const rvector& v2) {
  return spf_vv_sub<srvector,rvector,rvector>(v1,v2);
}

inline rvector operator-(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<rvector_slice,srvector,rvector>(v1,v2);
}

inline rvector operator-(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_sub<srvector,rvector_slice,rvector>(v1,v2);
}

inline srvector operator-(const srvector& v1, const srvector& v2) {
  return spsp_vv_sub<srvector,srvector,srvector,real>(v1,v2);
}

inline rvector& rvector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline rvector& rvector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline bool operator==(const srvector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator<(const srvector& v1, const srvector& v2) {
  return spsp_vv_less<srvector,srvector,real>(v1,v2);
}

inline bool operator<(const srvector& v1, const rvector& v2) {
  return spf_vv_less<srvector,rvector,real>(v1,v2);
}

inline bool operator<(const rvector& v1, const srvector& v2) {
  return fsp_vv_less<rvector,srvector,real>(v1,v2);
}

inline bool operator<(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_less<srvector,rvector_slice,real>(v1,v2);
}

inline bool operator<(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_less<rvector_slice,srvector,real>(v1,v2);
}

inline bool operator<=(const srvector& v1, const srvector& v2) {
  return spsp_vv_leq<srvector,srvector,real>(v1,v2);
}

inline bool operator<=(const srvector& v1, const rvector& v2) {
  return spf_vv_leq<srvector,rvector,real>(v1,v2);
}

inline bool operator<=(const rvector& v1, const srvector& v2) {
  return fsp_vv_leq<rvector,srvector,real>(v1,v2);
}

inline bool operator<=(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_leq<srvector,rvector_slice,real>(v1,v2);
}

inline bool operator<=(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_leq<rvector_slice,srvector,real>(v1,v2);
}

inline bool operator>(const srvector& v1, const srvector& v2) {
  return spsp_vv_greater<srvector,srvector,real>(v1,v2);
}

inline bool operator>(const srvector& v1, const rvector& v2) {
  return spf_vv_greater<srvector,rvector,real>(v1,v2);
}

inline bool operator>(const rvector& v1, const srvector& v2) {
  return fsp_vv_greater<rvector,srvector,real>(v1,v2);
}

inline bool operator>(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_greater<srvector,rvector_slice,real>(v1,v2);
}

inline bool operator>(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_greater<rvector_slice,srvector,real>(v1,v2);
}

inline bool operator>=(const srvector& v1, const srvector& v2) {
  return spsp_vv_geq<srvector,srvector,real>(v1,v2);
}

inline bool operator>=(const srvector& v1, const rvector& v2) {
  return spf_vv_geq<srvector,rvector,real>(v1,v2);
}

inline bool operator>=(const rvector& v1, const srvector& v2) {
  return fsp_vv_geq<rvector,srvector,real>(v1,v2);
}

inline bool operator>=(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_geq<srvector,rvector_slice,real>(v1,v2);
}

inline bool operator>=(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_geq<rvector_slice,srvector,real>(v1,v2);
}

inline bool operator!(const srvector& x) {
  return sp_v_not(x);
}

inline std::ostream& operator<<(std::ostream& os, const srvector& v) {
  return sp_v_output<srvector,real>(os,v);
}

inline std::istream& operator>>(std::istream& is, srvector& v) {
  return sp_v_input<srvector,real>(is,v);
}

class srvector_slice {
  private:
    std::vector<int>& p;
    std::vector<real>& x;
    srvector& orig;
    int start,end;
    int lb;
    int ub;
    int n;
    int nnz;
    int offset;

    srvector_slice(srvector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1)  {
      int i;

      for(i=0 ; i<v.get_nnz() && p[i]<lb-v.lb ; i++);

      start = i;

      for(i=start ; i<v.get_nnz() && p[i]<=ub-v.lb ; i++);

      end = i-1;

      nnz = end-start+1;
      offset = lb-v.lb;
    }

  public:

    int get_nnz() const {
      return nnz;
    }

    real density() const {
      return (double)nnz/n;
    }

    real& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, 0.0);
      end++;

      return x[k];
    }

    real operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    real operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator()(const int)"));
#endif
      real r = 0.0;

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    srvector_slice& operator=(const real& v) {
      return sl_vs_assign<srvector_slice,real,real,std::vector<real>::iterator>(*this,v);
    }

    srvector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<srvector_slice,srvector_slice,real,std::vector<real>::iterator>(*this,v);
    }

    srvector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<srvector_slice,srvector,real,std::vector<real>::iterator>(*this,v);
    }

    srvector_slice& operator=(const rvector& v) {
      return slf_vv_assign<srvector_slice,rvector,real,std::vector<real>::iterator>(*this,v);
    }

    srvector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<srvector_slice,rvector,real,std::vector<real>::iterator>(*this,v);
    }

    srvector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    srvector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    srvector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<srvector_slice,rvector,real>(*this,v);
    }

    srvector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<srvector_slice,rvector_slice,real>(*this,v);
    }

    srvector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    srvector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    srvector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<srvector_slice,rvector,real>(*this,v);
    }

    srvector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<srvector_slice,rvector_slice,real>(*this,v);
    }

    srvector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    srvector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }


    friend int Lb(const srvector_slice&);
    friend int Ub(const srvector_slice&);
    friend int VecLen(const srvector_slice&);

    friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
    friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class scvector;
    friend class sivector;
    friend class scivector;
    friend class srmatrix_subv;
    friend class rvector;
    friend class rvector_slice;
    friend class ivector;
    friend class ivector_slice;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline rvector::rvector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new real[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline rvector_slice& rvector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline rvector_slice& rvector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline srvector::srvector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline srvector& srvector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<srvector,srvector_slice,real>(*this,v);
}

inline srvector_slice srvector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator()(const int,const int)"));
#endif
  return srvector_slice(*this,i,j);
}

inline srvector operator-(const srvector_slice& v) {
  return sl_v_negative<srvector_slice,srvector>(v);
}

inline int Lb(const srvector_slice& v) {
  return v.lb;
}

inline int Ub(const srvector_slice& v) {
  return v.ub;
}

inline int VecLen(const srvector_slice& v) {
  return v.n;
}

inline real operator*(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_mult<srvector_slice,rvector,real,sparse_dot>(v1,v2);
}

inline real operator*(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<rvector,srvector_slice,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<srvector_slice,rvector_slice,real,sparse_dot>(v1,v2);
}

inline real operator*(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<rvector_slice,srvector_slice,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<srvector,srvector_slice,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<srvector_slice,srvector,real,sparse_dot>(v1,v2);
}

inline real operator*(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<srvector_slice,srvector_slice,real,sparse_dot>(v1,v2);
}

inline srvector operator*(const srvector_slice& v, const real& s) {
  return sp_vs_mult<srvector_slice,real,srvector>(v,s);
}

inline srvector operator/(const srvector_slice& v, const real& s) {
  return sp_vs_div<srvector_slice,real,srvector>(v,s);
}

inline srvector operator*(const real& s, const srvector_slice& v) {
  return sp_sv_mult<real,srvector_slice,srvector>(s,v);
}

inline rvector operator+(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_add<rvector,srvector_slice,rvector>(v1,v2);
}

inline rvector operator+(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_add<srvector_slice,rvector,rvector>(v1,v2);
}

inline rvector operator+(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<rvector_slice,srvector_slice,rvector>(v1,v2);
}

inline rvector operator+(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<srvector_slice,rvector_slice,rvector>(v1,v2);
}

inline srvector operator+(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<srvector_slice,srvector_slice,srvector,real>(v1,v2);
}

inline srvector operator+(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_add<srvector,srvector_slice,srvector,real>(v1,v2);
}

inline srvector operator+(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_add<srvector_slice,srvector,srvector,real>(v1,v2);
}

inline rvector operator-(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<rvector,srvector_slice,rvector>(v1,v2);
}

inline rvector operator-(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_sub<srvector_slice,rvector,rvector>(v1,v2);
}

inline rvector operator-(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<rvector_slice,srvector_slice,rvector>(v1,v2);
}

inline rvector operator-(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<srvector_slice,rvector_slice,rvector>(v1,v2);
}

inline srvector operator-(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<srvector_slice,srvector_slice,srvector,real>(v1,v2);
}

inline srvector operator-(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<srvector,srvector_slice,srvector,real>(v1,v2);
}

inline srvector operator-(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<srvector_slice,srvector,srvector,real>(v1,v2);
}

inline rvector& rvector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline srvector& srvector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline rvector& rvector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline srvector& srvector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline bool operator==(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_less<srvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_less<srvector_slice,srvector,real>(v1,v2);
}

inline bool operator<(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_less<srvector,srvector_slice,real>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_less<srvector_slice,rvector,real>(v1,v2);
}

inline bool operator<(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_less<rvector,srvector_slice,real>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_less<srvector_slice,rvector_slice,real>(v1,v2);
}

inline bool operator<(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_less<rvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_leq<srvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_leq<srvector_slice,srvector,real>(v1,v2);
}

inline bool operator<=(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_leq<srvector,srvector_slice,real>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_leq<srvector_slice,rvector,real>(v1,v2);
}

inline bool operator<=(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_leq<rvector,srvector_slice,real>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_leq<srvector_slice,rvector_slice,real>(v1,v2);
}

inline bool operator<=(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_leq<rvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator>(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_greater<srvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator>(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_greater<srvector_slice,srvector,real>(v1,v2);
}

inline bool operator>(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_greater<srvector,srvector_slice,real>(v1,v2);
}

inline bool operator>(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_greater<srvector_slice,rvector,real>(v1,v2);
}

inline bool operator>(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_greater<rvector,srvector_slice,real>(v1,v2);
}

inline bool operator>(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_greater<srvector_slice,rvector_slice,real>(v1,v2);
}

inline bool operator>(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_greater<rvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator>=(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_geq<srvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator>=(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_geq<srvector_slice,srvector,real>(v1,v2);
}

inline bool operator>=(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_geq<srvector,srvector_slice,real>(v1,v2);
}

inline bool operator>=(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_geq<srvector_slice,rvector,real>(v1,v2);
}

inline bool operator>=(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_geq<rvector,srvector_slice,real>(v1,v2);
}

inline bool operator>=(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_geq<srvector_slice,rvector_slice,real>(v1,v2);
}

inline bool operator>=(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_geq<rvector_slice,srvector_slice,real>(v1,v2);
}

inline bool operator!(const srvector_slice& x) {
  return sl_v_not(x);
}

inline std::ostream& operator<<(std::ostream& os, const srvector_slice& v) {
  return sl_v_output<srvector_slice, real>(os,v);
}

inline std::istream& operator>>(std::istream& is, srvector_slice& v) {
  return sl_v_input<srvector_slice, real>(is,v);
}

inline void accumulate(dotprecision& dot, const srvector& x, const srvector& y) {
  spsp_vv_accu<dotprecision,srvector,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector& x, const rvector& y) {
  spf_vv_accu<dotprecision,srvector,rvector,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector& x, const rvector_slice& y) {
  spf_vv_accu<dotprecision,srvector,rvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const rvector& x, const srvector& y) {
  fsp_vv_accu<dotprecision,rvector,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const rvector_slice& x, const srvector& y) {
  fsp_vv_accu<dotprecision,rvector_slice,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector_slice& x, const rvector& y) {
  slf_vv_accu<dotprecision,srvector_slice,rvector,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  slf_vv_accu<dotprecision,srvector_slice,rvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const rvector& x, const srvector_slice& y) {
  fsl_vv_accu<dotprecision,rvector,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<dotprecision,rvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<dotprecision,srvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector& x, const srvector_slice& y) {
  spsl_vv_accu<dotprecision,srvector,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate(dotprecision& dot, const srvector_slice& x, const srvector& y) {
  slsp_vv_accu<dotprecision,srvector_slice,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector& x, const srvector& y) {
  spsp_vv_accuapprox<dotprecision,srvector,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector& x, const rvector& y) {
  spf_vv_accuapprox<dotprecision,srvector,rvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector& x, const rvector_slice& y) {
  spf_vv_accuapprox<dotprecision,srvector,rvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const rvector& x, const srvector& y) {
  fsp_vv_accuapprox<dotprecision,rvector,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const rvector_slice& x, const srvector& y) {
  fsp_vv_accuapprox<dotprecision,rvector_slice,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const rvector& y) {
  slf_vv_accuapprox<dotprecision,srvector_slice,rvector,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  slf_vv_accuapprox<dotprecision,srvector_slice,rvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const rvector& x, const srvector_slice& y) {
  fsl_vv_accuapprox<dotprecision,rvector,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  fsl_vv_accuapprox<dotprecision,rvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  slsl_vv_accuapprox<dotprecision,srvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector& x, const srvector_slice& y) {
  spsl_vv_accuapprox<dotprecision,srvector,srvector_slice,sparse_dot>(dot,x,y);
}

inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const srvector& y) {
  slsp_vv_accuapprox<dotprecision,srvector_slice,srvector,sparse_dot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}


} //namespace cxsc

#include "sparsevector.inl"

#endif
