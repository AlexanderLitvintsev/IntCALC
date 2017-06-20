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

/* CVS $Id: scvector.hpp,v 1.8 2010/12/21 09:54:02 cxsc Exp $ */

#ifndef _CXSC_SCVECTOR_HPP_INCLUDED
#define _CXSC_SCVECTOR_HPP_INCLUDED

#include <complex.hpp>
#include <cvector.hpp>
#include <vector>
#include <iostream>
#include <srvector.hpp>
#include <sparsecdot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;
class scvector_slice;
class scmatrix;
class scmatrix_slice;
class scmatrix_subv;
class scivector;
class scivector_slice;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;

class scvector {
  private:
    std::vector<int> p;
    std::vector<complex> x;
    int lb;
    int ub;
    int n; 

  public:
    scvector() : lb(0), ub(-1) , n(0) {
    }

    explicit scvector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    scvector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    scvector(const cvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    scvector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(complex(v[i]));
          }
        }
    }

    scvector(const int n, const int nnz, const intvector& index, const cvector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i+Lb(values)] != 0.0) {
          p.push_back(index[i+Lb(index)]);
          x.push_back(values[i+Lb(values)]);
        }
      }
      
    }

    scvector(const srvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(complex(v.x[i]));
    }

    scvector(const srvector_slice&);
    scvector(const scvector_slice&);
    scvector(const srmatrix_subv& A);
    scvector(const scmatrix_subv& A);

    std::vector<int>& row_indices() {
      return p;
    }

    std::vector<complex>& values() {
      return x;
    }

    const std::vector<int>& row_indices() const {
      return p;
    }

    const std::vector<complex>& values() const {
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

    scvector& operator=(const srvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = complex(v.x[i]);
      return *this;
    } 

    scvector& operator=(const real& v) {
      return sp_vs_assign<scvector,real,complex>(*this,v);
    }

    scvector& operator=(const complex& v) {
      return sp_vs_assign<scvector,complex,complex>(*this,v);
    }

    scvector& operator=(const rvector& v) {
      return spf_vv_assign<scvector,rvector,complex>(*this,v);
    }

    scvector& operator=(const cvector& v) {
      return spf_vv_assign<scvector,cvector,complex>(*this,v);
    }

    scvector& operator=(const rvector_slice& v) {
      return spf_vv_assign<scvector,rvector_slice,complex>(*this,v);
    }

    scvector& operator=(const scvector_slice&);
    scvector& operator=(const srvector_slice&);

    complex& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, complex(0.0));

      return x[k];
    }

    complex operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    complex operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator()(const int)"));
#endif
      complex r(0.0);

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    scvector operator()(const intvector& per) {
      scvector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,complex> work;
      for(int i=0 ; i<get_nnz() ; i++)
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
 
      for(std::map<int,complex>::iterator it=work.begin() ; it!=work.end() ; it++) {
         v.p.push_back(it->first);
         v.x.push_back(it->second);
      }

      return v;
    }

    scvector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    scvector_slice operator()(const int, const int);

    scvector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    scvector& operator*=(const complex& s) {
      return sp_vs_multassign(*this,s);
    }

    scvector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    scvector& operator/=(const complex& s) {
      return sp_vs_divassign(*this,s);
    }

    scvector& operator+=(const rvector& v) {
      return spf_vv_addassign(*this,v);
    }

    scvector& operator+=(const cvector& v) {
      return spf_vv_addassign(*this,v);
    }

    scvector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scvector& operator+=(const cvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scvector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scvector& operator+=(const scvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scvector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    scvector& operator-=(const cvector& v) {
      return spf_vv_subassign(*this,v);
    }

    scvector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scvector& operator-=(const cvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scvector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scvector& operator-=(const scvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scvector& operator+=(const srvector_slice&);
    scvector& operator+=(const scvector_slice&);
    scvector& operator-=(const srvector_slice&);
    scvector& operator-=(const scvector_slice&);

    friend void SetLb(scvector&, const int);
    friend void SetUb(scvector&, const int);
    friend int Lb(const scvector&);
    friend int Ub(const scvector&);
    friend srvector Re(const scvector&);
    friend srvector Im (const scvector&);
    friend scvector Inf(const scivector&);
    friend scvector Sup (const scivector&);
    friend scvector mid(const scivector&);
    friend scvector diam(const scivector&);
    friend scvector mid(const scivector_slice&);
    friend scvector diam(const scivector_slice&);
    friend int VecLen(const scvector&);

    friend class srvector_slice;
    friend class scvector_slice;
    friend class scivector_slice;
    friend class scivector;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline cvector::cvector(const scvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector::cvector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector& cvector::operator=(const scvector& v) {
  return fsp_vv_assign<cvector,scvector,complex>(*this,v);
}

inline cvector& cvector::operator=(const scvector_slice& v) {
  return fsl_vv_assign<cvector,scvector_slice,complex>(*this,v);
}

inline cvector& cvector::operator=(const srvector& v) {
  return fsp_vv_assign<cvector,srvector,complex>(*this,v);
}

inline cvector& cvector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<cvector,srvector_slice,complex>(*this,v);
}

inline void SetLb(scvector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

inline void SetUb(scvector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

inline int Lb(const scvector& v) {
  return v.lb;
}

inline int Ub(const scvector& v) {
  return v.ub;
}

inline srvector Re(const scvector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Re(v.x[i]);
  return res;
}

inline srvector Im(const scvector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Im(v.x[i]);
  return res;
}

inline int VecLen(const scvector& v) {
  return v.n;
}

inline void Resize(scvector& v) {
  sp_v_resize(v);
}

inline void Resize(scvector& v, const int n) {
  sp_v_resize(v,n);
}

inline void Resize(scvector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

inline scvector operator-(const scvector& v) {
  return sp_v_negative(v);
}

inline complex operator*(const scvector& v1, const cvector& v2) {
  return spf_vv_mult<scvector,cvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const rvector& v2) {
  return spf_vv_mult<scvector,rvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector& v1, const cvector& v2) {
  return spf_vv_mult<srvector,cvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const rvector& v1, const scvector& v2) {
  return fsp_vv_mult<rvector,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector& v1, const srvector& v2) {
  return fsp_vv_mult<cvector,srvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector& v1, const scvector& v2) {
  return fsp_vv_mult<cvector,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_mult<scvector,rvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_mult<scvector,cvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_mult<srvector,cvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<cvector_slice,srvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<cvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<rvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const srvector& v2) {
  return spsp_vv_mult<scvector,srvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector& v1, const scvector& v2) {
  return spsp_vv_mult<srvector,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const scvector& v2) {
  return spsp_vv_mult<scvector,scvector,complex,sparse_cdot>(v1,v2);
}

inline scvector operator*(const scvector& v, const real& s) {
  return sp_vs_mult<scvector,real,scvector>(v,s);
}

inline scvector operator*(const scvector& v, const complex& s) {
  return sp_vs_mult<scvector,complex,scvector>(v,s);
}

inline scvector operator*(const srvector& v, const complex& s) {
  return sp_vs_mult<srvector,complex,scvector>(v,s);
}

inline scvector operator/(const scvector& v, const real& s) {
  return sp_vs_div<scvector,real,scvector>(v,s);
}

inline scvector operator/(const scvector& v, const complex& s) {
  return sp_vs_div<scvector,complex,scvector>(v,s);
}

inline scvector operator/(const srvector& v, const complex& s) {
  return sp_vs_div<srvector,complex,scvector>(v,s);
}

inline scvector operator*(const real& s, const scvector& v) {
  return sp_sv_mult<real,scvector,scvector>(s,v);
}

inline scvector operator*(const complex& s, const scvector& v) {
  return sp_sv_mult<complex,scvector,scvector>(s,v);
}

inline scvector operator*(const complex& s, const srvector& v) {
  return sp_sv_mult<complex,srvector,scvector>(s,v);
}

inline cvector operator+(const cvector& v1, const srvector& v2) {
  return fsp_vv_add<cvector,srvector,cvector>(v1,v2);
}

inline cvector operator+(const rvector& v1, const scvector& v2) {
  return fsp_vv_add<rvector,scvector,cvector>(v1,v2);
}

inline cvector operator+(const cvector& v1, const scvector& v2) {
  return fsp_vv_add<cvector,scvector,cvector>(v1,v2);
}

inline cvector operator+(const scvector& v1, const rvector& v2) {
  return spf_vv_add<scvector,rvector,cvector>(v1,v2);
}

inline cvector operator+(const srvector& v1, const cvector& v2) {
  return spf_vv_add<srvector,cvector,cvector>(v1,v2);
}

inline cvector operator+(const scvector& v1, const cvector& v2) {
  return spf_vv_add<scvector,cvector,cvector>(v1,v2);
}

inline cvector operator+(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_add<cvector_slice,srvector,cvector>(v1,v2);
}

inline cvector operator+(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_add<rvector_slice,scvector,cvector>(v1,v2);
}

inline cvector operator+(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_add<cvector_slice,scvector,cvector>(v1,v2);
}

inline cvector operator+(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_add<scvector,rvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_add<srvector,cvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_add<scvector,cvector_slice,cvector>(v1,v2);
}

inline scvector operator+(const scvector& v1, const srvector& v2) {
  return spsp_vv_add<scvector,srvector,scvector,complex>(v1,v2);
}

inline scvector operator+(const srvector& v1, const scvector& v2) {
  return spsp_vv_add<srvector,scvector,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector& v1, const scvector& v2) {
  return spsp_vv_add<scvector,scvector,scvector,complex>(v1,v2);
}

inline cvector operator-(const cvector& v1, const srvector& v2) {
  return fsp_vv_sub<cvector,srvector,cvector>(v1,v2);
}

inline cvector operator-(const rvector& v1, const scvector& v2) {
  return fsp_vv_sub<rvector,scvector,cvector>(v1,v2);
}

inline cvector operator-(const cvector& v1, const scvector& v2) {
  return fsp_vv_sub<cvector,scvector,cvector>(v1,v2);
}

inline cvector operator-(const scvector& v1, const rvector& v2) {
  return spf_vv_sub<scvector,rvector,cvector>(v1,v2);
}

inline cvector operator-(const srvector& v1, const cvector& v2) {
  return spf_vv_sub<srvector,cvector,cvector>(v1,v2);
}

inline cvector operator-(const scvector& v1, const cvector& v2) {
  return spf_vv_sub<scvector,cvector,cvector>(v1,v2);
}

inline cvector operator-(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<cvector_slice,srvector,cvector>(v1,v2);
}

inline cvector operator-(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<rvector_slice,scvector,cvector>(v1,v2);
}

inline cvector operator-(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<cvector_slice,scvector,cvector>(v1,v2);
}

inline cvector operator-(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_sub<scvector,rvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_sub<srvector,cvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_sub<scvector,cvector_slice,cvector>(v1,v2);
}

inline scvector operator-(const scvector& v1, const srvector& v2) {
  return spsp_vv_sub<scvector,srvector,scvector,complex>(v1,v2);
}

inline scvector operator-(const srvector& v1, const scvector& v2) {
  return spsp_vv_sub<srvector,scvector,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector& v1, const scvector& v2) {
  return spsp_vv_sub<scvector,scvector,scvector,complex>(v1,v2);
}

inline cvector& cvector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector& cvector::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline cvector& cvector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector& cvector::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline bool operator==(const scvector& v1, const scvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const scvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const cvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const cvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const scvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const scvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const cvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const cvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const cvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const cvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const scvector& v) {
  return sp_v_output<scvector,complex>(os,v);
}

inline std::istream& operator>>(std::istream& is, scvector& v) {
  return sp_v_input<scvector,complex>(is,v);
}


class scvector_slice {
  private:
    std::vector<int>& p;
    std::vector<complex>& x;
    scvector& orig;
    int start,end;
    int lb;
    int ub;
    int n;
    int nnz;
    int offset;

    scvector_slice(scvector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1)  {
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

    complex& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, complex(0.0));
      end++;

      return x[k];
    }

    complex operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    complex operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator()(const int)"));
#endif
      complex r(0.0);

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    scvector_slice& operator=(const real& v) {
      return sl_vs_assign<scvector_slice,real,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const complex& v) {
      return sl_vs_assign<scvector_slice,complex,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<scvector_slice,srvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const scvector_slice& v) {
      return slsl_vv_assign<scvector_slice,scvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<scvector_slice,srvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const scvector& v) {
      return slsp_vv_assign<scvector_slice,scvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const rvector& v) {
      return slf_vv_assign<scvector_slice,rvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const cvector& v) {
      return slf_vv_assign<scvector_slice,cvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<scvector_slice,rvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator=(const cvector_slice& v) {
      return slf_vv_assign<scvector_slice,cvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    scvector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    scvector_slice& operator*=(const complex& s) {
      return sl_vs_multassign(*this,s);
    }

    scvector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    scvector_slice& operator/=(const complex& s) {
      return sl_vs_divassign(*this,s);
    }

    scvector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<scvector_slice,rvector,complex>(*this,v);
    }

    scvector_slice& operator+=(const cvector& v) {
      return slf_vv_addassign<scvector_slice,cvector,complex>(*this,v);
    }

    scvector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<scvector_slice,rvector_slice,complex>(*this,v);
    }

    scvector_slice& operator+=(const cvector_slice& v) {
      return slf_vv_addassign<scvector_slice,cvector_slice,complex>(*this,v);
    }

    scvector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scvector_slice& operator+=(const scvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scvector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scvector_slice& operator+=(const scvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scvector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<scvector_slice,rvector,complex>(*this,v);
    }

    scvector_slice& operator-=(const cvector& v) {
      return slf_vv_subassign<scvector_slice,cvector,complex>(*this,v);
    }

    scvector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<scvector_slice,rvector_slice,complex>(*this,v);
    }

    scvector_slice& operator-=(const cvector_slice& v) {
      return slf_vv_subassign<scvector_slice,cvector_slice,complex>(*this,v);
    }

    scvector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scvector_slice& operator-=(const scvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scvector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    scvector_slice& operator-=(const scvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    friend int Lb(const scvector_slice&);
    friend int Ub(const scvector_slice&);
    friend srvector Re(const scvector_slice&);
    friend srvector Im(const scvector_slice&);
    friend int VecLen(const scvector_slice&);

//     friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
//     friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class scvector;
    friend class scivector;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline cvector::cvector(const scvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector::cvector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector_slice& cvector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const scvector& v) {
  *this = cvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const scvector_slice& v) {
  *this = cvector(v);
  return *this;
}

inline scvector::scvector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(complex(s.x[i]));
  }

}

inline scvector::scvector(const scvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n) {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline scvector& scvector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<scvector,srvector_slice,complex>(*this,v);
}

inline scvector& scvector::operator=(const scvector_slice& v) {
  return spsl_vv_assign<scvector,scvector_slice,complex>(*this,v);
}

inline scvector_slice scvector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
  if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator()(const int, const int)"));
#endif
  return scvector_slice(*this,i,j);
}

inline scvector operator-(const scvector_slice& v) {
  return sl_v_negative<scvector_slice,scvector>(v);
}

inline int Lb(const scvector_slice& v) {
  return v.lb;
}

inline int Ub(const scvector_slice& v) {
  return v.ub;
}

inline srvector Re(const scvector_slice& v) {
  return Re(scvector(v));
}

inline srvector Im(const scvector_slice& v) {
  return Im(scvector(v));
}

inline int VecLen(const scvector_slice& v) {
  return v.n;
}

inline complex operator*(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_mult<scvector_slice,rvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_mult<srvector_slice,cvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_mult<scvector_slice,cvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<cvector,srvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<rvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<cvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<scvector_slice,rvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<srvector_slice,cvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<scvector_slice,cvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<cvector_slice,srvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<rvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<cvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<scvector,srvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<srvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<scvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<scvector_slice,srvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<srvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<scvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<scvector_slice,srvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<srvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline complex operator*(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<scvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

inline scvector operator*(const scvector_slice& v, const real& s) {
  return sp_vs_mult<scvector_slice,real,scvector>(v,s);
}

inline scvector operator*(const scvector_slice& v, const complex& s) {
  return sp_vs_mult<scvector_slice,complex,scvector>(v,s);
}

inline scvector operator*(const srvector_slice& v, const complex& s) {
  return sp_vs_mult<srvector_slice,complex,scvector>(v,s);
}

inline scvector operator/(const scvector_slice& v, const real& s) {
  return sp_vs_div<scvector_slice,real,scvector>(v,s);
}

inline scvector operator/(const scvector_slice& v, const complex& s) {
  return sp_vs_div<scvector_slice,complex,scvector>(v,s);
}

inline scvector operator/(const srvector_slice& v, const complex& s) {
  return sp_vs_div<srvector_slice,complex,scvector>(v,s);
}

inline scvector operator*(const real& s, const scvector_slice& v) {
  return sp_sv_mult<real,scvector_slice,scvector>(s,v);
}

inline scvector operator*(const complex& s, const scvector_slice& v) {
  return sp_sv_mult<complex,scvector_slice,scvector>(s,v);
}

inline scvector operator*(const complex& s, const srvector_slice& v) {
  return sp_sv_mult<complex,srvector_slice,scvector>(s,v);
}

inline cvector operator+(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_add<cvector,srvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_add<rvector,scvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_add<cvector,scvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_add<scvector_slice,rvector,cvector>(v1,v2);
}

inline cvector operator+(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_add<srvector_slice,cvector,cvector>(v1,v2);
}

inline cvector operator+(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_add<scvector_slice,cvector,cvector>(v1,v2);
}

inline cvector operator+(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<cvector_slice,srvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<rvector_slice,scvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<cvector_slice,scvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<scvector_slice,rvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<srvector_slice,cvector_slice,cvector>(v1,v2);
}

inline cvector operator+(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<scvector_slice,cvector_slice,cvector>(v1,v2);
}

inline scvector operator+(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<scvector_slice,srvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<srvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<scvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_add<scvector,srvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_add<srvector,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_add<scvector,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_add<scvector_slice,srvector,scvector,complex>(v1,v2);
}

inline scvector operator+(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_add<srvector_slice,scvector,scvector,complex>(v1,v2);
}

inline scvector operator+(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_add<scvector_slice,scvector,scvector,complex>(v1,v2);
}

inline cvector operator-(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<cvector,srvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<rvector,scvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<cvector,scvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_sub<scvector_slice,rvector,cvector>(v1,v2);
}

inline cvector operator-(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_sub<srvector_slice,cvector,cvector>(v1,v2);
}

inline cvector operator-(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_sub<scvector_slice,cvector,cvector>(v1,v2);
}

inline cvector operator-(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<cvector_slice,srvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<rvector_slice,scvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<cvector_slice,scvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<scvector_slice,rvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<srvector_slice,cvector_slice,cvector>(v1,v2);
}

inline cvector operator-(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<scvector_slice,cvector_slice,cvector>(v1,v2);
}

inline scvector operator-(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<scvector_slice,srvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<srvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<scvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<scvector,srvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<srvector,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<scvector,scvector_slice,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<scvector_slice,srvector,scvector,complex>(v1,v2);
}

inline scvector operator-(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<srvector_slice,scvector,scvector,complex>(v1,v2);
}

inline scvector operator-(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<scvector_slice,scvector,scvector,complex>(v1,v2);
}

inline cvector& cvector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector& cvector::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline scvector& scvector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline scvector& scvector::operator+=(const scvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline cvector& cvector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector& cvector::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline scvector& scvector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline scvector& scvector::operator-=(const scvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline bool operator==(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const scvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const scvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const cvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const cvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const scvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const scvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const scvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const scvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const cvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const cvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const scvector_slice& v) {
  return sl_v_output<scvector_slice,complex>(os,v);
}

inline std::istream& operator>>(std::istream& is, scvector_slice& v) {
  return sl_v_input<scvector_slice,complex>(is,v);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const scvector& y) {
  spsp_vv_accu<cdotprecision,scvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const srvector& y) {
  spsp_vv_accu<cdotprecision,scvector,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector& x, const scvector& y) {
  spsp_vv_accu<cdotprecision,srvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const cvector& y) {
  spf_vv_accu<cdotprecision,scvector,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const rvector& y) {
  spf_vv_accu<cdotprecision,scvector,rvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector& x, const cvector& y) {
  spf_vv_accu<cdotprecision,srvector,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const cvector_slice& y) {
  spf_vv_accu<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const rvector_slice& y) {
  spf_vv_accu<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector& x, const cvector_slice& y) {
  spf_vv_accu<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,cvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector& x, const srvector& y) {
  fsp_vv_accu<cdotprecision,cvector,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const rvector& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,rvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector_slice& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector_slice& x, const srvector& y) {
  fsp_vv_accu<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const rvector_slice& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const cvector& y) {
  slf_vv_accu<cdotprecision,scvector_slice,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const rvector& y) {
  slf_vv_accu<cdotprecision,scvector_slice,rvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const cvector& y) {
  slf_vv_accu<cdotprecision,srvector_slice,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  slf_vv_accu<cdotprecision,scvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  slf_vv_accu<cdotprecision,scvector_slice,rvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  slf_vv_accu<cdotprecision,srvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector& x, const srvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const rvector& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,rvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,rvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  slsl_vv_accu<cdotprecision,scvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<cdotprecision,scvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  slsl_vv_accu<cdotprecision,srvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const scvector_slice& y) {
  spsl_vv_accu<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector& x, const srvector_slice& y) {
  spsl_vv_accu<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector& x, const scvector_slice& y) {
  spsl_vv_accu<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const scvector& y) {
  slsp_vv_accu<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const scvector_slice& x, const srvector& y) {
  slsp_vv_accu<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cdotprecision& dot, const srvector_slice& x, const scvector& y) {
  slsp_vv_accu<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const scvector& y) {
  spsp_vv_accuapprox<cdotprecision,scvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const srvector& y) {
  spsp_vv_accuapprox<cdotprecision,scvector,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const scvector& y) {
  spsp_vv_accuapprox<cdotprecision,srvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const cvector& y) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const rvector& y) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const cvector& y) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const cvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const rvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const cvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector& x, const srvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const rvector& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,rvector,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const srvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const cvector& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const rvector& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,rvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const cvector& y) {
  slf_vv_accuapprox<cdotprecision,srvector_slice,cvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,rvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,srvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector& x, const srvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const rvector& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,rvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,rvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,scvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,scvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,srvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const scvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector& x, const srvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector& x, const scvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const scvector& y) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const srvector& y) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const scvector& y) {
  slsp_vv_accuapprox<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot,x,y);
}

inline void accumulate(cidotprecision& dot, const scvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const rvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const rvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const rvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const scvector_slice& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

} //namespace cxsc

#include "sparsevector.inl"

#endif
