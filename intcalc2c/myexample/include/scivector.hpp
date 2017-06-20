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

/* CVS $Id: scivector.hpp,v 1.8 2010/12/21 09:54:01 cxsc Exp $ */

#ifndef _CXSC_SCIVECTOR_HPP_INCLUDED
#define _CXSC_SCIVECTOR_HPP_INCLUDED

#include <cinterval.hpp>
#include <civector.hpp>
#include <vector>
#include <iostream>
#include <srvector.hpp>
#include <scvector.hpp>
#include <sivector.hpp>
#include <sparsecidot.hpp>
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
class sivector_slice;
class simatrix;
class simatrix_slice;
class simatrix_subv;
class scivector_slice;
class scimatrix_slice;
class scimatrix_subv;

class scivector {
  private:
    std::vector<int> p;
    std::vector<cinterval> x;
    int lb;
    int ub;
    int n; 

  public:
    scivector() : lb(0), ub(-1) , n(0) {
    }

    explicit scivector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    scivector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    scivector(const civector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    scivector(const cvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(cinterval(v[i]));
          }
        }
    }

    scivector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(cinterval(v[i]));
          }
        }
    }

    scivector(const int n, const int nnz, const intvector& index, const civector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i+Lb(values)] != 0.0) {
          p.push_back(index[i+Lb(index)]);
          x.push_back(values[i+Lb(values)]);
        }
      }
      
    }

    scivector(const srvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(cinterval(v.x[i]));
    }

    scivector(const scvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz()); 
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(cinterval(v.x[i]));
    }

    scivector(const sivector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(cinterval(v.x[i]));
    }

    scivector(const srvector_slice&);
    scivector(const scvector_slice&);
    scivector(const sivector_slice&);
    scivector(const scivector_slice&);
    scivector(const srmatrix_subv& A);
    scivector(const scmatrix_subv& A);
    scivector(const simatrix_subv& A);
    scivector(const scimatrix_subv& A);

    std::vector<int>& row_indices() {
      return p;
    }

    std::vector<cinterval>& values() {
      return x;
    }

    const std::vector<int>& row_indices() const {
      return p;
    }

    const std::vector<cinterval>& values() const {
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

    scivector& operator=(const srvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = cinterval(v.x[i]);
      return *this;
    } 

    scivector& operator=(const sivector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = cinterval(v.x[i]);
      return *this;
    } 

    scivector& operator=(const scvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = cinterval(v.x[i]);
      return *this;
    } 

    scivector& operator=(const real& v) {
      return sp_vs_assign<scivector,real,cinterval>(*this,v);
    }

    scivector& operator=(const complex& v) {
      return sp_vs_assign<scivector,complex,cinterval>(*this,v);
    }

    scivector& operator=(const interval& v) {
      return sp_vs_assign<scivector,interval,cinterval>(*this,v);
    }

    scivector& operator=(const cinterval& v) {
      return sp_vs_assign<scivector,cinterval,cinterval>(*this,v);
    }

    scivector& operator=(const rvector& v) {
      return spf_vv_assign<scivector,rvector,cinterval>(*this,v);
    }

    scivector& operator=(const cvector& v) {
      return spf_vv_assign<scivector,cvector,cinterval>(*this,v);
    }

    scivector& operator=(const ivector& v) {
      return spf_vv_assign<scivector,ivector,cinterval>(*this,v);
    }

    scivector& operator=(const civector& v) {
      return spf_vv_assign<scivector,civector,cinterval>(*this,v);
    }

    scivector& operator=(const rvector_slice& v) {
      return spf_vv_assign<scivector,rvector_slice,cinterval>(*this,v);
    }

    scivector& operator=(const cvector_slice& v) {
      return spf_vv_assign<scivector,cvector_slice,cinterval>(*this,v);
    }

    scivector& operator=(const ivector_slice& v) {
      return spf_vv_assign<scivector,ivector_slice,cinterval>(*this,v);
    }

    scivector& operator=(const civector_slice& v) {
      return spf_vv_assign<scivector,civector_slice,cinterval>(*this,v);
    }

    scivector& operator=(const srvector_slice&);
    scivector& operator=(const scvector_slice&);
    scivector& operator=(const sivector_slice&);
    scivector& operator=(const scivector_slice&);

    cinterval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, cinterval(0.0));

      return x[k];
    }

    cinterval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    cinterval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector::operator()(const int)"));
#endif
      cinterval r(0.0);

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    scivector operator()(const intvector& per) {
      scivector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,cinterval> work;
      for(int i=0 ; i<get_nnz() ; i++)
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
 
      for(std::map<int,cinterval>::iterator it=work.begin() ; it!=work.end() ; it++) {
         v.p.push_back(it->first);
         v.x.push_back(it->second);
      }

      return v;
    }

    scivector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }


    scivector_slice operator()(const int, const int);

    scivector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    scivector& operator*=(const complex& s) {
      return sp_vs_multassign(*this,s);
    }

    scivector& operator*=(const interval& s) {
      return sp_vs_multassign(*this,s);
    }

    scivector& operator*=(const cinterval& s) {
      return sp_vs_multassign(*this,s);
    }

    scivector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    scivector& operator/=(const complex& s) {
      return sp_vs_divassign(*this,s);
    }

    scivector& operator/=(const interval& s) {
      return sp_vs_divassign(*this,s);
    }

    scivector& operator/=(const cinterval& s) {
      return sp_vs_divassign(*this,s);
    }

    scivector& operator+=(const rvector& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const cvector& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const ivector& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const civector& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const cvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const ivector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const civector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    scivector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scivector& operator+=(const scvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scivector& operator+=(const sivector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scivector& operator+=(const scivector& v) {
      return spsp_vv_addassign(*this,v);
    }

    scivector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const cvector& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const ivector& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const civector& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const cvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const ivector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const civector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    scivector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scivector& operator-=(const scvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scivector& operator-=(const sivector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scivector& operator-=(const scivector& v) {
      return spsp_vv_subassign(*this,v);
    }

    scivector& operator|=(const rvector& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const cvector& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const ivector& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const civector& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const rvector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const cvector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const ivector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const civector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    scivector& operator|=(const srvector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    scivector& operator|=(const scvector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    scivector& operator|=(const sivector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    scivector& operator|=(const scivector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    scivector& operator&=(const ivector& v) {
      return spf_vv_intersectassign(*this,v);
    }

    scivector& operator&=(const civector& v) {
      return spf_vv_intersectassign(*this,v);
    }

    scivector& operator&=(const ivector_slice& v) {
      return spf_vv_intersectassign(*this,v);
    }

    scivector& operator&=(const civector_slice& v) {
      return spf_vv_intersectassign(*this,v);
    }

    scivector& operator&=(const sivector& v) {
      return spsp_vv_intersectassign(*this,v);
    }

    scivector& operator&=(const scivector& v) {
      return spsp_vv_intersectassign(*this,v);
    }

    scivector& operator+=(const srvector_slice&);
    scivector& operator+=(const scvector_slice&);
    scivector& operator+=(const sivector_slice&);
    scivector& operator+=(const scivector_slice&);
    scivector& operator-=(const srvector_slice&);
    scivector& operator-=(const scvector_slice&);
    scivector& operator-=(const sivector_slice&);
    scivector& operator-=(const scivector_slice&);
    scivector& operator|=(const srvector_slice&);
    scivector& operator|=(const scvector_slice&);
    scivector& operator|=(const sivector_slice&);
    scivector& operator|=(const scivector_slice&);
    scivector& operator&=(const sivector_slice&);
    scivector& operator&=(const scivector_slice&);

    friend void SetLb(scivector&, const int);
    friend void SetUb(scivector&, const int);
    friend int Lb(const scivector&);
    friend int Ub(const scivector&);
    friend sivector Re(const scivector&);
    friend sivector Im(const scivector&);
    friend scvector Inf(const scivector&);
    friend scvector Sup(const scivector&);
    friend srvector InfRe(const scivector&);
    friend srvector InfIm(const scivector&);
    friend srvector SupRe(const scivector&);
    friend srvector SupIm(const scivector&);
    friend scivector conj(const scivector&);
    friend scivector conj(const scivector_slice&);
    friend sivector abs(const scivector&);
    friend scvector mid(const scivector&);
    friend scvector diam(const scivector&);
    friend int VecLen(const scivector&);

    friend class srvector_slice;
    friend class scvector_slice;
    friend class sivector_slice;
    friend class scivector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline civector::civector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const scvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const sivector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const scivector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector& civector::operator=(const scivector& v) {
  return fsp_vv_assign<civector,scivector,cinterval>(*this,v);
}

inline civector& civector::operator=(const scivector_slice& v) {
  return fsl_vv_assign<civector,scivector_slice,cinterval>(*this,v);
}

inline civector& civector::operator=(const srvector& v) {
  return fsp_vv_assign<civector,srvector,cinterval>(*this,v);
}

inline civector& civector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<civector,srvector_slice,cinterval>(*this,v);
}

inline civector& civector::operator=(const scvector& v) {
  return fsp_vv_assign<civector,scvector,cinterval>(*this,v);
}

inline civector& civector::operator=(const scvector_slice& v) {
  return fsl_vv_assign<civector,scvector_slice,cinterval>(*this,v);
}

inline civector& civector::operator=(const sivector& v) {
  return fsp_vv_assign<civector,sivector,cinterval>(*this,v);
}

inline civector& civector::operator=(const sivector_slice& v) {
  return fsl_vv_assign<civector,sivector_slice,cinterval>(*this,v);
}

inline void SetLb(scivector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

inline void SetUb(scivector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

inline int Lb(const scivector& v) {
  return v.lb;
}

inline int Ub(const scivector& v) {
  return v.ub;
}

inline void Resize(scivector& v) {
  sp_v_resize(v);
}

inline void Resize(scivector& v, const int n) {
  sp_v_resize(v,n);
}

inline void Resize(scivector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

inline sivector abs(const scivector& v) {
  sivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

inline scivector conj(const scivector& v) {
  scivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(conj(v.x[i]));
  return res;
}

inline scvector mid(const scivector& v) {
  scvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(mid(v.x[i]));
  return res;
}

inline scvector diam(const scivector& v) {
  scvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

inline sivector Re(const scivector& v) {
  sivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Re(v.x[i]);
  return res;
}

inline sivector Im(const scivector& v) {
  sivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Im(v.x[i]);
  return res;
}

inline scvector Inf(const scivector& v) {
  scvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Inf(v.x[i]);
  return res;
}

inline scvector Sup(const scivector& v) {
  scvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Sup(v.x[i]);
  return res;
}

inline srvector InfRe(const scivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = InfRe(v.x[i]);
  return res;
}

inline srvector InfIm(const scivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = InfIm(v.x[i]);
  return res;
}

inline srvector SupRe(const scivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = SupRe(v.x[i]);
  return res;
}

inline srvector SupIm(const scivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = SupIm(v.x[i]);
  return res;
}

inline int VecLen(const scivector& v) {
  return v.n;
}

inline bool in (const scivector& v1, const scivector& v2) {
  for(int i=0 ; i<VecLen(v1) ; i++)
    if(!in(v1(i+Lb(v1)), v2(i+Lb(v2)))) return false;
  return true;
}

inline scivector operator-(const scivector& v) {
  return sp_v_negative(v);
}

inline cinterval operator*(const scivector& v1, const cvector& v2) {
  return spf_vv_mult<scivector,cvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const rvector& v2) {
  return spf_vv_mult<scivector,rvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const ivector& v2) {
  return spf_vv_mult<scivector,ivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const civector& v2) {
  return spf_vv_mult<scivector,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const civector& v2) {
  return spf_vv_mult<scvector,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector& v1, const civector& v2) {
  return spf_vv_mult<srvector,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const civector& v2) {
  return spf_vv_mult<sivector,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const ivector& v2) {
  return spf_vv_mult<scvector,ivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const cvector& v2) {
  return spf_vv_mult<sivector,cvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const rvector& v1, const scivector& v2) {
  return fsp_vv_mult<rvector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector& v1, const scivector& v2) {
  return fsp_vv_mult<cvector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector& v1, const scivector& v2) {
  return fsp_vv_mult<ivector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const scivector& v2) {
  return fsp_vv_mult<civector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const srvector& v2) {
  return fsp_vv_mult<civector,srvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const scvector& v2) {
  return fsp_vv_mult<civector,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const sivector& v2) {
  return fsp_vv_mult<civector,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector& v1, const scvector& v2) {
  return fsp_vv_mult<ivector,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector& v1, const sivector& v2) {
  return fsp_vv_mult<cvector,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_mult<scivector,cvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_mult<scivector,rvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_mult<scivector,ivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const civector_slice& v2) {
  return spf_vv_mult<scivector,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const civector_slice& v2) {
  return spf_vv_mult<scvector,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector& v1, const civector_slice& v2) {
  return spf_vv_mult<srvector,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const civector_slice& v2) {
  return spf_vv_mult<sivector,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const ivector_slice& v2) {
  return spf_vv_mult<scvector,ivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const cvector_slice& v2) {
  return spf_vv_mult<sivector,cvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_mult<rvector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_mult<cvector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_mult<ivector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_mult<civector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<civector_slice,srvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<civector_slice,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<civector_slice,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<ivector_slice,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<cvector_slice,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const srvector& v2) {
  return spsp_vv_mult<scivector,srvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const scvector& v2) {
  return spsp_vv_mult<scivector,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const sivector& v2) {
  return spsp_vv_mult<scivector,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const scivector& v2) {
  return spsp_vv_mult<scivector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector& v1, const scivector& v2) {
  return spsp_vv_mult<srvector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const scivector& v2) {
  return spsp_vv_mult<scvector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const scivector& v2) {
  return spsp_vv_mult<sivector,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const sivector& v2) {
  return spsp_vv_mult<scvector,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const scvector& v2) {
  return spsp_vv_mult<sivector,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline scivector operator*(const scivector& v, const real& s) {
  return sp_vs_mult<scivector,real,scivector>(v,s);
}

inline scivector operator*(const scivector& v, const complex& s) {
  return sp_vs_mult<scivector,complex,scivector>(v,s);
}

inline scivector operator*(const scivector& v, const interval& s) {
  return sp_vs_mult<scivector,interval,scivector>(v,s);
}

inline scivector operator*(const scivector& v, const cinterval& s) {
  return sp_vs_mult<scivector,cinterval,scivector>(v,s);
}

inline scivector operator*(const scvector& v, const interval& s) {
  return sp_vs_mult<scvector,interval,scivector>(v,s);
}

inline scivector operator*(const sivector& v, const complex& s) {
  return sp_vs_mult<sivector,complex,scivector>(v,s);
}

inline scivector operator/(const scivector& v, const real& s) {
  return sp_vs_div<scivector,real,scivector>(v,s);
}

inline scivector operator/(const scivector& v, const complex& s) {
  return sp_vs_div<scivector,complex,scivector>(v,s);
}

inline scivector operator/(const scivector& v, const interval& s) {
  return sp_vs_div<scivector,interval,scivector>(v,s);
}

inline scivector operator/(const scivector& v, const cinterval& s) {
  return sp_vs_div<scivector,cinterval,scivector>(v,s);
}

inline scivector operator/(const srvector& v, const cinterval& s) {
  return sp_vs_div<srvector,cinterval,scivector>(v,s);
}

inline scivector operator/(const scvector& v, const cinterval& s) {
  return sp_vs_div<scvector,cinterval,scivector>(v,s);
}

inline scivector operator/(const sivector& v, const cinterval& s) {
  return sp_vs_div<sivector,cinterval,scivector>(v,s);
}

inline scivector operator/(const scvector& v, const interval& s) {
  return sp_vs_div<scvector,interval,scivector>(v,s);
}

inline scivector operator/(const sivector& v, const complex& s) {
  return sp_vs_div<sivector,complex,scivector>(v,s);
}

inline scivector operator*(const real& s, const scivector& v) {
  return sp_sv_mult<real,scivector,scivector>(s,v);
}

inline scivector operator*(const complex& s, const scivector& v) {
  return sp_sv_mult<complex,scivector,scivector>(s,v);
}

inline scivector operator*(const interval& s, const scivector& v) {
  return sp_sv_mult<interval,scivector,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const srvector& v) {
  return sp_sv_mult<cinterval,srvector,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const sivector& v) {
  return sp_sv_mult<cinterval,sivector,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const scvector& v) {
  return sp_sv_mult<cinterval,scvector,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const scivector& v) {
  return sp_sv_mult<cinterval,scivector,scivector>(s,v);
}

inline scivector operator*(const srvector& v, const cinterval& s) {
  return sp_sv_mult<cinterval,srvector,scivector>(s,v);
}

inline scivector operator*(const sivector& v, const cinterval& s) {
  return sp_sv_mult<cinterval,sivector,scivector>(s,v);
}

inline scivector operator*(const scvector& v, const cinterval& s) {
  return sp_sv_mult<cinterval,scvector,scivector>(s,v);
}

inline scivector operator*(const complex& s, const sivector& v) {
  return sp_sv_mult<complex,sivector,scivector>(s,v);
}

inline scivector operator*(const interval& s, const scvector& v) {
  return sp_sv_mult<interval,scvector,scivector>(s,v);
}

inline civector operator+(const civector& v1, const srvector& v2) {
  return fsp_vv_add<civector,srvector,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const scvector& v2) {
  return fsp_vv_add<civector,scvector,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const sivector& v2) {
  return fsp_vv_add<civector,sivector,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const scivector& v2) {
  return fsp_vv_add<civector,scivector,civector>(v1,v2);
}

inline civector operator+(const rvector& v1, const scivector& v2) {
  return fsp_vv_add<rvector,scivector,civector>(v1,v2);
}

inline civector operator+(const cvector& v1, const scivector& v2) {
  return fsp_vv_add<cvector,scivector,civector>(v1,v2);
}

inline civector operator+(const ivector& v1, const scivector& v2) {
  return fsp_vv_add<ivector,scivector,civector>(v1,v2);
}

inline civector operator+(const cvector& v1, const sivector& v2) {
  return fsp_vv_add<cvector,sivector,civector>(v1,v2);
}

inline civector operator+(const ivector& v1, const scvector& v2) {
  return fsp_vv_add<ivector,scvector,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const rvector& v2) {
  return spf_vv_add<scivector,rvector,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const cvector& v2) {
  return spf_vv_add<scivector,cvector,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const ivector& v2) {
  return spf_vv_add<scivector,ivector,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const civector& v2) {
  return spf_vv_add<scivector,civector,civector>(v1,v2);
}

inline civector operator+(const srvector& v1, const civector& v2) {
  return spf_vv_add<srvector,civector,civector>(v1,v2);
}

inline civector operator+(const scvector& v1, const civector& v2) {
  return spf_vv_add<scvector,civector,civector>(v1,v2);
}

inline civector operator+(const sivector& v1, const civector& v2) {
  return spf_vv_add<sivector,civector,civector>(v1,v2);
}

inline civector operator+(const scvector& v1, const ivector& v2) {
  return spf_vv_add<scvector,ivector,civector>(v1,v2);
}

inline civector operator+(const sivector& v1, const cvector& v2) {
  return spf_vv_add<sivector,cvector,civector>(v1,v2);
}


inline civector operator+(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_add<civector_slice,srvector,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_add<civector_slice,scvector,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_add<civector_slice,sivector,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_add<civector_slice,scivector,civector>(v1,v2);
}

inline civector operator+(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_add<rvector_slice,scivector,civector>(v1,v2);
}

inline civector operator+(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_add<cvector_slice,scivector,civector>(v1,v2);
}

inline civector operator+(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_add<ivector_slice,scivector,civector>(v1,v2);
}

inline civector operator+(const cvector_slice& v1, const sivector& v2) {
  return fsp_vv_add<cvector_slice,sivector,civector>(v1,v2);
}

inline civector operator+(const ivector_slice& v1, const scvector& v2) {
  return fsp_vv_add<ivector_slice,scvector,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_add<scivector,rvector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_add<scivector,cvector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_add<scivector,ivector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector& v1, const civector_slice& v2) {
  return spf_vv_add<scivector,civector_slice,civector>(v1,v2);
}

inline civector operator+(const srvector& v1, const civector_slice& v2) {
  return spf_vv_add<srvector,civector_slice,civector>(v1,v2);
}

inline civector operator+(const scvector& v1, const civector_slice& v2) {
  return spf_vv_add<scvector,civector_slice,civector>(v1,v2);
}

inline civector operator+(const sivector& v1, const civector_slice& v2) {
  return spf_vv_add<sivector,civector_slice,civector>(v1,v2);
}

inline civector operator+(const scvector& v1, const ivector_slice& v2) {
  return spf_vv_add<scvector,ivector_slice,civector>(v1,v2);
}

inline civector operator+(const sivector& v1, const cvector_slice& v2) {
  return spf_vv_add<sivector,cvector_slice,civector>(v1,v2);
}

inline scivector operator+(const scivector& v1, const srvector& v2) {
  return spsp_vv_add<scivector,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const scvector& v2) {
  return spsp_vv_add<scivector,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const sivector& v2) {
  return spsp_vv_add<scivector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const scivector& v2) {
  return spsp_vv_add<scivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const srvector& v1, const scivector& v2) {
  return spsp_vv_add<srvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector& v1, const scivector& v2) {
  return spsp_vv_add<scvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector& v1, const scivector& v2) {
  return spsp_vv_add<sivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector& v1, const sivector& v2) {
  return spsp_vv_add<scvector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector& v1, const scvector& v2) {
  return spsp_vv_add<sivector,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator-(const civector& v1, const srvector& v2) {
  return fsp_vv_sub<civector,srvector,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const sivector& v2) {
  return fsp_vv_sub<civector,sivector,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const scvector& v2) {
  return fsp_vv_sub<civector,scvector,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const scivector& v2) {
  return fsp_vv_sub<civector,scivector,civector>(v1,v2);
}

inline civector operator-(const rvector& v1, const scivector& v2) {
  return fsp_vv_sub<rvector,scivector,civector>(v1,v2);
}

inline civector operator-(const cvector& v1, const scivector& v2) {
  return fsp_vv_sub<cvector,scivector,civector>(v1,v2);
}

inline civector operator-(const ivector& v1, const scivector& v2) {
  return fsp_vv_sub<ivector,scivector,civector>(v1,v2);
}

inline civector operator-(const cvector& v1, const sivector& v2) {
  return fsp_vv_sub<cvector,sivector,civector>(v1,v2);
}

inline civector operator-(const ivector& v1, const scvector& v2) {
  return fsp_vv_sub<ivector,scvector,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const rvector& v2) {
  return spf_vv_sub<scivector,rvector,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const cvector& v2) {
  return spf_vv_sub<scivector,cvector,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const ivector& v2) {
  return spf_vv_sub<scivector,ivector,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const civector& v2) {
  return spf_vv_sub<scivector,civector,civector>(v1,v2);
}

inline civector operator-(const srvector& v1, const civector& v2) {
  return spf_vv_sub<srvector,civector,civector>(v1,v2);
}

inline civector operator-(const scvector& v1, const civector& v2) {
  return spf_vv_sub<scvector,civector,civector>(v1,v2);
}

inline civector operator-(const sivector& v1, const civector& v2) {
  return spf_vv_sub<sivector,civector,civector>(v1,v2);
}

inline civector operator-(const scvector& v1, const ivector& v2) {
  return spf_vv_sub<scvector,ivector,civector>(v1,v2);
}

inline civector operator-(const sivector& v1, const cvector& v2) {
  return spf_vv_sub<sivector,cvector,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<civector_slice,srvector,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<civector_slice,sivector,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<civector_slice,scvector,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_sub<civector_slice,scivector,civector>(v1,v2);
}

inline civector operator-(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_sub<rvector_slice,scivector,civector>(v1,v2);
}

inline civector operator-(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_sub<cvector_slice,scivector,civector>(v1,v2);
}

inline civector operator-(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_sub<ivector_slice,scivector,civector>(v1,v2);
}

inline civector operator-(const cvector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<cvector_slice,sivector,civector>(v1,v2);
}

inline civector operator-(const ivector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<ivector_slice,scvector,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_sub<scivector,rvector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_sub<scivector,cvector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_sub<scivector,ivector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector& v1, const civector_slice& v2) {
  return spf_vv_sub<scivector,civector_slice,civector>(v1,v2);
}

inline civector operator-(const srvector& v1, const civector_slice& v2) {
  return spf_vv_sub<srvector,civector_slice,civector>(v1,v2);
}

inline civector operator-(const scvector& v1, const civector_slice& v2) {
  return spf_vv_sub<scvector,civector_slice,civector>(v1,v2);
}

inline civector operator-(const sivector& v1, const civector_slice& v2) {
  return spf_vv_sub<sivector,civector_slice,civector>(v1,v2);
}

inline civector operator-(const scvector& v1, const ivector_slice& v2) {
  return spf_vv_sub<scvector,ivector_slice,civector>(v1,v2);
}

inline civector operator-(const sivector& v1, const cvector_slice& v2) {
  return spf_vv_sub<sivector,cvector_slice,civector>(v1,v2);
}

inline scivector operator-(const scivector& v1, const srvector& v2) {
  return spsp_vv_sub<scivector,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const scvector& v2) {
  return spsp_vv_sub<scivector,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const sivector& v2) {
  return spsp_vv_sub<scivector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const scivector& v2) {
  return spsp_vv_sub<scivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const srvector& v1, const scivector& v2) {
  return spsp_vv_sub<srvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector& v1, const scivector& v2) {
  return spsp_vv_sub<scvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector& v1, const scivector& v2) {
  return spsp_vv_sub<sivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector& v1, const sivector& v2) {
  return spsp_vv_sub<scvector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector& v1, const scvector& v2) {
  return spsp_vv_sub<sivector,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator|(const civector& v1, const srvector& v2) {
  return fsp_vv_hull<civector,srvector,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const scvector& v2) {
  return fsp_vv_hull<civector,scvector,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const sivector& v2) {
  return fsp_vv_hull<civector,sivector,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const scivector& v2) {
  return fsp_vv_hull<civector,scivector,civector>(v1,v2);
}

inline civector operator|(const rvector& v1, const scivector& v2) {
  return fsp_vv_hull<rvector,scivector,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const scivector& v2) {
  return fsp_vv_hull<cvector,scivector,civector>(v1,v2);
}

inline civector operator|(const ivector& v1, const scivector& v2) {
  return fsp_vv_hull<ivector,scivector,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const sivector& v2) {
  return fsp_vv_hull<cvector,sivector,civector>(v1,v2);
}

inline civector operator|(const ivector& v1, const scvector& v2) {
  return fsp_vv_hull<ivector,scvector,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const rvector& v2) {
  return spf_vv_hull<scivector,rvector,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const cvector& v2) {
  return spf_vv_hull<scivector,cvector,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const ivector& v2) {
  return spf_vv_hull<scivector,ivector,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const civector& v2) {
  return spf_vv_hull<scivector,civector,civector>(v1,v2);
}

inline civector operator|(const srvector& v1, const civector& v2) {
  return spf_vv_hull<srvector,civector,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const civector& v2) {
  return spf_vv_hull<scvector,civector,civector>(v1,v2);
}

inline civector operator|(const sivector& v1, const civector& v2) {
  return spf_vv_hull<sivector,civector,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const ivector& v2) {
  return spf_vv_hull<scvector,ivector,civector>(v1,v2);
}

inline civector operator|(const sivector& v1, const cvector& v2) {
  return spf_vv_hull<sivector,cvector,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<civector_slice,srvector,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_hull<civector_slice,scvector,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<civector_slice,sivector,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_hull<civector_slice,scivector,civector>(v1,v2);
}

inline civector operator|(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_hull<rvector_slice,scivector,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_hull<cvector_slice,scivector,civector>(v1,v2);
}

inline civector operator|(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_hull<ivector_slice,scivector,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<cvector_slice,sivector,civector>(v1,v2);
}

inline civector operator|(const ivector_slice& v1, const scvector& v2) {
  return fsp_vv_hull<ivector_slice,scvector,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_hull<scivector,rvector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_hull<scivector,cvector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_hull<scivector,ivector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector& v1, const civector_slice& v2) {
  return spf_vv_hull<scivector,civector_slice,civector>(v1,v2);
}

inline civector operator|(const srvector& v1, const civector_slice& v2) {
  return spf_vv_hull<srvector,civector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const civector_slice& v2) {
  return spf_vv_hull<scvector,civector_slice,civector>(v1,v2);
}

inline civector operator|(const sivector& v1, const civector_slice& v2) {
  return spf_vv_hull<sivector,civector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const ivector_slice& v2) {
  return spf_vv_hull<scvector,ivector_slice,civector>(v1,v2);
}

inline civector operator|(const sivector& v1, const cvector_slice& v2) {
  return spf_vv_hull<sivector,cvector_slice,civector>(v1,v2);
}

inline scivector operator|(const scivector& v1, const srvector& v2) {
  return spsp_vv_hull<scivector,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const scvector& v2) {
  return spsp_vv_hull<scivector,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const sivector& v2) {
  return spsp_vv_hull<scivector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const scivector& v2) {
  return spsp_vv_hull<scivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector& v1, const scivector& v2) {
  return spsp_vv_hull<srvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const scivector& v2) {
  return spsp_vv_hull<scvector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector& v1, const scivector& v2) {
  return spsp_vv_hull<sivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const sivector& v2) {
  return spsp_vv_hull<scvector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector& v1, const scvector& v2) {
  return spsp_vv_hull<sivector,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator|(const cvector& v1, const srvector& v2) {
  return fsp_vv_hull<cvector,srvector,civector>(v1,v2);
}

inline civector operator|(const rvector& v1, const scvector& v2) {
  return fsp_vv_hull<rvector,scvector,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const scvector& v2) {
  return fsp_vv_hull<cvector,scvector,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const rvector& v2) {
  return spf_vv_hull<scvector,rvector,civector>(v1,v2);
}

inline civector operator|(const srvector& v1, const cvector& v2) {
  return spf_vv_hull<srvector,cvector,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const cvector& v2) {
  return spf_vv_hull<scvector,cvector,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<cvector_slice,srvector,civector>(v1,v2);
}

inline civector operator|(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_hull<rvector_slice,scvector,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_hull<cvector_slice,scvector,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_hull<scvector,rvector_slice,civector>(v1,v2);
}

inline civector operator|(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_hull<srvector,cvector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_hull<scvector,cvector_slice,civector>(v1,v2);
}

inline scivector operator|(const scvector& v1, const srvector& v2) {
  return spsp_vv_hull<scvector,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector& v1, const scvector& v2) {
  return spsp_vv_hull<srvector,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const scvector& v2) {
  return spsp_vv_hull<scvector,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator|(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<cvector,srvector_slice,civector>(v1,v2);
}

inline civector operator|(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_hull<rvector,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_hull<cvector,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_hull<scvector_slice,rvector,civector>(v1,v2);
}

inline civector operator|(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_hull<srvector_slice,cvector,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_hull<scvector_slice,cvector,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<cvector_slice,srvector_slice,civector>(v1,v2);
}

inline civector operator|(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_hull<rvector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_hull<cvector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<scvector_slice,rvector_slice,civector>(v1,v2);
}

inline civector operator|(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_hull<srvector_slice,cvector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_hull<scvector_slice,cvector_slice,civector>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<scvector_slice,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_hull<srvector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_hull<scvector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<scvector,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_hull<srvector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_hull<scvector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<scvector_slice,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_hull<srvector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_hull<scvector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator&(const civector& v1, const sivector& v2) {
  return fsp_vv_intersect<civector,sivector,civector>(v1,v2);
}

inline civector operator&(const civector& v1, const scivector& v2) {
  return fsp_vv_intersect<civector,scivector,civector>(v1,v2);
}

inline civector operator&(const ivector& v1, const scivector& v2) {
  return fsp_vv_intersect<ivector,scivector,civector>(v1,v2);
}

inline civector operator&(const scivector& v1, const ivector& v2) {
  return spf_vv_intersect<scivector,ivector,civector>(v1,v2);
}

inline civector operator&(const scivector& v1, const civector& v2) {
  return spf_vv_intersect<scivector,civector,civector>(v1,v2);
}

inline civector operator&(const sivector& v1, const civector& v2) {
  return spf_vv_intersect<sivector,civector,civector>(v1,v2);
}

inline civector operator&(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_intersect<civector_slice,sivector,civector>(v1,v2);
}

inline civector operator&(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_intersect<civector_slice,scivector,civector>(v1,v2);
}

inline civector operator&(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_intersect<ivector_slice,scivector,civector>(v1,v2);
}

inline civector operator&(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_intersect<scivector,ivector_slice,civector>(v1,v2);
}

inline civector operator&(const scivector& v1, const civector_slice& v2) {
  return spf_vv_intersect<scivector,civector_slice,civector>(v1,v2);
}

inline civector operator&(const sivector& v1, const civector_slice& v2) {
  return spf_vv_intersect<sivector,civector_slice,civector>(v1,v2);
}

inline scivector operator&(const scivector& v1, const sivector& v2) {
  return spsp_vv_intersect<scivector,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector& v1, const scivector& v2) {
  return spsp_vv_intersect<scivector,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const sivector& v1, const scivector& v2) {
  return spsp_vv_intersect<sivector,scivector,scivector,cinterval>(v1,v2);
}

inline civector& civector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const scivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const scivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline civector& civector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const scivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const scivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline civector& civector::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const scvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const scivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const scvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const scivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline civector& civector::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline civector& civector::operator&=(const scivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline civector_slice& civector_slice::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline civector_slice& civector_slice::operator&=(const scivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline bool operator==(const scivector& v1, const scivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const sivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const scvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const scivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const scivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const scivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const cvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const ivector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const civector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const civector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const civector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const civector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const civector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const civector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const civector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const civector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const scivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const sivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const scvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const scivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const scivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const scivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const cvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const ivector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const civector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const civector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const civector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const civector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const ivector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const cvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const civector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const civector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const civector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const civector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const scivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator<(const scivector& v1, const scivector& v2) {
  return spsp_vv_less<scivector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const sivector& v2) {
  return spsp_vv_less<scivector,sivector,cinterval>(v1,v2);
}

inline bool operator<(const srvector& v1, const scivector& v2) {
  return spsp_vv_less<srvector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scvector& v1, const scivector& v2) {
  return spsp_vv_less<scvector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const sivector& v1, const scivector& v2) {
  return spsp_vv_less<sivector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const ivector& v2) {
  return spf_vv_less<scivector,ivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const civector& v2) {
  return spf_vv_less<scivector,civector,cinterval>(v1,v2);
}

inline bool operator<(const srvector& v1, const civector& v2) {
  return spf_vv_less<srvector,civector,cinterval>(v1,v2);
}

inline bool operator<(const scvector& v1, const civector& v2) {
  return spf_vv_less<scvector,civector,cinterval>(v1,v2);
}

inline bool operator<(const sivector& v1, const civector& v2) {
  return spf_vv_less<sivector,civector,cinterval>(v1,v2);
}

inline bool operator<(const civector& v1, const sivector& v2) {
  return fsp_vv_less<civector,sivector,cinterval>(v1,v2);
}

inline bool operator<(const civector& v1, const scivector& v2) {
  return fsp_vv_less<civector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const rvector& v1, const scivector& v2) {
  return fsp_vv_less<rvector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const cvector& v1, const scivector& v2) {
  return fsp_vv_less<cvector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const ivector& v1, const scivector& v2) {
  return fsp_vv_less<ivector,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_less<scivector,ivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const srvector& v1, const civector_slice& v2) {
  return spf_vv_less<srvector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scvector& v1, const civector_slice& v2) {
  return spf_vv_less<scvector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<(const sivector& v1, const civector_slice& v2) {
  return spf_vv_less<sivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const civector_slice& v2) {
  return spf_vv_less<scivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_less<civector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator<(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_less<rvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_less<ivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_less<cvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_less<civector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const scivector& v2) {
  return spsp_vv_leq<scivector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const sivector& v2) {
  return spsp_vv_leq<scivector,sivector,cinterval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const scivector& v2) {
  return spsp_vv_leq<srvector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scvector& v1, const scivector& v2) {
  return spsp_vv_leq<scvector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const scivector& v2) {
  return spsp_vv_leq<sivector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const ivector& v2) {
  return spf_vv_leq<scivector,ivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const civector& v2) {
  return spf_vv_leq<scivector,civector,cinterval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const civector& v2) {
  return spf_vv_leq<srvector,civector,cinterval>(v1,v2);
}

inline bool operator<=(const scvector& v1, const civector& v2) {
  return spf_vv_leq<scvector,civector,cinterval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const civector& v2) {
  return spf_vv_leq<sivector,civector,cinterval>(v1,v2);
}

inline bool operator<=(const civector& v1, const sivector& v2) {
  return fsp_vv_leq<civector,sivector,cinterval>(v1,v2);
}

inline bool operator<=(const civector& v1, const scivector& v2) {
  return fsp_vv_leq<civector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const rvector& v1, const scivector& v2) {
  return fsp_vv_leq<rvector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const cvector& v1, const scivector& v2) {
  return fsp_vv_leq<cvector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const ivector& v1, const scivector& v2) {
  return fsp_vv_leq<ivector,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_leq<scivector,ivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const civector_slice& v2) {
  return spf_vv_leq<srvector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scvector& v1, const civector_slice& v2) {
  return spf_vv_leq<scvector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const civector_slice& v2) {
  return spf_vv_leq<sivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const civector_slice& v2) {
  return spf_vv_leq<scivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_leq<civector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator<=(const rvector_slice& v1, const scivector& v2) {
  return fsp_vv_leq<rvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_leq<ivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const cvector_slice& v1, const scivector& v2) {
  return fsp_vv_leq<cvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_leq<civector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const scivector& v2) {
  return spsp_vv_greater<scivector,scivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const srvector& v2) {
  return spsp_vv_greater<scivector,srvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const sivector& v2) {
  return spsp_vv_greater<scivector,sivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const scvector& v2) {
  return spsp_vv_greater<scivector,scvector,cinterval>(v1,v2);
}

inline bool operator>(const sivector& v1, const scivector& v2) {
  return spsp_vv_greater<sivector,scivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const rvector& v2) {
  return spf_vv_greater<scivector,rvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const cvector& v2) {
  return spf_vv_greater<scivector,cvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const ivector& v2) {
  return spf_vv_greater<scivector,ivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const civector& v2) {
  return spf_vv_greater<scivector,civector,cinterval>(v1,v2);
}

inline bool operator>(const sivector& v1, const civector& v2) {
  return spf_vv_greater<sivector,civector,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const srvector& v2) {
  return fsp_vv_greater<civector,srvector,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const scvector& v2) {
  return fsp_vv_greater<civector,scvector,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const sivector& v2) {
  return fsp_vv_greater<civector,sivector,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const scivector& v2) {
  return fsp_vv_greater<civector,scivector,cinterval>(v1,v2);
}

inline bool operator>(const ivector& v1, const scivector& v2) {
  return fsp_vv_greater<ivector,scivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_greater<scivector,rvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_greater<scivector,ivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_greater<scivector,cvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const sivector& v1, const civector_slice& v2) {
  return spf_vv_greater<sivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const civector_slice& v2) {
  return spf_vv_greater<scivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_greater<civector_slice,srvector,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_greater<civector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_greater<civector_slice,scvector,cinterval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_greater<ivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_greater<civector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const scivector& v2) {
  return spsp_vv_geq<scivector,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const srvector& v2) {
  return spsp_vv_geq<scivector,srvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const sivector& v2) {
  return spsp_vv_geq<scivector,sivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const scvector& v2) {
  return spsp_vv_geq<scivector,scvector,cinterval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const scivector& v2) {
  return spsp_vv_geq<sivector,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const rvector& v2) {
  return spf_vv_geq<scivector,rvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const cvector& v2) {
  return spf_vv_geq<scivector,cvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const ivector& v2) {
  return spf_vv_geq<scivector,ivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const civector& v2) {
  return spf_vv_geq<scivector,civector,cinterval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const civector& v2) {
  return spf_vv_geq<sivector,civector,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const srvector& v2) {
  return fsp_vv_geq<civector,srvector,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const scvector& v2) {
  return fsp_vv_geq<civector,scvector,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const sivector& v2) {
  return fsp_vv_geq<civector,sivector,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const scivector& v2) {
  return fsp_vv_geq<civector,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const scivector& v2) {
  return fsp_vv_geq<ivector,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const rvector_slice& v2) {
  return spf_vv_geq<scivector,rvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const ivector_slice& v2) {
  return spf_vv_geq<scivector,ivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const cvector_slice& v2) {
  return spf_vv_geq<scivector,cvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const civector_slice& v2) {
  return spf_vv_geq<sivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const civector_slice& v2) {
  return spf_vv_geq<scivector,civector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const srvector& v2) {
  return fsp_vv_geq<civector_slice,srvector,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const sivector& v2) {
  return fsp_vv_geq<civector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const scvector& v2) {
  return fsp_vv_geq<civector,scvector,cinterval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const scivector& v2) {
  return fsp_vv_geq<ivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const scivector& v2) {
  return fsp_vv_geq<civector_slice,scivector,cinterval>(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const scivector& v) {
  return sp_v_output<scivector,cinterval>(os,v);
}

inline std::istream& operator>>(std::istream& is, scivector& v) {
  return sp_v_input<scivector,cinterval>(is,v);
}


class scivector_slice {
  private:
    std::vector<int>& p;
    std::vector<cinterval>& x;
    scivector& orig;
    int start,end;
    int lb;
    int ub;
    int n;
    int nnz;
    int offset;

    scivector_slice(scivector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1) {
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

    cinterval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, cinterval(0.0));
      end++;

      return x[k];
    }

    cinterval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    cinterval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector_slice::operator()(const int)"));
#endif
      cinterval r(0.0);

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    scivector_slice& operator=(const real& v) {
      return sl_vs_assign<scivector_slice,real,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const complex& v) {
      return sl_vs_assign<scivector_slice,complex,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const interval& v) {
      return sl_vs_assign<scivector_slice,interval,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const cinterval& v) {
      return sl_vs_assign<scivector_slice,cinterval,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<scivector_slice,srvector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const scvector_slice& v) {
      return slsl_vv_assign<scivector_slice,scvector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const sivector_slice& v) {
      return slsl_vv_assign<scivector_slice,sivector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const scivector_slice& v) {
      return slsl_vv_assign<scivector_slice,scivector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<scivector_slice,srvector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const scvector& v) {
      return slsp_vv_assign<scivector_slice,scvector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const sivector& v) {
      return slsp_vv_assign<scivector_slice,sivector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const scivector& v) {
      return slsp_vv_assign<scivector_slice,scivector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const rvector& v) {
      return slf_vv_assign<scivector_slice,rvector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const cvector& v) {
      return slf_vv_assign<scivector_slice,cvector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const ivector& v) {
      return slf_vv_assign<scivector_slice,ivector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const civector& v) {
      return slf_vv_assign<scivector_slice,civector,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<scivector_slice,rvector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const cvector_slice& v) {
      return slf_vv_assign<scivector_slice,cvector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const ivector_slice& v) {
      return slf_vv_assign<scivector_slice,ivector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator=(const civector_slice& v) {
      return slf_vv_assign<scivector_slice,civector_slice,cinterval,std::vector<cinterval>::iterator>(*this,v);
    }

    scivector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    scivector_slice& operator*=(const complex& s) {
      return sl_vs_multassign(*this,s);
    }

    scivector_slice& operator*=(const interval& s) {
      return sl_vs_multassign(*this,s);
    }

    scivector_slice& operator*=(const cinterval& s) {
      return sl_vs_multassign(*this,s);
    }

    scivector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    scivector_slice& operator/=(const complex& s) {
      return sl_vs_divassign(*this,s);
    }

    scivector_slice& operator/=(const interval& s) {
      return sl_vs_divassign(*this,s);
    }

    scivector_slice& operator/=(const cinterval& s) {
      return sl_vs_divassign(*this,s);
    }

    scivector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<scivector_slice,rvector,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const ivector& v) {
      return slf_vv_addassign<scivector_slice,ivector,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const cvector& v) {
      return slf_vv_addassign<scivector_slice,cvector,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const civector& v) {
      return slf_vv_addassign<scivector_slice,civector,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<scivector_slice,rvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const cvector_slice& v) {
      return slf_vv_addassign<scivector_slice,cvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const ivector_slice& v) {
      return slf_vv_addassign<scivector_slice,ivector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const civector_slice& v) {
      return slf_vv_addassign<scivector_slice,civector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const scvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const sivector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const scivector& v) {
      return slsp_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const scvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const sivector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scivector_slice& operator+=(const scivector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    scivector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<scivector_slice,rvector,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const ivector& v) {
      return slf_vv_subassign<scivector_slice,ivector,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const cvector& v) {
      return slf_vv_subassign<scivector_slice,cvector,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const civector& v) {
      return slf_vv_subassign<scivector_slice,civector,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<scivector_slice,rvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const cvector_slice& v) {
      return slf_vv_subassign<scivector_slice,cvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const ivector_slice& v) {
      return slf_vv_subassign<scivector_slice,ivector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const civector_slice& v) {
      return slf_vv_subassign<scivector_slice,civector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const scvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const sivector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const scivector& v) {
      return slsp_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const scvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const sivector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    scivector_slice& operator-=(const scivector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    scivector_slice& operator|=(const rvector& v) {
      return slf_vv_hullassign<scivector_slice,rvector,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const ivector& v) {
      return slf_vv_hullassign<scivector_slice,ivector,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const cvector& v) {
      return slf_vv_hullassign<scivector_slice,cvector,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const civector& v) {
      return slf_vv_hullassign<scivector_slice,civector,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const rvector_slice& v) {
      return slf_vv_hullassign<scivector_slice,rvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const cvector_slice& v) {
      return slf_vv_hullassign<scivector_slice,cvector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const ivector_slice& v) {
      return slf_vv_hullassign<scivector_slice,ivector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const civector_slice& v) {
      return slf_vv_hullassign<scivector_slice,civector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator|=(const srvector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const scvector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const sivector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const scivector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const srvector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const scvector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const sivector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    scivector_slice& operator|=(const scivector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    scivector_slice& operator&=(const ivector& v) {
      return slf_vv_intersectassign<scivector_slice,ivector,cinterval>(*this,v);
    }

    scivector_slice& operator&=(const civector& v) {
      return slf_vv_intersectassign<scivector_slice,civector,cinterval>(*this,v);
    }

    scivector_slice& operator&=(const ivector_slice& v) {
      return slf_vv_intersectassign<scivector_slice,ivector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator&=(const civector_slice& v) {
      return slf_vv_intersectassign<scivector_slice,civector_slice,cinterval>(*this,v);
    }

    scivector_slice& operator&=(const sivector& v) {
      return slsp_vv_intersectassign(*this,v);
    }

    scivector_slice& operator&=(const scivector& v) {
      return slsp_vv_intersectassign(*this,v);
    }

    scivector_slice& operator&=(const sivector_slice& v) {
      return slsl_vv_intersectassign(*this,v);
    }

    scivector_slice& operator&=(const scivector_slice& v) {
      return slsl_vv_intersectassign(*this,v);
    }

    friend int Lb(const scivector_slice&);
    friend int Ub(const scivector_slice&);
    friend sivector Re(const scivector_slice&);
    friend sivector Im(const scivector_slice&);
    friend scvector Inf(const scivector_slice&);
    friend scvector Sup(const scivector_slice&);
    friend srvector InfRe(const scivector_slice&);
    friend srvector InfIm(const scivector_slice&);
    friend srvector SupRe(const scivector_slice&);
    friend srvector SupIm(const scivector_slice&);
    friend sivector abs(const scivector_slice&);
    friend scivector conj(const scivector_slice&);
    friend scvector mid(const scivector_slice&);
    friend scvector diam(const scivector_slice&);
    friend int VecLen(const scivector_slice&);

//     friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
//     friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class scvector;
    friend class sivector;
    friend class scivector;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline civector::civector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const scvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const sivector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector::civector(const scivector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new cinterval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline civector_slice& civector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const scvector& v) {
  *this = cvector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const scvector_slice& v) {
  *this = cvector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const scivector& v) {
  *this = civector(v);
  return *this;
}

inline civector_slice& civector_slice::operator=(const scivector_slice& v) {
  *this = civector(v);
  return *this;
}

inline scivector::scivector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(cinterval(s.x[i]));
  }

}

inline scivector::scivector(const scvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(cinterval(s.x[i]));
  }

}

inline scivector::scivector(const sivector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(cinterval(s.x[i]));
  }

}

inline scivector::scivector(const scivector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline scivector& scivector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<scivector,srvector_slice,cinterval>(*this,v);
}

inline scivector& scivector::operator=(const scvector_slice& v) {
  return spsl_vv_assign<scivector,scvector_slice,cinterval>(*this,v);
}

inline scivector& scivector::operator=(const sivector_slice& v) {
  return spsl_vv_assign<scivector,sivector_slice,cinterval>(*this,v);
}

inline scivector& scivector::operator=(const scivector_slice& v) {
  return spsl_vv_assign<scivector,scivector_slice,cinterval>(*this,v);
}

inline scivector_slice scivector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
  if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scivector_slice::operator()(const int, const int)"));
#endif
  return scivector_slice(*this,i,j);
}

inline scivector operator-(const scivector_slice& v) {
  return sl_v_negative<scivector_slice,scivector>(v);
}

inline int Lb(const scivector_slice& v) {
  return v.lb;
}

inline int Ub(const scivector_slice& v) {
  return v.ub;
}

inline sivector Re(const scivector_slice& v) {
  return Re(scivector(v));
}

inline sivector Im(const scivector_slice& v) {
  return Im(scivector(v));
}

inline scvector Inf(const scivector_slice& v) {
  return Inf(scivector(v));
}

inline scvector Sup(const scivector_slice& v) {
  return Sup(scivector(v));
}

inline srvector InfRe(const scivector_slice& v) {
  return InfRe(scivector(v));
}

inline srvector InfIm(const scivector_slice& v) {
  return InfIm(scivector(v));
}

inline srvector SupRe(const scivector_slice& v) {
  return SupRe(scivector(v));
}

inline srvector SupIm(const scivector_slice& v) {
  return SupIm(scivector(v));
}

inline scivector conj(const scivector_slice& v) {
  scivector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(conj(v.x[i]));
  return res;
}

inline sivector abs(const scivector_slice& v) {
  sivector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

inline scvector mid(const scivector_slice& v) {
  scvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(mid(v.x[i]));
  return res;
}

inline scvector diam(const scivector_slice& v) {
  scvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<v.end ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

inline int VecLen(const scivector_slice& v) {
  return v.n;
}

inline cinterval operator*(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_mult<scivector_slice,rvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_mult<scivector_slice,cvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_mult<scivector_slice,ivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const civector& v2) {
  return slf_vv_mult<scivector_slice,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector_slice& v1, const civector& v2) {
  return slf_vv_mult<srvector_slice,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const civector& v2) {
  return slf_vv_mult<sivector_slice,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const civector& v2) {
  return slf_vv_mult<scvector_slice,civector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const ivector& v2) {
  return slf_vv_mult<scvector_slice,ivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const cvector& v2) {
  return slf_vv_mult<sivector_slice,cvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<civector,srvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<civector,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<civector,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_mult<civector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_mult<rvector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_mult<cvector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_mult<ivector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<cvector,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<ivector,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<scivector_slice,rvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<scivector_slice,ivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<scivector_slice,cvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_mult<scivector_slice,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector_slice& v1, const civector_slice& v2) {
  return slf_vv_mult<srvector_slice,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const civector_slice& v2) {
  return slf_vv_mult<scvector_slice,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_mult<sivector_slice,civector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<scvector_slice,ivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<sivector_slice,cvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<civector_slice,srvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<civector_slice,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<civector_slice,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_mult<civector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const rvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_mult<rvector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_mult<ivector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_mult<cvector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const cvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<cvector_slice,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const ivector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<ivector_slice,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<scivector,srvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<scivector,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<scivector,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_mult<scivector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_mult<srvector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_mult<scvector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_mult<sivector,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<scvector,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<sivector,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<scivector_slice,srvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<scivector_slice,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<scivector_slice,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_mult<scivector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_mult<srvector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_mult<sivector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_mult<scvector_slice,scivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<scvector_slice,sivector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<sivector_slice,scvector,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<scivector_slice,srvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<scivector_slice,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<scivector_slice,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_mult<scivector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_mult<srvector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_mult<scvector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_mult<sivector_slice,scivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const sivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<sivector_slice,scvector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline cinterval operator*(const scvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<scvector_slice,sivector_slice,cinterval,sparse_cidot>(v1,v2);
}

inline scivector operator*(const scivector_slice& v, const real& s) {
  return sp_vs_mult<scivector_slice,real,scivector>(v,s);
}

inline scivector operator*(const scivector_slice& v, const complex& s) {
  return sp_vs_mult<scivector_slice,complex,scivector>(v,s);
}

inline scivector operator*(const scivector_slice& v, const interval& s) {
  return sp_vs_mult<scivector_slice,interval,scivector>(v,s);
}

inline scivector operator*(const scivector_slice& v, const cinterval& s) {
  return sp_vs_mult<scivector_slice,cinterval,scivector>(v,s);
}

inline scivector operator*(const srvector_slice& v, const cinterval& s) {
  return sp_vs_mult<srvector_slice,cinterval,scivector>(v,s);
}

inline scivector operator*(const scvector_slice& v, const cinterval& s) {
  return sp_vs_mult<scvector_slice,cinterval,scivector>(v,s);
}

inline scivector operator*(const sivector_slice& v, const cinterval& s) {
  return sp_vs_mult<sivector_slice,cinterval,scivector>(v,s);
}

inline scivector operator*(const scvector_slice& v, const interval& s) {
  return sp_vs_mult<scvector_slice,interval,scivector>(v,s);
}

inline scivector operator*(const sivector_slice& v, const complex& s) {
  return sp_vs_mult<sivector_slice,complex,scivector>(v,s);
}

inline scivector operator/(const scivector_slice& v, const real& s) {
  return sp_vs_div<scivector_slice,real,scivector>(v,s);
}

inline scivector operator/(const scivector_slice& v, const complex& s) {
  return sp_vs_div<scivector_slice,complex,scivector>(v,s);
}

inline scivector operator/(const scivector_slice& v, const interval& s) {
  return sp_vs_div<scivector_slice,interval,scivector>(v,s);
}

inline scivector operator/(const scivector_slice& v, const cinterval& s) {
  return sp_vs_div<scivector_slice,cinterval,scivector>(v,s);
}

inline scivector operator/(const srvector_slice& v, const cinterval& s) {
  return sp_vs_div<srvector_slice,cinterval,scivector>(v,s);
}

inline scivector operator/(const scvector_slice& v, const cinterval& s) {
  return sp_vs_div<scvector_slice,cinterval,scivector>(v,s);
}

inline scivector operator/(const sivector_slice& v, const cinterval& s) {
  return sp_vs_div<sivector_slice,cinterval,scivector>(v,s);
}

inline scivector operator/(const scvector_slice& v, const interval& s) {
  return sp_vs_div<scvector_slice,interval,scivector>(v,s);
}

inline scivector operator/(const sivector_slice& v, const complex& s) {
  return sp_vs_div<sivector_slice,complex,scivector>(v,s);
}


inline scivector operator*(const real& s, const scivector_slice& v) {
  return sp_sv_mult<real,scivector_slice,scivector>(s,v);
}

inline scivector operator*(const complex& s, const scivector_slice& v) {
  return sp_sv_mult<complex,scivector_slice,scivector>(s,v);
}

inline scivector operator*(const interval& s, const scivector_slice& v) {
  return sp_sv_mult<interval,scivector_slice,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const scivector_slice& v) {
  return sp_sv_mult<cinterval,scivector_slice,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const srvector_slice& v) {
  return sp_sv_mult<cinterval,srvector_slice,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const scvector_slice& v) {
  return sp_sv_mult<cinterval,scvector_slice,scivector>(s,v);
}

inline scivector operator*(const cinterval& s, const sivector_slice& v) {
  return sp_sv_mult<cinterval,sivector_slice,scivector>(s,v);
}

inline scivector operator*(const complex& s, const sivector_slice& v) {
  return sp_sv_mult<complex,sivector_slice,scivector>(s,v);
}

inline scivector operator*(const interval& s, const scvector_slice& v) {
  return sp_sv_mult<interval,scvector_slice,scivector>(s,v);
}

inline civector operator+(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_add<civector,srvector_slice,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_add<civector,scvector_slice,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_add<civector,sivector_slice,civector>(v1,v2);
}

inline civector operator+(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_add<civector,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_add<rvector,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_add<cvector,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_add<ivector,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const cvector& v1, const sivector_slice& v2) {
  return fsl_vv_add<cvector,sivector_slice,civector>(v1,v2);
}

inline civector operator+(const ivector& v1, const scvector_slice& v2) {
  return fsl_vv_add<ivector,scvector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_add<scivector_slice,rvector,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_add<scivector_slice,cvector,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_add<scivector_slice,ivector,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const civector& v2) {
  return slf_vv_add<scivector_slice,civector,civector>(v1,v2);
}

inline civector operator+(const srvector_slice& v1, const civector& v2) {
  return slf_vv_add<srvector_slice,civector,civector>(v1,v2);
}

inline civector operator+(const scvector_slice& v1, const civector& v2) {
  return slf_vv_add<scvector_slice,civector,civector>(v1,v2);
}

inline civector operator+(const sivector_slice& v1, const civector& v2) {
  return slf_vv_add<sivector_slice,civector,civector>(v1,v2);
}

inline civector operator+(const scvector_slice& v1, const ivector& v2) {
  return slf_vv_add<scvector_slice,ivector,civector>(v1,v2);
}

inline civector operator+(const sivector_slice& v1, const cvector& v2) {
  return slf_vv_add<sivector_slice,cvector,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<civector_slice,srvector_slice,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<civector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<civector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator+(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_add<civector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const rvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_add<rvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const cvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_add<cvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_add<ivector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator+(const ivector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<ivector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator+(const cvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<cvector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<scivector_slice,rvector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<scivector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<scivector_slice,cvector_slice,civector>(v1,v2);
}

inline civector operator+(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_add<scivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator+(const srvector_slice& v1, const civector_slice& v2) {
  return slf_vv_add<srvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator+(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_add<sivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator+(const scvector_slice& v1, const civector_slice& v2) {
  return slf_vv_add<scvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator+(const scvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<scvector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator+(const sivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<sivector_slice,cvector_slice,civector>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<scivector_slice,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<scivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<scivector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_add<scivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_add<srvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_add<scvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_add<sivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<scvector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<sivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_add<scivector,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_add<scivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_add<scivector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_add<scivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_add<srvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_add<scvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_add<sivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector& v1, const sivector_slice& v2) {
  return spsl_vv_add<scvector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector& v1, const scvector_slice& v2) {
  return spsl_vv_add<sivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_add<scivector_slice,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_add<scivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_add<scivector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_add<scivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_add<srvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_add<scvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_add<sivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const scvector_slice& v1, const sivector& v2) {
  return slsp_vv_add<scvector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator+(const sivector_slice& v1, const scvector& v2) {
  return slsp_vv_add<sivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator-(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<civector,srvector_slice,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<civector,scvector_slice,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<civector,sivector_slice,civector>(v1,v2);
}

inline civector operator-(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_sub<civector,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_sub<rvector,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_sub<cvector,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_sub<ivector,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const cvector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<cvector,sivector_slice,civector>(v1,v2);
}

inline civector operator-(const ivector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<ivector,scvector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_sub<scivector_slice,rvector,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_sub<scivector_slice,cvector,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_sub<scivector_slice,ivector,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const civector& v2) {
  return slf_vv_sub<scivector_slice,civector,civector>(v1,v2);
}

inline civector operator-(const srvector_slice& v1, const civector& v2) {
  return slf_vv_sub<srvector_slice,civector,civector>(v1,v2);
}

inline civector operator-(const scvector_slice& v1, const civector& v2) {
  return slf_vv_sub<scvector_slice,civector,civector>(v1,v2);
}

inline civector operator-(const sivector_slice& v1, const civector& v2) {
  return slf_vv_sub<sivector_slice,civector,civector>(v1,v2);
}

inline civector operator-(const scvector_slice& v1, const ivector& v2) {
  return slf_vv_sub<scvector_slice,ivector,civector>(v1,v2);
}

inline civector operator-(const sivector_slice& v1, const cvector& v2) {
  return slf_vv_sub<sivector_slice,cvector,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<civector_slice,srvector_slice,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<civector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<civector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator-(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_sub<civector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const rvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_sub<rvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const cvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_sub<cvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_sub<ivector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator-(const ivector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<ivector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator-(const cvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<cvector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<scivector_slice,rvector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<scivector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<scivector_slice,cvector_slice,civector>(v1,v2);
}

inline civector operator-(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_sub<scivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator-(const srvector_slice& v1, const civector_slice& v2) {
  return slf_vv_sub<srvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator-(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_sub<sivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator-(const scvector_slice& v1, const civector_slice& v2) {
  return slf_vv_sub<scvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator-(const scvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<scvector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator-(const sivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<sivector_slice,cvector_slice,civector>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<scivector_slice,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<scivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<scivector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_sub<scivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_sub<srvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_sub<scvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_sub<sivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<scvector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<sivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<scivector,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<scivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<scivector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_sub<scivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_sub<srvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_sub<scvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_sub<sivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<scvector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<sivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<scivector_slice,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<scivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<scivector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_sub<scivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_sub<srvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_sub<scvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_sub<sivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const scvector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<scvector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator-(const sivector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<sivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator|(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<civector,srvector_slice,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_hull<civector,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<civector,sivector_slice,civector>(v1,v2);
}

inline civector operator|(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_hull<civector,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_hull<rvector,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_hull<cvector,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_hull<ivector,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<cvector,sivector_slice,civector>(v1,v2);
}

inline civector operator|(const ivector& v1, const scvector_slice& v2) {
  return fsl_vv_hull<ivector,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_hull<scivector_slice,rvector,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_hull<scivector_slice,cvector,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_hull<scivector_slice,ivector,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const civector& v2) {
  return slf_vv_hull<scivector_slice,civector,civector>(v1,v2);
}

inline civector operator|(const srvector_slice& v1, const civector& v2) {
  return slf_vv_hull<srvector_slice,civector,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const civector& v2) {
  return slf_vv_hull<scvector_slice,civector,civector>(v1,v2);
}

inline civector operator|(const sivector_slice& v1, const civector& v2) {
  return slf_vv_hull<sivector_slice,civector,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const ivector& v2) {
  return slf_vv_hull<scvector_slice,ivector,civector>(v1,v2);
}

inline civector operator|(const sivector_slice& v1, const cvector& v2) {
  return slf_vv_hull<sivector_slice,cvector,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<civector_slice,srvector_slice,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_hull<civector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<civector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator|(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_hull<civector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const rvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_hull<rvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_hull<cvector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_hull<ivector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator|(const ivector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_hull<ivector_slice,scvector_slice,civector>(v1,v2);
}

inline civector operator|(const cvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<cvector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<scivector_slice,rvector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<scivector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_hull<scivector_slice,cvector_slice,civector>(v1,v2);
}

inline civector operator|(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_hull<scivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator|(const srvector_slice& v1, const civector_slice& v2) {
  return slf_vv_hull<srvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator|(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_hull<sivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const civector_slice& v2) {
  return slf_vv_hull<scvector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator|(const scvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<scvector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator|(const sivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_hull<sivector_slice,cvector_slice,civector>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<scivector_slice,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_hull<scivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<scivector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_hull<scivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_hull<srvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_hull<scvector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_hull<sivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<scvector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_hull<sivector_slice,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<scivector,srvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_hull<scivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<scivector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_hull<scivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_hull<srvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_hull<scvector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_hull<sivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<scvector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector& v1, const scvector_slice& v2) {
  return spsl_vv_hull<sivector,scvector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<scivector_slice,srvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_hull<scivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<scivector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_hull<scivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_hull<srvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_hull<scvector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_hull<sivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const scvector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<scvector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator|(const sivector_slice& v1, const scvector& v2) {
  return slsp_vv_hull<sivector_slice,scvector,scivector,cinterval>(v1,v2);
}

inline civector operator&(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<civector,sivector_slice,civector>(v1,v2);
}

inline civector operator&(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_intersect<civector,scivector_slice,civector>(v1,v2);
}

inline civector operator&(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_intersect<ivector,scivector_slice,civector>(v1,v2);
}

inline civector operator&(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_intersect<scivector_slice,ivector,civector>(v1,v2);
}

inline civector operator&(const scivector_slice& v1, const civector& v2) {
  return slf_vv_intersect<scivector_slice,civector,civector>(v1,v2);
}

inline civector operator&(const sivector_slice& v1, const civector& v2) {
  return slf_vv_intersect<sivector_slice,civector,civector>(v1,v2);
}

inline civector operator&(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<civector_slice,sivector_slice,civector>(v1,v2);
}

inline civector operator&(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_intersect<civector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator&(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_intersect<ivector_slice,scivector_slice,civector>(v1,v2);
}

inline civector operator&(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_intersect<scivector_slice,ivector_slice,civector>(v1,v2);
}

inline civector operator&(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_intersect<scivector_slice,civector_slice,civector>(v1,v2);
}

inline civector operator&(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_intersect<sivector_slice,civector_slice,civector>(v1,v2);
}

inline scivector operator&(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_intersect<scivector_slice,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_intersect<scivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_intersect<sivector_slice,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_intersect<scivector,sivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_intersect<scivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_intersect<sivector,scivector_slice,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_intersect<scivector_slice,sivector,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_intersect<scivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline scivector operator&(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_intersect<sivector_slice,scivector,scivector,cinterval>(v1,v2);
}

inline civector& civector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector& civector::operator+=(const scivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline civector_slice& civector_slice::operator+=(const scivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline scivector& scivector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline scivector& scivector::operator+=(const scvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline scivector& scivector::operator+=(const sivector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline scivector& scivector::operator+=(const scivector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline civector& civector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector& civector::operator-=(const scivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline civector_slice& civector_slice::operator-=(const scivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline scivector& scivector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline scivector& scivector::operator-=(const scvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline scivector& scivector::operator-=(const sivector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline scivector& scivector::operator-=(const scivector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline civector& civector::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const scvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector& civector::operator|=(const scivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const scvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline civector_slice& civector_slice::operator|=(const scivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline scivector& scivector::operator|=(const srvector_slice& v2) {
  return spsl_vv_hullassign(*this,v2);
}

inline scivector& scivector::operator|=(const scvector_slice& v2) {
  return spsl_vv_hullassign(*this,v2);
}

inline scivector& scivector::operator|=(const sivector_slice& v2) {
  return spsl_vv_hullassign(*this,v2);
}

inline scivector& scivector::operator|=(const scivector_slice& v2) {
  return spsl_vv_hullassign(*this,v2);
}

inline civector& civector::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline civector& civector::operator&=(const scivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline civector_slice& civector_slice::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline civector_slice& civector_slice::operator&=(const scivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline scivector& scivector::operator&=(const sivector_slice& v2) {
  return spsl_vv_intersectassign(*this,v2);
}

inline scivector& scivector::operator&=(const scivector_slice& v2) {
  return spsl_vv_intersectassign(*this,v2);
}

inline bool operator==(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const civector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const civector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const civector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const civector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const civector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const scvector_slice& v1, const civector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const cvector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const scvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const sivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const scivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const scivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const scivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const scivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const scvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const sivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const scivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const scivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const scivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const scivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const scvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const sivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector& v1, const scivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const scivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scvector& v1, const scivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const scivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const cvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const ivector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const civector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const civector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const civector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const civector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const cvector& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const ivector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const cvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scivector_slice& v1, const civector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const civector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const civector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const scvector_slice& v1, const civector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const civector_slice& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const cvector_slice& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const scivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_less<scivector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_less<scivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_less<srvector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_less<scvector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_less<sivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_less<scivector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_less<scivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_less<srvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_less<scvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_less<sivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_less<scivector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_less<scivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_less<srvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_less<scvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_less<sivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_less<scivector_slice,ivector,cinterval>(v1,v2);
}

inline bool operator<(const scivector_slice& v1, const civector& v2) {
  return slf_vv_less<scivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const civector& v2) {
  return slf_vv_less<srvector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const civector& v2) {
  return slf_vv_less<sivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<(const scvector_slice& v1, const civector& v2) {
  return slf_vv_less<scvector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_less<civector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_less<civector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_less<rvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_less<cvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_less<ivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_leq<scivector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_leq<scivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_leq<srvector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scvector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_leq<scvector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_leq<sivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_leq<scivector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_leq<scivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const scivector& v2) {
  return slsp_vv_leq<srvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scvector_slice& v1, const scivector& v2) {
  return slsp_vv_leq<scvector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_leq<sivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_leq<scivector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_leq<scivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const scivector_slice& v2) {
  return spsl_vv_leq<srvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scvector& v1, const scivector_slice& v2) {
  return spsl_vv_leq<scvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_leq<sivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_leq<scivector_slice,ivector,cinterval>(v1,v2);
}

inline bool operator<=(const scivector_slice& v1, const civector& v2) {
  return slf_vv_leq<scivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const civector& v2) {
  return slf_vv_leq<srvector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const civector& v2) {
  return slf_vv_leq<sivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<=(const scvector_slice& v1, const civector& v2) {
  return slf_vv_leq<scvector_slice,civector,cinterval>(v1,v2);
}

inline bool operator<=(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_leq<civector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_leq<civector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const rvector& v1, const scivector_slice& v2) {
  return fsl_vv_leq<rvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const cvector& v1, const scivector_slice& v2) {
  return fsl_vv_leq<cvector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator<=(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_leq<ivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_greater<scivector_slice,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_greater<scivector_slice,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_greater<scivector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_greater<scivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_greater<sivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_greater<scivector_slice,srvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_greater<scivector_slice,scvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_greater<scivector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_greater<scivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_greater<sivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_greater<scivector,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_greater<scivector,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_greater<scivector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_greater<scivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_greater<sivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_greater<scivector_slice,rvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_greater<scivector_slice,cvector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_greater<scivector_slice,ivector,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const civector& v2) {
  return slf_vv_greater<scivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const civector& v2) {
  return slf_vv_greater<sivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_greater<civector,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_greater<civector,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_greater<civector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_greater<civector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_greater<ivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_greater<scivector_slice,rvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_greater<scivector_slice,ivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_greater<scivector_slice,cvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_greater<scivector_slice,civector_slice,cinterval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_greater<sivector_slice,civector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_greater<civector_slice,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_greater<civector_slice,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_greater<civector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_greater<civector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_greater<ivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_geq<scivector_slice,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_geq<scivector_slice,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_geq<scivector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_geq<scivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const scivector_slice& v2) {
  return slsl_vv_geq<sivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const srvector& v2) {
  return slsp_vv_geq<scivector_slice,srvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const scvector& v2) {
  return slsp_vv_geq<scivector_slice,scvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const sivector& v2) {
  return slsp_vv_geq<scivector_slice,sivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const scivector& v2) {
  return slsp_vv_geq<scivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const scivector& v2) {
  return slsp_vv_geq<sivector_slice,scivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const srvector_slice& v2) {
  return spsl_vv_geq<scivector,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const scvector_slice& v2) {
  return spsl_vv_geq<scivector,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const sivector_slice& v2) {
  return spsl_vv_geq<scivector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector& v1, const scivector_slice& v2) {
  return spsl_vv_geq<scivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const scivector_slice& v2) {
  return spsl_vv_geq<sivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const rvector& v2) {
  return slf_vv_geq<scivector_slice,rvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const cvector& v2) {
  return slf_vv_geq<scivector_slice,cvector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const ivector& v2) {
  return slf_vv_geq<scivector_slice,ivector,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const civector& v2) {
  return slf_vv_geq<scivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const civector& v2) {
  return slf_vv_geq<sivector_slice,civector,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const srvector_slice& v2) {
  return fsl_vv_geq<civector,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const scvector_slice& v2) {
  return fsl_vv_geq<civector,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const sivector_slice& v2) {
  return fsl_vv_geq<civector,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector& v1, const scivector_slice& v2) {
  return fsl_vv_geq<civector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const scivector_slice& v2) {
  return fsl_vv_geq<ivector,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_geq<scivector_slice,rvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_geq<scivector_slice,ivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const cvector_slice& v2) {
  return slf_vv_geq<scivector_slice,cvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const scivector_slice& v1, const civector_slice& v2) {
  return slf_vv_geq<scivector_slice,civector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const civector_slice& v2) {
  return slf_vv_geq<sivector_slice,civector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_geq<civector_slice,srvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_geq<civector_slice,scvector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_geq<civector_slice,sivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const civector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_geq<civector_slice,scivector_slice,cinterval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const scivector_slice& v2) {
  return fsl_vv_geq<ivector_slice,scivector_slice,cinterval>(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const scivector_slice& v) {
  return sl_v_output<scivector_slice,cinterval>(os,v);
}

inline std::istream& operator>>(std::istream& is, scivector_slice& v) {
  return sl_v_input<scivector_slice,cinterval>(is,v);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const rvector& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,rvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const cvector& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,cvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const ivector& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,ivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const civector& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const civector& v2) {
  slf_vv_accu<cidotprecision,srvector_slice,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const civector& v2) {
  slf_vv_accu<cidotprecision,sivector_slice,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const civector& v2) {
  slf_vv_accu<cidotprecision,scvector_slice,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const ivector& v2) {
  slf_vv_accu<cidotprecision,scvector_slice,ivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const cvector& v2) {
  slf_vv_accu<cidotprecision,sivector_slice,cvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const srvector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector,srvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const scvector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const sivector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const rvector& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,rvector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,cvector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,ivector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const sivector_slice& v2) {
  fsl_vv_accu<cidotprecision,cvector,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const scvector_slice& v2) {
  fsl_vv_accu<cidotprecision,ivector,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const rvector_slice& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,rvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const ivector_slice& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,ivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const cvector_slice& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,cvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const civector_slice& v2) {
  slf_vv_accu<cidotprecision,scivector_slice,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const civector_slice& v2) {
  slf_vv_accu<cidotprecision,srvector_slice,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const civector_slice& v2) {
  slf_vv_accu<cidotprecision,scvector_slice,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const civector_slice& v2) {
  slf_vv_accu<cidotprecision,sivector_slice,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const ivector_slice& v2) {
  slf_vv_accu<cidotprecision,scvector_slice,ivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const cvector_slice& v2) {
  slf_vv_accu<cidotprecision,sivector_slice,cvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const srvector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector_slice,srvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const scvector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector_slice,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const sivector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector_slice,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,civector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,rvector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,ivector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const scivector_slice& v2) {
  fsl_vv_accu<cidotprecision,cvector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const sivector_slice& v2) {
  fsl_vv_accu<cidotprecision,cvector_slice,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const scvector_slice& v2) {
  fsl_vv_accu<cidotprecision,ivector_slice,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const srvector_slice& v2) {
  spsl_vv_accu<cidotprecision,scivector,srvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const scvector_slice& v2) {
  spsl_vv_accu<cidotprecision,scivector,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const sivector_slice& v2) {
  spsl_vv_accu<cidotprecision,scivector,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const scivector_slice& v2) {
  spsl_vv_accu<cidotprecision,scivector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const scivector_slice& v2) {
  spsl_vv_accu<cidotprecision,srvector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const scivector_slice& v2) {
  spsl_vv_accu<cidotprecision,scvector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const scivector_slice& v2) {
  spsl_vv_accu<cidotprecision,sivector,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const sivector_slice& v2) {
  spsl_vv_accu<cidotprecision,scvector,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const scvector_slice& v2) {
  spsl_vv_accu<cidotprecision,sivector,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const srvector& v2) {
  slsp_vv_accu<cidotprecision,scivector_slice,srvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const scvector& v2) {
  slsp_vv_accu<cidotprecision,scivector_slice,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const sivector& v2) {
  slsp_vv_accu<cidotprecision,scivector_slice,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const scivector& v2) {
  slsp_vv_accu<cidotprecision,scivector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const scivector& v2) {
  slsp_vv_accu<cidotprecision,srvector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const scivector& v2) {
  slsp_vv_accu<cidotprecision,sivector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const scivector& v2) {
  slsp_vv_accu<cidotprecision,scvector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const sivector& v2) {
  slsp_vv_accu<cidotprecision,scvector_slice,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const scvector& v2) {
  slsp_vv_accu<cidotprecision,sivector_slice,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const srvector_slice& v2) {
  slsl_vv_accu<cidotprecision,scivector_slice,srvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const scvector_slice& v2) {
  slsl_vv_accu<cidotprecision,scivector_slice,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const sivector_slice& v2) {
  slsl_vv_accu<cidotprecision,scivector_slice,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector_slice& v1, const scivector_slice& v2) {
  slsl_vv_accu<cidotprecision,scivector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const scivector_slice& v2) {
  slsl_vv_accu<cidotprecision,srvector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const scivector_slice& v2) {
  slsl_vv_accu<cidotprecision,scvector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const scivector_slice& v2) {
  slsl_vv_accu<cidotprecision,sivector_slice,scivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const scvector_slice& v2) {
  slsl_vv_accu<cidotprecision,sivector_slice,scvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector_slice& v1, const sivector_slice& v2) {
  slsl_vv_accu<cidotprecision,scvector_slice,sivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const cvector& v2) {
  spf_vv_accu<cidotprecision,scivector,cvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const rvector& v2) {
  spf_vv_accu<cidotprecision,scivector,rvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const ivector& v2) {
  spf_vv_accu<cidotprecision,scivector,ivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const civector& v2) {
  spf_vv_accu<cidotprecision,scivector,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const civector& v2) {
  spf_vv_accu<cidotprecision,scvector,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const civector& v2) {
  spf_vv_accu<cidotprecision,srvector,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const civector& v2) {
  spf_vv_accu<cidotprecision,sivector,civector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const ivector& v2) {
  spf_vv_accu<cidotprecision,scvector,ivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const cvector& v2) {
  spf_vv_accu<cidotprecision,sivector,cvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const rvector& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,rvector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,cvector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,ivector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,civector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const srvector& v2) {
  fsp_vv_accu<cidotprecision,civector,srvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const scvector& v2) {
  fsp_vv_accu<cidotprecision,civector,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector& v1, const sivector& v2) {
  fsp_vv_accu<cidotprecision,civector,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector& v1, const scvector& v2) {
  fsp_vv_accu<cidotprecision,ivector,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector& v1, const sivector& v2) {
  fsp_vv_accu<cidotprecision,cvector,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const cvector_slice& v2) {
  spf_vv_accu<cidotprecision,scivector,cvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const rvector_slice& v2) {
  spf_vv_accu<cidotprecision,scivector,rvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const ivector_slice& v2) {
  spf_vv_accu<cidotprecision,scivector,ivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const civector_slice& v2) {
  spf_vv_accu<cidotprecision,scivector,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const civector_slice& v2) {
  spf_vv_accu<cidotprecision,scvector,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const civector_slice& v2) {
  spf_vv_accu<cidotprecision,srvector,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const civector_slice& v2) {
  spf_vv_accu<cidotprecision,sivector,civector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const ivector_slice& v2) {
  spf_vv_accu<cidotprecision,scvector,ivector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const cvector_slice& v2) {
  spf_vv_accu<cidotprecision,sivector,cvector_slice,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,rvector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,cvector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,ivector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const scivector& v2) {
  fsp_vv_accu<cidotprecision,civector_slice,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const srvector& v2) {
  fsp_vv_accu<cidotprecision,civector_slice,srvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const scvector& v2) {
  fsp_vv_accu<cidotprecision,civector_slice,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const civector_slice& v1, const sivector& v2) {
  fsp_vv_accu<cidotprecision,civector_slice,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const scvector& v2) {
  fsp_vv_accu<cidotprecision,ivector_slice,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const cvector_slice& v1, const sivector& v2) {
  fsp_vv_accu<cidotprecision,cvector_slice,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const srvector& v2) {
  spsp_vv_accu<cidotprecision,scivector,srvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const scvector& v2) {
  spsp_vv_accu<cidotprecision,scivector,scvector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const sivector& v2) {
  spsp_vv_accu<cidotprecision,scivector,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scivector& v1, const scivector& v2) {
  spsp_vv_accu<cidotprecision,scivector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const srvector& v1, const scivector& v2) {
  spsp_vv_accu<cidotprecision,srvector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const scivector& v2) {
  spsp_vv_accu<cidotprecision,scvector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const scivector& v2) {
  spsp_vv_accu<cidotprecision,sivector,scivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const scvector& v1, const sivector& v2) {
  spsp_vv_accu<cidotprecision,scvector,sivector,sparse_cidot>(dot,v1,v2);
}

inline void accumulate(cidotprecision& dot, const sivector& v1, const scvector& v2) {
  spsp_vv_accu<cidotprecision,sivector,scvector,sparse_cidot>(dot,v1,v2);
}

} //namespace cxsc

#include "sparsevector.inl"

#endif
