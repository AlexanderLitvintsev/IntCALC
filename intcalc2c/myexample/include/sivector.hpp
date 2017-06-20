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

/* CVS $Id: sivector.hpp,v 1.8 2010/12/21 09:54:02 cxsc Exp $ */

#ifndef _CXSC_SIVECTOR_HPP_INCLUDED
#define _CXSC_SIVECTOR_HPP_INCLUDED

#include <interval.hpp>
#include <ivector.hpp>
#include <vector>
#include <iostream>
#include <cidot.hpp>
#include <srvector.hpp>
#include <sparseidot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;
class sivector_slice;
class simatrix;
class simatrix_slice;
class simatrix_subv;
class scivector;
class scivector_slice;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;


class sivector {
  private:
    std::vector<int> p;
    std::vector<interval> x;
    int lb;
    int ub;
    int n; 

  public:
    sivector() : lb(0), ub(-1) , n(0) {
    }

    explicit sivector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    sivector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    sivector(const ivector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    sivector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(interval(v[i]));
          }
        }
    }

    sivector(const int n, const int nnz, const intvector& index, const ivector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[Lb(values)+i] != 0.0) {
          p.push_back(index[Lb(index)+i]);
          x.push_back(values[Lb(values)+i]);
        }
      }
      
    }

    sivector(const srvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(interval(v.x[i]));
    }

    sivector(const srvector_slice&);
    sivector(const sivector_slice&);
    sivector(const srmatrix_subv& A);
    sivector(const simatrix_subv& A);

    std::vector<int>& row_indices() {
      return p;
    }

    std::vector<interval>& values() {
      return x;
    }

    const std::vector<int>& row_indices() const {
      return p;
    }

    const std::vector<interval>& values() const {
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


    /* sivector& operator=(const sivector& v) {
      p = v.p;
      x = v.x;
      return *this;
    } */

    sivector& operator=(const srvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = interval(v.x[i]);
      return *this;
    } 

    sivector& operator=(const real& v) {
      return sp_vs_assign<sivector,real,interval>(*this,v);
    }

    sivector& operator=(const interval& v) {
      return sp_vs_assign<sivector,interval,interval>(*this,v);
    }

    sivector& operator=(const rvector& v) {
      return spf_vv_assign<sivector,rvector,interval>(*this,v);
    }

    sivector& operator=(const ivector& v) {
      return spf_vv_assign<sivector,ivector,interval>(*this,v);
    }

    sivector& operator=(const rvector_slice& v) {
      return spf_vv_assign<sivector,rvector_slice,interval>(*this,v);
    }

    sivector& operator=(const sivector_slice&);
    sivector& operator=(const srvector_slice&);

    interval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, interval(0.0));

      return x[k];
    }

    interval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    interval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator()(const int)"));
#endif
      interval r(0.0);

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    sivector operator()(const intvector& per) {
      sivector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,interval> work;
      for(int i=0 ; i<get_nnz() ; i++) {
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
      }
 
      for(std::map<int,interval>::iterator it=work.begin() ; it!=work.end() ; it++) {
         v.p.push_back(it->first);
         v.x.push_back(it->second);
      }

      return v;
    }

    sivector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    sivector_slice operator()(const int, const int);

    sivector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    sivector& operator*=(const interval& s) {
      return sp_vs_multassign(*this,s);
    }

    sivector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    sivector& operator/=(const interval& s) {
      return sp_vs_divassign(*this,s);
    }

    sivector& operator+=(const rvector& v) {
      return spf_vv_addassign(*this,v);
    }

    sivector& operator+=(const ivector& v) {
      return spf_vv_addassign(*this,v);
    }

    sivector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    sivector& operator+=(const ivector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    sivector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    sivector& operator+=(const sivector& v) {
      return spsp_vv_addassign(*this,v);
    }

    sivector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    sivector& operator-=(const ivector& v) {
      return spf_vv_subassign(*this,v);
    }

    sivector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    sivector& operator-=(const ivector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    sivector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    sivector& operator-=(const sivector& v) {
      return spsp_vv_subassign(*this,v);
    }

    sivector& operator|=(const rvector& v) {
      return spf_vv_hullassign(*this,v);
    }

    sivector& operator|=(const ivector& v) {
      return spf_vv_hullassign(*this,v);
    }

    sivector& operator|=(const rvector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    sivector& operator|=(const ivector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    sivector& operator|=(const srvector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    sivector& operator|=(const sivector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    sivector& operator&=(const ivector_slice& v) {
      return spf_vv_intersectassign(*this,v);
    }

    sivector& operator&=(const sivector& v) {
      return spsp_vv_intersectassign(*this,v);
    }

    sivector& operator+=(const srvector_slice&);
    sivector& operator+=(const sivector_slice&);
    sivector& operator-=(const srvector_slice&);
    sivector& operator-=(const sivector_slice&);

    friend void SetLb(sivector&, const int);
    friend void SetUb(sivector&, const int);
    friend int Lb(const sivector&);
    friend int Ub(const sivector&);
    friend srvector Inf(const sivector&);
    friend srvector Sup(const sivector&);
    friend sivector Re(const scivector&);
    friend sivector Im(const scivector&);
    friend sivector abs(const sivector&);
    friend sivector abs(const sivector_slice&);
    friend srvector mid(const sivector&);
    friend srvector diam(const sivector&);
    friend sivector abs(const scivector&);
    friend sivector abs(const scivector_slice&);
    friend srvector absmin(const sivector&);
    friend srvector absmax(const sivector&);
    friend int VecLen(const sivector&);
    friend sivector Blow(const sivector&, const real&);

    friend class srvector_slice;
    friend class sivector_slice;
    friend class scivector_slice;
    friend class scivector;
    friend class ivector;
    friend class ivector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline ivector::ivector(const sivector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector::ivector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector& ivector::operator=(const sivector& v) {
  return fsp_vv_assign<ivector,sivector,interval>(*this,v);
}

inline ivector& ivector::operator=(const sivector_slice& v) {
  return fsl_vv_assign<ivector,sivector_slice,interval>(*this,v);
}

inline ivector& ivector::operator=(const srvector& v) {
  return fsp_vv_assign<ivector,srvector,interval>(*this,v);
}

inline ivector& ivector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<ivector,srvector_slice,interval>(*this,v);
}

inline void SetLb(sivector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

inline void SetUb(sivector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

inline int Lb(const sivector& v) {
  return v.lb;
}

inline int Ub(const sivector& v) {
  return v.ub;
}

inline void Resize(sivector& v) {
  sp_v_resize(v);
}

inline void Resize(sivector& v, const int n) {
  sp_v_resize(v,n);
}

inline void Resize(sivector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

inline srvector Inf(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Inf(v.x[i]);
  return res;
}

inline srvector Sup(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Sup(v.x[i]);
  return res;
}

inline sivector abs(const sivector& v) {
  sivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

inline srvector absmin(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(AbsMin(v.x[i]));
  res.dropzeros();
  return res;
}
inline srvector absmax(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(AbsMax(v.x[i]));
  res.dropzeros();
  return res;
}

inline srvector mid(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++) {
    res.x.push_back(mid(v.x[i]));
  }
  return res;
}

inline srvector diam(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

inline int VecLen(const sivector& v) {
  return v.n;
}

inline sivector Blow(const sivector& v, const real& eps) {
  sivector res(v);
  for(unsigned int i=0 ; i<v.x.size() ; i++)
    res.x[i] = Blow(v.x[i],eps);
  return res;
}

inline bool in (const sivector& v1, const sivector& v2) {
  for(int i=0 ; i<VecLen(v1) ; i++)
    if(!in(v1(i+Lb(v1)), v2(i+Lb(v2)))) return false;
  return true;
}

inline bool Zero(const sivector& v1) {
  for(int i=0 ; i<VecLen(v1) ; i++)
    if(v1(i+Lb(v1)) != 0.0) return false;
  return true;
}

inline sivector operator-(const sivector& v) {
  return sp_v_negative(v);
}

inline interval operator*(const sivector& v1, const ivector& v2) {
  return spf_vv_mult<sivector,ivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const rvector& v2) {
  return spf_vv_mult<sivector,rvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector& v1, const ivector& v2) {
  return spf_vv_mult<srvector,ivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const rvector& v1, const sivector& v2) {
  return fsp_vv_mult<rvector,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector& v1, const srvector& v2) {
  return fsp_vv_mult<ivector,srvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector& v1, const sivector& v2) {
  return fsp_vv_mult<ivector,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_mult<sivector,rvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_mult<sivector,ivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_mult<srvector,ivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<ivector_slice,srvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<ivector_slice,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<rvector_slice,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const srvector& v2) {
  return spsp_vv_mult<sivector,srvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector& v1, const sivector& v2) {
  return spsp_vv_mult<srvector,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const sivector& v2) {
  return spsp_vv_mult<sivector,sivector,interval,sparse_idot>(v1,v2);
}

inline sivector operator*(const sivector& v, const real& s) {
  return sp_vs_mult<sivector,real,sivector>(v,s);
}

inline sivector operator*(const sivector& v, const interval& s) {
  return sp_vs_mult<sivector,interval,sivector>(v,s);
}

inline sivector operator*(const srvector& v, const interval& s) {
  return sp_vs_mult<srvector,interval,sivector>(v,s);
}

inline sivector operator/(const sivector& v, const real& s) {
  return sp_vs_div<sivector,real,sivector>(v,s);
}

inline sivector operator/(const sivector& v, const interval& s) {
  return sp_vs_div<sivector,interval,sivector>(v,s);
}

inline sivector operator/(const srvector& v, const interval& s) {
  return sp_vs_div<srvector,interval,sivector>(v,s);
}

inline sivector operator*(const real& s, const sivector& v) {
  return sp_sv_mult<real,sivector,sivector>(s,v);
}

inline sivector operator*(const interval& s, const sivector& v) {
  return sp_sv_mult<interval,sivector,sivector>(s,v);
}

inline sivector operator*(const interval& s, const srvector& v) {
  return sp_sv_mult<interval,srvector,sivector>(s,v);
}

inline ivector operator+(const ivector& v1, const srvector& v2) {
  return fsp_vv_add<ivector,srvector,ivector>(v1,v2);
}

inline ivector operator+(const rvector& v1, const sivector& v2) {
  return fsp_vv_add<rvector,sivector,ivector>(v1,v2);
}

inline ivector operator+(const ivector& v1, const sivector& v2) {
  return fsp_vv_add<ivector,sivector,ivector>(v1,v2);
}

inline ivector operator+(const sivector& v1, const rvector& v2) {
  return spf_vv_add<sivector,rvector,ivector>(v1,v2);
}

inline ivector operator+(const srvector& v1, const ivector& v2) {
  return spf_vv_add<srvector,ivector,ivector>(v1,v2);
}

inline ivector operator+(const sivector& v1, const ivector& v2) {
  return spf_vv_add<sivector,ivector,ivector>(v1,v2);
}

inline ivector operator+(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_add<ivector_slice,srvector,ivector>(v1,v2);
}

inline ivector operator+(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_add<rvector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator+(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_add<ivector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator+(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_add<sivector,rvector_slice,ivector>(v1,v2);
}

inline ivector operator+(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_add<srvector,ivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_add<sivector,ivector_slice,ivector>(v1,v2);
}

inline sivector operator+(const sivector& v1, const srvector& v2) {
  return spsp_vv_add<sivector,srvector,sivector,interval>(v1,v2);
}

inline sivector operator+(const srvector& v1, const sivector& v2) {
  return spsp_vv_add<srvector,sivector,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector& v1, const sivector& v2) {
  return spsp_vv_add<sivector,sivector,sivector,interval>(v1,v2);
}

inline ivector operator-(const ivector& v1, const srvector& v2) {
  return fsp_vv_sub<ivector,srvector,ivector>(v1,v2);
}

inline ivector operator-(const rvector& v1, const sivector& v2) {
  return fsp_vv_sub<rvector,sivector,ivector>(v1,v2);
}

inline ivector operator-(const ivector& v1, const sivector& v2) {
  return fsp_vv_sub<ivector,sivector,ivector>(v1,v2);
}

inline ivector operator-(const sivector& v1, const rvector& v2) {
  return spf_vv_sub<sivector,rvector,ivector>(v1,v2);
}

inline ivector operator-(const srvector& v1, const ivector& v2) {
  return spf_vv_sub<srvector,ivector,ivector>(v1,v2);
}

inline ivector operator-(const sivector& v1, const ivector& v2) {
  return spf_vv_sub<sivector,ivector,ivector>(v1,v2);
}

inline ivector operator-(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<ivector_slice,srvector,ivector>(v1,v2);
}

inline ivector operator-(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<rvector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator-(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<ivector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator-(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_sub<sivector,rvector_slice,ivector>(v1,v2);
}

inline ivector operator-(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_sub<srvector,ivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_sub<sivector,ivector_slice,ivector>(v1,v2);
}

inline sivector operator-(const sivector& v1, const srvector& v2) {
  return spsp_vv_sub<sivector,srvector,sivector,interval>(v1,v2);
}

inline sivector operator-(const srvector& v1, const sivector& v2) {
  return spsp_vv_sub<srvector,sivector,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector& v1, const sivector& v2) {
  return spsp_vv_sub<sivector,sivector,sivector,interval>(v1,v2);
}

inline ivector operator|(const rvector& v1, const srvector& v2) {
  return fsp_vv_hull<rvector,srvector,ivector>(v1,v2);
}

inline ivector operator|(const srvector& v1, const rvector& v2) {
  return spf_vv_hull<srvector,rvector,ivector>(v1,v2);
}

inline ivector operator|(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<rvector_slice,srvector,ivector>(v1,v2);
}

inline ivector operator|(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_hull<srvector,rvector_slice,ivector>(v1,v2);
}

inline sivector operator|(const srvector& v1, const srvector& v2) {
  return spsp_vv_hull<srvector,srvector,sivector,interval>(v1,v2);
}

inline ivector operator|(const ivector& v1, const srvector& v2) {
  return fsp_vv_hull<ivector,srvector,ivector>(v1,v2);
}

inline ivector operator|(const rvector& v1, const sivector& v2) {
  return fsp_vv_hull<rvector,sivector,ivector>(v1,v2);
}

inline ivector operator|(const ivector& v1, const sivector& v2) {
  return fsp_vv_hull<ivector,sivector,ivector>(v1,v2);
}

inline ivector operator|(const sivector& v1, const rvector& v2) {
  return spf_vv_hull<sivector,rvector,ivector>(v1,v2);
}

inline ivector operator|(const srvector& v1, const ivector& v2) {
  return spf_vv_hull<srvector,ivector,ivector>(v1,v2);
}

inline ivector operator|(const sivector& v1, const ivector& v2) {
  return spf_vv_hull<sivector,ivector,ivector>(v1,v2);
}

inline ivector operator|(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<ivector_slice,srvector,ivector>(v1,v2);
}

inline ivector operator|(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<rvector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator|(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<ivector_slice,sivector,ivector>(v1,v2);
}

inline ivector operator|(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_hull<sivector,rvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_hull<srvector,ivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_hull<sivector,ivector_slice,ivector>(v1,v2);
}

inline sivector operator|(const sivector& v1, const srvector& v2) {
  return spsp_vv_hull<sivector,srvector,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector& v1, const sivector& v2) {
  return spsp_vv_hull<srvector,sivector,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector& v1, const sivector& v2) {
  return spsp_vv_hull<sivector,sivector,sivector,interval>(v1,v2);
}

inline sivector operator&(const ivector& v1, const sivector& v2) {
  return fsp_vv_intersect<ivector,sivector,ivector>(v1,v2);
}

inline sivector operator&(const sivector& v1, const ivector& v2) {
  return spf_vv_intersect<sivector,ivector,ivector>(v1,v2);
}

inline sivector operator&(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_intersect<ivector_slice,sivector,ivector>(v1,v2);
}

inline sivector operator&(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_intersect<sivector,ivector_slice,ivector>(v1,v2);
}

inline sivector operator&(const sivector& v1, const sivector& v2) {
  return spsp_vv_intersect<sivector,sivector,sivector,interval>(v1,v2);
}

inline ivector& ivector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector& ivector::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline ivector& ivector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector& ivector::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector& ivector::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline bool operator==(const sivector& v1, const sivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const sivector& v2) {
  return spsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const ivector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const ivector& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const sivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const sivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const ivector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const ivector& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const ivector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const ivector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

inline bool operator<(const sivector& v1, const sivector& v2) {
  return spsp_vv_less<sivector,sivector,interval>(v1,v2);
}

inline bool operator<(const srvector& v1, const sivector& v2) {
  return spsp_vv_less<srvector,sivector,interval>(v1,v2);
}

inline bool operator<(const srvector& v1, const ivector& v2) {
  return spf_vv_less<srvector,ivector,interval>(v1,v2);
}

inline bool operator<(const sivector& v1, const ivector& v2) {
  return spf_vv_less<sivector,ivector,interval>(v1,v2);
}

inline bool operator<(const rvector& v1, const sivector& v2) {
  return fsp_vv_less<rvector,sivector,interval>(v1,v2);
}

inline bool operator<(const ivector& v1, const sivector& v2) {
  return fsp_vv_less<ivector,sivector,interval>(v1,v2);
}

inline bool operator<(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_less<srvector,ivector_slice,interval>(v1,v2);
}

inline bool operator<(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_less<sivector,ivector_slice,interval>(v1,v2);
}

inline bool operator<(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_less<rvector_slice,sivector,interval>(v1,v2);
}

inline bool operator<(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_less<ivector_slice,sivector,interval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const sivector& v2) {
  return spsp_vv_leq<sivector,sivector,interval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const sivector& v2) {
  return spsp_vv_leq<srvector,sivector,interval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const ivector& v2) {
  return spf_vv_leq<srvector,ivector,interval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const ivector& v2) {
  return spf_vv_leq<sivector,ivector,interval>(v1,v2);
}

inline bool operator<=(const rvector& v1, const sivector& v2) {
  return fsp_vv_leq<rvector,sivector,interval>(v1,v2);
}

inline bool operator<=(const ivector& v1, const sivector& v2) {
  return fsp_vv_leq<ivector,sivector,interval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_leq<srvector,ivector_slice,interval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_leq<sivector,ivector_slice,interval>(v1,v2);
}

inline bool operator<=(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_leq<rvector_slice,sivector,interval>(v1,v2);
}

inline bool operator<=(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_leq<ivector_slice,sivector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const sivector& v2) {
  return spsp_vv_greater<sivector,sivector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const srvector& v2) {
  return spsp_vv_greater<sivector,srvector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const rvector& v2) {
  return spf_vv_greater<sivector,rvector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const ivector& v2) {
  return spf_vv_greater<sivector,ivector,interval>(v1,v2);
}

inline bool operator>(const ivector& v1, const srvector& v2) {
  return fsp_vv_greater<ivector,srvector,interval>(v1,v2);
}

inline bool operator>(const ivector& v1, const sivector& v2) {
  return fsp_vv_greater<ivector,sivector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_greater<sivector,rvector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_greater<sivector,ivector_slice,interval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_greater<ivector_slice,srvector,interval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_greater<ivector_slice,sivector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const sivector& v2) {
  return spsp_vv_geq<sivector,sivector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const srvector& v2) {
  return spsp_vv_geq<sivector,srvector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const rvector& v2) {
  return spf_vv_geq<sivector,rvector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const ivector& v2) {
  return spf_vv_geq<sivector,ivector,interval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const srvector& v2) {
  return fsp_vv_geq<ivector,srvector,interval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const sivector& v2) {
  return fsp_vv_geq<ivector,sivector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_geq<sivector,rvector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_geq<sivector,ivector_slice,interval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_geq<ivector_slice,srvector,interval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_geq<ivector_slice,sivector,interval>(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const sivector& v) {
  return sp_v_output<sivector,interval>(os,v);
}

inline std::istream& operator>>(std::istream& is, sivector& v) {
  return sp_v_input<sivector,interval>(is,v);
}


class sivector_slice {
  private:
    std::vector<int>& p;
    std::vector<interval>& x;
    sivector& orig;
    int start,end;
    int lb;
    int ub;
    int n;
    int nnz;
    int offset;

    sivector_slice(sivector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1) {
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

    interval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, interval(0.0));
      end++;

      return x[k];
    }

    interval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    interval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator()(const int)"));
#endif
      interval r(0.0);

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    sivector_slice& operator=(const real& v) {
      return sl_vs_assign<sivector_slice,real,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const interval& v) {
      return sl_vs_assign<sivector_slice,interval,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<sivector_slice,srvector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const sivector_slice& v) {
      return slsl_vv_assign<sivector_slice,sivector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<sivector_slice,srvector,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const sivector& v) {
      return slsp_vv_assign<sivector_slice,sivector,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const rvector& v) {
      return slf_vv_assign<sivector_slice,rvector,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const ivector& v) {
      return slf_vv_assign<sivector_slice,ivector,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<sivector_slice,rvector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator=(const ivector_slice& v) {
      return slf_vv_assign<sivector_slice,ivector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    sivector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    sivector_slice& operator*=(const interval& s) {
      return sl_vs_multassign(*this,s);
    }

    sivector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    sivector_slice& operator/=(const interval& s) {
      return sl_vs_divassign(*this,s);
    }

    sivector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<sivector_slice,rvector,interval>(*this,v);
    }

    sivector_slice& operator+=(const ivector& v) {
      return slf_vv_addassign<sivector_slice,ivector,interval>(*this,v);
    }

    sivector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    sivector_slice& operator+=(const ivector_slice& v) {
      return slf_vv_addassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    sivector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    sivector_slice& operator+=(const sivector& v) {
      return slsp_vv_addassign(*this,v);
    }

    sivector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    sivector_slice& operator+=(const sivector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    sivector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<sivector_slice,rvector,interval>(*this,v);
    }

    sivector_slice& operator-=(const ivector& v) {
      return slf_vv_subassign<sivector_slice,ivector,interval>(*this,v);
    }

    sivector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    sivector_slice& operator-=(const ivector_slice& v) {
      return slf_vv_subassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    sivector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    sivector_slice& operator-=(const sivector& v) {
      return slsp_vv_subassign(*this,v);
    }

    sivector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    sivector_slice& operator-=(const sivector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    sivector_slice& operator|=(const rvector& v) {
      return slf_vv_hullassign<sivector_slice,rvector,interval>(*this,v);
    }

    sivector_slice& operator|=(const ivector& v) {
      return slf_vv_hullassign<sivector_slice,ivector,interval>(*this,v);
    }

    sivector_slice& operator|=(const rvector_slice& v) {
      return slf_vv_hullassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    sivector_slice& operator|=(const ivector_slice& v) {
      return slf_vv_hullassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    sivector_slice& operator|=(const srvector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    sivector_slice& operator|=(const sivector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    sivector_slice& operator|=(const srvector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    sivector_slice& operator|=(const sivector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    sivector_slice& operator&=(const ivector& v) {
      return slf_vv_intersectassign<sivector_slice,ivector,interval>(*this,v);
    }

    sivector_slice& operator&=(const ivector_slice& v) {
      return slf_vv_intersectassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    sivector_slice& operator&=(const sivector& v) {
      return slsp_vv_intersectassign(*this,v);
    }

    sivector_slice& operator&=(const sivector_slice& v) {
      return slsl_vv_intersectassign(*this,v);
    }

    friend int Lb(const sivector_slice&);
    friend int Ub(const sivector_slice&);
    friend srvector Inf(const sivector_slice&);
    friend srvector Sup(const sivector_slice&);
    friend sivector abs(const sivector_slice&);
    friend srvector mid(const sivector_slice&);
    friend srvector diam(const sivector_slice&);
    friend int VecLen(const sivector_slice&);

//     friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
//     friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class sivector;
    friend class scivector;
    friend class ivector;
    friend class ivector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline ivector::ivector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector::ivector(const sivector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector_slice& ivector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline sivector::sivector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(interval(s.x[i]));
  }

}

inline sivector::sivector(const sivector_slice& s) : lb(s.lb), ub(s.ub), n(s.n) {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline sivector& sivector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<sivector,srvector_slice,interval>(*this,v);
}

inline sivector& sivector::operator=(const sivector_slice& v) {
  return spsl_vv_assign<sivector,sivector_slice,interval>(*this,v);
}

inline sivector_slice sivector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
  if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator()(const int,const int)"));
#endif
  return sivector_slice(*this,i,j);
}

inline sivector operator-(const sivector_slice& v) {
  return sl_v_negative<sivector_slice,sivector>(v);
}

inline int Lb(const sivector_slice& v) {
  return v.lb;
}

inline int Ub(const sivector_slice& v) {
  return v.ub;
}

inline srvector Inf(const sivector_slice& v) {
  return Inf(sivector(v));
}

inline srvector Sup(const sivector_slice& v) {
  return Sup(sivector(v));
}

inline sivector abs(const sivector_slice& v) {
  sivector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

inline srvector mid(const sivector_slice& v) {
  srvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(mid(v.x[i]));
  return res;
}

inline srvector diam(const sivector_slice& v) {
  srvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<v.end ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

inline int VecLen(const sivector_slice& v) {
  return v.n;
}

inline interval operator*(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_mult<sivector_slice,rvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_mult<srvector_slice,ivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_mult<sivector_slice,ivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<ivector,srvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<rvector,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<ivector,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<sivector_slice,rvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<srvector_slice,ivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<sivector_slice,ivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<ivector_slice,srvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<rvector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<ivector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<sivector,srvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<srvector,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<sivector,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<sivector_slice,srvector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<srvector_slice,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<sivector_slice,sivector,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<sivector_slice,srvector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<srvector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline interval operator*(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<sivector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

inline sivector operator*(const sivector_slice& v, const real& s) {
  return sp_vs_mult<sivector_slice,real,sivector>(v,s);
}

inline sivector operator*(const sivector_slice& v, const interval& s) {
  return sp_vs_mult<sivector_slice,interval,sivector>(v,s);
}

inline sivector operator*(const srvector_slice& v, const interval& s) {
  return sp_vs_mult<srvector_slice,interval,sivector>(v,s);
}

inline sivector operator/(const sivector_slice& v, const real& s) {
  return sp_vs_div<sivector_slice,real,sivector>(v,s);
}

inline sivector operator/(const sivector_slice& v, const interval& s) {
  return sp_vs_div<sivector_slice,interval,sivector>(v,s);
}

inline sivector operator/(const srvector_slice& v, const interval& s) {
  return sp_vs_div<srvector_slice,interval,sivector>(v,s);
}

inline sivector operator*(const real& s, const sivector_slice& v) {
  return sp_sv_mult<real,sivector_slice,sivector>(s,v);
}

inline sivector operator*(const interval& s, const sivector_slice& v) {
  return sp_sv_mult<interval,sivector_slice,sivector>(s,v);
}

inline sivector operator*(const interval& s, const srvector_slice& v) {
  return sp_sv_mult<interval,srvector_slice,sivector>(s,v);
}

inline ivector operator+(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_add<ivector,srvector_slice,ivector>(v1,v2);
}

inline ivector operator+(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_add<rvector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_add<ivector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_add<sivector_slice,rvector,ivector>(v1,v2);
}

inline ivector operator+(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_add<srvector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator+(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_add<sivector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator+(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<ivector_slice,srvector_slice,ivector>(v1,v2);
}

inline ivector operator+(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<rvector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<ivector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<sivector_slice,rvector_slice,ivector>(v1,v2);
}

inline ivector operator+(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<srvector_slice,ivector_slice,ivector>(v1,v2);
}

inline ivector operator+(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<sivector_slice,ivector_slice,ivector>(v1,v2);
}

inline sivector operator+(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_add<sivector,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_add<srvector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_add<sivector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_add<sivector_slice,srvector,sivector,interval>(v1,v2);
}

inline sivector operator+(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_add<srvector_slice,sivector,sivector,interval>(v1,v2);
}

inline sivector operator+(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_add<sivector_slice,sivector,sivector,interval>(v1,v2);
}

inline ivector operator-(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<ivector,srvector_slice,ivector>(v1,v2);
}

inline ivector operator-(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<rvector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<ivector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_sub<sivector_slice,rvector,ivector>(v1,v2);
}

inline ivector operator-(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_sub<srvector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator-(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_sub<sivector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator-(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<ivector_slice,srvector_slice,ivector>(v1,v2);
}

inline ivector operator-(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<rvector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<ivector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<sivector_slice,rvector_slice,ivector>(v1,v2);
}

inline ivector operator-(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<srvector_slice,ivector_slice,ivector>(v1,v2);
}

inline ivector operator-(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<sivector_slice,ivector_slice,ivector>(v1,v2);
}

inline sivector operator-(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<sivector,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<srvector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<sivector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<sivector_slice,srvector,sivector,interval>(v1,v2);
}

inline sivector operator-(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<srvector_slice,sivector,sivector,interval>(v1,v2);
}

inline sivector operator-(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<sivector_slice,sivector,sivector,interval>(v1,v2);
}

inline ivector operator|(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<rvector,srvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_hull<srvector_slice,rvector,ivector>(v1,v2);
}

inline ivector operator|(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<rvector_slice,srvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<srvector_slice,rvector_slice,ivector>(v1,v2);
}

inline sivector operator|(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<srvector_slice,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<srvector,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<srvector_slice,srvector,sivector,interval>(v1,v2);
}

inline ivector operator|(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<ivector,srvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<rvector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<ivector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_hull<sivector_slice,rvector,ivector>(v1,v2);
}

inline ivector operator|(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_hull<srvector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator|(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_hull<sivector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator|(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<ivector_slice,srvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<rvector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<ivector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<sivector_slice,rvector_slice,ivector>(v1,v2);
}

inline ivector operator|(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<srvector_slice,ivector_slice,ivector>(v1,v2);
}

inline ivector operator|(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<sivector_slice,ivector_slice,ivector>(v1,v2);
}

inline sivector operator|(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<sivector,srvector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<srvector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<sivector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<sivector_slice,srvector,sivector,interval>(v1,v2);
}

inline sivector operator|(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<srvector_slice,sivector,sivector,interval>(v1,v2);
}

inline sivector operator|(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<sivector_slice,sivector,sivector,interval>(v1,v2);
}

inline ivector operator&(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<ivector,sivector_slice,ivector>(v1,v2);
}

inline ivector operator&(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_intersect<sivector_slice,ivector,ivector>(v1,v2);
}

inline ivector operator&(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<ivector_slice,sivector_slice,ivector>(v1,v2);
}

inline ivector operator&(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_intersect<sivector_slice,ivector_slice,ivector>(v1,v2);
}

inline sivector operator&(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_intersect<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator&(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_intersect<sivector,sivector_slice,sivector,interval>(v1,v2);
}

inline sivector operator&(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_intersect<sivector_slice,sivector,sivector,interval>(v1,v2);
}

inline ivector& ivector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector& ivector::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline sivector& sivector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline sivector& sivector::operator+=(const sivector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline ivector& ivector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector& ivector::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline sivector& sivector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline sivector& sivector::operator-=(const sivector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline ivector& ivector::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline bool operator==(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator==(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const sivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const sivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const ivector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const ivector& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const ivector& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const sivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const sivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const srvector& v1, const sivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector& v1, const sivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const srvector_slice& v1, const ivector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const sivector_slice& v1, const ivector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const rvector_slice& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator!=(const ivector_slice& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_less<srvector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_less<sivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_less<srvector_slice,sivector,interval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_less<sivector_slice,sivector,interval>(v1,v2);
}

inline bool operator<(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_less<srvector,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_less<sivector,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_less<srvector_slice,ivector,interval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_less<sivector_slice,ivector,interval>(v1,v2);
}

inline bool operator<(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_less<rvector,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_less<ivector,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_less<srvector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator<(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_less<sivector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator<(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_less<rvector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator<(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_less<ivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_leq<sivector_slice,sivector_slice,interval>(v1,v2);
}
inline bool operator<=(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_leq<srvector_slice,sivector,interval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_leq<sivector_slice,sivector,interval>(v1,v2);
}

inline bool operator<=(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_leq<srvector,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_leq<sivector,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_leq<srvector_slice,ivector,interval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_leq<sivector_slice,ivector,interval>(v1,v2);
}

inline bool operator<=(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_leq<rvector,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_leq<ivector,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_leq<srvector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator<=(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_leq<sivector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator<=(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_leq<rvector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator<=(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_leq<ivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_greater<sivector_slice,srvector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_greater<sivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_greater<sivector_slice,srvector,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_greater<sivector_slice,sivector,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_greater<sivector,srvector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_greater<sivector,sivector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_greater<sivector_slice,rvector,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_greater<sivector_slice,ivector,interval>(v1,v2);
}

inline bool operator>(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_greater<ivector,srvector_slice,interval>(v1,v2);
}

inline bool operator>(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_greater<ivector,sivector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_greater<sivector_slice,rvector_slice,interval>(v1,v2);
}

inline bool operator>(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_greater<sivector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_greater<ivector_slice,srvector_slice,interval>(v1,v2);
}

inline bool operator>(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_greater<ivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_geq<sivector_slice,srvector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_geq<sivector_slice,sivector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_geq<sivector_slice,srvector,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_geq<sivector_slice,sivector,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_geq<sivector,srvector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_geq<sivector,sivector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_geq<sivector_slice,rvector,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_geq<sivector_slice,ivector,interval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_geq<ivector,srvector_slice,interval>(v1,v2);
}

inline bool operator>=(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_geq<ivector,sivector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_geq<sivector_slice,rvector_slice,interval>(v1,v2);
}

inline bool operator>=(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_geq<sivector_slice,ivector_slice,interval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_geq<ivector_slice,srvector_slice,interval>(v1,v2);
}

inline bool operator>=(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_geq<ivector_slice,sivector_slice,interval>(v1,v2);
}

inline std::ostream& operator<<(std::ostream& os, const sivector_slice& v) {
  return sl_v_output<sivector_slice,interval>(os,v);
}

inline std::istream& operator>>(std::istream& is, sivector_slice& v) {
  return sl_v_input<sivector_slice,interval>(is,v);
}

inline void accumulate(idotprecision& dot, const sivector& x, const sivector& y) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const srvector& y) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector& x, const sivector& y) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const ivector& y) {
  spf_vv_accu<idotprecision,sivector,ivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const rvector& y) {
  spf_vv_accu<idotprecision,sivector,rvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector& x, const ivector& y) {
  spf_vv_accu<idotprecision,srvector,ivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const ivector_slice& y) {
  spf_vv_accu<idotprecision,sivector,ivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const rvector_slice& y) {
  spf_vv_accu<idotprecision,sivector,rvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector& x, const ivector_slice& y) {
  spf_vv_accu<idotprecision,srvector,ivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector& x, const sivector& y) {
  fsp_vv_accu<idotprecision,ivector,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector& x, const srvector& y) {
  fsp_vv_accu<idotprecision,ivector,srvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const rvector& x, const sivector& y) {
  fsp_vv_accu<idotprecision,rvector,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector_slice& x, const sivector& y) {
  fsp_vv_accu<idotprecision,ivector_slice,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector_slice& x, const srvector& y) {
  fsp_vv_accu<idotprecision,ivector_slice,srvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const rvector_slice& x, const sivector& y) {
  fsp_vv_accu<idotprecision,rvector_slice,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const ivector& y) {
  slf_vv_accu<idotprecision,sivector_slice,ivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const rvector& y) {
  slf_vv_accu<idotprecision,sivector_slice,rvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const ivector& y) {
  slf_vv_accu<idotprecision,srvector_slice,ivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const ivector_slice& y) {
  slf_vv_accu<idotprecision,sivector_slice,ivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const rvector_slice& y) {
  slf_vv_accu<idotprecision,sivector_slice,rvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const ivector_slice& y) {
  slf_vv_accu<idotprecision,srvector_slice,ivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,ivector,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector& x, const srvector_slice& y) {
  fsl_vv_accu<idotprecision,ivector,srvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const rvector& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,rvector,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector_slice& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,ivector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const ivector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<idotprecision,ivector_slice,srvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const rvector_slice& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,rvector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const sivector_slice& y) {
  slsl_vv_accu<idotprecision,sivector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<idotprecision,sivector_slice,srvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const sivector_slice& y) {
  slsl_vv_accu<idotprecision,srvector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const sivector_slice& y) {
  spsl_vv_accu<idotprecision,sivector,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector& x, const srvector_slice& y) {
  spsl_vv_accu<idotprecision,sivector,srvector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector& x, const sivector_slice& y) {
  spsl_vv_accu<idotprecision,srvector,sivector_slice,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const sivector& y) {
  slsp_vv_accu<idotprecision,sivector_slice,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const sivector_slice& x, const srvector& y) {
  slsp_vv_accu<idotprecision,sivector_slice,srvector,sparse_idot>(dot,x,y);
}

inline void accumulate(idotprecision& dot, const srvector_slice& x, const sivector& y) {
  slsp_vv_accu<idotprecision,srvector_slice,sivector,sparse_idot>(dot,x,y);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const rvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const rvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const rvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const rvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const ivector_slice& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const rvector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const sivector_slice& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

inline void accumulate(cidotprecision& dot, const srvector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

} //namespace cxsc

#include "sparsevector.inl"

#endif 
