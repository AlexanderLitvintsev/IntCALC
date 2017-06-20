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

/* CVS $Id: clss.hpp,v 1.8 2010/12/21 09:53:47 cxsc Exp $ */

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

#ifndef _CXSC_CLSS_HPP
#define _CXSC_CLSS_HPP

#ifndef _CXSC_LSS_CONSTANTS_DEFINED
#define _CXSC_LSS_CONSTANTS_DEFINED
//Control constants
static const int
  LSS_ONLY_PART_ONE = 0,
  LSS_ONLY_PART_TWO = 1,
  LSS_BOTH_PARTS    = 2;
#endif

#include <cimatrix.hpp>

namespace cxsc {

//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cmatrix&, cxsc::cimatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);

//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cvector&, cxsc::civector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);

//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);

//! Entry function of complex linear sytems solver
static inline void  clss(cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, int&, int=2, bool=false, int=LSS_BOTH_PARTS, int=-1);

} //namespace cxsc

#include <fastlss.hpp> 

#endif
 
