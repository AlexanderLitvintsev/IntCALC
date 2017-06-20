// Borland C++ Builder
// Copyright (c) 1995, 1999 by Borland International
// All rights reserved

// (DO NOT EDIT: machine generated header) 'drag.pas' rev: 5.00

#ifndef dragHPP
#define dragHPP

#pragma delphiheader begin
#pragma option push -w-
#pragma option push -Vx
#include <SysUtils.hpp>	// Pascal unit
#include <Math.hpp>	// Pascal unit
#include <Dialogs.hpp>	// Pascal unit
#include <SysInit.hpp>	// Pascal unit
#include <System.hpp>	// Pascal unit

//-- user supplied -----------------------------------------------------------

namespace Drag
{
//-- type declarations -------------------------------------------------------
typedef double nla[100][100];

typedef double nly[100];

typedef double nlaux[100][100];

//-- var, const, procedure ---------------------------------------------------
static const Byte korn = 0xc8;
extern PACKAGE double aprl[100][200];
extern PACKAGE double aprr[100][200];
extern PACKAGE double a[100][100];
extern PACKAGE double aap[100][100];
extern PACKAGE double y[100];
extern PACKAGE double yy[100];
extern PACKAGE double y0[100];
extern PACKAGE double prmt[100];
extern PACKAGE double eql[100];
extern PACKAGE double dery[100];
extern PACKAGE double ytrue[100];
extern PACKAGE double R[100];
extern PACKAGE double aux[100][100];
extern PACKAGE int ii;
extern PACKAGE int i;
extern PACKAGE int n;
extern PACKAGE int ihlf;
extern PACKAGE int j0;
extern PACKAGE int jk;
extern PACKAGE int jc;
extern PACKAGE int ntrue;
extern PACKAGE int ng;
extern PACKAGE int nb;
extern PACKAGE double d;
extern PACKAGE double direk;
extern PACKAGE double Kn;
extern PACKAGE bool gg;
extern PACKAGE int mhu[100];
extern PACKAGE int Zu[100];
extern PACKAGE int Zf[100];
extern PACKAGE double sng_real[100];
extern PACKAGE double sng_imag[100];
extern PACKAGE double ys_real[100][100];
extern PACKAGE double ys_imag[100][100];
extern PACKAGE double yr_baz;
extern PACKAGE double yi_baz;
extern PACKAGE void __fastcall f1(double * y);
extern PACKAGE void __fastcall f2(double * aa, double * y);
extern PACKAGE void __fastcall nldet(int n, double * a, double &d);
extern PACKAGE void __fastcall nleqc(int ndim, double &x, double &direk, double * y, double * yy, double 
	* y0, double * dery, double * prmt, double * aa, double * aap);
extern PACKAGE void __fastcall nleqd(int ndim, int ihlf, double &x, double &direk, double * y, double 
	* yy, double * prmt, double * aa, double * aap);
extern PACKAGE void __fastcall nlmain(int ndim, int &ihlf, double * y0, double * dery, double * eql, 
	double * y, double * yy, double * prmt, double * aa, double * aap, double * aux, double &direk);
extern PACKAGE void __fastcall dragmain(void);
extern PACKAGE double __fastcall fi(int i, double y1, double y2, double f1, double f2);

}	/* namespace Drag */
#if !defined(NO_IMPLICIT_NAMESPACE_USE)
using namespace Drag;
#endif
#pragma option pop	// -w-
#pragma option pop	// -Vx

#pragma delphiheader end.
//-- end unit ----------------------------------------------------------------
#endif	// drag
