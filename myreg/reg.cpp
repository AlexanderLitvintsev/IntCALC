#include "reg.h"
#include "ap.h"
#include "sobstvqri.h"
#include "graph_cpp.h"
#include "det.h"
#include "svd.h"
#include "svd_old.h"
Reg *Rasch;
Reg *RRasch;
//---------------------------------------------------------------------------

void Reg::ResultsCalc()
{
int i,j,k;
Double tmpDouble;
complex<double> tmpCompl;

//Считаем перетоки мощностей и токи
ParNumber=0;
for(i=0;i<n-1;i++)
{
 for(j=i+1;j<n;j++)
   {

    if(msw[i][j]==1)
      {
       Params[ParNumber].Node1=i+1;
       Params[ParNumber].Node2=j+1;
       Params[ParNumber].SFlow1=(urab[i]*cktr[i][j]-urab[j])*y[i][j]*conj(urab[i])*cktr[i][j];
//       Params[ParNumber].SFlow2=(urab[i]*cktr[i][j]-urab[j])*y[i][j]*conj(urab[j]);
       Params[ParNumber].SFlow2=(urab[j]-urab[i]*cktr[i][j])*y[i][j]*conj(urab[j]);
       Params[ParNumber].I=(urab[i]*cktr[i][j]-urab[j])*y[i][j]/sqrt(3.);
       ParNumber++;
      }
    else if(msw[i][j]==-1)
      {
       Params[ParNumber].Node1=i+1;
       Params[ParNumber].Node2=j+1;
       Params[ParNumber].SFlow1=(urab[i]-urab[j]/cktr[i][j])*y[j][i]*conj(urab[i]);
//       Params[ParNumber].SFlow2=(urab[i]-urab[j]/cktr[i][j])*y[j][i]*conj(urab[j])/cktr[i][j];
       Params[ParNumber].SFlow2=(urab[j]/cktr[i][j]-urab[i])*y[j][i]*conj(urab[j])/cktr[i][j];
       Params[ParNumber].I=(urab[i]-urab[j]/cktr[i][j])*y[j][i]/sqrt(3.);
       ParNumber++;
      }
    else if(msw[i][j]==2)
      {
       Params[ParNumber].Node1=i+1;
       Params[ParNumber].Node2=j+1;
/*
       if(real(urab[i])>real(urab[j]))
         tmpCompl=urab[i]-urab[j];
       else
         tmpCompl=urab[j]-urab[i];
       Params[ParNumber].SFlow1=tmpCompl*y[min(i,j)][max(i,j)]*conj(urab[i]);
       Params[ParNumber].SFlow2=tmpCompl*y[min(i,j)][max(i,j)]*conj(urab[j]);
*/
       Params[ParNumber].Ql=abs(urab[i])*abs(urab[i])*imag(y[max(i,j)][min(i,j)])+abs(urab[j])*abs(urab[j])*imag(y[max(i,j)][min(i,j)]);
       Params[ParNumber].SFlow1=(urab[i]-urab[j])*y[min(i,j)][max(i,j)]*conj(urab[i])-complex<double>(0,-abs(urab[i])*abs(urab[i])*imag(y[max(i,j)][min(i,j)]));
//     Params[ParNumber].SFlow2=(urab[i]-urab[j])*y[min(i,j)][max(i,j)]*conj(urab[j]);
       Params[ParNumber].SFlow2=(urab[j]-urab[i])*y[min(i,j)][max(i,j)]*conj(urab[j])-complex<double>(0,-abs(urab[j])*abs(urab[j])*imag(y[max(i,j)][min(i,j)]));
       Params[ParNumber].I=(urab[i]-urab[j])*y[min(i,j)][max(i,j)]/sqrt(3.);
//       Params[ParNumber].I=(Params[ParNumber].SFlow2-Params[ParNumber].SFlow1)/conj(urab[i])/sqrt(3.);

       ParNumber++;
      }
   }
}
//здесь потери мощности и напруги в линиях и потери мощностей в трансформаторах
for(i=0;i<n-1;i++)
{
 for(j=i+1;j<n;j++)
   {
    k=GetMsw(i+1,j+1,ParNumber);
    if(msw[i][j]==1)
     {
      Params[k].DSTM=-complex<double>(-real(y[i][j]),imag(y[i][j]))*abs(urab[i]*cktr[i][j]-urab[j])*abs(urab[i]*cktr[i][j]-urab[j]);
//      Params[k].DSTC=-complex<double>(-real(y[j][i]),imag(y[j][i]))*abs(urab[i]*cktr[i][j]-urab[j])*abs(urab[i]*cktr[i][j]-urab[j]);
      Params[k].DSTC=-complex<double>(-real(y[j][i]),imag(y[j][i]))*max(abs(urab[j]),abs(urab[i]))*max(abs(urab[j]),abs(urab[i]));
      Params[k].IsTrf=true;

      }
    else if(msw[i][j]==-1)
     {
      Params[k].DSTM=-complex<double>(-real(y[j][i]),imag(y[j][i]))*abs(urab[i]-urab[j]/cktr[i][j])*abs(urab[i]-urab[j]/cktr[i][j]);
//      Params[k].DSTC=-complex<double>(-real(y[i][j]),imag(y[i][j]))*abs(urab[i]-urab[j]/cktr[i][j])*abs(urab[i]-urab[j]/cktr[i][j]);
      Params[k].DSTC=-complex<double>(-real(y[i][j]),imag(y[i][j]))*max(abs(urab[j]),abs(urab[i]))*max(abs(urab[j]),abs(urab[i]));
      Params[k].IsTrf=true;
     }
    else if(msw[i][j]==2)
     {
//      Params[k].DS=Params[k].SFlow1-Params[k].SFlow2;
      Params[k].DS=-abs((urab[i]-urab[j])*(urab[i]-urab[j]))*complex<double>(-real(y[i][j]),imag(y[i][j]));
      Params[k].DU=-real(urab[i]-urab[j]);
      Params[k].IsLine=true;
     }
   }
}

}
//---------------------------------------------------------------------------
int Reg::GetMsw(int i,int j,int k)
{
int value=-1;
for(int l=0;l<k;l++)
  if(Params[l].Node1==i&&Params[l].Node2==j)
    {
     value=l;
     break;
    }
return value;
}
//---------------------------------------------------------------------------
void Reg::Obnul()
{
ubaz=kt=kl=0;
for (int i=0;i<n;i++)
  {
   mhu[i]=0;
   a0[i]=a1[i]=a2[i]=b0[i]=b1[i]=b2[i]=0;
   sng[i]=unom[i]=urab[i]=irab[i]=complex<double>(0,0);
   for(int j=0;j<n;j++)
     {
      cktr[i][j]=y[i][j]=ys[i][j]=z[i][j]=complex<double>(0,0);
      msw[i][j]=0;
     }
//   for(int j=0;j<n*2;j++)
//      Jac[i][j]=Jac[i+MaxNodesNumber][j]=0;
  }
for (int i=0;i<n*n;i++)
  {
   Params[i].Node1=Params[i].Node2=0;
   Params[i].SFlow1=Params[i].SFlow2=Params[i].DS=Params[i].DSTM=Params[i].DSTC=Params[i].I=complex<double>(0,0);
   Params[i].DU=0;
   Params[i].IsLine=Params[i].IsTrf=false;
  }
}
//---------------------------------------------------------------------------

void Reg::JacobiCalc()
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;

int i,j,k,nb=0,kJac=0,iJac=0;

//double obuslovl;

//for(i=0;i<n;i++)
//   urab[i]=urab[i];
for(i=0;i<n;i++)
 if(mhu[i]==0)
   nb++;
nJac=2*(n-nb);
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdr=dwPdi=dwQdr=dwQdi=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*real(urab[k]);
  dwQdr+=2.*imag(ys[k][k])*real(urab[k]);
  dwPdi+=2.*real(ys[k][k])*imag(urab[k]);
  dwQdi+=2.*imag(ys[k][k])*imag(urab[k]);

  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=real(ys[k][j])*real(urab[j])-imag(ys[k][j])*imag(urab[j]);
//          dwQdr+=-real(ys[k][j])*imag(urab[j])+imag(ys[k][j])*real(urab[j]);
          dwQdr+=real(ys[k][j])*imag(urab[j])+imag(ys[k][j])*real(urab[j]);
          dwPdi+=real(ys[k][j])*imag(urab[j])+imag(ys[k][j])*real(urab[j]);
          dwQdi+=-real(ys[k][j])*real(urab[j])+imag(ys[k][j])*imag(urab[j]);
         }
    }
  }
  else
  {
   dwPdr+=real(ys[k][i])*real(urab[k])+imag(ys[k][i])*imag(urab[k]);
   dwQdr+=-real(ys[k][i])*imag(urab[k])+imag(ys[k][i])*real(urab[k]);
   dwPdi+=real(ys[k][i])*imag(urab[k])-imag(ys[k][i])*real(urab[k]);
   dwQdi+=real(ys[k][i])*real(urab[k])+imag(ys[k][i])*imag(urab[k]);

  }
 Jac[kJac][iJac]=dwPdr;
 Jac[kJac][iJac+n-nb]=dwPdi;
 Jac[kJac+n-nb][iJac]=dwQdr;
 Jac[kJac+n-nb][iJac+n-nb]=dwQdi;
 iJac++;
 }
 kJac++;
 }
kJac=0;
for(k=0;k<n;k++)
if(mhu[k]!=0)
 {
  Nebal[kJac]=0;
  Nebal[kJac+n-nb]=0;
//  Nebal[kJac]+=-real(sng[k])+real(ys[k][k])*(real(urab[k])*real(urab[k])+imag(urab[k])*imag(urab[k]));
// Nebal[kJac+n-nb]+=-imag(sng[k])+imag(ys[k][k])*(real(urab[k])*real(urab[k])+imag(urab[k])*imag(urab[k]));

  tmp1=0;
  tmp2=0;
  tmp1+=real(-sng[k]+ys[k][k]*urab[k]*conj(urab[k]));
  tmp2+=imag(-sng[k]+ys[k][k]*urab[k]*conj(urab[k]));
  for(j=0;j<n;j++)
    if(j!=k)
    {
//     Nebal[kJac]+=real(ys[k][j])*(real(urab[j])*real(urab[k])+imag(urab[j])*imag(urab[k]))-imag(ys[k][j])*(imag(urab[j])*real(urab[k])-real(urab[j])*imag(urab[k]));
//     Nebal[kJac+n-nb]+=real(ys[k][j])*(imag(urab[j])*real(urab[k])-real(urab[j])*imag(urab[k]))+imag(ys[k][j])*(real(urab[j])*real(urab[k])+imag(urab[j])*imag(urab[k]));
     tmp1+=real(ys[k][j]*urab[j]*conj(urab[k]));
     tmp2+=imag(ys[k][j]*urab[j]*conj(urab[k]));
    }
  Nebal[kJac]=tmp1;
  Nebal[kJac+n-nb]=tmp2;
  kJac++;
 }
/*
tmpMatr=new double*[MaxNodesNumber*2];
 for(i=0;i<MaxNodesNumber*2;i++)
  {
   tmpMatr[i]=new double[MaxNodesNumber*4];
  }

ObrMatr(Jac,tmpMatr,nJac);
obuslovl=Norma(Jac,nJac,-1);
obuslovl*=Norma(tmpMatr,nJac,-1);
Form1->Memo1->Lines->Add("Число обусловленности:");
Form1->Memo1->Lines->Add(obuslovl);

//JacobiShow(Form1->Memo1);
 for(i=0;i<MaxNodesNumber*2;i++)
  {
   delete tmpMatr[i];

  }
delete tmpMatr;
tmpMatr=NULL;
*/
}
//---------------------------------------------------------------------------
void Reg::Jacobi2Calc()
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;
bool ogran=true;
int i,j,k,kJac=0,iJac=0;
//double obuslovl;

/*
if (ogran==true)
{
for(k=0;k<n;k++)
 if(mhu[k]==4)
 {
  tmp1=imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
   if(j!=k)
     tmp1+=abs(urab[k])*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
  sng[k];
  if(tmp1<QminQmax[0][k])
   {
   }
 }

}
*/
//for(i=0;i<n;i++)
//   urab[i]=urab[i];

/*
nJac=2*n;
for(i=0;i<n;i++)
 {
 if(mhu[i]==0)
   {
    nJac-=2;
    nb++;
   }

 else if (mhu[i]==4)
   {
    nJac-=1;
    ng++;
   }
 }
*/
//цикл для dP/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwPdr=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwPdr+=abs(urab[k])*(real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 Jac[kJac][iJac]=dwPdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;

//цикл для dQ/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwQdr=0;
  if(i==k)
  {
  dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdr+=abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdr+=abs(urab[k])*(imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 Jac[kJac+n-nb][iJac]=dwQdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;


//цикл для dP/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdi+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwPdi;
         }
    }
  }
  else
  {
   dwPdi+=abs(urab[k])*abs(urab[i])*(real(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
   dwPdi;
  }
 Jac[kJac][iJac+n-nb-ng]=dwPdi;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//цикл для  dQ/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwQdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdi+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdi+=abs(urab[k])*abs(urab[i])*(imag(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  }
 Jac[kJac+n-nb][iJac+n-nb-ng]=dwQdi;
 iJac++;
 }
 kJac++;
 }
//цикл для P
//sng_r[0]=complex<double>(-33.082442385,0);
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
  Nebal[kJac]=0;
  tmp1=0;
  tmp1+=-real(sng_r[k])+real(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp1+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac]=tmp1;
  kJac++;
 }

//цикл для Q
kJac=0;
for(k=0;k<n;k++)
if(mhu[k]!=0&&mhu[k]!=4)
 {
  Nebal[kJac+n]=0;
  tmp2=0;
  tmp2+=-imag(sng_r[k])+imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
    tmp2+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac+n-nb]=tmp2;
  kJac++;
 }


//ObrMatr(Jac,tmpMatr,nJac);
//obuslovl=Norma(Jac,nJac,-1);
//obuslovl*=Norma(tmpMatr,nJac,-1);
//Form1->Memo1->Lines->Add("Число обусловленности:");
//Form1->Memo1->Lines->Add(obuslovl);

//JacobiShow(Form1->Memo1);
}
//---------------------------------------------------------------------------
void Reg::Jacobi_KnCalc()
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;

int i,j,k,kJac=0,iJac=0;

//double obuslovl;

//for(i=0;i<n;i++)
//   urab[i]=urab[i];

/*
nJac=2*n;
for(i=0;i<n;i++)
 {
 if(mhu[i]==0)
   {
    nJac-=2;
    nb++;
   }

 else if (mhu[i]==4)
   {
    nJac-=1;
    ng++;
   }
 }
*/
//цикл для dP/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwPdr=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwPdr+=abs(urab[k])*(real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 Jac[kJac][iJac]=dwPdr*Kn;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;

//цикл для dQ/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwQdr=0;
  if(i==k)
  {
  dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdr+=abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdr+=abs(urab[k])*(imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 Jac[kJac+n-nb][iJac]=dwQdr*Kn;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;


//цикл для dP/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdi+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwPdi;
         }
    }
  }
  else
  {
   dwPdi+=abs(urab[k])*abs(urab[i])*(real(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
   dwPdi;
  }
 Jac[kJac][iJac+n-nb-ng]=dwPdi*Kn;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//цикл для  dQ/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwQdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdi+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdi+=abs(urab[k])*abs(urab[i])*(imag(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  }
 Jac[kJac+n-nb][iJac+n-nb-ng]=dwQdi*Kn;
 iJac++;
 }
 kJac++;
 }
//цикл для P
//sng_r[0]=complex<double>(-33.082442385,0);
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
  Nebal[kJac]=0;
  tmp1=0;
  tmp1+=-real(sng_r[k])+real(ys[k][k])*abs(urab[k])*abs(urab[k])*Kn;
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp1+=Kn*abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac]=tmp1;
  kJac++;
 }

//цикл для Q
kJac=0;
for(k=0;k<n;k++)
if(mhu[k]!=0&&mhu[k]!=4)
 {
  Nebal[kJac+n]=0;
  tmp2=0;
  tmp2+=-imag(sng_r[k])+Kn*imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
    tmp2+=Kn*abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac+n-nb]=tmp2;
  kJac++;
 }
for(k=0;k<nJac;k++)
Jac[k][nJac]=Nebal[k];

//ObrMatr(Jac,tmpMatr,nJac);
//obuslovl=Norma(Jac,nJac,-1);
//obuslovl*=Norma(tmpMatr,nJac,-1);
//Form1->Memo1->Lines->Add("Число обусловленности:");
//Form1->Memo1->Lines->Add(obuslovl);
//JacobiShow(Form1->Memo1);
}
//---------------------------------------------------------------------------



void Reg::JacobiX2Calc()
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;

int i,j,k,nb=0,ng=0,kJac=0,iJac=0;
//int Rb,Rg,Rn;
//double obuslovl;
//double Kn;
//for(i=0;i<n;i++)
//   urab[i]=urab[i];


nJac=2*n;
for(i=0;i<n;i++)
 {
 if(mhu[i]==0)
   {
    nJac-=2;
    nb++;
   }

 else if (mhu[i]==4)
   {
    nJac-=1;
    ng++;
   }
 for(j=nJac/2;j<nJac;j++)
   Jac[i][j]=0;
 }
nJac*=2;

//цикл для dP/dU

for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwPdr=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwPdr+=abs(urab[k])*(real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 dwPdr*=Kn;
 Jac[kJac][iJac]=dwPdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;

//цикл для dQ/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwQdr=0;
  if(i==k)
  {
  dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdr+=abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdr+=abs(urab[k])*(imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 dwQdr*=Kn;
 Jac[kJac+n-nb][iJac]=dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;

//***********************
//цикл для dVm/dU
//Rg=Rb=Rn=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   dwPdr+=2.*real(ys[k][k])*R[k];
   dwQdr+=2.*imag(ys[k][k])*R[k+n];
  }
 else
  {
   dwPdr+=R[k]*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]))+R[i]*fi(1,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
   dwQdr+=R[k+n]*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]))+R[i+n]*fi(1,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;


//dVf/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
      {
       dwPdr+=R[k]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
       dwQdr+=R[k+n]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
      }
    }
//   RCount=0;
   for(j=0;j<n;j++)
    {
     if(j!=k&&mhu[j]!=0)
      {
       dwPdr+=R[j]*abs(urab[j])*fi(3,real(ys[j][k]),-imag(ys[j][k]),arg(urab[j]),arg(urab[k]));
       if(mhu[j]!=4)
         dwQdr+=R[j+n]*abs(urab[j])*fi(3,imag(ys[j][k]),real(ys[j][k]),arg(urab[j]),arg(urab[k]));
      }
//     else if(mhu[j]==0)
//       RCount++;
    }
  }
 else
  {
   dwPdr+=R[k]*abs(urab[k])*fi(2,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[i]*abs(urab[k])*fi(3,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
   dwQdr+=R[k+n]*abs(urab[k])*fi(2,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[i+n]*abs(urab[k])*fi(3,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb-ng][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
 }
//***************************
kJac=0;
iJac=0;

//цикл для dP/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdi+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwPdi+=abs(urab[k])*abs(urab[i])*(real(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  }
 dwPdi*=Kn;
 Jac[kJac][iJac+n-nb-ng]=dwPdi;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//цикл для  dQ/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwQdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdi+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
   dwQdi+=abs(urab[k])*abs(urab[i])*(imag(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  }
 dwQdi*=Kn;
 Jac[kJac+n-nb][iJac+n-nb-ng]=dwQdi;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;
//***********************
//цикл для dVm/df
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
      {
       dwPdr+=R[k]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
       dwQdr+=R[k+n]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
      }
    }
   for(j=0;j<n;j++)
    {
     if(j!=k&&mhu[j]!=0)
      {
       dwPdr+=R[j]*abs(urab[j])*fi(3,real(ys[j][k]),-imag(ys[j][k]),arg(urab[j]),arg(urab[k]));
       if(mhu[j]!=4)
          dwQdr+=R[j+n]*abs(urab[j])*fi(3,imag(ys[j][k]),real(ys[j][k]),arg(urab[j]),arg(urab[k]));
      }
    }
  }
 else
  {
   dwPdr+=R[k]*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[i]*abs(urab[i])*fi(2,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
   dwQdr+=R[k+n]*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[i+n]*abs(urab[i])*fi(2,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;


//dVf/df
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
      {
       dwPdr+=R[k]*abs(urab[k])*abs(urab[j])*fi(4,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
       dwQdr+=R[k+n]*abs(urab[k])*abs(urab[j])*fi(4,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
      }
    }
   for(j=0;j<n;j++)
    {
     if(j!=k&&mhu[j]!=0)
      {
       dwPdr+=R[j]*abs(urab[k])*abs(urab[j])*fi(1,real(ys[j][k]),-imag(ys[j][k]),arg(urab[j]),arg(urab[k]));
       if(mhu[j]!=4)
          dwQdr+=R[j+n]*abs(urab[k])*abs(urab[j])*fi(1,imag(ys[j][k]),real(ys[j][k]),arg(urab[j]),arg(urab[k]));
      }
    }
  }
 else
  {
   dwPdr+=R[k]*abs(urab[k])*abs(urab[i])*fi(4,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[i]*abs(urab[k])*abs(urab[i])*fi(1,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
   dwQdr+=R[k+n]*abs(urab[k])*abs(urab[i])*fi(4,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[i+n]*abs(urab[k])*abs(urab[i])*fi(1,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb-ng][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
 }
//***************************
kJac=0;
iJac=0;

//dVm/dRp
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 if(k==i)
  {
   dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwPdr+=abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwPdr+=abs(urab[i])*fi(1,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
 }
 dwPdr*=Kn;
 Jac[kJac+nJac/2][iJac+nJac/2]=dwPdr;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;

//dVf/dRp
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 if(k==i)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwPdr+=abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwPdr+=abs(urab[k])*abs(urab[i])*fi(3,real(ys[i][k]),-imag(ys[i][k]),arg(urab[i]),arg(urab[k]));
 }
 dwPdr*=Kn;
 Jac[kJac+nJac/2+n-nb-ng][iJac+nJac/2]=dwPdr;
 iJac++;
 }
 kJac++;
 }


iJac=0;
kJac=0;
//dVm/dRq
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[i]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwQdr=0;
 if(k==i)
  {
   dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwQdr+=abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwQdr+=abs(urab[i])*fi(1,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
 }
 dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac+nJac/2+n-nb]=dwQdr;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;

//dVf/dRq
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwQdr=0;
 if(k==i)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwQdr+=abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwQdr+=abs(urab[k])*abs(urab[i])*fi(3,imag(ys[i][k]),real(ys[i][k]),arg(urab[i]),arg(urab[k]));
 }
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb-ng][iJac+nJac/2+n-nb]=dwQdr;
 iJac++;
 }
 kJac++;
 }

//цикл для P
//sng_r[0]=complex<double>(-33.082442385,0);
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
//  Nebal[kJac]=0;
  tmp1=0;
  tmp2=0;
  tmp1=-real(sng_r[k]);
  tmp2+=real(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp2+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Jac[kJac][nJac-1]=tmp2;
  Nebal[kJac]=tmp1+tmp2*Kn;
 // Nebal[kJac]=tmp1+tmp2;
  kJac++;
 }

//цикл для Q
kJac=0;
for(k=0;k<n;k++)
if(mhu[k]!=0&&mhu[k]!=4)
 {
//  Nebal[kJac+n-nb]=0;
  tmp2=0;
  tmp1=0;
  tmp1=-imag(sng_r[k]);
  tmp2+=imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
    tmp2+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Jac[kJac+n-nb][nJac-1]=tmp2;
  Nebal[kJac+n-nb]=tmp1+tmp2*Kn;
//  Nebal[kJac+n-nb]=tmp1+tmp2;
  kJac++;
 }

//цикл для Vm
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 tmp1=2.*real(ys[k][k])*abs(urab[k])*R[k];
 tmp2=2.*imag(ys[k][k])*abs(urab[k])*R[k+n-nb];
 for(j=0;j<n;j++)
  {
   if(j!=k)
    {
     tmp1+=R[k]*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     tmp2+=R[k+n-nb]*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 for(j=0;j<n;j++)
  {
   if(j!=k&&mhu[j]!=0)
    {
     tmp1+=R[j]*abs(urab[j])*fi(1,real(ys[j][k]),-imag(ys[j][k]),arg(urab[j]),arg(urab[k]));
     if(mhu[j]!=4)
        tmp2+=R[j+n-nb]*abs(urab[j])*fi(1,imag(ys[j][k]),real(ys[j][k]),arg(urab[j]),arg(urab[k]));
    }
  }
 Jac[nJac/2+kJac][nJac-1]=tmp1+tmp2;
 Nebal[nJac/2+kJac]=(tmp1*Kn+tmp2*Kn);
 kJac++;
 }

//Vf
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 tmp1=0;
 tmp2=0;
 for(j=0;j<n;j++)
  {
   if(j!=k)
    {
     tmp1+=R[k]*abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     tmp2+=R[k+n-nb]*abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 for(j=0;j<n;j++)
  {
   if(j!=k&&mhu[j]!=0)
    {
     tmp1+=R[j]*abs(urab[k])*abs(urab[j])*fi(3,real(ys[j][k]),-imag(ys[j][k]),arg(urab[j]),arg(urab[k]));
     if(mhu[j]!=4)
        tmp2+=R[j+n-nb]*abs(urab[k])*abs(urab[j])*fi(3,imag(ys[j][k]),real(ys[j][k]),arg(urab[j]),arg(urab[k]));
    }
  }
 Jac[kJac+nJac/2+n-nb-ng][nJac-1]=tmp1+tmp2;
 Nebal[kJac+nJac/2+n-nb-ng]=(Kn*tmp1+Kn*tmp2);
 kJac++;
 }
for(i=0;i<nJac;i++)
{
 Jac[i][nJac]=Nebal[i];
}



//ObrMatr(Jac,tmpMatr,nJac);
//obuslovl=Norma(Jac,nJac,-1);
//obuslovl*=Norma(tmpMatr,nJac,-1);
//Form1->Memo1->Lines->Add("Число обусловленности:");
//Form1->Memo1->Lines->Add(obuslovl);
//JacobiShow(Form1->Memo1);
}
//---------------------------------------------------------------------------
void Reg::JacobiX2TCalc()//расширенная матрица якоби - с параметром Т
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;

int i,j,k,l,kJac=0,iJac=0;
//double obuslovl;


for(i=0;i<nJac/2;i++)
 for(j=nJac/2;j<nJac+1;j++)
  Jac[i][j]=0;

//цикл для dP/dU

for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwPdr=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
   dwPdr+=abs(urab[k])*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
//   dwPdr+=abs(urab[k])*(real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
// dwPdr*=Kn;
 Jac[kJac][iJac]=dwPdr;
 Jac[kJac+nJac/2][iJac+nJac/2]=dwPdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;

//цикл для dQ/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwQdr=0;
  if(i==k)
  {
  dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdr+=abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
//          dwQdr+=abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
//   dwQdr+=abs(urab[k])*(imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
   dwQdr+=abs(urab[k])*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwQdr*=Kn;
 Jac[kJac+n-nb][iJac]=dwQdr;
 Jac[kJac+n-nb+nJac/2][iJac+nJac/2]=dwQdr;
 iJac++;
 }
 kJac++;
}
iJac=0;
kJac=0;

//***********************
//цикл для dVp/dU

for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   if(mhu[k]!=4)
   {
    dwPdr+=2.*real(ys[k][k])*R[Zu[k]];
    for(j=0;j<n;j++)
     if(j!=k)
      dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    for(l=0;l<n;l++)
     if(l!=k&&mhu[l]!=0)
      {
       dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
       if(mhu[l]!=4)
        dwPdr+= R[Zu[l]]*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
      }
   }
   else //для генераторного узла
   {
    for(j=0;j<n;j++)
     if(j!=k)
      dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    for(l=0;l<n;l++)
    {
     if(mhu[l]!=0&&mhu[l]!=4)
       dwPdr+= R[Zu[l]]*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
     if(l!=k&&mhu[l]!=0)
       dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
    }
   }
  }
 else
  {
   if(mhu[k]!=4)
     dwPdr+= R[Zu[k]]*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*fi(2,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwPdr*=Kn;
// dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//dVq/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   dwPdr+=2.*imag(ys[k][k])*R[Zu[k]];
   for(j=0;j<n;j++)
    if(j!=k)
     dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
   for(l=0;l<n;l++)
    if(l!=k&&mhu[l]!=0)
     {
      dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
      if(mhu[l]!=4)
       dwPdr+= R[Zu[l]]*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
     }
  }
 else
  {
   dwPdr+= R[Zu[k]]*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*fi(2,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwPdr*=Kn;
// dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;

//цикл для dP/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
//          dwPdi+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwPdi+=abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
//dwPdi+=abs(urab[k])*abs(urab[i])*(real(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  dwPdi+=abs(urab[k])*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwPdi*=Kn;
 Jac[kJac][iJac+n-nb-ng]=dwPdi;
 Jac[kJac+nJac/2][iJac+n-nb-ng+nJac/2]=dwPdi;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//цикл для  dQ/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwQdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
//          dwQdi+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwQdi+=abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
   //dwQdi+=abs(urab[k])*abs(urab[i])*(imag(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
   dwQdi+=abs(urab[k])*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwQdi*=Kn;
 Jac[kJac+n-nb][iJac+n-nb-ng]=dwQdi;
 Jac[kJac+n-nb+nJac/2][iJac+n-nb-ng+nJac/2]=dwQdi;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;
//***********************
//цикл для dVp/df
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   if(mhu[k]!=4)
   {
    for(j=0;j<n;j++)
     {
      if(j!=k)
       {
        dwPdr+=R[Zu[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
        dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
       }
     }
    for(l=0;l<n;l++)
     {
      if(l!=k&&mhu[l]!=0)
       {
        dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
        if(mhu[l]!=4)
           dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
       }
     }
   }
   else //для генераторных узлов
   {
    for(j=0;j<n;j++)
     {
      if(j!=k)
        dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    for(l=0;l<n;l++)
     {
      if(mhu[l]!=0&&mhu[l]!=4)
        dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
      if(l!=k&&mhu[l]!=0)
        dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
     }
   }

  }
 else
  {
   if(mhu[k]!=4)
     dwPdr+=R[Zu[k]]*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[i])*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[Zu[i]]*abs(urab[k])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*abs(urab[i])*fi(4,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwPdr*=Kn;
// dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;


//dVq/df
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
      {
       dwPdr+=R[Zu[k]]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
       dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
      }
    }
   for(l=0;l<n;l++)
    {
     if(l!=k&&mhu[l]!=0)
      {
       dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
       if(mhu[l]!=4)
          dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
      }
    }
  }
 else
  {
   dwPdr+=R[Zu[k]]*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[i])*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[Zu[i]]*abs(urab[k])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*abs(urab[i])*fi(4,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
// dwPdr*=Kn;
// dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}
//***************************
kJac=0;
iJac=0;

//цикл для P
//sng_r[0]=complex<double>(-33.082442385,0);
kJac=0;
//Form1->Memo1->Lines->Add("dY");
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
//  Nebal[kJac]=0;
  tmp1=0;
  tmp2=0;
  tmp1+=-dPn[k]+dPg[k];
//  Form1->Memo1->Lines->Add(AnsiString(dPn[k]-dPg[k]));
  tmp2+=-real(sng_r[k])+real(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp2+=abs(urab[k])*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
//     tmp2+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac]=tmp1*T+tmp2;
  Jac[kJac][nJac-1]=tmp1;
  kJac++;
 }

//цикл для Q
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
//  Nebal[kJac]=0;
  tmp1=0;
  tmp2=0;
  tmp1+=dQn[k];
//  Form1->Memo1->Lines->Add(AnsiString(dQn[k]));
  tmp2+=-imag(sng_r[k])+imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp2+=abs(urab[k])*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
//     tmp2+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac+n-nb]=tmp1*T+tmp2;
  Jac[kJac+n-nb][nJac-1]=tmp1;
  kJac++;
 }

//цикл для Vp
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
  tmp1=0;
  tmp2=0;
  for(l=0;l<n;l++)
   if(mhu[l]!=0&&mhu[l]!=4)
   {
    if(k==l)
     {
      tmp1+=2.*real(ys[k][k])*abs(urab[k])*R[Zu[k]];
      for(j=0;j<n;j++)
       if(j!=k)
         tmp1+=R[Zu[l]]*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
  for(l=0;l<n;l++)
   if(mhu[l]!=0)
   {
    if(k==l)
     {
      for(j=0;j<n;j++)
       if(j!=k)
         tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
// Jac[nJac/2+kJac][nJac-1]=tmp1+tmp2;
 Jac[nJac/2+kJac][nJac-1]=0;
 Nebal[nJac/2+kJac]=(tmp1+tmp2);
 kJac++;
 }


//Vq
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
  tmp1=0;
  tmp2=0;
  for(l=0;l<n;l++)
   if(mhu[l]!=0&&mhu[l]!=4)
   {
    if(k==l)
     {
      tmp1+=2.*imag(ys[k][k])*abs(urab[k])*R[Zu[k]];
      for(j=0;j<n;j++)
       if(j!=k)
         tmp1+=R[Zu[l]]*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
  for(l=0;l<n;l++)
   if(mhu[l]!=0)
   {
    if(k==l)
     {
      for(j=0;j<n;j++)
       if(j!=k)
         tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
// Jac[nJac/2+kJac+n-nb][nJac-1]=tmp1+tmp2;
 Jac[nJac/2+kJac+n-nb][nJac-1]=0;
 Nebal[nJac/2+kJac+n-nb]=(tmp1+tmp2);
 kJac++;
 }
for(i=0;i<nJac;i++)
{
 Jac[i][nJac]=Nebal[i];
}
//JacobiShow(Form1->Memo1);
}
//---------------------------------------------------------------------------
void Reg::JacobiX22Calc()//расширенная матрица якоби - моя
{
double dwPdr,dwPdi,dwQdr,dwQdi;
dwPdr=dwPdi=dwQdr=dwQdi=0;
double tmp1,tmp2;

int i,j,k,l,kJac=0,iJac=0;
//double obuslovl;

for(i=0;i<nJac/2;i++)
 for(j=nJac/2;j<nJac+1;j++)
  Jac[i][j]=0;

//цикл для dP/dU

for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwPdr=0;
  if(i==k)
  {
  dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwPdr+=abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
   dwPdr+=abs(urab[k])*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
//   dwPdr+=abs(urab[k])*(real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
  }
 dwPdr*=Kn;
 Jac[kJac][iJac]=dwPdr;
 Jac[kJac+nJac/2][iJac+nJac/2]=dwPdr;
 iJac++;
 }
 kJac++;
 }

iJac=0;
kJac=0;

//цикл для dQ/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
  dwQdr=0;
  if(i==k)
  {
  dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
          dwQdr+=abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
//          dwQdr+=abs(urab[j])*(imag(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
         }
    }
  }
  else
  {
//   dwQdr+=abs(urab[k])*(imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*sin(arg(urab[i]))-sin(arg(urab[k]))*cos(arg(urab[i]))));
   dwQdr+=abs(urab[k])*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwQdr*=Kn;
 Jac[kJac+n-nb][iJac]=dwQdr;
 Jac[kJac+n-nb+nJac/2][iJac+nJac/2]=dwQdr;
 iJac++;
 }
 kJac++;
}
iJac=0;
kJac=0;

//***********************
//цикл для dVp/dU

for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   if(mhu[k]!=4)
   {
    dwPdr+=2.*real(ys[k][k])*R[Zu[k]];
    for(j=0;j<n;j++)
     if(j!=k)
      dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    for(l=0;l<n;l++)
     if(l!=k&&mhu[l]!=0)
      {
       dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
       if(mhu[l]!=4)
        dwPdr+= R[Zu[l]]*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
      }
   }
   else //для генераторного узла
   {
    for(j=0;j<n;j++)
     if(j!=k)
      dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    for(l=0;l<n;l++)
    {
     if(mhu[l]!=0&&mhu[l]!=4)
       dwPdr+= R[Zu[l]]*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
     if(l!=k&&mhu[l]!=0)
       dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
    }
   }
  }
 else
  {
   if(mhu[k]!=4)
     dwPdr+= R[Zu[k]]*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*fi(2,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//dVq/dU
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   dwPdr+=2.*imag(ys[k][k])*R[Zu[k]];
   for(j=0;j<n;j++)
    if(j!=k)
     dwQdr+=R[Zf[k]]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
   for(l=0;l<n;l++)
    if(l!=k&&mhu[l]!=0)
     {
      dwQdr+=R[Zf[l]]*abs(urab[l])*fi(3,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
      if(mhu[l]!=4)
       dwPdr+= R[Zu[l]]*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
     }
  }
 else
  {
   dwPdr+= R[Zu[k]]*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*fi(2,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;

//цикл для dP/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwPdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
//          dwPdi+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwPdi+=abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
//dwPdi+=abs(urab[k])*abs(urab[i])*(real(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))-imag(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
  dwPdi+=abs(urab[k])*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwPdi*=Kn;
 Jac[kJac][iJac+n-nb-ng]=dwPdi;
 Jac[kJac+nJac/2][iJac+n-nb-ng+nJac/2]=dwPdi;
 iJac++;
 }
 kJac++;
 }
iJac=0;
kJac=0;


//цикл для  dQ/dD
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
  dwQdi=0;
  if(i==k)
  {
  for(j=0;j<n;j++)
    {
        if(j!=k)
         {
//          dwQdi+=abs(urab[k])*abs(urab[j])*(imag(ys[k][j])*(-sin(arg(urab[k]))*cos(arg(urab[j]))+cos(arg(urab[k]))*sin(arg(urab[j])))+real(ys[k][j])*(-sin(arg(urab[k]))*sin(arg(urab[j]))-cos(arg(urab[k]))*cos(arg(urab[j]))));
          dwQdi+=abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
         }
    }
  }
  else
  {
   //dwQdi+=abs(urab[k])*abs(urab[i])*(imag(ys[k][i])*(-cos(arg(urab[k]))*sin(arg(urab[i]))+sin(arg(urab[k]))*cos(arg(urab[i])))+real(ys[k][i])*(cos(arg(urab[k]))*cos(arg(urab[i]))+sin(arg(urab[k]))*sin(arg(urab[i]))));
   dwQdi+=abs(urab[k])*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwQdi*=Kn;
 Jac[kJac+n-nb][iJac+n-nb-ng]=dwQdi;
 Jac[kJac+n-nb+nJac/2][iJac+n-nb-ng+nJac/2]=dwQdi;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;
//***********************
//цикл для dVp/df
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   if(mhu[k]!=4)
   {
    for(j=0;j<n;j++)
     {
      if(j!=k)
       {
        dwPdr+=R[Zu[k]]*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
        dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
       }
     }
    for(l=0;l<n;l++)
     {
      if(l!=k&&mhu[l]!=0)
       {
        dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
        if(mhu[l]!=4)
           dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
       }
     }
   }
   else //для генераторных узлов
   {
    for(j=0;j<n;j++)
     {
      if(j!=k)
        dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    for(l=0;l<n;l++)
     {
      if(mhu[l]!=0&&mhu[l]!=4)
        dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
      if(l!=k&&mhu[l]!=0)
        dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
     }
   }

  }
 else
  {
   if(mhu[k]!=4)
     dwPdr+=R[Zu[k]]*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[i])*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[Zu[i]]*abs(urab[k])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*abs(urab[i])*fi(4,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}

iJac=0;
kJac=0;


//dVq/df
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 dwQdr=0;
 if(i==k)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
      {
       dwPdr+=R[Zu[k]]*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
       dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(4,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
      }
    }
   for(l=0;l<n;l++)
    {
     if(l!=k&&mhu[l]!=0)
      {
       dwQdr+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
       if(mhu[l]!=4)
          dwPdr+=R[Zu[l]]*abs(urab[k])*fi(2,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
      }
    }
  }
 else
  {
   dwPdr+=R[Zu[k]]*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[k]]*abs(urab[k])*abs(urab[i])*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwPdr+=R[Zu[i]]*abs(urab[k])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
   dwQdr+=R[Zf[i]]*abs(urab[k])*abs(urab[i])*fi(4,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
  }
 dwPdr*=Kn;
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac+n-nb-ng]=dwPdr+dwQdr;
 iJac++;
 }
 kJac++;
}
//***************************
kJac=0;
iJac=0;

/*
//dVp/dRu
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwPdr=0;
 if(k==i)
  {
   dwPdr+=2.*real(ys[k][k])*abs(urab[k]);
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwPdr+=abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwPdr+=abs(urab[k])*fi(1,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
 }
 dwPdr*=Kn;
 Jac[kJac+nJac/2][iJac+nJac/2]=dwPdr;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;

//dVQ/dRu
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 dwQdr=0;
 if(k==i)
  {
   dwQdr+=2.*imag(ys[k][k])*abs(urab[k]);
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwQdr+=abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwQdr+=abs(urab[k])*fi(1,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
 }
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac+nJac/2]=dwQdr;
 iJac++;
 }
 kJac++;
 }


iJac=0;
kJac=0;
//dVp/dRf
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwPdr=0;
 if(k==i)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwPdr+=abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwPdr+=abs(urab[k])*abs(urab[i])*fi(3,real(ys[k][i]),-imag(ys[k][i]),arg(urab[k]),arg(urab[i]));
 }
 dwPdr*=Kn;
 Jac[kJac+nJac/2][iJac+nJac/2+n-nb-ng]=dwPdr;
 iJac++;
 }
 kJac++;
 }

kJac=0;
iJac=0;

//dVq/dRf
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
 iJac=0;
 for(i=0;i<n;i++)
 if(mhu[i]!=0)
 {
 dwQdr=0;
 if(k==i)
  {
   for(j=0;j<n;j++)
    {
     if(j!=k)
       dwQdr+=abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 else
 {
  dwQdr+=abs(urab[k])*abs(urab[i])*fi(3,imag(ys[k][i]),real(ys[k][i]),arg(urab[k]),arg(urab[i]));
 }
 dwQdr*=Kn;
 Jac[kJac+nJac/2+n-nb][iJac+nJac/2+n-nb-ng]=dwQdr;
 iJac++;
 }
 kJac++;
 }
*/
//цикл для P
//sng_r[0]=complex<double>(-33.082442385,0);
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
//  Nebal[kJac]=0;
  tmp1=0;
  tmp2=0;
  tmp1+=-real(sng_r[k]);
  tmp2+=real(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp2+=abs(urab[k])*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
//     tmp2+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac]=tmp1+tmp2*Kn;
  Jac[kJac][nJac-1]=tmp2;
  kJac++;
 }

//цикл для Q
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
//  Nebal[kJac]=0;
  tmp1=0;
  tmp2=0;
  tmp1+=-imag(sng_r[k]);
  tmp2+=imag(ys[k][k])*abs(urab[k])*abs(urab[k]);
  for(j=0;j<n;j++)
    if(j!=k)
    {
     tmp2+=abs(urab[k])*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
//     tmp2+=abs(urab[k])*abs(urab[j])*(real(ys[k][j])*(cos(arg(urab[k]))*cos(arg(urab[j]))+sin(arg(urab[k]))*sin(arg(urab[j])))-imag(ys[k][j])*(cos(arg(urab[k]))*sin(arg(urab[j]))-sin(arg(urab[k]))*cos(arg(urab[j]))));
    }
  Nebal[kJac+n-nb]=tmp1+tmp2*Kn;
  Jac[kJac+n-nb][nJac-1]=tmp2;
  kJac++;
 }

//цикл для Vp
/*
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
 tmp1=2.*real(ys[k][k])*abs(urab[k])*R[Zu[k]];
 tmp2=0;
  for(j=0;j<n;j++)
  {
   if(j!=k)
    {
     tmp1+=R[Zu[k]]*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     tmp2+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 for(l=0;l<n;l++)
  {
   if(l!=k&&mhu[l]!=0)
    {
     tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
     if(mhu[l]!=4)
        tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
    }
  }
 Jac[nJac/2+kJac][nJac-1]=tmp1+tmp2;
 Nebal[nJac/2+kJac]=(tmp1*Kn+tmp2*Kn);
 kJac++;
 }
*/
//цикл для Vp
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0)
 {
  tmp1=0;
  tmp2=0;
  for(l=0;l<n;l++)
   if(mhu[l]!=0&&mhu[l]!=4)
   {
    if(k==l)
     {
      tmp1+=2.*real(ys[k][k])*abs(urab[k])*R[Zu[k]];
      for(j=0;j<n;j++)
       if(j!=k)
         tmp1+=R[Zu[l]]*abs(urab[j])*fi(1,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
  for(l=0;l<n;l++)
   if(mhu[l]!=0)
   {
    if(k==l)
     {
      for(j=0;j<n;j++)
       if(j!=k)
         tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[j])*fi(2,real(ys[k][j]),-imag(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,real(ys[k][l]),-imag(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
 Jac[nJac/2+kJac][nJac-1]=tmp1+tmp2;
 Nebal[nJac/2+kJac]=(tmp1*Kn+tmp2*Kn);
 kJac++;
 }


//Vq
/*
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
  tmp1=2.*imag(ys[k][k])*abs(urab[k])*R[Zu[k]];
  tmp2=0;
  for(j=0;j<n;j++)
  {
   if(j!=k)
    {
     tmp1+=R[Zu[k]]*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
     tmp2+=R[Zf[k]]*abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
    }
  }
 for(l=0;l<n;l++)
  {
   if(l!=k&&mhu[l]!=0)
    {
     tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
     if(mhu[j]!=4)
        tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
    }
  }
 Jac[nJac/2+kJac+n-nb][nJac-1]=tmp1+tmp2;
 Nebal[nJac/2+kJac+n-nb]=(tmp1*Kn+tmp2*Kn);
 kJac++;
 }
*/
//Vq
kJac=0;
for(k=0;k<n;k++)
 if(mhu[k]!=0&&mhu[k]!=4)
 {
  tmp1=0;
  tmp2=0;
  for(l=0;l<n;l++)
   if(mhu[l]!=0&&mhu[l]!=4)
   {
    if(k==l)
     {
      tmp1+=2.*imag(ys[k][k])*abs(urab[k])*R[Zu[k]];
      for(j=0;j<n;j++)
       if(j!=k)
         tmp1+=R[Zu[l]]*abs(urab[j])*fi(1,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp1+=R[Zu[l]]*abs(urab[k])*fi(1,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
  for(l=0;l<n;l++)
   if(mhu[l]!=0)
   {
    if(k==l)
     {
      for(j=0;j<n;j++)
       if(j!=k)
         tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[j])*fi(2,imag(ys[k][j]),real(ys[k][j]),arg(urab[k]),arg(urab[j]));
     }
    else
      tmp2+=R[Zf[l]]*abs(urab[k])*abs(urab[l])*fi(3,imag(ys[k][l]),real(ys[k][l]),arg(urab[k]),arg(urab[l]));
   }
 Jac[nJac/2+kJac+n-nb][nJac-1]=tmp1+tmp2;
 Nebal[nJac/2+kJac+n-nb]=(tmp1*Kn+tmp2*Kn);
 kJac++;
 }
for(i=0;i<nJac;i++)
{
 Jac[i][nJac]=Nebal[i];
}
//JacobiShow(Form1->Memo1);
}
//---------------------------------------------------------------------------






int Reg::Gauss()
{
int i,j,k,state;
double tmp,tmp2[MaxNodesNumber];
double tmpJac[MaxNodesNumber][MaxNodesNumber];
complex<double>tmpCompl;
double errsumm;
bool ok;
bool RegOk;
AnsiString tmpStr;
for(i=0;i<nJac;i++)
{
 Jac[i][nJac]=Nebal[i];
}
//JacobiShow(Form1->Memo1);
/*
Form1->Memo1->Lines->Add("***********");
for(i=0;i<nJac;i++)
{
tmpStr="";
for(j=0;j<nJac;j++)
  {
   tmpJac[i][j]=Jac[i][j];
   tmpStr+=FloatToStrF(tmpJac[i][j],ffFixed,3,7)+"; ";

  }
Form1->Memo1->Lines->Add(tmpStr);
}
*/
ok=ReorganizeJacobi();
if(ok==false)
  {
   ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
   goto WrongJac;
  }

//здесь метод без обратного хода

for(i=0;i<nJac;i++)
{
 tmp=Jac[i][i];
 if(tmp==0)
   ok=ReorganizeJacobi();
 if(ok==false)
   {
    ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
    goto WrongJac;
   }
 tmp=Jac[i][i];
 for(j=0;j<nJac+1;j++)
   Jac[i][j]=Jac[i][j]/tmp;
 for(k=0;k<nJac;k++)
   {
    tmp=Jac[k][i];
    for(j=0;j<nJac+1;j++)
     {
      if(k!=i)
         {
          Jac[k][j]=Jac[k][j]-Jac[i][j]*tmp;
         }
     }
   }
}
tmp=0;
for(i=0;i<nJac;i++)
  tmp2[i]=0;


//здесь метод с обратным ходом
/*
for(i=0;i<nJac-1;i++)
{

 for(k=i+1;k<nJac;k++)
   {
    tmp=Jac[k][i];
    for(j=i;j<nJac;j++)
      {
      Jac[k][j]=Jac[k][j]-tmp*Jac[i][j]/Jac[i][i];
      Jac[k][j];
      }
   }
// else ShowMessage("Матрица Якоби не верна!");

}
*/
k=0;
errsumm=0;
RegOk=true;
for(i=0;i<nJac/2;i++)
  {
   if(mhu[i]==0)
     k++;
//   tmpCompl=polar(Jac[i][nJac],Jac[i+nJac/2][nJac]);
   tmpCompl=complex<double>(Jac[i][nJac],Jac[i+nJac/2][nJac]);
//   urab[i+1]=complex<double>(abs(urab[i+1])-Jac[i][nJac],arg(urab[i+1])-Jac[i+nJac/2][nJac]);
   urab[i+k]=urab[i+k]-tmpCompl;
//   errsumm+=abs(tmpCompl);
   errsumm=abs(tmpCompl);
   if(errsumm>=eps)
     RegOk=false;
   tmp=abs(urab[i+k]);
   tmp;
   tmp=arg(urab[i+k])*180./M_PI;
   tmp;
  }
errsumm=errsumm/double(nJac/2);
//if(errsumm<eps)
if(RegOk==true)
  state=0;//сошелся!
else
  state=1;

goto SuccesfulJac;
WrongJac:
state=2;
SuccesfulJac:
return state;
}
//---------------------------------------------------------------------------

int Reg::Gauss2()
{
int i,j,k,state,nb,ng;
double tmp,tmp2[MaxNodesNumber];
double tmpJac[MaxNodesNumber][MaxNodesNumber];
complex<double>tmpCompl;
double errsumm;
bool ok;
bool RegOk;
AnsiString tmpStr;
//double *dU;double *dD;

ok=true;
for(i=0;i<nJac;i++)
{
 Jac[i][nJac]=Nebal[i];
}

ok=ReorganizeJacobi();
if(ok==false)
  {
   ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
   goto WrongJac;
  }

//здесь метод без обратного хода

for(i=0;i<nJac;i++)
{
 tmp=Jac[i][i];
 if(tmp==0)
   ok=ReorganizeJacobi();
 if(ok==false)
   {
    ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
    goto WrongJac;
   }
 tmp=Jac[i][i];
 for(j=0;j<nJac+1;j++)
   Jac[i][j]=Jac[i][j]/tmp;
 for(k=0;k<nJac;k++)
   {
    tmp=Jac[k][i];
    for(j=0;j<nJac+1;j++)
     {
      if(k!=i)
         {
          Jac[k][j]=Jac[k][j]-Jac[i][j]*tmp;
         }
     }
   }
}
tmp=0;
for(i=0;i<nJac;i++)
  tmp2[i]=0;


//здесь метод с обратным ходом
/*
for(i=0;i<nJac2-1;i++)
{

 for(k=i+1;k<nJac2;k++)
   {
    tmp=Jac[k][i];
    for(j=i;j<nJac2;j++)
      {
      Jac[k][j]=Jac[k][j]-tmp*Jac[i][j]/Jac[i][i];
      Jac[k][j];
      }
   }
// else ShowMessage("Матрица Якоби не верна!");

}
*/
k=0;
errsumm=0;
RegOk=true;

//цикл для U


dU=new double [MaxNodesNumber];
dD=new double [MaxNodesNumber];

nb=ng=0;
for(i=0;i<n;i++)
  if(mhu[i]==0)
    nb++;
for(i=0;i<n;i++)
  if(mhu[i]==4)
    ng++;


j=0;
for(i=0;i<n;i++)
 {
  if(mhu[i]==0||mhu[i]==4)
    {
     dU[i]=0;
     j++;
    }
  else
     dU[i]=Jac[i-j][nJac];
 }
j=0;
for(i=0;i<n;i++)
 {
  if(mhu[i]==0)
    {
     dD[i]=0;
     j++;
    }
  else
     dD[i]=Jac[i+n-nb-ng-j][nJac];
 }

for(i=0;i<n;i++)
  {
   tmpCompl=complex<double>(dU[i],dD[i]);
   urab[i]=complex<double>(abs(urab[i]),arg(urab[i]))-tmpCompl;
   urab[i]=polar(real(urab[i]),imag(urab[i]));
   errsumm=fabs(dU[i])+fabs(dD[i]);
   if(errsumm>=eps)
     RegOk=false;
  }
//if(errsumm<eps)

for(i=0;i<nJac;i++)
  if(fabs(Nebal[i])>eps)
    RegOk=false;

if(RegOk==true)
  state=0;//сошелся!
else
  state=1;
delete dD;
delete dU;
goto SuccesfulJac;
WrongJac:
state=2;
SuccesfulJac:
return state;
}

//---------------------------------------------------------------------------
void Reg::Gauss(complex<double> **Matr,complex<double> *Rez,int razm,int *uzlNo)
{
int i,j,k;
complex<double> tmp;
bool ok=true;
ok=ReorganizeMatr(Matr,razm);
if(ok==false)
 {
  ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
  goto WrongMatr;
 }
for(i=0;i<razm;i++)
{
 tmp=Matr[uzlNo[i]][uzlNo[i]];
 if(tmp==complex<double>(0,0))
   ok=ReorganizeMatr(Matr,razm);
 if(ok==false)
   {
    ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
    goto WrongMatr;
   }
 tmp=Matr[uzlNo[i]][uzlNo[i]];
 for(j=0;j<razm+1;j++)
   Matr[uzlNo[i]][uzlNo[j]]=Matr[uzlNo[i]][uzlNo[j]]/tmp;
 for(k=0;k<razm;k++)
   {
    tmp=Matr[uzlNo[k]][uzlNo[i]];
    for(j=0;j<razm+1;j++)
     {
      if(uzlNo[k]!=uzlNo[i])
         {
          Matr[uzlNo[k]][uzlNo[j]]=Matr[uzlNo[k]][uzlNo[j]]-Matr[uzlNo[i]][uzlNo[j]]*tmp;
         }
     }
   }
}
for(i=0;i<razm;i++)
Rez[uzlNo[i]]=Matr[uzlNo[i]][uzlNo[razm]]/Matr[uzlNo[i]][uzlNo[i]];
WrongMatr:
}
//---------------------------------------------------------------------------
bool Reg::Gauss(double **Matr,double *Rez,int razm)
{
int i,j,k;
double tmp;
bool ok;

bool state=false;
ok=true;
ok=ReorganizeMatr(Matr,razm);
if(ok==false)
  {
   ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
   goto WrongMatr1;
  }

//здесь метод без обратного хода

for(i=0;i<razm;i++)
{
 tmp=Matr[i][i];
 if(tmp==0)
   ok=ReorganizeMatr(Matr,razm);
 if(ok==false)
   {
    ShowMessage("Ошибка в вычислении корней - проверьте исх. данные");
    goto WrongMatr1;
   }
 tmp=Matr[i][i];
 for(j=0;j<razm+1;j++)
   Matr[i][j]=Matr[i][j]/tmp;
 for(k=0;k<razm;k++)
   {
    tmp=Matr[k][i];
      if(k!=i)
      {
       for(j=0;j<razm+1;j++)
         {
          Matr[k][j]=Matr[k][j]-Matr[i][j]*tmp;
         }
      }
   }
}
for(i=0;i<razm;i++)
  Rez[i]=Matr[i][razm]/Matr[i][i];
state=true;
WrongMatr1:
return state;
}
//---------------------------------------------------------------------------
bool Reg::RaschCheck(double *Rez,complex<double> *urab,int razm,int NoKn,int Raschtype,bool LimitNebals,double dUmax,double dfmax,double dKnmax,double dTmax,double dVri)
{
double err;
int NoKn1;
dU=new double [n];
dD=new double [n];

int i,k;
bool RegOk=true;

k=0;
for(i=0;i<n;i++)
 {
  if(mhu[i]==0||mhu[i]==4)
    dU[i]=0;
  else
   {
    if(LimitNebals==true)
    {
    if(fabs(Rez[k])<=dUmax*abs(unom[i]))
      dU[i]=Rez[k];
    else if(Rez[k]<0)
      dU[i]=-dUmax*abs(unom[i]);
    else
      dU[i]=dUmax*abs(unom[i]);
    }
    else
     dU[i]=Rez[k];
    if(fabs(Rez[k])>eps)
     RegOk=false;
    k++;
   }
 }
k=0;
for(i=0;i<n;i++)
 {
  if(mhu[i]==0)
    dD[i]=0;
  else
   {
    if(LimitNebals==true)
    {
    if(fabs(Rez[k+n-nb-ng])<=dfmax*M_PI/180.)
      dD[i]=Rez[k+n-nb-ng];
    else if(Rez[k+n-nb-ng]<0)
      dD[i]=-dfmax*M_PI/180.;
    else
      dD[i]=dfmax*M_PI/180.;
    }
    else
      dD[i]=Rez[k+n-nb-ng];
    if(fabs(Rez[k+n-nb-ng])>eps)
     RegOk=false;
    k++;
   }
 }
k=0;
for(i=0;i<n;i++)
 {
  urab[i]=polar(abs(urab[i])-dU[i],arg(urab[i])-dD[i]);
 }

if(Raschtype!=0)
{
if(NoKn!=-1)
{
for(i=0;i<nJac/2;i++)
 {
  if(NoKn!=i)
   {
    if(LimitNebals==true&&R[i]!=0)
    {
    if(fabs((R[i]-Rez[i+nJac/2])/R[i])<=dVri)
      R[i]-=Rez[i+nJac/2];
    else if(Rez[i+nJac/2]<0)
//    R[i]-=-dVri;
      R[i]-=-dVri*fabs(R[i]);
    else
//      R[i]-=dVri;
      R[i]-=dVri*fabs(R[i]);
    }
    else
      R[i]-=Rez[i+nJac/2];
   }
  else if(Raschtype==1)
    {
    if(LimitNebals==true)
     {
     if(fabs(Rez[i+nJac/2])<dKnmax)
       Kn-=Rez[i+nJac/2];
     else if (Rez[i+nJac/2]<0)
       Kn-=-dKnmax;
     else
       Kn-=dKnmax;
     }
    else
       Kn-=Rez[i+nJac/2];
/*
     if(Kn<0.1&&Raschtype==1)
      Kn=0.1;
     if(Kn>2.&&Raschtype==1)
      Kn=2.;
*/
    }
  else if(Raschtype==2)
    {
    if(LimitNebals==true&&T!=0)
     {
     if(fabs((T-Rez[i+nJac/2])/T)<dTmax)
       T-=Rez[i+nJac/2];
     else if (Rez[i+nJac/2]<0)
//       T-=-dTmax;
       T-=-dTmax*fabs(T);
     else
//       T-=dTmax;
       T-=dTmax*fabs(T);
     }
    else
       T-=Rez[i+nJac/2];
    }
  if (fabs(Rez[i+nJac/2])>eps)
    RegOk=false;
 }
}
}

for(i=0;i<nJac;i++)
 if(fabs(Nebal[i])>eps)
  RegOk=false;
delete dU;
delete dD;

/*
double NormV=0;

for(i=0;i<nJac;i++)
 {
  NormV+=Nebal[i]*Nebal[i];
 }
NormV=sqrt(NormV);
Form1->Memo1->Lines->Add(FloatToStrF(NormV,ffFixed,7,4));
*/
return RegOk;
 }
//---------------------------------------------------------------------------
void Reg::JacobiShow(TMemo *Memo)
{
int i,j;
AnsiString tmpStr;
Memo->Lines->Add("*******Jacobi******************************");
for(i=0;i<nJac;i++)
  {
//   tmpStr=AnsiString(i+1)+" ";
   tmpStr="";
   for(j=0;j<nJac+1;j++)
     tmpStr+=FloatToStrF(Jac[i][j],ffFixed,16,16)+";";
   Memo->Lines->Add(tmpStr);
  }
Memo->Lines->Add("*************************************");
}
//---------------------------------------------------------------------------
void Reg::YsShow(TMemo *Memo)
{
int i,j;
AnsiString tmpStr;
Memo->Lines->Add("*************************************");
Memo->Lines->Add("Real:");
for(i=0;i<n;i++)
  {
   tmpStr=AnsiString(i)+" ";
   for(j=0;j<n;j++)
     tmpStr+=FloatToStrF(real(ys[i][j]),ffFixed,7,7)+"; ";
   Memo->Lines->Add(tmpStr);
  }
Memo->Lines->Add("Imag:");
for(i=0;i<n;i++)
  {
   tmpStr=AnsiString(i)+" ";
   for(j=0;j<n;j++)
     tmpStr+=FloatToStrF(imag(ys[i][j]),ffFixed,7,7)+"; ";
   Memo->Lines->Add(tmpStr);
  }
Memo->Lines->Add("Z-Real:");
for(i=0;i<n;i++)
  {
   tmpStr=AnsiString(i)+" ";
   for(j=0;j<n;j++)
     if(ys[i][j]!=complex<double>(0,0))
       tmpStr+=FloatToStrF(real(complex<double>(1,0)/ys[i][j]),ffFixed,7,7)+"; ";
     else
       tmpStr+="$";
   Memo->Lines->Add(tmpStr);
  }
Memo->Lines->Add("Z-Imag:");
for(i=0;i<n;i++)
  {
   tmpStr=AnsiString(i)+" ";
   for(j=0;j<n;j++)
     if(ys[i][j]!=complex<double>(0,0))
       tmpStr+=FloatToStrF(imag(complex<double>(1,0)/ys[i][j]),ffFixed,7,7)+"; ";
     else
       tmpStr+="$";
   Memo->Lines->Add(tmpStr);
  }
Memo->Lines->Add ("Коэфф. Тр-ии");
Memo->Lines->Add ("i    j    cktr    cktm");
for(i=0;i<n;i++)
  for(j=0;j<n;j++)
//    if(cktr[i][j]!=complex<double>(1,0)&&cktr[i][j]!=complex<double>(0,0))
    if(y[i][j]!=complex<double>(0,0)&&((msw[i][j]==1)||(msw[i][j]==2&&i<j)))

//       Memo->Lines->Add (AnsiString (i+1)+";"+AnsiString(j+1)+";"+AnsiString(real(cktr[i][j]))+";"+AnsiString(imag(cktr[i][j]))+";"+AnsiString(real(1./y[i][j]))+";"+AnsiString(imag(1./y[i][j]))+";");
       Memo->Lines->Add (AnsiString (Nu[i])+";"+AnsiString(Nu[j])+";"+AnsiString(real(cktr[i][j]))+";"+AnsiString(imag(cktr[i][j]))+";"+AnsiString(real(1./y[i][j]))+";"+AnsiString(imag(1./y[i][j]))+";");
Memo->Lines->Add("*************************************");
}
//---------------------------------------------------------------------------
bool Reg::ReorganizeJacobi()
{
double tmp;
bool ok=true;
for(int i=0;i<nJac;i++)
  if(Jac[i][i]==0)
   {
    ok=false;
    for(int k=i+1;k<nJac;k++)
     {
      if(Jac[k][i]!=0)
       {
        for(int j=0;j<nJac+1;j++)
          {

           tmp=Jac[i][j];
           Jac[i][j]=Jac[k][j];
           Jac[k][j]=tmp;
          }
        ok=true;
        break;
       }
     }
    if(ok==false)
      for(int k=0;k<nJac;k++)
        if(Jac[k][i]!=0&&Jac[i][k]!=0)
          {
           for(int j=0;j<nJac+1;j++)
             {
              tmp=Jac[i][j];
              Jac[i][j]=Jac[k][j];
              Jac[k][j]=tmp;
             }
           ok=true;
           break;
          }
    if(ok==false)
      for(int k=0;k<nJac;k++)
        if(Jac[k][i]!=0)
          {
           for(int j=0;j<nJac+1;j++)
             {
              Jac[i][j]+=Jac[k][j];
             }
           ok=true;
           break;
          }

   }
return ok;
}
//---------------------------------------------------------------------------
bool Reg::ReorganizeMatr(complex<double> **Matr,int razm)
{
complex<double> tmp,zero=complex<double>(0,0);
bool ok=true;
for(int i=0;i<razm;i++)
  if(Matr[i][i]==zero)
   {
    ok=false;
    for(int k=i+1;k<razm;k++)
     {
      if(Matr[k][i]!=zero)
       {
        for(int j=0;j<razm+1;j++)
          {
           tmp=Matr[i][j];
           Matr[i][j]=Matr[k][j];
           Matr[k][j]=tmp;
          }
        ok=true;
        break;
       }
     }
    if(ok==false)
      for(int k=0;k<razm;k++)
        if(Matr[k][i]!=zero&&Matr[i][k]!=zero)
          {
           for(int j=0;j<razm+1;j++)
             {
              tmp=Matr[i][j];
              Matr[i][j]=Matr[k][j];
              Matr[k][j]=tmp;
             }
           ok=true;
           break;
          }
   }
return ok;

}

//---------------------------------------------------------------------------
bool Reg::ReorganizeMatr(double **Matr,int razm)
{
double tmp;
bool ok=true;
for(int i=0;i<razm;i++)
  if(Matr[i][i]==0)
   {
    ok=false;
    for(int k=i+1;k<razm;k++)
     {
      if(Matr[k][i]!=0)
       {
        for(int j=0;j<razm+1;j++)
          {

           tmp=Matr[i][j];
           Matr[i][j]=Matr[k][j];
           Matr[k][j]=tmp;
          }
        ok=true;
        break;
       }
     }
    if(ok==false)
      for(int k=0;k<razm;k++)
        if(Matr[k][i]!=0&&Matr[i][k]!=0)
          {
           for(int j=0;j<razm+1;j++)
             {
              tmp=Matr[i][j];
              Matr[i][j]=Matr[k][j];
              Matr[k][j]=tmp;
             }
           ok=true;
           break;
          }
    if(ok==false)
      for(int k=0;k<razm;k++)
        if(Matr[k][i]!=0)
          {
           for(int j=0;j<razm+1;j++)
             {
              Matr[i][j]+=Matr[k][j];
             }
           ok=true;
           break;
          }

   }
return ok;
}
//---------------------------------------------------------------------------


void Reg::ReCreateYs()
{
int i,j;
double X0,R0,G0,B0;
complex<double> Un;
for(i=0;i<kt;i++)
  {
  if(Form1->TrfList[i].ktr>1)
    Un=urab[Form1->TrfList[i].NodejNumber-1];
  else
    Un=urab[Form1->TrfList[i].NodeiNumber-1];
  X0=(Form1->TrfList[i].uk*abs(Un)*abs(Un)/(100.*Form1->TrfList[i].sn))/Form1->TrfList[i].n;
  R0=((Form1->TrfList[i].pk*abs(Un)*abs(Un)/1000.)/(Form1->TrfList[i].sn*Form1->TrfList[i].sn))/Form1->TrfList[i].n;
  G0=(Form1->TrfList[i].px/(abs(Un)*abs(Un)*1000.))*Form1->TrfList[i].n;
  B0=(Form1->TrfList[i].ix*Form1->TrfList[i].sn/(100.*abs(Un)*abs(Un)))*Form1->TrfList[i].n;
  if(Form1->TrfList[i].ktr>1)
    {
     y[Form1->TrfList[i].NodeiNumber-1][Form1->TrfList[i].NodejNumber-1]=1./complex<double>(R0,X0);
     y[Form1->TrfList[i].NodejNumber-1][Form1->TrfList[i].NodeiNumber-1]=complex<double>(G0,B0);
    }
  else
    {
     y[Form1->TrfList[i].NodejNumber-1][Form1->TrfList[i].NodeiNumber-1]=1./complex<double>(R0,X0);
     y[Form1->TrfList[i].NodeiNumber-1][Form1->TrfList[i].NodejNumber-1]=complex<double>(G0,B0);
    }

  }
}
//---------------------------------------------------------------------------

double Reg::Norma(double **Matr,int razm,int skip)
{
int i,j;
double maxmod=0,summod=0;
for(i=0;i<razm;i++)
  {
   if(i==skip)
     continue;
   summod=0;
   for(j=0;j<razm;j++)
     {
      if(j==skip)
        continue;
      summod+=fabs(Matr[i][j]);
     }
   if(summod>maxmod)
     maxmod=summod;
  }
return maxmod;
}
//---------------------------------------------------------------------------
double Reg::Norma(complex<double> **Matr,int razm,int skip)
{
int i,j;
double maxmod=0,summod=0;
for(i=0;i<razm;i++)
  {
   if(i==skip)
     continue;
   summod=0;
   for(j=0;j<razm;j++)
     {
      if(j==skip)
        continue;
      summod+=fabs(real(Matr[i][j]));
     }
   if(summod>maxmod)
     maxmod=summod;
   summod=0;
   for(j=0;j<razm;j++)
     {
      if(j==skip)
        continue;
      summod+=fabs(imag(Matr[i][j]));
     }
   if(summod>maxmod)
     maxmod=summod;

  }
return maxmod;
}
//---------------------------------------------------------------------------
void Reg::CheckRes(TMemo *Memo)
{
int i,j;
complex<double> tmpCompl;
for(i=0;i<n;i++)
  {
   tmpCompl=complex<double>(0,0);
   for(j=0;j<n;j++)
     tmpCompl+=ys[i][j]*urab[j];//*conj(urab[i]);
   Memo->Lines->Add(AnsiString(real(tmpCompl))+"; "+AnsiString(imag(tmpCompl)));
  }
}
//---------------------------------------------------------------------------

Reg::Reg(int nuzl)
{
int i;
RegimState=false;
mhu=new int [nuzl+1];
msw=new int *[nuzl+1];
a0=new double [nuzl+1];
a1=new double [nuzl+1];
a2=new double [nuzl+1];
b0=new double [nuzl+1];
b1=new double [nuzl+1];
b2=new double [nuzl+1];
sng=new complex<double>[nuzl+1];
unom=new complex<double>[nuzl+1];
urab=new complex<double>[nuzl+1];
sng_r=new complex<double>[nuzl+1];
Pg=new double[nuzl+1];
Qg=new double[nuzl+1];
dPg=new double[nuzl+1];
dPn=new double[nuzl+1];
dQn=new double[nuzl+1];
irab=new complex<double>[nuzl+1];
cktr=new complex<double>*[nuzl+1];
y=new complex<double>*[nuzl+1];
ys=new complex<double>*[nuzl+1];
z=new complex<double>*[nuzl+1];
Params=new struct Prms[nuzl*nuzl+1];
Nu=new int [nuzl+1];
for(i=0;i<nuzl+1;i++)
  {
   msw[i]=new int[nuzl+1];
   cktr[i]=new complex<double>[nuzl+1];
   y[i]=new complex<double>[nuzl+1];
   ys[i]=new complex<double>[nuzl+1];
   z[i]=new complex<double>[nuzl+1];
   Nu[i]=-1;
  }
QminQmax=new double*[2];
for(i=0;i<2;i++)
  QminQmax[i]=new double[nuzl+1];
Zu=new int[nuzl+1];
Zf=new int[nuzl+1];

R=new double [nuzl*2];
Jac=new double *[nuzl*4];
Nebal=new double[nuzl*4];
for(i=0;i<nuzl*4;i++)
 Jac[i]=new double[nuzl*4];

/*
mhu=new int [MaxNodesNumber+5];
msw=new int *[MaxNodesNumber+5];
a0=new double [MaxNodesNumber+5];
a1=new double [MaxNodesNumber+5];
a2=new double [MaxNodesNumber+5];
b0=new double [MaxNodesNumber+5];
b1=new double [MaxNodesNumber+5];
b2=new double [MaxNodesNumber+5];
sng=new complex<double>[MaxNodesNumber+5];
unom=new complex<double>[MaxNodesNumber+5];
urab=new complex<double>[MaxNodesNumber+5];
sng_r=new complex<double>[MaxNodesNumber+5];
Pg=new double[MaxNodesNumber+5];
Qg=new double[MaxNodesNumber+5];
dPg=new double[MaxNodesNumber+5];
dPn=new double[MaxNodesNumber+5];
dQn=new double[MaxNodesNumber+5];
irab=new complex<double>[MaxNodesNumber+5];
cktr=new complex<double>*[MaxNodesNumber+5];
y=new complex<double>*[MaxNodesNumber+5];
ys=new complex<double>*[MaxNodesNumber+5];
z=new complex<double>*[MaxNodesNumber+5];
Params=new struct Prms[MaxNodesNumber*MaxNodesNumber+5];
Nu=new int [MaxNodesNumber+5];
for(i=0;i<MaxNodesNumber+5;i++)
  {
   msw[i]=new int[MaxNodesNumber+5];
   cktr[i]=new complex<double>[MaxNodesNumber+5];
   y[i]=new complex<double>[MaxNodesNumber+5];
   ys[i]=new complex<double>[MaxNodesNumber+5];
   z[i]=new complex<double>[MaxNodesNumber+5];
   Nu[i]=-1;
  }
QminQmax=new double*[2];
for(i=0;i<2;i++)
  QminQmax[i]=new double[MaxNodesNumber+5];
Zu=new int[MaxNodesNumber+5];
Zf=new int[MaxNodesNumber+5];

R=new double [MaxNodesNumber*2+5];
Jac=new double *[MaxNodesNumber*4+5];
Nebal=new double[MaxNodesNumber*4+5];
for(i=0;i<MaxNodesNumber*4+5;i++)
 Jac[i]=new double[MaxNodesNumber*4+5];
*/
n=nuzl;

oldnJac=-1;
nJac=0;
}
//---------------------------------------------------------------------------
void Reg::Initialize(int Nodes)
{
/*
int i;
delete[] mhu;
delete[] a0;
delete[] a1;
delete[] a2;
delete[] b0;
delete[] b1;
delete[] b2;
delete[] sng;
delete[] unom;
delete[] urab;
delete[] sng_r;
delete[] Pg;
delete[] Qg;
delete[] dPg;
delete[] dPn;
delete[] dQn;
delete[] irab;
delete[] Params;
delete[] Nu;
for(i=0;i<n+5;i++)
  {
   delete[] msw[i];
   delete[] cktr[i];
   delete[] y[i];
   delete[] ys[i];
   delete[] z[i];
  }
delete[] msw;
delete[] cktr;
delete[] y;
delete[] ys;
delete[] z;
delete[] Zu;
delete[] Zf;

for(i=0;i<2;i++)
  delete[] QminQmax[i];
delete[] QminQmax;

for(i=0;i<n*4+5;i++)
 delete[] Jac[i];
delete[] Jac;
delete[] Nebal;
delete []R;


mhu=new int [Nodes+5];
a0=new double [Nodes+5];
a1=new double [Nodes+5];
a2=new double [Nodes+5];
b0=new double [Nodes+5];
b1=new double [Nodes+5];
b2=new double [Nodes+5];
sng=new complex<double>[Nodes+5];
unom=new complex<double>[Nodes+5];
urab=new complex<double>[Nodes+5];
sng_r=new complex<double>[Nodes+5];
Pg=new double[Nodes+5];
Qg=new double[Nodes+5];
dPg=new double[Nodes+5];
dPn=new double[Nodes+5];
dQn=new double[Nodes+5];
irab=new complex<double>[Nodes+5];
msw=new int *[Nodes+5];
cktr=new complex<double>*[Nodes+5];
y=new complex<double>*[Nodes+5];
ys=new complex<double>*[Nodes+5];
z=new complex<double>*[Nodes+5];
Params=new struct Prms[Nodes*Nodes+5];
Nu=new int [Nodes+5];
for(i=0;i<Nodes+5;i++)
  {
   msw[i]=new int[Nodes+5];
   cktr[i]=new complex<double>[Nodes+5];
   y[i]=new complex<double>[Nodes+5];
   ys[i]=new complex<double>[Nodes+5];
   z[i]=new complex<double>[Nodes+5];
   Nu[i]=-1;
  }
QminQmax=new double*[2];
for(i=0;i<2;i++)
  QminQmax[i]=new double[Nodes+5];
Zu=new int[Nodes+5];
Zf=new int[Nodes+5];

R=new double [Nodes*2+5];
Jac=new double *[Nodes*4+5];
Nebal=new double[Nodes*4+5];
for(i=0;i<Nodes*4+5;i++)
 Jac[i]=new double[Nodes*4+5];
n=Nodes;
*/
//---------------------------------------------------------------------------
}
Reg::~Reg()
{
int i;
delete []mhu;
delete []a0;delete []a1;delete []a2;
delete []b0;delete []b1;delete []b2;
delete []sng;delete []unom;delete []urab;
delete []sng_r;
delete []Pg;
delete []Qg;
delete []dPg;
delete []dPn;
delete []dQn;
delete []irab;
delete []Params;
delete []Nu;
for(i=0;i<n+1;i++)
  {
   delete []msw[i];
   delete []cktr[i];
   delete []y[i];
   delete []ys[i];
   delete []z[i];
  }
delete []msw;
delete []cktr;
delete []y;
delete []ys;
delete []z;
delete []Zu;
delete []Zf;

for(i=0;i<2;i++)
  delete QminQmax[i];
delete []QminQmax;

delete []R;

for(i=0;i<n*4;i++)
 delete []Jac[i];
delete []Jac;
delete []Nebal;
}

//---------------------------------------------------------------------------
void Reg::InitializeJac(int Nodes)
{
/*
int i;
for(i=0;i<oldnJac+2;i++)
 delete Jac[i];
delete []Jac;
delete []Nebal;
//delete []R;

Jac=new double *[Nodes+2];
Nebal=new double[Nodes+2];
for(i=0;i<Nodes+2;i++)
 Jac[i]=new double[Nodes+2];
oldnJac=Nodes;
*/
}
//---------------------------------------------------------------------------



int Reg::YsCreate()
{
int i,j;
int krm=0;
//complex<double>tmpCompl;
for(i=0;i<n;i++)
for(j=0;j<n;j++)
ys[i][j]=complex<double>(0,0);

for(i=0;i<n;i++)
{
// if(mhu[i]==4) mhu[i]=1;
// if(mhu[i]==3||mhu[i]==0) continue;
 if(mhu[i]!=3&&mhu[i]!=0)
   krm++;
 for(j=0;j<n;j++)
   {

    if(i==j)continue;
 //   tmpCompl=y[i][j]*cktr[i][j]*cktr[i][j];
    if(msw[i][j]==1)//для повышающего трансформатора
      ys[i][i]-=y[i][j]*cktr[i][j]*cktr[i][j];
      //ys[i][i]-=y[i][j]*cktr[i][j]*cktr[i][j]+y[j][i]*cktr[i][j]*cktr[i][j];
    else if(msw[i][j]==-1||msw[i][j]==2)//для понижающего трансформатора и линий
      ys[i][i]-=y[i][j]+y[j][i];
   }
 for(j=i+1;j<n;j++)
   {
//    if(mhu[j]==3||mhu[j]==0) continue;
    if(msw[i][j]==1)
      ys[i][j]=ys[j][i]=y[i][j]*cktr[i][j];
    else if(msw[i][j]==-1)
      ys[i][j]=ys[j][i]=y[j][i]/cktr[i][j];
    else if(msw[i][j]==2)
      ys[i][j]=ys[j][i]=y[i][j];

   }
//если есть, то учитываем шунтирующие реакторы
ys[i][i]-=y[i][i];
}
//YsShow(Form1->Memo1);
return krm;
}
//---------------------------------------------------------------------------
void Reg::ZCreate(int *ki,int *kj,int &ii,int &jj)
{
int i,j,k;
complex<double> tmpCompl;
ii=0;
for(i=0;i<n;i++)
{
 if(mhu[i]==0||mhu[i]==3) continue;
 jj=0;

 for(j=0;j<n;j++)
   {
    if(mhu[j]==0||mhu[j]==3) continue;

    z[ii][jj]=ys[i][j];
    ki[ii]=i;
    kj[jj]=j;
    jj++;
   }
 ii++;
}
for(i=0;i<ii;i++)
 for(j=0;j<jj;j++)
   {
    if(i==j) z[i][j+jj]=complex<double>(1,0);
    else     z[i][j+jj]=complex<double>(0,0);
   }
for(i=0;i<ii;i++)
{
 tmpCompl=z[i][i];
 for(j=0;j<jj*2;j++)
   {

    if(tmpCompl!=complex<double>(0,0))
      z[i][j]=z[i][j]/tmpCompl;
    else ShowMessage("Матрица проводимостей не верна!");

   }

 for(k=0;k<ii;k++)
   {
    tmpCompl=z[k][i];
    for(j=0;j<jj*2;j++)
     {
      if(k!=i)
         {
          z[k][j]=z[k][j]-z[i][j]*tmpCompl;
         }
     }
   }
}
for(i=0;i<ii;i++)
  for(j=0;j<jj;j++)
   z[i][j]=z[i][j+jj];
}

//---------------------------------------------------------------------------
//void Reg::Zsimple()
void Reg::ObrMatr(complex<double> **Matr,complex<double> **ObrMatr,int razm)
{
int i,j,k;
complex<double> tmpCompl;
for(i=0;i<razm;i++)
{
 for(j=0;j<razm;j++)
   {
    ObrMatr[i][j]=Matr[i][j];
   }
}
for(i=0;i<razm;i++)
 for(j=0;j<razm;j++)
   {
    if(i==j) ObrMatr[i][j+razm]=complex<double>(1,0);
    else     ObrMatr[i][j+razm]=complex<double>(0,0);
   }
for(i=0;i<razm;i++)
{
 tmpCompl=ObrMatr[i][i];
 for(j=0;j<razm*2;j++)
   {

    if(tmpCompl!=complex<double>(0,0))
      ObrMatr[i][j]=ObrMatr[i][j]/tmpCompl;
    else ShowMessage("Матрица проводимостей не верна!");

   }

 for(k=0;k<razm;k++)
   {
    tmpCompl=ObrMatr[k][i];
    for(j=0;j<razm*2;j++)
     {
      if(k!=i)
         {
          ObrMatr[k][j]=ObrMatr[k][j]-ObrMatr[i][j]*tmpCompl;
         }
     }
   }
}
for(i=0;i<razm;i++)
  for(j=0;j<razm;j++)
   ObrMatr[i][j]=ObrMatr[i][j+razm];

}
//---------------------------------------------------------------------------
void Reg::ObrMatr(double **Matr,double **ObrMatr,int razm)
{
int i,j,k;
double tmpCompl;
for(i=0;i<razm;i++)
{
 for(j=0;j<razm;j++)
   {
    ObrMatr[i][j]=Matr[i][j];
   }
}
for(i=0;i<razm;i++)
 for(j=0;j<razm;j++)
   {
    if(i==j) ObrMatr[i][j+razm]=1;
    else     ObrMatr[i][j+razm]=0;
   }
for(i=0;i<razm;i++)
{
 tmpCompl=ObrMatr[i][i];
 for(j=0;j<razm*2;j++)
   {

    if(tmpCompl!=0)
      ObrMatr[i][j]=ObrMatr[i][j]/tmpCompl;
    else ShowMessage("Матрица проводимостей не верна!");

   }

 for(k=0;k<razm;k++)
   {
    tmpCompl=ObrMatr[k][i];
    for(j=0;j<razm*2;j++)
     {
      if(k!=i)
         {
          ObrMatr[k][j]=ObrMatr[k][j]-ObrMatr[i][j]*tmpCompl;
         }
     }
   }
}
for(i=0;i<razm;i++)
  for(j=0;j<razm;j++)
   ObrMatr[i][j]=ObrMatr[i][j+razm];

}
//---------------------------------------------------------------------------
double Reg::fi(int i,double y1,double y2,double f1,double f2)
{
//double tmp11;
switch (i)
{
case 1:
  {
   return y1*(cos(f1)*cos(f2)+sin(f1)*sin(f2))+y2*(cos(f1)*sin(f2)-sin(f1)*cos(f2));
  }
case 2:
  {
//   tmp11=y1*(-sin(f1)*cos(f2)+cos(f1)*sin(f2))+y2*(-sin(f1)*sin(f2)-cos(f1)*cos(f2));
   return y1*(-sin(f1)*cos(f2)+cos(f1)*sin(f2))+y2*(-sin(f1)*sin(f2)-cos(f1)*cos(f2));
  }
case 3:
  {
   return y1*(-cos(f1)*sin(f2)+sin(f1)*cos(f2))+y2*(cos(f1)*cos(f2)+sin(f1)*sin(f2));
  }
case 4:
  {
   return y1*(-cos(f1)*cos(f2)-sin(f1)*sin(f2))+y2*(-cos(f1)*sin(f2)+sin(f1)*cos(f2));
  }
case 5:
  {
   return y1*(sin(f1)*sin(f2)+cos(f1)*cos(f2))+y2*(-sin(f1)*cos(f2)+cos(f1)*sin(f2));
  }
}
}
//---------------------------------------------------------------------------
bool Reg::SobstvQri(double **matr,double *evr,double *evi,double **v,int nrazm, bool isShow,TMemo *ShowMemo)
{
int i,j;
bool isOk;
AnsiString tmpStr;
ap::real_2d_array newmatr;
ap::real_1d_array my_evr;
ap::real_1d_array my_evi;
ap::real_2d_array my_v;
newmatr.setbounds(1, nrazm, 1, nrazm);
my_evr.setbounds(1, nrazm);
my_evi.setbounds(1, nrazm);
my_v.setbounds(1, nrazm, 1, nrazm);
for (i=0;i<nrazm;i++)
 for(j=0;j<nrazm;j++)
  newmatr(i+1,j+1)=matr[i][j];
isOk=hessenbergqrieigenvaluesandvectors(newmatr,nrazm,my_evr,my_evi,my_v);
for(i=0;i<nrazm;i++)
 {
  evr[i]=my_evr(i+1);
  evi[i]=my_evi(i+1);
  for(j=0;j<nrazm;j++)
    v[i][j]=my_v(i+1,j+1);
 }
if(isShow==true)
{
tmpStr="";
ShowMemo->Lines->Add("****Evr****");
for(i=0;i<nrazm;i++)
  tmpStr+=AnsiString(evr[i])+";";
ShowMemo->Lines->Add(tmpStr);
tmpStr="";
ShowMemo->Lines->Add("****Evi****");
for(i=0;i<nrazm;i++)
  tmpStr+=AnsiString(evi[i])+";";
ShowMemo->Lines->Add(tmpStr);
ShowMemo->Lines->Add("****Vi****");
for(i=0;i<nrazm;i++)
 {
  tmpStr="";
  for(j=0;j<nrazm;j++)
   tmpStr+=AnsiString(v[j][i])+";";
  ShowMemo->Lines->Add(tmpStr);
 }
}
return isOk;
}
//---------------------------------------------------------------------------
bool Reg::QR(double **a,double *b,double *x,int n,bool z)
{
int i,j,k,l,m;
double r1,r2,r,p,t,c,s,a1,a2;
double d,d1;

for(i=n-1;i>=0;i--)
 {
  for(j=n-1;j>=0;j--)
   a[i+1][j+1]=a[i][j];
  b[i+1]=b[i];
 }

m=n;
for (i=1;i<=n;i++)
 x[i]=b[i];
if(z==true)
{
 for (k=1;k<=n;k++)
 {
  l=k+1;
  if(l>m)
    continue;
  for (j=l;j<=m;j++)
  {
   r2=a[k][j];
   if(r2==0)
    continue;
   r1=a[k][k];
   a1=fabs(r1);
   a2=fabs(r2);
   if(a1!=0)
     goto L2;
   t=1;
   a[k][k]=a[k][j];
   goto L6;
L2:
   if(a1<=a2)
     goto L3;
   p=r2/r1;
   r=sqrt(1.+p*p);
   t=p/(r+1.);
   a[k][k]=r*r1;
   goto L6;
L3:
   p=r1/r2;
   r=sqrt(1.+p*p);
   t=1./(r+fabs(p));
   if(p<0)
    t=-t;
   if(r1>=0)
     goto L5;
   else
     goto L4;
L4:
   a[k][k]=-r*a2;
     goto L6;
L5:
   a[k][k]=r*a2;
L6:
   a[k][j]=t;
   r=t*t;
   r1=1./(1.+r);
   c=(1.-r)*r1;
   s=(t+t)*r1;
   if(l>n)
     continue;
//     goto L8;
   for (i=l;i<=n;i++)
    {
     r1=a[i][k];
     r2=a[i][j];
     a[i][k]=r1*c+r2*s;
     a[i][j]=r2*c-r1*s;
    }
  }
 }
}
x[1]=x[1]/a[1][1];
if(n==1)
 goto L13;
for(k=2;k<=n;k++)
{
 l=k-1;
 d=0;
 for(j=1;j<=l;j++)
  {
   d1=a[k][j];
   d=d+d1*x[j];
  }
 d1=x[k];
 d=d1-d;
 r2=d;
 x[k]=r2/a[k][k];
}
L13:
k=n;
goto L16;
L14:
j=m;
L15:
t=a[k][j];
r=t*t;
r1=1./(1.+r);
c=(1.-r)*r1;
s=(t+t)*r1;
r1=x[k];
r2=x[j];
x[k]=c*r1-s*r2;
x[j]=s*r1+c*r2;
j=j-1;
if(j>k)
 goto L15;
L16:
k=k-1;
if(k>0)
  goto L14;

for(i=0;i<n;i++)
 {
  x[i]=x[i+1];
 }


return true;
}
//---------------------------------------------------------------------------
void Reg::NjNbNg()
{
nJac=2*n;
nb=ng=0;
for(int i=0;i<n;i++)
 {
 if(mhu[i]==0)
   {
    nJac-=2;
    nb++;
   }

 else if (mhu[i]==4)
   {
    nJac-=1;
    ng++;
   }
 }

}
//---------------------------------------------------------------------------
void Reg::SHN(bool stattrue)
{
int i;
double p,q;
if(stattrue==true)
{
for(i=0;i<n;i++)
// if(mhu[i]==1)
 {
  p=real(sng[i])*(a0[i]+a1[i]*(abs(urab[i])/abs(unom[i]))+a2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
//  if(mhu[i]==4)
   p-=Pg[i];
  q=imag(sng[i])*(b0[i]+b1[i]*(abs(urab[i])/abs(unom[i]))+b2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
  q+=Qg[i];
  sng_r[i]=complex<double>(p,q);
 }
}
else
for(i=0;i<n;i++)
 sng_r[i]=sng[i]+complex<double>(-Pg[i],Qg[i]);
}

//---------------------------------------------------------------------------
void Reg::SHN_Hard(bool stattrue,double T)
{
int i;
double p,q;
if(stattrue==true)
{
for(i=0;i<n;i++)
// if(mhu[i]==1)
 {
  p=real(sng[i])*(a0[i]+a1[i]*(abs(urab[i])/abs(unom[i]))+a2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
  p-=Pg[i];
  if(fabs(dPn[i])>0)
   {
    if(real(sng[i])>=0)
     p+=T*dPn[i];
    else
     p-=T*dPn[i];
   }
  if(fabs(dPg[i])>0)
   p-=T*dPg[i];
 //  if(mhu[i]==4)


  q=imag(sng[i])*(b0[i]+b1[i]*(abs(urab[i])/abs(unom[i]))+b2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
  q+=Qg[i];
  if(fabs(dQn[i])>0)
   {
    if(imag(sng[i])<=0)
     q-=T*dQn[i];
    else
     q+=T*dQn[i];
   }
  sng_r[i]=complex<double>(p,q);
 }
}
else
for(i=0;i<n;i++)
 sng_r[i]=sng[i]+complex<double>(-Pg[i],Qg[i]);
}

//---------------------------------------------------------------------------

bool Reg::Limits(bool limtype,bool useT,double T)
{
int i,j;
double p,q;
bool Limit=false;
//if(limtype==0)//
 for(i=0;i<n;i++)
  if(mhu_n[i]==4)
  {
//   p=real(sng[i])*(a0[i]+a1[i]*(abs(urab[i])/abs(unom[i]))+a2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
  p=real(sng[i])*(a0[i]+a1[i]*(abs(urab[i])/abs(unom[i]))+a2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
  p-=Pg[i];
  if(useT==true)
  {
  if(fabs(dPn[i])>0)
   {
    if(real(sng[i])>=0)
     p+=T*dPn[i];
    else
     p-=T*dPn[i];
   }
  if(fabs(dPg[i])>0)
   p-=T*dPg[i];
  }
   q=-imag(sng[i])*(b0[i]+b1[i]*(abs(urab[i])/abs(unom[i]))+b2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
//   p-=Pg[i];

  if(useT==true)
  {
  if(fabs(dQn[i])>0)
   {
    if(imag(sng[i])<=0)
     q+=T*dQn[i];
    else
     q-=T*dQn[i];
   }
   }
   q+=imag(ys[i][i])*abs(urab[i])*abs(urab[i]);

  for(j=0;j<n;j++)
   if(j!=i)
     q+=abs(urab[i])*abs(urab[j])*fi(1,imag(ys[i][j]),real(ys[i][j]),arg(urab[i]),arg(urab[j]));
//  q+=
   if(q<QminQmax[0][i]&&limtype==1)//ограничения ведущие к избытку реакт. мощности
    {
     q=QminQmax[0][i];
     Qg[i]=q;
     q+=imag(sng[i])*(b0[i]+b1[i]*(abs(urab[i])/abs(unom[i]))+b2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
     if(useT==true)
     if(fabs(dQn[i])>0)
      {
       if(imag(sng[i])<=0)
        q-=T*dQn[i];
       else
        q+=T*dQn[i];
      }
     Limit=true;
     sng_r[i]=complex<double>(p,q);
     mhu[i]=1;
    }
   else if(q>QminQmax[1][i]&&limtype==0)//ограничения ведущие к недостатку реакт. мощности
    {
     q=QminQmax[1][i];
     Qg[i]=q;
     q+=imag(sng[i])*(b0[i]+b1[i]*(abs(urab[i])/abs(unom[i]))+b2[i]*(abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
     if(useT==true)
     if(fabs(dQn[i])>0)
      {
       if(imag(sng[i])<=0)
        q-=T*dQn[i];
       else
        q+=T*dQn[i];
      }
     Limit=true;
     sng_r[i]=complex<double>(p,q);
     mhu[i]=1;
    }
  }
return Limit;
}

//---------------------------------------------------------------------------
void Reg::RaschR0(double **Matr,int nMatr,double *R)
{
double **v=new double *[nMatr+1];
for(int i=0;i<nMatr+1;i++)
 v[i]=new double [nMatr+1];
double *evr=new double [nMatr+1];
double *evi=new double [nMatr+1];
TMemo *tmpMemo;
SobstvQri(Matr,evr,evi,v,nMatr,false,tmpMemo);
//delete tmpMemo;
double minev=1000000000000000000000;
int minevNo;
bool flag=false;
Form1->Memo1->Lines->Add("SobstvZn");
for(int i=0;i<nMatr;i++)
 if(sqrt(evr[i]*evr[i]+evi[i]*evi[i])<minev)
  {
   minevNo=i;

   minev=sqrt(evr[i]*evr[i]+evi[i]*evi[i]);
   if (evi[i]<0)
    {minevNo--;flag=true;}
   else
     flag=false;
  }
for(int i=0;i<nMatr;i++)
  Form1->Memo1->Lines->Add(AnsiString(evr[i]));
Form1->Memo1->Lines->Add("***");
for(int i=0;i<nMatr;i++)
  Form1->Memo1->Lines->Add(AnsiString(evi[i]));

//if(evi[minevNo]!=0)
 double NormV=0;
 Form1->Memo1->Lines->Add("VectorMinSobstv");
 for(int i=0;i<nMatr;i++)
  {
   if(evi[minevNo]!=0)
     R[i]=sqrt(v[i][minevNo]*v[i][minevNo]+v[i][minevNo+1]*v[i][minevNo+1]);
   else
     R[i]=v[i][minevNo];

   NormV+=R[i]*R[i];
   if(evi[minevNo]!=0)
     R[i]*=Sign(v[i][minevNo])*Sign(v[i][minevNo+1]);
   Form1->Memo1->Lines->Add(AnsiString(R[i]));
  }
 NormV=sqrt(NormV);
 for(int i=0;i<nMatr;i++)
  {
   R[i]/=NormV;
//   Form1->Memo1->Lines->Add(AnsiString(R[i]));
  }
delete []evr;
delete []evi;
for(int i=0;i<nMatr+1;i++)
 delete v[i];
delete []v;
}
//---------------------------------------------------------------------------
void Reg::RaschSv(double **Matr,int nMatr,double *Rr,double *Ri)
{
double **v=new double *[nMatr+1];
for(int i=0;i<nMatr+1;i++)
 v[i]=new double [nMatr+1];
double *evr=new double [nMatr+1];
double *evi=new double [nMatr+1];
TMemo *tmpMemo;
SobstvQri(Matr,evr,evi,v,nMatr,false,tmpMemo);
//delete tmpMemo;
double minev=1000000000000000000000;
int minevNo;
for(int i=0;i<nMatr;i++)
 if(sqrt(evr[i]*evr[i]+evi[i]*evi[i])<minev)
  {
   minevNo=i;
   minev=sqrt(evr[i]*evr[i]+evi[i]*evi[i]);
//   if (evi[i]<0)
//    minevNo--;
  }
if(evi[minevNo]>0)
 for(int i=0;i<nMatr;i++)
  {
   Rr[i]=v[i][minevNo];
   Ri[i]=v[i][minevNo+1];
  }
else if(evi[minevNo]<0)
 for(int i=0;i<nMatr;i++)
  {
   Rr[i]=v[i][minevNo];
   Ri[i]=v[i][minevNo-1];
  }
else
 for(int i=0;i<nMatr;i++)
  {
   Rr[i]=v[i][minevNo];
   Ri[i]=0;
  }
delete []evr;
delete []evi;
for(int i=0;i<nMatr+1;i++)
 delete v[i];
delete []v;
}


//---------------------------------------------------------------------------
void Reg::EkvivalentToStar()
{
int i,j,k,l,ubazNo;
double p,q;
complex<double> tmpCompl,dU;

for(i=0;i<Rasch->n;i++)
  if(Rasch->mhu[i]==0)
    {
     ubazNo=i;
     break;
    }

for(i=0;i<n;i++)
  for(j=0;j<n;j++)
    {
     ys[i][j]=complex<double>(0,0);
     y[i][j]=complex<double>(0,0);
     msw[i][j]=0;
    }

k=ubazNo;

l=0;
j=0;
SHN(true);
for(i=0;i<n;i++)
  if(mhu[i]!=0&&sng_r[i]==complex<double>(0,0))
    j++;
for(i=0;i<n;i++)
  {
  if(mhu[i]!=0&&sng_r[i]!=complex<double>(0,0))
     {

//      dU=complex<double>(real(urab[i]*0.05),fabs(imag(urab[i]*0.05)));


      if(mhu[i]==4)
       dU=complex<double>(min(abs(urab[i]),abs(urab[k]))*0.05,min(abs(urab[i]),abs(urab[k]))*0.05);
      else
       dU=complex<double>(1,0.1);      
      p=real(sng_r[i]);
      q=imag(sng_r[i]);

//      p=real(sng_r[i])*(a0[i]+(a1[i]*abs(urab[i])/abs(unom[i]))+(a2[i]*abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));
//      q=imag(sng_r[i])*(b0[i]+(b1[i]*abs(urab[i])/abs(unom[i]))+(b2[i]*abs(urab[i])*abs(urab[i])/(abs(unom[i])*abs(unom[i]))));

      if(abs(urab[i])>abs(urab[k]))
        {
            if(p>0)
             {
              cktr[k][i]=(urab[i]+dU)/urab[k];
//              cktr[k][i]=abs(urab[i]+dU)/abs(urab[k]);
              tmpCompl=(cktr[k][i]*urab[k]-urab[i]);
              tmpCompl=polar(real(tmpCompl),imag(tmpCompl));
              tmpCompl=complex<double>(p,q)/conj(urab[i]);
              tmpCompl=(complex<double>(p,q)/conj(urab[i]))/(cktr[k][i]*urab[k]-urab[i]);

             }
            else if((abs(urab[i])-abs(urab[k]))/abs(urab[i])>0.4)
             {
              dU=dU*urab[k]/urab[i];
              cktr[k][i]=urab[i]/(urab[k]+dU);
//              cktr[k][i]=abs(urab[i])/abs(urab[k]+dU);
              tmpCompl=(urab[i]-urab[k]*cktr[k][i]);
              tmpCompl=-complex<double>(p,q)/conj(urab[i]);
              tmpCompl=(-complex<double>(p,q)/conj(urab[i]))/(urab[i]-urab[k]*cktr[k][i]);
             }

            else
             {
              cktr[k][i]=complex<double>(1,0);
              tmpCompl=-complex<double>(p,q)/conj(urab[i]);
              tmpCompl=(-complex<double>(p,q)/conj(urab[i]))/(urab[i]-urab[k]);

             }

        }
      else if((abs(urab[k])-abs(urab[i]))/abs(urab[k])>0.4)
        {
            trf1:
            if(p>0)
             {
              dU=dU*urab[i]/urab[k];
              cktr[k][i]=(urab[i]+dU)/(urab[k]);
//              cktr[k][i]=abs(urab[i]+dU)/abs(urab[k]);
              tmpCompl=complex<double>(p,q)/conj(urab[i]);
              tmpCompl=urab[k]-urab[i]/cktr[k][i];
              tmpCompl=(cktr[k][i]*complex<double>(p,q)/conj(urab[i]))/(urab[k]-urab[i]/cktr[k][i]);

//              tmpCompl=cktr[k][i]*urab[k]-urab[i];
//              tmpCompl=(complex<double>(p,q)/conj(urab[i]))/(cktr[k][i]*urab[k]-urab[i]);
             }

            else
             {
//              dU=dU*urab[i]/urab[k];
              cktr[k][i]=(urab[i])/(urab[k]+dU);
//              cktr[k][i]=abs(urab[i]+dU)/abs(urab[k]);
              tmpCompl=complex<double>(p,q)/conj(urab[i]);
              tmpCompl=urab[k]-urab[i]/cktr[k][i];
              tmpCompl=(-cktr[k][i]*complex<double>(p,q)/conj(urab[i]))/(urab[i]/cktr[k][i]-urab[k]);
             }

        }
      else
        {
//         dU=urab[k]*0.1;
         cktr[k][i]=complex<double>(1,0);
            if(p>0)
             {
              tmpCompl=complex<double>(p,q)/conj(urab[i]);
              tmpCompl=urab[k]-urab[i];
              tmpCompl=(complex<double>(p,q)/conj(urab[i]))/(urab[k]-urab[i]);
             }
            else
             {
              tmpCompl=-complex<double>(p,q)/conj(urab[i]);
              tmpCompl=urab[k]-urab[i];
              tmpCompl=(-complex<double>(p,q)/conj(urab[i]))/(urab[i]-urab[k]);
             }
            if (real(tmpCompl)<0)
              goto trf1;
        }
      cktr[i][k]=1./cktr[k][i];
      if(abs(cktr[k][i])>1)
        {
         ys[k][k]-=tmpCompl*cktr[k][i]*cktr[k][i];
         ys[i][i]-=tmpCompl;
         ys[i][k]=ys[k][i]=tmpCompl*cktr[k][i];
        }
      else
        {
/*
         ys[k][k]-=tmpCompl*cktr[k][i]*cktr[k][i];
         ys[i][i]-=tmpCompl;
         ys[i][k]=ys[k][i]=tmpCompl*cktr[k][i];
*/


         ys[i][i]-=tmpCompl/(cktr[k][i]*cktr[k][i]);
         ys[k][k]-=tmpCompl;
         ys[i][k]=ys[k][i]=tmpCompl/cktr[k][i];

        }

      if(abs(cktr[k][i])>1)
        y[k][i]=tmpCompl;
      else if (abs(cktr[k][i])<1)
        y[i][k]=tmpCompl;
      else
        y[min(k,i)][max(k,i)]=tmpCompl;

      tmpCompl=1./ys[i][k];

//      Memo1->Lines->Add(AnsiString(k+1)+ ";"+ AnsiString(i+1)+ ";"+AnsiString(real(cktr[k][i]))+";"+AnsiString(imag(cktr[k][i]))+";"+AnsiString(real(tmpCompl))+";"+AnsiString(imag(tmpCompl)));
/*
cktr[k][i]=cktr[i][k]=1;
/*
            if(p>0)
             {
              tmpCompl=complex<double>(p,q)/conj(urab[i]);
//              tmpCompl=urab[k]-urab[i];
              tmpCompl=(complex<double>(p,q)/conj(urab[i]))/(urab[k]-urab[i]);
             }


      ys[k][k]-=tmpCompl*cktr[k][i]*cktr[k][i];
      ys[i][i]-=tmpCompl;
      ys[i][k]=ys[k][i]=tmpCompl*cktr[k][i];
      tmpCompl=1./ys[i][k];
      Memo1->Lines->Add(AnsiString(k+1)+ ";"+ AnsiString(i+1)+ ";"+AnsiString(real(cktr[k][i]))+";"+AnsiString(imag(cktr[k][i]))+";"+AnsiString(real(tmpCompl))+";"+AnsiString(imag(tmpCompl)));
*/
     sng[l]=sng[i];
     sng_r[l]=sng_r[i];
     Pg[l]=Pg[i];
     Qg[l]=Qg[i];
     dPn[l]=dPn[i];
     dPg[l]=dPg[i];
     dQn[l]=dQn[i];
     urab[l]=urab[i];
     unom[l]=unom[i];
     a0[l]=a0[i];
     a1[l]=a1[i];
     a2[l]=a2[i];
     b0[l]=b0[i];
     b1[l]=b1[i];
     b2[l]=b2[i];
     mhu[l]=mhu[i];
     Nu[l]=Nu[i];
     cktr[l][k]=cktr[i][k];
     cktr[k][l]=cktr[k][i];
     if(abs(cktr[l][k])>1)
        {
         msw[l][k]=1;
         msw[k][l]=-1;
        }
     else if(abs(cktr[l][k])<1)
        {
         msw[l][k]=-1;
         msw[k][l]=1;
        }
     else
         msw[l][k]=msw[k][l]=2;

     ys[l][l]=ys[i][i];
     ys[l][k]=ys[k][l]=ys[i][k];
     y[l][k]=y[i][k];
     y[k][l]=y[k][i];

//     y[min(l,k)][max(l,k)]=y[min(i,k)][max(i,k)];
     l++;
     }
  else if(mhu[i]==0)
    {
     l++;
    }
  }
//Memo1->Lines->Add("*******************");
n-=j;

}
//---------------------------------------------------------------------------
int Reg::Sign(int digit)
{
if (digit<0)return -1;
else return 1;
}

double Reg::Sign(double digit)
{
if (digit<0)return -1.;
else return 1.;
}


//---------------------------------------------------------------------------
double Reg::ObuslMatr(double **matr,int nmatr)
{

ap::real_2d_array tmp;
tmp.setbounds(1, nmatr, 1, nmatr);
for(int i=0;i<nmatr;i++)
 for(int j=0;j<nmatr;j++)
  tmp(i+1,j+1)=matr[i][j];
double Norma=rcond1(tmp,nmatr);
return Norma;
}
//---------------------------------------------------------------------------
bool Reg::singular(double **matr,int nmatr,double *Vr)
{
bool result;
ap::real_2d_array tmp1;
tmp1.setbounds(1, nmatr, 1, nmatr);
for(int i=0;i<nmatr;i++)
 for(int j=0;j<nmatr;j++)
  tmp1(i+1,j+1)=matr[i][j];
//  tmp1(i+1,j+1)=double(random(15));
ap::real_1d_array w1;
ap::real_2d_array u1;
ap::real_2d_array vt1;

w1.setbounds(1,nmatr);
u1.setbounds(1, nmatr, 1, nmatr);
vt1.setbounds(1, nmatr, 1, nmatr);
//balancematrix(tmp1, nmatr, w1);
result=svddecomposition(tmp1,nmatr,nmatr,0,1,w1,u1,vt1);

//result=svddecomposition_old(tmp1,nmatr,nmatr,w1,vt1);
/*
double minsingular=9999999999999999999999;
int minsingularNo=0;
for(int i=1;i<=nmatr;i++)
 if(w1(i)<minsingular)
 {
  minsingular=w1(i);
  minsingularNo=i;
 }
*/

//Form1->Memo1->Lines->Add("Singular");
double NormV=0;
for(int i=1;i<=nmatr;i++)
 {
  NormV+=vt1(nmatr,i)*vt1(nmatr,i);
 }
NormV=sqrt(NormV);
if(NormV>=1)
 for(int i=1;i<=nmatr;i++)
  {
//   Form1->Memo1->Lines->Add(AnsiString(vt1(nmatr,i)));
   //Vr[i-1]=vt1(nmatr,i)/NormV;
   Vr[i-1]=vt1(nmatr,i);
  }
else
 for(int i=1;i<=nmatr;i++)
  {
   Vr[i-1]=vt1(nmatr,i)/NormV;
  }

//Form1->Memo1->Lines->Add("Singular-End");

return result;
}
//---------------------------------------------------------------------------

bool Reg::LU(double **a,double *b,double *x,int nmatr)
{
ap::real_2d_array tmp_a;
ap::real_1d_array tmp_b;
ap::real_1d_array tmp_x;
tmp_b.setbounds(1,nmatr);
tmp_x.setbounds(1,nmatr);
tmp_a.setbounds(1, nmatr, 1, nmatr);
for (int i=0;i<nmatr;i++)
{
tmp_b(i+1)=b[i];
for (int j=0;j<nmatr;j++)
 tmp_a(i+1,j+1)=a[i][j];
}
bool result=solvesystem(tmp_a,tmp_b,nmatr,tmp_x);
if(result==true)
for(int i=0;i<nmatr;i++)
 x[i]=tmp_x(i+1);
return result;
}
