#include "graph_cpp.h"
#include "sobstv.h"

//---------------------------------------------------------------------------
//int SobstvZn(int n,complex<double> A[][maxArrayNo],complex<double> *Ev,complex<double> V[][maxArrayNo])
int SobstvZn(int n,complex<double> **A,complex<double> *Ev,complex<double> **V)
{
int ir1,ir2,a,b,c,d,e,f,g,h,r,s,t,u,v,w,x,y,z,ba,bb,bc,bd,be,bf,min0,na,nb,nc,nr,nt;
double rab[maxArrayNo];
double k,l,m,o,p,q,d_sqrt,d_cabs,d_abs,bg,bh,bj,bk,bl,bm,bn,bo,bp,bq,d_real,sys051;
double d_aimag,bi,f0,f1,f2,f3;
double tmp1,tmp2;
complex<double> cmplx,bs,c_sqrt;
int ierr;
sys051=1.192093e-07;
ir1=1;
ir2=n;
f=ir2-1;
h=ir1+1;
if(f<h)
  goto lb12;
for(c=h;c<=f;c++)
  {
   nc=n+c;
   m=0;
   rab[c]=0;
   rab[nc]=0;
   q=0;
   for(a=c;a<=ir2;a++)
     q+=fabs(real(A[a][c-1]))+fabs(imag(A[a][c-1]));

   if(q==0)
     goto lb11;
   g=c+ir2;
   for(d=c;d<=ir2;d++)
     {
      a=g-d;
      na=n+a;
      rab[a]=real(A[a][c-1])/q;
      rab[na]=imag(A[a][c-1])/q;
      m+=rab[a]*rab[a]+rab[na]*rab[na];
     }
   l=sqrt(m);
   k=abs(complex<double>(rab[c],rab[nc]));
   if(k==0)
     goto lb3;
   m+=k*l;
   l=l/k;
   rab[c]=(1.+l)*rab[c];
   rab[nc]=(1.+l)*rab[nc];
   goto lb4;
   lb3:
   rab[c]=l;
   A[c][c-1]=complex<double>(q,imag(A[c][c-1]));//real A[c][c-1]=q;
   lb4:
   for(b=c;b<=n;b++)
     {
      p=0;
      o=0;
      for(d=c;d<=ir2;d++)
        {
         a=g-d;
         na=n+a;
         p+=rab[a]*real(A[a][b])+rab[na]*imag(A[a][b]);
         o+=rab[a]*imag(A[a][b])-rab[na]*real(A[a][b]);
        }
      p=p/m;
      o=o/m;
      for(a=c;a<=ir2;a++)
        {
         na=n+a;
         A[a][b]=complex<double>(real(A[a][b])-p*rab[a]+o*rab[na],imag(A[a][b])-p*rab[na]-o*rab[a]);
        }
     }
      for(a=1;a<=ir2;a++)
        {
         p=0;
         o=0;
         for(e=c;e<=ir2;e++)
           {
            b=g-e;
            nb=n+b;
            p+=rab[b]*real(A[a][b])-rab[nb]*imag(A[a][b]);
            o+=rab[b]*imag(A[a][b])+rab[nb]*real(A[a][b]);
           }
         p=p/m;
         o=o/m;
         for(b=c;b<=ir2;b++)
           {
            nb=n+b;
            A[a][b]=complex<double>(real(A[a][b])-p*rab[b]-o*rab[nb],imag(A[a][b])+p*rab[nb]-o*rab[b]);
           }
        }
        rab[c]=q*rab[c];
        rab[nc]=q*rab[nc];
        A[c][c-1]=complex<double>(real(A[c][c-1])*-l,imag(A[c][c-1])*-l);
        lb11:
  }
lb12:
ierr=0;
for(r=1;r<=n;r++)
  for(s=1;s<=n;s++)
    {
     if(r!=s)
       V[r][s]=complex<double>(0,0);
     else
       V[r][s]=complex<double>(1,0);
    }
bf=ir2-ir1-1;
if(bf<0)
  goto lb25;
else if(bf==0)
  goto lb20;
else
  goto lb14;
lb14:
for(x=1;x<=bf;x++)
  {
   r=ir2-x;
   nr=n+r;
   if((rab[r]==0&&rab[nr]==0)||(real(A[r][r-1])==0&&imag(A[r][r-1])==0))
     continue;
   bq=real(A[r][r-1])*rab[r]+imag(A[r][r-1])*rab[nr];
   bb=r+1;
   for(t=bb;t<=ir2;t++)
     {
      nt=n+t;
      rab[t]=real(A[t][r-1]);
      rab[nt]=imag(A[t][r-1]);
     }
   for(s=r;s<=ir2;s++)
     {
      bh=0;
      bg=0;
      for(t=r;t<=ir2;t++)
        {
         nt=n+t;
         bh+=rab[t]*real(V[t][s])+rab[nt]*imag(V[t][s]);
         bg+=rab[t]*imag(V[t][s])-rab[nt]*real(V[t][s]);
        }
      bh=bh/bq;
      bg=bg/bq;
      for(t=r;t<=ir2;t++)
        {
         nt=n+t;
         V[t][s]+=complex<double>(bh*rab[t]-bg*rab[nt],bh*rab[nt]+bg*rab[t]);
        }
     }
  }
lb20:
u=ir1+1;
for(r=u;r<=ir2;r++)
  {
   z=min(r+1,ir2);
   if(imag(A[r][r-1])==0)
      continue;
   bq=abs(A[r][r-1]);
   bn=real(A[r][r-1])/bq;
   bm=imag(A[r][r-1])/bq;
   A[r][r-1]=complex<double>(bq,0);
   for(s=r;s<=n;s++)
     {
      bg=bn*imag(A[r][s])-bm*real(A[r][s]);
      A[r][s]=complex<double>(bn*real(A[r][s])+bm*imag(A[r][s]),bg);
     }
   for(s=1;s<=z;s++)
     {
      bg=bn*imag(A[s][r])+bm*real(A[s][r]);
      A[s][r]=complex<double>(bn*real(A[s][r])-bm*imag(A[s][r]),bg);
     }
   for(s=ir1;s<=ir2;s++)
     {
      bg=bn*imag(V[s][r])+bm*real(V[s][r]);
      V[s][r]=complex<double>(bn*real(V[s][r])-bm*imag(V[s][r]),bg);
     }
  }
lb25:
for(r=1;r<=n;r++)
  {
   if(r>=ir1&&r<=ir2)
     continue;
   Ev[r]=A[r][r];
  }
w=ir2;
bj=0;
bi=0;
lb27:
if(w<ir1)
  goto lb46;
bc=0;
be=w-1;
lb28:
for(z=ir1;z<=w;z++)
  {
   u=w+ir1-z;
   if(u==ir1)
     goto lb30;
   if(fabs(real(A[u][u-1]))<sys051*(fabs(real(A[u-1][u-1]))+fabs(imag(A[u-1][u-1]))+fabs(real(A[u][u]))+fabs(imag(A[u][u]))))
     goto lb30;
  }
lb30:
if(u==w)
  goto lb45;
if(bc==30)
  goto lb56;
if(bc==10||bc==20)
  goto lb32;
bh=real(A[w][w]);
bg=imag(A[w][w]);
bl=real(A[be][w])*real(A[w][be]);
bk=imag(A[be][w])*real(A[w][be]);
if(bl==0&&bk==0)
  goto lb33;
//bn=real(A[be][be]-complex<double>(bh,bg))/2.;
//bm=imag(A[be][be]-complex<double>(bh,bg))/2.;
bn=(real(A[be][be])-bh)/2.;
bm=(imag(A[be][be])-bg)/2.;

bs=sqrt(complex<double>(bn*bn-bm*bm+bl,2.*bn*bm+bk));
bp=real(bs);
bo=imag(bs);
if((bn*bp+bm*bo)>=0)
  goto lb31;
bp=-bp;
bo=-bo;

lb31:
aa02c(bl,bk,(bn+bp),(bm+bo),&f1,&f2);
bh-=f1;
bg-=f2;
goto lb33;

lb32:
bh=fabs(real(A[w][be]))+fabs(real(A[be][w-2]));
bg=0;

lb33:
for(r=ir1;r<=w;r++)
  A[r][r]-=complex<double>(bh,bg);
bj+=bh;
bi+=bg;
bc+=1;
bd=u+1;
for(r=bd;r<=w;r++)
  {
   bh=real(A[r][r-1]);
   A[r][r-1]=complex<double>(0,imag(A[r][r-1]));
   f1=fabs(real(A[r-1][r-1]));
   f2=fabs(imag(A[r-1][r-1]));
   f3=fabs(bh);
   f0=f1;
   if(f0<f2)
     f0=f2;
   if(f0<f3)
     f0=f3;
   f1=f1/f0;
   f2=f2/f0;
   f3=f3/f0;
   bq=f0*sqrt(f1*f1+f2*f2+f3*f3);
   bl=real(A[r-1][r-1])/bq;
   bk=imag(A[r-1][r-1])/bq;
   Ev[r-1]=complex<double>(bl,bk);
   A[r-1][r-1]=complex<double>(bq,0);
   A[r][r-1]=complex<double>(real(A[r][r-1]),bh/bq);
   for(s=r;s<=n;s++)
     {
      bn=real(A[r-1][s]);
      bm=imag(A[r-1][s]);
      bp=real(A[r][s]);
      bo=imag(A[r][s]);
      A[r-1][s]=complex<double>(bl*bn+bk*bm+imag(A[r][r-1])*bp,bl*bm-bk*bn+imag(A[r][r-1])*bo);
      A[r][s]=complex<double>(bl*bp-bk*bo-imag(A[r][r-1])*bn,bl*bo+bk*bp-imag(A[r][r-1])*bm);
     }
  }
bg=imag(A[w][w]);
if(bg==0)
  goto lb38;
bq=abs(complex<double>(real(A[w][w]),bg));
bh=real(A[w][w])/bq;
bg=bg/bq;
A[w][w]=complex<double>(bq,0);
if(w==n)
  goto lb38;
bb=w+1;
for(s=bb;s<=n;s++)
  {
   bn=real(A[w][s]);
   bm=imag(A[w][s]);
   A[w][s]=complex<double>(bh*bn+bg*bm,bh*bm-bg*bn);
  }
lb38:
for(s=bd;s<=w;s++)
  {
   bl=real(Ev[s-1]);
   bk=imag(Ev[s-1]);
   for(r=1;r<=s;r++)
     {
      bn=real(A[r][s-1]);
      bm=0;
      bp=real(A[r][s]);
      bo=imag(A[r][s]);
      if(r==s)
        goto lb39;
      bm=imag(A[r][s-1]);
      A[r][s-1]=complex<double>(real(A[r][s-1]),bl*bm+bk*bn+imag(A[s][s-1])*bo);
      lb39:
      A[r][s-1]=complex<double>(bl*bn-bk*bm+imag(A[s][s-1])*bp,imag(A[r][s-1]));
      A[r][s]=complex<double>(bl*bp+bk*bo-imag(A[s][s-1])*bn,bl*bo-bk*bp-imag(A[s][s-1])*bm);
     }
   for(r=ir1;r<=ir2;r++)
     {
      bn=real(V[r][s-1]);
      bm=imag(V[r][s-1]);
      bp=real(V[r][s]);
      bo=imag(V[r][s]);
      V[r][s-1]=complex<double>(bl*bn-bk*bm+imag(A[s][s-1])*bp,bl*bm+bk*bn+imag(A[s][s-1])*bo);
      V[r][s]=complex<double>(bl*bp+bk*bo-imag(A[s][s-1])*bn,bl*bo-bk*bp-imag(A[s][s-1])*bm);
     }
  }
if(bg==0)
  goto lb28;
for(r=1;r<=w;r++)
  {
   bn=real(A[r][w]);
   bm=imag(A[r][s]);
   A[r][w]=complex<double>(bh*bn-bg*bm,bh*bm+bg*bn);
  }
for(r=ir1;r<=ir2;r++)
  {
   bn=real(V[r][w]);
   bm=imag(V[r][w]);
   V[r][w]=complex<double>(bh*bn-bg*bm,bg*bm+bg*bn);
  }
goto lb28;

lb45:
A[w][w]+=complex<double>(bj,bi);
Ev[w]=A[w][w];
w=be;
goto lb27;

lb46:
bq=0;
for(r=1;r<=n;r++)
  for(s=r;s<=n;s++)
    {
     bq+=fabs(real(A[r][s]))+fabs(imag(A[r][s]));
    }
if(n==1||bq==0)
  goto lb57;
for(ba=2;ba<=n;ba++)
  {
   w=n+2-ba;
   bl=real(Ev[w]);
   bk=imag(Ev[w]);
   be=w-1;
   for(x=1;x<=be;x++)
     {
      r=w-x;
      bp=real(A[r][w]);
      bo=imag(A[r][w]);
      if(r==be)
        goto lb49;
      bb=r+1;
      for(s=bb;s<=be;s++)
        {
         bp+=real(A[r][s])*real(A[s][w])-imag(A[r][s])*imag(A[s][w]);
         bo+=real(A[r][s])*imag(A[s][w])+imag(A[r][s])*real(A[s][w]);
        }
      lb49:
      bn=bl-real(Ev[r]);
      bm=bk-imag(Ev[r]);
      if(bn==0&&bm==0)
        bn=sys051*bq;
      aa02c(bp,bo,bn,bm,&tmp1,&tmp2);
      A[r][w]=complex<double>(tmp1,tmp2);
     }
  }
be=n-1;
for(r=1;r<=be;r++)
  {
   if(r>=ir1&&r<=ir2)
     continue;
   bb=r+1;
   for(s=bb;s<=n;s++)
     {
      V[r][s]=A[r][s];
     }
  }
for(y=ir1;y<=be;y++)
  {
   s=n+ir1-y;
   v=min(s-1,ir2);
   for(r=ir1;r<=ir2;r++)
     {
      bp=real(V[r][s]);
      bo=imag(V[r][s]);
      for(t=ir1;t<=v;t++)
        {
         bp+=real(V[r][t])*real(A[t][s])-imag(V[r][t])*imag(A[t][s]);
         bo+=real(V[r][t])*imag(A[t][s])+imag(V[r][t])*real(A[t][s]);
        }
      V[r][s]=complex<double>(bp,bo);
     }
  }
goto lb57;

lb56:
ierr=w;
lb57:
return ierr;
}
//---------------------------------------------------------------------------


void aa02c(double xr,double xi, double yr,double yi,double *zr,double *zi)
{
double f,s;

if(fabs(yr)>fabs(yi))
  {
   f=yi/yr;
   s=f*yi+yr;
   *zr=(f*xi+xr)/s;
   *zi=(xi-f*xr)/s;
  }
else
  {
   f=yr/yi;
   s=f*yr+yi;
   *zr=(f*xr+xi)/s;
   *zi=(f*xi-xr)/s;
  }

}

//---------------------------------------------------------------------------