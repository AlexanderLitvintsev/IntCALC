//#include "graph_cpp.h"
//#include "sobstvqri.h"
#include <complex.h>
#include <Controls.hpp>
#include <ExtCtrls.hpp>

//---------------------------------------------------------------------------
using namespace std;

struct Prms
{
int Node1,Node2;
complex<double> SFlow1,SFlow2,DS,DSTM,DSTC,I;
double Ql;
double DU;
bool IsLine,IsTrf;
};

//---------------------------------------------------------------------------
class GesseParams
{
public:
        int RaschType;
        int nJac;
        int n;
        int n_p,n_q,n_u,n_f;
        int kG;
        int *Zp;
        int *Zq;
        int *Zu,*Zf;
        int *Zp_r;
        int *Zq_r;
        int *Zu_r;
        int *Zf_r;
        GesseParams(int SetRaschtype, int SetnJac, int Setn,int *mhu);
        ~GesseParams();
        int MaxNebalNo;
        int MaxSummNo;

};
//---------------------------------------------------------------------------

class Reg
{
public:
        Reg(int nuzl);
        ~Reg();
        int kt,kl,n,ParNumber,nb,ng;
        int nl;
        int nJac,oldnJac;
        int kShReact;
        double f;
        double s;
        double Kn;
        int *mhu;
        int *mhu_n;
        int **msw;
        double *a0,*a1,*a2,*b0,*b1,*b2,ubaz,eps;
        complex<double> *sng, *sng_r, *unom,*urab,*irab,**cktr,**y,**ys,**z,**ysTemp;
        double *Pg,*Qg,*dPg,*dQn,*dPn;
        int *Nu;
        double **QminQmax;
        double T;

        double **Jac;
        double **Ges;
        double **tmpMatr;
        double *Nebal;
        double *R;
        double *dU;
        double *dD;
        double *dR;
        int *Zu,*Zf;
        struct Prms *Params;
        bool RegimState;
        double NormU;
        bool BkUse;
        bool CurrReg;

        int YsCreate();
        void ZCreate(int *ki,int *kj,int &ii,int &jj);
        void ObrMatr(double **Matr,double **ObrMatr,int razm);
        void ObrMatr(complex<double> **Matr,complex<double> **ObrMatr,int razm);

        void Obnul();
        void ResultsCalc();
        void JacobiCalc();
        void Jacobi_KnCalc();
        void Jacobi2Calc();//Расчет матр. якоби для модулей и фаз
        void JacobiX2Calc();//Расчет расширенной матрицы Якоби
        void JacobiX22Calc();//Расчет расширенной матрицы Якоби - моя
        void JacobiX2TCalc();//Расчет расширенной матрицы Якоби с параметром Т
        void JacobiShow(TMemo *Memo);
        void YsShow(TMemo *Memo);
        void NjNbNg();
        void SHN(bool stattrue);
        void SHN_Hard(bool stattrue,double T);
        bool Limits(bool limtype,bool useT,double T);
        double Norma(double **Matr,int razm,int skip);
        double Norma(complex<double> **Matr,int razm,int skip);
        double fi(int i,double y1,double y2,double f1,double f2);
        double fin(int i,double y1,double y2,double f1,double f2);

        void GesseCalc(GesseParams *GP);
        double BkCalc(GesseParams *GP,double *DX,double *Nebal);

        void CheckRes(TMemo *Memo);
        int Gauss();
        int Gauss2();
        void Gauss(complex<double> **Matr,complex<double> *Rez,int razm,int *uzlNo);
        bool Gauss(double **Matr,double *Rez,int razm);
        bool QR(double **a,double *B,double *X,int N,bool Z);
        bool LU(double **a,double *b,double *x,int nmatr);
        bool RaschCheck(double *Rez,complex<double> *urab,int razm,int NoKn,int Raschtype,bool LimitNebals,double dUmax,double dfmax,double dKnmax,double dTmax,double dVri);
        bool RaschCheck_new(bool BkUse, double Bk,double *Rez,complex<double> *urab,int razm,int NoKn,int Raschtype,bool LimitNebals,double dUmax,double dfmax,double dKnmax,double dTmax,double dVri);
        void RaschR0(double **Matr,int nMatr,double *R);
        void RaschSv(double **Matr,int nMatr,double *S,bool traspMatr);

        bool ReorganizeJacobi();
        bool ReorganizeMatr(complex<double> **Matr,int razm);
        bool ReorganizeMatr(double **Matr,int razm);
        void ReCreateYs();
        bool SobstvQri(double **matr,double *evr,double *evi,double **v,int nrazm, bool isShow,TMemo *ShowMemo);
        int GetMsw(int i,int j,int k);

        void Initialize(int Nodes);
        void InitializeJac(int Nodes);

        void EkvivalentToStar();
        //знак числа
        int Sign(int digit);
        double Sign(double digit);
        double ObuslMatr(double **matr,int nmatr);
        bool singular(double **matr,int nmatr,double *Vr);
        

};

extern Reg *Rasch;
extern Reg *RRasch;
