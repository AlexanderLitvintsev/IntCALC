//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "RegimTable.h"
#include "graph_parm.h"
#include "VectorR0.h"
//#include "graph_cpp.h"
//#include "reg.h"
//#include "Rasch.h"
//#include "graph_cpp.h"
#include "reg.h"
#include "Values.h"
#include "selectdir.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFormR *FormR;

//---------------------------------------------------------------------------
__fastcall TFormR::TFormR(TComponent* Owner)
        : TForm(Owner)
{
Vr=new double [1];
VrTemp=new double [1];
//VrRes=new double [1];
dSdi=new double [1];
Svr=new double [1];
//Svi=new double [1];
VrResTrasp=new double [1];
mhu=new int [1];
StarFilesPath=PrmsFileName=GetCurrentDir();
PrmsFileName+="\\params_default.csv";
newSchema=true;
OldRegimType=-1;
}
//---------------------------------------------------------------------------
void __fastcall TFormR::Button2Click(TObject *Sender)
{
ClearGrid(StringGrid1);
ClearGrid(StringGrid2);
ClearGrid(StringGrid3);
ClearGrid(StringGrid4);
ClearGrid(StringGrid5);
ClearGrid(StringGridSensor);
ClearGrid(StringGridVr);
ClearGrid(StringGridGradLink);
LabelRegim->Caption="Режим не рассчитан!";
LabelRegim->Font->Color=clRed;
LabelIters->Caption=0;
LabelResSensor->Caption=0;
LabelResKn->Caption=0;
LabelResTpred->Caption=0;
LabelResGradLink->Caption=0;
Form1->CreateRasch();
FillNodes(StringGrid1,StringGrid3);
FillLinks(StringGrid2);
newVr=true;
newSchema=true;
}

//---------------------------------------------------------------------------

void TFormR::FillNodes(TStringGrid *SG1,TStringGrid *SG2)
{
int i;
//Rasch=new Reg;
//Form1->CreateRasch();
//Form1->CopyRasch();
SG1->RowCount=Rasch->n+1;
SG2->RowCount=Rasch->n+1;

SG1->Cells[0][0]="№";
SG1->Cells[1][0]="№ узла";
SG1->Cells[2][0]="Тип узла";
SG1->Cells[3][0]="Uном";
SG1->Cells[4][0]="Pн";
SG1->Cells[5][0]="Qн";
SG1->Cells[6][0]="Pг";
SG1->Cells[7][0]="Qг";
SG1->Cells[8][0]="Gш";
SG1->Cells[9][0]="Bш";
SG1->Cells[10][0]="Qmin";
SG1->Cells[11][0]="Qmax";

SG2->ColCount=13;
SG2->Cells[0][0]="№";
SG2->Cells[1][0]="№ узла";
SG2->Cells[2][0]="A0";
SG2->Cells[3][0]="A1";
SG2->Cells[4][0]="A2";
SG2->Cells[5][0]="B0";
SG2->Cells[6][0]="B1";
SG2->Cells[7][0]="B2";
//SG2->Cells[8][0]="Qmin";
//SG2->Cells[9][0]="Qmax";
SG2->Cells[8][0]="dPн";
SG2->Cells[9][0]="dQн";
SG2->Cells[10][0]="dPг";

for(i=0;i<Rasch->n;i++)
 {
  SG1->Cells[0][i+1]=i+1;
  if(Rasch->Nu[i]>0)
    SG1->Cells[1][i+1]=AnsiString (Rasch->Nu[i]);
//  else
  SG1->Cells[2][i+1]=AnsiString(Rasch->mhu[i]);
  SG1->Cells[3][i+1]=FloatToStrF(abs(Rasch->unom[i]),ffFixed,10,4);
  SG1->Cells[4][i+1]=FloatToStrF(real(Rasch->sng[i]),ffFixed,10,4);
  SG1->Cells[5][i+1]=FloatToStrF(-imag(Rasch->sng[i]),ffFixed,10,4);
  SG1->Cells[6][i+1]=FloatToStrF(Rasch->Pg[i],ffFixed,10,4);
  SG1->Cells[7][i+1]=FloatToStrF(Rasch->Qg[i],ffFixed,10,4);
  SG1->Cells[8][i+1]=FloatToStrF(real(Rasch->y[i][i]),ffFixed,10,8);
  SG1->Cells[9][i+1]=FloatToStrF(-imag(Rasch->y[i][i]),ffFixed,10,8);
  SG1->Cells[10][i+1]=FloatToStrF(Rasch->QminQmax[0][i],ffFixed,10,4);
  SG1->Cells[11][i+1]=FloatToStrF(Rasch->QminQmax[1][i],ffFixed,10,4);

  SG2->Cells[0][i+1]=i+1;
  if(Rasch->Nu[i]>0)
    SG2->Cells[1][i+1]=AnsiString (Rasch->Nu[i]);
  SG2->Cells[2][i+1]=FloatToStrF(Rasch->a0[i],ffFixed,10,4);
  SG2->Cells[3][i+1]=FloatToStrF(Rasch->a1[i],ffFixed,10,4);
  SG2->Cells[4][i+1]=FloatToStrF(Rasch->a2[i],ffFixed,10,4);
  SG2->Cells[5][i+1]=FloatToStrF(Rasch->b0[i],ffFixed,10,4);
  SG2->Cells[6][i+1]=FloatToStrF(Rasch->b1[i],ffFixed,10,4);
  SG2->Cells[7][i+1]=FloatToStrF(Rasch->b2[i],ffFixed,10,4);
//  SG2->Cells[8][i+1]=FloatToStrF(Rasch->QminQmax[0][i],ffFixed,10,4);
//  SG2->Cells[9][i+1]=FloatToStrF(Rasch->QminQmax[1][i],ffFixed,10,4);
  SG2->Cells[8][i+1]=FloatToStrF(Rasch->dPn[i],ffFixed,10,4);
  SG2->Cells[9][i+1]=FloatToStrF(Rasch->dQn[i],ffFixed,10,4);
  SG2->Cells[10][i+1]=FloatToStrF(Rasch->dPg[i],ffFixed,10,4);
 }
//Rasch=NULL;
//delete Rasch;
}
//---------------------------------------------------------------------------

void TFormR::FillLinks(TStringGrid *SG)
{
int i,j,k;
double tmpD;
complex<double> tmpCD;
//Rasch=new Reg;
//Form1->CreateRasch();
//Form1->CopyRasch();
SG->RowCount=1;

SG->Cells[0][0]="№i";
SG->Cells[1][0]="№j";
SG->Cells[2][0]="№ узла i";
SG->Cells[3][0]="№ узла j";
SG->Cells[4][0]="Rij";
SG->Cells[5][0]="Xij";
SG->Cells[6][0]="Gij";
SG->Cells[7][0]="Bij";
SG->Cells[8][0]="Kтр'";
SG->Cells[9][0]="Kтр''";

k=1;
for(i=0;i<Rasch->n;i++)
   for(j=i+1;j<Rasch->n;j++)
    {
     if(Rasch->msw[i][j]==2)
       {
        tmpCD=complex<double>(1,0)/Rasch->y[i][j];
        SG->Cells[4][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[5][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
        tmpCD=2.*Rasch->y[j][i];
        SG->Cells[6][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[7][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
       }
     else if(Rasch->msw[i][j]==1)
       {
        tmpCD=complex<double>(1,0)/Rasch->y[i][j];
        SG->Cells[4][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[5][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
        tmpCD=Rasch->y[j][i];
        SG->Cells[6][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[7][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
        tmpCD=complex<double>(1,0)/Rasch->cktr[i][j];
        SG->Cells[8][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[9][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
       }
     else if(Rasch->msw[i][j]==-1)
       {
        tmpCD=complex<double>(1,0)/Rasch->y[j][i];
        SG->Cells[4][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[5][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
        tmpCD=Rasch->y[i][j];
        SG->Cells[6][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[7][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
        tmpCD=complex<double>(1,0)/Rasch->cktr[i][j];
        SG->Cells[8][k]=FloatToStrF(real(tmpCD),ffFixed,10,8);
        SG->Cells[9][k]=FloatToStrF(imag(tmpCD),ffFixed,10,8);
       }
     if(Rasch->msw[i][j]!=0)
       {
       SG->Cells[0][k]=AnsiString(i+1);
       SG->Cells[1][k]=AnsiString(j+1);
       if(Rasch->Nu[i]>0)
         SG->Cells[2][k]=AnsiString(Rasch->Nu[i]);
       if(Rasch->Nu[j]>0)
         SG->Cells[3][k]=AnsiString(Rasch->Nu[j]);
       k++;
       SG->RowCount++;
       }
    }
//Rasch=NULL;
//delete Rasch;
}
//---------------------------------------------------------------------------

void __fastcall TFormR::PageControlChange(TObject *Sender)
{
//TPageControl
if (PageControl->ActivePageIndex==0)
 {
  NewB->Caption="Новый узел";
  DelB->Caption="Удалить узлы";
  NewB->Enabled=true;
  DelB->Enabled=true;
 }
else if (PageControl->ActivePageIndex==1)
 {
  NewB->Caption="Новая связь";
  DelB->Caption="Удалить связи";
  NewB->Enabled=true;
  DelB->Enabled=true;
 }
else
 {
  NewB->Enabled=false;
  DelB->Enabled=false;
 }

if(PageControl->ActivePageIndex!=3&&PageControl->ActivePageIndex!=4)
 {Button1->Enabled=true;Button3->Enabled=true;}
else
 {Button1->Enabled=false;Button3->Enabled=false;}
}
//---------------------------------------------------------------------------

void __fastcall TFormR::NewBClick(TObject *Sender)
{
if (PageControl->ActivePageIndex==0)
  StringGrid1->RowCount+=1;
else if (PageControl->ActivePageIndex==1)
  StringGrid2->RowCount+=1;
}
//---------------------------------------------------------------------------

void __fastcall TFormR::DelBClick(TObject *Sender)
{
int i;
TGridRect myRect;
if (PageControl->ActivePageIndex==0)
 {
  myRect=StringGrid1->Selection;
  for(i=myRect.Bottom+1;i<StringGrid1->RowCount;i++)
    StringGrid1->Rows[i-(myRect.Bottom-myRect.Top+1)]=StringGrid1->Rows[i];
  StringGrid1->RowCount-=myRect.Bottom-myRect.Top+1;
  Sort_i(StringGrid1);
 }
else if(PageControl->ActivePageIndex==1)
 {
  myRect=StringGrid2->Selection;
  for(i=myRect.Bottom+1;i<StringGrid2->RowCount;i++)
    StringGrid2->Rows[i-(myRect.Bottom-myRect.Top+1)]=StringGrid2->Rows[i];
  StringGrid2->RowCount-=myRect.Bottom-myRect.Top+1;
 }
}
//---------------------------------------------------------------------------
void TFormR::Sort_i(TStringGrid *StringGrid)
{
int i,j,k;
k=10000000000;
TStrings *tmp=new TStringList();
try
{
for(i=1;i<StringGrid->RowCount;i++)
{
k=10000000000;
for(j=i;j<StringGrid->RowCount;j++)
 if(StrToInt(StringGrid->Cells[0][j])<k)
 {
  k=StringGrid->Cells[0][j].ToInt();
  if(i!=j)
    {
    tmp->Assign(StringGrid->Rows[j]);
    StringGrid->Rows[j]=StringGrid->Rows[i];
    StringGrid->Rows[i]=tmp;
    }
 }
}
}
__finally
 {delete tmp;}

}
//---------------------------------------------------------------------------
void TFormR::Sort_ij(TStringGrid *StringGrid)
{
int i,j,k,l;
Sort_i(StringGrid);
k=10000000000;
TStrings *tmp=new TStringList();
try
{
for(i=1;i<StringGrid->RowCount;i++)
{
l=StrToInt(StringGrid->Cells[0][i]);
k=10000000000;
for(j=i;j<StringGrid->RowCount;j++)
 if(StrToInt(StringGrid->Cells[1][j])<k&&StrToInt(StringGrid->Cells[0][j])==l)
 {
  k=StringGrid->Cells[1][j].ToInt();
  if(i!=j)
    {
    tmp->Assign(StringGrid->Rows[j]);
    StringGrid->Rows[j]=StringGrid->Rows[i];
    StringGrid->Rows[i]=tmp;
    }
 }
}
}
__finally
 {delete tmp;}

}
//---------------------------------------------------------------------------


int TFormR::FindMaxNodeNumber(TStringGrid *StringGrid,int ColNumber,int CurrMax)
{
int i;
for(i=1;i<StringGrid->RowCount;i++)
 if(StringGrid->Cells[ColNumber][i].ToInt()>CurrMax)
  CurrMax=StringGrid->Cells[ColNumber][i].ToInt();
return CurrMax;
}
//---------------------------------------------------------------------------

void TFormR::WriteRRasch(bool CurrReg)
{
int i,j,k;
Sort_i(StringGrid1);
Sort_i(StringGrid3);
Sort_ij(StringGrid2);
int MaxNodeNumber=FindMaxNodeNumber(StringGrid1,0,0);
MaxNodeNumber=FindMaxNodeNumber(StringGrid3,0,MaxNodeNumber);
MaxNodeNumber=FindMaxNodeNumber(StringGrid2,0,MaxNodeNumber);
MaxNodeNumber=FindMaxNodeNumber(StringGrid2,1,MaxNodeNumber);

//Rasch=new Reg;
if(CurrReg==false)
{
delete Rasch;
Rasch=new Reg(MaxNodeNumber);
Rasch->Initialize(MaxNodeNumber);
Rasch->n=MaxNodeNumber;
Rasch->nl=StringGrid2->RowCount-1;
Rasch->Obnul();
Rasch->f=50;
Rasch->s=0;
}
Rasch->eps=Form1->precision;
Rasch->kt=0;
Rasch->kl=0;

for(i=1;i<StringGrid1->RowCount;i++)
 {
  if(StringGrid1->Cells[1][i]!="")
   Rasch->Nu[i-1]=StringGrid1->Cells[1][i].ToInt();
  Rasch->mhu[i-1]=StringGrid1->Cells[2][i].ToInt();
  Rasch->unom[i-1]=polar(double(MyStrToFloat(StringGrid1->Cells[3][i])),0.);
  Rasch->sng[i-1]=complex<double>(MyStrToFloat(StringGrid1->Cells[4][i]),-MyStrToFloat(StringGrid1->Cells[5][i]));
  Rasch->Pg[i-1]=MyStrToFloat(StringGrid1->Cells[6][i]);
  Rasch->Qg[i-1]=MyStrToFloat(StringGrid1->Cells[7][i]);
  Rasch->y[i-1][i-1]=complex<double>(MyStrToFloat(StringGrid1->Cells[8][i]),-MyStrToFloat(StringGrid1->Cells[9][i]));
  Rasch->QminQmax[0][i-1]=MyStrToFloat(StringGrid1->Cells[10][i]);
  Rasch->QminQmax[1][i-1]=MyStrToFloat(StringGrid1->Cells[11][i]);
 }
for(i=1;i<StringGrid3->RowCount;i++)
 {
  Rasch->a0[i-1]=MyStrToFloat(StringGrid3->Cells[2][i]);
  if(Rasch->a0[i-1]==0)
    Rasch->a0[i-1]=1;
  Rasch->a1[i-1]=MyStrToFloat(StringGrid3->Cells[3][i]);
  Rasch->a2[i-1]=MyStrToFloat(StringGrid3->Cells[4][i]);
  Rasch->b0[i-1]=MyStrToFloat(StringGrid3->Cells[5][i]);
  if(Rasch->b0[i-1]==0)
    Rasch->b0[i-1]=1;
  Rasch->b1[i-1]=MyStrToFloat(StringGrid3->Cells[6][i]);
  Rasch->b2[i-1]=MyStrToFloat(StringGrid3->Cells[7][i]);

  Rasch->dPn[i-1]=MyStrToFloat(StringGrid3->Cells[8][i]);
  Rasch->dQn[i-1]=MyStrToFloat(StringGrid3->Cells[9][i]);
  Rasch->dPg[i-1]=MyStrToFloat(StringGrid3->Cells[10][i]);

//  Rasch->QminQmax[0][i-1]=MyStrToFloat(StringGrid3->Cells[8][i]);
//  Rasch->QminQmax[1][i-1]=MyStrToFloat(StringGrid3->Cells[9][i]);
//  Rasch->dPn[i-1]=MyStrToFloat(StringGrid3->Cells[10][i]);
//  Rasch->dQn[i-1]=MyStrToFloat(StringGrid3->Cells[11][i]);
//  Rasch->dPg[i-1]=MyStrToFloat(StringGrid3->Cells[12][i]);
 }
for(k=1;k<StringGrid2->RowCount;k++)
 {
  i=StringGrid2->Cells[0][k].ToInt();
  j=StringGrid2->Cells[1][k].ToInt();
  Rasch->cktr[i-1][j-1]=complex<double>(MyStrToFloat(StringGrid2->Cells[8][k]),MyStrToFloat(StringGrid2->Cells[9][k]));
  if(Rasch->cktr[i-1][j-1]==complex<double>(0,0))
    Rasch->cktr[i-1][j-1]=complex<double>(1,0);
  Rasch->cktr[j-1][i-1]=Rasch->cktr[i-1][j-1];
  Rasch->cktr[i-1][j-1]=complex<double>(1,0)/Rasch->cktr[i-1][j-1];
  if(real(Rasch->cktr[i-1][j-1])>1)
   {
    Rasch->y[i-1][j-1]=complex<double>(1,0)/complex<double>(MyStrToFloat(StringGrid2->Cells[4][k]),MyStrToFloat(StringGrid2->Cells[5][k]));
    Rasch->y[j-1][i-1]=complex<double>(MyStrToFloat(StringGrid2->Cells[6][k]),MyStrToFloat(StringGrid2->Cells[7][k]));
//    Rasch->y[j-1][i-1]/=2.;
    Rasch->msw[i-1][j-1]=1;
    Rasch->msw[j-1][i-1]=-1;
    Rasch->kt++;
   }
  else if(real(Rasch->cktr[i-1][j-1])<1)
   {
    Rasch->y[i-1][j-1]=complex<double>(MyStrToFloat(StringGrid2->Cells[6][k]),MyStrToFloat(StringGrid2->Cells[7][k]));
//    Rasch->y[i-1][j-1]/=2.;
    Rasch->y[j-1][i-1]=complex<double>(1,0)/complex<double>(MyStrToFloat(StringGrid2->Cells[4][k]),MyStrToFloat(StringGrid2->Cells[5][k]));
    Rasch->msw[i-1][j-1]=-1;
    Rasch->msw[j-1][i-1]=1;
    Rasch->kt++;
   }
  else
   {
    Rasch->y[i-1][j-1]=complex<double>(1,0)/complex<double>(MyStrToFloat(StringGrid2->Cells[4][k]),MyStrToFloat(StringGrid2->Cells[5][k]));
    Rasch->y[j-1][i-1]=complex<double>(MyStrToFloat(StringGrid2->Cells[6][k]),MyStrToFloat(StringGrid2->Cells[7][k]));
    Rasch->y[j-1][i-1]/=2.;
    Rasch->msw[i-1][j-1]=2;
    Rasch->msw[j-1][i-1]=2;
    Rasch->kl++;
   }

 }
//Form1->SaveRasch();
//Rasch=NULL;
//delete Rasch;

}

//---------------------------------------------------------------------------
float TFormR::MyStrToFloat(AnsiString str)
{
float tmp;
if(str!=""&&str!=NULL)
{
for(int i=1;i<=str.Length();i++)
 if(str[i]=='.')
   str[i]=',';
 tmp=StrToFloat(str);
}
else
 tmp=0;
return tmp;
}

void TFormR::FillResults(int ItNumber)
{
int i;
//Rasch=new Reg;
//Form1->CopyRasch();
double *Pn,*Qn,*Pg,*Qg;
complex<double> tmp;
Pn=new double[Rasch->n];
Qn=new double[Rasch->n];
Pg=new double[Rasch->n];
Qg=new double[Rasch->n];

for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==1||Rasch->mhu[i]==4)
 {
  Pn[i]=real(Rasch->sng[i])*(Rasch->a0[i]+Rasch->a1[i]*(abs(Rasch->urab[i])/abs(Rasch->unom[i]))+Rasch->a2[i]*(abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
  Pg[i]=Rasch->Pg[i];
  Qn[i]=-imag(Rasch->sng[i])*(Rasch->b0[i]+Rasch->b1[i]*(abs(Rasch->urab[i])/abs(Rasch->unom[i]))+Rasch->b2[i]*(abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
  if(Rasch->mhu[i]==4)
   {
    tmp=complex<double>(0,0);
    for(int j=0;j<Rasch->n;j++)
    {
     tmp+=Rasch->ys[i][j]*Rasch->urab[j];
    }
    tmp=tmp*conj(Rasch->urab[i]);
    Qg[i]=imag(tmp)+(Qn[i]);
   }
  else
    Qg[i]=Rasch->Qg[i];
 }
 else if(Rasch->mhu[i]==0)
 {
  tmp=complex<double>(0,0);
  for(int j=0;j<Rasch->n;j++)
  {
   tmp+=Rasch->ys[i][j]*Rasch->urab[j];
  }
  tmp=tmp*conj(Rasch->urab[i]);
  Pg[i]=-real(tmp);
  Qg[i]=imag(tmp);
  Pn[i]=0;
  Qn[i]=0;
 }
 else
 {
  Pn[i]=Pg[i]=Qn[i]=Qg[i]=0;
 }
}

StringGrid5->RowCount=Rasch->n+1;
StringGrid5->Cells[0][0]="№";
StringGrid5->Cells[1][0]="№ узла";
StringGrid5->Cells[2][0]="Тип узла";
StringGrid5->Cells[3][0]="U";
StringGrid5->Cells[4][0]="f";
StringGrid5->Cells[5][0]="Pн";
StringGrid5->Cells[6][0]="Qн";
StringGrid5->Cells[7][0]="Pг";
StringGrid5->Cells[8][0]="Qг";
StringGrid5->Cells[9][0]="Qmin";
StringGrid5->Cells[10][0]="Qmax";


for(i=0;i<Rasch->n;i++)
 {
  StringGrid5->Cells[0][i+1]=i+1;
  if(Rasch->Nu[i]>0)
    StringGrid5->Cells[1][i+1]=AnsiString (Rasch->Nu[i]);
  StringGrid5->Cells[2][i+1]=AnsiString(Rasch->mhu[i]);
  StringGrid5->Cells[3][i+1]=FloatToStrF(abs(Rasch->urab[i]),ffFixed,10,4);
  StringGrid5->Cells[4][i+1]=FloatToStrF(arg(Rasch->urab[i])*180./M_PI,ffFixed,10,4);
  StringGrid5->Cells[5][i+1]=FloatToStrF(Pn[i],ffFixed,10,4);
  StringGrid5->Cells[6][i+1]=FloatToStrF(Qn[i],ffFixed,10,4);
  StringGrid5->Cells[7][i+1]=FloatToStrF(Pg[i],ffFixed,10,4);
  StringGrid5->Cells[8][i+1]=FloatToStrF(Qg[i],ffFixed,10,4);
  StringGrid5->Cells[9][i+1]=FloatToStrF(Rasch->QminQmax[0][i],ffFixed,10,4);
  StringGrid5->Cells[10][i+1]=FloatToStrF(Rasch->QminQmax[1][i],ffFixed,10,4);
 }


StringGrid4->RowCount=Vals[ItNumber].mswNumber+1;
StringGrid4->Cells[0][0]="№i";
StringGrid4->Cells[1][0]="№j";
StringGrid4->Cells[2][0]="№ узла i";
StringGrid4->Cells[3][0]="№ узла j";
StringGrid4->Cells[4][0]="Модуль тока";
StringGrid4->Cells[5][0]="Фаза тока";
StringGrid4->Cells[6][0]="Акт. ток";
StringGrid4->Cells[7][0]="Реакт. ток";
StringGrid4->Cells[8][0]="Pi";
StringGrid4->Cells[9][0]="Qi";
StringGrid4->Cells[10][0]="Pj";
StringGrid4->Cells[11][0]="Qj";
StringGrid4->Cells[12][0]="Зарядная мощность Ql";
StringGrid4->Cells[13][0]="Потери dP";
StringGrid4->Cells[14][0]="Потери dQ";
StringGrid4->Cells[15][0]="Потери dU";
StringGrid4->Cells[16][0]="Потери в меди dP";
StringGrid4->Cells[17][0]="Потери в меди dQ";
StringGrid4->Cells[18][0]="Потери в стали dP";
StringGrid4->Cells[19][0]="Потери в стали dQ";

for(i=0;i<Vals[ItNumber].mswNumber;i++)
{
StringGrid4->Cells[0][i+1]=AnsiString(Vals[ItNumber].Node1[i]);
StringGrid4->Cells[1][i+1]=AnsiString(Vals[ItNumber].Node2[i]);
if(Vals[ItNumber].RealNode[Vals[ItNumber].Node1[i]-1]>0)
 StringGrid4->Cells[2][i+1]=AnsiString(Vals[ItNumber].RealNode[Vals[ItNumber].Node1[i]-1]);
if(Vals[ItNumber].RealNode[Vals[ItNumber].Node2[i]-1]>0)
 StringGrid4->Cells[3][i+1]=AnsiString(Vals[ItNumber].RealNode[Vals[ItNumber].Node2[i]-1]);
StringGrid4->Cells[4][i+1]=FloatToStrF(abs(Vals[ItNumber].I[i]),ffFixed,10,4);
StringGrid4->Cells[5][i+1]=FloatToStrF(arg(Vals[ItNumber].I[i]),ffFixed,10,4);
StringGrid4->Cells[6][i+1]=FloatToStrF(real(Vals[ItNumber].I[i]),ffFixed,10,4);
StringGrid4->Cells[7][i+1]=FloatToStrF(imag(Vals[ItNumber].I[i]),ffFixed,10,4);
StringGrid4->Cells[8][i+1]=FloatToStrF(real(Vals[ItNumber].S1[i]),ffFixed,10,4);
StringGrid4->Cells[9][i+1]=FloatToStrF(-imag(Vals[ItNumber].S1[i]),ffFixed,10,4);
StringGrid4->Cells[10][i+1]=FloatToStrF(real(Vals[ItNumber].S2[i]),ffFixed,10,4);
StringGrid4->Cells[11][i+1]=FloatToStrF(-imag(Vals[ItNumber].S2[i]),ffFixed,10,4);
StringGrid4->Cells[12][i+1]=FloatToStrF(Vals[ItNumber].Ql[i],ffFixed,10,4);
StringGrid4->Cells[13][i+1]=FloatToStrF(real(Vals[ItNumber].DS[i]),ffFixed,10,4);
StringGrid4->Cells[14][i+1]=FloatToStrF(imag(Vals[ItNumber].DS[i]),ffFixed,10,4);
StringGrid4->Cells[15][i+1]=FloatToStrF(Vals[ItNumber].DU[i],ffFixed,10,4);
StringGrid4->Cells[16][i+1]=FloatToStrF(real(Vals[ItNumber].DSTM[i]),ffFixed,10,4);
StringGrid4->Cells[17][i+1]=FloatToStrF(imag(Vals[ItNumber].DSTM[i]),ffFixed,10,4);
StringGrid4->Cells[18][i+1]=FloatToStrF(real(Vals[ItNumber].DSTC[i]),ffFixed,10,4);
StringGrid4->Cells[19][i+1]=FloatToStrF(imag(Vals[ItNumber].DSTC[i]),ffFixed,10,4);
}




//Rasch=NULL;
//delete Rasch;
delete []Pn;
delete []Qn;
delete []Pg;
delete []Qg;
}


//---------------------------------------------------------------------------

void TFormR::ClearGrid (TStringGrid *StringGrid)
{
for (int i=0;i<StringGrid->RowCount;i++)
 StringGrid->Rows[i]->Clear();

}




void __fastcall TFormR::ComboBox1Change(TObject *Sender)
{
/*
if (ComboBox1->ItemIndex==4)
  CheckStarFiles->Enabled=true;
else
  CheckStarFiles->Enabled=false;
if (ComboBox1->ItemIndex==3)
  GroupPredel->Enabled=true;
else
  GroupPredel->Enabled=false;
*/  
}
//---------------------------------------------------------------------------

void __fastcall TFormR::FormActivate(TObject *Sender)
{
//ComboBox1->ItemIndex=0;
}
//---------------------------------------------------------------------------

void __fastcall TFormR::CheckNebalClick(TObject *Sender)
{
if(CheckNebal->Checked==true)
 GroupNebal->Enabled=true;
else
 GroupNebal->Enabled=false;



}
//---------------------------------------------------------------------------

bool TFormR::ReadRaschetType()
{
bool ParamsOk=true;
RaschetType=ComboBox1->ItemIndex;
MaxIters=MyStrToFloat(EditMaxIter->Text);
if(MaxIters<=0)
 {ParamsOk=false;MaxIters=50;}
prec=MyStrToFloat(EditPrec->Text)/100.;
if(prec<=0)
 {ParamsOk=false;prec=0.001;}
Showtxt=CheckTxtFile->Checked;
BkUse=CheckBkUse->Checked;
CurrReg=CheckCurrReg->Checked;

FindSensor=CheckSensor->Checked;
UseQminQmax=CheckQminQmax->Checked;
LimitNebals=CheckNebal->Checked;
dUmax=MyStrToFloat(EditdUmax->Text)/100.;
dfmax=MyStrToFloat(Editdfmax->Text);
//dfmax=MyStrToFloat(Editdfmax->Text)/100.;
dKnmax=MyStrToFloat(EditdKnmax->Text);
dTmax=MyStrToFloat(EditdTmax->Text)/100.;
dVrimax=MyStrToFloat(EditdVrimax->Text)/100.;
if(LimitNebals==true)
 {
  if(dUmax<=0||dfmax<=0)
   ParamsOk=false;
  if((RaschetType==3)&&(dKnmax<=0||dTmax<=0||dVrimax<=0))
   ParamsOk=false;
 }
PredelRaschetMetod=RadioGroupPredelType->ItemIndex;
if(RaschetType==3&&PredelRaschetMetod<0)
 ParamsOk=false;
Kn0=MyStrToFloat(EditKn0->Text);
if(RaschetType==3&&PredelRaschetMetod==0&&Kn0<=0)
 ParamsOk=false;
T0=MyStrToFloat(EditT0->Text);
if(RaschetType==3&&PredelRaschetMetod==0&&T0<=0)
 ParamsOk=false;
RaschetVr0=CheckAutoVr0->Checked;
TempRegims=CheckProm->Checked;
Gradients=CheckGrad->Checked;
ShowGraph=CheckGraph->Checked;
dT=MyStrToFloat(EditdT->Text);
dKn=MyStrToFloat(EditdKn->Text);
if(dT<=0)
 dT=0;
if(dKn<=0)
 dKn=0;
if(CheckStarFiles->Enabled==true&&CheckStarFiles->Checked==true)
 SaveStarFiles=true;
else
 SaveStarFiles=false;

if(CheckParamsToFile->Checked==true)
 SavePrmsToFile=true;
else
 SavePrmsToFile=false;

return ParamsOk;
}
//---------------------------------------------------------------------------
void TFormR::RegimRaschet()
{
bool Regim,CalcVr=false,SetCurrReg=true;
int RegType,tmptype;
int n;

if(ReadRaschetType()==true)
{
Form1FillParams();

n=FindMaxNodeNumber(StringGrid1,0,0);
n=FindMaxNodeNumber(StringGrid3,0,n);
n=FindMaxNodeNumber(StringGrid2,0,n);
n=FindMaxNodeNumber(StringGrid2,1,n);

if(RaschetType==0)
 RegType=1;
else if(RaschetType==1)
 RegType=0;
else if(RaschetType==2)
 RegType=3;
else if(RaschetType==4)
 RegType=6;

if(RaschetType==3&&PredelRaschetMetod==0)
 RegType=7;
else if(RaschetType==3&&PredelRaschetMetod==1)
 RegType=8;
else if(RaschetType==3&&PredelRaschetMetod==2)
 RegType=9;
else if(RaschetType==3&&PredelRaschetMetod==3)
 RegType=5;

if(OldRegimType!=RegType)
 CurrReg=false;
OldRegimType=RegType;

if((newSchema==false&&Rasch->RegimState==true&&CurrReg==true)!=true)
 {
  SetCurrReg=false;
  CalcVr=true;
 }
WriteRRasch(SetCurrReg);
newSchema=false;

//n=Form1->CreateRasch();


//Values::NodesNumber=n;
//Values::mswNumber=n*n;
//Vals=new Values[1];


if(CalcVr==true)
{
if(RaschetVr0==true&&RaschetType==3)
 {
  delete []Vr;
  Vr=new double[n*2];
  if(PredelRaschetMetod==0||PredelRaschetMetod==3)
   Form1->RaschetVr0(0,nVr,Vr);
  else
   Form1->RaschetVr0(1,nVr,Vr);
//  for(int ll=0;ll<nVr;ll++)
//    Vr[ll]=1;
 }
else if(RaschetType==3)
  nVr=FillVr(newVr);
 }
Regim=Form1->UstRegimRaschet(RegType,BkUse,SetCurrReg,Showtxt,Vr,nVr,SavePrmsToFile,PrmsFileName);
if(Regim==true)
{
 LabelRegim->Caption="Режим рассчитан";
 LabelRegim->Font->Color=clGreen;
 LabelIters->Caption=Iters;
 FillResults(0);
 if(FindSensor==true)
 {
  RaschetSv(nJ,nb,ng,nuzl,false);
  FillSensorGrid();
 }
 if(RaschetType==4)
 {
  TStringGrid *tmpSG1,*tmpSG2,*tmpSG3;
  tmpSG1=new TStringGrid(FormR);
  tmpSG2=new TStringGrid(FormR);
  tmpSG1->Parent=FormR;
  tmpSG2->Parent=FormR;
  tmpSG1->Visible=false;
  tmpSG2->Visible=false;
  FillNodes(tmpSG1,tmpSG2);
  SaveStringGrid(tmpSG1,0,StarFilesPath+"\\узлы.csv");
  SaveStringGrid(tmpSG2,2,StarFilesPath+"\\узлы-дополнительно.csv");
  FillLinks(tmpSG1);
  SaveStringGrid(tmpSG1,1,StarFilesPath+"\\ветви.csv");
  delete tmpSG1;
  delete tmpSG2;

 }

 if(RaschetType==3)
 {

  T=Rasch->T;
  Kn=Rasch->Kn;
  Rasch->NjNbNg();

  if(PredelRaschetMetod==0||PredelRaschetMetod==1||PredelRaschetMetod==3)
  {
   nVr=Rasch->nJac;
   VrRes=Rasch->R;
   RaschetVrResTrasp(true);
   FillVrResToGrid();
   VrRes=NULL;
  }

//  Rasch=NULL;
//  delete Rasch;
  if(PredelRaschetMetod==0||PredelRaschetMetod==3)
   LabelResKn->Caption=AnsiString(Kn);
  else
   LabelResTpred->Caption=AnsiString(T);
 }
}
else
{
 LabelRegim->Caption="Режим не рассчитан!";
 LabelRegim->Font->Color=clRed;
 LabelIters->Caption=Iters;
}
//delete []Vals;
if(Regim==true)
{
if(RaschetType==3&&TempRegims==true)
 {
  if(PredelRaschetMetod==1||PredelRaschetMetod==2)
  {
   if(dT<=0)
    {dT=T/double(49);nP=50;}
   else
    {nP=int(T/dT);dT=T/double(nP-1);}
   tmptype=0;
  }
  else
  {
   if(dKn<=0)
    {dKn=fabs(1.-Kn)/double(49);nP=50;}
   else
    {nP=int(fabs(1.-Kn)/dKn);dKn=fabs(1.-Kn)/double(nP-1);}
   tmptype=1;
  }
//  Values::NodesNumber=n;
//  Values::mswNumber=n*n;
//  Vals=new Values[nP];
  Form1->RaschetSeries(tmptype,T,Kn,nP,false);
  if(Gradients==true)
   {
    FinddSdiMax();
    FilldSdiToGrid();

   }
//  Form5->Vals=Vals;

  Form5->ItersNumber=nP;
  delete []ValsGraph;
  ValsGraph=Vals;
  Vals=NULL;
  delete []Vals;
  Vals=new Values[1];
  if(ShowGraph==true)
    Form5->Show();

//  Vals=NULL;
//  delete []Vals;
 }
}

//delete []Vr;
}
else
 ShowMessage("Расчет невозможен - проверьте параметры!");

}

//---------------------------------------------------------------------------
void TFormR::Form1FillParams()
{
Form1->NumberOfIters=MaxIters;
Form1->precision=prec;
Form1->UseQminQmax=UseQminQmax;
Form1->dUmax=dUmax;
Form1->dfmax=dfmax;
Form1->dTmax=dTmax;
Form1->dKnmax=dKnmax;
Form1->dVrimax=dVrimax;
Form1->Kn0=Kn0;
Form1->T0=T0;
Form1->LimitNebals=LimitNebals;
}
//---------------------------------------------------------------------------
void __fastcall TFormR::ButtonRaschetClick(TObject *Sender)
{
RegimRaschet();
PageControl->ActivePageIndex=4;
}
//---------------------------------------------------------------------------

void __fastcall TFormR::Button4Click(TObject *Sender)
{
nVr=FillVr(newVr);
FormVectorVr->Vr=Vr;
FormVectorVr->n=nVr;
if(FormVectorVr->ShowModal()==mrOk)
{
 FormVectorVr->Vr=NULL;
 delete []FormVectorVr->Vr;
}
}
//---------------------------------------------------------------------------
int TFormR::FillVr(bool &FillVr)
{
int nJ=-1;
int tmpint;

nJ=0;
for(int i=1;i<StringGrid1->RowCount;i++)
  {
   tmpint=StringGrid1->Cells[2][i].ToInt();
   if(tmpint==4)
     nJ+=1;
   else if(tmpint!=0)
     nJ+=2;
  }

tmpint=FindMaxNodeNumber(StringGrid1,0,0);
if(tmpint!=nU||FillVr==true)
{
nU=tmpint;
Vr=NULL;
delete []Vr;
Vr=new double [nJ];
for(int i=0;i<nJ;i++)
 Vr[i]=1;
FillVr=false;
}
return nJ;
}
//---------------------------------------------------------------------------

int TFormR::FinddSdiMax()
{
double maxdS=0;
int maxdSLinkNumber,i;
delete []dSdi;
dSdi=new double[Vals[nP-1].mswNumber];

for(i=0;i<Vals[nP-1].mswNumber;i++)
{
 if(Vals[nP-1].IsLine[i]==true)
  {
   dSdi[i]=(abs(Vals[nP-1].DS[i])-abs(Vals[nP-2].DS[i]));
   if(dSdi[i]>maxdS)
     {
      maxdS=dSdi[i];
      maxdSLinkNumber=i;
     }
  }
 else if(Vals[nP-1].IsTrf[i]==true)
  {
   dSdi[i]=(abs(Vals[nP-1].DSTM[i])-abs(Vals[nP-2].DSTM[i]));
   if(dSdi[i]>maxdS)
     {
      maxdS=dSdi[i];
      maxdSLinkNumber=i;
     }
  }
}
//ShowMessage(AnsiString(Vals[nP-1].RealNode[Vals[nP-1].Node1[maxdSLinkNumber]-1])+" - "+AnsiString(Vals[nP-1].RealNode[Vals[nP-1].Node2[maxdSLinkNumber]-1]));
return maxdSLinkNumber;
}
//---------------------------------------------------------------------------
void TFormR::RaschetSv(int &nJlocal,int &nblocal,int &nglocal,int &nlocal,bool transpMatr)
{
delete []Svr;
//delete []Svi;
delete []mhu;
//Rasch=new Reg;
//Form1->CopyRasch();
Rasch->YsCreate();
Rasch->NjNbNg();
Svr=new double[Rasch->nJac+1];
//Svi=new double[Rasch->nJac+1];
mhu=new int [Rasch->n];
nJlocal=Rasch->nJac;
nblocal=Rasch->nb;
nglocal=Rasch->ng;
nlocal=Rasch->n;
for(int i=0;i<Rasch->n;i++)
 mhu[i]=Rasch->mhu[i];
Rasch->Jacobi2Calc();
Rasch->RaschSv(Rasch->Jac,Rasch->nJac,Svr,transpMatr);

}
//---------------------------------------------------------------------------
void TFormR::RaschetVrResTrasp(bool transpMatr)
{
int oldnJac=Rasch->nJac;
delete []VrResTrasp;
Rasch->NjNbNg();
VrResTrasp=new double[Rasch->nJac+1];
Rasch->Jacobi2Calc();
Rasch->RaschSv(Rasch->Jac,Rasch->nJac,VrResTrasp,transpMatr);
Rasch->nJac=oldnJac;
}
//---------------------------------------------------------------------------

void TFormR::FillSensorGrid()
{
int i,k;
double max=0;
int maxi;
//int Ucount=(nJ+ng)/2-ng;

StringGridSensor->RowCount=nJ+1;
StringGridSensor->Cells[0][0]="Узел";
StringGridSensor->Cells[1][0]="S[i]";
k=0;
for(i=0;i<nuzl;i++)
{
 if(mhu[i]!=0&&mhu[i]!=4)
 {
 if(Vals[0].RealNode[i]>0)
  StringGridSensor->Cells[0][k+1]=AnsiString(Vals[0].RealNode[i])+"(модуль)";
 else
  StringGridSensor->Cells[0][k]=AnsiString(i+1)+"(модуль)";
 k++;
 }
}
for(i=0;i<nuzl;i++)
{
 if(mhu[i]!=0)
 {
 if(Vals[0].RealNode[i]>0)
  StringGridSensor->Cells[0][k+1]=AnsiString(Vals[0].RealNode[i])+"(фаза)";
 else
  StringGridSensor->Cells[0][k]=AnsiString(i+1)+"(фаза)";
 k++;
 }
}

for(i=1;i<nJ+1;i++)
{

 StringGridSensor->Cells[1][i]=FloatToStrF(Svr[i-1],ffFixed,10,8);
}

TStrings *tmp=new TStringList();
try
{
for(k=1;k<nJ+1;k++)
{
 max=0;
 for(i=k;i<nJ+1;i++)
 {
  if((MyStrToFloat(StringGridSensor->Cells[1][i])*MyStrToFloat(StringGridSensor->Cells[1][i]))>max)
    {maxi=i;max=MyStrToFloat(StringGridSensor->Cells[1][i])*MyStrToFloat(StringGridSensor->Cells[1][i]);}
 }
if(k!=maxi)
 {
  tmp->Assign(StringGridSensor->Rows[k]);
  StringGridSensor->Rows[k]=StringGridSensor->Rows[maxi];
  StringGridSensor->Rows[maxi]=tmp;
 }
}
}
__finally
{delete tmp;}
LabelResSensor->Caption=StringGridSensor->Cells[0][1];

}
//---------------------------------------------------------------------------
void TFormR::FillVrResToGrid()
{
StringGridVr->RowCount=nVr+1;
StringGridVr->Cells[0][0]="№ ";
StringGridVr->Cells[1][0]="S[i]";
StringGridVr->Cells[2][0]="R[i]";

int i;
for(i=1;i<nVr+1;i++)
 {
 StringGridVr->Cells[0][i]=AnsiString(i);
 }

for(i=1;i<nVr+1;i++)
{
 StringGridVr->Cells[1][i]=FloatToStrF(VrRes[i-1],ffFixed,10,8);
 StringGridVr->Cells[2][i]=FloatToStrF(VrResTrasp[i-1],ffFixed,10,8);
}


}
//---------------------------------------------------------------------------

void TFormR::FilldSdiToGrid()
{
int i,k;
double max;
int maxi;
StringGridGradLink->Cells[0][0]="Ветвь";
StringGridGradLink->Cells[0][1]="dS'";
StringGridGradLink->RowCount=Vals[nP-1].mswNumber+1;
for (i=1;i<Vals[nP-1].mswNumber+1;i++)
 {
  StringGridGradLink->Cells[0][i]=AnsiString(Vals[nP-1].RealNode[Vals[nP-1].Node1[i-1]-1])+" - "+AnsiString(Vals[nP-1].RealNode[Vals[nP-1].Node2[i-1]-1]);
  StringGridGradLink->Cells[1][i]=dSdi[i-1];
 }
TStrings *tmp=new TStringList();

try
{
for(k=1;k<Vals[nP-1].mswNumber+1;k++)
{
 max=-9999999999;
 for(i=k;i<Vals[nP-1].mswNumber+1;i++)
 {
  if(MyStrToFloat(StringGridGradLink->Cells[1][i])>max)
    {maxi=i;max=MyStrToFloat(StringGridGradLink->Cells[1][i]);}
 }
if(k!=maxi)
 {
  tmp->Assign(StringGridGradLink->Rows[k]);
  StringGridGradLink->Rows[k]=StringGridGradLink->Rows[maxi];
  StringGridGradLink->Rows[maxi]=tmp;
 }
}
}
__finally
{delete tmp;}
LabelResGradLink->Caption=StringGridGradLink->Cells[0][1];
}

void __fastcall TFormR::CheckPromClick(TObject *Sender)
{
if(CheckProm->Checked==true)
 {
  CheckGrad->Enabled=true;
  CheckGraph->Enabled=true;
 }
else
 {
  CheckGrad->Checked=false;
  CheckGraph->Checked=false;

  CheckGrad->Enabled=false;
  CheckGraph->Enabled=false;
 }

}
//---------------------------------------------------------------------------
void TFormR::SaveStringGrid(TStringGrid *CurrSG,int type,AnsiString FileName)
{
TStrings *tmp=new TStringList();
AnsiString tmpStr;

switch(type)
{
case 0:
 {
  tmpStr="Узлы";
  break;
 }
case 1:
 {
  tmpStr="Ветви";
  break;
 }
case 2:
 {
  tmpStr="Узлы-дополнительно";
  break;
 }
case 5:
 {
  tmpStr="Результаты-узлы";
  break;
 }
case 6:
 {
  tmpStr="Результаты-ветви";
  break;
 }

}


try
{
tmp->Add(tmpStr);
tmp->Add(AnsiString(CurrSG->RowCount));
tmp->Add(AnsiString(CurrSG->ColCount));
for(int i=0;i<CurrSG->RowCount;i++)
 {
  tmpStr="";
  for(int j=0;j<CurrSG->ColCount;j++)
   tmpStr+=CurrSG->Cells[j][i]+";";
  tmp->Add(tmpStr);
 }
if(FileName==NULL)
{
Form1->SaveDialog1->FileName="";
Form1->SaveDialog1->DefaultExt="csv";
Form1->SaveDialog1->Filter="Таблица с разделителями (*.csv)| *.csv";
if (Form1->SaveDialog1->Execute())
 tmp->SaveToFile(Form1->SaveDialog1->FileName);
Form1->SaveDialog1->DefaultExt="els";
Form1->SaveDialog1->Filter="Электрическая схема (*.els)| *.els";
}
else
 tmp->SaveToFile(FileName);
}
__finally
{delete tmp;}


}


//---------------------------------------------------------------------------
void TFormR::OpenStringGrid(TStringGrid *CurrSG,int type)
{
AnsiString tmpStr;
Form1->OpenDialog1->DefaultExt="csv";
Form1->OpenDialog1->Filter="Таблица с разделителями (*.csv)| *.csv";

switch(type)
{
case 0:
 {
  tmpStr="Узлы";
  break;
 }
case 1:
 {
  tmpStr="Ветви";
  break;
 }
case 2:
 {
  tmpStr="Узлы-дополнительно";
  break;
 }
case 5:
 {
  tmpStr="Результаты-узлы";
  break;
 }
case 6:
 {
  tmpStr="Результаты-ветви";
  break;
 }

}

if(Form1->OpenDialog1->Execute())
{
TStrings *tmp=new TStringList();
TStrings *tmpCurr=new TStringList();
try
{
tmp->LoadFromFile(Form1->OpenDialog1->FileName);
if(tmp->Strings[0]==tmpStr)
{
CurrSG->RowCount=tmp->Strings[1].ToInt();
CurrSG->ColCount=tmp->Strings[2].ToInt();
for(int i=0;i<CurrSG->RowCount;i++)
 {
  tmpStr="";
  tmpCurr->Clear();
  for(int j=1;j<=tmp->Strings[i+3].Length();j++)
   {
    if(tmp->Strings[i+3][j]==';')
     {tmpCurr->Add(tmpStr);tmpStr="";}
    else
     tmpStr+=tmp->Strings[i+3][j];
   }
  CurrSG->Rows[i]=tmpCurr;
 }
newVr=true;
newSchema=true;
}
else
 ShowMessage("Несоответствие типа, таблица не открыта.");
}
__finally
{delete tmp;delete tmpCurr;}

}

Form1->OpenDialog1->DefaultExt="els";
Form1->OpenDialog1->Filter="Электрическая схема (*.els)| *.els";

}
//---------------------------------------------------------------------------
void __fastcall TFormR::Button1Click(TObject *Sender)
{
TStringGrid *tmpSG;
switch(PageControl->ActivePageIndex)
{
case 0:
 {
 tmpSG=StringGrid1;
 break;
 }
case 1:
 {
 tmpSG=StringGrid2;
 break;
 }
case 2:
 {
 tmpSG=StringGrid3;
 break;
 }
case 5:
 {
 tmpSG=StringGrid5;
 break;
 }
case 6:
 {
 tmpSG=StringGrid4;
 break;
 }
}
SaveStringGrid(tmpSG,PageControl->ActivePageIndex,NULL);
tmpSG=NULL;
delete tmpSG;
}
//---------------------------------------------------------------------------

void __fastcall TFormR::Button3Click(TObject *Sender)
{
TStringGrid *tmpSG;
switch(PageControl->ActivePageIndex)
{
case 0:
 {
 tmpSG=StringGrid1;
 break;
 }
case 1:
 {
 tmpSG=StringGrid2;
 break;
 }
case 2:
 {
 tmpSG=StringGrid3;
 break;
 }
case 5:
 {
 tmpSG=StringGrid5;
 break;
 }
case 6:
 {
 tmpSG=StringGrid4;
 break;
 }
}
OpenStringGrid(tmpSG,PageControl->ActivePageIndex);
tmpSG=NULL;
delete tmpSG;

}
//---------------------------------------------------------------------------

void __fastcall TFormR::Button5Click(TObject *Sender)
{
  int n;
  n=FindMaxNodeNumber(StringGrid1,0,0);
  n=FindMaxNodeNumber(StringGrid3,0,n);
  n=FindMaxNodeNumber(StringGrid2,0,n);
  n=FindMaxNodeNumber(StringGrid2,1,n);
  delete []VrTemp;
  VrTemp=new double[n*2];
  if(PredelRaschetMetod==0||PredelRaschetMetod==3)
   Form1->RaschetVr0(0,nVrTemp,VrTemp);
  else
   Form1->RaschetVr0(1,nVrTemp,VrTemp);
  double summProizv=0;
  double summSquarePredel=0;
  double summSquareRegim=0;
  int i;
  Form1->Memo1->Lines->Add("Векторы:");
  double NormV1=0,NormV2=0;
  for(i=0;i<nVrTemp;i++)
   {
    NormV1+=VrTemp[i]*VrTemp[i];
    NormV2+=Rasch->R[i]*Rasch->R[i];
   }
  NormV1=sqrt(NormV1);
  NormV2=sqrt(NormV2);
  Form1->Memo1->Lines->Add("Предельный:");
  for(i=0;i<nVrTemp;i++)
   {
    Form1->Memo1->Lines->Add(FloatToStrF(VrTemp[i]/NormV1,ffFixed,10,8)+"; "+FloatToStrF(Rasch->R[i]/NormV2,ffFixed,10,8)+";");
    Form1->Memo1->Lines->Add(FloatToStrF(Rasch->R[i]/NormV2,ffFixed,10,8));
    summProizv+=Rasch->R[i]*VrTemp[i];
    summSquarePredel+=Rasch->R[i]*Rasch->R[i];
    summSquareRegim+=VrTemp[i]*VrTemp[i];
   }
  Form1->Memo1->Lines->Add("S'");
  for(i=0;i<nVrTemp;i++)
   {
    Form1->Memo1->Lines->Add(FloatToStrF(VrTemp[i]/NormV1,ffFixed,10,8));
   }
  Form1->Memo1->Lines->Add("**********");
  double cosfi=summProizv/(sqrt(summSquarePredel)*sqrt(summSquareRegim));
  ShowMessage(acos(cosfi)*180./M_PI);
}
//---------------------------------------------------------------------------





void __fastcall TFormR::ButtStarPathClick(TObject *Sender)
{
FormSelectDir->Caption="Выберите папку";
if(FormSelectDir->ShowModal()==mrOk)
 StarFilesPath=FormSelectDir->DirName;

}
//---------------------------------------------------------------------------


void __fastcall TFormR::Button6Click(TObject *Sender)
{
int n,i;
n=FindMaxNodeNumber(StringGrid1,0,0);
n=FindMaxNodeNumber(StringGrid3,0,n);
n=FindMaxNodeNumber(StringGrid2,0,n);
n=FindMaxNodeNumber(StringGrid2,1,n);
delete []VrTemp;
VrTemp=new double[n*2];
if(PredelRaschetMetod==0||PredelRaschetMetod==3)
 Form1->RaschetTranspVr(0,0,nVrTemp,VrTemp);
else
 Form1->RaschetTranspVr(1,0,nVrTemp,VrTemp);

double NormV1=0;
for(i=0;i<nVrTemp;i++)
  NormV1+=VrTemp[i]*VrTemp[i];

NormV1=sqrt(NormV1);

Form1->Memo1->Lines->Add("Вектор R транспонированной матрицы Якоби");
for(i=0;i<nVrTemp;i++)
  Form1->Memo1->Lines->Add(FloatToStrF(VrTemp[i]/NormV1,ffFixed,10,8));

}
//---------------------------------------------------------------------------


void __fastcall TFormR::Button7Click(TObject *Sender)
{
int n,i;
n=FindMaxNodeNumber(StringGrid1,0,0);
n=FindMaxNodeNumber(StringGrid3,0,n);
n=FindMaxNodeNumber(StringGrid2,0,n);
n=FindMaxNodeNumber(StringGrid2,1,n);
delete []VrTemp;
VrTemp=new double[n*4];
if(PredelRaschetMetod==0||PredelRaschetMetod==3)
 Form1->RaschetTranspVr(0,1,nVrTemp,VrTemp);
else
 Form1->RaschetTranspVr(1,1,nVrTemp,VrTemp);

double NormV1=0;
for(i=0;i<nVrTemp;i++)
 NormV1+=VrTemp[i]*VrTemp[i];

NormV1=sqrt(NormV1);
Form1->Memo1->Lines->Add("Вектор R транспонированной матрицы Якоби");
for(i=0;i<nVrTemp;i++)
 Form1->Memo1->Lines->Add(FloatToStrF(VrTemp[i]/NormV1,ffFixed,10,8));

}
//---------------------------------------------------------------------------


void __fastcall TFormR::FormCreate(TObject *Sender)
{
ComboBox1->ItemIndex=0;
}
//---------------------------------------------------------------------------


void __fastcall TFormR::ButtPrmsFilePathClick(TObject *Sender)
{
Form1->SaveDialog1->FileName="";
Form1->SaveDialog1->DefaultExt="csv";
Form1->SaveDialog1->Filter="Таблица с разделителями (*.csv)| *.csv";
if (Form1->SaveDialog1->Execute())
 PrmsFileName=Form1->SaveDialog1->FileName;
Form1->SaveDialog1->DefaultExt="els";
Form1->SaveDialog1->Filter="Электрическая схема (*.els)| *.els";
}
//---------------------------------------------------------------------------



