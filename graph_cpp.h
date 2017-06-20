//---------------------------------------------------------------------------

#ifndef graph_cppH
#define graph_cppH
//---------------------------------------------------------------------------
//#include "all.h"
#include <ActnList.hpp>
#include <Classes.hpp>
#include <Controls.hpp>
#include <Dialogs.hpp>
#include <ImgList.hpp>
#include <Menus.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Graphics.hpp>
#include <ComCtrls.hpp>
#include <ToolWin.hpp>
#include <math.h>
#include <complex.h>
#include <process.h>
#include "ap.h"
#include <iostream>
#include <complex>
//#include "Values.h"

//#include "consts.h"
//#include "Rasch.h"

//#include <stdio.h>
//---------------------------------------------------------------------------
using namespace std;
static int c__10 = 10;
const int UP=0;
const int DOWN=1;
const int LEFT=2;
const int RIGHT=3;
const int MaxNodesInElem=6;
const int MaxNumberOfElem=200;
//
const int MaxNodesNumber=200;
const int MaxLinksNumber=MaxNodesNumber*10;
const int MaxCycles=50;
//const int MaxElementParams=40;
const int MaxElementParams=40;
const int VerNo=1;

/*
const int UP=0;
const int DOWN=1;
const int LEFT=2;
const int RIGHT=3;
const int MaxNodesInElem=6;
const int MaxNumberOfElem=200;
//
const int MaxNodesNumber=1;
const int MaxLinksNumber=MaxNodesNumber*10;
const int MaxCycles=50;
//const int MaxElementParams=40;
const int MaxElementParams=40;
*/
//int NumberOfIters=500;
//double f0=50;


//---------------------------------------------------------------------------
struct Stru
   {
    int i,j,type,u4,u5,u6;
   };
//---------------------------------------------------------------------------
struct Pqu
   {
    int NodeNumber,Hui,TipN;
    double P,Pg,dPg,dPn;
    double Q,Qg,dQn;
    double Un;
    double a0,a1,a2,b0,b1,b2;
    double Qmin;
    double Qmax;
   };
//---------------------------------------------------------------------------
struct Uzl
   {
    int NodeNumber;
    int Hui;
    int a[3],b[3];
   };
//---------------------------------------------------------------------------
struct Str
   {
	int NodeiNumber,
        NodejNumber,
        Mswij,Mhui,Mhuj;
    double cktr,cktm,yr,ym;
    bool flag;
   };
//---------------------------------------------------------------------------
struct Lin
   {
    int NodeiNumber,
        NodejNumber,TipN;
    double r0,x0,b0,Length,n;
   };
//---------------------------------------------------------------------------
struct Trf
   {
    int NodeiNumber,
        NodejNumber,TipN;
    double ktr,un,sn,pk,uk,px,ix,xt,n;
    double ktr_i,r,x,g,b;
    bool paramtype;
   };
//---------------------------------------------------------------------------
struct Kol
   {
       int NumberOfTrf,NumberOfLines,NumberOfNodes;
	   double Ubaz,Precision;
   };
//---------------------------------------------------------------------------
struct ShReact
   {
       int NodeNumber;
       double G,B;
   };
//---------------------------------------------------------------------------


class MyImage : public TImage
{

public:
        MyImage(TForm*Form);
        int NumberOfPoints;
        int ElementNumber;
        TRect RectCoord[MaxNodesInElem];
        int Direction[MaxNodesInElem];
        bool PointIsLinked[MaxNodesInElem];
        bool IsClicked;
        int LinkObjectNumber[MaxNodesInElem];
        int ObjectNumber;
        double Params[MaxElementParams];
        AnsiString ParamNames[MaxElementParams];
		int NumberOfParams;



        int MainRectNo;
        int BalanceRectNo;
        int type;
        AnsiString ElementName;


        bool PointIsForwarded[MaxNodesInElem];
        int PointNumber[MaxNodesInElem];
        int RealPointNumber[MaxNodesInElem];
        int InPointNo[2];
        int InIndex;
        bool ParFromFile;
        bool OtherParams;
        AnsiString ParFileName;

};

//---------------------------------------------------------------------------

class TemplateObject
{
public:
	   AnsiString ElementName;
       int NumberOfPoints;
       TRect RectCoord[MaxNodesInElem];
       int Direction[MaxNodesInElem];
       int NumberOfParams;
       int type;
       AnsiString ParamsNames[MaxElementParams];
       AnsiString PictureFileName;
       AnsiString IcoPictureFileName;

};
//---------------------------------------------------------------------------
class tmpLineCoords
{
public:
        POINT Start;
        POINT End;
        int Obj1Number;
        int Node1Number;
        int Obj2Number;
        int Node2Number;
};
//---------------------------------------------------------------------------
class LineCoords:public tmpLineCoords
{
public:
		static int LineNumber;
        LineCoords();
        ~LineCoords();
};
int LineCoords::LineNumber;
//---------------------------------------------------------------------------
class NodesNum
{
public:
        NodesNum();
        ~NodesNum();
        int ObjNumber;
        int NumberOfNodes;
        bool NodeIsForwarded[MaxNodesInElem];
        int NodeNumber[MaxNodesInElem];
        int InNodeNo;
        int OutNode1No;
        int OutNode2No;
        int Number;
        //        static
        static int Quantity;
};
int NodesNum::Quantity=-1;
//---------------------------------------------------------------------------
/*
struct Prms
{
int Node1,Node2;
complex<double> SFlow1,SFlow2,DS,DSTM,DSTC,I;
double Ql;
double DU;
bool IsLine,IsTrf;
};
*/
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
class ParamForCycle
{
public:
       ParamForCycle();
       ~ParamForCycle();
       static int Quantity;
       int ElementNo;
       int NumberOfIters;
       int NumberOfParams;
       double Params[MaxElementParams][50];
       AnsiString FileName;
       bool Ok;

       void Add(int ElNo,AnsiString FName);
       bool ReadFromFile();

};
int ParamForCycle::Quantity=0;
//---------------------------------------------------------------------------

class MyForm : public TForm
{
public:
      __fastcall MyForm(TComponent* Owner);
};



class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TActionList *ActionList1;
        TAction *Obj1;
        TImageList *ImageList1;
        TAction *Obj2;
        TPopupMenu *PopMenu;
        TMenuItem *N1;
        TMenuItem *N2;
        TMenuItem *N3;
        TMainMenu *MainMenu1;
        TMenuItem *N4;
		TMenuItem *N5;
        TMenuItem *N6;
        TMenuItem *N7;
        TMenuItem *N9;
        TMenuItem *N8;
        TSaveDialog *SaveDialog1;
        TMenuItem *N12;
        TMenuItem *N13;
        TMenuItem *N14;
        TMenuItem *N15;
        TOpenDialog *OpenDialog1;
        TMenuItem *N16;
        TMenuItem *N22;
        TMenuItem *N23;
        TMenuItem *N25;
        TMenuItem *N32;
        TMenuItem *N33;
        TMenuItem *N34;
        TMenuItem *N35;
        TMenuItem *N39;
        TMenuItem *N40;
        TMenuItem *N11;
        TMenuItem *N10;
        TMemo *Memo1;
        void __fastcall Image1MouseDown(TObject *Sender,
          TMouseButton Button, TShiftState Shift, int X, int Y);
		void __fastcall Image1MouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
        void __fastcall Image1MouseUp(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
        void __fastcall Obj1Execute(TObject *Sender);
        void __fastcall FormClick(TObject *Sender);
        void __fastcall FormMouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
        void __fastcall Image1Click(TObject *Sender);
        void __fastcall N2Click(TObject *Sender);
        void __fastcall N3Click(TObject *Sender);
        void __fastcall N1Click(TObject *Sender);
        void __fastcall N5Click(TObject *Sender);
        void __fastcall N6Click(TObject *Sender);
        void __fastcall N7Click(TObject *Sender);
        void __fastcall MainMenu1Change(TObject *Sender, TMenuItem *Source,
          bool Rebuild);
        void __fastcall N8Click(TObject *Sender);
        void __fastcall N9Click(TObject *Sender);
        void __fastcall N11Click(TObject *Sender);
        void __fastcall N13Click(TObject *Sender);
        void __fastcall N14Click(TObject *Sender);
        void __fastcall N15Click(TObject *Sender);
        void __fastcall N16Click(TObject *Sender);
        void __fastcall FormPaint(TObject *Sender);
        void __fastcall N19Click(TObject *Sender);
		void __fastcall N20Click(TObject *Sender);
        void __fastcall N21Click(TObject *Sender);
        void __fastcall N22Click(TObject *Sender);
        void __fastcall N24Click(TObject *Sender);
        void __fastcall N25Click(TObject *Sender);
        void __fastcall N32Click(TObject *Sender);
        void __fastcall N35Click(TObject *Sender);
        void __fastcall N34Click(TObject *Sender);
        void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
        void __fastcall N39Click(TObject *Sender);
        void __fastcall N40Click(TObject *Sender);
        void __fastcall N10Click(TObject *Sender);
//        void __fastcall Button2Click(TObject *Sender);
private:	// User declarations
public:		// User declarations

        //���������� �������� �� �����
        void AddElement(int ElemNo,AnsiString ElementName,int &NumberOfElem,TPicture *Pic,POINT point,int NumberOfPoints,TRect *Rects, int NumberOfParams,AnsiString *ParamNames,double *Params, int *Dirs);
        //�������� ��������
        void DeleteElement(int ElemNo,int &NumberOfElem);
        //��������� ����� � �������� ��������
        bool PointInRect(int x,int y,TRect Rect);
        //������ �����������
        void Redraw();
        //����������� �����
        void RedrawWithoutRefresh();
		//��������� ��� ����� ������� ������
        void Link2Points(POINT P1,POINT P2,int Dir1,int Dir2,int MaxWidth,int MaxHeight,bool is1Sheena,bool is2Sheena);
        //��������, ��������������� ����������� Dir
        int ReverseDir(int Dir);
        //��������� ����� ��������
        void BreakLinks();
        //������� �� ����� ������ ����� ���� ������� �������
        void PrintNumber(MyImage *Object,TRect Rect,int Number,int RectNumber);
        //������������ ��������� �����
        void Numerate();
        //���� ���� ����� ������ ���� ������� � ������ ����� ��. �������
        bool FindLink(int ObjStart,int NodeStart,int &ObjFinish,int &NodeFinish);
        //��������� ����� � ����
        void Save(AnsiString FName);
        //��������� ����� �� �����
        void Open(AnsiString FName);
        //������� �����
        void CleanSchema();

        void MenuOpen();
        void MenuSave();
        void MenuSaveAs();
        void MenuClose();
        void MenuNew();
        bool MenuCloseSchema();
        void MenuExit();


        //���������� � �������

        void Prepare();
        //��� ���� �������
        int PointType(int CurrObjNumber,int PointNumber);
        //������ ������� Z
        bool Raschet(bool showresult);

        //������� ������ � ����������� �����.
        bool RaschetIskl(bool showresult);

        //������ ������� �������
        bool RaschetJ(bool showresult,bool &continueJ);


        bool RaschetJEkviv(bool showresult,bool &continueJ);

        bool RaschetJProba(bool showresult,bool &continueJ);

        bool RaschetJX2(bool showresult,bool &continueJ);
//        bool RaschetJKn(bool showresult,bool &continueJ,bool askKn,bool &RegimExt,bool ogran);
        bool RaschetJKn(double Kn,bool UseUnom);
        bool RaschetJX22(bool showresult,bool &continueJ,bool WritePrmsToFile,AnsiString PrmsFileName);
        bool RaschetJX2T(bool showresult,bool &continueJ,bool WritePrmsToFile,AnsiString PrmsFileName);
		bool RaschetJSimpleT(bool showresult,bool &continueJ);
        bool RaschetJT(double T,bool UseUnom);
        void RaschetVr0(bool type,int &nVr,double *Vr);
        void RaschetSingular(bool type,int &nVr,double *Vr);
		void RaschetTranspVr(bool type,bool MatrixType,int &nVr,double *Vr);
		void IntervalRaschet();
        //������ ������� �������������
        void CalculateMain(AnsiString FName);
        //����� ����������� (���� � �����. ����)
        void ResultsOut(int kit,bool Regim,bool Showtxt);
        //������ ������ ������, � �.�. � ����������
        bool UstRegimRaschet(int type,bool BkUse,bool CurrReg,bool Showtxt,double *Vr0,int nVr0, bool WritePrmsToFile,AnsiString PrmsFileName);

        void FillRealNodesNumbers();

        //������� ����������
        int MyRound(double val);

        void HardRaschet(int type);
        bool RaschetSeries(int type,double T,double Kn,int &n,bool Showtxt);

        void CopyRasch();
        int CreateRasch();
        int CreateRasch63Proba();
        void SaveRasch();
        void CopyResultsParams(int sourceind, int destind);
//        RECT RectToLink(int ObjNumber,int NodeNumber);
        //int ObjectNumberToRedraw;

        bool NewFile;//���� ���� ��� �� ����������
        bool HasChanged;//���� ����� ��� ����� ���������
        AnsiString FileName;//��� ����� ������� ����� (els)
        int NumberOfElem;//����� ��������� �� �����
        int MouseX,MouseY;//���������� ������� ��� mouseclick �� ��������
        int CurrMouseX,CurrMouseY;//������� ���������� ������� ��� ���������� ��. ���� ��� ��������������
        POINT OldMousePos;//������ ������� �������� (�� ������ ��������������)
        TPicture *PicToCreate;//�������� ��������, ������� ����� ���������� �� �����
        bool CreatePicFlag;//������� � ���, ��� ���� ������ ������ �� ������ ���������, �.�. ��� �������� �� ����� � ���� ����� �������� ���. ����-�
        bool MouseBtnLeftIsDown;//����� ������ ���� ������� (�� ��������)

        bool OnObject;//��� �������� ����
        bool ObjectIsClicked;//���� �������� � �������� ����� - �.�. ���� ����� �����
//        bool ObjectIsSelected;

        MyImage *tmpImage;//��������� ������
        MyImage *ClickedObject;//�������, �� �������� ����� �����
        MyImage *SelectedObject;//��������� �������
        MyImage *ImageForRightBtn;//�������, �� ������� �������� ������ �������
        LineCoords *LineArray;//������ - ������ � ������������ ����� � ��������� �� � ��������
        tmpLineCoords *tmpLineArray;//�� ��, ������ ���������
        TImage **Picture;//������ ������� �������� ��������� (�� ����� ���. ���������)
        MyImage **Element;//���� ��������

        NodesNum *NodeNum;//�������� ������ ����� ��������, ����� � �.�.

        TStringList *List;//����� �����,
        AnsiString ElementName;//��� ��������, ���. ����� ��������� �� �����
        int MainElementNo;//����� ��������� ��������
        int BalanceElementNo;//����� ����������� ��������
        int LinesNumber;//������� ����� ����� �����
        int NumberOfNodes;//����� �����, ������� ����� ��������� ����� �������
//        int NumberOfStr;
        int type;//��� �������� ��������
        TRect NodesCoord[MaxNodesInElem];//���������� ��� �����
        int NumberOfElemParams;//����� ��� ����������
        double Params[MaxElementParams];//���� ���������
        AnsiString ParamNames[MaxElementParams];//�������� ����������
        int Dirs[MaxNodesInElem];//�����������
        int NumberOfElements;//����� ��������� �� ������ ���������
        int ElementNo;//����� ��������
        AnsiString UstOutDir;//���������� ��� ������ ������ ���-��� �������
        int CurrIt;//������� �������� ����� ����������
        int MaxNumberOfIt;//������������ ����� �������� ��� ����� ���������� (���������� �� ���� ���������)
        int HardItersNumber;//����� �������� � ����� ����������
        double precision;//�������� �������
        bool IsTxtsOut;//��������� �� ��������� ����� � �������������� �������
        //������� �����, ������������ � �������
//        complex<double> Jac[MaxNodesNumber][MaxNodesNumber];
		complex<double> **Jac;



        int NumberOfIters;
        double f0;
        bool UseQminQmax;
        bool LimitNebals;
        double dUmax;
        double dfmax;
        double dKnmax;
        double dTmax;
        double dVrimax;
        double Kn0;
        double T0;
        int Iters;

//        int RectNo;
        TImage **Image;//�������� ������������ �������� ���������
        TToolButton **ToolBtn;//������ ������ ���������
        TemplateObject *TemplObj;//��������� ������

//����� ����������� ��������� �����
        struct Pqu *PquList;
        struct Trf *TrfList;
        struct Uzl *UzlList;
        struct Str *StrList;
        struct Lin *LinList;
        struct Kol KolList;
        struct Stru *StruList;
        struct ShReact *ShReactList;


        int MaxNodeNumber;//����. ����� ����
        bool kzMode;//����� ������� ��
        int NumberOfTrf,NumberOfLines;//����� ������, ����� ���.
//        Reg Rasch;
//        Reg *Rasch;
        ParamForCycle CyclePar[MaxNumberOfElem];//���������, �����. ��� ������������ ������� ���������� ������
        AnsiString RezFileName;//�������, ���������� notepad.exe  � ��� ����� ��� ������� � ��������
//        Values Vals[MaxCycles];//������������ ���������
//        Values *Vals;
//������������ ��������� ������
        complex<double> MaxGradS[MaxNodesNumber];//��������
        complex<double> MaxGradSTrC[MaxNodesNumber];//�������� � ����� ��-���
        complex<double> MinGradS[MaxNodesNumber];//���. ����.
        complex<double> MinGradSTrC[MaxNodesNumber];
        int Node1,Node2;//���� �����, ���. � ������ ���. ���������

        bool RealNodes;
        double Mk;
        __fastcall TForm1(TComponent* Owner);

};


//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;

//---------------------------------------------------------------------------
#endif
