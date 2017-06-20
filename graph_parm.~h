//---------------------------------------------------------------------------

#ifndef graph_parmH
#define graph_parmH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <CheckLst.hpp>
#include <ExtCtrls.hpp>
//#include "Values.h"
//#include "graph_cpp.h"
//---------------------------------------------------------------------------
class TForm5 : public TForm
{
__published:	// IDE-managed Components
        TGroupBox *GroupBox6;
        TPanel *Panel3;
        TComboBox *comboy1;
        TLabel *Label1;
        TPanel *Panel4;
        TLabel *Label2;
        TComboBox *comboy2;
        TListBox *ListGraph;
        TButton *Button1;
        TButton *Button2;
        TButton *Button3;
        TGroupBox *GroupBox5;
        TPanel *Panel5;
        TLabel *Label3;
        TComboBox *combox1;
        TPanel *Panel6;
        TLabel *Label4;
        TComboBox *combox2;
        TButton *Button4;
        TButton *Button5;
        TButton *Button6;
        TButton *Button7;
        void __fastcall FormShow(TObject *Sender);
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall combox1Change(TObject *Sender);
        void __fastcall comboy1Change(TObject *Sender);
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button4Click(TObject *Sender);
        void __fastcall Button3Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall Button5Click(TObject *Sender);
        void __fastcall Button7Click(TObject *Sender);
        void __fastcall Button6Click(TObject *Sender);

private:	// User declarations
public:		// User declarations
        bool IsMsw;
        int BusyGraph[10];
        int CreateGraph();
        bool ParamsFromList(AnsiString Str,int &Node1,int &Node2,int &type, int &Paramtype,AnsiString &AxeTitle);
        void AddSeriesToGraph(int Node1X,int Node2X,int Node1Y,int Node2Y,int typeX,int ParamtypeX,int typeY,int ParamtypeY,int GraphFormNo,int SerNo,int Count,AnsiString AxeXTitle,AnsiString AxeYTitle);
        void FillSeries(int type,int Paramtype,int ValsIndex,int Count,double *SerValue);
        TColor SerColor[15];
        bool combox1State;
        bool combox2State;
        bool comboy1State;
        bool comboy2State;
        TStringList *ListUzlParams;
        TStringList *ListLinkTrfParams;
        TStringList *ListLinkLineParams;
//        Values *Vals;
//        void *Vals;
        int ItersNumber;
        bool DeleteGraph(int GraphNo);

        __fastcall TForm5(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm5 *Form5;
//---------------------------------------------------------------------------
#endif
