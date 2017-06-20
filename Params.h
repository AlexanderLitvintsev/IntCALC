//---------------------------------------------------------------------------

#ifndef ParamsH
#define ParamsH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "graph_cpp.h"
#include <ExtCtrls.hpp>
//#include "Consts.h"

//---------------------------------------------------------------------------
//const int MaxElementParams=40;


class TForm2 : public TForm
{
__published:	// IDE-managed Components
        TPanel *Panel2;
        TPanel *Panel1;
        TButton *Button1;
        TButton *Button2;
        TLabel *Label1;
        TCheckBox *ParamsCheck;
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall FormActivate(TObject *Sender);
        void __fastcall FClick(TObject *Sender);
        void __fastcall DbClick(TObject *Sender);
        void __fastcall ParamsCheckClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
        int ParamsNumber;
//        AnsiString ParamNames[MaxElementParams];
        AnsiString *ParamNames;
        AnsiString ElementName;
        AnsiString ParFileName;
//        float Params[MaxElementParams];
        float *Params;
        TLabel **Label;
        TEdit **Edit;
        void ClearForm();
        TButton *FBut;
        TButton *DbBut;
        bool ParFromFile;
        bool ShowDB;
        bool UseParamsCheck;
        bool ParamsCheckEnabled;
        TPanel *CurrPanel;


        __fastcall TForm2(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm2 *Form2;
//---------------------------------------------------------------------------
#endif