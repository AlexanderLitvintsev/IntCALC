//---------------------------------------------------------------------------

#ifndef VectorR0H
#define VectorR0H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Grids.hpp>
//---------------------------------------------------------------------------
class TFormVectorVr : public TForm
{
__published:	// IDE-managed Components
        TStringGrid *StringGrid1;
        TButton *Button1;
        TButton *Button2;
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall FormActivate(TObject *Sender);
private:	// User declarations
public:		// User declarations
        int n;
        double *Vr;
        void FillTable();
        void SaveTable();
        __fastcall TFormVectorVr(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFormVectorVr *FormVectorVr;
//---------------------------------------------------------------------------
#endif
