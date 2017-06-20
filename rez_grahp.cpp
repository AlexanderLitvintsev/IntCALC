//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "rez_grahp.h"
#include "graph_parm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm4 *Form4;

//---------------------------------------------------------------------------
__fastcall TForm4::TForm4(TComponent* Owner)
        : TForm(Owner)
{
Ser=new TLineSeries*[15];
}
//---------------------------------------------------------------------------

void __fastcall TForm4::FormClose(TObject *Sender, TCloseAction &Action)
{
Form5->DeleteGraph(FormNo);
}
//---------------------------------------------------------------------------


void __fastcall TForm4::FormCreate(TObject *Sender)
{
Chart1->Left=0;
Chart1->Top=0;
Chart1->Width=this->Width;
Chart1->Height=this->Height;
}
//---------------------------------------------------------------------------

void __fastcall TForm4::FormCanResize(TObject *Sender, int &NewWidth,
      int &NewHeight, bool &Resize)
{
Chart1->Left=0;
Chart1->Top=0;
Chart1->Width=this->Width;
Chart1->Height=this->Height-20;
}
//---------------------------------------------------------------------------






