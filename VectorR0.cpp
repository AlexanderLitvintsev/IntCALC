//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "VectorR0.h"
#include "RegimTable.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFormVectorVr *FormVectorVr;
//---------------------------------------------------------------------------
__fastcall TFormVectorVr::TFormVectorVr(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void TFormVectorVr::FillTable()
{
StringGrid1->RowCount=n+1;
StringGrid1->Cells[0][0]="�";
StringGrid1->Cells[1][0]="VRi";
for(int i=1;i<n+1;i++)
  {
   StringGrid1->Cells[0][i]=AnsiString(i);
   StringGrid1->Cells[1][i]=FloatToStrF(Vr[i-1],ffFixed,10,8);
  }

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void TFormVectorVr::SaveTable()
{
//n=StringGrid1->RowCount;
for(int i=1;i<n+1;i++)
 Vr[i-1]=FormR->MyStrToFloat(StringGrid1->Cells[1][i]);
}

//---------------------------------------------------------------------------
void __fastcall TFormVectorVr::Button1Click(TObject *Sender)
{
SaveTable();
ModalResult = mrOk;        
}
//---------------------------------------------------------------------------

void __fastcall TFormVectorVr::Button2Click(TObject *Sender)
{
ModalResult = mrCancel;
}
//---------------------------------------------------------------------------


void __fastcall TFormVectorVr::FormActivate(TObject *Sender)
{
FillTable();
}
//---------------------------------------------------------------------------

