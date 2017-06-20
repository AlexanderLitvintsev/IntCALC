//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "elembase.h"
//#include "graph_cpp.h"
#include "Params.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm6 *Form6;
//---------------------------------------------------------------------------
__fastcall TForm6::TForm6(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm6::FormClose(TObject *Sender, TCloseAction &Action)
{

Table1->Active=false;
}
//---------------------------------------------------------------------------
void __fastcall TForm6::FormCreate(TObject *Sender)
{
type=-1;        
}
//---------------------------------------------------------------------------


void __fastcall TForm6::FormShow(TObject *Sender)
{
int i;
if(type!=-1)
{
 Table1->Filter="ElemType="+AnsiString(type);
 Table1->Filtered=true;
}
type=-1;

for(i=0;i<DBGrid1->Columns->Count-1;i++)
 {
  if(i<ParamsNumber)
   {
    DBGrid1->Columns->Items[i+1]->Visible=true;
    DBGrid1->Columns->Items[i+1]->Title->Caption=ParamNames[i];
    DBGrid1->Columns->Items[i+1]->Width=60;
   }
  else
    DBGrid1->Columns->Items[i+1]->Visible=false;
 }

Table1->Active=true;
//Table1->Edit();

}
//---------------------------------------------------------------------------

void __fastcall TForm6::Button2Click(TObject *Sender)
{
ModalResult=mrCancel;


}
//---------------------------------------------------------------------------


void __fastcall TForm6::Button1Click(TObject *Sender)
{
float tmp;
if(Table1->RecordCount>0)
  for(int i=0;i<ParamsNumber;i++)
    Form2->Params[i]=DBGrid1->Fields[i+1]->AsFloat;
ModalResult=mrOk;

//ShowMessage(tmp);

}
//---------------------------------------------------------------------------



