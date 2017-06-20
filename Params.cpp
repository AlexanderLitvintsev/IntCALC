//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Params.h"
#include "elembase.h"
#include "graph_cpp.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

TForm2 *Form2;
const int MaxElementsParams=40;
//---------------------------------------------------------------------------
__fastcall TForm2::TForm2(TComponent* Owner)
        : TForm(Owner)
{
Params=new float [MaxElementsParams];
ParamNames=new AnsiString [MaxElementsParams];
}
//---------------------------------------------------------------------------
void __fastcall TForm2::Button1Click(TObject *Sender)
{
bool OkFlag=true;
AnsiString Field="";
int i;
for (i=0;i<ParamsNumber;i++)
  {
//   count=0;
   Field=Edit[i]->Text;
   for (int j=0;j<Edit[i]->Text.Length();j++)
     if(Edit[i]->Text[j+1]=='.')
       {
        Field[j+1]=',';
        Edit[i]->Text=Field;

       }
     if(Edit[i]->Text.Length()==0)
       Edit[i]->Text='0';
     OkFlag=true;
     try
       {
       Params[i]=StrToFloat(Edit[i]->Text);
       }
     catch (...)
       {
        ShowMessage("�������� ���������!");
        OkFlag=false;
        break;
       }
  }
if(OkFlag==true)
  {
   Form1->ImageForRightBtn->NumberOfParams=ParamsNumber;
   for(i=0;i<ParamsNumber;i++)
     {
      Form1->ImageForRightBtn->Params[i]=Params[i];
//      Form1->ImageForRightBtn->ParamNames[i]=ParamNames[i];
     }
   Form1->ImageForRightBtn->ParFromFile=ParFromFile;
   Form1->ImageForRightBtn->ParFileName=ParFileName;
   Form1->ImageForRightBtn->OtherParams=ParamsCheckEnabled;
   Form1->HasChanged=true;
   ClearForm();

   ModalResult=mrOk;
  }
else
  {
   Edit[i]->Text=FloatToStr(Params[i]);
   //   ShowMessage("������� ��������� ���������");
  }
}
//---------------------------------------------------------------------------
void __fastcall TForm2::Button2Click(TObject *Sender)
{
ClearForm();
ModalResult=mrCancel;

//Form2->Close();
}
//---------------------------------------------------------------------------
void __fastcall TForm2::FormActivate(TObject *Sender)
{
//int minheight=Button2->Top+Button2->Height+10;
int minheight=140;
//ParamsNumber=10;
int i;

if(ElementName!="����"&&ElementName!="�������"&&ElementName!="��������")
  ShowDB=true; //���������� ������ "����������"
else
  ShowDB=false;


if(ElementName.UpperCase().Pos("�������������")>0)
  UseParamsCheck=true; //���������� ����� ��� �������������� ����������
else
  UseParamsCheck=false;




this->Visible=false;
if(ParamsNumber>0)
{
Label=new TLabel*[ParamsNumber+1];
Edit=new TEdit*[ParamsNumber+1];

Panel1->AutoSize=false;
Panel2->AutoSize=false;
ParamsCheck->Visible=false;

int MaxWidth=0;
int MaxHeight;
bool OtherParamsVisible=false,flag=false;

int LabelTop=24;
CurrPanel=Panel1;
for (i=0;i<ParamsNumber;i++)
   {
    Label[i]=new TLabel(CurrPanel);
    Label[i]->Parent=CurrPanel;
//    Label[i]->Anchors=Label1->Anchors;
    Label[i]->Top=LabelTop;
    LabelTop+=24;
    Label[i]->Height=24;
//    Label[i]->Width=130;

    Label[i]->AutoSize=true;
    Label[i]->Left=1;
    //    Label[i]->Caption=i;
    Label[i]->Caption=ParamNames[i];
    Label[i]->Alignment=taRightJustify;
    if(Label[i]->Width>MaxWidth)
       MaxWidth=Label[i]->Width;
    if(UseParamsCheck==true&&ParamNames[i]=="n, ��")
      {
       CurrPanel=Panel2;
       LabelTop=ParamsCheck->Top+ParamsCheck->Height+5;

      }
   }

CurrPanel=Panel1;
LabelTop=24;
for(i=0;i<ParamsNumber;i++)
   {
    if(flag==true)
      OtherParamsVisible=true;
    Label[i]->AutoSize=false;
    Label[i]->Width=MaxWidth;
    Label[i]->Left=1;
    Edit[i]=new TEdit(CurrPanel);
    Edit[i]->Parent=CurrPanel;
//    Edit[i]->Anchors=Label1->Anchors;
    Edit[i]->Left=Label[i]->Left+Label[i]->Width+5;;
    Edit[i]->Top=LabelTop;


    LabelTop+=24;
    Edit[i]->Height=24;
    Edit[i]->Width=80;
    Edit[i]->Text=FloatToStrF(Params[i],ffFixed,7,7);
    if(UseParamsCheck==true&&ParamNames[i]=="n, ��")
      {
       LabelTop=24;
       CurrPanel=Panel2;
       LabelTop=ParamsCheck->Top+ParamsCheck->Height+5;
       flag=true;
      }
   }
if(ShowDB==true)
  {
   DbBut=new TButton(Form2);
   DbBut->Parent=Form2;
   DbBut->Caption="����������";
   DbBut->OnClick=DbClick;
  }
if(OtherParamsVisible==true)
  ParamsCheck->Visible=true;
if(ParamsCheckEnabled==true)
 ParamsCheck->Checked=true;
else
 ParamsCheck->Checked=false;
for(i=0;i<ParamsNumber;i++)
 {
  if(Edit[i]->Parent==Panel2)
    Edit[i]->Enabled=ParamsCheckEnabled;
  else
    Edit[i]->Enabled=!ParamsCheckEnabled;
 }



Panel1->AutoSize=true;
Panel2->AutoSize=true;
Panel1->Left=0;
Panel1->Top=0;
Panel2->Top=0;
Panel2->Left=Panel1->Left+Panel1->Width;
MaxHeight=0;
MaxHeight=max(Panel1->Height,Panel2->Height);
Panel1->AutoSize=false;
Panel2->AutoSize=false;
Panel1->Height=MaxHeight;
Panel2->Height=MaxHeight;

Button1->Top=5;
Button2->Top=Button1->Top+Button1->Height+10;
Button1->Left=Panel2->Left+Panel2->Width+5;
Button2->Left=Button1->Left;
MaxHeight=max(MaxHeight,Button2->Top+Button2->Height+5);
MaxWidth=Button1->Left+Button1->Width+5;

if(ShowDB==true)
 {
  DbBut->Left=Button1->Left;
  DbBut->Top=Button2->Top+Button2->Height+10;
  MaxWidth=max(MaxWidth,DbBut->Left+DbBut->Width+5);
  MaxHeight=max(MaxHeight,DbBut->Top+DbBut->Height+5);
 }



this->Width=MaxWidth;
this->Height=MaxHeight+20;
this->Visible=true;
}
}

//---------------------------------------------------------------------------
void TForm2::ClearForm()
{
int i;
if(ParamsNumber>0)
{
//if(ParFromFile==true)
  {
  for (i=0;i<ParamsNumber;i++)
    {
     delete Label[i];
     delete Edit[i];

    }
/*
   if(ElementName!="����")
     {
      delete Label[i];
      delete Edit[i];
      delete FBut;
      delete DbBut;
     }
*/
   if(ShowDB==true)
     delete DbBut;
   }
/*
else
  for (int i=0;i<ParamsNumber;i++)
    {
     delete Label[i];
     delete Edit[i];
    }
*/
delete []Label;
delete []Edit;
Label=NULL;
Edit=NULL;
ParFromFile=false;
}
}

//---------------------------------------------------------------------------
void __fastcall TForm2::FClick(TObject *Sender)
{

AnsiString DefExt=Form1->OpenDialog1->DefaultExt;
AnsiString Filter=Form1->OpenDialog1->Filter;
Form1->OpenDialog1->DefaultExt="txt";
Form1->OpenDialog1->Filter="��������� ���� � ����������� (*.txt)| *.txt";
Form1->OpenDialog1->FileName="";
Form1->OpenDialog1->Execute();
if(Form1->OpenDialog1->FileName!="")
  {
   ParFileName=Form1->OpenDialog1->FileName;
   ParFromFile=true;
   Edit[ParamsNumber]->Text=ParFileName;
   Edit[ParamsNumber]->Enabled=true;
   Edit[ParamsNumber]->ReadOnly=true;
   Form1->OpenDialog1->FileName="";
  }
else
  {
   ParFileName="";
   ParFromFile=false;
   Edit[ParamsNumber]->Text="�� ������";
   Edit[ParamsNumber]->Enabled=false;
   Edit[ParamsNumber]->ReadOnly=true;
  }

//Form1->OpenDialog1-
Form1->OpenDialog1->DefaultExt=DefExt;
Form1->OpenDialog1->Filter=Filter;

}
//---------------------------------------------------------------------------
void __fastcall TForm2::DbClick(TObject *Sender)
{


Form6->type=Form1->ImageForRightBtn->type;
Form6->ParamsNumber=ParamsNumber-1;
for(int i=0;i<ParamsNumber-1;i++)
  Form6->ParamNames[i]=ParamNames[i];
Form6->ShowModal();
if(Form6->ModalResult==mrOk)
  {
   for(int i=0;i<ParamsNumber;i++)
     Edit[i]->Text=FloatToStrF(Params[i],ffFixed,7,7);
  }
}

void __fastcall TForm2::ParamsCheckClick(TObject *Sender)
{
if(ParamsCheck->Checked==true)
 {
  ParamsCheckEnabled=true;
  for(int i=0;i<ParamsNumber;i++)
   {
    if(Edit[i]->Parent==Panel2)
     Edit[i]->Enabled=true;
    else
     Edit[i]->Enabled=false;
   }
 }
else
 {
  ParamsCheckEnabled=false;
  for(int i=0;i<ParamsNumber;i++)
   {
    if(Edit[i]->Parent==Panel2)
     Edit[i]->Enabled=false;
    else
     Edit[i]->Enabled=true;
   }
 }

}
//---------------------------------------------------------------------------


