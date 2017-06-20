//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "graph_parm.h"
#include "graph_cpp.h"
#include "rez_grahp.h"
#include "Values.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm5 *Form5;
TForm4 **GraphForm;
//---------------------------------------------------------------------------
__fastcall TForm5::TForm5(TComponent* Owner)
        : TForm(Owner)
{
for(int i=0;i<10;i++)
  BusyGraph[i]=false;
GraphForm=new TForm4*[10];
SerColor[0]=clRed;
SerColor[1]=clBlack;
SerColor[2]=clGreen;
SerColor[3]=clYellow;
SerColor[4]=clBlue;
SerColor[5]=clWhite;
SerColor[6]=clPurple;
SerColor[7]=clDkGray;
SerColor[8]=clGray;
SerColor[9]=clTeal;
SerColor[10]=clAqua;
SerColor[11]=clOlive;
SerColor[12]=clSilver;
SerColor[13]=clNavy;
SerColor[14]=clFuchsia;

}
//---------------------------------------------------------------------------
int TForm5::CreateGraph()
{
int GraphNo=-1;
for(int i=0;i<10;i++)
  if(BusyGraph[i]==false)
   {
    GraphForm[i]=new TForm4(Application);
    GraphForm[i]->Visible=true;
    GraphForm[i]->FormNo=i;
    BusyGraph[i]=true;
    GraphNo=i;
    break;
   }
return GraphNo;
}
//---------------------------------------------------------------------------
bool TForm5::DeleteGraph(int GraphNo)
{
bool Ok=false;
if(BusyGraph[GraphNo]==true)
{
 delete GraphForm[GraphNo];
 BusyGraph[GraphNo]=false;
 Ok=true;
}
return Ok;
}
//---------------------------------------------------------------------------

void __fastcall TForm5::FormShow(TObject *Sender)
{
int i;
combox1->Items->Clear();
//if(Form1->RealNodes==true)
//combox1->Items->Add("Номер итерации утяжеления");
for(i=0;i<ValsGraph[0].NodesNumber;i++)
   combox1->Items->Add(AnsiString(ValsGraph[0].RealNode[i]));
for(i=0;i<ValsGraph[0].mswNumber;i++)
  combox1->Items->Add(AnsiString(ValsGraph[0].RealNode[ValsGraph[0].Node1[i]-1])+"-"+AnsiString(ValsGraph[0].RealNode[ValsGraph[0].Node2[i]-1]));
combox1->Text=combox1->Items->Strings[0];
combox2->Items=ListUzlParams;
combox2->Text=combox2->Items->Strings[0];
comboy1->Items=combox1->Items;
comboy2->Items=combox2->Items;
comboy2->Items->Delete(0);
comboy1->Text=comboy1->Items->Strings[0];
comboy2->Text=comboy2->Items->Strings[0];
ListGraph->Clear();

}
//---------------------------------------------------------------------------

void __fastcall TForm5::FormCreate(TObject *Sender)
{
ListUzlParams=new TStringList;
try
{
ListUzlParams->Clear();
ListUzlParams->Add("Номер итерации утяжеления");
ListUzlParams->Add("Модуль напряжения");
ListUzlParams->Add("Фаза напряжения");
ListUzlParams->Add("Активная составляющая напряжения");
ListUzlParams->Add("Реактивная составляющая напряжения");
}
__finally
{
}
ListLinkTrfParams=new TStringList;
try
{
ListLinkTrfParams->Clear();
ListLinkTrfParams->Add("Номер итерации утяжеления");
ListLinkTrfParams->Add("Модуль тока");
ListLinkTrfParams->Add("Фаза тока");
ListLinkTrfParams->Add("Активная составляющая тока");
ListLinkTrfParams->Add("Реактивная составляющая тока");
ListLinkTrfParams->Add("Активные потери в меди трансформатора");
ListLinkTrfParams->Add("Реактивные потери в меди трансформатора");
ListLinkTrfParams->Add("Активные потери в стали трансформатора");
ListLinkTrfParams->Add("Реактивные потери в стали трансформатора");
}
__finally
{
}
ListLinkLineParams=new TStringList;
try
{
ListLinkLineParams->Clear();
ListLinkLineParams->Add("Номер итерации утяжеления");
ListLinkLineParams->Add("Модуль тока");
ListLinkLineParams->Add("Фаза тока");
ListLinkLineParams->Add("Активная составляющая тока");
ListLinkLineParams->Add("Реактивная составляющая тока");
ListLinkLineParams->Add("Модуль потерь напряжения");
ListLinkLineParams->Add("Активные потери мощности в линии");
ListLinkLineParams->Add("Реактивные потери мощности в линии");
}
__finally
{
}
combox2->Items=ListUzlParams;
comboy2->Items=ListUzlParams;
}
//---------------------------------------------------------------------------

void __fastcall TForm5::combox1Change(TObject *Sender)
{
if(combox1->Text.Pos("-")>0)
  {
   if(ValsGraph[0].IsLine[combox1->ItemIndex-ValsGraph[0].NodesNumber]==true)
     combox2->Items=ListLinkLineParams;
   else
     combox2->Items=ListLinkTrfParams;
  }
else
  combox2->Items=ListUzlParams;
combox2->Text=combox2->Items->Strings[0];
}
//---------------------------------------------------------------------------

void __fastcall TForm5::comboy1Change(TObject *Sender)
{
static bool Link;
static bool prevLine;
bool currLine=false,currLink=false;

if(comboy1->Text.Pos("-")>0)
  {
   currLink=true;
   if(ValsGraph[0].IsLine[comboy1->ItemIndex-ValsGraph[0].NodesNumber]==true)
     {comboy2->Items=ListLinkLineParams;currLine=true;}
   else
     {comboy2->Items=ListLinkTrfParams;;currLine=false;}
  }
else
  {
   currLink=false;
   comboy2->Items=ListUzlParams;
  }

comboy2->Items->Delete(0);
if((prevLine!=currLine)||(Link!=currLink))
{
comboy2->Text=comboy2->Items->Strings[0];
}
Link=currLink;
prevLine=currLine;
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button1Click(TObject *Sender)
{
AnsiString GraphName="";
if(comboy1->Text.Pos("-")>0)
 GraphName+="Ветвь "+comboy1->Text+": ";
else
 GraphName+="Узел "+comboy1->Text+": ";
GraphName+=comboy2->Text;
ListGraph->Items->Add(GraphName);
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button4Click(TObject *Sender)
{

ListGraph->Clear();
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button3Click(TObject *Sender)
{
this->Close();
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button2Click(TObject *Sender)
{
AnsiString AxeXStr="";
AnsiString AxeXTitle,AxeYTitle;
int i,Node1X,Node2X,typeX,ParamtypeX,Node1Y,Node2Y,typeY,ParamtypeY;
int k=CreateGraph();
int SerNo=0;
if(k==-1)
 {
  ShowMessage("Слишком много открытых графиков!");
  goto Label_TooManyGraphs;
 }

AxeXStr="";//combox1->Text+;
if(combox1->Text.Pos("-")>0)
 AxeXStr+="Ветвь "+combox1->Text+": ";
else
 AxeXStr+="Узел "+combox1->Text+": ";
AxeXStr+=combox2->Text;

if(combox2->Text!="Номер итерации утяжеления")
 ParamsFromList(AxeXStr,Node1X,Node2X,typeX,ParamtypeX,AxeXTitle);
else
 {
  Node1X=Node2X=-1;typeX=4,ParamtypeX=-1;
  AxeXTitle="Номер итерации утяжеления";
 }
for(i=0;i<ListGraph->Items->Count;i++)
{
if(i==15)
 {ShowMessage("Больше 15 графиков - многовато!");break;}
if(ParamsFromList(ListGraph->Items->Strings[i],Node1Y,Node2Y,typeY,ParamtypeY,AxeYTitle))
 {
//  AddSeriesToGraph(Node1X,Node2X,Node1Y,Node2Y,typeX,ParamtypeX,typeY,ParamtypeY,k,SerNo,Form1->HardItersNumber,AxeXTitle,AxeYTitle);
  AddSeriesToGraph(Node1X,Node2X,Node1Y,Node2Y,typeX,ParamtypeX,typeY,ParamtypeY,k,SerNo,ItersNumber,AxeXTitle,AxeYTitle);
  SerNo++;
 }

}

Label_TooManyGraphs:
}
//---------------------------------------------------------------------------

bool TForm5::ParamsFromList(AnsiString Str,int &Node1,int &Node2,int &type, int &Paramtype, AnsiString &AxeTitle)
{
bool retValue;
int a,i,j;
AnsiString tmpStr;


//0-узел
//1-линия
//2-трансформатор
//3-линия либо трансформатор
if((a=Str.AnsiPos("Узел"))>0)
{
 type=0;
 tmpStr="";
 for(i=a+5;i<=Str.Length();i++)
  {
  if(Str[i]!=':')
   tmpStr+=Str[i];
  else
   break;
  }
 Node1=tmpStr.ToIntDef(-1);
 if(Str.AnsiPos("Модуль напряжения")>0)
  {Paramtype=0;AxeTitle="U, кВ";}
 else if(Str.AnsiPos("Фаза напряжения")>0)
  {Paramtype=1;AxeTitle="f, Град";}
 else if(Str.AnsiPos("Активная составляющая напряжения")>0)
  {Paramtype=2;AxeTitle="U', кВ";}
 else if(Str.AnsiPos("Реактивная составляющая напряжения")>0)
  {Paramtype=3;AxeTitle="U'', кВ";}
 else
  Paramtype=-1;
retValue=true;
}
else if((a=Str.AnsiPos("Ветвь"))>0)
{
type=3;
tmpStr="";
j=Str.Length();
for(i=a+6;i<=Str.Length();i++)
  {
  if(Str[i]!='-')
   tmpStr+=Str[i];
  else
   {
    a=i;
    break;
   }
  }
Node1=tmpStr.ToIntDef(-1);
tmpStr="";
for(i=a+1;i<=Str.Length();i++)
  {
  if(Str[i]!=':')
   tmpStr+=Str[i];
  else
   break;
  }
Node2=tmpStr.ToIntDef(-1);
if(Str.AnsiPos("Модуль тока")>0)
 {Paramtype=0;AxeTitle="I, A";}
else if(Str.AnsiPos("Фаза тока")>0)
 {Paramtype=1;AxeTitle="f, Град";}
else if(Str.AnsiPos("Активная составляющая тока")>0)
 {Paramtype=2;AxeTitle="I', A";}
else if(Str.AnsiPos("Реактивная составляющая тока")>0)
 {Paramtype=3;AxeTitle="I'', A";}
//для трансформаторов
else if(Str.AnsiPos("Активные потери в меди")>0)
 {Paramtype=4;type=2;AxeTitle="dP, МВт";}
else if(Str.AnsiPos("Реактивные потери в меди")>0)
 {Paramtype=5;type=2;AxeTitle="dQ, МВт";}
else if(Str.AnsiPos("Активные потери в стали")>0)
 {Paramtype=6;type=2;AxeTitle="dP, МВт";}
else if(Str.AnsiPos("Реактивные потери в стали")>0)
 {Paramtype=7;type=2;AxeTitle="dQ, МВт";}
//для линий
else if(Str.AnsiPos("Модуль потерь напряжения")>0)
 {Paramtype=4;type=1;AxeTitle="|U|, кВ";}
else if(Str.AnsiPos("Активные потери мощности в линии")>0)
 {Paramtype=5;type=1;AxeTitle="dP, МВт";}
else if(Str.AnsiPos("Реактивные потери мощности в линии")>0)
 {Paramtype=6;type=1;AxeTitle="dQ, МВт";}
else
 Paramtype=-1;
retValue=true;
}
return retValue;
}
//---------------------------------------------------------------------------
void TForm5::AddSeriesToGraph(int Node1X,int Node2X,int Node1Y,int Node2Y,int typeX,int ParamtypeX,int typeY,int ParamtypeY,int GraphFormNo,int SerNo,int Count,AnsiString AxeXTitle,AnsiString AxeYTitle)
{
int i,k;
AnsiString SerName;
double *SerValueX=new double [Count];
double *SerValueY=new double [Count];

if(typeX==0)
 {
 for(i=0;i<ValsGraph[0].NodesNumber;i++)
  if(Node1X==ValsGraph[0].RealNode[i])
   {
    k=i;
    break;
   }
  }
else if(typeX!=4)
 {
  for(i=0;i<ValsGraph[0].NodesNumber;i++)
  if(Node1X==ValsGraph[0].RealNode[ValsGraph[0].Node1[i]-1]&&Node2X==ValsGraph[0].RealNode[ValsGraph[0].Node2[i]-1])
   {
    k=i;
    break;
   }
 }
if(typeX!=4)
  FillSeries(typeX,ParamtypeX,k,Count,SerValueX);
else  //итерация цикла утяжеления
 for(i=0;i<Count;i++)
  SerValueX[i]=i;


if(typeY==0)
 {
  for(i=0;i<ValsGraph[0].NodesNumber;i++)
  if(Node1Y==ValsGraph[0].RealNode[i])
   {
    k=i;
    break;
   }
 }
else
 {
  for(i=0;i<ValsGraph[0].mswNumber;i++)
  if(Node1Y==ValsGraph[0].RealNode[ValsGraph[0].Node1[i]-1]&&Node2Y==ValsGraph[0].RealNode[ValsGraph[0].Node2[i]-1])
   {
    k=i;
    break;
   }
 }
FillSeries(typeY,ParamtypeY,k,Count,SerValueY);

if(typeY==0)//для узлов
{
if(ParamtypeY==0)
 {SerName="Модуль напряжения";}
else if(ParamtypeY==1)
 {SerName="Фаза напряжения";}
else if(ParamtypeY==2)
 {SerName="Активная составляющая напряжения";}
else if(ParamtypeY==3)
 {SerName="Реактивная составляющая напряжения";}
}
else if(typeY==1)//линия
{
if(ParamtypeY==4)
 {SerName="Модуль потерь напряжения";}
else if(ParamtypeY==5)
 {SerName="Активные потери мощности в линии";}
else if(ParamtypeY==6)
 {SerName="Реактивные потери мощности в линии";}
}
else if(typeY==2)//трансформатор
{
if(ParamtypeY==4)
 {SerName="Активные потери в меди";}
else if(ParamtypeY==5)
 {SerName="Реактивные потери в меди";}
else if(ParamtypeY==6)
 {SerName="Активные потери в стали";}
else if(ParamtypeY==7)
 {SerName="Реактивные потери в стали";}

}
else if(typeY==3)//трансформатор или линия
{
if(ParamtypeY==0)
 {SerName="Модуль тока";}
else if(ParamtypeY==1)
 {SerName="Фаза тока";}
else if(ParamtypeY==2)
 {SerName="Активная составляющая тока";}
else if(ParamtypeY==3)
 {SerName="Реактивная составляющая тока";}
}


GraphForm[GraphFormNo]->Ser[SerNo]=new TLineSeries(GraphForm[GraphFormNo]->Chart1);
if(typeY==0)
 GraphForm[GraphFormNo]->Ser[SerNo]->Title=SerName+": узел "+AnsiString(Node1Y);
else
 GraphForm[GraphFormNo]->Ser[SerNo]->Title=SerName+": ветвь "+AnsiString(Node1Y)+"-"+AnsiString(Node2Y);
GraphForm[GraphFormNo]->Ser[SerNo]->Name="MySereis"+AnsiString(SerNo);
GraphForm[GraphFormNo]->Ser[SerNo]->SeriesColor=SerColor[SerNo];
GraphForm[GraphFormNo]->Chart1->LeftAxis->Title->Caption=AxeYTitle;
GraphForm[GraphFormNo]->Chart1->BottomAxis->Title->Caption=AxeXTitle;
for(i=0;i<Count;i++)
 {
  GraphForm[GraphFormNo]->Ser[SerNo]->AddXY(SerValueX[i],SerValueY[i],"",SerColor[SerNo]);
 }
GraphForm[GraphFormNo]->Chart1->AddSeries(GraphForm[GraphFormNo]->Ser[SerNo]);

}
//---------------------------------------------------------------------------
void TForm5::FillSeries(int type,int Paramtype,int ValsIndex,int Count,double *SerValue)
{
int i;
if(type==0)//для узлов
{
if(Paramtype==0)
 for(i=0;i<Count;i++)
  SerValue[i]=abs(ValsGraph[i].u[ValsIndex]);
else if(Paramtype==1)
 for(i=0;i<Count;i++)
  SerValue[i]=arg(ValsGraph[i].u[ValsIndex]);
else if(Paramtype==2)
 for(i=0;i<Count;i++)
  SerValue[i]=real(ValsGraph[i].u[ValsIndex]);
else if(Paramtype==3)
 for(i=0;i<Count;i++)
  SerValue[i]=imag(ValsGraph[i].u[ValsIndex]);
}
else if(type==3)//для линии или трансформаторов
{
if(Paramtype==0)
 for(i=0;i<Count;i++)
  SerValue[i]=abs(ValsGraph[i].I[ValsIndex]);
else if(Paramtype==1)
 for(i=0;i<Count;i++)
  SerValue[i]=arg(ValsGraph[i].I[ValsIndex]);
else if(Paramtype==2)
 for(i=0;i<Count;i++)
  SerValue[i]=real(ValsGraph[i].I[ValsIndex]);
else if(Paramtype==3)
 for(i=0;i<Count;i++)
  SerValue[i]=imag(ValsGraph[i].I[ValsIndex]);
}
else if(type==1)//для линии
{
if(Paramtype==4)
 for(i=0;i<Count;i++)
  SerValue[i]=ValsGraph[i].DU[ValsIndex];
else if(Paramtype==5)
 for(i=0;i<Count;i++)
  SerValue[i]=real(ValsGraph[i].DS[ValsIndex]);
else if(Paramtype==6)
 for(i=0;i<Count;i++)
  SerValue[i]=imag(ValsGraph[i].DS[ValsIndex]);
}
else if(type==2)//для трансформаторов
{
if(Paramtype==4)
 for(i=0;i<Count;i++)
  SerValue[i]=real(ValsGraph[i].DSTM[ValsIndex]);
else if(Paramtype==5)
 for(i=0;i<Count;i++)
  SerValue[i]=imag(ValsGraph[i].DSTM[ValsIndex]);
else if(Paramtype==6)
 for(i=0;i<Count;i++)
  SerValue[i]=real(ValsGraph[i].DSTC[ValsIndex]);
else if(Paramtype==7)
 for(i=0;i<Count;i++)
  SerValue[i]=imag(ValsGraph[i].DSTC[ValsIndex]);
}
}


void __fastcall TForm5::Button5Click(TObject *Sender)
{


for(int i=0;i<ListGraph->Items->Count;i++)
 if(ListGraph->Selected[i]==true)
  {
   ListGraph->Items->Delete(i);
   i--;
  }
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button7Click(TObject *Sender)
{
AnsiString GraphName;
if(comboy1->Text.Pos("-")==0)
{
for(int i=0;i<ValsGraph[0].NodesNumber;i++)
 {
  GraphName="Узел "+AnsiString(ValsGraph[0].RealNode[i])+": "+comboy2->Text;
  ListGraph->Items->Add(GraphName);
 }
}
else
 ShowMessage("Необходимо выбрать параметр узла");
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button6Click(TObject *Sender)
{
bool IsLine;
AnsiString GraphName;
if(comboy1->Text.Pos("-")>0)
  {
  if(ValsGraph[0].IsLine[comboy1->ItemIndex-ValsGraph[0].NodesNumber]==true)
   IsLine=true;
  else
   IsLine=false;

  for (int i=0;i<comboy1->Items->Count-ValsGraph[0].NodesNumber;i++)
   if(ValsGraph[0].IsLine[i]==IsLine)
    {
     GraphName="Ветвь "+AnsiString(ValsGraph[0].RealNode[ValsGraph[0].Node1[i]-1])+"-"+AnsiString(ValsGraph[0].RealNode[ValsGraph[0].Node2[i]-1])+": "+comboy2->Text;
     ListGraph->Items->Add(GraphName);
    }
  }
else
 ShowMessage("Необходимо выбрать параметр связи");
}
//---------------------------------------------------------------------------

