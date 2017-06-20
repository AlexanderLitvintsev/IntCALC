//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "graph_cpp.h"
#include "Params.h"
#include "rez_grahp.h"
#include "graph_parm.h"
#include "Tools.h"
#include "sobstv.h"
#include "elembase.h"
//#include "reg.h"
#include "ap.h"
#include "det.h"
#include "sobstvqri.h"
#include "RegimTable.h"
#include "Reg.h"
#include "Values.h"
#include "svd.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

//extern Reg *Rascht;


//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{


AnsiString Dir,tmpStr;
int ListCount,CurrElement,Count;
bool ParamsRead=false;
Rasch=new Reg(2);
Vals=new Values[1];
ValsGraph=new Values[1];
RealNodes=false;

Picture=new TImage*[MaxNumberOfElem];
Element=new MyImage*[MaxNumberOfElem];
StruList=new struct Stru[1]; 


NewFile=true;
HasChanged=false;
SelectedObject=NULL;
NumberOfElem=0;
MouseBtnLeftIsDown=false;
LinesNumber=0;
ObjectIsClicked=false;
MainElementNo=-1;
BalanceElementNo=-1;
MaxNumberOfIt=1;
NumberOfIters=500;
f0=50;
dUmax=0.05;
dfmax=10;
dKnmax=0.1;
dTmax=5;
dVrimax=200;
LimitNebals=true;

precision=0.001;
IsTxtsOut=false;
Mk=1;

UseQminQmax=true;

Dir=GetCurrentDir();




TStringList *ParamsList=new TStringList;
TStringList *OneStr=new TStringList;
try
   {
    OneStr->LoadFromFile(Dir+"\\params.ini");
//    NumberOfElements=StrToInt();
    Count=18;
    tmpStr="";
    for (int i=18;i<=OneStr->Strings[OneStr->IndexOfName("NumberOfElements")].Length();i++)
       {
        tmpStr+=OneStr->Strings[OneStr->IndexOfName("NumberOfElements")][Count];
        Count++;
       }
   }
   __finally
   {
   delete OneStr;
   }
NumberOfElements=StrToInt(tmpStr);
if(NumberOfElements>0)
  TemplObj=new TemplateObject[NumberOfElements];
else
  goto finish;


try
   {
    ParamsList->LoadFromFile(Dir+"\\params.ini");

    ListCount=ParamsList->Count;
    CurrElement=-1;
    for (int i=0;i<ListCount;i++)
      {
       if(ParamsList->Strings[i]=="[object]"&&ParamsRead==false)
         {
          ParamsRead=true;
          CurrElement++;
         }
	   else if(ParamsList->Strings[i].Pos("Name=")>0)
		   TemplObj[CurrElement].ElementName=ParamsList->Strings[i].SubString(6,ParamsList->Strings[i].Length());
	   else if(ParamsList->Strings[i].Pos("type=")>0)
		   TemplObj[CurrElement].type=StrToInt(ParamsList->Strings[i].SubString(6,ParamsList->Strings[i].Length()));
	   else if(ParamsList->Strings[i].Pos("Picture=")>0)
		  {
		   TemplObj[CurrElement].PictureFileName=Dir+"\\pics\\"+ParamsList->Strings[i].SubString(9,ParamsList->Strings[i].Length());
		   TemplObj[CurrElement].IcoPictureFileName=Dir+"\\pics\\icos\\"+ParamsList->Strings[i].SubString(9,ParamsList->Strings[i].Length());
		  }
	   else if(ParamsList->Strings[i].Pos("NumberOfPoints=")>0)
          TemplObj[CurrElement].NumberOfPoints=StrToInt(ParamsList->Strings[i].SubString(16,ParamsList->Strings[i].Length()));

       for(int k=0;k<MaxNodesInElem;k++)
          {
		   if(ParamsList->Strings[i].Pos("Point"+AnsiString(k+1)+'=')>0)
           {
           tmpStr="";
           Count=0;
           for (int j=8;j<=ParamsList->Strings[i].Length();j++)
             {
              if(ParamsList->Strings[i][j]!=';')
                {
                 tmpStr+=ParamsList->Strings[i][j];

                }

              else if(ParamsList->Strings[i][j]==';')
                {
                 Count++;
                 switch(Count)
                   {
                    case 1:
                        TemplObj[CurrElement].RectCoord[k].Left=StrToInt(tmpStr);
                        tmpStr="";
                        break;
                    case 2:
                        TemplObj[CurrElement].RectCoord[k].Top=StrToInt(tmpStr);
                        tmpStr="";
                        break;
                    case 3:
                        TemplObj[CurrElement].RectCoord[k].Right=StrToInt(tmpStr);
                        tmpStr="";
                        break;
                    case 4:
                        TemplObj[CurrElement].RectCoord[k].Bottom=StrToInt(tmpStr);
                        tmpStr="";
                        break;

                   }
                }

             }

           }
          }

       for(int k=0;k<MaxNodesInElem;k++)
         {
		 if(ParamsList->Strings[i].Pos("Dir"+AnsiString(k+1)+'=')>0)
            TemplObj[CurrElement].Direction[k]=StrToInt(ParamsList->Strings[i].SubString(6,ParamsList->Strings[i].Length()));
         }
	   if(ParamsList->Strings[i].Pos("Parameters=")>0)
          {
           Count=0;
           tmpStr="";
           for(int j=12;j<=ParamsList->Strings[i].Length();j++)
             {
              if(ParamsList->Strings[i][j]!=';')
                {
                 tmpStr+=ParamsList->Strings[i][j];
                }
              else if(ParamsList->Strings[i][j]==';')
                {
                 Count++;
                 switch(Count)
                   {
                    case 1:
                          TemplObj[CurrElement].NumberOfParams=StrToInt(tmpStr);
                          tmpStr="";
                          break;
                    default:
                          TemplObj[CurrElement].ParamsNames[Count-2]=tmpStr;
                          tmpStr="";
                          break;
                   }
                }
             }
          }
	   else if(ParamsList->Strings[i].Pos("[end]")>0)
          {
           ParamsRead=false;
          }
      }
   }
   __finally
   {
    delete ParamsList;
   }
Image=new TImage*[NumberOfElements];
for (int i=0;i<NumberOfElements;i++)
  {
   Image[i]=new TImage(Form1);
   Image[i]->Parent=Form1;
   Image[i]->Visible=false;
   Image[i]->Enabled=false;
   Image[i]->AutoSize=true;
   Image[i]->Picture->LoadFromFile(TemplObj[i].PictureFileName);
  }
  goto well;
finish:
   ShowMessage("������ ������ ����� ���������!");
   Application->Terminate();

well:

}
//---------------------------------------------------------------------------

MyImage::MyImage(TForm *Form1):TImage(Form1)
{
NumberOfPoints=0;
for (int i=0;i<MaxNodesInElem;i++)
   {
    RectCoord[i]=TRect(0,0,0,0);

    //    RectCoord[3].y=-1;
    PointIsLinked[i]=false;
    IsClicked=false;
    LinkObjectNumber[i]=-1;
    PointIsForwarded[i]=false;
    PointNumber[i]=-1;
    RealPointNumber[i]=-1;

   }
InPointNo[0]=-1;
InPointNo[1]=-1;
InIndex=-1;
BalanceRectNo=-1;

}
//---------------------------------------------------------------------------
void TForm1::AddElement(int ElemNo,AnsiString ElementName,int &NumberOfElem,TPicture *Pic,POINT point,int NumberOfPoints,TRect *Rect,int NumberOfParams,AnsiString *ParamNames,double *Params,int *Dirs)
{
int i;
HasChanged=true;
NumberOfElem++;

Picture[NumberOfElem-1]=new TImage(Form1);
Picture[NumberOfElem-1]->Parent=Form1;
Picture[NumberOfElem-1]->Visible=false;
Element[NumberOfElem-1]=new MyImage(Form1);
Element[NumberOfElem-1]->Parent=Form1;

Element[NumberOfElem-1]->OnMouseMove=Image1MouseMove;
Element[NumberOfElem-1]->OnMouseDown=Image1MouseDown;
Element[NumberOfElem-1]->OnMouseUp=Image1MouseUp;
Element[NumberOfElem-1]->OnClick=Image1Click;

Element[NumberOfElem-1]->Picture=Pic;
Picture[NumberOfElem-1]->Picture=Pic;

Element[NumberOfElem-1]->Left=point.x;
Element[NumberOfElem-1]->Top=point.y;

Element[NumberOfElem-1]->Height=Pic->Height;

Element[NumberOfElem-1]->Width=Pic->Width;
Element[NumberOfElem-1]->ObjectNumber=NumberOfElem-1;
Element[NumberOfElem-1]->NumberOfPoints=NumberOfPoints;
Element[NumberOfElem-1]->Transparent=false;
Element[NumberOfElem-1]->PopupMenu=PopMenu;

for(i=0;i<NumberOfPoints;i++)
  {
   Element[NumberOfElem-1]->PointIsLinked[i]=false;
   Element[NumberOfElem-1]->RectCoord[i]=Rect[i];
   Element[NumberOfElem-1]->Direction[i]=Dirs[i];
  }
Element[NumberOfElem-1]->NumberOfParams=NumberOfParams;
Element[NumberOfElem-1]->type=type;
for (i=0;i<NumberOfParams;i++)
   {
    Element[NumberOfElem-1]->ParamNames[i]=ParamNames[i];
    Element[NumberOfElem-1]->Params[i]=Params[i];
   }
Element[NumberOfElem-1]->ElementNumber=ElemNo;
Element[NumberOfElem-1]->ElementName=ElementName;
Element[NumberOfElem-1]->MainRectNo=-1;
Element[NumberOfElem-1]->ParFromFile=false;
Element[NumberOfElem-1]->ParFileName="";
Element[NumberOfElem-1]->Stretch=true;
Element[NumberOfElem-1]->Height=MyRound(double(Element[NumberOfElem-1]->Height)*Mk);
Element[NumberOfElem-1]->Width=MyRound(double(Element[NumberOfElem-1]->Width)*Mk);
}
//---------------------------------------------------------------------------
void TForm1::DeleteElement(int ElemNo,int &NumberOfElem)
{
//MyImage *tmpElement;
int i,k=-1;

for(i=0;i<CyclePar[0].Quantity;i++)
  if(CyclePar[i].ElementNo==ElemNo)
   k=i;
if(k>=0)
  for(i=k;i<CyclePar[0].Quantity-1;i++)
  {
   CyclePar[i]=CyclePar[i+1];
  }
if(k>=0)
  {
   CyclePar[k].ElementNo=-1;
   CyclePar[k].Ok=false;
   CyclePar[0].Quantity-=1;
   MaxNumberOfIt=1;
  }
for(i=ElemNo;i<NumberOfElem-1;i++)
  {
   delete Element[i];
   Element[i]=new MyImage(Form1);
   Element[i]->Parent=Form1;
   Element[i]->Stretch=true;
   Element[i]->Visible=false;
   Element[i]->OnMouseMove=Image1MouseMove;
   Element[i]->OnMouseDown=Image1MouseDown;
   Element[i]->OnMouseUp=Image1MouseUp;
   Element[i]->OnClick=Image1Click;
   Element[i]->NumberOfPoints=Element[i+1]->NumberOfPoints;
   Element[i]->NumberOfParams=Element[i+1]->NumberOfParams;
   Element[i]->type=Element[i+1]->type;
   Element[i]->ParFromFile=Element[i+1]->ParFromFile;
   Element[i]->ParFileName=Element[i+1]->ParFileName;
   if(MainElementNo==i+1)
      MainElementNo=i;
   if(BalanceElementNo==i+1)
      BalanceElementNo=i;

   for (int j=0;j<MaxNodesInElem;j++)
    {
     Element[i]->RectCoord[j]=Element[i+1]->RectCoord[j];
     Element[i]->PointIsLinked[j]=Element[i+1]->PointIsLinked[j];
     Element[i]->LinkObjectNumber[j]=Element[i+1]->LinkObjectNumber[j];
     Element[i]->Direction[j]=Element[i+1]->Direction[j];
     Element[i]->RealPointNumber[j]=Element[i+1]->RealPointNumber[j];
    }
   for (int j=0;j<Element[i+1]->NumberOfParams;j++)
    {
     Element[i]->Params[j]=Element[i+1]->Params[j];
     Element[i]->ParamNames[j]=Element[i+1]->ParamNames[j];
    }
   Element[i]->IsClicked=Element[i+1]->IsClicked;
   Element[i]->Left=Element[i+1]->Left;
   Element[i]->Top=Element[i+1]->Top;
   Element[i]->Height=Element[i+1]->Height;
   Element[i]->Width=Element[i+1]->Width;
   Element[i]->Picture=Element[i+1]->Picture;
   Picture[i]->Picture=Picture[i+1]->Picture;
   Element[i]->ObjectNumber=i;
   Element[i]->Transparent=false;
   Element[i]->PopupMenu=PopMenu;
   Element[i]->ElementNumber=Element[i+1]->ElementNumber;
   Element[i]->MainRectNo=Element[i+1]->MainRectNo;
   Element[i]->ElementName=Element[i+1]->ElementName;
   Element[i]->Visible=true;
   for (int j=0;j<LineCoords::LineNumber;j++)
    {
     if(Element[i+1]->ObjectNumber==LineArray[j].Obj1Number)
       LineArray[j].Obj1Number=i;
     else if(Element[i+1]->ObjectNumber==LineArray[j].Obj2Number)
       LineArray[j].Obj2Number=i;
    }
  }
delete Element[NumberOfElem-1];
//delete tmpElement;
if(MainElementNo==ElemNo)
  MainElementNo=-1;
if(BalanceElementNo==ElemNo)
  BalanceElementNo=-1;

NumberOfElem--;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Image1Click(TObject *Sender)
{

static int Node1Number;
POINT CurrPos;
CurrPos=ScreenToClient(Mouse->CursorPos);
MyImage * ptr=dynamic_cast<MyImage*>(Sender);
bool flag=false;
if(SelectedObject!=NULL)
  SelectedObject->Picture=Picture[SelectedObject->ObjectNumber]->Picture;

//if(SelectedObject!=ptr)
{
ptr->Picture=Picture[ptr->ObjectNumber]->Picture;
ptr->Canvas->Brush->Style=bsSolid;
ptr->Canvas->Brush->Color=clBlack;
ptr->Canvas->Pen->Color=clBlack;

for (int i=0;i<ptr->NumberOfPoints;i++)
   if(ptr->MainRectNo!=i)
   ptr->Canvas->Rectangle(ptr->RectCoord[i]);


}
flag=false;
if(ClickedObject!=ptr&&ObjectIsClicked==false)
  {
  for (int j=0;j<ptr->NumberOfPoints;j++)
     if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true&&ptr->PointIsLinked[j]==false)
        flag=true;

  if(flag==true)
   {
    ObjectIsClicked=true;
    MouseX=CurrPos.x;
    MouseY=CurrPos.y;
    Form1->Canvas->MoveTo(MouseX,MouseY);
    ClickedObject=ptr;
    flag=false;
    for (int j=0;j<ptr->NumberOfPoints;j++)
       if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true&&ptr->PointIsLinked[j]==false)
         Node1Number=j;
    flag=false;

   }
  }
else if(ObjectIsClicked==true&&ClickedObject!=ptr)
  {
  for (int j=0;j<ptr->NumberOfPoints;j++)
     if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true&&ptr->PointIsLinked[j]==false)
        flag=true;
    if(flag==true)
      {
       for(int i=0;i<ClickedObject->NumberOfPoints;i++)
           if(PointInRect(MouseX-ClickedObject->Left,MouseY-ClickedObject->Top,ClickedObject->RectCoord[i])==true)
             {
              if((ClickedObject->Direction[i]==UP||ClickedObject->Direction[i]==DOWN)&&ClickedObject->ElementName=="����")
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2;
                 MouseX=ClickedObject->Left+MyRound(double(ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2)*Mk);
                 MouseY=ClickedObject->Top+MyRound(double(ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2)*Mk);

                }
              else if(ClickedObject->Direction[i]==UP)
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Top;
                 MouseX=ClickedObject->Left+MyRound(double(ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2)*Mk);
                 MouseY=ClickedObject->Top+MyRound(double(ClickedObject->RectCoord[i].Top)*Mk);

                }
              else if(ClickedObject->Direction[i]==DOWN)
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Bottom;
                 MouseX=ClickedObject->Left+MyRound(double(ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2)*Mk);
                 MouseY=ClickedObject->Top+int(double(ClickedObject->RectCoord[i].Bottom)*Mk);

                }
              else if((ClickedObject->Direction[i]==LEFT||ClickedObject->Direction[i]==RIGHT)&&ClickedObject->ElementName=="����")
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2;
                 MouseX=ClickedObject->Left+MyRound(double(ClickedObject->RectCoord[i].Left+(ClickedObject->RectCoord[i].Right-ClickedObject->RectCoord[i].Left)/2)*Mk);
                 MouseY=ClickedObject->Top+MyRound(double(ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2)*Mk);

                }

              else if(ClickedObject->Direction[i]==LEFT)
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Left;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2;
                 MouseX=ClickedObject->Left+MyRound(double(ClickedObject->RectCoord[i].Left)*Mk);
                 MouseY=ClickedObject->Top+MyRound(double(ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2)*Mk);

                }
              else if(ClickedObject->Direction[i]==RIGHT)
                {
//                 MouseX=ClickedObject->Left+ClickedObject->RectCoord[i].Right;
//                 MouseY=ClickedObject->Top+ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2;
                 MouseX=ClickedObject->Left+int(double(ClickedObject->RectCoord[i].Right)*Mk);
                 MouseY=ClickedObject->Top+MyRound(double(ClickedObject->RectCoord[i].Top+(ClickedObject->RectCoord[i].Bottom-ClickedObject->RectCoord[i].Top)/2)*Mk);

                }
             }
       for(int i=0;i<ptr->NumberOfPoints;i++)
           if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[i])==true)
             {
              if((ptr->Direction[i]==UP||ptr->Direction[i]==DOWN)&&ptr->ElementName=="����")
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2;
                 CurrPos.x=ptr->Left+MyRound(double(ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2)*Mk);
                 CurrPos.y=ptr->Top+MyRound(double(ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2)*Mk);
                }

              else if(ptr->Direction[i]==UP)
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Top;
                 CurrPos.x=ptr->Left+MyRound(double(ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2)*Mk);
                 CurrPos.y=ptr->Top+MyRound(double(ptr->RectCoord[i].Top)*Mk);

                }
              else if(ptr->Direction[i]==DOWN)
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Bottom;
                 CurrPos.x=ptr->Left+MyRound(double(ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2)*Mk);
                 CurrPos.y=ptr->Top+int(double(ptr->RectCoord[i].Bottom)*Mk);

                }
              else if((ptr->Direction[i]==LEFT||ptr->Direction[i]==RIGHT)&&ptr->ElementName=="����")
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2;
                 CurrPos.x=ptr->Left+MyRound(double(ptr->RectCoord[i].Left+(ptr->RectCoord[i].Right-ptr->RectCoord[i].Left)/2)*Mk);
                 CurrPos.y=ptr->Top+MyRound(double(ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2)*Mk);

                }

              else if(ptr->Direction[i]==LEFT)
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Left;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2;
                 CurrPos.x=ptr->Left+MyRound(double(ptr->RectCoord[i].Left)*Mk);
                 CurrPos.y=ptr->Top+MyRound(double(ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2)*Mk);

                }
              else if(ptr->Direction[i]==RIGHT)
                {
//                 CurrPos.x=ptr->Left+ptr->RectCoord[i].Right;
//                 CurrPos.y=ptr->Top+ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2;
                 CurrPos.x=ptr->Left+int(double(ptr->RectCoord[i].Right)*Mk);
                 CurrPos.y=ptr->Top+MyRound(double(ptr->RectCoord[i].Top+(ptr->RectCoord[i].Bottom-ptr->RectCoord[i].Top)/2)*Mk);

                }
             }

       Form1->Canvas->Pen->Color=clWhite;
       Form1->Canvas->MoveTo(MouseX,MouseY);
       Form1->Canvas->LineTo(OldMousePos.x,OldMousePos.y);
       Form1->Canvas->Pen->Color=clWhite;
       Form1->Canvas->MoveTo(MouseX,MouseY);
       Form1->Canvas->LineTo(CurrPos.x,CurrPos.y);
       Form1->Canvas->Pen->Color=clBlack;
       ObjectIsClicked=false;
       for(int i=0;i<ptr->NumberOfPoints;i++)
          {
           if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[i])==true)
               ptr->PointIsLinked[i]=true;
          }
       for(int i=0;i<ClickedObject->NumberOfPoints;i++)
          {
           if(PointInRect(MouseX-ClickedObject->Left,MouseY-ClickedObject->Top,ClickedObject->RectCoord[i])==true)
               ClickedObject->PointIsLinked[i]=true;

          }

       tmpLineArray=new tmpLineCoords[LineCoords::LineNumber];
       for (int i=0;i<LineCoords::LineNumber;i++)
          {
           tmpLineArray[i].Start=LineArray[i].Start;
           tmpLineArray[i].End=LineArray[i].End;
           tmpLineArray[i].Obj1Number=LineArray[i].Obj1Number;
           tmpLineArray[i].Obj2Number=LineArray[i].Obj2Number;
           tmpLineArray[i].Node1Number=LineArray[i].Node1Number;
           tmpLineArray[i].Node2Number=LineArray[i].Node2Number;

          }
       delete []LineArray;
       LineArray=new LineCoords[LineCoords::LineNumber+1];
       for (int i=0;i<LineCoords::LineNumber;i++)
          {
           LineArray[i].Start=tmpLineArray[i].Start;
           LineArray[i].End=tmpLineArray[i].End;
           LineArray[i].Obj1Number=tmpLineArray[i].Obj1Number;
           LineArray[i].Obj2Number=tmpLineArray[i].Obj2Number;
           LineArray[i].Node1Number=tmpLineArray[i].Node1Number;
           LineArray[i].Node2Number=tmpLineArray[i].Node2Number;

          }
       delete []tmpLineArray;
       LineCoords::LineNumber++;//�.�. ���������� ���� �������
       LineArray[LineCoords::LineNumber-1].Start.x=MouseX;
       LineArray[LineCoords::LineNumber-1].Start.y=MouseY;
       LineArray[LineCoords::LineNumber-1].End=CurrPos;
       for(int i=0;i<ptr->NumberOfPoints;i++)
         if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[i])==true&&ptr->PointIsLinked[i]==true)
            LineArray[LineCoords::LineNumber-1].Node2Number=i;
       LineArray[LineCoords::LineNumber-1].Node1Number=Node1Number;
       LineArray[LineCoords::LineNumber-1].Obj2Number=ptr->ObjectNumber;
       LineArray[LineCoords::LineNumber-1].Obj1Number=ClickedObject->ObjectNumber;

       ClickedObject=NULL;

       Redraw();
      }
  }

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Image1MouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
POINT CurrPos;
CurrPos=ScreenToClient(Mouse->CursorPos);
bool flag=false;

if(Button==mbLeft)
{

MouseBtnLeftIsDown=true;


    CurrMouseX=CurrPos.x;
    CurrMouseY=CurrPos.y;
}

else if(Button==mbRight)
{

MyImage *ptr=dynamic_cast<MyImage*>(Sender);
ImageForRightBtn=ptr;
  for(int j=0;j<ptr->NumberOfPoints;j++)
     if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true)
       flag=true;

     if(flag==true)
     {
      flag=false;
      N9->Enabled=true;
      N9->Visible=true;
      N16->Enabled=true;
      N16->Visible=true;
      N32->Enabled=true;
      N32->Visible=true;
      for(int j=0;j<ptr->NumberOfPoints;j++)
         if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true)
            {
            ptr->MainRectNo=j;
            ptr->BalanceRectNo=j;

            }
     }
  else
     {
      N9->Enabled=false;
      N9->Visible=false;
      N16->Enabled=false;
      N16->Visible=false;
      N32->Enabled=false;
      N32->Visible=false;
      ptr->MainRectNo=-1;
     }
PopMenu->Popup(Mouse->CursorPos.x,Mouse->CursorPos.y);

}

}
//---------------------------------------------------------------------------


void __fastcall TForm1::Image1MouseMove(TObject *Sender, TShiftState Shift,
      int X, int Y)
{
static bool NeedToChange=false;
MyImage * ptr=dynamic_cast<MyImage*>(Sender);

OnObject=true;
POINT CurrPos=ScreenToClient(Mouse->CursorPos);;
tmpImage=ptr;

ptr->Hint=ptr->ElementName;
if(ptr!=SelectedObject)

   {

    ptr->Canvas->Brush->Style=bsSolid;
    ptr->Canvas->Brush->Color=clRed;
    ptr->Canvas->Pen->Color=clRed;

    for (int i=0;i<ptr->NumberOfPoints;i++)
       if(ptr->MainRectNo!=i)
         ptr->Canvas->FillRect(ptr->RectCoord[i]);
    ptr->IsClicked=true;

   }
    if(NeedToChange==true)
    for(int i=0;i<ptr->NumberOfPoints;i++)
        if(ptr->MainRectNo!=i)
         {
          ptr->Canvas->Pen->Color=clBlack;
          ptr->Canvas->FillRect(ptr->RectCoord[i]);
          NeedToChange=false;
         }

    for(int j=0;j<ptr->NumberOfPoints;j++)
       if(PointInRect(CurrPos.x-ptr->Left,CurrPos.y-ptr->Top,ptr->RectCoord[j])==true)
//       if(PointInRect(X,Y,ptr->RectCoord[j])==true)
         {
          ptr->Canvas->Pen->Color=clGreen;
          ptr->Canvas->Brush->Color=clGreen;
          ptr->Canvas->FillRect(ptr->RectCoord[j]);
          ptr->Canvas->Pen->Color=clBlack;
          ptr->Canvas->Brush->Color=clBlack;
          NeedToChange=true;
         }

if(MouseBtnLeftIsDown==true&&ClickedObject==NULL)
   {
    HasChanged=true;
    if(ClickedObject==NULL)
       {
        ObjectIsClicked=false;
        ClickedObject=NULL;
       }
    CurrPos=ScreenToClient(Mouse->CursorPos);
    ptr->Left+=CurrPos.x-CurrMouseX;
    ptr->Top+=CurrPos.y-CurrMouseY;
    CurrPos=ScreenToClient(Mouse->CursorPos);
    CurrMouseX=CurrPos.x;
    CurrMouseY=CurrPos.y;
   }
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Image1MouseUp(TObject *Sender, TMouseButton Button,
      TShiftState Shift, int X, int Y)
{
MyImage * ptr=dynamic_cast<MyImage*>(Sender);
MouseBtnLeftIsDown=false;

for(int i=0;i<LineCoords::LineNumber;i++)
  {
   if(LineArray[i].Obj1Number==ptr->ObjectNumber)
     {
      switch (ptr->Direction[LineArray[i].Node1Number])
        {
         case LEFT:
             {
              LineArray[i].Start.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Left)*Mk);
              LineArray[i].Start.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Top+abs(ptr->RectCoord[LineArray[i].Node1Number].Top-ptr->RectCoord[LineArray[i].Node1Number].Bottom)/2)*Mk);
              break;
             }
         case RIGHT:
             {
//              LineArray[i].Start.x=ptr->Left+ptr->RectCoord[LineArray[i].Node1Number].Right;
//              LineArray[i].Start.y=ptr->Top+ptr->RectCoord[LineArray[i].Node1Number].Top+abs(ptr->RectCoord[LineArray[i].Node1Number].Top-ptr->RectCoord[LineArray[i].Node1Number].Bottom)/2;
              LineArray[i].Start.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Right)*Mk);
              LineArray[i].Start.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Top+abs(ptr->RectCoord[LineArray[i].Node1Number].Top-ptr->RectCoord[LineArray[i].Node1Number].Bottom)/2)*Mk);

              break;
             }
         case UP:
             {
              LineArray[i].Start.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Left+(ptr->RectCoord[LineArray[i].Node1Number].Right-ptr->RectCoord[LineArray[i].Node1Number].Left)/2)*Mk);
              LineArray[i].Start.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Top)*Mk);
              break;
             }
         case DOWN:
             {
              LineArray[i].Start.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Left+(ptr->RectCoord[LineArray[i].Node1Number].Right-ptr->RectCoord[LineArray[i].Node1Number].Left)/2)*Mk);
              LineArray[i].Start.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Bottom)*Mk);
              break;
             }

        }
      if(ptr->ElementName=="����")
        {
         LineArray[i].Start.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Left+(ptr->RectCoord[LineArray[i].Node1Number].Right-ptr->RectCoord[LineArray[i].Node1Number].Left)/2)*Mk);
         LineArray[i].Start.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node1Number].Top+abs(ptr->RectCoord[LineArray[i].Node1Number].Top-ptr->RectCoord[LineArray[i].Node1Number].Bottom)/2)*Mk);
        }

     }
   else if(LineArray[i].Obj2Number==ptr->ObjectNumber)
     {
      switch (ptr->Direction[LineArray[i].Node2Number])
        {
         case LEFT:
             {
              LineArray[i].End.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Left)*Mk);
              LineArray[i].End.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Top+abs(ptr->RectCoord[LineArray[i].Node2Number].Top-ptr->RectCoord[LineArray[i].Node2Number].Bottom)/2)*Mk);

              break;
             }
         case RIGHT:
             {
              LineArray[i].End.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Right)*Mk);
              LineArray[i].End.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Top+abs(ptr->RectCoord[LineArray[i].Node2Number].Top-ptr->RectCoord[LineArray[i].Node2Number].Bottom)/2)*Mk);
              break;
             }
         case UP:
             {
              LineArray[i].End.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Left+(ptr->RectCoord[LineArray[i].Node2Number].Right-ptr->RectCoord[LineArray[i].Node2Number].Left)/2)*Mk);
              LineArray[i].End.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Top)*Mk);
              break;
             }
         case DOWN:
             {
              LineArray[i].End.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Left+(ptr->RectCoord[LineArray[i].Node2Number].Right-ptr->RectCoord[LineArray[i].Node2Number].Left)/2)*Mk);
              LineArray[i].End.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Bottom)*Mk);
              break;
             }

        }
      if(ptr->ElementName=="����")
        {
         LineArray[i].End.x=ptr->Left+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Left+(ptr->RectCoord[LineArray[i].Node2Number].Right-ptr->RectCoord[LineArray[i].Node2Number].Left)/2)*Mk);
         LineArray[i].End.y=ptr->Top+MyRound(double(ptr->RectCoord[LineArray[i].Node2Number].Top+abs(ptr->RectCoord[LineArray[i].Node2Number].Top-ptr->RectCoord[LineArray[i].Node2Number].Bottom)/2)*Mk);
        }

     }
  }

Redraw();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Obj1Execute(TObject *Sender)
{

TToolButton * ptr=dynamic_cast<TToolButton*>(Sender);
int ind=ptr->Index;
PicToCreate=Image[ind]->Picture;
NumberOfNodes=TemplObj[ind].NumberOfPoints;
for (int i=0;i<NumberOfNodes;i++)
{
NodesCoord[i]=TemplObj[ind].RectCoord[i];
Dirs[i]=TemplObj[ind].Direction[i];
}
NumberOfElemParams=TemplObj[ind].NumberOfParams;
for (int i=0;i<NumberOfElemParams;i++)
  {
   ParamNames[i]=TemplObj[ind].ParamsNames[i];
   Params[i]=0;
  }
ElementName=TemplObj[ind].ElementName;
type=TemplObj[ind].type;
ElementNo=ind;
CreatePicFlag=true;

}
//---------------------------------------------------------------------------




void __fastcall TForm1::FormClick(TObject *Sender)
{
POINT ImagePos;
if (CreatePicFlag==true)
   {
    if(SelectedObject!=NULL)
       {
        SelectedObject->Picture=Image[SelectedObject->ElementNumber]->Picture;
        SelectedObject=NULL;
       }
    ImagePos=ScreenToClient(Mouse->CursorPos);
    AddElement(ElementNo,ElementName,NumberOfElem,PicToCreate,ImagePos,NumberOfNodes,NodesCoord,NumberOfElemParams,ParamNames,Params,Dirs);
    CreatePicFlag=false;
   }
else if(SelectedObject!=NULL&&ObjectIsClicked==false)
   {
    SelectedObject->Picture=Image[SelectedObject->ElementNumber]->Picture;
    SelectedObject==NULL;
   }
if(ObjectIsClicked==true)
   {
    ObjectIsClicked=false;
    Form1->Canvas->Pen->Color=clWhite;
    ImagePos=ScreenToClient(Mouse->CursorPos);
    Form1->Canvas->MoveTo(MouseX,MouseY);
    Form1->Canvas->LineTo(ImagePos.x,ImagePos.y);
    Form1->Canvas->Pen->Color=clBlack;
    ClickedObject=NULL;
   }

}
//---------------------------------------------------------------------------



bool TForm1::PointInRect(int x,int y,TRect Rect)
{
int x1,y1;
x1=MyRound(double(x)/Mk);
y1=MyRound(double(y)/Mk);
if(x1>=Rect.Left&&x1<=Rect.Right&&y1>=Rect.Top&&y1<=Rect.Bottom)
   return true;
else
   return false;
}

//---------------------------------------------------------------------------
void __fastcall TForm1::FormMouseMove(TObject *Sender, TShiftState Shift,
      int X, int Y)
{

POINT CurrPos;
if(OnObject==true&&tmpImage!=SelectedObject)
  {
   tmpImage->Transparent=false;
   tmpImage->Picture=Picture[tmpImage->ObjectNumber]->Picture;
   OnObject=false;

  }

if(ObjectIsClicked==true)
   {
    Form1->Canvas->Pen->Color=clWhite;
    Form1->Canvas->MoveTo(MouseX,MouseY);
    Form1->Canvas->LineTo(OldMousePos.x,OldMousePos.y);
    Form1->Canvas->MoveTo(MouseX,MouseY);
    CurrPos=ScreenToClient(Mouse->CursorPos);
    Form1->Canvas->Pen->Color=clBlack;
    Form1->Canvas->LineTo(CurrPos.x,CurrPos.y);
    OldMousePos=CurrPos;
    //ClickedObject=NULL;
   }
}
//---------------------------------------------------------------------------
LineCoords::LineCoords()
{
//LineNumber=0;
//LineNumber++;
}
//---------------------------------------------------------------------------
LineCoords::~LineCoords()
{
//LineNumber--;
}
//---------------------------------------------------------------------------
void TForm1::RedrawWithoutRefresh()
{
int Obj1,Obj2,Node1,Node2;
POINT P1,P2;
bool is1Sheena,is2Sheena;


for (int i=0;i<LineCoords::LineNumber;i++)
   {
    is1Sheena=false;
    is2Sheena=false;
    Obj1=LineArray[i].Obj1Number;
    Obj2=LineArray[i].Obj2Number;
    Node1=LineArray[i].Node1Number;
    Node2=LineArray[i].Node2Number;

    if(Element[Obj1]->ElementName=="����")
       is1Sheena=true;
    if(Element[Obj2]->ElementName=="����")
       is2Sheena=true;
    P1.x=Element[Obj1]->Left+MyRound(double(Element[Obj1]->RectCoord[Node1].Left)*Mk);
    P1.y=Element[Obj1]->Top+MyRound(double(Element[Obj1]->RectCoord[Node1].Top)*Mk);
    P2.x=Element[Obj2]->Left+MyRound(double(Element[Obj2]->RectCoord[Node2].Left)*Mk);
    P2.y=Element[Obj2]->Top+MyRound(double(Element[Obj2]->RectCoord[Node2].Top)*Mk);

    if(is1Sheena==true)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;

      }
    else if(Element[Obj1]->Direction[Node1]==UP)
      {

       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x,LineArray[i].Start.y-1);
       P1.y-=1;
      }
    else if(Element[Obj1]->Direction[Node1]==DOWN)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top))*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x,LineArray[i].Start.y+1);
       P1.y+=1;
      }
    else if(Element[Obj1]->Direction[Node1]==LEFT)
      {
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x-1,LineArray[i].Start.y);
       P1.x-=1;
      }
    else if(Element[Obj1]->Direction[Node1]==RIGHT)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left))*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x+1,LineArray[i].Start.y);
       P1.x+=1;
      }

    if(is2Sheena==true)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;

      }
    else if(Element[Obj2]->Direction[Node2]==UP)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x,LineArray[i].End.y-1);
       P2.y-=1;
      }
    else if(Element[Obj2]->Direction[Node2]==DOWN)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top))*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x,LineArray[i].End.y+1);
       P2.y+=1;
      }
    else if(Element[Obj2]->Direction[Node2]==LEFT)
      {
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x-1,LineArray[i].End.y);
       P2.x-=1;
      }
    else if(Element[Obj2]->Direction[Node2]==RIGHT)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left))*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x+1,LineArray[i].End.y);
       P2.x+=1;

      }
    Link2Points(P1,P2,Element[Obj1]->Direction[Node1],Element[Obj2]->Direction[Node2],max(Element[Obj1]->Width,Element[Obj2]->Width),max(Element[Obj1]->Height,Element[Obj2]->Height),is1Sheena,is2Sheena);


   }
}
//---------------------------------------------------------------------------
void TForm1::Redraw()
{
int Obj1,Obj2,Node1,Node2;
POINT P1,P2;
Form1->Refresh();
bool is1Sheena,is2Sheena;


for (int i=0;i<LineCoords::LineNumber;i++)
   {
    is1Sheena=false;
    is2Sheena=false;
    Obj1=LineArray[i].Obj1Number;
    Obj2=LineArray[i].Obj2Number;
    Node1=LineArray[i].Node1Number;
    Node2=LineArray[i].Node2Number;

    if(Element[Obj1]->ElementName=="����")
       is1Sheena=true;
    if(Element[Obj2]->ElementName=="����")
       is2Sheena=true;
    P1.x=Element[Obj1]->Left+MyRound(double(Element[Obj1]->RectCoord[Node1].Left)*Mk);
    P1.y=Element[Obj1]->Top+MyRound(double(Element[Obj1]->RectCoord[Node1].Top)*Mk);
    P2.x=Element[Obj2]->Left+MyRound(double(Element[Obj2]->RectCoord[Node2].Left)*Mk);
    P2.y=Element[Obj2]->Top+MyRound(double(Element[Obj2]->RectCoord[Node2].Top)*Mk);

    if(is1Sheena==true)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;

      }
    else if(Element[Obj1]->Direction[Node1]==UP)
      {

       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x,LineArray[i].Start.y-1);
       P1.y-=1;
      }
    else if(Element[Obj1]->Direction[Node1]==DOWN)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left)/2)*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top))*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x,LineArray[i].Start.y+1);
       P1.y+=1;
      }
    else if(Element[Obj1]->Direction[Node1]==LEFT)
      {
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x-1,LineArray[i].Start.y);
       P1.x-=1;
      }
    else if(Element[Obj1]->Direction[Node1]==RIGHT)
      {
       P1.x+=MyRound(double((Element[Obj1]->RectCoord[Node1].Right-Element[Obj1]->RectCoord[Node1].Left))*Mk);
       P1.y+=MyRound(double((Element[Obj1]->RectCoord[Node1].Bottom-Element[Obj1]->RectCoord[Node1].Top)/2)*Mk);
       LineArray[i].Start.x=P1.x;
       LineArray[i].Start.y=P1.y;
       Form1->Canvas->MoveTo(LineArray[i].Start.x,LineArray[i].Start.y);
       Form1->Canvas->LineTo(LineArray[i].Start.x+1,LineArray[i].Start.y);
       P1.x+=1;
      }

    if(is2Sheena==true)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;

      }
    else if(Element[Obj2]->Direction[Node2]==UP)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x,LineArray[i].End.y-1);
       P2.y-=1;
      }
    else if(Element[Obj2]->Direction[Node2]==DOWN)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left)/2)*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top))*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x,LineArray[i].End.y+1);
       P2.y+=1;
      }
    else if(Element[Obj2]->Direction[Node2]==LEFT)
      {
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x-1,LineArray[i].End.y);
       P2.x-=1;
      }
    else if(Element[Obj2]->Direction[Node2]==RIGHT)
      {
       P2.x+=MyRound(double((Element[Obj2]->RectCoord[Node2].Right-Element[Obj2]->RectCoord[Node2].Left))*Mk);
       P2.y+=MyRound(double((Element[Obj2]->RectCoord[Node2].Bottom-Element[Obj2]->RectCoord[Node2].Top)/2)*Mk);
       LineArray[i].End.x=P2.x;
       LineArray[i].End.y=P2.y;
       Form1->Canvas->MoveTo(LineArray[i].End.x,LineArray[i].End.y);
       Form1->Canvas->LineTo(LineArray[i].End.x+1,LineArray[i].End.y);
       P2.x+=1;

      }
    Link2Points(P1,P2,Element[Obj1]->Direction[Node1],Element[Obj2]->Direction[Node2],max(Element[Obj1]->Width,Element[Obj2]->Width),max(Element[Obj1]->Height,Element[Obj2]->Height),is1Sheena,is2Sheena);


   }


for (int i=0;i<NumberOfElem;i++)
 if(RealNodes==false)
   {
    for(int j=0;j<Element[i]->NumberOfPoints;j++)
      if(Element[i]->PointNumber[j]!=-1)
        PrintNumber(Element[i],Element[i]->RectCoord[j],Element[i]->PointNumber[j],j);
   }
 else
   {
    for(int j=0;j<Element[i]->NumberOfPoints;j++)
      if(Element[i]->RealPointNumber[j]!=-1)
        PrintNumber(Element[i],Element[i]->RectCoord[j],Element[i]->RealPointNumber[j],j);
   }


}
//---------------------------------------------------------------------------
void TForm1::Link2Points(POINT P1,POINT P2,int Dir1,int Dir2,int MaxWidth,int MaxHeight,bool is1Sheena,bool is2Sheena)
{
POINT PA,PB,tmpP;
int tmpDir;
bool tmpFlag;
if(P1.x>P2.x)
  {
   tmpP=P1;
   P1=P2;
   P2=tmpP;
   tmpDir=Dir1;
   Dir1=Dir2;
   Dir2=tmpDir;
   tmpFlag=is1Sheena;
   is1Sheena=is2Sheena;
   is2Sheena=tmpFlag;
  }

if(((Dir1==LEFT||Dir1==RIGHT)&&is1Sheena==true)||(Dir1!=LEFT&&is1Sheena==false))
  {
   if((Dir2!=RIGHT&&is2Sheena==false)||((Dir2==RIGHT||Dir2==LEFT)&&is2Sheena==true))
     {
      Form1->Canvas->MoveTo(P1.x,P1.y);
      Form1->Canvas->LineTo(P1.x+(P2.x-P1.x)/2,P1.y);
      Form1->Canvas->LineTo(P1.x+(P2.x-P1.x)/2,P2.y);
      Form1->Canvas->LineTo(P2.x,P2.y);
     }
     else
     {
      Form1->Canvas->MoveTo(P1.x,P1.y);
      Form1->Canvas->LineTo(P2.x,P1.y);
      Form1->Canvas->LineTo(P2.x,P2.y);
     }
  }
else if(((Dir1==UP||Dir1==DOWN)&&is1Sheena==true)||(Dir1==LEFT&&is1Sheena==false))
  {
   if((Dir2!=RIGHT&&is2Sheena==false)||((Dir2==RIGHT||Dir2==LEFT)&&is2Sheena==true))
     {
      Form1->Canvas->MoveTo(P1.x,P1.y);
      Form1->Canvas->LineTo(P1.x,P2.y);
      Form1->Canvas->LineTo(P2.x,P2.y);
     }
   else
     {
      Form1->Canvas->MoveTo(P1.x,P1.y);
      Form1->Canvas->LineTo(P1.x,P1.y+(P2.y-P1.y)/2);
      Form1->Canvas->LineTo(P2.x,P1.y+(P2.y-P1.y)/2);
      Form1->Canvas->LineTo(P2.x,P2.y);
     }
  }

}

//---------------------------------------------------------------------------
void TForm1::BreakLinks()
{
int i=0;
while(i<LineCoords::LineNumber)
   {
    if(LineArray[i].Obj1Number==ImageForRightBtn->ObjectNumber)
       {
        ImageForRightBtn->PointIsLinked[LineArray[i].Node1Number]=false;
        Element[LineArray[i].Obj2Number]->PointIsLinked[LineArray[i].Node2Number]=false;
        for (int j=i;j<LineCoords::LineNumber-1;j++)
           {
            LineArray[j]=LineArray[j+1];
           }
        LineCoords::LineNumber--;
       }
    else if(LineArray[i].Obj2Number==ImageForRightBtn->ObjectNumber)
       {
        ImageForRightBtn->PointIsLinked[LineArray[i].Node2Number]=false;
        Element[LineArray[i].Obj1Number]->PointIsLinked[LineArray[i].Node1Number]=false;
        for (int j=i;j<LineCoords::LineNumber-1;j++)
           {
            LineArray[j]=LineArray[j+1];
           }
        LineCoords::LineNumber--;

       }
    else i++;
   }


}
//---------------------------------------------------------------------------

void __fastcall TForm1::N2Click(TObject *Sender)
{

Form2->ParamsNumber=ImageForRightBtn->NumberOfParams;
for (int i=0;i<MaxElementParams;i++)
  {
   Form2->ParamNames[i]=0;
   Form2->Params[i]=0;
  }
//   ShowMessage(ImageForRightBtn->Name);
  for (int i=0;i<ImageForRightBtn->NumberOfParams;i++)
  {
   Form2->ParamNames[i]=ImageForRightBtn->ParamNames[i];
   Form2->Params[i]=ImageForRightBtn->Params[i];

  }
Form2->ElementName=ImageForRightBtn->ElementName;
Form2->ParFromFile=ImageForRightBtn->ParFromFile;
Form2->ParFileName=ImageForRightBtn->ParFileName;
Form2->ParamsCheckEnabled=ImageForRightBtn->OtherParams;
Form2->ShowModal();
if(ImageForRightBtn->ParFromFile==true)
  {
   CyclePar[CyclePar[0].Quantity].Add(ImageForRightBtn->ObjectNumber,ImageForRightBtn->ParFileName);
   if(CyclePar[CyclePar[0].Quantity-1].NumberOfIters>MaxNumberOfIt)
      MaxNumberOfIt=CyclePar[CyclePar[0].Quantity-1].NumberOfIters;
  }

}
//---------------------------------------------------------------------------
int TForm1::ReverseDir(int Dir)
{
int rev;
switch (Dir)
{
case UP:

   rev=DOWN;
   break;

case DOWN:

//   return UP;
   rev=UP;
   break;

case RIGHT:

//   return LEFT;
   rev=LEFT;
   break;

case LEFT:

//   return RIGHT;
   rev=RIGHT;
   break;
}
return rev;
}
//---------------------------------------------------------------------------





void __fastcall TForm1::N3Click(TObject *Sender)
{
BreakLinks();
int ObjNo=ImageForRightBtn->ObjectNumber;
    ClickedObject=NULL;
    SelectedObject=NULL;
    tmpImage=NULL;
ImageForRightBtn=NULL;
DeleteElement(ObjNo,NumberOfElem);
Redraw();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N1Click(TObject *Sender)
{
BreakLinks();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N5Click(TObject *Sender)
{
ImageForRightBtn=SelectedObject;
N3Click(Sender);

}
//---------------------------------------------------------------------------

void __fastcall TForm1::N6Click(TObject *Sender)
{
ImageForRightBtn=SelectedObject;
BreakLinks();

}
//---------------------------------------------------------------------------

void __fastcall TForm1::N7Click(TObject *Sender)
{
ImageForRightBtn=SelectedObject;
N2Click(Sender);
}
//---------------------------------------------------------------------------

void __fastcall TForm1::MainMenu1Change(TObject *Sender, TMenuItem *Source,
      bool Rebuild)
{
if(SelectedObject==NULL)
  {
   N5->Enabled=false;
   N6->Enabled=false;
   N7->Enabled=false;
//   N11->Enabled=false;
  }
else
  {
   N5->Enabled=true;
   N6->Enabled=true;
   N7->Enabled=true;

   //N11->Enabled=true;
  }
if(MainElementNo!=-1)
  N8->Enabled=true;
else
  N8->Enabled=false;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N8Click(TObject *Sender)
{
int Quantity=NodeNum[0].Quantity;
if(Quantity>0)
{
for (int i=0;i<Quantity;i++)
  {
//   delete NodeNum[i];
  }
delete []NodeNum;
NodeNum=NULL;
}
for(int i=0;i<NumberOfElem;i++)
   {
    Element[i]->Picture=Image[Element[i]->ElementNumber]->Picture;
    Picture[i]->Picture=Image[Element[i]->ElementNumber]->Picture;
   }

Numerate();
Redraw();
}
//---------------------------------------------------------------------------



void __fastcall TForm1::N9Click(TObject *Sender)
{
if(MainElementNo!=-1&&MainElementNo!=ImageForRightBtn->ObjectNumber)
   {
    Element[MainElementNo]->Picture=Image[Element[MainElementNo]->ElementNumber]->Picture;
    Picture[MainElementNo]->Picture=Image[Element[MainElementNo]->ElementNumber]->Picture;
    Element[MainElementNo]->MainRectNo=-1;
   }
ImageForRightBtn->Picture=Image[ImageForRightBtn->ElementNumber]->Picture;
ImageForRightBtn->Canvas->Brush->Color=clWhite;
ImageForRightBtn->Canvas->Font->Color=clTeal;
ImageForRightBtn->Canvas->Font->Name="Times New Roman";
ImageForRightBtn->Canvas->Font->Size=6;
int x,y;

if(ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top<5)
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Bottom;
else if((ImageForRightBtn->Height-ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Bottom)<5)
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top-abs(ImageForRightBtn->Canvas->Font->Height);
else
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top;

if(ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left<5)
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Right;
else if((ImageForRightBtn->Width-ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Right)<5)
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left-6;
else
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left;
ImageForRightBtn->Canvas->Pen->Color=clTeal;
ImageForRightBtn->Canvas->TextOut(x,y,"�");


Picture[ImageForRightBtn->ObjectNumber]->Picture=ImageForRightBtn->Picture;
ImageForRightBtn->Canvas->Brush->Color=clBlack;
MainElementNo=ImageForRightBtn->ObjectNumber;
}
//---------------------------------------------------------------------------
void TForm1::PrintNumber(MyImage *Object,TRect Rect,int Number,int RectNumber)
{

//Object->Canvas->Font->Name=
int x,y;
Object->Picture=Picture[Object->ObjectNumber]->Picture;
//Object->Picture=Image[Object->ElementNumber]->Picture;
Object->Canvas->Font->Name="Times New Roman";
Object->Canvas->Font->Height=12;
Object->Canvas->Font->Color=clPurple;
Object->Canvas->Font->Color=clRed;
Object->Canvas->Brush->Style=bsClear;
Object->Canvas->Pen->Color=clBlack;
Form1->Canvas->Font->Name="Times New Roman";
Form1->Canvas->Font->Height=12;
Form1->Canvas->Font->Color=clPurple;
Form1->Canvas->Font->Color=clRed;
Form1->Canvas->Brush->Style=bsClear;
Form1->Canvas->Pen->Color=clBlack;

if(Object->Direction[RectNumber]==UP)
  {
   x=Object->Left+MyRound(double(Rect.Right)*Mk);
   y=Object->Top-Form1->Canvas->Font->Height;
  }
else if(Object->Direction[RectNumber]==DOWN)
  {
   x=Object->Left+MyRound(double(Rect.Right)*Mk);
   y=Object->Top+Object->Height;
  }
else if(Object->Direction[RectNumber]==LEFT)
  {
   x=Object->Left-Form1->Canvas->Font->Height;
   y=Object->Top+MyRound(double(Rect.Top)*Mk)-Form1->Canvas->Font->Height/2;
  }
else if(Object->Direction[RectNumber]==RIGHT)
  {
   x=Object->Left+Object->Width;
   y=Object->Top+MyRound(double(Rect.Top)*Mk)-Form1->Canvas->Font->Height/2;
  }


Form1->Canvas->TextOut(x,y,Number);

Object->Canvas->Brush->Color=clBlack;
Object->Canvas->Pen->Color=clBlack;
Picture[Object->ObjectNumber]->Picture=Object->Picture;

}
//---------------------------------------------------------------------------
void TForm1::Numerate()
{
int i=1,j=0;
int CurrNo=0;
int StartObj,FinishObj,FinishPoint;
bool ForwardFlag=true;
bool finish=false;

for(i=0;i<NumberOfElem;i++)
  {
   for(j=0;j<Element[i]->NumberOfPoints;j++)
      {
       Element[i]->PointNumber[j]=-1;
       Element[i]->PointIsForwarded[j]=false;
      }
  }
StartObj=MainElementNo;


while (finish==false)
  {
   if(ForwardFlag==true)
   {
    ForwardFlag=false;
//    for (j=0;j<NodeNum[StepNo].NumberOfNodes;j++)
    for (j=0;j<Element[StartObj]->NumberOfPoints;j++)
      if(Element[StartObj]->PointNumber[j]==-1||Element[StartObj]->ElementName=="����")
        {
         if(Element[StartObj]->ElementName!="����")
            {
             CurrNo++;
             Element[StartObj]->PointNumber[j]=CurrNo;
             if(RealNodes==false)
               PrintNumber(Element[StartObj],Element[StartObj]->RectCoord[j],CurrNo,j);
            }
         //NodesNum::Number++;
         if(FindLink(StartObj,j,FinishObj,FinishPoint)==true)
           {
            if(Element[FinishObj]->PointNumber[FinishPoint]==-1)
              {
               if(Element[FinishObj]->ElementName=="����")
                  for(int k=0;k<Element[FinishObj]->NumberOfPoints;k++)
                    {
                     Element[FinishObj]->PointNumber[k]=Element[StartObj]->PointNumber[j];
                     if(RealNodes==false)
                       PrintNumber(Element[FinishObj],Element[FinishObj]->RectCoord[k],Element[StartObj]->PointNumber[j],k);
                    }
               else
                  {
                   Element[FinishObj]->PointNumber[FinishPoint]=Element[StartObj]->PointNumber[j];
                   if(RealNodes==false)
                     PrintNumber(Element[FinishObj],Element[FinishObj]->RectCoord[FinishPoint],Element[StartObj]->PointNumber[j],FinishPoint);
                  }

              }
           }
         else
           {
            Element[StartObj]->PointIsForwarded[j]=true;
           }
        }

    for (i=0;i<Element[StartObj]->NumberOfPoints;i++)
     {
      if(Element[StartObj]->PointIsForwarded[i]==false&&FindLink(StartObj,i,FinishObj,FinishPoint)==true)
        {
         Element[StartObj]->PointIsForwarded[i]=true;
         Element[FinishObj]->PointIsForwarded[FinishPoint]=true;
         Element[FinishObj]->InIndex++;
         Element[FinishObj]->InPointNo[Element[FinishObj]->InIndex]=FinishPoint;

         ForwardFlag=true;
         StartObj=FinishObj;
//         StartPoint=FinishPoint;
         goto cyc1End;
        }

     }
     cyc1End:
   }//forward end
   else if(ForwardFlag==false)//back
   {
//    flag=0;

    if(FindLink(StartObj,Element[StartObj]->InPointNo[Element[StartObj]->InIndex],FinishObj,FinishPoint)==true)
        {
         Element[StartObj]->InIndex--;
         StartObj=FinishObj;
//         StartPoint=FinishPoint;
         ForwardFlag=true;
        }
    else
        finish=true;
   }

  }
if(BalanceElementNo!=-1)
  {
   ImageForRightBtn=Element[BalanceElementNo];
   N16Click(N16);
   for(i=0;i<ImageForRightBtn->NumberOfPoints;i++)
//     if(ImageForRightBtn->BalanceRectNo!=i)
        if(RealNodes==false)
          PrintNumber(ImageForRightBtn,ImageForRightBtn->RectCoord[i],ImageForRightBtn->PointNumber[i],i);

  }
}
//---------------------------------------------------------------------------
bool TForm1::FindLink(int ObjStart,int NodeStart,int &ObjFinish,int &NodeFinish)
{
for (int i=0;i<LineCoords::LineNumber;i++)
  {
   if (ObjStart==LineArray[i].Obj1Number&&NodeStart==LineArray[i].Node1Number)
      {
       ObjFinish=LineArray[i].Obj2Number;
       NodeFinish=LineArray[i].Node2Number;
       return true;
      }
   else if (ObjStart==LineArray[i].Obj2Number&&NodeStart==LineArray[i].Node2Number)
      {
       ObjFinish=LineArray[i].Obj1Number;
       NodeFinish=LineArray[i].Node1Number;
       return true;
      }
  }
return false;
}
//---------------------------------------------------------------------------
NodesNum::NodesNum()
{
for (int i=0;i<MaxNodesInElem;i++)
  {
   NodeIsForwarded[i]=false;
   NodeNumber[i]=-1;

  }
ObjNumber=-1;
NumberOfNodes-1;
InNodeNo=-1;
OutNode1No=-1;
OutNode2No=-1;
Number=Quantity+1;
Quantity++;
}
//---------------------------------------------------------------------------
NodesNum::~NodesNum()
{
Quantity--;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::N11Click(TObject *Sender)
{
FormR->Show();
}
//---------------------------------------------------------------------------
void TForm1::Save(AnsiString FName)
{
int tmpInt;
int i,k;
int FHandle;
char tmpChar[100];
char *ver="Ver";

DeleteFile(FName);
FHandle=FileCreate(FName);

FileWrite(FHandle,&("Ver"),sizeof("Ver"));
FileWrite(FHandle,&(VerNo),sizeof(VerNo));
FileWrite(FHandle,&(NumberOfElem),sizeof(NumberOfElem));
FileWrite(FHandle,&(LineCoords::LineNumber),sizeof(LineCoords::LineNumber));
FileWrite(FHandle,&(MainElementNo),sizeof(MainElementNo));
FileWrite(FHandle,&(BalanceElementNo),sizeof(BalanceElementNo));
FileWrite(FHandle,&(precision),sizeof(precision));
FileWrite(FHandle,&(RealNodes),sizeof(RealNodes));
FileWrite(FHandle,&(Mk),sizeof(Mk));

for (int j=0;j<NumberOfElem;j++)
{
FileWrite(FHandle,&(Element[j]->Left),sizeof(Element[j]->Left));
FileWrite(FHandle,&(Element[j]->Top),sizeof(Element[j]->Top));
FileWrite(FHandle,&(Element[j]->NumberOfPoints),sizeof(Element[j]->NumberOfPoints));
FileWrite(FHandle,&(Element[j]->ElementNumber),sizeof(Element[j]->ElementNumber));
FileWrite(FHandle,&(Element[j]->type),sizeof(Element[j]->type));


for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->RectCoord[i]),sizeof(Element[j]->RectCoord[i]));
for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->Direction[i]),sizeof(Element[j]->Direction[i]));
for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->PointIsLinked[i]),sizeof(Element[j]->PointIsLinked[i]));
for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->PointIsForwarded[i]),sizeof(Element[j]->PointIsForwarded[i]));
for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->PointNumber[i]),sizeof(Element[j]->PointNumber[i]));
for(i=0;i<Element[j]->NumberOfPoints;i++)
  FileWrite(FHandle,&(Element[j]->RealPointNumber[i]),sizeof(Element[j]->RealPointNumber[i]));

//FileWrite(FHandle,&(Element[j]->ObjectNumber),sizeof(Element[j]->ObjectNumber));
FileWrite(FHandle,&(Element[j]->NumberOfParams),sizeof(Element[j]->NumberOfParams));
for(i=0;i<Element[j]->NumberOfParams;i++)
  FileWrite(FHandle,&(Element[j]->Params[i]),sizeof(Element[j]->Params[i]));

for(i=0;i<Element[j]->NumberOfParams;i++)
  {
   for(k=0;k<Element[j]->ParamNames[i].Length();k++)
     {
      tmpChar[k]=Element[j]->ParamNames[i][k+1];

     }
//   tmpChar[k]='\n';
//   tmpInt=Element[j]->ParamNames[i].Length()+1;
   tmpInt=Element[j]->ParamNames[i].Length();
   FileWrite(FHandle,&(tmpInt),sizeof(tmpInt));
   FileWrite(FHandle,&(tmpChar),tmpInt);

  }
FileWrite(FHandle,&(Element[j]->MainRectNo),sizeof(Element[j]->MainRectNo));
FileWrite(FHandle,&(Element[j]->BalanceRectNo),sizeof(Element[j]->BalanceRectNo));
//tmpChar=Element[j]->ElementName.c_str();
for(k=0;k<Element[j]->ElementName.Length();k++)
     {
      tmpChar[k]=Element[j]->ElementName[k+1];

     }
//tmpChar[k]='\n';
//tmpInt=Element[j]->ElementName.Length()+1;
tmpInt=Element[j]->ElementName.Length();
FileWrite(FHandle,&(tmpInt),sizeof(tmpInt));

FileWrite(FHandle,&(tmpChar),tmpInt);
FileWrite(FHandle,&(Element[j]->InPointNo[0]),sizeof(Element[j]->InPointNo[0]));
FileWrite(FHandle,&(Element[j]->InPointNo[1]),sizeof(Element[j]->InPointNo[1]));
FileWrite(FHandle,&(Element[j]->InIndex),sizeof(Element[j]->InIndex));
FileWrite(FHandle,&(Element[j]->OtherParams),sizeof(Element[j]->OtherParams));
}
for (int j=0;j<LineCoords::LineNumber;j++)
{
FileWrite(FHandle,&(LineArray[j].Start),sizeof(LineArray[j].Start));
FileWrite(FHandle,&(LineArray[j].End),sizeof(LineArray[j].End));
FileWrite(FHandle,&(LineArray[j].Obj1Number),sizeof(LineArray[j].Obj1Number));
FileWrite(FHandle,&(LineArray[j].Obj2Number),sizeof(LineArray[j].Obj2Number));
FileWrite(FHandle,&(LineArray[j].Node1Number),sizeof(LineArray[j].Node1Number));
FileWrite(FHandle,&(LineArray[j].Node2Number),sizeof(LineArray[j].Node2Number));
}
FileClose(FHandle);




}
//---------------------------------------------------------------------------




void __fastcall TForm1::N13Click(TObject *Sender)
{
MenuSave();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::N14Click(TObject *Sender)
{
MenuSaveAs();

}
//---------------------------------------------------------------------------
void TForm1::Open(AnsiString FName)
{
int i;
int FHandle=FileOpen(FName,fmOpenRead);
int NumberOfPoints,ElementNumber,NumberOfLinks,ElementsRead,
    tmpInt,PointNumber[MaxNodesInElem],RealPointNumber[MaxNodesInElem],InPointNo[2],MainRectNo,InIndex,BalanceRectNo;
int VerNo;
bool UseVerNo=false;
bool ParamsType;
char tmpChar[100];
char *tmphar;
bool PointIsLinked[MaxNodesInElem],PointIsForwarded[MaxNodesInElem];
POINT Pos;
for(i=0;i<NumberOfElem;i++)
  {
   delete Element[i];
   Element[i]=NULL;
  }

SelectedObject=NULL;
NumberOfElem=0;
MouseBtnLeftIsDown=false;
ObjectIsClicked=false;
MainElementNo=-1;
BalanceElementNo=-1;
MaxNumberOfIt=1;

if(LineCoords::LineNumber>0)
{
  delete []LineArray;
  LineArray=NULL;
}
LinesNumber=0;

tmpInt=3;
FileRead(FHandle,&(tmpChar),sizeof("Ver"));
if(tmpChar[0]=='V'&&tmpChar[1]=='e'&&tmpChar[2]=='r')
 {
  UseVerNo=true;
  FileRead(FHandle,&(VerNo),sizeof(VerNo));
 }
else
 FileSeek(FHandle,0,0);


FileRead(FHandle,&(ElementsRead),sizeof(ElementsRead));
FileRead(FHandle,&(NumberOfLinks),sizeof(NumberOfLinks));
FileRead(FHandle,&(MainElementNo),sizeof(MainElementNo));
FileRead(FHandle,&(BalanceElementNo),sizeof(BalanceElementNo));
FileRead(FHandle,&(precision),sizeof(precision));
FileRead(FHandle,&(RealNodes),sizeof(RealNodes));
FileRead(FHandle,&(Mk),sizeof(Mk));

NumberOfElem=0;
for (int j=0;j<ElementsRead;j++)
{
FileRead(FHandle,&(Pos.x),sizeof(Pos.x));
FileRead(FHandle,&(Pos.y),sizeof(Pos.y));
FileRead(FHandle,&NumberOfPoints,sizeof(NumberOfPoints));
FileRead(FHandle,&ElementNumber,sizeof(ElementNumber));
FileRead(FHandle,&(type),sizeof(type));


for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(NodesCoord[i]),sizeof(NodesCoord[i]));
for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(Dirs[i]),sizeof(Dirs[i]));
for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(PointIsLinked[i]),sizeof(PointIsLinked[i]));
for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(PointIsForwarded[i]),sizeof(PointIsForwarded[i]));
for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(PointNumber[i]),sizeof(PointNumber[i]));
for(i=0;i<NumberOfPoints;i++)
  FileRead(FHandle,&(RealPointNumber[i]),sizeof(RealPointNumber[i]));


//FileWrite(FHandle,&(Element[j]->ObjectNumber),sizeof(Element[j]->ObjectNumber));
FileRead(FHandle,&(NumberOfElemParams),sizeof(NumberOfElemParams));
for(i=0;i<NumberOfElemParams;i++)
  FileRead(FHandle,&(Params[i]),sizeof(Params[i]));

for(i=0;i<NumberOfElemParams;i++)
  {
//   tmpChar=Element[j]->ParamNames[i].c_str();
//   tmpInt=Element[j]->ParamNames[i].Length();
   FileRead(FHandle,&(tmpInt),sizeof(tmpInt));
   FileRead(FHandle,&(tmpChar),tmpInt);
   ParamNames[i]="";
   for(int k=0;k<tmpInt;k++)
      ParamNames[i]+=tmpChar[k];
      //ParamNames[i][k+1]=tmpChar[k];
  }

FileRead(FHandle,&(MainRectNo),sizeof(MainRectNo));
FileRead(FHandle,&(BalanceRectNo),sizeof(BalanceRectNo));
//tmpChar=Element[j]->ElementName.c_str();
//tmpInt=Element[j]->ElementName.Length();
FileRead(FHandle,&(tmpInt),sizeof(tmpInt));
FileRead(FHandle,&(tmpChar),tmpInt);
ElementName="";
for(int k=0;k<tmpInt;k++)
  ElementName+=tmpChar[k];

//ElementName=AnsiString(tmpChar);
FileRead(FHandle,&(InPointNo[0]),sizeof(InPointNo[0]));
FileRead(FHandle,&(InPointNo[1]),sizeof(InPointNo[1]));
FileRead(FHandle,&(InIndex),sizeof(InIndex));

if(UseVerNo==true)
 FileRead(FHandle,&(ParamsType),sizeof(ParamsType));

PicToCreate=Image[ElementNumber]->Picture;

AddElement(ElementNumber,ElementName,NumberOfElem,
PicToCreate,Pos,NumberOfPoints,
NodesCoord,NumberOfElemParams,ParamNames,Params,Dirs);
Element[j]->InIndex=InIndex;
Element[j]->InPointNo[0]=InPointNo[0];
Element[j]->InPointNo[1]=InPointNo[1];
Element[j]->MainRectNo=MainRectNo;
Element[j]->BalanceRectNo=BalanceRectNo;
if(UseVerNo==true)
 Element[j]->OtherParams=ParamsType;

for(i=0;i<NumberOfPoints;i++)
   {
    Element[j]->PointNumber[i]=PointNumber[i];
    Element[j]->RealPointNumber[i]=RealPointNumber[i];
    Element[j]->PointIsLinked[i]=PointIsLinked[i];
    Element[j]->PointIsForwarded[i]=PointIsForwarded[i];
   }
}
LineArray=new LineCoords[NumberOfLinks];
LineCoords::LineNumber=NumberOfLinks;
for (int j=0;j<NumberOfLinks;j++)
{
FileRead(FHandle,&(LineArray[j].Start),sizeof(LineArray[j].Start));
FileRead(FHandle,&(LineArray[j].End),sizeof(LineArray[j].End));
FileRead(FHandle,&(LineArray[j].Obj1Number),sizeof(LineArray[j].Obj1Number));
FileRead(FHandle,&(LineArray[j].Obj2Number),sizeof(LineArray[j].Obj2Number));
FileRead(FHandle,&(LineArray[j].Node1Number),sizeof(LineArray[j].Node1Number));
FileRead(FHandle,&(LineArray[j].Node2Number),sizeof(LineArray[j].Node2Number));
}

int minTop=99999999999,minLeft=99999999999;

for (int j=0;j<NumberOfElem;j++)
{
 if(Element[j]->Top<minTop)
  minTop=Element[j]->Top;
 if(Element[j]->Left<minLeft)
  minLeft=Element[j]->Left;
}

for (int j=0;j<NumberOfElem;j++)
{
 if(minTop<0)
  Element[j]->Top-=minTop;
 if(minLeft<0)
  Element[j]->Left-=minLeft;
}
if(MainElementNo>=0)
   {
    ImageForRightBtn=Element[MainElementNo];
    N9Click(N9);
   }
if(RealNodes==true)
  N34->Checked=true;
else
  N35->Checked=true;

Redraw();
FileClose(FHandle);
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N15Click(TObject *Sender)
{
MenuOpen();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N16Click(TObject *Sender)
{
if(BalanceElementNo!=-1&&BalanceElementNo!=ImageForRightBtn->ObjectNumber)
   {
    Element[BalanceElementNo]->Picture=Image[Element[BalanceElementNo]->ElementNumber]->Picture;
    Picture[BalanceElementNo]->Picture=Image[Element[BalanceElementNo]->ElementNumber]->Picture;
    Element[BalanceElementNo]->BalanceRectNo=-1;


   }
int x,y;
TRect Rect;
ImageForRightBtn->Picture=Image[ImageForRightBtn->ElementNumber]->Picture;
ImageForRightBtn->Canvas->Brush->Color=clYellow;
ImageForRightBtn->Canvas->Font->Name="Times New Roman";
ImageForRightBtn->Canvas->Font->Size=6;
ImageForRightBtn->Canvas->Font->Color=clYellow;

if(ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top<5)
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Bottom;
//  y=2;
else if((ImageForRightBtn->Height-ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Bottom)<5)
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top-abs(ImageForRightBtn->Canvas->Font->Height);
//  y=ImageForRightBtn->Height-abs(ImageForRightBtn->Canvas->Font->Height)-2;
else
  y=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Top;
if(ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left<5)
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Right;
else if((ImageForRightBtn->Width-ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Right)<5)
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left-6;
else
  x=ImageForRightBtn->RectCoord[ImageForRightBtn->MainRectNo].Left;
ImageForRightBtn->Canvas->Pen->Color=clYellow;

ImageForRightBtn->Canvas->RoundRect(x-1,y-1,x+8,y+abs(ImageForRightBtn->Canvas->Font->Height),2,2);
ImageForRightBtn->Canvas->Font->Color=clBlack;
ImageForRightBtn->Canvas->TextOut(x,y,"�");

Picture[ImageForRightBtn->ObjectNumber]->Picture=ImageForRightBtn->Picture;
ImageForRightBtn->Canvas->Brush->Color=clBlack;
BalanceElementNo=ImageForRightBtn->ObjectNumber;
}
//---------------------------------------------------------------------------



void __fastcall TForm1::FormPaint(TObject *Sender)
{
RedrawWithoutRefresh();
}
//---------------------------------------------------------------------------

void TForm1::Prepare()//���������� �����, � ���, ������� ��� �������
{
int i,j,k;
MaxNodeNumber=0;
NumberOfTrf=0,NumberOfNodes=0;

int CurrTrf=0,Node1,Node2,mswij,mhui,mhuj,CurrStr=0,CurrPqu=0,CurrLine=0,CurrShReact=0;
int OldCurrStr=0;
int TipN,n;
double Uk,Un,Sn,Xt,X0,R0,Ktr,Ktr_i,Pk,Px,Ix,Pn,Qn,Ed,Cos,Sin,a0,a1,a2,b0,b1,b2,KPD;
double Xd,Xg,Pr,Qr,PZag,Delt,Eg;
double ie,Ubaz,Cktr,Cktm,Yr,Ym,G0,B0,Yr2,Ym2;
bool flag=true;
bool paramtype=false;
/*
for(i=0;i<CyclePar[0].Quantity;i++)
{
  for(j=0;j<MaxElementParams;j++)
   Element[CyclePar[i].ElementNo]->Params[j]=CyclePar[i].Params[j][CurrIt];
}
*/
for(i=0;i<NumberOfElem;i++)
  {
   for (j=0;j<Element[i]->NumberOfPoints;j++)
    if(Element[i]->PointNumber[j]>MaxNodeNumber)
       MaxNodeNumber=Element[i]->PointNumber[j];
   if((Element[i]->type==1||Element[i]->type==2||Element[i]->type==3)&&(Element[i]->PointNumber[0]!=-1))
    {
     LinesNumber+=3;
     NumberOfTrf+=3;
    }
   if(Element[i]->NumberOfPoints==2&&Element[i]->PointNumber[0]!=-1)
     LinesNumber++;
   if(Element[i]->NumberOfPoints==3&&Element[i]->PointNumber[0]!=-1)
    NumberOfNodes+=4;
   else if(Element[i]->type==9||Element[i]->type==11)
    {
     NumberOfNodes+=2;
     LinesNumber++;
    }
   else
    NumberOfNodes+=Element[i]->NumberOfPoints;
   if((Element[i]->type==4)&&(Element[i]->PointNumber[0]!=-1))
    NumberOfTrf++;
  }
//LinesNumber=NumberOfNodes;
TrfList=new struct Trf[NumberOfTrf+1];
PquList=new struct Pqu[NumberOfNodes+1];
StrList=new struct Str[NumberOfNodes*NumberOfNodes+1];
LinList=new struct Lin[LinesNumber+1];
ShReactList= new struct ShReact[NumberOfNodes+1];
for(i=0;i<NumberOfTrf+1;i++)
  {TrfList[i].NodeiNumber=TrfList[i].NodejNumber=TrfList[i].TipN=
  TrfList[i].NodeiNumber=TrfList[i].ktr=TrfList[i].un=TrfList[i].sn=
  TrfList[i].pk=TrfList[i].uk=TrfList[i].px=TrfList[i].ix=TrfList[i].xt=TrfList[i].n=0;
  TrfList[i].ktr_i=TrfList[i].r=TrfList[i].x=TrfList[i].g=TrfList[i].b=0;
  TrfList[i].paramtype=false;
  }
for(i=0;i<NumberOfNodes+1;i++)
 {
  PquList[i].NodeNumber=PquList[i].Hui=PquList[i].TipN=PquList[i].P=PquList[i].Pg=PquList[i].dPg=PquList[i].dPn=PquList[i].Q
  =PquList[i].Qg=PquList[i].dQn=PquList[i].Un=PquList[i].a0=PquList[i].a1=PquList[i].a2=PquList[i].b0=PquList[i].b1
  =PquList[i].b2=PquList[i].Qmin=PquList[i].Qmax=0;
  ShReactList[i].NodeNumber=-1;
  ShReactList[i].B=ShReactList[i].G=0;
 }
for(i=0;i<NumberOfNodes*NumberOfNodes+1;i++)
  StrList[i].NodeiNumber=StrList[i].NodejNumber=StrList[i].Mswij=StrList[i].Mhui
  =StrList[i].Mhuj=StrList[i].cktr=StrList[i].cktm=StrList[i].yr=StrList[i].ym=StrList[i].flag=0;
for(i=0;i<LinesNumber+1;i++)
  LinList[i].NodeiNumber=LinList[i].NodejNumber=LinList[i].TipN=LinList[i].r0=LinList[i].x0
  =LinList[i].b0=LinList[i].Length=LinList[i].n=0;
struct Str tmpStrList;
struct Lin tmpLinList;
struct Pqu tmpPquList;
struct Trf tmpTrfList;
struct ShReact tmpShReactList;
//�������������� ��-��
for(i=0;i<NumberOfElem;i++)
  if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==1||Element[i]->type==2))
    {
     MaxNodeNumber++;
     paramtype=Element[i]->OtherParams;
     for(j=0;j<3;j++)
     {
     switch (j)
     {
     case 0:

//        Node1=Element[i]->PointNumber[0];
        Node2=Element[i]->PointNumber[0];
        if(paramtype==false)
        {Uk=(Element[i]->Params[6]+Element[i]->Params[7]-Element[i]->Params[8])/2.;
        Xt=Uk*Element[i]->Params[1]*Element[i]->Params[1]/(100.*Element[i]->Params[0]);
        Px=Element[i]->Params[9];
        Ix=Element[i]->Params[10];
        }
        else
        {
        Yr=Element[i]->Params[12];
        Ym=Element[i]->Params[13];
        G0=Element[i]->Params[22];
        B0=Element[i]->Params[23];
        }
        Ktr=1.000001;Ktr_i=0;
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[1])
            {

             if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                mhui=0;
             else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
                mhui=0;
             else if(StruList[k].u5<10)
                mhui=StruList[k].u5;
             else
                mhui=1;
            }
        break;
     case 1:
        //Node1=Element[i]->PointNumber[1];
        Node2=Element[i]->PointNumber[1];
        if(paramtype==false)
        {Uk=(Element[i]->Params[6]-Element[i]->Params[7]+Element[i]->Params[8])/2.;
        Xt=Uk*max(Element[i]->Params[2],Element[i]->Params[1])*max(Element[i]->Params[2],Element[i]->Params[1])/(100.*Element[i]->Params[0]);
        Ktr=Element[i]->Params[2]/Element[i]->Params[1];
        Ktr_i=0;
        Px=0;
        Ix=0;}
        else
        {
        Yr=Element[i]->Params[14];
        Ym=Element[i]->Params[15];
        G0=0;
        B0=0;
        Ktr=real(complex<double>(1,0)/complex<double>(Element[i]->Params[16],Element[i]->Params[17]));
        Ktr_i=imag(complex<double>(1,0)/complex<double>(Element[i]->Params[16],Element[i]->Params[17]));
        }
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[1]&&StruList[k].j==Element[i]->PointNumber[2])
            {

             if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                mhui=0;
             else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
                mhui=0;
             else if(StruList[k].u5<10)
                mhui=StruList[k].u5;
             else
                mhui=1;
            }

        break;
     case 2:
        //Node1=Element[i]->PointNumber[2];
        Node2=Element[i]->PointNumber[2];
        if(paramtype==false)
        {
        Uk=(-Element[i]->Params[6]+Element[i]->Params[7]+Element[i]->Params[8])/2.;
        Xt=Uk*max(Element[i]->Params[3],Element[i]->Params[1])*max(Element[i]->Params[3],Element[i]->Params[1])/(100.*Element[i]->Params[0]);
        Ktr=Element[i]->Params[3]/Element[i]->Params[1];
        Ktr_i=0;
        Px=0;
        Ix=0;}
        else
        {
        Yr=Element[i]->Params[18];
        Ym=Element[i]->Params[19];
        G0=0;
        B0=0;
//        Ktr=Element[i]->Params[16];
//        Ktr_i=Element[i]->Params[17];
        Ktr=real(complex<double>(1,0)/complex<double>(Element[i]->Params[20],Element[i]->Params[21]));
        Ktr_i=imag(complex<double>(1,0)/complex<double>(Element[i]->Params[20],Element[i]->Params[21]));
        }
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[2])
            {

             if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
                mhui=0;
             else if((StruList[k].u6>1&&StruList[k].u6!=5)&&kzMode==true)
                mhui=0;
             else if(StruList[k].u6<10)
                mhui=StruList[k].u6;
             else
                mhui=1;
            }

        break;
     }
     mhuj=1;
     mswij=Element[i]->type*10+j;
     //Node2=MaxNodeNumber;
     Node1=MaxNodeNumber;
     Pk=Element[i]->Params[4]/2.0;//Pk
     //uk=ukws+ukwn-uksn
     //Xt=Uk*Uwn/100*Sn

     TrfList[CurrTrf].NodeiNumber=Node1;
     TrfList[CurrTrf].NodejNumber=Node2;
    if(Ktr>1.000001)
       {
        TrfList[CurrTrf].NodeiNumber=Node2;
        TrfList[CurrTrf].NodejNumber=Node1;
        Ktr=real(complex<double>(1,0)/complex<double>(Ktr,Ktr_i));
        Ktr_i=imag(complex<double>(1,0)/complex<double>(Ktr,Ktr_i));
//        TrfList[CurrTrf].un=Element[i]->Params[2];
        //Xt=Element[i]->Params[4]*Element[i]->Params[2]*Element[i]->Params[2]/(100.*Element[i]->Params[0]);
        StrList[CurrStr].NodeiNumber=Node2;
        StrList[CurrStr].NodejNumber=Node1;
        StrList[CurrStr].Mhui=mhuj;
        StrList[CurrStr].Mhuj=mhui;

       }

     TrfList[CurrTrf].ktr=Ktr;
     TrfList[CurrTrf].ktr_i=Ktr_i;
/*     for (int mm=0;mm<3;mm++)
        if(Element[i]->Params[mm+1]>TrfList[CurrTrf].un)
           TrfList[CurrTrf].un=Element[i]->Params[mm+1];
*/
     if(paramtype==false)
     {TrfList[CurrTrf].un=max(Element[i]->Params[1],max(Element[i]->Params[2],Element[i]->Params[3]));
     TrfList[CurrTrf].sn=Element[i]->Params[0];
     TrfList[CurrTrf].pk=Pk;
     TrfList[CurrTrf].uk=Uk;
     TrfList[CurrTrf].px=Px;
     TrfList[CurrTrf].ix=Ix;
     TrfList[CurrTrf].xt=Xt;
     TrfList[CurrTrf].n=Element[i]->Params[11];
     }
     else
     {
     TrfList[CurrTrf].r=Yr;
     TrfList[CurrTrf].x=Ym;
     TrfList[CurrTrf].b=B0;
     TrfList[CurrTrf].g=G0;
     }
     TrfList[CurrTrf].TipN=mswij;
     TrfList[CurrTrf].paramtype=paramtype;


     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=mswij;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=mhuj;

     PquList[CurrPqu].NodeNumber=Node2;
     PquList[CurrPqu].Hui=mhui;
     if(paramtype==false)
       PquList[CurrPqu].Un=Element[i]->Params[1]*Ktr;
     PquList[CurrPqu].TipN=0;
     CurrTrf++;
     CurrStr++;
     CurrPqu++;
     }
    PquList[CurrPqu].NodeNumber=Node1;
    PquList[CurrPqu].Hui=mhuj;
    if(paramtype==false)
      PquList[CurrPqu].Un=Element[i]->Params[1];
    PquList[CurrPqu].TipN=0;
    CurrPqu++;
     //if(
    }
//�������������� ��-�� � �������. ��������
for(i=0;i<NumberOfElem;i++)
  if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==3))
    {
     MaxNodeNumber++;
     paramtype=Element[i]->OtherParams;
     for(j=0;j<3;j++)
     {
     switch (j)
     {
     case 0:

        //Node1=Element[i]->PointNumber[0];
        Node2=Element[i]->PointNumber[0];
        if(paramtype==false)
        {
        Uk=0.125*(Element[i]->Params[4]);
        Xt=Uk*Element[i]->Params[1]*Element[i]->Params[1]/(100.*Element[i]->Params[0]);
        Px=Element[i]->Params[5];
        Ix=Element[i]->Params[6];
        }
        else
        {
        Yr=Element[i]->Params[8];
        Ym=Element[i]->Params[9];
        G0=Element[i]->Params[14];
        B0=Element[i]->Params[15];
        }
        Ktr=1.000001;
        Ktr_i=0;
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[1])
            {

             if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                mhui=0;
           //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
           //     mhui=0;
             else if(StruList[k].u5<10)
                mhui=StruList[k].u5;
             else
                mhui=1;
            }
        break;
     case 1:
        //Node1=Element[i]->PointNumber[1];
        Node2=Element[i]->PointNumber[1];
        if(paramtype==false)
        {
        Uk=1.75*(Element[i]->Params[4]);
        Xt=Uk*max(Element[i]->Params[1],Element[i]->Params[2])*max(Element[i]->Params[1],Element[i]->Params[2])/(100.*Element[i]->Params[0]);
        Ktr=Element[i]->Params[2]/Element[i]->Params[1];
        Ktr_i=0;
        Px=0;
        Ix=0;
        }
        else
        {
        Yr=Element[i]->Params[10];
        Ym=Element[i]->Params[11];
        G0=0;
        B0=0;
        Ktr=real(complex<double>(1,0)/complex<double>(Element[i]->Params[12],Element[i]->Params[13]));
        Ktr_i=imag(complex<double>(1,0)/complex<double>(Element[i]->Params[12],Element[i]->Params[13]));
        }
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[1]&&StruList[k].j==Element[i]->PointNumber[2])
            {

             if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                mhui=2;
           //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
           //     mhui=0;
             else if(StruList[k].u5<10)
                mhui=StruList[k].u5;
             else
                mhui=1;
            }
        break;
     case 2:
        //Node1=Element[i]->PointNumber[2];
        Node2=Element[i]->PointNumber[2];
        //Uk=(-Element[i]->Params[6]+Element[i]->Params[7]+Element[i]->Params[8])/2.;
        //Xt=Uk*Element[i]->Params[3]*Element[i]->Params[3]/(100.*Element[i]->Params[0]);
        //Ktr=Element[i]->Params[3]/Element[i]->Params[1];
        //Px=0;
        //Ix=0;
        for(k=0;k<NumberOfLines;k++)
          if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[2])
            {

             if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
                mhui=2;
           //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
           //     mhui=0;
             else if(StruList[k].u6<10)
                mhui=StruList[k].u6;
             else
                mhui=1;
            }

        break;
     }
     mhuj=1;
     mswij=Element[i]->type*10+j;
     //Node2=MaxNodeNumber;
     Node1=MaxNodeNumber;
     Pk=Element[i]->Params[3]/2.0;//Pk
     //uk=ukws+ukwn-uksn
     //Xt=Uk*Uwn/100*Sn
     TrfList[CurrTrf].NodeiNumber=Node1;
     TrfList[CurrTrf].NodejNumber=Node2;
    if(Ktr>1.000001)
       {
        TrfList[CurrTrf].NodeiNumber=Node2;
        TrfList[CurrTrf].NodejNumber=Node1;
        Ktr=real(complex<double>(1,0)/complex<double>(Ktr,Ktr_i));
        Ktr=imag(complex<double>(1,0)/complex<double>(Ktr,Ktr_i));        
//        TrfList[CurrTrf].un=Element[i]->Params[2];
        //Xt=Element[i]->Params[4]*Element[i]->Params[2]*Element[i]->Params[2]/(100.*Element[i]->Params[0]);
        StrList[CurrStr].NodeiNumber=Node2;
        StrList[CurrStr].NodejNumber=Node1;
        StrList[CurrStr].Mhui=mhuj;
        StrList[CurrStr].Mhuj=mhui;

       }
     TrfList[CurrTrf].ktr=Ktr;
     TrfList[CurrTrf].ktr_i=Ktr_i;

//     TrfList[CurrTrf].un=max(Element[i]->Params[2],Element[i]->Params[1]);
     if(paramtype==false)
     {
     TrfList[CurrTrf].un=max(Element[i]->Params[1],Element[i]->Params[2]);
     TrfList[CurrTrf].sn=Element[i]->Params[0];
     TrfList[CurrTrf].pk=Pk;
     TrfList[CurrTrf].uk=Uk;
     TrfList[CurrTrf].px=Px;
     TrfList[CurrTrf].ix=Ix;
     TrfList[CurrTrf].xt=Xt;
     TrfList[CurrTrf].n=Element[i]->Params[7];
     TrfList[CurrTrf].TipN=mswij;
     }
     else
     {
      TrfList[CurrTrf].r=Yr;
      TrfList[CurrTrf].x=Ym;
      TrfList[CurrTrf].g=G0;
      TrfList[CurrTrf].b=B0;
     }
     TrfList[CurrTrf].paramtype=paramtype;
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=mswij;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=mhuj;
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node2)
         flag=false;
     if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node2;
       PquList[CurrPqu].Hui=mhui;
       if(paramtype==false)
         PquList[CurrPqu].Un=Element[i]->Params[1]*Ktr;
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }
     CurrTrf++;
     CurrStr++;

     }
    flag=true;
    for(k=0;k<NumberOfNodes;k++)
      if(PquList[k].NodeNumber==Node1)
        flag=false;
    if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node1;
       PquList[CurrPqu].Hui=mhuj;
       if(paramtype==false)
         PquList[CurrPqu].Un=Element[i]->Params[1];
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }

    }
//�������������� ��-��
for(i=0;i<NumberOfElem;i++)
  if(Element[i]->PointNumber[0]!=-1&&Element[i]->type==4)
    {
     for(k=0;k<NumberOfLines;k++)
       if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[1])
         {
          if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
             mhui=0;
           //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
           //     mhui=0;
          else if(StruList[k].u5<10)
            mhui=StruList[k].u5;
          else
            mhui=1;
          if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
             mhuj=2;
           //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
           //     mhui=0;
          else if(StruList[k].u6<10)
            mhuj=StruList[k].u6;
          else
            mhuj=1;
         }
     Node1=Element[i]->PointNumber[0];
     Node2=Element[i]->PointNumber[1];
     mswij=Element[i]->type;
     paramtype=Element[i]->OtherParams;
     if(paramtype==false)
     {
     Ktr=Element[i]->Params[2]/Element[i]->Params[1];
     }
     else
     {
      TrfList[CurrTrf].r=Element[i]->Params[8];
      TrfList[CurrTrf].x=Element[i]->Params[9];
      Ktr=TrfList[CurrTrf].ktr=real(complex<double>(1.,0)/complex<double>(Element[i]->Params[10],Element[i]->Params[11]));
      TrfList[CurrTrf].ktr_i=imag(complex<double>(1.,0)/complex<double>(Element[i]->Params[10],Element[i]->Params[11]));
      TrfList[CurrTrf].g=Element[i]->Params[12];
      TrfList[CurrTrf].b=Element[i]->Params[13];
     }
     if(Ktr>=1)
       {
        TrfList[CurrTrf].NodeiNumber=Node2;
        TrfList[CurrTrf].NodejNumber=Node1;
//        Ktr=1./Ktr;
        if(paramtype==false)
        {TrfList[CurrTrf].un=Element[i]->Params[2];
        Xt=Element[i]->Params[4]*Element[i]->Params[2]*Element[i]->Params[2]/(100.*Element[i]->Params[0]);
        }
        StrList[CurrStr].NodeiNumber=Node2;
        StrList[CurrStr].NodejNumber=Node1;
        StrList[CurrStr].Mhui=mhuj;
        StrList[CurrStr].Mhuj=mhui;

       }
     else
       {
        TrfList[CurrTrf].NodeiNumber=Node1;
        TrfList[CurrTrf].NodejNumber=Node2;
        if(paramtype==false)
        {TrfList[CurrTrf].un=Element[i]->Params[1];
        Xt=Element[i]->Params[4]*Element[i]->Params[1]*Element[i]->Params[1]/(100.*Element[i]->Params[0]);
        }
        StrList[CurrStr].NodeiNumber=Node1;
        StrList[CurrStr].NodejNumber=Node2;
        StrList[CurrStr].Mswij=mswij;
        StrList[CurrStr].Mhui=mhui;
        StrList[CurrStr].Mhuj=mhuj;

       }
//     Xt=Element[i]->Params[4]*Element[i]->Params[1]*Element[i]->Params[1]/(100.*Element[i]->Params[0]);

     TrfList[CurrTrf].ktr=min(Ktr,1./Ktr);
     if(paramtype==false)
     {TrfList[CurrTrf].sn=Element[i]->Params[0];
     TrfList[CurrTrf].pk=Element[i]->Params[3];
     TrfList[CurrTrf].uk=Element[i]->Params[4];
     TrfList[CurrTrf].px=Element[i]->Params[5];
     TrfList[CurrTrf].ix=Element[i]->Params[6];
     TrfList[CurrTrf].xt=Xt;
     TrfList[CurrTrf].n=Element[i]->Params[7];
     }
     TrfList[CurrTrf].TipN=mswij;

     StrList[CurrStr].Mswij=mswij;
     TrfList[CurrStr].paramtype=paramtype;
     CurrTrf++;
     CurrStr++;

     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node1)
         flag=false;
     if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node1;
       PquList[CurrPqu].Hui=mhui;
       if(paramtype==false)
         PquList[CurrPqu].Un=Element[i]->Params[1];
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node2)
         flag=false;
     if(flag==true)
       {
        PquList[CurrPqu].NodeNumber=Node2;
        PquList[CurrPqu].Hui=mhuj;
        if(paramtype==false)
          PquList[CurrPqu].Un=Element[i]->Params[1]*Ktr;   //������
        PquList[CurrPqu].TipN=0;
        CurrPqu++;
       }
    }
//��������� ��������
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&Element[i]->type==5)
    {
     for(j=0;j<2;j++)
       {
//        switch(j)
 //        {
//          case 0:
          Node1=Element[i]->PointNumber[0];
          Node2=Element[i]->PointNumber[j+1];
//          break;


         for(k=0;k<NumberOfLines;k++)
           if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[j+1])
             {
              if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                 mhui=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u5<10)
                mhui=StruList[k].u5;
              else
                mhui=1;
              if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
                 mhuj=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u6<10)
                mhuj=StruList[k].u6;
              else
                mhuj=1;
             }

   //      }
     mswij=Element[i]->type*10.+j+1;
     X0=Element[i]->Params[2]*(1.-Element[i]->Params[3]);
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=mswij;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=mhuj;
     CurrStr++;
     LinList[CurrLine].NodeiNumber=Node1;
     LinList[CurrLine].NodejNumber=Node2;
     LinList[CurrLine].r0=0;
     LinList[CurrLine].x0=X0;
     LinList[CurrLine].b0=0;
     LinList[CurrLine].Length=1;
     LinList[CurrLine].n=Element[i]->Params[4];
     LinList[CurrLine].TipN=mswij;
     CurrLine++;
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node1)
         flag=false;
     if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node1;
       PquList[CurrPqu].Hui=mhui;
       PquList[CurrPqu].Un=Element[i]->Params[0];
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }

     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node2)
         flag=false;
     if(flag==true)
       {
        PquList[CurrPqu].NodeNumber=Node2;
        PquList[CurrPqu].Hui=mhuj;
        PquList[CurrPqu].Un=Element[i]->Params[0];
        PquList[CurrPqu].TipN=0;
        CurrPqu++;
       }

       }
    }
//��������
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&Element[i]->type==6)
    {
         Node1=Element[i]->PointNumber[0];
         Node2=Element[i]->PointNumber[1];
         mswij=Element[i]->type;
         for(k=0;k<NumberOfLines;k++)
           if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[1])
             {
              if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                 mhui=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u5<10)
                 mhui=StruList[k].u5;
              else
                 mhui=1;
              if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
                 mhuj=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u6<10)
                mhuj=StruList[k].u6;
              else
                mhuj=1;
             }
     X0=Element[i]->Params[2];
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=mswij;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=mhuj;
     CurrStr++;
     LinList[CurrLine].NodeiNumber=Node1;
     LinList[CurrLine].NodejNumber=Node2;
     LinList[CurrLine].r0=0;
     LinList[CurrLine].x0=X0;
     LinList[CurrLine].b0=0;
     LinList[CurrLine].Length=1;
     LinList[CurrLine].n=Element[i]->Params[3];
     LinList[CurrLine].TipN=mswij;
     CurrLine++;
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node1)
         flag=false;
     if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node1;
       PquList[CurrPqu].Hui=mhui;
       PquList[CurrPqu].Un=Element[i]->Params[0];
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node2)
         flag=false;
     if(flag==true)
       {
        PquList[CurrPqu].NodeNumber=Node2;
        PquList[CurrPqu].Hui=mhuj;
        PquList[CurrPqu].Un=Element[i]->Params[0];
        PquList[CurrPqu].TipN=0;
        CurrPqu++;
       }

    }

//����������� ��������
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&Element[i]->type==14)
    {
     ShReactList[CurrShReact].NodeNumber=Element[i]->PointNumber[0];
     ShReactList[CurrShReact].G=Element[i]->Params[1]/1000000.;
     ShReactList[CurrShReact].B=(1.+Rasch->s)*Element[i]->Params[2]/1000000.;
     CurrShReact++;
    }

//�����
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&Element[i]->type==7)
    {
         Node1=Element[i]->PointNumber[0];
         Node2=Element[i]->PointNumber[1];
         mswij=Element[i]->type;
         for(k=0;k<NumberOfLines;k++)
           if(StruList[k].i==Element[i]->PointNumber[0]&&StruList[k].j==Element[i]->PointNumber[1])
             {
              if((StruList[k].u5==3||StruList[k].u5==6)&&kzMode==false)
                 mhui=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u5<10)
                 mhui=StruList[k].u5;
              else
                 mhui=1;
              if((StruList[k].u6==3||StruList[k].u6==6)&&kzMode==false)
                 mhuj=2;
               //  else if((StruList[k].u5>1&&StruList[k].u5!=5)&&kzMode==true)
               //     mhui=0;
              else if(StruList[k].u6<10)
                mhuj=StruList[k].u6;
              else
                mhuj=1;
             }
//     X0=Element[i]->Params[2];
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=mswij;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=mhuj;
     CurrStr++;
     LinList[CurrLine].NodeiNumber=Node1;
     LinList[CurrLine].NodejNumber=Node2;
     LinList[CurrLine].r0=Element[i]->Params[1];
     LinList[CurrLine].x0=Element[i]->Params[2];
     LinList[CurrLine].b0=Element[i]->Params[3];
     LinList[CurrLine].Length=Element[i]->Params[4];
     LinList[CurrLine].n=Element[i]->Params[5];
     LinList[CurrLine].TipN=mswij;
     CurrLine++;
     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node1)
         flag=false;
     if(flag==true)
      {
       PquList[CurrPqu].NodeNumber=Node1;
       PquList[CurrPqu].Hui=mhui;
       PquList[CurrPqu].Un=Element[i]->Params[0];
       PquList[CurrPqu].TipN=0;
       CurrPqu++;
      }

     flag=true;
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node2)
         flag=false;
     if(flag==true)
       {
        PquList[CurrPqu].NodeNumber=Node2;
        PquList[CurrPqu].Hui=mhuj;
        PquList[CurrPqu].Un=Element[i]->Params[0];
        PquList[CurrPqu].TipN=0;
        CurrPqu++;
       }
    }
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==10||Element[i]->type==12||Element[i]->type==13))
    {
//     for(k=0;k<NumberOfNodes;k++)
//       if(PquList[k].NodeNumber==Element[i]->PointNumber[0])
//         {
          Node1=Element[i]->PointNumber[0];
          mswij=Element[i]->type*10+Element[i]->PointNumber[0];
          //TipN=Element[i]->type;
          mhui=2;
          if(Element[i]->type==12||Element[i]->type==13)
            {
             n=Element[i]->Params[5];
             if(n==0)
               n=1;
             Un=Element[i]->Params[1];
             Pn=Element[i]->Params[0]/1000.*n;
             Cos=Element[i]->Params[2];
             KPD=Element[i]->Params[4];
             if(KPD>1||KPD<0.8)
               KPD=1.;
             X0=Element[i]->Params[3];//I�
             if(X0==0)
               X0=0.2;
             Ed=1.1;
             if(Element[i]->type==12)
               {
                X0=1./X0;
                Ed=0.9;
                if(X0>1)
                  X0=0.2;
               }
             Sn=Pn/KPD/Cos;
             Qn=-Sn*sqrt(1-Cos*Cos);
             a0=Element[i]->Params[6];
             if(a0==0)
               a0=1;
             a1=Element[i]->Params[7];
             a2=Element[i]->Params[8];
             b0=Element[i]->Params[9];
             if(b0==0)
               b0=1;
             b1=Element[i]->Params[10];
             b2=Element[i]->Params[11];
            }
          else if(Element[i]->type==10)
            {
             Un=Element[i]->Params[2];
             Pn=Element[i]->Params[0];
             Qn=Element[i]->Params[1];
             X0=0.35;
             Ed=0.85;
             a0=Element[i]->Params[3];
             if(a0==0)
              a0=1;
             a1=Element[i]->Params[4];
             a2=Element[i]->Params[5];
             b0=Element[i]->Params[6];
             if(b0==0)
               b0=1;
             b1=Element[i]->Params[7];
             b2=Element[i]->Params[8];
             Pr=Element[i]->Params[9];
             Qr=Element[i]->Params[10];
            }
//         }
           for(k=0;k<NumberOfLines;k++)
            {
             if(StrList[k].NodeiNumber==Element[i]->PointNumber[0])
                 StrList[k].Mhui=mhui;
             else if(StrList[k].NodejNumber==Element[i]->PointNumber[0])
                 StrList[k].Mhuj=mhui;
            }
           for(k=0;k<NumberOfNodes;k++)
             if(PquList[k].NodeNumber==Element[i]->PointNumber[0])
              if(kzMode==false)
               {

                a0=(a0*Pn+PquList[k].a0*PquList[k].P)/(Pn+PquList[k].P+.001);
                a1=(a1*Pn+PquList[k].a1*PquList[k].P)/(Pn+PquList[k].P+.001);
                a2=(a2*Pn+PquList[k].a2*PquList[k].P)/(Pn+PquList[k].P+.001);
                b0=(b0*Qn+PquList[k].b0*PquList[k].Q)/(Qn+PquList[k].Q+.001);
                b1=(b1*Qn+PquList[k].b1*PquList[k].Q)/(Qn+PquList[k].Q+.001);
                b2=(b2*Qn+PquList[k].b2*PquList[k].Q)/(Qn+PquList[k].Q+.001);

                Pn=PquList[k].P+Pn;
                Qn=PquList[k].Q+Qn;
                PquList[k].Un=Un;
                PquList[k].TipN=Element[i]->type;
                PquList[k].P=Pn;
                PquList[k].Q=Qn;
                PquList[k].a0=a0;
                PquList[k].a1=a1;
                PquList[k].a2=a2;
                PquList[k].b0=b0;
                PquList[k].b1=b1;
                PquList[k].b2=b2;
                PquList[k].dPn+=Pr;
                PquList[k].dQn+=Qr;
               }
           //����� �� �������� ��� ��
    }
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==8))
    {
     Node1=Element[i]->PointNumber[0];
     mhui=0;
     TipN=8;
     Un=Element[i]->Params[1];
     for(j=0;j<CurrStr;j++)
       {
        if(StrList[j].NodeiNumber==Node1)
           StrList[j].Mhui=mhui;
        else if(StrList[j].NodejNumber==Node1)
           StrList[j].Mhuj=mhui;
       }
     for(j=0;j<CurrPqu;j++)
       {
        if(PquList[j].NodeNumber==Node1)
           {
            PquList[j].Un=Un;
            PquList[j].Hui=mhui;
            PquList[j].TipN=TipN;
           }
       }

    }
//���������
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==9))
    {
     Node1=Element[i]->PointNumber[0];
     Xd=Element[i]->Params[4];
//     Xd2=Element[i]->Params[5];
     Un=Element[i]->Params[3];
     for(j=0;j<NumberOfNodes;j++)
       if(PquList[j].NodeNumber==Node1)
          {
           if(PquList[j].Hui==4)
            {
             PquList[j].Hui=1;
             mhui=1;
            }
           else
             mhui=PquList[j].Hui;
          }
/*
     MaxNodeNumber++;
     Node2=MaxNodeNumber;
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=7;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=4;
     if(kzMode==true)
        StrList[CurrStr].Mhuj=3;
     CurrStr++;
     LinList[CurrLine].NodeiNumber=Node1;
     LinList[CurrLine].NodejNumber=Node2;
     LinList[CurrLine].r0=0;
     LinList[CurrLine].x0=Xd/100.;
     LinList[CurrLine].b0=0;
     LinList[CurrLine].Length=1;
     LinList[CurrLine].TipN=7;
     LinList[CurrLine].n=1;
     CurrLine++;

     Node1=Node2;
*/
     mhui=4;
     TipN=9;
     Un=Element[i]->Params[3];
     Pn=Element[i]->Params[0];
     Pr=Element[i]->Params[1];
     Px=Element[i]->Params[10];//dP
     if(Pr==0)
       Pr=Pn;
//     Cos=Element[i]->Params[2];
//     Xd=Element[i]->Params[4];
//     Xg=Xd*(Un*Un)/(Pn/Cos);
     n=Element[i]->Params[7];
//     PZag=Pr/Pn;
/*
     for(ie=1;ie<=1.07;ie+=0.01)
       {
        Sin=Pr*Xg/(Un*ie)/Un;
        if(Sin>1)
          Sin=1;
        Delt=asin(Sin)*360./(2.*M_PI);
        Qr=((Un*ie)*Un*cos(Delt*(2.*M_PI))-Un*Un)/Xg;
        Cos=Pr/sqrt(Pr*Pr+Qr*Qr);
        if(Cos>0.8&&Cos<0.9)
          goto CosCycleEnd;
        Sin=Pr*Xg/Un/Un;
        if(Sin>1)
          Sin=1;
        Delt=asin(Sin)*360./(2.*M_PI);
        Qr=(Un*Un*cos(Delt*(2.*M_PI)/360.)-Un*Un)/Xg;
        Cos=Pr/sqrt(Pr*Pr+Qr*Qr);
       }
     CosCycleEnd:
     if(kzMode==true)
       {
       //����� �������� ��� ������ kz
       }
     for(j=0;j<LinesNumber;j++)
       {
        if(StrList[j].NodeiNumber==Node1)
          StrList[j].Mhui=mhui;
        if(StrList[j].NodejNumber==Node1)
          StrList[j].Mhuj=mhui;
       }
*/
     Qn=Element[i]->Params[8];//Qmin
     Qr=Element[i]->Params[9];//Qmax
     for(k=0;k<NumberOfNodes;k++)
       if(PquList[k].NodeNumber==Node1)
       {
         PquList[k].NodeNumber=Node1;
         PquList[k].Hui=mhui;
         PquList[k].Un=Un;
         PquList[k].TipN=TipN;
//         PquList[k].P+=-Pr*n;
         PquList[k].Pg+=Pr*n;
//         PquList[k].Q=Qr*n;
         PquList[k].a0=1.;
         PquList[k].a1=0;
         PquList[k].a2=0;
         PquList[k].b0=1;
         PquList[k].b1=0;
         PquList[k].b2=0;
         PquList[k].Qmin=Qn;
         PquList[k].Qmax=Qr;
         PquList[k].dPg=Px;
       }

    }

for(i=0;i<NumberOfElem;i++)
   if(Element[i]->PointNumber[0]!=-1&&(Element[i]->type==11))
    {
     Node1=Element[i]->PointNumber[0];
     Xd=Element[i]->Params[3];
//     Xd2=Element[i]->Params[4];
     Un=Element[i]->Params[1];
     for(j=0;j<NumberOfNodes;j++)
       if(PquList[j].NodeNumber==Node1)
          {
           if(PquList[j].Hui==7)
            {
             PquList[j].Hui=1;
             mhui=1;
            }
           else
             mhui=PquList[j].Hui;
          }
     MaxNodeNumber++;
     Node2=MaxNodeNumber;
     StrList[CurrStr].NodeiNumber=Node1;
     StrList[CurrStr].NodejNumber=Node2;
     StrList[CurrStr].Mswij=7;
     StrList[CurrStr].Mhui=mhui;
     StrList[CurrStr].Mhuj=4;
     if(kzMode==true)
       StrList[CurrStr].Mhuj=3;
     CurrStr++;
     LinList[CurrLine].NodeiNumber=Node1;
     LinList[CurrLine].NodejNumber=Node2;
     LinList[CurrLine].r0=0;
     LinList[CurrLine].x0=Xd/100.;
     LinList[CurrLine].b0=0;
     LinList[CurrLine].Length=1;
     LinList[CurrLine].n=1;
     LinList[CurrLine].TipN=7;
     CurrLine++;
     Node2=Node1;
     mhui=4;
     TipN=11;
     Un=Element[i]->Params[1];
     Sn=Element[i]->Params[0];
     Eg=Element[i]->Params[2];
     Xd=Element[i]->Params[3];
     if(Sn==0)
       Xg=0;
     else
       Xg=Xd*(Un*Un)/Sn;
     n=Element[i]->Params[6];
     if(Eg==0)
       Eg=1;
     if(Sn==0)
       {
        Qr=0;
        PZag=0;
       }
     else
       {
        Qr=(Eg*Un*Un-Un*Un)/Xg;
		PZag=abs((long)Qr)/Sn;
       }
     if(PZag>1)
       {
        Qr=-Sn;
        Eg=(Qr*Xg+Un*Un)/(Un*Un);
       }
     for(j=0;j<LinesNumber;j++)
       {
        if(StrList[j].NodeiNumber==Node1)
          StrList[j].Mhui=mhui;
        if(StrList[j].NodejNumber==Node2)
          StrList[j].Mhuj=mhui;
       }
     PquList[CurrPqu].NodeNumber=Node1;
     PquList[CurrPqu].Un=Un;
     PquList[CurrPqu].Hui=mhui;
     PquList[CurrPqu].TipN=TipN;
     PquList[CurrPqu].P=0;
     PquList[CurrPqu].Q=-Qr*n;
     PquList[CurrPqu].a1=0;
     PquList[CurrPqu].a2=Eg;
     PquList[CurrPqu].b0=Xg;
     PquList[CurrPqu].b1=0;
     PquList[CurrPqu].b2=0;
     if(Sn==0)
       PquList[CurrPqu].a0=0;
     else
       PquList[CurrPqu].a0=Qr/Sn;
     CurrPqu++;
     if(kzMode==true)
       {
        //kz �����
       }
    }
//���������� ��� �������� ��� �������
for(i=0;i<CurrPqu;i++)
  {
   if(PquList[i].Hui==5)
     {
      PquList[i].Hui=0;
      PquList[i].Un=0;
     }
   if(PquList[i].Hui==0)
     Ubaz=PquList[i].Un;
   if(PquList[i].Hui==2)
     PquList[i].Hui=20;
  }
for(i=0;i<CurrPqu;i++)
  {
   if(PquList[i].Hui==1)
     PquList[i].Hui=2;
   if(PquList[i].Hui==20)
     PquList[i].Hui=1;
  }
for(i=0;i<CurrPqu;i++)
   if((PquList[i].Hui==2&&PquList[i].TipN==0)||
      (PquList[i].TipN==8&&PquList[i].Hui==0&&kzMode==false)||
      (PquList[i].TipN==0&&PquList[i].Hui==0&&kzMode==true)||
      (PquList[i].TipN==8&&PquList[i].Hui==3&&kzMode==true))
       {
        PquList[i].P=0;
        PquList[i].Q=0;
        PquList[i].a0=0;
        PquList[i].a1=0;
        PquList[i].a2=0;
        PquList[i].b0=0;
        PquList[i].b1=0;
        PquList[i].b2=0;
       }
for(i=0;i<LinesNumber;i++)
  {
   if(StrList[i].Mswij<5)
     StrList[i].Mswij=1;
   else if(StrList[i].Mswij>9&&int(StrList[i].Mswij/10)<5)
     StrList[i].Mswij=1;
   else
     StrList[i].Mswij=2;
  }
OldCurrStr=CurrStr;
for(i=0;i<CurrPqu;i++)
  for(j=0;j<CurrPqu;j++)
    {
     flag=false;
     for(k=0;k<CurrStr;k++)
      {
       if(StrList[k].NodeiNumber==i+1&&StrList[k].NodejNumber==j+1)
          {
           flag=true;
           if(StrList[k].flag==true)
              goto cycfindend;
           if(StrList[k].Mswij==2)
             {
              for(int l=0;l<CurrLine;l++)
               if((LinList[l].NodeiNumber==i+1&&LinList[l].NodejNumber==j+1)||
                 (LinList[l].NodeiNumber==j+1&&LinList[l].NodejNumber==i+1))
                 {
                  Cktr=1;
                  Cktm=0;
                  R0=LinList[l].r0*LinList[l].Length/LinList[l].n;
                  X0=(1.+Rasch->s)*LinList[l].x0*LinList[l].Length/LinList[l].n;
                  B0=(1.+Rasch->s)*LinList[l].b0*LinList[l].Length*0.5*LinList[l].n;
                  if(i<j)
                   {
                    Yr=R0/(R0*R0+X0*X0);
                    Ym=-X0/(R0*R0+X0*X0);
                    Yr2=0;
                    Ym2=B0;
                   }
                  else
                   {
                    Yr=0;
                    Ym=B0;
                    Yr2=R0/(R0*R0+X0*X0);
                    Ym2=-X0/(R0*R0+X0*X0);
                   }
                 }

             }
           else if(StrList[k].Mswij==1)
             {
              for(int l=0;l<CurrTrf;l++)
               if((TrfList[l].NodeiNumber==i+1&&TrfList[l].NodejNumber==j+1)||
                 (TrfList[l].NodeiNumber==j+1&&TrfList[l].NodejNumber==i+1))
                 {
                  Cktr=TrfList[l].ktr;
                  if(TrfList[l].paramtype==false)
                   Cktm=0;
                  else
                   Cktm=TrfList[l].ktr_i;
                  if(Cktr<1)
                     StrList[k].Mswij=-1;


//                  TrfList[l].un*=TrfList[l].ktr;
                  if(TrfList[l].paramtype==false)
                  {
                  X0=(1.+Rasch->s)*(TrfList[l].uk*TrfList[l].un*TrfList[l].un/(100.*TrfList[l].sn))/TrfList[l].n;
                  R0=((TrfList[l].pk*TrfList[l].un*TrfList[l].un/1000.)/(TrfList[l].sn*TrfList[l].sn))/TrfList[l].n;
                  G0=(TrfList[l].px/(TrfList[l].un*TrfList[l].un*1000.))*TrfList[l].n;
                  B0=(1.+Rasch->s)*(TrfList[l].ix*TrfList[l].sn/(100.*TrfList[l].un*TrfList[l].un))*TrfList[l].n;
                  }
                  if(Cktr>=1)
                    {
                     if(TrfList[l].paramtype==false)
                     {Yr=R0/(R0*R0+X0*X0);
                     Ym=-X0/(R0*R0+X0*X0);
                     Yr2=G0;
                     Ym2=-B0;}
                     else
                     {
                      Yr=real(complex<double>(1.,0)/complex<double>(TrfList[l].r,TrfList[l].x));
                      Ym=imag(complex<double>(1.,0)/complex<double>(TrfList[l].r,TrfList[l].x));
                      Yr2=TrfList[l].g;
                      Ym2=-TrfList[l].b;
                     }
                    }
                  else
                    {
                     if(TrfList[l].paramtype==false)
                     {Yr2=R0/(R0*R0+X0*X0);
                     Ym2=-X0/(R0*R0+X0*X0);
                     Yr=G0;
                     Ym=-B0;
                     }
                     else
                     {
                      Yr2=real(complex<double>(1.,0)/complex<double>(TrfList[l].r,TrfList[l].x));
                      Ym2=imag(complex<double>(1.,0)/complex<double>(TrfList[l].r,TrfList[l].x));
                      Yr=TrfList[l].g;
                      Ym=-TrfList[l].b;
                     }

                    }
                 }
             }
           StrList[k].cktr=Cktr;
           StrList[k].cktm=Cktm;
           StrList[k].yr=Yr;
           StrList[k].ym=Ym;
//           StrList[k].flag=true;
           StrList[CurrStr].NodeiNumber=j+1;
           StrList[CurrStr].NodejNumber=i+1;
           if(StrList[k].Mswij==2)
              StrList[CurrStr].Mswij=2;
           else
              StrList[CurrStr].Mswij=-1*StrList[k].Mswij;
//              StrList[CurrStr].Mswij=-1;
           StrList[CurrStr].cktr=real(complex<double>(1.,0)/complex<double>(Cktr,Cktm));
           StrList[CurrStr].cktm=imag(complex<double>(1.,0)/complex<double>(Cktr,Cktm));
//           StrList[CurrStr].cktm=0;
           StrList[CurrStr].yr=Yr2;
           StrList[CurrStr].ym=Ym2;
           StrList[CurrStr].flag=true;
           CurrStr++;
          }
      }
       if(flag==false)
          {
           for(k=0;k<OldCurrStr;k++)
             if(StrList[k].NodeiNumber==j+1&&StrList[k].NodejNumber==i+1)
               flag=true;
           if(flag==false)
           {
            StrList[CurrStr].NodeiNumber=i+1;
            StrList[CurrStr].NodejNumber=j+1;
            StrList[CurrStr].cktr=1;
            StrList[CurrStr].cktm=0;
            StrList[CurrStr].yr=0;
            StrList[CurrStr].ym=0;
//           StrList[CurrStr].flag=true;
            StrList[CurrStr].Mswij=0;
            if(i==j)
              for(int l=0;l<CurrPqu;l++)
                 if(PquList[l].NodeNumber==i+1)
                    StrList[CurrStr].Mswij=PquList[l].Hui;
            CurrStr++;
           }
           flag=false;
          }

     cycfindend:
    }

k=1000000000;
for (i=0;i<CurrStr;i++)
   {
//   int l=StrList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrStr;j++)
      if(StrList[j].NodeiNumber<k)
         {
         k=StrList[j].NodeiNumber;
         tmpStrList=StrList[j];
         StrList[j]=StrList[i];
         StrList[i]=tmpStrList;
         }
   }
k=1000000000;
for (i=0;i<CurrStr;i++)
   {
   int l=StrList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrStr;j++)
      if(StrList[j].NodejNumber<k&&StrList[j].NodeiNumber==l)
         {
         k=StrList[j].NodejNumber;
         tmpStrList=StrList[j];
         StrList[j]=StrList[i];
         StrList[i]=tmpStrList;
         }
   }

k=1000000000;
for (i=0;i<CurrLine;i++)
   {
//   int l=LinList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrLine;j++)
      if(LinList[j].NodeiNumber<k)
         {
         k=LinList[j].NodeiNumber;
         tmpLinList=LinList[j];
         LinList[j]=LinList[i];
         LinList[i]=tmpLinList;
         }
   }
k=1000000000;
for (i=0;i<CurrLine;i++)
   {
   int l=LinList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrLine;j++)
      if(LinList[j].NodejNumber<k&&LinList[j].NodeiNumber==l)
         {
         k=LinList[j].NodejNumber;
         tmpLinList=LinList[j];
         LinList[j]=LinList[i];
         LinList[i]=tmpLinList;
         }
   }

k=1000000000;
for (i=0;i<CurrTrf;i++)
   {
//   int l=TrfList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrTrf;j++)
      if(TrfList[j].NodeiNumber<k)
         {
         k=TrfList[j].NodeiNumber;
         tmpTrfList=TrfList[j];
         TrfList[j]=TrfList[i];
         TrfList[i]=tmpTrfList;
         }
   }
k=1000000000;
for (i=0;i<CurrTrf;i++)
   {
   int l=TrfList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrTrf;j++)
      if(TrfList[j].NodejNumber<k&&TrfList[j].NodeiNumber==l)
         {
         k=TrfList[j].NodejNumber;
         tmpTrfList=TrfList[j];
         TrfList[j]=TrfList[i];
         TrfList[i]=tmpTrfList;
         }
   }

k=1000000000;
for (i=0;i<CurrPqu;i++)
   {
//   int l=TrfList[i].NodeiNumber;
   k=1000000000;
   for(j=i;j<CurrPqu;j++)
      if(PquList[j].NodeNumber<k)
         {
         k=PquList[j].NodeNumber;
         tmpPquList=PquList[j];
         PquList[j]=PquList[i];
         PquList[i]=tmpPquList;
         }
   }
k=1000000000;
for (i=0;i<CurrShReact;i++)
   {
   k=1000000000;
   for(j=i;j<CurrShReact;j++)
      if(ShReactList[j].NodeNumber<k)
         {
         k=ShReactList[j].NodeNumber;
         tmpShReactList=ShReactList[j];
         ShReactList[j]=ShReactList[i];
         ShReactList[i]=tmpShReactList;
         }
   }




TStringList *tmpList=new TStringList;
AnsiString tmpStr,tmpDigit;
;

try
{
tmpList->Clear();
for(i=0;i<CurrStr;i++)
  {
   tmpStr="";

   tmpStr+=AnsiString(StrList[i].NodeiNumber)+", "+AnsiString(StrList[i].NodejNumber)+", "+AnsiString(StrList[i].Mswij)+", ";
   tmpDigit=FloatToStrF(StrList[i].cktr,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(StrList[i].cktm,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(StrList[i].yr,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(StrList[i].ym,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit;
   tmpList->Add(tmpStr);
  }

UstOutDir=GetCurrentDir();
if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regikstr.txt");
tmpList->Clear();
for(i=0;i<CurrPqu;i++)
  {
   tmpStr="";
   tmpStr+=AnsiString(PquList[i].NodeNumber)+", ";

   tmpDigit=FloatToStrF(PquList[i].P,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].Q,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].Un,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit;

   tmpList->Add(tmpStr);
  }

if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regikpqu.txt");
tmpList->Clear();
for(i=0;i<CurrPqu;i++)
  {
   tmpStr="";
   tmpStr+=AnsiString(PquList[i].NodeNumber)+", "+AnsiString(PquList[i].Hui)+", ";

   tmpDigit=FloatToStrF(PquList[i].a0,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].a1,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].a2,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].b0,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].b1,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(PquList[i].b2,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit;

   tmpList->Add(tmpStr);
  }
if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regikuzl.txt");
tmpList->Clear();
for(i=0;i<CurrLine;i++)
  {
   tmpStr="";
   tmpStr+=AnsiString(LinList[i].NodeiNumber)+", "+AnsiString(LinList[i].NodejNumber)+", ";

   tmpDigit=FloatToStrF(LinList[i].r0,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(LinList[i].x0,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(LinList[i].b0,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(LinList[i].Length,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";
   tmpStr+=AnsiString(LinList[i].n);

   tmpList->Add(tmpStr);
  }
if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regiklin.txt");
tmpList->Clear();
for(i=0;i<CurrTrf;i++)
  {
   tmpStr="";
   tmpStr+=AnsiString(TrfList[i].NodeiNumber)+", "+AnsiString(TrfList[i].NodejNumber)+", ";

   tmpDigit=FloatToStrF(TrfList[i].ktr,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].un,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].sn,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].pk,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].uk,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].px,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";

   tmpDigit=FloatToStrF(TrfList[i].ix,ffFixed,7,7);
   for(j=1;j<=tmpDigit.Length();j++)
    if(tmpDigit[j]==',')
       tmpDigit[j]='.';
   tmpStr+=tmpDigit+", ";
   tmpStr+=AnsiString(TrfList[i].n);

   tmpList->Add(tmpStr);
  }
if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regiktrf.txt");
tmpList->Clear();
tmpStr="";
tmpStr+=AnsiString(CurrTrf)+", "+AnsiString(CurrLine)+", "+AnsiString(CurrPqu)+", ";

tmpDigit=FloatToStrF(Ubaz,ffFixed,7,7);
for(j=1;j<=tmpDigit.Length();j++)
 if(tmpDigit[j]==',')
    tmpDigit[j]='.';
tmpStr+=tmpDigit+", 0.001";
tmpList->Add(tmpStr);
//tmpList->SaveToFile("c:\\Temp\\fxpr\\regikkol.txt");
if(IsTxtsOut==true)
tmpList->SaveToFile(UstOutDir+"\\regikkol.txt");
}
__finally
{
delete tmpList;
}
delete Rasch;
Rasch=new Reg(CurrPqu);
Rasch->Initialize(CurrPqu);
Rasch->n=CurrPqu;
Rasch->Obnul();
Rasch->ubaz=Ubaz;
Rasch->kt=CurrTrf;
Rasch->kl=CurrLine;
Rasch->kShReact=CurrShReact;
for(i=0;i<CurrPqu;i++)
  {
   Rasch->sng[i]=complex<double>(PquList[i].P,PquList[i].Q);
   Rasch->Pg[i]=PquList[i].Pg;
   Rasch->Qg[i]=PquList[i].Qg;
   Rasch->dPg[i]=PquList[i].dPg;
   Rasch->dPn[i]=PquList[i].dPn;
   Rasch->dQn[i]=PquList[i].dQn;
   Rasch->unom[i]=PquList[i].Un;
   if(Rasch->unom[i]==complex<double>(0,0))
     for(j=0;j<CurrTrf;j++)
      if((TrfList[j].NodeiNumber==i+1||TrfList[j].NodejNumber==i+1)&&(TrfList[j].ktr==1.000001))
       for(int jj=0;jj<CurrPqu;jj++)
         if(jj!=i&&(PquList[jj].NodeNumber==TrfList[j].NodeiNumber||PquList[jj].NodeNumber==TrfList[j].NodejNumber))
          Rasch->unom[i]=PquList[jj].Un;
   Rasch->a0[i]=PquList[i].a0;
   Rasch->a1[i]=PquList[i].a1;
   Rasch->a2[i]=PquList[i].a2;
   Rasch->b0[i]=PquList[i].b0;
   Rasch->b1[i]=PquList[i].b1;
   Rasch->b2[i]=PquList[i].b2;
   Rasch->mhu[i]=PquList[i].Hui;
   if(PquList[i].Qmax==0&&PquList[i].Qmin==0)
     {
      PquList[i].Qmax=999999999999999;
      PquList[i].Qmin=-999999999999999;
     }
   Rasch->QminQmax[0][i]=PquList[i].Qmin;
   Rasch->QminQmax[1][i]=PquList[i].Qmax;
   for(j=0;j<CurrPqu;j++)
     {
      Rasch->msw[i][j]=StrList[i*CurrPqu+j].Mswij;
      Rasch->cktr[i][j]=complex<double>(StrList[i*CurrPqu+j].cktr,StrList[i*CurrPqu+j].cktm);
      Rasch->y[i][j]=complex<double>(StrList[i*CurrPqu+j].yr,StrList[i*CurrPqu+j].ym);
     }
  }
for(i=0;i<CurrShReact;i++)
  {
   Rasch->y[ShReactList[i].NodeNumber-1][ShReactList[i].NodeNumber-1]=complex<double>(ShReactList[i].G,-ShReactList[i].B);
  }
Rasch->eps=precision;
FillRealNodesNumbers();
delete []TrfList;
delete []PquList;
delete []StrList;
delete []LinList;
delete []ShReactList;

}

//---------------------------------------------------------------------------
int TForm1::PointType(int CurrObjNumber,int PointNumber)
{
int i,j,ptype=-1,oldptype=-1;

for (i=0;i<NumberOfElem;i++)
{
 for(j=0;j<Element[i]->NumberOfPoints;j++)
  if(Element[i]->PointNumber[j]!=-1&&Element[i]->PointNumber[j]==PointNumber&&Element[i]->ObjectNumber!=CurrObjNumber&&Element[i]->ElementName!="����")
   {
    if(Element[i]->ObjectNumber==BalanceElementNo&&Element[i]->BalanceRectNo==j)
      ptype=0;
    else if(Element[i]->ElementName=="��������")
      ptype=2;
    else if(Element[i]->ElementName=="���������")
      ptype=4;
    else if(Element[i]->ElementName=="�������")
      ptype=6;
    else if(Element[i]->ElementName=="���������� �����������")
      ptype=7;
    else if(Element[i]->ElementName=="���������� ���������")
      ptype=8;
    else if(Element[i]->ElementName=="����������� ���������")
      ptype=9;
    else
      ptype=1;
    if((oldptype>ptype&&ptype!=0)||oldptype==0)
      ptype=oldptype;
    oldptype=ptype;
   }
  else if(ptype==-1&&Element[i]->ElementName!="����")
    ptype=1;
}
return ptype;
}
//---------------------------------------------------------------------------

bool TForm1::Raschet(bool showresult)
{
int i,j,k;
int *ki=new int [Rasch->n];
int *kj=new int [Rasch->n];
int ii,jj;
complex<double> tmpCompl,tmpUYs,tmpU,tmpUYb;
complex<double> *tngr=new complex<double>[Rasch->n];

int krm;//���������� �����, �� ������ ��������� � ��
int kit;//���-�� ��������
double p,q;
double osh;//������
krm=0;
AnsiString tmpStr,tmpStr2;
bool Regim=false;
bool RegOk;

//double obuslovl;
//int skip=-1;


//��������� ������� ������������� Ys
krm=Rasch->YsCreate();
//�������� �������� ������� - z[][]
Rasch->ZCreate(ki,kj,ii,jj);

/*
for(i=0;i<Rasch->n;i++)
  if(Rasch->mhu[i]==0)
    skip=i;
obuslovl=Rasch->Norma(Rasch->ys,Rasch->n,skip);
obuslovl*=Rasch->Norma(Rasch->z,Rasch->n-1,-1);
Memo1->Lines->Add("����� ���������������:");
Memo1->Lines->Add(obuslovl);
*/


//������ ������� ����������
for(i=0;i<Rasch->n;i++)
{
 Rasch->urab[i]=Rasch->unom[i];
 if(abs(Rasch->unom[i])!=0)
   {
    if(Rasch->mhu[i]==4)
      Rasch->sng_r[i]=complex<double>(-Rasch->Pg[i],0)+Rasch->sng[i];
    else
      Rasch->sng_r[i]=Rasch->sng[i];
    tngr[i]=Rasch->sng_r[i]/conj(Rasch->unom[i]);
    tngr[i];
   }

}
kit=0;
goto label1;
label2:
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==1)
   {
    p=real(Rasch->sng_r[i])*(Rasch->a0[i]+(Rasch->a1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->a2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
    q=imag(Rasch->sng_r[i])*(Rasch->b0[i]+(Rasch->b1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->b2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
    tngr[i]=complex<double>(p,q)/conj(Rasch->urab[i]);
//    tngr[i];
   }
}
label1:
for(ii=0;ii<krm;ii++)
{
 tmpUYs=0;
 tmpUYb=0;
// Rasch->JacobiCalc();
 for(jj=0;jj<krm;jj++)
   tmpUYs+=Rasch->ys[ki[ii]][kj[jj]]*Rasch->urab[kj[jj]];

/*
 for(j=0;j<Rasch->n;j++)
   if(Rasch->mhu[j]==0||Rasch->mhu[j]==3)
     tmpUYb+=Rasch->y[min(j,ki[ii])][max(j,ki[ii])]*Rasch->urab[j];//����� ����� �������� ����������� �������������??
*/
 for(j=0;j<Rasch->n;j++)
   if(Rasch->mhu[j]==0||Rasch->mhu[j]==3)
     tmpUYb+=Rasch->ys[j][ki[ii]]*Rasch->urab[j];


 Rasch->irab[ki[ii]]=tmpUYs-tngr[ki[ii]]+tmpUYb;
 Rasch->irab[ki[ii]];
}
osh=0;
RegOk=true;
for(ii=0;ii<krm;ii++)
{
 tmpU=complex<double>(0,0);
 for(jj=0;jj<krm;jj++)
   {
    tmpU+=Rasch->z[ii][jj]*Rasch->irab[kj[jj]];
   }
 Rasch->urab[ki[ii]]-=tmpU;
// osh+=abs(tmpU);
 osh=abs(tmpU);
 if(osh>=Rasch->eps)
   RegOk=false;
}
kit++;
if(kit>NumberOfIters)
  goto Raskhod;//����� ����������
osh=osh/(Rasch->n-1);
//if(osh>Rasch->eps)
if(RegOk==false)
  goto label2;//�� ����� �� ��������� �������� - ���������� ������
else
  {
   Regim=true;
   goto Ok;//����� ���������!
  }

Raskhod:
ShowMessage("����� �� �������� - ���������� �����������!");
Regim=false;

Ok:
//����� ������� ��������� ������������ ��������� ������
tmpU=complex<double>(0,0);
for(i=0;i<Rasch->n;i++)
  tmpU+=Rasch->ys[0][i]*Rasch->urab[i];
tmpU=tmpU*conj(Rasch->urab[0]);
Rasch->ResultsCalc();
//� �������, ���� ����
ResultsOut(kit,Regim,showresult);
//delete []Vals;
delete []ki;
delete []kj;
delete []tngr;
return Regim;
}

//---------------------------------------------------------------------------
bool TForm1::RaschetIskl(bool showresult)
{
int i,j,k,l,ubazNo,newubazNo,uzlDelCount;

complex<double> tmpCompl;
complex<double> *tngr=new complex<double>[Rasch->n];
complex<double> *uPrev;
int krm;//���������� �����, �� ������ ��������� � ��
int kit;//���-�� ��������
double p,q;
double osh;//������
krm=0;
AnsiString tmpStr,tmpStr2;
bool Regim=false;
bool RegOk;

Rasch->n++;
int *uzlNo,*uzlDel;
uzlNo=new int [Rasch->n+1];
uzlDel=new int [Rasch->n];

uPrev=new complex<double> [Rasch->n];
Rasch->ysTemp=new complex<double> *[Rasch->n+1];
for(i=0;i<Rasch->n+1;i++)
  Rasch->ysTemp[i]=new complex<double> [Rasch->n+2];


//��������� ������� ������������� Ys
Rasch->n--;//
krm=Rasch->YsCreate();

/*
complex<double> r1,r2,r3;
complex<double> one=complex<double>(1,0);
r1=one/(Rasch->ys[2][0]*Rasch->ys[0][1])/(one/Rasch->ys[0][1]+one/Rasch->ys[1][2]+one/Rasch->ys[2][0]);
r2=one/(Rasch->ys[0][1]*Rasch->ys[1][2])/(one/Rasch->ys[0][1]+one/Rasch->ys[1][2]+one/Rasch->ys[2][0]);
r3=one/(Rasch->ys[1][2]*Rasch->ys[2][0])/(one/Rasch->ys[0][1]+one/Rasch->ys[1][2]+one/Rasch->ys[2][0]);
Rasch->n++;
for(i=0;i<4;i++)
  for(j=0;j<4;j++)
    Rasch->ys[i][j]=complex<double>(0,0);
Rasch->ys[0][0]-=Rasch->ys[0][3]=Rasch->ys[3][0]=one/r1;
Rasch->ys[1][1]-=Rasch->ys[1][3]=Rasch->ys[3][1]=one/r2;
Rasch->ys[2][2]-=Rasch->ys[2][3]=Rasch->ys[3][2]=one/r3;
Rasch->ys[3][3]-=Rasch->ys[0][3]+Rasch->ys[1][3]+Rasch->ys[2][3];
Rasch->sng[3]=complex<double>(0,0);
Rasch->mhu[3]=2;
Rasch->unom[3]=Rasch->urab[3]=Rasch->unom[2];
Rasch->a0[3]=Rasch->b0[3]=1;
Rasch->a1[3]=Rasch->a2[3]=Rasch->b1[3]=Rasch->b2[3]=0;
*/

int m=0;
proba111:
for(i=0;i<Rasch->n+1;i++)
   uzlNo[i]=i;

uzlDelCount=0;

for(i=0;i<Rasch->n;i++)
  {
   if(Rasch->sng[i]==complex<double>(0,0)&&Rasch->mhu[i]!=0)
    {
     uzlDel[uzlDelCount]=i;
     uzlDelCount++;
    }
  }

//uzlDelCount=0;//��� ������ �� ��� ������� ������!!!
for(i=0;i<Rasch->n;i++)
  {
   Rasch->urab[i]=Rasch->unom[i];
   if(Rasch->mhu[i]==0)
     {
      Rasch->ubaz=real(Rasch->unom[i]);
      ubazNo=i;
     }
  }

for(i=0;i<Rasch->n;i++)
  for(j=0;j<Rasch->n;j++)
    Rasch->ysTemp[i][j]=Rasch->ys[i][j];

for(i=0;i<Rasch->n;i++)
  if(Rasch->mhu[i]!=0)
    tngr[i]=-Rasch->ubaz*Rasch->ysTemp[i][ubazNo];

for(l=0;l<uzlDelCount;l++)
{
for (i=0;i<Rasch->n;i++)
//  if(i!=1&&i!=3)
  if(i!=uzlDel[l])
  {
  for (j=0;j<Rasch->n;j++)
//   if(j!=1&&j!=3)
   if(j!=uzlDel[l])
    {
     if(i==j)
       Rasch->ys[i][j]=Rasch->ysTemp[i][j]-=Rasch->ys[i][uzlDel[l]]*Rasch->ys[uzlDel[l]][i]/Rasch->ys[uzlDel[l]][uzlDel[l]];
     else
       Rasch->ys[i][j]=Rasch->ysTemp[i][j]-=Rasch->ys[i][uzlDel[l]]*Rasch->ys[uzlDel[l]][j]/Rasch->ys[uzlDel[l]][uzlDel[l]];
//     Rasch->ysTemp[i][j];
    }
   tngr[i]-=tngr[uzlDel[l]]*Rasch->ys[uzlDel[l]][i]/Rasch->ys[uzlDel[l]][uzlDel[l]];
   tngr[i];
  }
for(i=uzlDel[l]-l;i<Rasch->n;i++)
  uzlNo[i]=uzlNo[i+1];

}
for(i=0;i<Rasch->n-uzlDelCount;i++)
  if(uzlNo[i]==ubazNo)
    newubazNo=i;
for(i=newubazNo;i<Rasch->n-uzlDelCount;i++)
  uzlNo[i]=uzlNo[i+1];

kit=0;
back:
kit++;
Regim=true;
for(i=0;i<Rasch->n;i++)
 {
  for(j=0;j<Rasch->n;j++)
//     Rasch->ysTemp[i][j]=Rasch->ys[i][j];
     Rasch->ys[i][j]=Rasch->ysTemp[i][j];
 //  Rasch->ysTemp[i][Rasch->n-1]=tngr[i]+Rasch->sng[i]/conj(Rasch->urab[i]);

  if(Rasch->mhu[i]!=0)
   {
    p=real(Rasch->sng[i])*(Rasch->a0[i]+(Rasch->a1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->a2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
    q=imag(Rasch->sng[i])*(Rasch->b0[i]+(Rasch->b1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->b2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
   }
  else
   {
    p=0;
    q=0;
   }
//  Rasch->ysTemp[i][Rasch->n]=tngr[i]+Rasch->sng[i]/conj(Rasch->urab[i]);
  Rasch->ys[i][Rasch->n]=tngr[i]+complex<double>(p,q)/conj(Rasch->urab[i]);
  uPrev[i]=Rasch->urab[i];
 }

Rasch->Gauss(Rasch->ys,Rasch->urab,Rasch->n-uzlDelCount-1,uzlNo);
for(i=0;i<Rasch->n-uzlDelCount;i++)
  {
   osh=abs(uPrev[uzlNo[i]])-abs(Rasch->urab[uzlNo[i]]);
   if(fabs(osh)>Rasch->eps)
     Regim=false;
  }
if(Regim==false&&kit>NumberOfIters)
 goto RegimRaskhod;
else if(Regim==false)
 goto back;
else
 goto RegimOk;

/*
//����� ������� ��������� ������������ ��������� ������
Rasch->ResultsCalc();
//� �������, ���� ����
*/

RegimRaskhod:
ShowMessage("����� �� �������� - ���������� �����������!");
Regim=false;

RegimOk:
bool Yes=true;

//����� tngr[]- ������ ���������� ������������ �� ����� ��� �����
for(i=0;i<Rasch->n;i++)
  {
  for(j=0;j<uzlDelCount;j++)
     if(uzlDel[j]==i)
       Yes=false;
  if(Yes==true)
    {
     tngr[i]=-Rasch->ysTemp[i][i];
     for(k=0;k<Rasch->n;k++)
       {
        for(j=0;j<uzlDelCount;j++)
          if(uzlDel[j]==i)
            Yes=false;
        if(Yes==true&&k!=i)
          tngr[i]-=Rasch->ysTemp[i][k];
        Yes=true;
       }
    }
  Yes=true;
  tngr[i]=complex<double>(0,0);//��� ������ �� ��� ������� ������!!!!
  }


//for(i=0;i<Rasch->n-uzlDelCount;i++)
for(i=0;i<Rasch->n;i++)
  {
  for(j=0;j<uzlDelCount;j++)
     if(uzlDel[j]==i)
       Yes=false;
  if(Rasch->mhu[i]!=0&&Yes==true)
     {
      p=real(Rasch->sng[i])*(Rasch->a0[i]+(Rasch->a1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->a2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
      q=imag(Rasch->sng[i])*(Rasch->b0[i]+(Rasch->b1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->b2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
      Rasch->ys[i][newubazNo]=Rasch->ys[newubazNo][i]=(Rasch->urab[i]*tngr[i]+complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[ubazNo]-Rasch->urab[i]);
//      Rasch->ysTemp[i][i]=-Rasch->ysTemp[i][newubazNo];
      tmpCompl=complex<double>(1,0)/Rasch->ys[i][newubazNo];
      tmpCompl=complex<double>(1,0)/(-Rasch->ys[i][newubazNo]);
      Rasch->ys[i][i]=-Rasch->ys[i][newubazNo]-tngr[i];
     }
  Yes=true;
  }
//Rasch->ysTemp[newubazNo][newubazNo]=0;
Rasch->ys[newubazNo][newubazNo]=-tngr[newubazNo];

//Yes=true;
//for(i=0;j<Rasch->n-uzlDelCount;i++)
for(i=0;i<Rasch->n;i++)
  {
   for(j=0;j<uzlDelCount;j++)
     if(uzlDel[j]==i)
       Yes=false;
  if(Rasch->mhu[i]!=0&&Yes==true)
    Rasch->ys[newubazNo][newubazNo]-=Rasch->ys[newubazNo][i];
  Yes=true;
  }
m++;

/*
for(i=0;i<Rasch->n;i++)
  for(j=0;j<Rasch->n;j++)
    Rasch->ys[i][j]=Rasch->ysTemp[i][j];
if(m==1)
  goto proba111;
*/
Rasch->ysTemp[newubazNo][newubazNo];
uPrev[0]=complex<double>(0,0);
uPrev[1]=complex<double>(0,0);

for(i=0;i<Rasch->n;i++)
  {
   uPrev[0]+=Rasch->urab[i]*Rasch->ysTemp[0][i];
  }
Yes=true;
krm=Rasch->YsCreate();
for(i=0;i<Rasch->n;i++)
  {
   for(j=0;j<uzlDelCount;j++)
     if(uzlDel[j]==i)
       Yes=false;
   if(Yes==true)
//   uPrev[1]+=Rasch->urab[uzlNo[i]]*Rasch->ys[0][uzlNo[i]];
   uPrev[1]+=Rasch->urab[i]*Rasch->ys[0][i];
   Yes=true;
  }
uPrev[0];
uPrev[1];





//Vals=new Values[1];
if(showresult==true)
  ResultsOut(kit,Regim,true);
//delete []Vals;
for(i=0;i<Rasch->n+1;i++)
  {
   delete Rasch->ysTemp[i];
  }
delete uPrev;
delete []uzlNo;
delete []uzlDel;
delete []tngr;
return Regim;
}
//---------------------------------------------------------------------------

bool TForm1::RaschetJ(bool showresult,bool &continueJ)
{
int i,state,kit=0;
bool Regim,limits=false;

Rasch->YsCreate();
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
double *dX=new double [Rasch->n*2];
//Rasch->mhu_n=new int [Rasch->n];
if(Rasch->CurrReg==false)
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
//Rasch->SHN(false);

/*
Memo1->Lines->Add("Z0");
for(i=0;i<Rasch->n;i++)
 if(Rasch->mhu[i]!=0&&Rasch->mhu[i]!=4)
  Memo1->Lines->Add(AnsiString(abs(Rasch->urab[i])));
for(i=0;i<Rasch->n;i++)
 if(Rasch->mhu[i]!=0)
  Memo1->Lines->Add(AnsiString(arg(Rasch->urab[i])));
Memo1->Lines->Add("End-Z");
*/
GesseParams *GP=new GesseParams(0,Rasch->nJac,Rasch->n,Rasch->mhu);
double Bk;
bool BkUse=false;
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
//   Memo1->Lines->Add(AnsiString(i)+" obusl: "+AnsiString (int(1./Rasch->ObuslMatr(Rasch->Jac,Rasch->nJac))));
   Rasch->LU(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac);
//   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   if(Rasch->BkUse==true)
   {
    Bk=Rasch->BkCalc(GP,dX,Rasch->Nebal);
    if(Bk>=1)
     BkUse=true;
    else
     BkUse=false;
   }
   state=!Rasch->RaschCheck_new(BkUse,Bk,dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
delete GP;
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacRasch;
  }

for(i=0;i<Rasch->n;i++)
 Rasch->mhu_n[i]=Rasch->mhu[i];

if(UseQminQmax==true)
if(Rasch->Limits(0,false,0)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
if(UseQminQmax==true)
if(Rasch->Limits(1,false,0)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}

if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacRasch;
  }
for(i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];

Rasch->ResultsCalc();
ResultsOut(kit,Regim,showresult);

continueJ=true;
WrongJacRasch:
delete []dX;
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------
bool TForm1::RaschetJKn(double Kn,bool UseUnom)
{
int i,state,kit=0;
bool Regim,limits=false;
//AnsiString InputKn;
bool sobstvOk;

/*
if(askKn==true)
{
Rasch->Kn=1;
InputKn=AnsiString(Kn);
InputKn = InputBox("������� ��", "����������� �������������", InputKn);
Rasch->Kn=StrToFloat(InputKn);
}
*/
//Rasch->YsCreate();

Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
double *dX=new double [Rasch->n*2];

if(UseUnom==true)
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];

for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Kn=Kn;
   Rasch->Jacobi_KnCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }

if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;

// for(i=0;i<Rasch->n;i++)
//  Rasch->mhu_n[i]=Rasch->mhu[i];

//if(ogran==true)
if(UseQminQmax==true)
if(Rasch->Limits(0,false,0)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi_KnCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
//if(ogran==true)
if(UseQminQmax==true)
if(Rasch->Limits(1,false,0)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi_KnCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}

//for(i=0;i<Rasch->n;i++)
//  Rasch->mhu[i]=Rasch->mhu_n[i];

Rasch->ResultsCalc();
ResultsOut(kit,Regim,false);



delete []dX;
Rasch->RegimState=Regim;
return Regim;
}

//---------------------------------------------------------------------------
bool TForm1::RaschetJProba(bool showresult,bool &continueJ)
{
int i,j,state;
bool Regim;
Rasch->n=63;
for(i=0;i<Rasch->n;i++)
  for(j=0;j<Rasch->n;j++)
    Rasch->y[i][j]=complex<double>(0,0);

Rasch->mhu[0]=0;Rasch->Nu[0]=1;
Rasch->mhu[1]=4;Rasch->Nu[1]=2;
Rasch->mhu[2]=4;Rasch->Nu[2]=3;
Rasch->mhu[3]=4;Rasch->Nu[3]=4;
Rasch->mhu[4]=4;Rasch->Nu[4]=5;
Rasch->mhu[5]=4;Rasch->Nu[5]=6;
Rasch->mhu[6]=4;Rasch->Nu[6]=7;
Rasch->mhu[7]=4;Rasch->Nu[7]=8;
Rasch->mhu[8]=4;Rasch->Nu[8]=9;
Rasch->mhu[9]=4;Rasch->Nu[9]=10;
Rasch->mhu[10]=4;Rasch->Nu[10]=11;
Rasch->mhu[11]=4;Rasch->Nu[11]=12;
Rasch->mhu[12]=1;Rasch->Nu[12]=101;
Rasch->mhu[13]=1;Rasch->Nu[13]=102;
Rasch->mhu[14]=1;Rasch->Nu[14]=103;
Rasch->mhu[15]=1;Rasch->Nu[15]=500;
Rasch->mhu[16]=1;Rasch->Nu[16]=703;
Rasch->mhu[17]=1;Rasch->Nu[17]=704;
Rasch->mhu[18]=1;Rasch->Nu[18]=705;
Rasch->mhu[19]=4;Rasch->Nu[19]=706;
Rasch->mhu[20]=1;Rasch->Nu[20]=707;
Rasch->mhu[21]=1;Rasch->Nu[21]=708;
Rasch->mhu[22]=1;Rasch->Nu[22]=710;
Rasch->mhu[23]=1;Rasch->Nu[23]=712;
Rasch->mhu[24]=1;Rasch->Nu[24]=714;
Rasch->mhu[25]=1;Rasch->Nu[25]=715;
Rasch->mhu[26]=1;Rasch->Nu[26]=716;
Rasch->mhu[27]=1;Rasch->Nu[27]=717;
Rasch->mhu[28]=1;Rasch->Nu[28]=718;
Rasch->mhu[29]=1;Rasch->Nu[29]=719;
Rasch->mhu[30]=1;Rasch->Nu[30]=722;
Rasch->mhu[31]=1;Rasch->Nu[31]=731;
Rasch->mhu[32]=1;Rasch->Nu[32]=760;
Rasch->mhu[33]=1;Rasch->Nu[33]=810;
Rasch->mhu[34]=1;Rasch->Nu[34]=811;
Rasch->mhu[35]=1;Rasch->Nu[35]=812;
Rasch->mhu[36]=1;Rasch->Nu[36]=814;
Rasch->mhu[37]=1;Rasch->Nu[37]=818;
Rasch->mhu[38]=1;Rasch->Nu[38]=819;
Rasch->mhu[39]=1;Rasch->Nu[39]=821;
Rasch->mhu[40]=1;Rasch->Nu[40]=822;
Rasch->mhu[41]=1;Rasch->Nu[41]=824;
Rasch->mhu[42]=1;Rasch->Nu[42]=826;
Rasch->mhu[43]=1;Rasch->Nu[43]=901;
Rasch->mhu[44]=1;Rasch->Nu[44]=902;
Rasch->mhu[45]=1;Rasch->Nu[45]=903;
Rasch->mhu[46]=1;Rasch->Nu[46]=904;
Rasch->mhu[47]=1;Rasch->Nu[47]=905;
Rasch->mhu[48]=1;Rasch->Nu[48]=906;
Rasch->mhu[49]=1;Rasch->Nu[49]=907;
Rasch->mhu[50]=1;Rasch->Nu[50]=908;
Rasch->mhu[51]=1;Rasch->Nu[51]=912;
Rasch->mhu[52]=1;Rasch->Nu[52]=913;
Rasch->mhu[53]=1;Rasch->Nu[53]=914;
Rasch->mhu[54]=1;Rasch->Nu[54]=917;
Rasch->mhu[55]=1;Rasch->Nu[55]=918;
Rasch->mhu[56]=1;Rasch->Nu[56]=920;
Rasch->mhu[57]=1;Rasch->Nu[57]=921;
Rasch->mhu[58]=1;Rasch->Nu[58]=923;
Rasch->mhu[59]=1;Rasch->Nu[59]=924;
Rasch->mhu[60]=1;Rasch->Nu[60]=929;
Rasch->mhu[61]=1;Rasch->Nu[61]=937;
Rasch->mhu[62]=1;Rasch->Nu[62]=944;
Rasch->unom[0]=complex<double>(530,0);
Rasch->unom[1]=complex<double>(540,0);
Rasch->unom[2]=complex<double>(525,0);
Rasch->unom[3]=complex<double>(500,0);
Rasch->unom[4]=complex<double>(525,0);
Rasch->unom[5]=complex<double>(525,0);
Rasch->unom[6]=complex<double>(525,0);
Rasch->unom[7]=complex<double>(525,0);
Rasch->unom[8]=complex<double>(542,0);
Rasch->unom[9]=complex<double>(542,0);
Rasch->unom[10]=complex<double>(510,0);
Rasch->unom[11]=complex<double>(525,0);
Rasch->unom[12]=complex<double>(537,0);
Rasch->unom[13]=complex<double>(536,0);
Rasch->unom[14]=complex<double>(538,0);
Rasch->unom[15]=complex<double>(524.01,0);
Rasch->unom[16]=complex<double>(520.34,0);
Rasch->unom[17]=complex<double>(518.43,0);
Rasch->unom[18]=complex<double>(535.43,0);
Rasch->unom[19]=complex<double>(242,0);
Rasch->unom[20]=complex<double>(1143.4,0);
Rasch->unom[21]=complex<double>(231.96,0);
Rasch->unom[22]=complex<double>(520.37,0);
Rasch->unom[23]=complex<double>(517.51,0);
Rasch->unom[24]=complex<double>(505.13,0);
Rasch->unom[25]=complex<double>(500.52,0);
Rasch->unom[26]=complex<double>(501.06,0);
Rasch->unom[27]=complex<double>(512.41,0);
Rasch->unom[28]=complex<double>(1124.8,0);
Rasch->unom[29]=complex<double>(510.63,0);
Rasch->unom[30]=complex<double>(519.77,0);
Rasch->unom[31]=complex<double>(522.13,0);
Rasch->unom[32]=complex<double>(1181.4,0);
Rasch->unom[33]=complex<double>(524.12,0);
Rasch->unom[34]=complex<double>(512.62,0);
Rasch->unom[35]=complex<double>(540.88,0);
Rasch->unom[36]=complex<double>(512.09,0);
Rasch->unom[37]=complex<double>(511.78,0);
Rasch->unom[38]=complex<double>(508.39,0);
Rasch->unom[39]=complex<double>(504,0);
Rasch->unom[40]=complex<double>(519.37,0);
Rasch->unom[41]=complex<double>(518.83,0);
Rasch->unom[42]=complex<double>(1170.2,0);
Rasch->unom[43]=complex<double>(529.77,0);
Rasch->unom[44]=complex<double>(1149.4,0);
Rasch->unom[45]=complex<double>(1138.6,0);
Rasch->unom[46]=complex<double>(531.2,0);
Rasch->unom[47]=complex<double>(527.5,0);
Rasch->unom[48]=complex<double>(536.90,0);
Rasch->unom[49]=complex<double>(533.3,0);
Rasch->unom[50]=complex<double>(519.91,0);
Rasch->unom[51]=complex<double>(519.55,0);
Rasch->unom[52]=complex<double>(520.49,0);
Rasch->unom[53]=complex<double>(1191.2,0);
Rasch->unom[54]=complex<double>(538.46,0);
Rasch->unom[55]=complex<double>(536.88,0);
Rasch->unom[56]=complex<double>(528.31,0);
Rasch->unom[57]=complex<double>(511.7,0);
Rasch->unom[58]=complex<double>(520.46,0);
Rasch->unom[59]=complex<double>(516.15,0);
Rasch->unom[60]=complex<double>(236.03,0);
Rasch->unom[61]=complex<double>(236.61,0);
Rasch->unom[62]=complex<double>(521.54,0);
Rasch->sng[0]=complex<double>(0,0);
Rasch->sng[1]=complex<double>(-3000,0);
Rasch->sng[2]=complex<double>(-1955,0);
Rasch->sng[3]=complex<double>(-2400,0);
Rasch->sng[4]=complex<double>(-500,0);
Rasch->sng[5]=complex<double>(-2400,0);
Rasch->sng[6]=complex<double>(-500,0);
Rasch->sng[7]=complex<double>(-1000,0);
Rasch->sng[8]=complex<double>(-1200,0);
Rasch->sng[9]=complex<double>(-500,0);
Rasch->sng[10]=complex<double>(-4480,0);
Rasch->sng[11]=complex<double>(1170,-925);
Rasch->sng[12]=complex<double>(0,0);
Rasch->sng[13]=complex<double>(0,0);
Rasch->sng[14]=complex<double>(0,0);
Rasch->sng[15]=complex<double>(0,0);
Rasch->sng[16]=complex<double>(241,-114.5);
Rasch->sng[17]=complex<double>(541,-232.8);
Rasch->sng[18]=complex<double>(0,0);
Rasch->sng[19]=complex<double>(-429,-660);
Rasch->sng[20]=complex<double>(0,0);Rasch->y[20][20]=complex<double>(0,-0.001875);
Rasch->sng[21]=complex<double>(477,-213);
Rasch->sng[22]=complex<double>(32,0);
Rasch->sng[23]=complex<double>(295,-196.70);
Rasch->sng[24]=complex<double>(0,0);
Rasch->sng[25]=complex<double>(1005,-400.4);
Rasch->sng[26]=complex<double>(535,-150.30);
Rasch->sng[27]=complex<double>(120,-61.50);
Rasch->sng[28]=complex<double>(0,0);Rasch->y[28][28]=complex<double>(0,-0.000625);
Rasch->sng[29]=complex<double>(642,-249.2);
Rasch->sng[30]=complex<double>(0,0);
Rasch->sng[31]=complex<double>(272,-197.40);
Rasch->sng[32]=complex<double>(0,0);Rasch->y[32][32]=complex<double>(0,-0.001250);
Rasch->sng[33]=complex<double>(2460,-733.80);
Rasch->sng[34]=complex<double>(1000,-242.70);
Rasch->sng[35]=complex<double>(487,-26.30);
Rasch->sng[36]=complex<double>(658,-273.80);
Rasch->sng[37]=complex<double>(850,0);
Rasch->sng[38]=complex<double>(1059,-405.90);
Rasch->sng[39]=complex<double>(850,-516.90);
Rasch->sng[40]=complex<double>(567.00,0);
Rasch->sng[41]=complex<double>(1062,-427.30);
Rasch->sng[42]=complex<double>(0,0);Rasch->y[42][42]=complex<double>(0,-0.000625);
Rasch->sng[43]=complex<double>(707,-503.30);Rasch->y[43][43]=complex<double>(0,-0.000653);
Rasch->sng[44]=complex<double>(0,0);Rasch->y[44][44]=complex<double>(0,-0.001250);
Rasch->sng[45]=complex<double>(0,0);Rasch->y[45][45]=complex<double>(0,-0.001875);
Rasch->sng[46]=complex<double>(382,-263.2);
Rasch->sng[47]=complex<double>(402,-133);
Rasch->sng[48]=complex<double>(321.0,-172.0);
Rasch->sng[49]=complex<double>(240,-142.1);
Rasch->sng[50]=complex<double>(776,-183.0);Rasch->y[50][50]=complex<double>(0,-0.001959);
Rasch->sng[51]=complex<double>(0,0);
Rasch->sng[52]=complex<double>(637,-187.0);
Rasch->sng[53]=complex<double>(0,0);Rasch->y[53][53]=complex<double>(0,-0.001250);
Rasch->sng[54]=complex<double>(432,-170.1);
Rasch->sng[55]=complex<double>(217,-105.2);
Rasch->sng[56]=complex<double>(166.0,-61.5);
Rasch->sng[57]=complex<double>(0,-497.5);Rasch->y[57][57]=complex<double>(0,-0.001959);
Rasch->sng[58]=complex<double>(398,-270.60);
Rasch->sng[59]=complex<double>(200,-117.70);
Rasch->sng[60]=complex<double>(600,-275.70);
Rasch->sng[61]=complex<double>(320,-215.10);
Rasch->sng[62]=complex<double>(0,0);
for(i=0;i<Rasch->n;i++)
 {
  Rasch->a0[i]=Rasch->b0[i]=1;
  Rasch->a1[i]=Rasch->a2[i]=Rasch->b1[i]=Rasch->b2[i]=0;
  for(j=0;j<Rasch->n;j++)
   if(i!=j)
   {
    Rasch->y[i][j]=complex<double>(0,0);
    Rasch->msw[i][j]=0;
    Rasch->cktr[i][j]=complex<double>(1,0);
   }
 }
Rasch->msw[0][15]=Rasch->msw[15][0]=2;Rasch->y[0][15]=1./complex<double>(0,10);
Rasch->msw[1][33]=Rasch->msw[33][1]=2;Rasch->y[1][33]=1./complex<double>(0,7.3);
Rasch->msw[2][41]=Rasch->msw[41][2]=2;Rasch->y[2][41]=1./complex<double>(0,12.2);
Rasch->msw[3][27]=Rasch->msw[27][3]=2;Rasch->y[3][27]=1./complex<double>(0,10.4);
Rasch->msw[4][52]=Rasch->msw[52][4]=2;Rasch->y[4][52]=1./complex<double>(0,61.1);
Rasch->msw[5][50]=Rasch->msw[50][5]=2;Rasch->y[5][50]=1./complex<double>(0,10.2);
Rasch->msw[6][51]=Rasch->msw[51][6]=2;Rasch->y[6][51]=1./complex<double>(0,61.1);
Rasch->msw[7][57]=Rasch->msw[57][7]=2;Rasch->y[7][57]=1./complex<double>(0,17.7);
Rasch->msw[8][60]=Rasch->msw[60][8]=2;Rasch->y[8][60]=1./complex<double>(0,4.6);
Rasch->msw[9][61]=Rasch->msw[61][9]=2;Rasch->y[9][61]=1./complex<double>(0,5.8);
Rasch->msw[10][24]=Rasch->msw[24][10]=2;Rasch->y[10][24]=1./complex<double>(0,7.7);
Rasch->msw[11][30]=Rasch->msw[30][11]=2;Rasch->y[11][30]=1./complex<double>(0,17);
Rasch->msw[12][47]=Rasch->msw[47][12]=2;Rasch->y[12][47]=1./complex<double>(0,79.9);
Rasch->msw[13][48]=Rasch->msw[48][13]=2;Rasch->y[13][48]=1./complex<double>(0,79.9);
Rasch->msw[14][54]=Rasch->msw[54][14]=2;Rasch->y[14][54]=1./complex<double>(0,79.9);
Rasch->msw[15][33]=Rasch->msw[33][15]=2;Rasch->y[15][33]=1./complex<double>(0,70);
Rasch->msw[15][36]=Rasch->msw[36][15]=2;Rasch->y[15][36]=1./complex<double>(5.2,69.5);Rasch->y[36][15]=complex<double>(0,.000893);
Rasch->msw[16][17]=Rasch->msw[17][16]=2;Rasch->y[16][17]=1./complex<double>(3.1,39);Rasch->y[17][16]=complex<double>(0,.000460);
Rasch->msw[16][31]=Rasch->msw[31][16]=2;Rasch->y[16][31]=1./complex<double>(3.1,39);Rasch->y[31][16]=complex<double>(0,.000460);
Rasch->msw[17][22]=Rasch->msw[22][17]=2;Rasch->y[17][22]=1./complex<double>(5.4,68.7);Rasch->y[22][17]=complex<double>(0,.000797);
Rasch->msw[17][27]=Rasch->msw[27][17]=2;Rasch->y[17][27]=1./complex<double>(5.5,70);Rasch->y[27][17]=complex<double>(0,.000825);
Rasch->msw[17][29]=Rasch->msw[29][17]=2;Rasch->y[17][29]=1./complex<double>(2.3,25);Rasch->y[29][17]=complex<double>(0,.000290);
Rasch->msw[17][30]=Rasch->msw[30][17]=2;Rasch->y[17][30]=1./complex<double>(6.4,85.1);Rasch->y[30][17]=complex<double>(0,.001015);

//Rasch->msw[20][18]=-1;Rasch->msw[18][20]=1;Rasch->cktr[18][20]=complex<double>(2.1299,0);Rasch->cktr[20][12]=1./Rasch->cktr[18][20];Rasch->y[18][20]=1./complex<double>(0,40.6)/Rasch->cktr[18][20]/Rasch->cktr[18][20];
//Rasch->msw[18][21]=-1;Rasch->msw[21][18]=1;Rasch->cktr[21][18]=complex<double>(2.1299,0);Rasch->cktr[18][21]=1./Rasch->cktr[21][18];Rasch->y[21][18]=1./complex<double>(0,27.0)/Rasch->cktr[21][18]/Rasch->cktr[21][18];
Rasch->msw[20][18]=-1;Rasch->msw[18][20]=1;Rasch->cktr[18][20]=complex<double>(2.1299,0);Rasch->cktr[20][18]=1./Rasch->cktr[18][20];Rasch->y[18][20]=1./complex<double>(0,40.6);
Rasch->msw[18][21]=-1;Rasch->msw[21][18]=1;Rasch->cktr[21][18]=complex<double>(2.1299,0);Rasch->cktr[18][21]=1./Rasch->cktr[21][18];Rasch->y[21][18]=1./complex<double>(0,27.0);

Rasch->msw[18][23]=Rasch->msw[23][18]=2;Rasch->y[18][23]=1./complex<double>(8.2,84.7);Rasch->y[23][18]=complex<double>(0,.000991);
Rasch->msw[18][31]=Rasch->msw[31][18]=2;Rasch->y[18][31]=1./complex<double>(5.2,55);Rasch->y[31][18]=complex<double>(0,.000650);
Rasch->msw[18][58]=Rasch->msw[58][18]=2;Rasch->y[18][58]=1./complex<double>(10.3,109);Rasch->y[58][18]=complex<double>(0,.001265);
Rasch->msw[19][21]=Rasch->msw[21][19]=2;Rasch->y[19][21]=1./complex<double>(0,20.6);
Rasch->msw[20][28]=Rasch->msw[28][20]=2;Rasch->y[20][28]=1./complex<double>(4.02,114);Rasch->y[28][20]=complex<double>(0,.001988);
Rasch->msw[20][32]=Rasch->msw[32][20]=2;Rasch->y[20][32]=1./complex<double>(3.75,92.2);Rasch->y[32][20]=complex<double>(0,.001644);
Rasch->msw[22][23]=Rasch->msw[23][22]=2;Rasch->y[22][23]=1./complex<double>(0,19.5);
Rasch->msw[22][31]=Rasch->msw[31][22]=2;Rasch->y[22][31]=1./complex<double>(0,75);
Rasch->msw[23][24]=Rasch->msw[24][23]=2;Rasch->y[23][24]=1./complex<double>(6.6,70.1);Rasch->y[24][23]=complex<double>(0,.003260);
Rasch->msw[24][25]=Rasch->msw[25][24]=2;Rasch->y[24][25]=1./complex<double>(.3,4.8);Rasch->y[25][24]=complex<double>(0,.000230);
Rasch->msw[25][26]=Rasch->msw[26][25]=2;Rasch->y[25][26]=1./complex<double>(1,12.5);Rasch->y[26][25]=complex<double>(0,.000620);
Rasch->msw[26][27]=Rasch->msw[27][26]=2;Rasch->y[26][27]=1./complex<double>(4,42);Rasch->y[27][26]=complex<double>(0,.001950);

//Rasch->msw[28][27]=-1;Rasch->msw[27][28]=1;Rasch->cktr[27][28]=complex<double>(2.1701,0);Rasch->cktr[28][27]=1./Rasch->cktr[27][28];Rasch->y[27][28]=1./complex<double>(0,27.1)/Rasch->cktr[27][28]/Rasch->cktr[27][28];
Rasch->msw[28][27]=-1;Rasch->msw[27][28]=1;Rasch->cktr[27][28]=complex<double>(2.1701,0);Rasch->cktr[28][27]=1./Rasch->cktr[27][28];Rasch->y[27][28]=1./complex<double>(0,27.1);

Rasch->msw[27][29]=Rasch->msw[29][27]=2;Rasch->y[27][29]=1./complex<double>(10.2,108);Rasch->y[29][27]=complex<double>(0,.001260);
Rasch->msw[27][30]=Rasch->msw[30][27]=2;Rasch->y[27][30]=1./complex<double>(2.9,35.5);Rasch->y[30][27]=complex<double>(0,.000430);
Rasch->msw[32][53]=Rasch->msw[53][32]=2;Rasch->y[32][53]=1./complex<double>(3.6,88.6);Rasch->y[53][32]=complex<double>(0,.001572);
Rasch->msw[33][34]=Rasch->msw[34][33]=2;Rasch->y[33][34]=1./complex<double>(6.9,66.4);Rasch->y[34][33]=complex<double>(0,.000842);
Rasch->msw[33][38]=Rasch->msw[38][33]=2;Rasch->y[33][38]=1./complex<double>(6,75);Rasch->y[38][33]=complex<double>(0,.000072);
Rasch->msw[34][35]=Rasch->msw[35][34]=2;Rasch->y[34][35]=1./complex<double>(7.6,89);Rasch->y[35][34]=complex<double>(0,.000911);
Rasch->msw[34][38]=Rasch->msw[38][34]=2;Rasch->y[34][38]=1./complex<double>(2,19.2);Rasch->y[38][34]=complex<double>(0,.000250);
Rasch->msw[35][54]=Rasch->msw[54][35]=2;Rasch->y[35][54]=1./complex<double>(9,93);Rasch->y[54][35]=complex<double>(0,.001160);
Rasch->msw[36][37]=Rasch->msw[37][36]=2;Rasch->y[36][37]=1./complex<double>(7,98.4);Rasch->y[37][36]=complex<double>(0,.001862);
Rasch->msw[36][39]=Rasch->msw[39][36]=2;Rasch->y[36][39]=1./complex<double>(8.8,97.6);Rasch->y[39][36]=complex<double>(0,.001141);
Rasch->msw[37][38]=Rasch->msw[38][37]=2;Rasch->y[37][38]=1./complex<double>(1.2,15);Rasch->y[38][37]=complex<double>(0,.000220);

//Rasch->msw[42][37]=-1;Rasch->msw[37][42]=1;Rasch->cktr[37][42]=complex<double>(2.2999,0);Rasch->cktr[42][37]=1./Rasch->cktr[37][42];Rasch->y[37][42]=1./complex<double>(0,40.6)/Rasch->cktr[37][42]/Rasch->cktr[37][42];
Rasch->msw[42][37]=-1;Rasch->msw[37][42]=1;Rasch->cktr[37][42]=complex<double>(2.2999,0);Rasch->cktr[42][37]=1./Rasch->cktr[37][42];Rasch->y[37][42]=1./complex<double>(0,40.6);

Rasch->msw[38][41]=Rasch->msw[41][38]=2;Rasch->y[38][41]=1./complex<double>(3.37,45);Rasch->y[41][38]=complex<double>(0,.000605);
Rasch->msw[39][40]=Rasch->msw[40][39]=2;Rasch->y[39][40]=1./complex<double>(5.75,62.5);Rasch->y[40][39]=complex<double>(0,.000835);
Rasch->msw[39][41]=Rasch->msw[41][39]=2;Rasch->y[39][41]=1./complex<double>(5,57.4);Rasch->y[41][39]=complex<double>(0,.000765);
Rasch->msw[40][43]=Rasch->msw[43][40]=2;Rasch->y[40][43]=1./complex<double>(10,10.6);Rasch->y[43][40]=complex<double>(0,.001265);
Rasch->msw[41][49]=Rasch->msw[49][41]=2;Rasch->y[41][49]=1./complex<double>(4.08,46.8);Rasch->y[49][41]=complex<double>(0,.000625);
Rasch->msw[42][44]=Rasch->msw[44][42]=2;Rasch->y[42][44]=1./complex<double>(4,88);Rasch->y[44][42]=complex<double>(0,.001528);
Rasch->msw[43][48]=Rasch->msw[48][43]=2;Rasch->y[43][48]=1./complex<double>(8.5,80);Rasch->y[48][43]=complex<double>(0,.001050);
Rasch->msw[43][49]=Rasch->msw[49][43]=2;Rasch->y[43][49]=1./complex<double>(.7,6.5);Rasch->y[49][43]=complex<double>(0,.000352);
Rasch->msw[44][45]=Rasch->msw[45][44]=2;Rasch->y[44][45]=1./complex<double>(4.28,101.8);Rasch->y[45][44]=complex<double>(0,.001751);

//Rasch->msw[44][49]=-1;Rasch->msw[49][44]=1;Rasch->cktr[49][44]=complex<double>(2.1299,0);Rasch->cktr[44][49]=1./Rasch->cktr[49][44];Rasch->y[49][44]=1./complex<double>(0,40.6)/Rasch->cktr[49][44]/Rasch->cktr[49][44];
//Rasch->msw[45][46]=-1;Rasch->msw[46][45]=1;Rasch->cktr[46][45]=complex<double>(2.1299,0);Rasch->cktr[45][46]=1./Rasch->cktr[46][45];Rasch->y[46][45]=1./complex<double>(0,81.2)/Rasch->cktr[46][45]/Rasch->cktr[46][45];
Rasch->msw[44][49]=-1;Rasch->msw[49][44]=1;Rasch->cktr[49][44]=complex<double>(2.1299,0);Rasch->cktr[44][49]=1./Rasch->cktr[49][44];Rasch->y[49][44]=1./complex<double>(0,40.6);
Rasch->msw[45][46]=-1;Rasch->msw[46][45]=1;Rasch->cktr[46][45]=complex<double>(2.1299,0);Rasch->cktr[45][46]=1./Rasch->cktr[46][45];Rasch->y[46][45]=1./complex<double>(0,81.2);

Rasch->msw[45][53]=Rasch->msw[53][45]=2;Rasch->y[45][53]=1./complex<double>(5.2,126.6);Rasch->y[53][45]=complex<double>(0,.002165);
Rasch->msw[46][54]=Rasch->msw[54][46]=2;Rasch->y[46][54]=1./complex<double>(5.6,51);Rasch->y[54][46]=complex<double>(0,.000680);
Rasch->msw[47][48]=Rasch->msw[48][47]=2;Rasch->y[47][48]=1./complex<double>(11.5,110);Rasch->y[48][47]=complex<double>(0,.001440);
Rasch->msw[47][50]=Rasch->msw[50][47]=2;Rasch->y[47][50]=1./complex<double>(8.1,84.4);Rasch->y[50][47]=complex<double>(0,.001178);
Rasch->msw[50][51]=Rasch->msw[51][50]=2;Rasch->y[50][51]=1./complex<double>(.4,4.6);Rasch->y[51][50]=complex<double>(0,.000574);
Rasch->msw[50][52]=Rasch->msw[52][50]=2;Rasch->y[50][52]=1./complex<double>(.4,4.6);Rasch->y[52][50]=complex<double>(0,.000574);
Rasch->msw[50][55]=Rasch->msw[55][50]=2;Rasch->y[50][55]=1./complex<double>(10.7,114);Rasch->y[55][50]=complex<double>(0,.001330);
Rasch->msw[50][57]=Rasch->msw[57][50]=2;Rasch->y[50][57]=1./complex<double>(.32,36);Rasch->y[57][50]=complex<double>(0,.000490);

//Rasch->msw[50][61]=-1;Rasch->msw[61][50]=1;Rasch->cktr[61][50]=complex<double>(2.1701,0);Rasch->cktr[50][61]=1./Rasch->cktr[61][50];Rasch->y[61][50]=1./complex<double>(0,29.6)/Rasch->cktr[61][50]/Rasch->cktr[61][50];
//Rasch->msw[53][51]=-1;Rasch->msw[51][53]=1;Rasch->cktr[51][53]=complex<double>(2.2999,0);Rasch->cktr[53][51]=1./Rasch->cktr[51][53];Rasch->y[51][53]=1./complex<double>(0,81.2)/Rasch->cktr[51][53]/Rasch->cktr[51][53];
//Rasch->msw[53][52]=-1;Rasch->msw[52][53]=1;Rasch->cktr[52][53]=complex<double>(2.2999,0);Rasch->cktr[53][52]=1./Rasch->cktr[52][53];Rasch->y[52][53]=1./complex<double>(0,81.2)/Rasch->cktr[52][53]/Rasch->cktr[52][53];
Rasch->msw[50][61]=-1;Rasch->msw[61][50]=1;Rasch->cktr[61][50]=complex<double>(2.1701,0);Rasch->cktr[50][61]=1./Rasch->cktr[61][50];Rasch->y[61][50]=1./complex<double>(0,29.6);
Rasch->msw[53][51]=-1;Rasch->msw[51][53]=1;Rasch->cktr[51][53]=complex<double>(2.2999,0);Rasch->cktr[53][51]=1./Rasch->cktr[51][53];Rasch->y[51][53]=1./complex<double>(0,81.2);
Rasch->msw[53][52]=-1;Rasch->msw[52][53]=1;Rasch->cktr[52][53]=complex<double>(2.2999,0);Rasch->cktr[53][52]=1./Rasch->cktr[52][53];Rasch->y[52][53]=1./complex<double>(0,81.2);


Rasch->msw[54][55]=Rasch->msw[55][54]=2;Rasch->y[54][55]=1./complex<double>(8.8,83);Rasch->y[55][54]=complex<double>(0,.001110);
Rasch->msw[55][56]=Rasch->msw[56][55]=2;Rasch->y[55][56]=1./complex<double>(4.06,41.2);Rasch->y[56][55]=complex<double>(0,.000550);
Rasch->msw[56][57]=Rasch->msw[57][56]=2;Rasch->y[56][57]=1./complex<double>(6.6,66.6);Rasch->y[57][56]=complex<double>(0,.000890);
Rasch->msw[57][58]=Rasch->msw[58][57]=2;Rasch->y[57][58]=1./complex<double>(10.9,99.3);Rasch->y[58][57]=complex<double>(0,.001260);

//Rasch->msw[57][60]=-1;Rasch->msw[60][57]=1;Rasch->cktr[60][57]=complex<double>(2.1701,0);Rasch->cktr[57][60]=1./Rasch->cktr[60][57];Rasch->y[60][57]=1./complex<double>(0,61.3)/Rasch->cktr[60][57]/Rasch->cktr[60][57];
Rasch->msw[57][60]=-1;Rasch->msw[60][57]=1;Rasch->cktr[60][57]=complex<double>(2.1701,0);Rasch->cktr[57][60]=1./Rasch->cktr[60][57];Rasch->y[60][57]=1./complex<double>(0,61.3);

Rasch->msw[58][59]=Rasch->msw[59][58]=2;Rasch->y[58][59]=1./complex<double>(4.3,46);Rasch->y[59][58]=complex<double>(0,.000540);
Rasch->msw[58][62]=Rasch->msw[62][58]=2;Rasch->y[58][62]=1./complex<double>(0,67.3);
Rasch->msw[59][62]=Rasch->msw[62][59]=2;Rasch->y[59][62]=1./complex<double>(0,86.8);
Rasch->msw[60][61]=Rasch->msw[61][60]=2;Rasch->y[60][61]=1./complex<double>(8,31);Rasch->y[61][60]=complex<double>(0,.000822);
Rasch->YsCreate();
//for(i=0;i<Rasch->n;i++)
//  Rasch->urab[i]=Rasch->unom[i];
//Rasch->NjNbNg();
//double *dX=new double[Rasch->nJac*2+1];
/*
for(i=0;i<NumberOfIters;i++)
  {
//   Rasch->JacobiCalc();
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
  }
*/
//for(i=0;i<Rasch->nJac;i++)
// Rasch->R[i]=-1;
//Rasch->nJac*=2;
/*
for(i=0;i<NumberOfIters;i++)
  {
//   Rasch->JacobiCalc();
   Rasch->JacobiX22Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,Rasch->nJac/2-1,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
  }

if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacRasch;
  }

Rasch->ResultsCalc();
//� �������, ���� ����
//Vals=new Values[1];
if(showresult==true)
  ResultsOut(i,Regim,true);
//delete []Vals;
*/


continueJ=true;
WrongJacRasch:
//delete []dX;
return Regim;
}

//---------------------------------------------------------------------------


bool TForm1::RaschetJEkviv(bool showresult,bool &continueJ)
{
int i,j,k,l,ubazNo;

complex<double> tmpCompl,dU;
complex<double> *tngr=new complex<double> [Rasch->n];
//int krm;//���������� �����, �� ������ ��������� � ��
int kit;//���-�� ��������
double p,q;
//double osh;//������
//krm=0;
//AnsiString tmpStr,tmpStr2;
bool Regim=false;
bool RegOk,Yes;
int state;

Rasch->ysTemp=new complex<double> *[Rasch->n+1];
for(i=0;i<Rasch->n+1;i++)
  Rasch->ysTemp[i]=new complex<double> [Rasch->n+1];
for(i=0;i<Rasch->n;i++)
  if(Rasch->mhu[i]==0)
    {
     ubazNo=i;
     break;
    }

//��������� ������� ������������� Ys
if(Rasch->CurrReg==false)
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
Rasch->YsCreate();
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
kit=0;
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   state=Rasch->Gauss2();
   kit++;
   if(state==0||state==2)
     break;
   state=4;
  }
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacRasch2;
  }
continueJ=true;
for (i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==4||Rasch->mhu[i]==0)
  {
   Rasch->sng[i]=complex<double>(0,0);
   for(j=0;j<Rasch->n;j++)
     Rasch->sng[i]+=Rasch->ys[i][j]*Rasch->urab[j]*conj(Rasch->urab[i]);
   Rasch->sng[i];
  }
}

/*
Rasch->sng[1]=complex<double>(-477,346.5);
Rasch->sng[6]=complex<double>(-440,88);
Rasch->sng[11]=complex<double>(-70,28.3);

Rasch->urab[0]=polar(229.503,-.36*M_PI/180.);
Rasch->urab[2]=polar(227.057,-1.39*M_PI/180.);
Rasch->urab[3]=polar(244.895,14.1*M_PI/180.);
Rasch->urab[4]=polar(39.772,11.86*M_PI/180.);
Rasch->urab[5]=polar(37.861,11.89*M_PI/180.);
Rasch->urab[6]=polar(20.50,18.19*M_PI/180.);
Rasch->urab[7]=polar(34.648,-5.92*M_PI/180.);
Rasch->urab[8]=polar(35.794,-2.36*M_PI/180.);
Rasch->urab[9]=polar(37.720,-1.62*M_PI/180.);
Rasch->urab[10]=polar(221.168,-4.16*M_PI/180.);
Rasch->urab[11]=polar(6.3,1.35*M_PI/180.);
*/

for(i=0;i<Rasch->n;i++)
  for(j=0;j<Rasch->n;j++)
    {
     Rasch->ys[i][j]=complex<double>(0,0);
     Rasch->y[i][j]=complex<double>(0,0);
     Rasch->msw[i][j]=0;
    }

k=ubazNo;

//Memo1->Lines->Add("************************************");
//Memo1->Lines->Add("��������� ����� �����:");

l=0;
j=0;
for(i=0;i<Rasch->n;i++)
  if(Rasch->mhu[i]!=0&&Rasch->sng[i]==complex<double>(0,0))
    j++;
for(i=0;i<Rasch->n;i++)
  {
  if(Rasch->mhu[i]!=0&&Rasch->sng[i]!=complex<double>(0,0))
     {

//      dU=complex<double>(real(Rasch->urab[i]*0.05),fabs(imag(Rasch->urab[i]*0.05)));
      dU=complex<double>(1,0.1);
      p=real(Rasch->sng[i])*(Rasch->a0[i]+(Rasch->a1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->a2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
      q=imag(Rasch->sng[i])*(Rasch->b0[i]+(Rasch->b1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->b2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));

      if(abs(Rasch->urab[i])>abs(Rasch->urab[k]))
        {
            if(p>0)
             {
              Rasch->cktr[k][i]=(Rasch->urab[i]+dU)/Rasch->urab[k];
//              Rasch->cktr[k][i]=abs(Rasch->urab[i]+dU)/abs(Rasch->urab[k]);
              tmpCompl=(Rasch->cktr[k][i]*Rasch->urab[k]-Rasch->urab[i]);
              tmpCompl=polar(real(tmpCompl),imag(tmpCompl));
              tmpCompl=complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=(complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->cktr[k][i]*Rasch->urab[k]-Rasch->urab[i]);

             }
            else if((abs(Rasch->urab[i])-abs(Rasch->urab[k]))/abs(Rasch->urab[i])>0.4)
             {
              dU=dU*Rasch->urab[k]/Rasch->urab[i];
              Rasch->cktr[k][i]=Rasch->urab[i]/(Rasch->urab[k]+dU);
//              Rasch->cktr[k][i]=abs(Rasch->urab[i])/abs(Rasch->urab[k]+dU);
              tmpCompl=(Rasch->urab[i]-Rasch->urab[k]*Rasch->cktr[k][i]);
              tmpCompl=-complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=(-complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[i]-Rasch->urab[k]*Rasch->cktr[k][i]);
             }

            else
             {
              Rasch->cktr[k][i]=complex<double>(1,0);
              tmpCompl=-complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=(-complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[i]-Rasch->urab[k]);

             }

        }
      else if((abs(Rasch->urab[k])-abs(Rasch->urab[i]))/abs(Rasch->urab[k])>0.4)
        {
            trf1:
            if(p>0)
             {
              dU=dU*Rasch->urab[i]/Rasch->urab[k];
              Rasch->cktr[k][i]=(Rasch->urab[i]+dU)/(Rasch->urab[k]);
//              Rasch->cktr[k][i]=abs(Rasch->urab[i]+dU)/abs(Rasch->urab[k]);
              tmpCompl=complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=Rasch->urab[k]-Rasch->urab[i]/Rasch->cktr[k][i];
              tmpCompl=(Rasch->cktr[k][i]*complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[k]-Rasch->urab[i]/Rasch->cktr[k][i]);

//              tmpCompl=Rasch->cktr[k][i]*Rasch->urab[k]-Rasch->urab[i];
//              tmpCompl=(complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->cktr[k][i]*Rasch->urab[k]-Rasch->urab[i]);
             }

            else
             {
//              dU=dU*Rasch->urab[i]/Rasch->urab[k];
              Rasch->cktr[k][i]=(Rasch->urab[i])/(Rasch->urab[k]+dU);
//              Rasch->cktr[k][i]=abs(Rasch->urab[i]+dU)/abs(Rasch->urab[k]);
              tmpCompl=complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=Rasch->urab[k]-Rasch->urab[i]/Rasch->cktr[k][i];
              tmpCompl=(-Rasch->cktr[k][i]*complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[i]/Rasch->cktr[k][i]-Rasch->urab[k]);
             }

        }
      else
        {
//         dU=Rasch->urab[k]*0.1;
         Rasch->cktr[k][i]=complex<double>(1,0);
            if(p>0)
             {
              tmpCompl=complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=Rasch->urab[k]-Rasch->urab[i];
              tmpCompl=(complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[k]-Rasch->urab[i]);
             }
            else
             {
              tmpCompl=-complex<double>(p,q)/conj(Rasch->urab[i]);
              tmpCompl=Rasch->urab[k]-Rasch->urab[i];
              tmpCompl=(-complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[i]-Rasch->urab[k]);
             }
//            if (real(tmpCompl)<0)
//              goto trf1;
        }
      Rasch->cktr[i][k]=1./Rasch->cktr[k][i];
      if(abs(Rasch->cktr[k][i])>1)
        {
         Rasch->ys[k][k]-=tmpCompl*Rasch->cktr[k][i]*Rasch->cktr[k][i];
         Rasch->ys[i][i]-=tmpCompl;
         Rasch->ys[i][k]=Rasch->ys[k][i]=tmpCompl*Rasch->cktr[k][i];
        }
      else
        {
/*
         Rasch->ys[k][k]-=tmpCompl*Rasch->cktr[k][i]*Rasch->cktr[k][i];
         Rasch->ys[i][i]-=tmpCompl;
         Rasch->ys[i][k]=Rasch->ys[k][i]=tmpCompl*Rasch->cktr[k][i];
*/


         Rasch->ys[i][i]-=tmpCompl/(Rasch->cktr[k][i]*Rasch->cktr[k][i]);
         Rasch->ys[k][k]-=tmpCompl;
         Rasch->ys[i][k]=Rasch->ys[k][i]=tmpCompl/Rasch->cktr[k][i];

        }

      if(abs(Rasch->cktr[k][i])>1)
        Rasch->y[k][i]=tmpCompl;
      else if (abs(Rasch->cktr[k][i])<1)
        Rasch->y[i][k]=tmpCompl;
      else
        Rasch->y[min(k,i)][max(k,i)]=tmpCompl;

      tmpCompl=1./Rasch->ys[i][k];
//      Memo1->Lines->Add(AnsiString(k+1)+ ";"+ AnsiString(i+1)+ ";"+AnsiString(real(Rasch->cktr[k][i]))+";"+AnsiString(imag(Rasch->cktr[k][i]))+";"+AnsiString(real(tmpCompl))+";"+AnsiString(imag(tmpCompl)));
/*
Rasch->cktr[k][i]=Rasch->cktr[i][k]=1;
/*
            if(p>0)
             {
              tmpCompl=complex<double>(p,q)/conj(Rasch->urab[i]);
//              tmpCompl=Rasch->urab[k]-Rasch->urab[i];
              tmpCompl=(complex<double>(p,q)/conj(Rasch->urab[i]))/(Rasch->urab[k]-Rasch->urab[i]);
             }


      Rasch->ys[k][k]-=tmpCompl*Rasch->cktr[k][i]*Rasch->cktr[k][i];
      Rasch->ys[i][i]-=tmpCompl;
      Rasch->ys[i][k]=Rasch->ys[k][i]=tmpCompl*Rasch->cktr[k][i];
      tmpCompl=1./Rasch->ys[i][k];
      Memo1->Lines->Add(AnsiString(k+1)+ ";"+ AnsiString(i+1)+ ";"+AnsiString(real(Rasch->cktr[k][i]))+";"+AnsiString(imag(Rasch->cktr[k][i]))+";"+AnsiString(real(tmpCompl))+";"+AnsiString(imag(tmpCompl)));
*/
     Rasch->sng[l]=Rasch->sng[i];
     Rasch->urab[l]=Rasch->urab[i];
     Rasch->unom[l]=Rasch->unom[i];
     Rasch->a0[l]=Rasch->a0[i];
     Rasch->a1[l]=Rasch->a1[i];
     Rasch->a2[l]=Rasch->a2[i];
     Rasch->b0[l]=Rasch->b0[i];
     Rasch->b1[l]=Rasch->b1[i];
     Rasch->b2[l]=Rasch->b2[i];
     Rasch->mhu[l]=Rasch->mhu[i];
     Rasch->cktr[l][k]=Rasch->cktr[i][k];
     Rasch->cktr[k][l]=Rasch->cktr[k][i];
     if(abs(Rasch->cktr[l][k])>1)
        {
         Rasch->msw[l][k]=1;
         Rasch->msw[k][l]=-1;
        }
     else if(abs(Rasch->cktr[l][k])<1)
        {
         Rasch->msw[l][k]=-1;
         Rasch->msw[k][l]=1;
        }
     else
         Rasch->msw[l][k]=Rasch->msw[k][l]=2;

     Rasch->ys[l][l]=Rasch->ys[i][i];
     Rasch->ys[l][k]=Rasch->ys[k][l]=Rasch->ys[i][k];
     Rasch->y[l][k]=Rasch->y[i][k];
     Rasch->y[k][l]=Rasch->y[k][i];

//     Rasch->y[min(l,k)][max(l,k)]=Rasch->y[min(i,k)][max(i,k)];
     l++;




     }
  else if(Rasch->mhu[i]==0)
    {
     l++;
    }
  }
//Memo1->Lines->Add("*******************");
Rasch->n-=j;


for(i=0;i<Rasch->n;i++)
  {
   Rasch->urab[i]=Rasch->unom[i];
  }





//if(showresult==true)
//  ResultsOut(kit,Regim);

for(i=0;i<Rasch->n+1;i++)
  {
   delete []Rasch->ysTemp[i];
  }

kit=0;

Rasch->NjNbNg();
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   state=Rasch->Gauss2();
   kit++;
   if(state==0||state==2)
     break;
   state=4;
  }
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacRasch2;
  }
continueJ=true;
Rasch->ResultsCalc();
//Vals=new Values[1];
if(showresult==true)
  ResultsOut(kit,Regim,true);
//delete []Vals;

WrongJacRasch2:
//Rasch->YsShow(Memo1);
delete []tngr;
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------
bool TForm1::RaschetJX2(bool showresult,bool &continueJ)
{
int i,state,kit,k=0;
bool Regim,limits=false;
double minJ,minJNo,Kn;
Rasch->YsCreate();
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);

double *dX=new double [Rasch->n*4];
for(i=0;i<Rasch->n;i++)
  Rasch->mhu_n[i]=Rasch->mhu[i];


k=0;
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
   Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zu[i]=k;
   k++;
  }
}
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0)
   Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zf[i]=k;
   k++;
  }
}
kit=0;
minJNo=Rasch->nJac-1;

Rasch->NjNbNg();
Rasch->Jacobi2Calc();
Rasch->RaschR0(Rasch->Jac,Rasch->nJac,Rasch->R);
Rasch->nJac*=2;
Rasch->InitializeJac(Rasch->nJac);
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();
   if(Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1)==0)//������ � ���������� ������
     {
      state=2;
      break;
     }
   if(Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax)==1)//����� �������!
     {
      state=0;
      break;
     }
   state=4;
   kit++;
  }

if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX2Rasch;
  }


if(UseQminQmax==true)
if(Rasch->Limits(0,false,0)==true)
 {
  Rasch->NjNbNg();
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }
/*
 minJ=10000000000000000;
 minJNo=0;

 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }


*/
 minJNo=Rasch->nJac-1;
// Rasch->NjNbNg();
// Rasch->Jacobi2Calc();
// Rasch->RaschR0(Rasch->Jac,Rasch->nJac,Rasch->R);
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac);
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
if(UseQminQmax==true)
if(Rasch->Limits(1,false,0)==true)
 {
  Rasch->NjNbNg();
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }
/*
 minJ=10000000000000000;
 minJNo=0;
 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }
 Rasch->nJac*=2;
*/
 minJNo=Rasch->nJac-1;
// Rasch->NjNbNg();
// Rasch->Jacobi2Calc();
// Rasch->RaschR0(Rasch->Jac,Rasch->nJac,Rasch->R);
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac);
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX2Rasch;
  }

continueJ=true;
ShowMessage(Rasch->Kn);
Rasch->ResultsCalc();
ResultsOut(kit,Regim,showresult);



WrongJacX2Rasch:

for(i=0;i<Rasch->n;i++)
  Rasch->mhu_n[i]=Rasch->mhu[i];


delete []dX;
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------
bool TForm1::RaschetJX22(bool showresult,bool &continueJ,bool WritePrmsToFile,AnsiString PrmsFileName)  //������ � ����������� �������� �����
{
int i,k,state,kit=0;
bool Regim,sobstvOk,limits=false;
double minJ,minJNo;
Rasch->YsCreate();
/*
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
*/  
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
Rasch->Jacobi2Calc();
double *dX=new double [Rasch->n*4];
/*
for(i=0;i<NumberOfIters;i++)
  {

   Rasch->SHN(true);
   Rasch->Jacobi_KnCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
  }
*/
//if(state!=0)
Regim=Rasch->RegimState;
if(Regim==true)
 state=0;
else
 state=4;

if(Regim==false&&Rasch->CurrReg==false)
{
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
for(i=0;i<NumberOfIters;i++)
  {

   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   if(state==0||state==2)
     break;
   state=4;
  }
}
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
//   goto WrongJacX22Rasch;
  }
k=0;
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
   Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zu[i]=k;
   k++;
  }
}
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0)
   Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zf[i]=k;
   k++;
  }
}

minJ=10000000000000000;
minJNo=0;
Rasch->Jacobi2Calc();
for(i=0;i<Rasch->nJac;i++)
  if(Rasch->Jac[i][i]<minJ)
    {
     minJ=Rasch->Jac[i][i];
     minJNo=i;
    }

if(Regim==false)
 for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];

minJNo=Rasch->nJac-1;
Rasch->nJac*=2;
Rasch->nJac+=1;
Rasch->InitializeJac(Rasch->nJac);

ap::real_2d_array tmpm;
tmpm.setbounds(1,Rasch->nJac,1,Rasch->nJac);
double determUR;
double determUPR;
long double det2nUR;
long double det2nUPR;
long double determUPRFull;
long double determURFull;
double obuslovl;
double NormZ;
double NormH=0;
double NormF=0;
double NormV=0;
double NormUs=0;
double Kntemp=Rasch->Kn;
AnsiString MaxNebalNo,MaxSummNo;
AnsiString tmpString;
for(i=0;i<Rasch->nJac;i++)
 dX[i]=0;

TStringList *tmpPrmsList=new TStringList;
try
{
 tmpPrmsList->Add("� ����;||H||;||F||;||V||;U;det(dF/dX)*;10^x;det(dH/dZ)*;10^y;||dZ||;Kn;Obusl;");
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();

   //���� ����� ������������ ������������� ��������� ��� ������ � ����
   if(WritePrmsToFile==true)
   {
    for(int ii=0;ii<Rasch->nJac;ii++)
    for(int jj=0;jj<Rasch->nJac;jj++)
      tmpm(ii+1,jj+1)=Rasch->Jac[ii][jj];
    determUPR=determinant(tmpm,Rasch->nJac);
    det2nUPR=powl(100,Rasch->nJac);
    determUPRFull=(long double) (det2nUPR)*(long double) (determUPR);

    Rasch->nJac=(Rasch->nJac-1)/2;
    Rasch->Jacobi2Calc();
    for(int ii=0;ii<Rasch->nJac;ii++)
    for(int jj=0;jj<Rasch->nJac;jj++)
      tmpm(ii+1,jj+1)=Rasch->Jac[ii][jj];
    determUR=determinant(tmpm,Rasch->nJac);
    det2nUR=powl(100,Rasch->nJac);
    determURFull=(long double) (det2nUR)*(long double) (determUR);
    Rasch->nJac*=2;
    Rasch->nJac+=1;
    Rasch->JacobiX22Calc();
    obuslovl=Rasch->ObuslMatr(Rasch->Jac,Rasch->nJac);

    if(obuslovl!=0)
     obuslovl=1./obuslovl;
    else
     obuslovl=-1;

    NormH=0;
    NormF=0;
    NormV=0;
    for(int ii=0;ii<Rasch->nJac;ii++)
     {
      if(ii<(Rasch->nJac-1)/2)
       NormF+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
      else if(ii<Rasch->nJac-1)
       NormV+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
      NormH+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
     }

    NormH=sqrt(NormH);
    NormF=sqrt(NormF);
    NormV=sqrt(NormV);
    NormUs=sqrt(Rasch->Nebal[Rasch->nJac-1]*Rasch->Nebal[Rasch->nJac-1]);

    NormZ=0;
    for(int ii=0;ii<Rasch->nJac;ii++)
      {
       NormZ+=dX[ii]*dX[ii];
      }
     NormZ=sqrt(NormZ);

   }
   //***********************

   if(Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1)==0)//������ � ���������� ������
     {
      state=2;
      break;
     }
   //���� ����� ������������ ������������� ��������� ��� ������ � ����
   if(WritePrmsToFile==true)
    {
    Kntemp=Rasch->Kn;
    }
    //***************************************


   if(Rasch->RaschCheck_new(false,0,dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax)==1)//����� �������!
     {
      state=0;
//      break;
     }
   else
     state=4;

  //���� ����� ������������ ������������� ��������� ��� ������ � ����
  if(WritePrmsToFile==true)
   {
   tmpString=AnsiString(i)+";";
   tmpString+=FloatToStrF(NormH,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormF,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormV,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormUs,ffFixed,11,7)+";";

   tmpString+=FloatToStrF(determUR,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(det2nUR,ffFixed,7,4)+";";

   tmpString+=FloatToStrF(determUPR,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(det2nUPR,ffFixed,7,4)+";";

   tmpString+=FloatToStrF(NormZ,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(Kntemp,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(obuslovl,ffFixed,7,4)+";";
   tmpPrmsList->Add(tmpString);
   }
  //**********************************
   if(state==0)
    break;
   kit++;
  }
  if(WritePrmsToFile==true)
   tmpPrmsList->SaveToFile(PrmsFileName);
  //***********************************
}
__finally
{
delete tmpPrmsList;
}



if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX22Rasch;
  }

for(i=0;i<Rasch->n;i++)
  Rasch->mhu_n[i]=Rasch->mhu[i];
if(UseQminQmax==true)
if(Rasch->Limits(0,false,0)==true)
 {
  Rasch->NjNbNg();
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }

 minJ=10000000000000000;
 minJNo=0;
 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }
 minJNo=Rasch->nJac-1;
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac);
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
 }



if(UseQminQmax==true)
if(Rasch->Limits(1,false,0)==true)
 {
  Rasch->NjNbNg();
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }


 minJ=10000000000000000;
 minJNo=0;
 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }
 minJNo=Rasch->nJac-1;
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac); 
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX22Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,1,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX22Rasch;
  }
continueJ=true;

for(i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];

 Rasch->ResultsCalc();
 ResultsOut(kit,Regim,showresult);




WrongJacX22Rasch:
delete []dX;
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------
bool TForm1::RaschetJX2T(bool showresult,bool &continueJ,bool WritePrmsToFile,AnsiString PrmsFileName)  //������ � ����������� �������� ����� � ���������� �
{
int i,k,state,kit=0;
bool Regim,sobstvOk,limits=false;
double minJ,minJNo;
Rasch->YsCreate();
/*
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
*/
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
Rasch->Jacobi2Calc();
double *dX=new double [Rasch->n*4];

Regim=Rasch->RegimState;
if(Regim==true)
 state=0;
else
 state=4;
if(Regim==false)
{
for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck_new(false,0,dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
//                                BkUse,Bk,dX,Rasch->urab,Rasch->nJac,minJNo,2,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax)==1)
   if(state==0||state==2)
     break;
   state=4;
  }
}
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
//   goto WrongJacX22Rasch;
  }
k=0;
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
   Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zu[i]=k;
   k++;
  }
}
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==0)
   Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
 else
  {
   Rasch->Zf[i]=k;
   k++;
  }
}

minJ=10000000000000000;
minJNo=0;
Rasch->Jacobi2Calc();
for(i=0;i<Rasch->nJac;i++)
  if(Rasch->Jac[i][i]<minJ)
    {
     minJ=Rasch->Jac[i][i];
     minJNo=i;
    }

if(Regim==false)
 for(i=0;i<Rasch->n;i++)
  Rasch->urab[i]=Rasch->unom[i];


minJNo=Rasch->nJac-1;

Rasch->nJac*=2;
Rasch->nJac+=1;
Rasch->InitializeJac(Rasch->nJac);
/*
//Memo1->Lines->Add("����� ��������");

Memo1->Lines->Add("Z0");

for(i=0;i<Rasch->n;i++)
 if(Rasch->mhu[i]!=0&&Rasch->mhu[i]!=4)
  Memo1->Lines->Add(AnsiString(abs(Rasch->urab[i])));
for(i=0;i<Rasch->n;i++)
 if(Rasch->mhu[i]!=0)
  Memo1->Lines->Add(AnsiString(arg(Rasch->urab[i])));

for(i=0;i<(Rasch->nJac-1)/2;i++)
  Memo1->Lines->Add(AnsiString(Rasch->R[i]));

Memo1->Lines->Add(AnsiString(Rasch->T));

Memo1->Lines->Add("End-Z");
*/
GesseParams *GP=new GesseParams(2,Rasch->nJac,Rasch->n,Rasch->mhu);
double Bk=0;
bool BkUse=false;
AnsiString BkStr;

ap::real_2d_array tmpm;
tmpm.setbounds(1,Rasch->nJac,1,Rasch->nJac);
double determUR;
double determUPR;
long double det2nUR;
long double det2nUPR;
long double determUPRFull;
long double determURFull;
double obuslovl;
double NormZ;
double dY;
double NormH=0;
double NormF=0;
double NormV=0;
double NormUs=0;
double Ttemp=Rasch->T;
AnsiString MaxNebalNo,MaxSummNo;
AnsiString tmpString;
//AnsiString PnStr,PgStr,QnStr,QgStr,UStr,fStr;
for(i=0;i<Rasch->nJac;i++)
 dX[i]=0;

TStringList *tmpPrmsList=new TStringList;
try
{
tmpPrmsList->Add("� ����;||H||;||F||;||V||;U;det(dF/dX)*;10^x;det(dH/dZ)*;10^y;||dZ||;T;||dY||;Obusl;Bk;MaxNebNo;MaxGesseNo;");
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX2TCalc();

   //���� ����� ������������ ������������� ��������� ��� ������ � ����
   if(WritePrmsToFile==true)
   {
    for(int ii=0;ii<Rasch->nJac;ii++)
    for(int jj=0;jj<Rasch->nJac;jj++)
      tmpm(ii+1,jj+1)=Rasch->Jac[ii][jj];
    determUPR=determinant(tmpm,Rasch->nJac);
    det2nUPR=powl(100,Rasch->nJac);
    determUPRFull=(long double) (det2nUPR)*(long double) (determUPR);
    determUR=determinant(tmpm,(Rasch->nJac-1)/2);
    det2nUR=powl(100,(Rasch->nJac-1)/2);
    determURFull=(long double) (det2nUR)*(long double) (determUR);
    obuslovl=Rasch->ObuslMatr(Rasch->Jac,Rasch->nJac);
    if(obuslovl!=0)
     obuslovl=1./obuslovl;
    else
     obuslovl=-1;

    NormH=0;
    NormF=0;
    NormV=0;
    for(int ii=0;ii<Rasch->nJac;ii++)
     {
      if(ii<(Rasch->nJac-1)/2)
       NormF+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
      else if(ii<Rasch->nJac-1)
       NormV+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
      NormH+=Rasch->Nebal[ii]*Rasch->Nebal[ii];
     }

    NormH=sqrt(NormH);
    NormF=sqrt(NormF);
    NormV=sqrt(NormV);
    NormUs=sqrt(Rasch->Nebal[Rasch->nJac-1]*Rasch->Nebal[Rasch->nJac-1]);

    NormZ=0;
    for(int ii=0;ii<Rasch->nJac;ii++)
      {
       NormZ+=dX[ii]*dX[ii];
      }
     NormZ=sqrt(NormZ);

    BkStr=FloatToStrF(Bk,ffFixed,7,4)+";";
   MaxNebalNo=FloatToStrF(GP->MaxNebalNo,ffFixed,7,4)+";";
   MaxSummNo=FloatToStrF(GP->MaxSummNo,ffFixed,7,4)+";";
   }
   //***********************
   if(Rasch->LU(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac)==false)//������ � ���������� ������
     {
      state=2;
      break;
     }
   if(Rasch->BkUse==true)
   {
    Bk=Rasch->BkCalc(GP,dX,Rasch->Nebal);
    if(Bk>=1)
     BkUse=true;
    else
     BkUse=false;
   }

  //���� ����� ������������ ������������� ��������� ��� ������ � ����
  if(WritePrmsToFile==true)
   {
   dY=0;
   for(int ii=0;ii<Rasch->n;ii++)
    {
      dY+=(Rasch->dPn[ii]-Rasch->dPg[ii])*(Rasch->dPn[ii]-Rasch->dPg[ii])*(Rasch->T-Ttemp)*(Rasch->T-Ttemp);
    }
   dY=sqrt(dY);

   Ttemp=Rasch->T;
   }
   //***************************************
   if(Rasch->RaschCheck_new(BkUse,Bk,dX,Rasch->urab,Rasch->nJac,minJNo,2,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax)==1)//����� �������!
     {
      state=0;
     }
   else
     state=4;

  //���� ����� ������������ ������������� ��������� ��� ������ � ����
  if(WritePrmsToFile==true)
   {
   tmpString=AnsiString(i)+";";
   tmpString+=FloatToStrF(NormH,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormF,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormV,ffFixed,11,7)+";";
   tmpString+=FloatToStrF(NormUs,ffFixed,11,7)+";";

   tmpString+=FloatToStrF(determUR,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(det2nUR,ffFixed,7,4)+";";

   tmpString+=FloatToStrF(determUPR,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(det2nUPR,ffFixed,7,4)+";";

   tmpString+=FloatToStrF(NormZ,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(Ttemp,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(dY,ffFixed,7,4)+";";
   tmpString+=FloatToStrF(obuslovl,ffFixed,7,4)+";";
   tmpString+=BkStr;
   tmpString+=MaxNebalNo;
   tmpString+=MaxSummNo;
   tmpPrmsList->Add(tmpString);
   }
  //**********************************
   if(state==0)
    break;
   kit++;
  }
  //���� ����� ���������� ������������� ��������� � ����
  if(WritePrmsToFile==true)
   tmpPrmsList->SaveToFile(PrmsFileName);
  //***********************************
}

__finally
{
delete tmpPrmsList;
}

delete GP;
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX22Rasch;
  }

for(i=0;i<Rasch->n;i++)
  Rasch->mhu_n[i]=Rasch->mhu[i];
if(UseQminQmax==true)
if(Rasch->Limits(0,true,Rasch->T)==true)
 {
  Rasch->NjNbNg();
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }

 minJ=10000000000000000;
 minJNo=0;
 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }
 minJNo=Rasch->nJac-1;
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac);
 //Rasch->Kn=Rasch->T=15;
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX2TCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);

   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,2,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}

if(UseQminQmax==true)
if(Rasch->Limits(1,true,Rasch->T)==true)
 {
  Rasch->NjNbNg();
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 k=0;
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
    Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zu[i]=k;
    k++;
   }
 }
 for(i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]==0)
    Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
  else
   {
    Rasch->Zf[i]=k;
    k++;
   }
 }


 minJ=10000000000000000;
 minJNo=0;
 Rasch->Jacobi2Calc();
 for(i=0;i<Rasch->nJac;i++)
   if(Rasch->Jac[i][i]<minJ)
     {
      minJ=Rasch->Jac[i][i];
      minJNo=i;
     }
 minJNo=Rasch->nJac-1;
 Rasch->nJac*=2;
 Rasch->InitializeJac(Rasch->nJac);
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN(true);
   Rasch->JacobiX2TCalc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,minJNo,2,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
//Memo1->Lines->Add("********");
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto WrongJacX22Rasch;
  }
continueJ=true;

for(i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];

 Rasch->ResultsCalc();
 ResultsOut(kit,Regim,showresult);



//ShowMessage(Rasch->T);

WrongJacX22Rasch:
delete []dX;
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------
bool TForm1::RaschetJSimpleT(bool showresult,bool &continueJ)
{
int i,j,k,state;
float d,T,T0,Tpred;
bool flag,Regim;
bool ok;
int kit=0;

double *dX=new double [Rasch->n*2];
Rasch->YsCreate();
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
//Rasch->mhu_n=new int [Rasch->n];
for(i=0;i<Rasch->n;i++)
 Rasch->mhu_n[i]=Rasch->mhu[i];
for(i=0;i<Rasch->n;i++)
   Rasch->urab[i]=Rasch->unom[i];
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
if(state==0)
 flag=true;
else
 flag=false;
float maxdPQ=0;
int maxdPQNumber=-1;
int maxIs=-1;
for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dPg[i])>0&&fabs(Rasch->Pg[i])>maxdPQ)
  {maxdPQ=fabs(Rasch->Pg[i]);
   maxdPQNumber=i;
   maxIs=0;
  }
for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dPn[i])>0&&fabs(real(Rasch->sng[i]))>maxdPQ)
  {maxdPQ=fabs(real(Rasch->sng[i]));
   maxdPQNumber=i;
   maxIs=1;
  }

for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dQn[i])>0&&fabs(imag(Rasch->sng[i]))>maxdPQ)
  {maxdPQ=fabs(imag(Rasch->sng[i]));
   maxdPQNumber=i;
   maxIs=2;
  }
if(maxdPQNumber!=-1)
{
if(maxIs==0)
 T0=fabs(Rasch->Pg[maxdPQNumber])/fabs(Rasch->dPg[maxdPQNumber]);
else if(maxIs==1)
 T0=fabs(real(Rasch->sng[maxdPQNumber]))/fabs(Rasch->dPn[maxdPQNumber]);
else
 T0=fabs(imag(Rasch->sng[maxdPQNumber]))/fabs(Rasch->dQn[maxdPQNumber]);
}

ok=false;
d=2.;
T=T0;
//flag=true;
float CurrPrec=0.01;
kit=0;
bool limits=false;
while(ok==false&&maxdPQNumber!=-1)
{
if(flag==false)
  for(i=0;i<Rasch->n;i++)
   Rasch->urab[i]=Rasch->unom[i];
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
if(UseQminQmax==true)
if(Rasch->Limits(0,true,T)==true)
 {
  limits=true;
  if(state!=0)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
//    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
 }
if(UseQminQmax==true)
if(Rasch->Limits(1,true,T)==true)
 {
  limits=true;
  if(state!=0)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
//    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
 }
if(state==0)//������� �����
 {
  T+=T0/d;
  if(fabs(Tpred-T)<CurrPrec)
    {ok=true;T-=T0/d;Regim=true;}
  else
  {
  if(flag==false)
    d*=2;
  Tpred=T;
  flag=true;
  }
 }
else
 {
  T-=T0/d;
  if(flag==true)
    d*=2;
  Tpred=T;
  flag=false;
 }
}
/*
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
if(state==4)
  Regim=false;
else if(state==0)
  Regim=true;
else if(state==2)
  {
   continueJ=false;
   goto Wrong;
  }
continueJ=true;
*/
Rasch->ResultsCalc();
ResultsOut(kit,Regim,showresult);
HardItersNumber=T;
Rasch->T=T;
for(i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];


delete []dX;
//delete []Rasch->mhu_n;
ShowMessage(T);
Wrong:
Rasch->RegimState=Regim;
return Regim;
}
//---------------------------------------------------------------------------





//---------------------------------------------------------------------------
void TForm1::ResultsOut(int kit,bool Regim,bool Showtxt)
{
AnsiString tmpStr;
int tmpNode1,tmpNode2;
int i;
TStringList *TmpList=new TStringList;
TStringList *TmpList2=new TStringList;
try
{
TmpList->Add("");
TmpList->Add("������� ����� �����: "+AnsiString(CurrIt));
TmpList->Add("���������� � �����");
TmpList->Add("i    U������      U����");
Vals[CurrIt].NodesNumber=Rasch->n;
Vals[CurrIt].mswNumber=Rasch->ParNumber;
Vals[CurrIt].RegimOk=Regim;
Iters=kit;

for(i=0;i<Rasch->n;i++)
 {
  if(RealNodes==false)
    {
     tmpStr=AnsiString(i+1)+" "+FloatToStrF(abs(Rasch->urab[i]),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->urab[i])*180./M_PI,ffFixed,10,4);
     if(Rasch->Nu[i]!=-1)
       Vals[CurrIt].RealNode[i]=Rasch->Nu[i];
     else
       Vals[CurrIt].RealNode[i]=i+1;
    }
  else
    {
     if(Rasch->Nu[i]!=-1)
      {
      tmpStr=AnsiString(Rasch->Nu[i])+" "+FloatToStrF(abs(Rasch->urab[i]),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->urab[i])*180./M_PI,ffFixed,10,4);
      Vals[CurrIt].RealNode[i]=Rasch->Nu[i];
      }
     else
       {
       tmpStr=AnsiString(i+1)+" "+FloatToStrF(abs(Rasch->urab[i]),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->urab[i])*180./M_PI,ffFixed,10,4);
       Vals[CurrIt].RealNode[i]=i+1;
       }
    }
  Vals[CurrIt].u[i]=Rasch->urab[i];
  TmpList->Add(tmpStr);
}

TmpList->Add("");
TmpList->Add("���������� ��������: "+AnsiString(kit));
if(Regim==false)
   TmpList->Add("����� �� �������!!");
TmpList->Add("");
TmpList->Add("���� � ������ �����");
TmpList->Add("i,j,     I���,     I�����,    I���,    I���");

for(i=0;i<Rasch->ParNumber;i++)
   {
    if(RealNodes==false)
     {
      tmpStr=AnsiString(Rasch->Params[i].Node1)+","+AnsiString(Rasch->Params[i].Node2)+" "+FloatToStrF(real(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(abs(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->Params[i].I)*180./M_PI,ffFixed,10,4);
      Vals[CurrIt].Node1[i]=Rasch->Params[i].Node1;
      Vals[CurrIt].Node2[i]=Rasch->Params[i].Node2;
     }
    else
     {

      if(Rasch->Nu[Rasch->Params[i].Node1-1]!=-1)
       tmpNode1=Rasch->Nu[Rasch->Params[i].Node1-1];
      else
       tmpNode1=Rasch->Params[i].Node1;
      if(Rasch->Nu[Rasch->Params[i].Node2-1]!=-1)
       tmpNode2=Rasch->Nu[Rasch->Params[i].Node2-1];
      else
       tmpNode2=Rasch->Params[i].Node2;

      //tmpStr=AnsiString(Rasch->Nu[Rasch->Params[i].Node1-1])+","+AnsiString(Rasch->Nu[Rasch->Params[i].Node2-1])+" "+FloatToStrF(real(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(abs(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->Params[i].I)*180./M_PI,ffFixed,10,4);
      tmpStr=AnsiString(tmpNode1)+","+AnsiString(tmpNode2)+" "+FloatToStrF(real(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(abs(Rasch->Params[i].I),ffFixed,10,4)+'\t'+FloatToStrF(arg(Rasch->Params[i].I)*180./M_PI,ffFixed,10,4);
//      Vals[CurrIt].Node1[i]=tmpNode1;
//      Vals[CurrIt].Node2[i]=tmpNode2;
      Vals[CurrIt].Node1[i]=Rasch->Params[i].Node1;
      Vals[CurrIt].Node2[i]=Rasch->Params[i].Node2;

     }
    Vals[CurrIt].I[i]=Rasch->Params[i].I;
    TmpList->Add(tmpStr);
   }

TmpList->Add("");
TmpList->Add("�������� �������� � ������ �����");
TmpList->Add("i,j,     Pi,     Qi,    Pj,    Qj");
for(i=0;i<Rasch->ParNumber;i++)
   {
    if(RealNodes==false)
      {
       tmpStr=AnsiString(Rasch->Params[i].Node1)+","+AnsiString(Rasch->Params[i].Node2)+" "+FloatToStrF(real(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].SFlow2),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow2),ffFixed,10,4);
      }
    else
      {
      if(Rasch->Nu[Rasch->Params[i].Node1-1]!=-1)
       tmpNode1=Rasch->Nu[Rasch->Params[i].Node1-1];
      else
       tmpNode1=Rasch->Params[i].Node1;
      if(Rasch->Nu[Rasch->Params[i].Node2-1]!=-1)
       tmpNode2=Rasch->Nu[Rasch->Params[i].Node2-1];
      else
       tmpNode2=Rasch->Params[i].Node2;
//       tmpStr=AnsiString(Rasch->Nu[Rasch->Params[i].Node1-1])+","+AnsiString(Rasch->Nu[Rasch->Params[i].Node2-1])+" "+FloatToStrF(real(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].SFlow2),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow2),ffFixed,10,4);
       tmpStr=AnsiString(tmpNode1)+","+AnsiString(tmpNode2)+" "+FloatToStrF(real(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow1),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].SFlow2),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].SFlow2),ffFixed,10,4);
      }
    Vals[CurrIt].S1[i]=Rasch->Params[i].SFlow1;
    Vals[CurrIt].S2[i]=Rasch->Params[i].SFlow2;
    TmpList->Add(tmpStr);
   }
TmpList->Add("");
TmpList->Add("������ �������� � ���������� � ���");
TmpList->Add("i,j,     dP,     dQ,    dU");
for(i=0;i<Rasch->ParNumber;i++)
   if(Rasch->Params[i].IsLine==true)
   {
    if(RealNodes==false)
      tmpStr=AnsiString(Rasch->Params[i].Node1)+","+AnsiString(Rasch->Params[i].Node2)+" "+FloatToStrF(real(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(Rasch->Params[i].DU,ffFixed,10,4);
    else
      {
      if(Rasch->Nu[Rasch->Params[i].Node1-1]!=-1)
       tmpNode1=Rasch->Nu[Rasch->Params[i].Node1-1];
      else
       tmpNode1=Rasch->Params[i].Node1;
      if(Rasch->Nu[Rasch->Params[i].Node2-1]!=-1)
       tmpNode2=Rasch->Nu[Rasch->Params[i].Node2-1];
      else
       tmpNode2=Rasch->Params[i].Node2;

//      tmpStr=AnsiString(Rasch->Nu[Rasch->Params[i].Node1-1])+","+AnsiString(Rasch->Nu[Rasch->Params[i].Node2-1])+" "+FloatToStrF(real(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(Rasch->Params[i].DU,ffFixed,10,4);
      tmpStr=AnsiString(tmpNode1)+","+AnsiString(tmpNode2)+" "+FloatToStrF(real(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DS),ffFixed,10,4)+'\t'+FloatToStrF(Rasch->Params[i].DU,ffFixed,10,4);
      }
    Vals[CurrIt].DS[i]=Rasch->Params[i].DS;
    Vals[CurrIt].DU[i]=Rasch->Params[i].DU;
    Vals[CurrIt].Ql[i]=Rasch->Params[i].Ql;
    Vals[CurrIt].IsLine[i]=true;
    TmpList->Add(tmpStr);
   }
TmpList->Add("");
TmpList->Add("������ �������� � ���������������");
TmpList->Add("i,j,     dPTM,     dQTM,    dPTC,    dPQC");
for(i=0;i<Rasch->ParNumber;i++)
   if(Rasch->Params[i].IsTrf==true)
   {
    if(RealNodes==false)
      tmpStr=AnsiString(Rasch->Params[i].Node1)+","+AnsiString(Rasch->Params[i].Node2)+" "+FloatToStrF(real(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].DSTC),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTC),ffFixed,10,4);
    else
      {
      if(Rasch->Nu[Rasch->Params[i].Node1-1]!=-1)
       tmpNode1=Rasch->Nu[Rasch->Params[i].Node1-1];
      else
       tmpNode1=Rasch->Params[i].Node1;
      if(Rasch->Nu[Rasch->Params[i].Node2-1]!=-1)
       tmpNode2=Rasch->Nu[Rasch->Params[i].Node2-1];
      else
       tmpNode2=Rasch->Params[i].Node2;
//      tmpStr=AnsiString(Rasch->Nu[Rasch->Params[i].Node1-1])+","+AnsiString(Rasch->Nu[Rasch->Params[i].Node2-1])+" "+FloatToStrF(real(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].DSTC),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTC),ffFixed,10,4);
      tmpStr=AnsiString(tmpNode1)+","+AnsiString(tmpNode2)+" "+FloatToStrF(real(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTM),ffFixed,10,4)+'\t'+FloatToStrF(real(Rasch->Params[i].DSTC),ffFixed,10,4)+'\t'+FloatToStrF(imag(Rasch->Params[i].DSTC),ffFixed,10,4);
      }
    Vals[CurrIt].DSTM[i]=Rasch->Params[i].DSTM;
    Vals[CurrIt].DSTC[i]=Rasch->Params[i].DSTC;
    Vals[CurrIt].IsTrf[i]=true;
    TmpList->Add(tmpStr);
   }

tmpStr=FileName.SubString(1,FileName.Length()-4)+"res.txt";

if(CurrIt!=0&&Showtxt==true)
 {
  TmpList2->LoadFromFile(tmpStr);
 }
else
 {
  TmpList2->Clear();
 }
for(i=0;i<TmpList->Count;i++)
 TmpList2->Add(TmpList->Strings[i]);
if(Showtxt==true)
 TmpList2->SaveToFile(tmpStr);
RezFileName="notepad.exe "+tmpStr;
//WinExec(tmpStr.c_str(),SW_RESTORE);

}
__finally
{
delete TmpList;
delete TmpList2;
}

}



ParamForCycle::ParamForCycle()
{
Ok=false;
NumberOfIters=0;
NumberOfParams=0;
FileName="";
ElementNo=-1;
}
//---------------------------------------------------------------------------
ParamForCycle::~ParamForCycle()
{
Quantity--;
}
//---------------------------------------------------------------------------
void ParamForCycle::Add(int ElNo,AnsiString FName)
{
FileName=FName;
if(ReadFromFile()==true)
  {
   Ok=true;
   Quantity++;
   ElementNo=ElNo;
  }

else
  {
   Ok=false;
   FileName="";
   ShowMessage("�� ���� ��������� ����!");
  }

}
//---------------------------------------------------------------------------
bool ParamForCycle::ReadFromFile()
{
TStringList *tmpList=new TStringList();
AnsiString tmpStr,tmpStr2;
int i,j,k,n,max,maxk[MaxElementParams];
double tmpDouble;
try
{
tmpList->LoadFromFile(FileName);
n=tmpList->Count;
max=0;
for(i=0;i<n;i++)
  {
   tmpStr=tmpList->Strings[i];
   tmpStr2="";
   k=-1;
   for(j=1;j<=tmpStr.Length();j++)
     {
      if(tmpStr[j]==','||tmpStr[j]=='.')
        {
         tmpStr[j]=',';
         tmpStr2=tmpStr2+tmpStr[j];
        }
      else if(tmpStr[j]==';')
        {
         tmpDouble=StrToFloat(tmpStr2);
         tmpStr2="";
         k++;
         if(k>max)
            max=k;
         maxk[i]=k;
         Params[i][k]=tmpDouble;
        }
      else
        {
         tmpStr2=tmpStr2+tmpStr[j];
        }
     }

  }
for(i=0;i<n;i++)
 if(maxk[i]<max)
   for(j=maxk[i];j<max;j++)
     Params[i][j+1]=Params[i][j];
NumberOfIters=max+1;
NumberOfParams=n;
Ok=true;
}
catch(...)
{
Ok=false;
}
return Ok;
}
//---------------------------------------------------------------------------
void TForm1::CalculateMain(AnsiString FName)
{
AnsiString tmpStr;
List=new TStringList;
bool Finish=false;
int i,j,k,l;
i=0;j=0,k=0,l=0;
int FindObject[MaxNodesInElem],FindPoint[MaxNodesInElem],FindPointType[2];
NumberOfLines=0;int CurrLine=0;
//AnsiString FName;
for (i=0;i<NumberOfElem;i++)
  {
   if(Element[i]->NumberOfPoints==3&&Element[i]->Name!="����")
     NumberOfLines+=3;
   else if(Element[i]->NumberOfPoints==2)
     NumberOfLines+=1;
  }
delete []StruList;
StruList=new struct Stru[NumberOfLines];
i=0;
try
   {
      while(i<NumberOfElem)
      {
       tmpStr="";
       if(Element[i]->NumberOfPoints==2)
          for (j=0;j<Element[i]->NumberOfPoints;j++)
           {
            if(Element[i]->PointNumber[j]!=-1)
              {
               tmpStr+=AnsiString(Element[i]->PointNumber[j])+", ";
               if(j==0)
                  StruList[CurrLine].i=Element[i]->PointNumber[j];
               else
                  StruList[CurrLine].j=Element[i]->PointNumber[j];

               FindPointType[j]=PointType(i,Element[i]->PointNumber[j]);
              }
            else
              goto cyc2end;
           }
      else if(Element[i]->NumberOfPoints==3&&Element[i]->ElementName!="����")
         {
          k++;
          if(Element[i]->PointNumber[j]!=-1)
           {
            switch(k)
            {
            case 1:
               tmpStr+=AnsiString(Element[i]->PointNumber[0])+", ";
               tmpStr+=AnsiString(Element[i]->PointNumber[1])+", ";
               StruList[CurrLine].i=Element[i]->PointNumber[0];
               StruList[CurrLine].j=Element[i]->PointNumber[1];
               break;
            case 2:
               tmpStr+=AnsiString(Element[i]->PointNumber[1])+", ";
               tmpStr+=AnsiString(Element[i]->PointNumber[2])+", ";
               StruList[CurrLine].i=Element[i]->PointNumber[1];
               StruList[CurrLine].j=Element[i]->PointNumber[2];

               break;
            case 3:
               tmpStr+=AnsiString(Element[i]->PointNumber[0])+", ";
               tmpStr+=AnsiString(Element[i]->PointNumber[2])+", ";
               StruList[CurrLine].i=Element[i]->PointNumber[0];
               StruList[CurrLine].j=Element[i]->PointNumber[2];

//               k=0;
               break;
            }
            for(l=0;l<2;l++)
            {
            if(k==3&&l==1)
              j=2;
            else if(k==3&&l==0)
              j=0;
            else
              j=l+k-1;
            FindPointType[l]=PointType(i,Element[i]->PointNumber[j]);

            }//while
           }
         }
       else
         {
          goto cyc2end;
         }
       tmpStr+=AnsiString(Element[i]->type)+", ";
       StruList[CurrLine].type=Element[i]->type;


       if(k!=0)
         {
          tmpStr+=AnsiString(k)+", ";
          StruList[CurrLine].u4=k;
         }
       else
         StruList[CurrLine].u4=0;
       if(k==3)
         k=0;


       tmpStr+=AnsiString(FindPointType[0])+", ";
       tmpStr+=AnsiString(FindPointType[1]);
       StruList[CurrLine].u5=FindPointType[0];
       StruList[CurrLine].u6=FindPointType[1];
//       tmpStr+='\n';
       List->Add(tmpStr);
       CurrLine++;
       cyc2end:
       if(k==0)
          i++;
      }
NumberOfLines=CurrLine;
List->Add("");
List->Add("");


AnsiString tmpNumber;
for(i=0;i<NumberOfElem;i++)
   if(Element[i]->ElementName!="����")
   {
    tmpStr="";
    tmpStr+=AnsiString(Element[i]->type)+", ";
    for(j=0;j<Element[i]->NumberOfPoints;j++)
      if(Element[i]->PointNumber[j]!=-1)
        tmpStr+=AnsiString(Element[i]->PointNumber[j])+", ";
      else
        goto cyc3end;
    for(j=0;j<Element[i]->NumberOfParams-1;j++)
      {
       tmpNumber=FloatToStr(Element[i]->Params[j]);
       for(k=1;k<=tmpNumber.Length();k++)
         if(tmpNumber[k]==',')
           tmpNumber[k]='.';
       tmpStr+=tmpNumber+", ";
      }
    tmpNumber=FloatToStr(Element[i]->Params[Element[i]->NumberOfParams-1]);
    for(k=1;k<=tmpNumber.Length();k++)
      if(tmpNumber[k]==',')
        tmpNumber[k]='.';
    tmpStr+=tmpNumber;

//    tmpStr+='\n';
    List->Add(tmpStr);
    cyc3end:
   }
if(FName!=NULL)
  List->SaveToFile(FName);
//    List->Add(
   }
__finally
   {
    delete List;
   }
}
//---------------------------------------------------------------------------

__fastcall MyForm::MyForm(TComponent* Owner)
        : TForm(Owner)
{

}
//---------------------------------------------------------------------------


void __fastcall TForm1::N19Click(TObject *Sender)
{
Form5->IsMsw=false;
Form5->Show();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N20Click(TObject *Sender)
{
Form5->IsMsw=true;
Form5->Show();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::N21Click(TObject *Sender)
{
Form3->Show();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N22Click(TObject *Sender)
{
MenuExit();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N24Click(TObject *Sender)
{
AnsiString prec;

prec=InputBox("������� �������� �������","������� �������� ������� �����:",AnsiString(precision));
for(int i=1;i<=prec.Length();i++)
  if(prec[i]=='.')
    prec[i]=',';
try
{
precision=StrToFloat(prec);
}
catch(...)
{
ShowMessage("������������ ������ �����!");
}

}
//---------------------------------------------------------------------------


void __fastcall TForm1::N25Click(TObject *Sender)
{
if(N25->Checked==false)
  {
   N25->Checked=true;
   IsTxtsOut=true;
  }
else
  {
   N25->Checked=false;
   IsTxtsOut=false;
  }
}
//---------------------------------------------------------------------------

bool TForm1::UstRegimRaschet(int type,bool BkUse,bool CurrReg,bool Showtxt,double *Vr0,int nVr0, bool WritePrmsToFile,AnsiString PrmsFileName)
{
bool continueJ;
bool Regim;
CurrIt=0;
Rasch->mhu_n=new int[Rasch->n];
//Values::NodesNumber=Rasch->n;
//Values::mswNumber=Rasch->n*Rasch->n;
delete []Vals;
Vals=new Values[1];
//Vals[0].Initialize(Rasch->n,Rasch->n*Rasch->n);
Vals[0].Initialize(Rasch->n,Rasch->nl);
for(int i=0;i<Rasch->n;i++)
 Rasch->mhu_n[i]=Rasch->mhu[i];
double normtmp=0;
double uMax=0;
Rasch->BkUse=BkUse;
Rasch->CurrReg=CurrReg;
if(type==7||type==8||type==5)
 {
  Rasch->NormU=0;

  for(int i=0;i<Rasch->n;i++)
   {
    if(Rasch->mhu[i]!=0&&Rasch->mhu[i]!=4)
     Rasch->NormU+=abs(Rasch->unom[i])*abs(Rasch->unom[i]);
   }
//  Rasch->NormU=sqrt(Rasch->NormU)/(Rasch->n-Rasch->nb-Rasch->ng);

  Rasch->NormU=1;
  if(CurrReg==false)
   {
   for(int i=0;i<nVr0;i++)
    {
     Rasch->R[i]=Vr0[i];
     normtmp+=Vr0[i]*Vr0[i];
    }
   normtmp=sqrt(normtmp);
   }
  else
   {
    normtmp=1;
   }
  if(normtmp!=0)
  for(int i=0;i<nVr0;i++)
    Rasch->R[i]=Rasch->R[i]*Rasch->NormU/normtmp;
  normtmp=0;
  for(int i=0;i<nVr0;i++)
    normtmp+=Rasch->R[i]*Rasch->R[i];
  normtmp=sqrt(normtmp);
 }
Rasch->T=T0;
Rasch->Kn=Kn0;
if(type==0)
 Regim=Raschet(Showtxt);
else if(type==1)
 Regim=RaschetJ(Showtxt,continueJ);
else if(type==2)
 Regim=RaschetIskl(Showtxt);
else if(type==3)
 Regim=RaschetJEkviv(Showtxt,continueJ);
else if(type==4)
 Regim=RaschetJX2(Showtxt,continueJ);
//else if(type==5)
// Regim=RaschetDrag(Showtxt,continueJ);//���������� ����� - ����� ���������, ��
else if(type==6)
 {
  Regim=RaschetJ(false,continueJ);
  if(Regim==true)
   {
   Rasch->EkvivalentToStar();
   Regim=RaschetJ(Showtxt,continueJ);
   }
 }
else if(type==7)
 Regim=RaschetJX22(Showtxt,continueJ,WritePrmsToFile,PrmsFileName);//���������� ����� - ����� �������, ��
else if(type==8)
 Regim=RaschetJX2T(Showtxt,continueJ,WritePrmsToFile,PrmsFileName);//���������� ����� - ����� �������, �
else if(type==9)
 Regim=RaschetJSimpleT(Showtxt,continueJ);//���������� ����� - ���������� ����������, �
else if(type==10)
 {
  Regim=RaschetJ(false,continueJ);

 }
if(Showtxt==true)
 WinExec(RezFileName.c_str(),SW_RESTORE);
if(Regim==true)
 CopyResultsParams(0,0);
FormR->Iters=Iters;
//Vals=NULL;
//delete []Vals;
if(type!=6)
for(int i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];
delete []Rasch->mhu_n;
return Regim;
}

//---------------------------------------------------------------------------
void TForm1::RaschetVr0(bool type,int &nVr,double *Vr)
{
/*
 bool Regim,flag=false;
 Rasch->YsCreate();
 Rasch->NjNbNg();
 Rasch->InitializeJac(Rasch->nJac);
 Rasch->Kn=Kn0;
 Rasch->T=T0;
 Rasch->mhu_n=new int [Rasch->n];
 for(int i=0;i<Rasch->n;i++)
  {
   Rasch->urab[i]=Rasch->unom[i];
   Rasch->mhu_n[i]=Rasch->mhu[i];
  }
 CurrIt=0;
// Values::NodesNumber=Rasch->n;
// Values::mswNumber=Rasch->n*Rasch->n;
 delete []Vals;
 Vals=new Values[1];
 Vals[0].Initialize(Rasch->n,Rasch->n*Rasch->n);
 if(type==0)
  Regim=RaschetJKn(Kn0,false);
 else
  Regim=RaschetJT(T0,false);
 if(Regim==false)
 {
  for(int i=0;i<Rasch->n;i++)
   Rasch->urab[i]=Rasch->unom[i];
  Regim=RaschetJ(false,false);
  flag=true;
 }
 for(int i=0;i<Rasch->n;i++)
   Rasch->mhu[i]=Rasch->mhu_n[i];
 delete []Rasch->mhu_n;
 if(Regim==false)
  {
   for(int i=0;i<Rasch->n;i++)
    Rasch->urab[i]=Rasch->unom[i];
  }
 if(Regim==true&&type==0)
   Rasch->Jacobi_KnCalc();
 else
   Rasch->Jacobi2Calc();
 //Rasch->RaschR0(Rasch->Jac,Rasch->nJac,Vr);
 double NormForSingular=0;
if(Regim==true)
{
 for(int i=0;i<Rasch->n;i++)
 {
  if(Rasch->mhu[i]!=0&&Rasch->mhu[i]!=4)
    NormForSingular+=abs(Rasch->urab[i])*abs(Rasch->urab[i]);
  if(Rasch->mhu[i]!=0)
    NormForSingular+=arg(Rasch->urab[i])*arg(Rasch->urab[i]);
 }
 NormForSingular=sqrt(NormForSingular);
}
 Rasch->singular(Rasch->Jac,Rasch->nJac,Vr);
 //Rasch->RaschR0(Rasch->Jac,Rasch->nJac,Vr);


 nVr=Rasch->nJac;
 */
}
//---------------------------------------------------------------------------
void TForm1::RaschetTranspVr(bool type,bool MatrixType,int &nVr,double *Vr)
{
ap::real_2d_array tmpm;
if(Rasch->RegimState==true&&type==0)
  Rasch->Jacobi_KnCalc();
else if(Rasch->RegimState==true)
  Rasch->Jacobi2Calc();
nVr=Rasch->nJac;
int k;
if(MatrixType==1)
 {
  Rasch->NjNbNg();
  k=0;
  for(int i=0;i<Rasch->n;i++)
  {
   if(Rasch->mhu[i]==0||Rasch->mhu[i]==4)
     Rasch->Zu[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
   else
    {
     Rasch->Zu[i]=k;
     k++;
    }
  }
  for(int i=0;i<Rasch->n;i++)
  {
   if(Rasch->mhu[i]==0)
     Rasch->Zf[i]=-1;//���� � �������� ����� ������ �� ��� ��������, ������ �������� �������!!!
   else
    {
     Rasch->Zf[i]=k;
     k++;
    }
  }
  Rasch->nJac*=2;
  if(Rasch->RegimState==true&&type==0)
    Rasch->JacobiX22Calc();
  else if(Rasch->RegimState==true)
    Rasch->JacobiX2TCalc();
//  Rasch->JacobiShow(Form1->Memo1);
  nVr=Rasch->nJac;
  Rasch->NjNbNg();
}
tmpm.setbounds(1,nVr,1,nVr);
for(int i=0;i<nVr;i++)
 for(int j=0;j<nVr;j++)
  tmpm(i+1,j+1)=Rasch->Jac[i][j];
double determ=determinant(tmpm,nVr);
//Memo1->Lines->Add("Determinant");
//Memo1->Lines->Add(determ);
//Memo1->Lines->Add("***********");
double **transpJac=new double *[nVr];
for(int i=0;i<nVr;i++)
 {
  transpJac[i]=new double [nVr];
 }
for(int i=0;i<nVr;i++)
  for(int j=0;j<nVr;j++)
//     transpJac[i][j]=Rasch->Jac[i][j];
   transpJac[i][j]=Rasch->Jac[j][i];

 Rasch->RaschR0(transpJac,nVr,Vr);
for(int i=0;i<nVr;i++)
 delete transpJac[i];
delete []transpJac;

}
//---------------------------------------------------------------------------

void TForm1::RaschetSingular(bool type,int &nVr,double *Vr)
{
/*
 bool Regim,flag=false;
 Rasch->YsCreate();
 Rasch->NjNbNg();
 Rasch->InitializeJac(Rasch->nJac);
 Rasch->Kn=Kn0;
 Rasch->T=T0;
 Rasch->mhu_n=new int [Rasch->n];
 for(int i=0;i<Rasch->n;i++)
  {
   Rasch->urab[i]=Rasch->unom[i];
   Rasch->mhu_n[i]=Rasch->mhu[i];
  }
 CurrIt=0;
// Values::NodesNumber=Rasch->n;
// Values::mswNumber=Rasch->n*Rasch->n;
 delete []Vals;
 Vals=new Values[1];
 Vals[0].Initialize(Rasch->n,Rasch->n*Rasch->n);
 if(type==0)
  Regim=RaschetJKn(Kn0,false);
 else
  Regim=RaschetJT(T0,false);
 if(Regim==false)
 {
  for(int i=0;i<Rasch->n;i++)
   Rasch->urab[i]=Rasch->unom[i];
  Regim=RaschetJ(false,false);
  flag=true;
 }
 for(int i=0;i<Rasch->n;i++)
   Rasch->mhu[i]=Rasch->mhu_n[i];
 delete []Rasch->mhu_n;
 if(Regim==false)
  {
   for(int i=0;i<Rasch->n;i++)
    Rasch->urab[i]=Rasch->unom[i];
  }
 if(Regim==true&&type==0)
   Rasch->Jacobi_KnCalc();
 else
   Rasch->Jacobi2Calc();
 Rasch->singular(Rasch->Jac,Rasch->nJac,Vr);
*/
}

//---------------------------------------------------------------------------
void TForm1::HardRaschet(int type)
{
int i,j,k,state;
float d,T,T0,Tpred;
bool flag;
bool ok;
bool continueJ;
if(MainElementNo==-1)
  {
   for(i=0;i<NumberOfElem;i++)
     if(Element[i]->NumberOfPoints==1)
       {
        MainElementNo=i;
        break;
       }
   Numerate();
  }
CalculateMain(NULL);
Prepare();
switch (type)
{
case 0:
{Raschet(false);
break;}
case 1:
{RaschetJ(false,continueJ);
break;}
case 2:
{RaschetIskl(false);
break;}
case 3:
{RaschetJEkviv(false,continueJ);
break;}
case 4:
{RaschetJX2(false,continueJ);
break;}
//case 6:
//{RaschetJKn(false,continueJ,true,false,true);
//break;}
case 7:
{RaschetJX22(false,continueJ,false,"");
break;}
}
float maxdPQ=0;
int maxdPQNumber=-1;
int maxIs=-1;
for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dPg[i])>0&&fabs(Rasch->Pg[i])>maxdPQ)
  {maxdPQ=fabs(Rasch->Pg[i]);
   maxdPQNumber=i;
   maxIs=0;
  }
for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dPn[i])>0&&fabs(real(Rasch->sng[i]))>maxdPQ)
  {maxdPQ=fabs(real(Rasch->sng[i]));
   maxdPQNumber=i;
   maxIs=1;
  }

for(i=0;i<Rasch->n;i++)
  if(fabs(Rasch->dQn[i])>0&&fabs(imag(Rasch->sng[i]))>maxdPQ)
  {maxdPQ=fabs(imag(Rasch->sng[i]));
   maxdPQNumber=i;
   maxIs=2;
  }
if(maxdPQNumber!=-1)
{
if(maxIs==0)
 T0=fabs(Rasch->Pg[maxdPQNumber])/fabs(Rasch->dPg[maxdPQNumber]);
else if(maxIs==1)
 T0=fabs(real(Rasch->sng[maxdPQNumber]))/fabs(Rasch->dPn[maxdPQNumber]);
else
 T0=fabs(imag(Rasch->sng[maxdPQNumber]))/fabs(Rasch->dQn[maxdPQNumber]);
}

double *dX=new double [Rasch->n*2];
Rasch->YsCreate();
Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
Rasch->mhu_n=new int [Rasch->n];
ok=false;
d=2.;
T=T0;
flag=true;
float CurrPrec=0.01;
int kit=0;
//T=0;
while(ok==false&&maxdPQNumber!=-1)
{
if(flag==false)
  for(i=0;i<Rasch->n;i++)
   Rasch->urab[i]=Rasch->unom[i];
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
if(state==0)//������� �����
 {
  T+=T0/d;
  if(flag==false)
    d*=2;
  if(fabs(Tpred-T)<CurrPrec)
    ok=true;
  Tpred=T;
  flag=true;
 }
else
 {
  T-=T0/d;
  if(flag==true)
    d*=2;
  Tpred=T;
  flag=false;
 }
}

Values::NodesNumber=Rasch->n;
Values::mswNumber=Rasch->n*Rasch->n;
Vals=new Values[int(T)+1];
CurrIt=0;

for(j=0;j<int(T);j++)
{
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,double(j));
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
Rasch->ResultsCalc();
ResultsOut(kit,state,true);
CurrIt++;
}
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
Rasch->ResultsCalc();
//if(showresult==true)
ResultsOut(kit,true,true);
CurrIt=0;
HardItersNumber=T;
delete []dX;
delete []Rasch->mhu_n;
WinExec(RezFileName.c_str(),SW_RESTORE);
ShowMessage(T);
}

//---------------------------------------------------------------------------
bool TForm1::RaschetSeries(int type,double T,double Kn,int &n,bool Showtxt)
{
int i;
bool Regim,continueJ;
bool UseUnom=true;
Rasch->mhu_n=new int [Rasch->n];
for(i=0;i<Rasch->n;i++)
 Rasch->mhu_n[i]=Rasch->mhu[i];

//Values::NodesNumber=Rasch->n;
//Values::mswNumber=Rasch->n*Rasch->n;
delete []Vals;
Vals=new Values[n];
for(i=0;i<n;i++)
 Vals[i].Initialize(Rasch->n,Rasch->nl);
// Vals[i].Initialize(Rasch->n,Rasch->n*Rasch->n);
double dP;
if(type==0)
 dP=T/double(n-1);
else
 {dP=fabs(1.-Kn)/double(n-1);dP-=dP*.000001;}

Rasch->YsCreate();
CurrIt=0;
if(type==0)
for(i=0;i<n;i++)
{
 Regim=RaschetJT(dP*double(i),UseUnom);
 if(Regim==true)
  UseUnom=false;
 else
  {n=i;break;}
 CurrIt++;
}
else
for(i=0;i<n;i++)
{
 Regim=RaschetJKn(1.-dP*double(i),UseUnom);
 if(Regim==true)
  UseUnom=false;
 else
  {n=i;break;}
 CurrIt++;
}
//if(Regim==true)
//if(n>0)
// FormR->Vals=Vals;
for(i=0;i<Rasch->n;i++)
 Rasch->mhu[i]=Rasch->mhu_n[i];
delete []Rasch->mhu_n;
//Vals=NULL;
//delete []Vals;
return Regim;
}


//---------------------------------------------------------------------------
bool TForm1::RaschetJT(double T,bool UseUnom)
{

int i,j,kit,state;
bool Regim,limits;
double *dX=new double[Rasch->n*2];
if(UseUnom==true)
for(i=0;i<Rasch->n;i++)
 Rasch->urab[i]=Rasch->unom[i];

Rasch->NjNbNg();
Rasch->InitializeJac(Rasch->nJac);
kit=0;
for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
if(UseQminQmax==true)
if(Rasch->Limits(0,true,T)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}
if(UseQminQmax==true)
if(Rasch->Limits(1,true,T)==true)
 {
  limits=true;
  if(Regim==false)
   {
    for(i=0;i<Rasch->n;i++)
     Rasch->urab[i]=Rasch->unom[i];
    kit=0;
   }
 Rasch->NjNbNg();
 for(i=0;i<NumberOfIters;i++)
  {
   Rasch->SHN_Hard(true,T);
   Rasch->Jacobi2Calc();
   Rasch->QR(Rasch->Jac,Rasch->Nebal,dX,Rasch->nJac,1);
   state=!Rasch->RaschCheck(dX,Rasch->urab,Rasch->nJac,-1,0,LimitNebals,dUmax,dfmax,dKnmax,dTmax,dVrimax);
   if(state==0||state==2)
     break;
   state=4;
   kit++;
  }
}

Rasch->ResultsCalc();
ResultsOut(kit,state,false);
delete dX;
if (state==0)
 Regim=true;
Rasch->RegimState=Regim;
return Regim;
}

//---------------------------------------------------------------------------







void __fastcall TForm1::N32Click(TObject *Sender)
{
//����� ������ ����� ����
int tmpN;
tmpN=ImageForRightBtn->RealPointNumber[ImageForRightBtn->MainRectNo];
ImageForRightBtn->RealPointNumber[ImageForRightBtn->MainRectNo]=StrToInt(
InputBox("������� ����� ����","������� ����� �����:",AnsiString(tmpN)));
Redraw();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::N35Click(TObject *Sender)
{
RealNodes=false;
N35->Checked=true;
Redraw();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::N34Click(TObject *Sender)
{
RealNodes=true;
N34->Checked=true;
Redraw();
}
//---------------------------------------------------------------------------
void TForm1::FillRealNodesNumbers()
{
for(int k=0;k<Rasch->n;k++)
Rasch->Nu[k]=-1;
for(int k=0;k<Rasch->n;k++)
  for (int i=0;i<NumberOfElem;i++)
    for(int j=0;j<Element[i]->NumberOfPoints;j++)
      if(Element[i]->RealPointNumber[j]!=-1&&Element[i]->PointNumber[j]==k+1)
        Rasch->Nu[k]=Element[i]->RealPointNumber[j];

}
//---------------------------------------------------------------------------









int TForm1::MyRound(double val)
{
double tmpval;
tmpval=val-double(int(val));
if(tmpval>=0.5)
 return int(val)+1;
else
 return int(val);
}

//---------------------------------------------------------------------------
void __fastcall TForm1::FormClose(TObject *Sender, TCloseAction &Action)
{
Action=caNone;
MenuCloseSchema();
//this->Hide();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::N39Click(TObject *Sender)
{

MenuNew();


}
//---------------------------------------------------------------------------
void TForm1::CleanSchema()
{
int i;


for(i=0;i<NumberOfElem;i++)
  {
   delete Element[i];
   Element[i]=NULL;
  }

SelectedObject=NULL;
NumberOfElem=0;
MouseBtnLeftIsDown=false;
ObjectIsClicked=false;
MainElementNo=-1;
BalanceElementNo=-1;
MaxNumberOfIt=1;

NumberOfElem=0;
if(LineCoords::LineNumber>0)
{
  delete []LineArray;
  LineArray=NULL;
}
LineCoords::LineNumber=0;
LinesNumber=0;
Redraw();

}
//---------------------------------------------------------------------------

void TForm1::MenuNew()
{
int MsgDlgResult;
if(HasChanged==true)
{
MsgDlgResult=MessageDlg("��������� ���������?", mtConfirmation, TMsgDlgButtons() << mbYes << mbNo << mbCancel, 0);
if (MsgDlgResult == mrYes)
   N13Click(N13);
else if(MsgDlgResult == mrCancel)
goto NothingToDo;
}
CleanSchema();
HasChanged=false;
NewFile=true;
FileName="";
Form1->Caption="����� �����";
if(Form1->Visible==false)
  Form1->Show();
NothingToDo:
}

//---------------------------------------------------------------------------
void TForm1::MenuOpen()
{
int MsgDlgResult;
if(NewFile==true||HasChanged==true)
{
MsgDlgResult=MessageDlg("��������� ���������?", mtConfirmation, TMsgDlgButtons() << mbYes << mbNo << mbCancel, 0);
if (MsgDlgResult == mrYes)
   MenuSave();
else if(MsgDlgResult == mrCancel)
goto NothingToDo;
}
if(OpenDialog1->Execute())
  {
    Open(OpenDialog1->FileName);
    NewFile=false;
    HasChanged=false;
    FileName=OpenDialog1->FileName;
    Form1->Caption="����� - "+OpenDialog1->FileName;
    if(Form1->Visible==false)
      Form1->Show();

   }
NothingToDo:
}


//---------------------------------------------------------------------------
void TForm1::MenuSave()
{


SaveDialog1->Filter="������������� ����� (*.els)| *.els";
SaveDialog1->DefaultExt="els";
if(NewFile==true&&HasChanged==true)
 {
  SaveDialog1->Execute();
  if(SaveDialog1->FileName!="")
     {

      Save(SaveDialog1->FileName);
      FileName=SaveDialog1->FileName;
      NewFile=false;
      HasChanged=false;
      Form1->Caption="����� - "+SaveDialog1->FileName;
      SaveDialog1->FileName="";

     }

 }
else if (HasChanged==true)
 {
  Save(FileName);
  Form1->Caption="����� - "+FileName;
  SaveDialog1->FileName="";
  NewFile=false;
  HasChanged=false;
 }

}


//---------------------------------------------------------------------------
void TForm1::MenuSaveAs()
{
  SaveDialog1->Execute();
  if(SaveDialog1->FileName!="")
     {
      Save(SaveDialog1->FileName);
      FileName=SaveDialog1->FileName;
      NewFile=false;
      Form1->Caption="����� - "+SaveDialog1->FileName;
      SaveDialog1->FileName="";
     }
}
//---------------------------------------------------------------------------
bool TForm1::MenuCloseSchema()
{
int MsgDlgResult;
bool state=false;
if(HasChanged==true)
{
MsgDlgResult=MessageDlg("��������� ���������?", mtConfirmation, TMsgDlgButtons() << mbYes << mbNo << mbCancel, 0);
if (MsgDlgResult == mrYes)
   N13Click(N13);
else if(MsgDlgResult == mrCancel)
  goto NothingToDo;
}
CleanSchema();
HasChanged=false;
FileName="";
Form1->Caption="����� �����";
Form1->Hide();
state=true; //��������� ��� ��������
NothingToDo:
return state;
}

void __fastcall TForm1::N40Click(TObject *Sender)
{
MenuCloseSchema();
}
//---------------------------------------------------------------------------
void TForm1::MenuExit()
{
if(MenuCloseSchema()==true)
 Application->Terminate();

}


void TForm1::CopyRasch()
{
RRasch=Rasch;
}

//---------------------------------------------------------------------------
int TForm1::CreateRasch()
{
CalculateMain(NULL);
Prepare();
FillRealNodesNumbers();
Rasch->YsCreate();
return Rasch->n;
}
//---------------------------------------------------------------------------
int TForm1::CreateRasch63Proba()
{
/*
RaschetJProba(false,false);
Rasch->YsCreate();
return Rasch->n;
*/
}
//---------------------------------------------------------------------------

void TForm1::SaveRasch()
{
//delete Rasch;
//Rasch=RRasch;
}
//---------------------------------------------------------------------------
void TForm1::CopyResultsParams(int sourceind, int destind)
{
//FormR->Vals[destind]=Vals[sourceind];
FormR->T=Rasch->T;
FormR->Kn=Rasch->Kn;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::N10Click(TObject *Sender)
{
Form5->Show();
}
//---------------------------------------------------------------------------
void TForm1::IntervalRaschet()
{
// ������ ������������ �������
int i,j,k;
int *ki=new int [Rasch->n];
int *kj=new int [Rasch->n];
int ii,jj;
complex<double> tmpCompl,tmpUYs,tmpU,tmpUYb;
complex<double> *tngr=new complex<double>[Rasch->n];

int krm;//���������� �����, �� ������ ��������� � ��
int kit;//���-�� ��������
double p,q;
double osh;//������
krm=0;
AnsiString tmpStr,tmpStr2;
bool Regim=false;
bool RegOk;

//double obuslovl;
//int skip=-1;


//��������� ������� ������������� Ys
krm=Rasch->YsCreate();
ShowMessage(Rasch->ys[0][0]);

/*
//�������� �������� ������� - z[][]
Rasch->ZCreate(ki,kj,ii,jj);



//������ ������� ����������
for(i=0;i<Rasch->n;i++)
{
 Rasch->urab[i]=Rasch->unom[i];
 if(abs(Rasch->unom[i])!=0)
   {
    if(Rasch->mhu[i]==4)
      Rasch->sng_r[i]=complex<double>(-Rasch->Pg[i],0)+Rasch->sng[i];
    else
      Rasch->sng_r[i]=Rasch->sng[i];
    tngr[i]=Rasch->sng_r[i]/conj(Rasch->unom[i]);
    tngr[i];
   }

}
kit=0;
goto label1;
label2:
for(i=0;i<Rasch->n;i++)
{
 if(Rasch->mhu[i]==1)
   {
    p=real(Rasch->sng_r[i])*(Rasch->a0[i]+(Rasch->a1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->a2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
    q=imag(Rasch->sng_r[i])*(Rasch->b0[i]+(Rasch->b1[i]*abs(Rasch->urab[i])/abs(Rasch->unom[i]))+(Rasch->b2[i]*abs(Rasch->urab[i])*abs(Rasch->urab[i])/(abs(Rasch->unom[i])*abs(Rasch->unom[i]))));
    tngr[i]=complex<double>(p,q)/conj(Rasch->urab[i]);
//    tngr[i];
   }
}
label1:
for(ii=0;ii<krm;ii++)
{
 tmpUYs=0;
 tmpUYb=0;
// Rasch->JacobiCalc();
 for(jj=0;jj<krm;jj++)
   tmpUYs+=Rasch->ys[ki[ii]][kj[jj]]*Rasch->urab[kj[jj]];


 for(j=0;j<Rasch->n;j++)
   if(Rasch->mhu[j]==0||Rasch->mhu[j]==3)
     tmpUYb+=Rasch->ys[j][ki[ii]]*Rasch->urab[j];


 Rasch->irab[ki[ii]]=tmpUYs-tngr[ki[ii]]+tmpUYb;
 Rasch->irab[ki[ii]];
}
osh=0;
RegOk=true;
for(ii=0;ii<krm;ii++)
{
 tmpU=complex<double>(0,0);
 for(jj=0;jj<krm;jj++)
   {
    tmpU+=Rasch->z[ii][jj]*Rasch->irab[kj[jj]];
   }
 Rasch->urab[ki[ii]]-=tmpU;
// osh+=abs(tmpU);
 osh=abs(tmpU);
 if(osh>=Rasch->eps)
   RegOk=false;
}
kit++;
if(kit>NumberOfIters)
  goto Raskhod;//����� ����������
osh=osh/(Rasch->n-1);
//if(osh>Rasch->eps)
if(RegOk==false)
  goto label2;//�� ����� �� ��������� �������� - ���������� ������
else
  {
   Regim=true;
   goto Ok;//����� ���������!
  }

Raskhod:
ShowMessage("����� �� �������� - ���������� �����������!");
Regim=false;

Ok:
//����� ������� ��������� ������������ ��������� ������
tmpU=complex<double>(0,0);
for(i=0;i<Rasch->n;i++)
  tmpU+=Rasch->ys[0][i]*Rasch->urab[i];
tmpU=tmpU*conj(Rasch->urab[0]);
Rasch->ResultsCalc();
//� �������, ���� ����
//ResultsOut(kit,Regim,showresult);      <-------------
//delete []Vals;
delete []ki;
delete []kj;
delete []tngr;
//return Regim;
*/
}
//---------------------------------------------------------------------------
