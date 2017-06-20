//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "selectdir.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFormSelectDir *FormSelectDir;
//---------------------------------------------------------------------------
__fastcall TFormSelectDir::TFormSelectDir(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TFormSelectDir::DriveComboBox1Change(TObject *Sender)
{
DirectoryListBox1->Drive=DriveComboBox1->Drive;
}
//---------------------------------------------------------------------------
void __fastcall TFormSelectDir::Button1Click(TObject *Sender)
{
DirName=DirectoryListBox1->Directory;
FileName=FileListBox1->FileName;
ModalResult = mrOk;

}
//---------------------------------------------------------------------------
void __fastcall TFormSelectDir::Button2Click(TObject *Sender)
{
ModalResult = mrCancel;        
}
//---------------------------------------------------------------------------
void __fastcall TFormSelectDir::DirectoryListBox1Change(TObject *Sender)
{
FileListBox1->Directory=DirectoryListBox1->Directory;
}
//---------------------------------------------------------------------------
