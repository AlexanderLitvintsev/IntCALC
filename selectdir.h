//---------------------------------------------------------------------------

#ifndef selectdirH
#define selectdirH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <FileCtrl.hpp>
//---------------------------------------------------------------------------
class TFormSelectDir : public TForm
{
__published:	// IDE-managed Components
        TDirectoryListBox *DirectoryListBox1;
        TDriveComboBox *DriveComboBox1;
        TButton *Button1;
        TButton *Button2;
        TFileListBox *FileListBox1;
        void __fastcall DriveComboBox1Change(TObject *Sender);
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall DirectoryListBox1Change(TObject *Sender);
private:	// User declarations
public:		// User declarations
        AnsiString DirName;
        AnsiString FileName;
        __fastcall TFormSelectDir(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFormSelectDir *FormSelectDir;
//---------------------------------------------------------------------------
#endif
