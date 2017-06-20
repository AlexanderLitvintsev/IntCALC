//---------------------------------------------------------------------------

#ifndef ToolsH
#define ToolsH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ToolWin.hpp>
#include <ImgList.hpp>
//---------------------------------------------------------------------------
class TForm3 : public TForm
{
__published:	// IDE-managed Components
        TImageList *ImageList1;
        TCoolBar *CoolBar1;
        TToolBar *ToolBar1;
        TToolBar *ToolBar2;
        TPageScroller *PageScroller1;
        TImageList *sysicons;
        TToolButton *ToolButton1;
        TComboBox *ComboBox1;
        TToolButton *ToolButton2;
        TToolButton *ToolButton3;
        TToolButton *ToolButton4;
	TToolButton *ToolButton5;
        void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
        void __fastcall ComboBox1KeyPress(TObject *Sender, char &Key);
        void __fastcall ToolButton1Click(TObject *Sender);
        void __fastcall ToolButton2Click(TObject *Sender);
        void __fastcall ToolButton3Click(TObject *Sender);
	void __fastcall ToolButton5Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
        TToolButton **ToolBtn;
        TFileStream **FS;
        __fastcall TForm3(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm3 *Form3;
//---------------------------------------------------------------------------
#endif
