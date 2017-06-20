//---------------------------------------------------------------------------

#ifndef RegimTableH
#define RegimTableH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <Grids.hpp>
#include "graph_cpp.h"
#include <Buttons.hpp>
//#include "reg.h"
//#include "graph_cpp.cpp"
//---------------------------------------------------------------------------
class TFormR : public TForm
{
__published:	// IDE-managed Components
        TPageControl *PageControl;
        TTabSheet *TabSheetNodes;
        TTabSheet *TabSheetLinks;
        TPanel *PanelButtons;
        TButton *Button2;
        TStringGrid *StringGrid1;
        TStringGrid *StringGrid2;
        TTabSheet *TabAdditional;
        TButton *NewB;
        TButton *DelB;
        TTabSheet *TabParams;
        TStringGrid *StringGrid3;
        TTabSheet *TabRusults;
        TTabSheet *TabRusultsNodes;
        TTabSheet *TabRusultsLinks;
        TStringGrid *StringGrid5;
        TStringGrid *StringGrid4;
        TLabel *Label2;
        TComboBox *ComboBox1;
        TGroupBox *GroupCommon;
        TLabel *Label3;
        TEdit *EditMaxIter;
        TLabel *Label4;
        TEdit *EditPrec;
        TGroupBox *GroupPredel;
        TGroupBox *GroupLimits;
        TCheckBox *CheckQminQmax;
        TRadioGroup *RadioGroupPredelType;
        TGroupBox *GroupBox2;
        TLabel *Label7;
        TButton *Button4;
        TCheckBox *CheckAutoVr0;
        TLabel *Label6;
        TLabel *Label5;
        TEdit *EditKn0;
        TEdit *EditT0;
        TGroupBox *GroupBox3;
        TCheckBox *CheckProm;
        TLabel *Label8;
        TEdit *EditdT;
        TGroupBox *GroupNebal;
        TLabel *Label11;
        TEdit *EditdUmax;
        TLabel *Label12;
        TEdit *Editdfmax;
        TLabel *Label13;
        TEdit *EditdKnmax;
        TLabel *Label14;
        TEdit *EditdTmax;
        TCheckBox *CheckNebal;
        TCheckBox *CheckTxtFile;
        TCheckBox *CheckGraph;
        TCheckBox *CheckGrad;
        TLabel *Label9;
        TEdit *EditdKn;
        TLabel *Label10;
        TCheckBox *CheckSensor;
        TGroupBox *GroupResCommon;
        TLabel *LabelRegim;
        TLabel *LabelIters;
        TLabel *LabelItersCaption;
        TGroupBox *GroupBox1;
        TLabel *Label15;
        TLabel *LabelResSensor;
        TLabel *Label18;
        TLabel *Label19;
        TStringGrid *StringGridSensor;
        TGroupBox *GroupBox4;
        TLabel *LabelResKn;
        TLabel *LabelTpred;
        TLabel *LabelResTpred;
        TStringGrid *StringGridVr;
        TLabel *Label24;
        TLabel *Label25;
        TStringGrid *StringGridGradLink;
        TLabel *Label26;
        TLabel *LabelResGradLink;
        TButton *ButtonRaschet;
        TLabel *Label1;
        TEdit *EditdVrimax;
        TButton *Button1;
        TButton *Button3;
        TButton *Button5;
        TButton *Button6;
        TCheckBox *CheckStarFiles;
        TBitBtn *ButtStarPath;
        TButton *Button7;
        TCheckBox *CheckBkUse;
        TCheckBox *CheckBox2;
        TLabel *Label16;
        TCheckBox *CheckCurrReg;
        TCheckBox *CheckParamsToFile;
        TBitBtn *ButtPrmsFilePath;
        TLabel *Label17;
        TLabel *LabelKn;
        TGroupBox *GroupBox5;
        TStringGrid *StringGrid6;
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall PageControlChange(TObject *Sender);
        void __fastcall NewBClick(TObject *Sender);
        void __fastcall DelBClick(TObject *Sender);
        void __fastcall ComboBox1Change(TObject *Sender);
        void __fastcall FormActivate(TObject *Sender);
        void __fastcall CheckNebalClick(TObject *Sender);
        void __fastcall ButtonRaschetClick(TObject *Sender);
        void __fastcall Button4Click(TObject *Sender);
        void __fastcall CheckPromClick(TObject *Sender);
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall Button3Click(TObject *Sender);
        void __fastcall Button5Click(TObject *Sender);
        void __fastcall ButtStarPathClick(TObject *Sender);
        void __fastcall Button6Click(TObject *Sender);
        void __fastcall Button7Click(TObject *Sender);
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall ButtPrmsFilePathClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
//        Reg *RRasch;
//        Values *Vals;

        int RaschetType;
        int Iters;
        int MaxIters;
        double prec;
        bool Showtxt;
        bool FindSensor;
        bool UseQminQmax;
        bool LimitNebals;
        double dUmax;
        double dfmax;
        double dKnmax;
        double dTmax;
        double dVrimax;
        int PredelRaschetMetod;
        double Kn0;
        double T0;
        bool RaschetVr0;
        double *Vr;
        double *VrRes;
        double *VrResTrasp;
        int nVr;
        bool TempRegims;
        bool Gradients;
        bool ShowGraph;
        bool SaveStarFiles;
        bool SavePrmsToFile;
        double dT;
        double dKn;
        int nP;
        int nU;
        int nJ,nb,ng,nuzl;
        int *mhu;

        double T;
        double Kn;
        bool newVr;
        double *dSdi;
        double *Svr;//,*Svi;

        double *VrTemp;
        int nVrTemp;
        bool BkUse;
        bool CurrReg;
        bool newSchema;
        int OldRegimType;

        AnsiString StarFilesPath;
        AnsiString PrmsFileName;        

        void FillNodes(TStringGrid *SG1,TStringGrid *SG2);
        void FillLinks(TStringGrid *SG);
        void Sort_i(TStringGrid *StringGrid);
        void Sort_ij(TStringGrid *StringGrid);
        void WriteRRasch(bool CurrReg);
        float MyStrToFloat(AnsiString str);
        int FindMaxNodeNumber(TStringGrid *StringGrid,int ColNumber,int CurrMax);
        void FillResults(int ItNumber);
        void Form1FillParams();
        void ClearGrid (TStringGrid *StringGrid);
        bool ReadRaschetType();
        int FillVr(bool &FillVr);
        int FinddSdiMax();
        void RaschetSv(int &nJlocal,int &nblocal,int &nglocal,int &nlocal,bool transpMatr);
        void RaschetVrResTrasp(bool transpMatr);
        void FillSensorGrid();
        void FillVrResToGrid();
        void FilldSdiToGrid();
        void RegimRaschet();

        void SaveStringGrid(TStringGrid *CurrSG,int type,AnsiString FileName);
        void OpenStringGrid(TStringGrid *CurrSG,int type);
        __fastcall TFormR(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFormR *FormR;
//---------------------------------------------------------------------------
#endif
