//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
USERES("graph_prj.res");
USEFORM("graph_cpp.cpp", Form1);
USEFORM("Params.cpp", Form2);
USEFORM("Tools.cpp", Form3);
USEFORM("rez_grahp.cpp", Form4);
USEFORM("graph_parm.cpp", Form5);
USEFORM("elembase.cpp", Form6);
USEUNIT("reg.cpp");
USEUNIT("det.cpp");
USEUNIT("sobstvqri.cpp");
USEFORM("RegimTable.cpp", FormR);
USEFORM("VectorR0.cpp", FormVectorVr);
USEUNIT("values.cpp");
USEFORM("selectdir.cpp", FormSelectDir);
USEFORM("ContactR.cpp", FormContactR);
USEUNIT("svd.cpp");
USEUNIT("svd_old.cpp");
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
        try
        {
                 Application->Initialize();
                 Application->Title = "Interval Calculation 2";
		Application->CreateForm(__classid(TForm1), &Form1);
                 Application->CreateForm(__classid(TForm2), &Form2);
                 Application->CreateForm(__classid(TForm3), &Form3);
                 Application->CreateForm(__classid(TForm5), &Form5);
                 Application->CreateForm(__classid(TForm6), &Form6);
                 Application->CreateForm(__classid(TFormR), &FormR);
                 Application->CreateForm(__classid(TFormVectorVr), &FormVectorVr);
                 Application->CreateForm(__classid(TFormSelectDir), &FormSelectDir);
                 Application->CreateForm(__classid(TFormContactR), &FormContactR);
                 Application->Run();
        }

        catch (Exception &exception)
        {

                 Application->ShowException(&exception);
        }
        return 0;
}
//---------------------------------------------------------------------------
