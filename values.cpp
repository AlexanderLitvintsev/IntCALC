#include "Values.h"
Values *Vals;
Values *ValsGraph;

Values::Values()
{
/*
mswNumber=0;
NodesNumber=0;
for(int i=0;i<MaxNodesNumber;i++)
{
RealNode[i]=-1;
u[i]=S1[i]=S2[i]=complex<double>(0,0);
}
for(int i=0;i<MaxLinksNumber;i++)
{
Node1[i]=Node2[i]=-1;
I[i]=DS[i]=DSTM[i]=DSTC[i]=complex<double>(0,0);
DU[i]=0;
Ql[i]=0;
IsLine[i]=IsTrf[i]=false;
}
RegimOk=false;
*/
//mswNumber=0;
//NodesNumber=0;
RealNode=new int[NodesNumber+2];
Node1=new int[mswNumber+2];
Node2=new int[mswNumber+2];
IsLine=new bool[mswNumber+2];
IsTrf=new bool[mswNumber+2];
u=new complex<double>[NodesNumber+2];
I=new complex<double>[mswNumber+2];
S1=new complex<double>[mswNumber+2];
S2=new complex<double>[mswNumber+2];
DS=new complex<double>[mswNumber+2];
DSTM=new complex<double>[mswNumber+2];
DSTC=new complex<double>[mswNumber+2];
DU=new double[mswNumber+2];
Ql=new double[mswNumber+2];
for(int i=0;i<NodesNumber+2;i++)
{
RealNode[i]=-1;
u[i]=complex<double>(0,0);
}
for(int i=0;i<mswNumber+2;i++)
{
Node1[i]=Node2[i]=-1;
I[i]=S1[i]=S2[i]=DS[i]=DSTM[i]=DSTC[i]=complex<double>(0,0);
DU[i]=0;
Ql[i]=0;
IsLine[i]=IsTrf[i]=false;
}
RegimOk=false;

}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Values::~Values()
{
delete []RealNode;
delete []Node1;
delete []Node2;
delete []IsLine;
delete []IsTrf;
delete []u;
delete []I;
delete []S1;
delete []S2;
delete []DS;
delete []DSTM;
delete []DSTC;
delete []DU;
delete []Ql;

}
//---------------------------------------------------------------------------
void Values::Initialize(int Nodes,int Links)
{
delete[] RealNode;
delete[] Node1;
delete[] Node2;
delete[] IsLine;
delete[] IsTrf;
delete[] u;
delete[] I;
delete[] S1;
delete[] S2;
delete[] DS;
delete[] DSTM;
delete[] DSTC;
delete[] DU;
delete[] Ql;
RealNode=new int[Nodes+1];
Node1=new int[Links+1];
Node2=new int[Links+1];
IsLine=new bool[Links+1];
IsTrf=new bool[Links+1];
u=new complex<double>[Nodes+1];
I=new complex<double>[Links+1];
S1=new complex<double>[Links+1];
S2=new complex<double>[Links+1];
DS=new complex<double>[Links+1];
DSTM=new complex<double>[Links+1];
DSTC=new complex<double>[Links+1];
DU=new double[Links+1];
Ql=new double[Links+1];
for(int i=0;i<Nodes+1;i++)
{
RealNode[i]=-1;
u[i]=complex<double>(0,0);
}
for(int i=0;i<Links+1;i++)
{
Node1[i]=Node2[i]=-1;
I[i]=S1[i]=S2[i]=DS[i]=DSTM[i]=DSTC[i]=complex<double>(0,0);
DU[i]=0;
Ql[i]=0;
IsLine[i]=IsTrf[i]=false;
}
RegimOk=false;
//mswNumber=Links;
//NodesNumber=Nodes;
}
//---------------------------------------------------------------------------
