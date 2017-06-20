
#include <iostream>
#include <complex>

//---------------------------------------------------------------------------
using namespace std;
class Values
{
public:
       static int mswNumber;
       static int NodesNumber;
/*
       int RealNode[MaxNodesNumber];
       int Node1[MaxLinksNumber],Node2[MaxLinksNumber];
       bool IsLine[MaxLinksNumber];
       bool IsTrf[MaxLinksNumber];

       complex<double>u[MaxNodesNumber];
       complex<double>I[MaxLinksNumber];
       complex<double>S1[MaxNodesNumber];
       complex<double>S2[MaxNodesNumber];
       complex<double>DS[MaxLinksNumber];
       complex<double>DSTM[MaxLinksNumber];
       complex<double>DSTC[MaxLinksNumber];
       double DU[MaxLinksNumber];
       double Ql[MaxLinksNumber];
*/
       int *RealNode;
       int *Node1,*Node2;
       bool *IsLine;
       bool *IsTrf;

       complex<double> *u;
       complex<double> *I;
       complex<double> *S1;
       complex<double> *S2;
       complex<double> *DS;
       complex<double> *DSTM;
       complex<double> *DSTC;
       double *DU;
       double *Ql;

       bool RegimOk;
       Values();
       ~Values();
       void Initialize(int Nodes,int Links);
};
int Values::NodesNumber=0;
int Values::mswNumber=0;

extern Values *Vals;
extern Values *ValsGraph;
