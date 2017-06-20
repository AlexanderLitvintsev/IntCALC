// Author: Alexander Litvintsev <alexanderlitvintsev@yahoo.com>

#include "intcalc2c.h"

char math_type[2];
int iter_max = 20;
ap::complex_2d_array Yc1_temp_1, Yc2_temp_1, Yc1_temp_2, Yc2_temp_2, Yc1_temp_3, Yc2_temp_3, Yc1_temp_4, Yc2_temp_4;
ap::complex_1d_array temp1_1, temp2_1, temp3_1, temp1_2, temp2_2, temp3_2;
ap::complex_2d_array TEMP_U_1, TEMP_U_2; // for get_K2U();
double U_a1_1_mod,U_b1_1_mod,U_c1_1_mod,U_a2_1_mod,U_b2_1_mod,U_c2_1_mod,U_a3_1_mod,U_b3_1_mod,U_c3_1_mod;
double U_a1_2_mod,U_b1_2_mod,U_c1_2_mod,U_a2_2_mod,U_b2_2_mod,U_c2_2_mod,U_a3_2_mod,U_b3_2_mod,U_c3_2_mod;
double U4f_1_mod, U4f_2_mod, U5f_1_mod, U5f_2_mod, U6f_1_mod, U6f_2_mod;
ap::complex S1, S2, S3;
ap::complex_1d_array MF_U_low, MF_U_high;
ap::complex_2d_array ISC_U_unsim, MF_ZM_low, MF_ZM_high;
char seps[] = " .,!?:;\n";
long double pi = acos((long double) -1);
double INIT_S = 20.0;
double IP_P = 20.0;
double IP_Q = 20.0;
double IP_P1 = 20.0;
double IP_Q1 = 20.0;
double IP_P2 = 20.0;
double IP_Q2 = 20.0;
double IP_P3 = 20.0;
double IP_Q3 = 20.0;

double IP_P1_A = 0.0;
double IP_Q1_A = 0.0;
double IP_P1_B = 0.0;
double IP_Q1_B = 0.0;
double IP_P1_C = 0.0;
double IP_Q1_C = 0.0;
double IP_P2_A = 0.0;
double IP_Q2_A = 0.0;
double IP_P2_B = 0.0;
double IP_Q2_B = 0.0;
double IP_P2_C = 0.0;
double IP_Q2_C = 0.0;
double IP_P3_A = 0.0;
double IP_Q3_A = 0.0;
double IP_P3_B = 0.0;
double IP_Q3_B = 0.0;
double IP_P3_C = 0.0;
double IP_Q3_C = 0.0;

double IP_P1_A_low = 0.0;
double IP_P1_A_high = 0.0;
double IP_Q1_A_low = 0.0;
double IP_Q1_A_high = 0.0;
double IP_P1_B_low = 0.0;
double IP_P1_B_high = 0.0;
double IP_Q1_B_low = 0.0;
double IP_Q1_B_high = 0.0;
double IP_P1_C_low = 0.0;
double IP_P1_C_high = 0.0;
double IP_Q1_C_low = 0.0;
double IP_Q1_C_high = 0.0;
double IP_P2_A_low = 0.0;
double IP_P2_A_high = 0.0;
double IP_Q2_A_low = 0.0;
double IP_Q2_A_high = 0.0;
double IP_P2_B_low = 0.0;
double IP_P2_B_high = 0.0;
double IP_Q2_B_low = 0.0;
double IP_Q2_B_high = 0.0;
double IP_P2_C_low = 0.0;
double IP_P2_C_high = 0.0;
double IP_Q2_C_low = 0.0;
double IP_Q2_C_high = 0.0;
double IP_P3_A_low = 0.0;
double IP_P3_A_high = 0.0;
double IP_Q3_A_low = 0.0;
double IP_Q3_A_high = 0.0;
double IP_P3_B_low = 0.0;
double IP_P3_B_high = 0.0;
double IP_Q3_B_low = 0.0;
double IP_Q3_B_high = 0.0;
double IP_P3_C_low = 0.0;
double IP_P3_C_high = 0.0;
double IP_Q3_C_low = 0.0;
double IP_Q3_C_high = 0.0;

double IP_P4_A = 0.0;
double IP_Q4_A = 0.0;
double IP_P4_B = 0.0;
double IP_Q4_B = 0.0;
double IP_P4_C = 0.0;
double IP_Q4_C = 0.0;

double IP_L_mistake = 0.00;
double IP_Ua_l = 112.58330249197702407928401219788;
double IP_Ub_l = 112.58330249197702407928401219788;
double IP_Uc_l = 112.58330249197702407928401219788;
double Length = 50.0;
double Height = 19.0;
double Yground = 0.01;
double X1 = - 2.0;
double X2 = 1.0;
double X3 = 3.0;
double Y1 = 19.0;
double Y2 = 21.5;
double Y3 = 19.0;
double temperature=0.0;
double RP_Umod_min[3], RP_Umod_max[3];
cinterval RP_Yc[6][6], RP_Z[3][3];
interval IP_L;
cinterval iGIG_A;
interval iGIG_GRAD;
int triangle_ptype;
int param_interval;
int param_unsim;
int phase_fault = 0; // обрыв фазы: 0 - нет, 1 - обрыв фазы A
//===============================================================================================================
double SC_UL;
double SC_L;
double SC_UL_1;
double SC_UL_2;
double SC_L_1;
double SC_L_2;
double SC_X0;
double SC_X0_1;
double SC_X0_2;
interval iSC_UL, iSC_UF, iSC_L, iSC_X0, iSC_X1;
cinterval iSC_IA, iSC_IA1, iSC_IA2, iSC_IB, iSC_IC, iSC_IS, iSC_temp_01, iSC_temp_02, iSC_a, iSC_UA1, iSC_UA2;
cinterval iSC_UA, iSC_UB, iSC_UC;

double BAR_UA[419904];
double BAR_UB[419904];
double BAR_UC[419904];

//===============================================================================================================
/*
void cimatrixinv(cimatrix M_A, cimatrix M_X)
{
cimatrix M_J(3,3);
Re(M_J(1,1)) = 1.0;
Im(M_J(1,1)) = 0.0;
Re(M_J(2,2)) = 1.0;
Im(M_J(2,2)) = 0.0;
Re(M_J(3,3)) = 1.0;
Im(M_J(3,3)) = 0.0;
cout << Re(M_A(0,0));

	for (int i=0; i<100; i++)
	{
	M_X = mid(M_X) + M_X * (M_J - M_A * mid(M_X));
	}

cout << Re(M_X(1,1)) << "\n";
}
*/

void readfile( char *FileNam )
{
	/*
		FILE *txt;
        int i=0, count = 0;
        char *text;
        char temp[255];
        char *words[255] = {NULL};

        if((txt=fopen("text.txt","r")) == NULL) 
        {
        cout << "Ошибка при открытии файла.";
        }
		else
        {
                cout << "Файл успешно открылся!";
                fscanf(txt, "%s", temp); 
                text = new char [strlen(temp) +1];
                strcpy (text, temp);
                do
                {
                        words[count] = strtok( text, seps );
 
                        while (words[count] != NULL)
                        {
                                count++;                                
                                words[count] = strtok( NULL, seps );
                        }
 
                        while (i<count)
                        {
                                puts(words[i]);
                                i++;
                        }
                }while ((fscanf(txt, "%s", text)) != EOF);
        }
        fclose(txt)
*/
}
////////////////////////////////////////////////////////////////////////////////////////
void get_interval_grad()
{


if (Inf(Re(iGIG_A)) > NULL && Sup(Re(iGIG_A)) > NULL)
{
iGIG_GRAD = atan(Im(iGIG_A) / Re(iGIG_A)) * (180.0 / Pi());
//cout << "a > 0" << "\n";
}
else if (Inf(Re(iGIG_A)) < NULL && Sup(Re(iGIG_A)) < NULL && Inf(Im(iGIG_A)) > NULL && Sup(Im(iGIG_A)) > NULL)
{
iGIG_GRAD = 180.0 + atan(Im(iGIG_A) / Re(iGIG_A)) * (180.0 / Pi());
//cout << "a < 0, b > 0" << "\n";
}
else if (Inf(Re(iGIG_A)) < NULL && Sup(Re(iGIG_A)) < NULL && Inf(Im(iGIG_A)) < NULL && Sup(Im(iGIG_A)) < NULL)
{
iGIG_GRAD = -180.0 + atan(Im(iGIG_A) / Re(iGIG_A)) * (180.0 / Pi());
//cout << "a < 0, b < 0" << "\n";
}	
else cout << "ERROR void get_interval_grad()" << "\n";
/*
interval iMF_Ey_grad = atan(Im(iMF_Ey_2) / Re(iMF_Ey_2)) * (180.0 / Pi());
interval iMF_Hx_grad = atan(Im(iMF_Hx_2) / Re(iMF_Hx_2)) * (180.0 / Pi());
interval iMF_Hy_grad = atan(Im(iMF_Hy_2) / Re(iMF_Hy_2)) * (180.0 / Pi());
*/
}
////////////////////////////////////////////////////////////////////////////////////////

void get_PQ_5()
{
double IP_MISTAKE_MIN = 0.95;
double IP_MISTAKE_MAX = 1.05;
/////////////////////////////
IP_P1_A_low = IP_P1_A_low * IP_MISTAKE_MIN;
IP_P1_A_high = IP_P1_A_high * IP_MISTAKE_MAX;
IP_Q1_A_low = IP_Q1_A_low * IP_MISTAKE_MIN;
IP_Q1_A_high = IP_Q1_A_high * IP_MISTAKE_MAX;

IP_P1_B_low = IP_P1_B_low * IP_MISTAKE_MIN;
IP_P1_B_high = IP_P1_B_high * IP_MISTAKE_MAX;
IP_Q1_B_low = IP_Q1_B_low * IP_MISTAKE_MIN;
IP_Q1_B_high = IP_Q1_B_high * IP_MISTAKE_MAX;

IP_P1_C_low = IP_P1_C_low * IP_MISTAKE_MIN;
IP_P1_C_high = IP_P1_C_high * IP_MISTAKE_MAX;
IP_Q1_C_low = IP_Q1_C_low * IP_MISTAKE_MIN;
IP_Q1_C_high = IP_Q1_C_high * IP_MISTAKE_MAX;

/////////////////////////////
IP_P2_A_low = IP_P2_A_low * IP_MISTAKE_MIN;
IP_P2_A_high = IP_P2_A_high * IP_MISTAKE_MAX;
IP_Q2_A_low = IP_Q2_A_low * IP_MISTAKE_MIN;
IP_Q2_A_high = IP_Q2_A_high * IP_MISTAKE_MAX;


IP_P2_B_low = IP_P2_B_low * IP_MISTAKE_MIN;
IP_P2_B_high = IP_P2_B_high * IP_MISTAKE_MAX;
IP_Q2_B_low = IP_Q2_B_low * IP_MISTAKE_MIN;
IP_Q2_B_high = IP_Q2_B_high * IP_MISTAKE_MAX;

IP_P2_C_low = IP_P2_C_low * IP_MISTAKE_MIN;
IP_P2_C_high = IP_P2_C_high * IP_MISTAKE_MAX;
IP_Q2_C_low = IP_Q2_C_low * IP_MISTAKE_MIN;
IP_Q2_C_high = IP_Q2_C_high * IP_MISTAKE_MAX;

/////////////////////////////

IP_P3_A_low = IP_P3_A_low * IP_MISTAKE_MIN;
IP_P3_A_high = IP_P3_A_high * IP_MISTAKE_MAX;
IP_Q3_A_low = IP_Q3_A_low * IP_MISTAKE_MIN;
IP_Q3_A_high = IP_Q3_A_high * IP_MISTAKE_MAX;

IP_P3_B_low = IP_P3_B_low * IP_MISTAKE_MIN;
IP_P3_B_high = IP_P3_B_high * IP_MISTAKE_MAX;
IP_Q3_B_low = IP_Q3_B_low * IP_MISTAKE_MIN;
IP_Q3_B_high = IP_Q3_B_high * IP_MISTAKE_MAX;

IP_P3_C_low = IP_P3_C_low * IP_MISTAKE_MIN;
IP_P3_C_high = IP_P3_C_high * IP_MISTAKE_MAX;
IP_Q3_C_low = IP_Q3_C_low * IP_MISTAKE_MIN;
IP_Q3_C_high = IP_Q3_C_high * IP_MISTAKE_MAX;

/////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_PQ_10()
{
double IP_MISTAKE_MIN = 0.90;
double IP_MISTAKE_MAX = 1.10;
/////////////////////////////
IP_P1_A_low = IP_P1_A_low * IP_MISTAKE_MIN;
IP_P1_A_high = IP_P1_A_high * IP_MISTAKE_MAX;
IP_Q1_A_low = IP_Q1_A_low * IP_MISTAKE_MIN;
IP_Q1_A_high = IP_Q1_A_high * IP_MISTAKE_MAX;

IP_P1_B_low = IP_P1_B_low * IP_MISTAKE_MIN;
IP_P1_B_high = IP_P1_B_high * IP_MISTAKE_MAX;
IP_Q1_B_low = IP_Q1_B_low * IP_MISTAKE_MIN;
IP_Q1_B_high = IP_Q1_B_high * IP_MISTAKE_MAX;

IP_P1_C_low = IP_P1_C_low * IP_MISTAKE_MIN;
IP_P1_C_high = IP_P1_C_high * IP_MISTAKE_MAX;
IP_Q1_C_low = IP_Q1_C_low * IP_MISTAKE_MIN;
IP_Q1_C_high = IP_Q1_C_high * IP_MISTAKE_MAX;

/////////////////////////////
IP_P2_A_low = IP_P2_A_low * IP_MISTAKE_MIN;
IP_P2_A_high = IP_P2_A_high * IP_MISTAKE_MAX;
IP_Q2_A_low = IP_Q2_A_low * IP_MISTAKE_MIN;
IP_Q2_A_high = IP_Q2_A_high * IP_MISTAKE_MAX;


IP_P2_B_low = IP_P2_B_low * IP_MISTAKE_MIN;
IP_P2_B_high = IP_P2_B_high * IP_MISTAKE_MAX;
IP_Q2_B_low = IP_Q2_B_low * IP_MISTAKE_MIN;
IP_Q2_B_high = IP_Q2_B_high * IP_MISTAKE_MAX;

IP_P2_C_low = IP_P2_C_low * IP_MISTAKE_MIN;
IP_P2_C_high = IP_P2_C_high * IP_MISTAKE_MAX;
IP_Q2_C_low = IP_Q2_C_low * IP_MISTAKE_MIN;
IP_Q2_C_high = IP_Q2_C_high * IP_MISTAKE_MAX;

/////////////////////////////

IP_P3_A_low = IP_P3_A_low * IP_MISTAKE_MIN;
IP_P3_A_high = IP_P3_A_high * IP_MISTAKE_MAX;
IP_Q3_A_low = IP_Q3_A_low * IP_MISTAKE_MIN;
IP_Q3_A_high = IP_Q3_A_high * IP_MISTAKE_MAX;

IP_P3_B_low = IP_P3_B_low * IP_MISTAKE_MIN;
IP_P3_B_high = IP_P3_B_high * IP_MISTAKE_MAX;
IP_Q3_B_low = IP_Q3_B_low * IP_MISTAKE_MIN;
IP_Q3_B_high = IP_Q3_B_high * IP_MISTAKE_MAX;

IP_P3_C_low = IP_P3_C_low * IP_MISTAKE_MIN;
IP_P3_C_high = IP_P3_C_high * IP_MISTAKE_MAX;
IP_Q3_C_low = IP_Q3_C_low * IP_MISTAKE_MIN;
IP_Q3_C_high = IP_Q3_C_high * IP_MISTAKE_MAX;

/////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_K2U()
{

complex cSC_1_3, cSC_a, cSC_aa, cSC_2_3_of_Pi;
Re(cSC_1_3) = double(1.0 / 3.0);
Im(cSC_1_3) = 0.0;
Re(cSC_2_3_of_Pi) = 0.0;
Im(cSC_2_3_of_Pi) = double((2.0 / 3.0) * pi);
cSC_a = exp(cSC_2_3_of_Pi);
cSC_aa = pow(cSC_a, 2.0);
ap::complex_2d_array SC_U_1, SC_U_2, ISC_S, ISC_U120_1, ISC_U120_2, ISC_TEMP_1, ISC_TEMP_2;
ap::complex ISC_1_3;
SC_U_1.setlength(3,1);
SC_U_2.setlength(3,1);
ISC_TEMP_1.setlength(3,1);
ISC_TEMP_2.setlength(3,1);
ISC_S.setlength(3,3);
ISC_U120_1.setlength(3,1);
ISC_U120_2.setlength(3,1);
ISC_1_3.x = _double(Re(cSC_1_3));
ISC_1_3.y = 0.0;

for (int i = 0; i < 3; i++)      
{
   for (int j = 0; j < 1; j++)  
   {
	SC_U_1(i,j) = TEMP_U_1(i,j);
	SC_U_2(i,j) = TEMP_U_2(i,j);
   }
}
/*
SC_U_1(0,0).x = 133.0;
SC_U_1(0,0).y = 0.0;
SC_U_1(1,0).x = 108.286;
SC_U_1(1,0).y = -77.221;
SC_U_1(2,0).x = 108.286;
SC_U_1(2,0).y = 77.221;
SC_U_2(0,0).x = 133.0;
SC_U_2(0,0).y = 0.0;
SC_U_2(1,0).x = 108.286;
SC_U_2(1,0).y = -77.221;
SC_U_2(2,0).x = 108.286;
SC_U_2(2,0).y = 77.221;
*/
ISC_S(0,0).x = 1.0;
ISC_S(0,0).y = 0.0;
ISC_S(0,1).x = _double(Re(cSC_a));
ISC_S(0,1).y = _double(Im(cSC_a));
ISC_S(0,2).x = _double(Re(cSC_aa));
ISC_S(0,2).y = _double(Im(cSC_aa));
ISC_S(1,0).x = 1.0;
ISC_S(1,0).y = 0.0;
ISC_S(1,1).x = _double(Re(cSC_aa));
ISC_S(1,1).y = _double(Im(cSC_aa));
ISC_S(1,2).x = _double(Re(cSC_a));
ISC_S(1,2).y = _double(Im(cSC_a));
ISC_S(2,0).x = 1.0;
ISC_S(2,0).y = 0.0;
ISC_S(2,1).x = 1.0;
ISC_S(2,1).y = 0.0;
ISC_S(2,2).x = 1.0;
ISC_S(2,2).y = 0.0;

// [S] * [U]
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 ISC_TEMP_1(i,j) += ISC_S(i,k) * SC_U_1(k,j);
		 ISC_TEMP_2(i,j) += ISC_S(i,k) * SC_U_2(k,j);
	  }
   }
}
// (1/3) * [S] * [U]
for (int i = 0; i < 3; i++)      
{
   for (int j = 0; j < 1; j++)  
   {
	ISC_U120_1(i,j) = ISC_1_3 * ISC_TEMP_1(i,j);
	ISC_U120_2(i,j) = ISC_1_3 * ISC_TEMP_2(i,j);
   }
}
cinterval iISC_U120[3][1], iISC_K2U;
interval iISC_100;

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_U120[i][j])) = _real(min(ISC_U120_1(i,j).x, ISC_U120_2(i,j).x));
	Inf(Im(iISC_U120[i][j])) = _real(min(ISC_U120_1(i,j).y, ISC_U120_2(i,j).y));
	Sup(Re(iISC_U120[i][j])) = _real(max(ISC_U120_1(i,j).x, ISC_U120_2(i,j).x));
	Sup(Im(iISC_U120[i][j])) = _real(max(ISC_U120_1(i,j).y, ISC_U120_2(i,j).y));
	}
}
/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_U120[i][j])) = real(ISC_U120_1(i,j).x);
	Inf(Im(iISC_U120[i][j])) = real(ISC_U120_1(i,j).y);
	Sup(Re(iISC_U120[i][j])) = real(ISC_U120_2(i,j).x);
	Sup(Im(iISC_U120[i][j])) = real(ISC_U120_2(i,j).y);
	}
}
*/

/*
Inf(Re(iISC_U120[0][0])) = 5.0;
Inf(Im(iISC_U120[0][0])) = 5.0;
Sup(Re(iISC_U120[0][0])) = 5.0;
Sup(Im(iISC_U120[0][0])) = 5.0;

Inf(Re(iISC_U120[1][0])) = 3.0;
Inf(Im(iISC_U120[1][0])) = 3.0;
Sup(Re(iISC_U120[1][0])) = 3.0;
Sup(Im(iISC_U120[1][0])) = 3.0;
*/
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "U120 [" << i << "] = " << iISC_U120[i][j] << "\n";
	}
}


iISC_100 = 100.0;

//iISC_K2U = (abs(iISC_U120[1][0]) / abs(iISC_U120[0][0])) * iISC_100;

interval iISC_K2U_ADD, iISC_K2U_U120_ABS0, iISC_K2U_U120_ABS1;

iISC_K2U_U120_ABS0 = abs(iISC_U120[0][0]);
iISC_K2U_U120_ABS1 = abs(iISC_U120[1][0]);
/*
iISC_K2U_U120_ABS0.x =  _double (max ( abs(Inf(Re(iISC_U120[0][0]))), abs(Sup(Re(iISC_U120[0][0]))) ));
iISC_K2U_U120_ABS0.y =  _double (max ( abs(Inf(Im(iISC_U120[0][0]))), abs(Sup(Im(iISC_U120[0][0]))) ));
iISC_K2U_U120_ABS1.x =  _double (max ( abs(Inf(Re(iISC_U120[1][0]))), abs(Sup(Re(iISC_U120[1][0]))) ));
iISC_K2U_U120_ABS1.y =  _double (max ( abs(Inf(Im(iISC_U120[1][0]))), abs(Sup(Im(iISC_U120[1][0]))) ));
*/
cout << "ABS U120 [0] = " << iISC_K2U_U120_ABS0 << "\n";
cout << "ABS U120 [1] = " << iISC_K2U_U120_ABS1 << "\n";

iISC_K2U_ADD = ( iISC_K2U_U120_ABS1 / iISC_K2U_U120_ABS0 ) * iISC_100;
cout << "K2U = " << iISC_K2U_ADD << "\n";

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void getintvloss()
{

// http://lvvas.co.cc/index.php?option=com_content&view=article&id=135&Itemid=157
// initial parameters

double MagConst;

interval IP_MagConst, IP_WC_temperature, IP_Yground, IP_freq, IP_cfreq, IP_WC_R0, IP_WC_S, IP_WC_R,
		 IP_WC_R_Equal, IP_X1, IP_X2, IP_X3, IP_Y1, IP_Y2, IP_Y3;

interval RP_d12, RP_d13, RP_d23, RP_D12, RP_D13, RP_D23;

interval IP_U_A_l, IP_U_B_l, IP_U_C_l, IP_U_A_f, IP_U_B_f, IP_U_C_f;

cinterval IP_S;

Re(IP_S) = 15.0 * 1000000.0;
Im(IP_S) = 7.5 * 1000000.0;

MagConst = 0.00000125663706;

IP_U_A_l = IP_Ua_l * 1000.0;
IP_U_B_l = IP_Ub_l * 1000.0;
IP_U_C_l = IP_Uc_l * 1000.0;
//cout << "U" << IP_U_A_l;
IP_U_A_f = IP_U_A_l / sqrt(3.0);
IP_U_B_f = IP_U_B_l / sqrt(3.0);
IP_U_C_f = IP_U_C_l / sqrt(3.0);


if (param_interval == 0)
{

IP_MagConst = MagConst;	// магнитная постоянная
IP_Yground = Yground;//(0.05,0.3); // проводимость земли
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
IP_L = Length;
IP_WC_R0 = 0.1197*(1.0 + 0.004*(temperature-20.0));
IP_WC_R0 = 0.2;		// сопротивление 1 км провода, См/м
IP_WC_S = 150.0;	// 240.0
IP_WC_R = 0.0085;	// 0.012, радиус провода
IP_X1 = X1;
IP_X2 = X2;
IP_X3 = X3;
IP_Y1 = Y1;
IP_Y2 = Y2;
IP_Y3 = Y3;
}

if (param_interval == 1)
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_WC_temperature) = -40.0;	// -40
Sup(IP_WC_temperature) = 40.0;	// 40
Inf(IP_Yground) = 0.05; // проводимость земли
Sup(IP_Yground) = 0.3;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.1197*(1.0 + 0.004*(IP_WC_temperature-20.0));
IP_WC_S = 150.0;	// 240.0
IP_WC_R = 0.0085;	// 0.012
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}

if (param_interval == 2)	// fazanord min & max parameters
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_Yground) = 0.01; // проводимость земли
Sup(IP_Yground) = 0.01;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.122;
IP_WC_S = 240.0;
IP_WC_R = 0.012;
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}
//cout << IP_X1 << "\n";
//cout << IP_X2 << "\n";

//cout << "L = " << IP_L << " км" << "\n";
//cout << "L = " << IP_L << " км" << "\n";
//cout << IP_X2 << "\n";
interval IP_X12;
Inf(IP_X12) = Inf(IP_X1) - Inf(IP_X2);
Sup(IP_X12) = Sup(IP_X1) - Sup(IP_X2);
Inf(IP_X12) = Inf(IP_X12) * Inf(IP_X12);
Sup(IP_X12) = Sup(IP_X12) * Sup(IP_X12);

RP_d12 = sqrt(IP_X12 + Power(IP_Y1 - IP_Y2, 2.0));
RP_d13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 - IP_Y3), 2.0));
RP_d23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 - IP_Y3), 2.0));

RP_D12 = sqrt(IP_X12 + Power(IP_Y1 + IP_Y2, 2.0));
RP_D13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 + IP_Y3), 2.0));
RP_D23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 + IP_Y3), 2.0));
//RP_D12 = RP_D23; 
//cout << RP_D12 << "\n";
cinterval RP_Z11, RP_Z12, RP_Z13, RP_Z21, RP_Z22, RP_Z23, RP_Z31, RP_Z32, RP_Z33, RP_Z11out, RP_Z11inn;
interval pow_0_755, pow_0_83;
pow_0_755 = 0.755;
pow_0_83 = 0.83;

if (phase_fault == 0)
{

Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;
}

if (phase_fault == 1)
{
Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;

// -----------------------------------------------------------------------------------------
Re(RP_Z11) = 1000000.0;
Im(RP_Z11) = 0.0;
}

cinterval RP_d_p[3][3], RP_d_n[3][3];
RP_Z[0][0] = RP_Z11;
RP_Z[0][1] = RP_Z12;
RP_Z[0][2] = RP_Z13;
RP_Z[1][0] = RP_Z21;
RP_Z[1][1] = RP_Z22;
RP_Z[1][2] = RP_Z23;
RP_Z[2][0] = RP_Z31;
RP_Z[2][1] = RP_Z32;
RP_Z[2][2] = RP_Z33;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	RP_Z[i][j] = IP_L * RP_Z[i][j];

	}
}

//cout << RP_Z[0][0] << "\n";

ap::complex_2d_array Z_inf;
ap::complex_2d_array Z_sup;
ap::complex_2d_array D_inf;
ap::complex_2d_array D_sup;
Z_inf.setlength(3,3);
Z_sup.setlength(3,3);
D_inf.setlength(3,3);
D_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Z_inf(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	Z_inf(i,j).y = _double(Inf(Im(RP_Z[i][j])));
	Z_sup(i,j).x = _double(Sup(Re(RP_Z[i][j])));
	Z_sup(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	}
}
//cout << Z_inf(0,0).x << "\n";
cmatrixinverse(Z_inf,3);
cmatrixinverse(Z_sup,3);
//cout << Z_inf(0,0).x << "\n";
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_d_p[i][j])) = Z_inf(i,j).x;
	Inf(Im(RP_d_p[i][j])) = Z_inf(i,j).y;
	Sup(Re(RP_d_p[i][j])) = Z_sup(i,j).x;
	Sup(Im(RP_d_p[i][j])) = Z_sup(i,j).y;
	Inf(Re(RP_d_n[i][j])) = - Z_inf(i,j).x;// -Z
	Inf(Im(RP_d_n[i][j])) = - Z_inf(i,j).y;
	Sup(Re(RP_d_n[i][j])) = - Z_sup(i,j).x;
	Sup(Im(RP_d_n[i][j])) = - Z_sup(i,j).y;
	}
}

/*
//inv
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_Z[i][j]) = Inf(RP_Z[i][j]) / 100.0;
	Sup(RP_Z[i][j]) = Sup(RP_Z[i][j]) / 100.0;
	}
}

cinterval M_A[3][3], M_B[3][3], M_X[3][3], M_J[3][3];

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
		M_A[i][j] = RP_Z[i][j]; 
		Re(M_J[i][j]) = 0.0;
		Im(M_J[i][j]) = 0.0;
	}
}

Re(M_J[0][0]) = 1.0;
Im(M_J[0][0]) = 0.0;
Re(M_J[1][1]) = 1.0;
Im(M_J[1][1]) = 0.0;
Re(M_J[2][2]) = 1.0;
Im(M_J[2][2]) = 0.0;

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	M_B[i][j] = M_J[i][j] - M_A[i][j]; 
	}
}

cout << abs(M_B[0][0]) << "\n";


// норма комплексной матрицы

double M_B_norm_row, M_B_norm_col;

interval M_B_row_1, M_B_row_2, M_B_row_3;

//M_B_row_1;
for (int j=0; j<3;j++)
{
M_B_row_1 += abs(M_B[0][j]);
M_B_row_2 += abs(M_B[1][j]);
M_B_row_3 += abs(M_B[2][j]);
}

 
//cout << M_B_row_2 << "\n";

const int len = 3;
const int start = 0;
double arr[len];
arr[0] = _double(Sup(M_B_row_1));
arr[1] = _double(Sup(M_B_row_2));
arr[2] = _double(Sup(M_B_row_3));
double* first = arr;
double* last = arr + len;
double* asd;
asd = max_element(first,last);


cout << "Matrix norm" << *asd << "\n";

M_B_norm_col = 1.0 / (1.0 - *asd);


Inf(M_X[0][0]) = - M_B_norm_col;
Sup(M_X[0][0]) = M_B_norm_col + 2.0;
Inf(M_X[0][1]) = - M_B_norm_col;
Sup(M_X[0][1]) = M_B_norm_col;
Inf(M_X[0][2]) = - M_B_norm_col;
Sup(M_X[0][2]) = M_B_norm_col;
Inf(M_X[1][0]) = - M_B_norm_col;
Sup(M_X[1][0]) = M_B_norm_col;
Inf(M_X[1][1]) = - M_B_norm_col;
Sup(M_X[1][1]) = M_B_norm_col + 2.0;
Inf(M_X[1][2]) = - M_B_norm_col;
Sup(M_X[1][2]) = M_B_norm_col;
Inf(M_X[2][0]) = - M_B_norm_col;
Sup(M_X[2][0]) = M_B_norm_col;
Inf(M_X[2][1]) = - M_B_norm_col;
Sup(M_X[2][1]) = M_B_norm_col;
Inf(M_X[2][2]) = - M_B_norm_col;
Sup(M_X[2][2]) = M_B_norm_col + 2.0;
cout << "M_X" << M_X[0][0] << "\n";

cimatrix iM_X(3,3), iM_J(3,3), iM_A(3,3);

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	iM_X(i+1,j+1) = M_X[i][j];
	}
}
cout << "M_X" << iM_X(1,1) << "\n";


	for (int i=0; i<10; i++)
	{
	M_X[i][j] = mid(M_X) + M_X * (M_J - M_A * mid(M_X));
	}
	
//cout << Re(M_X(1,1)) << "\n";
*/
// Yrc

cinterval RP_Yrc[6][6];
RP_Yrc[0][0] = RP_d_p[0][0];
RP_Yrc[0][1] = RP_d_p[0][1];
RP_Yrc[0][2] = RP_d_p[0][2];
RP_Yrc[0][3] = RP_d_n[0][3];
RP_Yrc[0][4] = RP_d_n[0][1];
RP_Yrc[0][5] = RP_d_n[0][2];
RP_Yrc[1][0] = RP_d_p[1][0];
RP_Yrc[1][1] = RP_d_p[1][1];
RP_Yrc[1][2] = RP_d_p[1][2];
RP_Yrc[1][3] = RP_d_n[1][0];
RP_Yrc[1][4] = RP_d_n[1][1];
RP_Yrc[1][5] = RP_d_n[1][2];
RP_Yrc[2][0] = RP_d_p[2][0];
RP_Yrc[2][1] = RP_d_p[2][1];
RP_Yrc[2][2] = RP_d_p[2][2];
RP_Yrc[2][3] = RP_d_n[2][0];
RP_Yrc[2][4] = RP_d_n[2][1];
RP_Yrc[2][5] = RP_d_n[2][2];
RP_Yrc[3][0] = RP_d_n[0][0];
RP_Yrc[3][1] = RP_d_n[0][1];
RP_Yrc[3][2] = RP_d_n[0][2];
RP_Yrc[3][3] = RP_d_p[0][0];
RP_Yrc[3][4] = RP_d_p[0][1];
RP_Yrc[3][5] = RP_d_p[0][2];
RP_Yrc[4][0] = RP_d_n[1][0];
RP_Yrc[4][1] = RP_d_n[1][1];
RP_Yrc[4][2] = RP_d_n[1][2];
RP_Yrc[4][3] = RP_d_p[1][0];
RP_Yrc[4][4] = RP_d_p[1][1];			
RP_Yrc[4][5] = RP_d_p[1][2];
RP_Yrc[5][0] = RP_d_n[2][0];
RP_Yrc[5][1] = RP_d_n[2][1];
RP_Yrc[5][2] = RP_d_n[2][2];
RP_Yrc[5][3] = RP_d_p[2][0];
RP_Yrc[5][4] = RP_d_p[2][1];
RP_Yrc[5][5] = RP_d_p[2][2];

const int n_ABC=3;
interval RP_A[n_ABC][n_ABC], RP_B[n_ABC][n_ABC], RP_C[n_ABC][n_ABC];

RP_A[0][0] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y1) / (IP_WC_R * 100.0));
RP_A[0][1] = 1.8 * power(10.0, 7.0) * ln(RP_D12 / RP_d12);
RP_A[0][2] = 1.8 * power(10.0, 7.0) * ln(RP_D13 / RP_d13);
RP_A[1][0] = RP_A[0][1];
RP_A[1][1] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y2) / (IP_WC_R * 100.0));
RP_A[1][2] = 1.8 * power(10.0, 7.0) * ln(RP_D23 / RP_d23);
RP_A[2][0] = RP_A[0][2];
RP_A[2][1] = RP_A[1][2];
RP_A[2][2] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y3) / (IP_WC_R * 100.0));

ap::real_2d_array a_inf;
ap::real_2d_array a_sup;
a_inf.setlength(3,3);
a_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	a_inf(i,j) = _double(Inf(RP_A[i][j]));
	a_sup(i,j) = _double(Sup(RP_A[i][j]));
	}
}
//cout << a_inf(0,0) << "\n";
rmatrixinverse(a_inf,3);
rmatrixinverse(a_sup,3);
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_B[i][j]) = a_inf(i,j);
	Sup(RP_B[i][j]) = a_sup(i,j);
	}
}

/*
RP_C[0][0] = RP_B[0][0] + RP_B[0][1] + RP_B[0][2];
RP_C[0][1] = - RP_B[0][1];
RP_C[0][2] = - RP_B[0][2];
RP_C[1][0] = - RP_B[1][0];
RP_C[1][1] = RP_B[1][0] + RP_B[1][1] + RP_B[1][2];
RP_C[1][2] = - RP_B[1][2];
RP_C[2][0] = - RP_B[2][0];
RP_C[2][1] = - RP_B[2][1];
RP_C[2][2] = RP_B[2][0] + RP_B[2][1] + RP_B[2][2];
*/

//const int n_CY=6;
interval RP_Cy[6][6];
cinterval RP_Yc_inv[6][6]; 
interval i_zero(0.0);

RP_Cy[0][0] = RP_B[0][0];
RP_Cy[0][1] = RP_B[0][1];
RP_Cy[0][2] = RP_B[0][2];
RP_Cy[0][3] = i_zero;
RP_Cy[0][4] = i_zero;
RP_Cy[0][5] = i_zero;
RP_Cy[1][0] = RP_B[1][0];
RP_Cy[1][1] = RP_B[1][1];
RP_Cy[1][2] = RP_B[1][2];
RP_Cy[1][3] = i_zero;
RP_Cy[1][4] = i_zero;
RP_Cy[1][5] = i_zero;
RP_Cy[2][0] = RP_B[2][0];
RP_Cy[2][1] = RP_B[2][1];
RP_Cy[2][2] = RP_B[2][2];
RP_Cy[2][3] = i_zero;
RP_Cy[2][4] = i_zero; 
RP_Cy[2][5] = i_zero;
RP_Cy[3][0] = i_zero;
RP_Cy[3][1] = i_zero;
RP_Cy[3][2] = i_zero;
RP_Cy[3][3] = RP_B[0][0];
RP_Cy[3][4] = RP_B[0][1];
RP_Cy[3][5] = RP_B[0][2];
RP_Cy[4][0] = i_zero;
RP_Cy[4][1] = i_zero;
RP_Cy[4][2] = i_zero;
RP_Cy[4][3] = RP_B[1][0];
RP_Cy[4][4] = RP_B[1][1];
RP_Cy[4][5] = RP_B[1][2];
RP_Cy[5][0] = i_zero;
RP_Cy[5][1] = i_zero;
RP_Cy[5][2] = i_zero;
RP_Cy[5][3] = RP_B[2][0];
RP_Cy[5][4] = RP_B[2][1];
RP_Cy[5][5] = RP_B[2][2];

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Cy[i][j] = RP_Cy[i][j] * 0.5 * IP_L * IP_cfreq;
	//cout << "Cy" << RP_Cy[i][j] << "\n";
	}
}

cinterval RP_m_temp[6][6];
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(RP_m_temp[i][j]) = i_zero;
	Im(RP_m_temp[i][j]) = i_zero;
	Re(RP_Yc[i][j]) = i_zero;
	Im(RP_Yc[i][j]) = i_zero;
	}
}


for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Im(RP_m_temp[i][j]) = RP_Cy[i][j];
	}
}

double RP_Yc_I_Re, RP_Yc_I_Im, RP_Yc_S_Re, RP_Yc_S_Im; 
double RP_Yrc_I_Re, RP_Yrc_I_Im, RP_Yrc_S_Re, RP_Yrc_S_Im; 
double RP_m_temp_I_Re, RP_m_temp_I_Im, RP_m_temp_S_Re, RP_m_temp_S_Im;

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Yrc_I_Re = _double(Inf(Re(RP_Yrc[i][j])));
	RP_Yrc_I_Im = _double(Inf(Im(RP_Yrc[i][j])));
	RP_Yrc_S_Re = _double(Sup(Re(RP_Yrc[i][j])));
	RP_Yrc_S_Im = _double(Sup(Im(RP_Yrc[i][j])));
	
	RP_m_temp_I_Re = _double(Inf(Re(RP_m_temp[i][j])));
	RP_m_temp_I_Im = _double(Inf(Im(RP_m_temp[i][j])));
	RP_m_temp_S_Re = _double(Sup(Re(RP_m_temp[i][j])));
	RP_m_temp_S_Im = _double(Sup(Im(RP_m_temp[i][j])));

	RP_Yc_I_Re = RP_Yrc_I_Re + RP_m_temp_I_Re;
	RP_Yc_I_Im = RP_Yrc_I_Im + RP_m_temp_I_Im;
	RP_Yc_S_Re = RP_Yrc_S_Re + RP_m_temp_S_Re;
	RP_Yc_S_Im = RP_Yrc_S_Im + RP_m_temp_S_Im;

	Inf(Re(RP_Yc[i][j])) = _real(RP_Yc_I_Re);
	Inf(Im(RP_Yc[i][j])) = _real(RP_Yc_I_Im);
	Sup(Re(RP_Yc[i][j])) = _real(RP_Yc_S_Re);
	Sup(Im(RP_Yc[i][j])) = _real(RP_Yc_S_Im);
	}
}

//cout << "Yrc = " << RP_Yrc[0][0] << "\n";
//cout << "w * Cy = " << RP_m_temp[0][0] << "\n";
//cout << "Yc = " << RP_Yc[0][0] << "\n";

cinterval RP_Yc22[3][3], RP_Yc12[3][3], RP_Yc22_inv[3][3];

RP_Yc22[0][0] = RP_Yc[3][3];
RP_Yc22[0][1] = RP_Yc[3][4];
RP_Yc22[0][2] = RP_Yc[3][5];
RP_Yc22[1][0] = RP_Yc[4][3];
RP_Yc22[1][1] = RP_Yc[4][4];
RP_Yc22[1][2] = RP_Yc[4][5];
RP_Yc22[2][0] = RP_Yc[5][3];
RP_Yc22[2][1] = RP_Yc[5][4];
RP_Yc22[2][2] = RP_Yc[5][5];

RP_Yc12[0][0] = RP_Yc[3][0];
RP_Yc12[0][1] = RP_Yc[3][1];
RP_Yc12[0][2] = RP_Yc[3][2];
RP_Yc12[1][0] = RP_Yc[4][0];
RP_Yc12[1][1] = RP_Yc[4][1];
RP_Yc12[1][2] = RP_Yc[4][2];
RP_Yc12[2][0] = RP_Yc[5][0];
RP_Yc12[2][1] = RP_Yc[5][1];
RP_Yc12[2][2] = RP_Yc[5][2];
//cout << "Yc22 00"<< RP_Yc[3][3];// <-----------[-0.01]

ap::complex_2d_array Yc12_1, Yc12_2, Yc22_1, Yc22_2;
Yc12_1.setlength(3,3);
Yc12_2.setlength(3,3);
Yc22_1.setlength(3,3);
Yc22_2.setlength(3,3);
//cout << "noninv1" << RP_Yc22[0][0] << endl;//<----[-0.01]
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Yc22_1(i,j).x = _double(Inf(Re(RP_Yc22[i][j])));
	Yc22_1(i,j).y = _double(Inf(Im(RP_Yc22[i][j])));
	Yc22_2(i,j).x = _double(Sup(Re(RP_Yc22[i][j])));
	Yc22_2(i,j).y = _double(Sup(Im(RP_Yc22[i][j])));
	}
}
//cout << "noninv" << Yc22_1(0,0).y << endl;//<----[-0.01]
cmatrixinverse(Yc22_1, 3);
cmatrixinverse(Yc22_2, 3);
//!!!!!!!!!!
//cout << "inv" << Yc22_1(0,0).x << endl;
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Inf(Im(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).y,Yc22_2(i,j).y));
	Sup(Re(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Sup(Im(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).y,Yc22_2(i,j).y));
	}
}

/*
cout << "00" << RP_Yc22[0][0] << endl;
cout << "01" << RP_Yc22[0][1] << endl;
cout << "02" << RP_Yc22[0][2] << endl;
cout << "10" << RP_Yc22[1][0] << endl;
cout << "11" << RP_Yc22[1][1] << endl;
cout << "12" << RP_Yc22[1][2] << endl;
cout << "20" << RP_Yc22[2][0] << endl;
cout << "21" << RP_Yc22[2][1] << endl;
cout << "22" << RP_Yc22[2][2] << endl;
*/

cinterval RP_U1[3], RP_U2[3];
interval IP_arg_A, IP_arg_B, IP_arg_C, IP_180;
interval RP1_cos_A, RP1_cos_B, RP1_cos_C, RP1_sin_A, RP1_sin_B, RP1_sin_C;
interval RP2_cos_A, RP2_cos_B, RP2_cos_C, RP2_sin_A, RP2_sin_B, RP2_sin_C;

IP_arg_A = 0.0;
IP_arg_B = - 120.0;
IP_arg_C = 120.0;
//cout << "120" << IP_arg_B << "\n" << endl;
Inf(IP_180) = 180.0;
Sup(IP_180) = 180.0;
RP1_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP1_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP1_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP1_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP1_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP1_sin_C = sin((IP_arg_C / IP_180) * Pi());

RP2_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP2_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP2_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP2_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP2_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP2_sin_C = sin((IP_arg_C / IP_180) * Pi());
/*
interval fd;
Inf(fd) = 0.5;
Sup(fd) = 0.5;
interval fdd = cos(fd);
*/
//cout << "cos" << fdd << "\n" << endl;

//cout << "grad" << RP2_cos_B * Pi() << "\n" << endl;
//cout << "PI" << Pi() << "\n" << endl;
Re(RP_U1[0]) = IP_U_A_f * RP1_cos_A;
Im(RP_U1[0]) = IP_U_A_f * RP1_sin_A;
Re(RP_U1[1]) = IP_U_B_f * RP1_cos_B;
Im(RP_U1[1]) = IP_U_B_f * RP1_sin_B;
Re(RP_U1[2]) = IP_U_C_f * RP1_cos_C;
Im(RP_U1[2]) = IP_U_C_f * RP1_sin_C;

Re(RP_U2[0]) = IP_U_A_f * RP2_cos_A;
Im(RP_U2[0]) = IP_U_A_f * RP2_sin_A;
Re(RP_U2[1]) = IP_U_B_f * RP2_cos_B;
Im(RP_U2[1]) = IP_U_B_f * RP2_sin_B;
Re(RP_U2[2]) = IP_U_C_f * RP2_cos_C;
Im(RP_U2[2]) = IP_U_C_f * RP2_sin_C;

cinterval RP_I2[3], TP_SU[3];
cinterval TP_TMP1[3], TP_TMP2[3], TP_TMP3[3], TP_TMP4[3], TP_TMP5[3];

//--------------------------------------------------------------------------------------
ap::complex_1d_array temp1_1, temp2_1, temp3_1;
ap::complex_1d_array temp1_2, temp2_2, temp3_2;
ap::complex_1d_array U1_1, U1_2, U2_1, U2_2, I2_1, I2_2;
U1_1.setlength(3);
temp1_1.setlength(3);
temp2_1.setlength(3);
temp3_1.setlength(3);
I2_1.setlength(3);

U1_2.setlength(3);
temp1_2.setlength(3);
temp2_2.setlength(3);
temp3_2.setlength(3);
I2_2.setlength(3);

ap::complex_1d_array U456_1, U456_2;
U456_1.setlength(3);
U456_2.setlength(3);
ap::complex I4_1, I5_1, I6_1, U1_1_c, U2_1_c, U3_1_c, U4_1_c, U5_1_c, U6_1_c, Imax;
ap::complex I4_2, I5_2, I6_2, U1_2_c, U2_2_c, U3_2_c, U4_2_c, U5_2_c, U6_2_c;
ap::complex S;
S.x = _double(Inf(Re(IP_S)));
S.y = _double(Inf(Im(IP_S)));

ap::complex_2d_array Yc22_1_inv, Yc22_2_inv;
Yc22_1_inv.setlength(3,3);
Yc22_2_inv.setlength(3,3);

for(int i = 0 ; i < 3; i++)
{
	for(int j = 0 ; j < 3; j++)
	{
	Yc22_1_inv(i,j).x = _double(Inf(Re(RP_Yc22_inv[i][j])));
	Yc22_1_inv(i,j).y = _double(Inf(Im(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).x = _double(Sup(Re(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).y = _double(Sup(Im(RP_Yc22_inv[i][j])));

	Yc12_1(i,j).x = _double(Inf(Re(RP_Yc12[i][j])));
	Yc12_1(i,j).y = _double(Inf(Im(RP_Yc12[i][j])));
	Yc12_2(i,j).x = _double(Sup(Re(RP_Yc12[i][j])));
	Yc12_2(i,j).y = _double(Sup(Im(RP_Yc12[i][j])));
	}
}

U1_1_c.x = _double(Inf(Re(RP_U1[0]))); 
U1_1_c.y = _double(Inf(Im(RP_U1[0]))); 
U2_1_c.x = _double(Inf(Re(RP_U1[1]))); 
U2_1_c.y = _double(Inf(Im(RP_U1[1])));
U3_1_c.x = _double(Inf(Re(RP_U1[2]))); 
U3_1_c.y = _double(Inf(Im(RP_U1[2])));
U1_2_c.x = _double(Sup(Re(RP_U1[0]))); 
U1_2_c.y = _double(Sup(Im(RP_U1[0]))); 
U2_2_c.x = _double(Sup(Re(RP_U1[1]))); 
U2_2_c.y = _double(Sup(Im(RP_U1[1])));
U3_2_c.x = _double(Sup(Re(RP_U1[2]))); 
U3_2_c.y = _double(Sup(Im(RP_U1[2])));

U4_1_c.x = _double(Inf(Re(RP_U2[0]))); 
U4_1_c.y = _double(Inf(Im(RP_U2[0]))); 
U5_1_c.x = _double(Inf(Re(RP_U2[1]))); 
U5_1_c.y = _double(Inf(Im(RP_U2[1])));
U6_1_c.x = _double(Inf(Re(RP_U2[2]))); 
U6_1_c.y = _double(Inf(Im(RP_U2[2])));
U4_2_c.x = _double(Sup(Re(RP_U2[0]))); 
U4_2_c.y = _double(Sup(Im(RP_U2[0]))); 
U5_2_c.x = _double(Sup(Re(RP_U2[1]))); 
U5_2_c.y = _double(Sup(Im(RP_U2[1])));
U6_2_c.x = _double(Sup(Re(RP_U2[2]))); 
U6_2_c.y = _double(Sup(Im(RP_U2[2])));


//cout << "U4 = " << U4_2_c.x << endl << "\n";
// Нижнее значение <-------------------------------------------------------- Run process 1
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

I4_1.x = - (S / U4_1_c).x;
I4_1.y = (S / U4_1_c).y;
I5_1.x = - (S / U5_1_c).x;
I5_1.y = (S / U5_1_c).y;
I6_1.x = - (S / U6_1_c).x;
I6_1.y = (S / U6_1_c).y;

I2_1(0) = I4_1;
I2_1(1) = I5_1;
I2_1(2) = I6_1;

U1_1(0) = U1_1_c;
U1_1(1) = U2_1_c;
U1_1(2) = U3_1_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_1(j).x = 0;
  temp1_1(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yc12_1(j,i) * U1_1(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_1(i) = I2_1(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_1(j) += Yc22_1_inv(j,i) * temp2_1(i);
		}
	}

U4_1_c = temp3_1(0);
U5_1_c = temp3_1(1);
U6_1_c = temp3_1(2);

// Верхнее значение <-------------------------------------------------------- Run process 2

I4_2.x = - (S / U4_2_c).x;
I4_2.y = (S / U4_2_c).y;
I5_2.x = - (S / U5_2_c).x;
I5_2.y = (S / U5_2_c).y;
I6_2.x = - (S / U6_2_c).x;
I6_2.y = (S / U6_2_c).y;

I2_2(0) = I4_2;
I2_2(1) = I5_2;
I2_2(2) = I6_2;

U1_2(0) = U1_2_c;
U1_2(1) = U2_2_c;
U1_2(2) = U3_2_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yc12_2(j,i) * U1_2(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_2(i) = I2_2(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_2(j) += Yc22_2_inv(j,i) * temp2_2(i);
		}
	}

U4_2_c = temp3_2(0);
U5_2_c = temp3_2(1);
U6_2_c = temp3_2(2);



//cout << "Iteration = " << iter_count << endl << "\n";

for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}
/*
for(int i = 0 ; i < 3; i++)
{
cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}
*/
};



for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}


for(int i = 0 ; i < 3; i++)
{
cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
ap::complex_1d_array apRP_U2_low, apRP_U2_high, apRP_SU_low, apRP_SU_high;
apRP_U2_low.setlength(3);
apRP_U2_high.setlength(3);
apRP_SU_low.setlength(3);
apRP_SU_high.setlength(3);

for (int iter_count = 0; iter_count < iter_max; iter_count++)							// run process
{

for (int i = 0; i < 3; i++)
{
apRP_U2_low(i).x = _double(Inf(Re(RP_U2[i])));
apRP_U2_low(i).y = _double(Inf(Im(RP_U2[i])));
apRP_U2_high(i).x = _double(Sup(Re(RP_U2[i])));
apRP_U2_high(i).y = _double(Sup(Im(RP_U2[i])));
}

apRP_SU_low(0) = S / apRP_U2_low(0);
apRP_SU_low(1) = S / apRP_U2_low(1);
apRP_SU_low(2) = S / apRP_U2_low(2);
apRP_SU_high(0) = S / apRP_U2_high(0);
apRP_SU_high(1) = S / apRP_U2_high(1);
apRP_SU_high(2) = S / apRP_U2_high(2);

for (int i = 0; i < 3; i++)
{
Inf(Re(TP_SU[i])) = _real(min(apRP_SU_low(i).x,apRP_SU_high(i).x));
Inf(Im(TP_SU[i])) = _real(min(apRP_SU_low(i).y,apRP_SU_high(i).y));
Sup(Re(TP_SU[i])) = _real(max(apRP_SU_low(i).x,apRP_SU_high(i).x));
Sup(Im(TP_SU[i])) = _real(max(apRP_SU_low(i).y,apRP_SU_high(i).y));
}

	
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_U2_low [" << i << "] = [" << apRP_U2_low(i).x << "]+j*[" << apRP_U2_low(i).y << "]" << endl << "\n";
//cout << "apRP_U2_high [" << i << "] = [" << apRP_U2_high(i).x << "]+j*[" << apRP_U2_high(i).y << "]" << endl << "\n";
//}
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_SU_low [" << i << "] = [" << apRP_SU_low(i).x << "]+j*[" << apRP_SU_low(i).y << "]" << endl << "\n";
//cout << "apRP_SU_high [" << i << "] = [" << apRP_SU_high(i).x << "]+j*[" << apRP_SU_high(i).y << "]" << endl << "\n";
//}



Re(RP_I2[0]) = -1.0 * Re(TP_SU[0]); // ERR
Im(RP_I2[0]) = Im(TP_SU[0]); // OK
Re(RP_I2[1]) = -1.0 * Re(TP_SU[1]); // ERR
Im(RP_I2[1]) = Im(TP_SU[1]);  // OK
Re(RP_I2[2]) = -1.0 * Re(TP_SU[2]);  // OK
Im(RP_I2[2]) = Im(TP_SU[2]); // OK


double TP_TMP1_S_Re, TP_TMP1_S_Im, TP_TMP1_I_Re, TP_TMP1_I_Im;
double TP_TMP2_S_Re, TP_TMP2_S_Im, TP_TMP2_I_Re, TP_TMP2_I_Im;
double TP_TMP3_S_Re, TP_TMP3_S_Im, TP_TMP3_I_Re, TP_TMP3_I_Im;
double TP_TMP_S_Re, TP_TMP_S_Im, TP_TMP_I_Re, TP_TMP_I_Im;

for(int j = 0 ; j < 3; j++)
 {
 Re(TP_TMP2[j]) = 0.0;
 Im(TP_TMP2[j]) = 0.0;
 TP_TMP1[0] = RP_Yc12[j][0] * RP_U1[0];
 TP_TMP1[1] = RP_Yc12[j][1] * RP_U1[1];
 TP_TMP1[2] = RP_Yc12[j][2] * RP_U1[2];
 TP_TMP1_I_Re = _double(Inf(Re(TP_TMP1[0])));
 TP_TMP1_I_Im = _double(Inf(Im(TP_TMP1[0])));
 TP_TMP1_S_Re = _double(Sup(Re(TP_TMP1[0])));
 TP_TMP1_S_Im = _double(Sup(Im(TP_TMP1[0])));

 TP_TMP2_I_Re = _double(Inf(Re(TP_TMP1[1])));
 TP_TMP2_I_Im = _double(Inf(Im(TP_TMP1[1])));
 TP_TMP2_S_Re = _double(Sup(Re(TP_TMP1[1])));
 TP_TMP2_S_Im = _double(Sup(Im(TP_TMP1[1])));

 TP_TMP3_I_Re = _double(Inf(Re(TP_TMP1[2])));
 TP_TMP3_I_Im = _double(Inf(Im(TP_TMP1[2])));
 TP_TMP3_S_Re = _double(Sup(Re(TP_TMP1[2])));
 TP_TMP3_S_Im = _double(Sup(Im(TP_TMP1[2])));

 TP_TMP_I_Re = TP_TMP1_I_Re + TP_TMP2_I_Re + TP_TMP3_I_Re;
 TP_TMP_I_Im = TP_TMP1_I_Im + TP_TMP2_I_Im + TP_TMP3_I_Im;
 TP_TMP_S_Re = TP_TMP1_S_Re + TP_TMP2_S_Re + TP_TMP3_S_Re;
 TP_TMP_S_Im = TP_TMP1_S_Im + TP_TMP2_S_Im + TP_TMP3_S_Im;

 Inf(Re(TP_TMP2[j])) = _real(min(TP_TMP_I_Re,TP_TMP_S_Re));
 Inf(Im(TP_TMP2[j])) = _real(min(TP_TMP_I_Im,TP_TMP_S_Im));
 Sup(Re(TP_TMP2[j])) = _real(max(TP_TMP_I_Re,TP_TMP_S_Re));
 Sup(Im(TP_TMP2[j])) = _real(max(TP_TMP_I_Im,TP_TMP_S_Im));
 }


// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 3; i++)
	{
	TP_TMP3[i] = RP_I2[i] - TP_TMP2[i];
	}	
	
	for(int j = 0 ; j < 3; j++)
	{
	Inf(Re(TP_TMP4[j])) = real(0.0);
	Inf(Im(TP_TMP4[j])) = real(0.0);
	Sup(Re(TP_TMP4[j])) = real(0.0);
	Sup(Im(TP_TMP4[j])) = real(0.0);
	for(int i = 0 ; i < 3; i++)
	{
	TP_TMP4[j] += RP_Yc22_inv[j][i] * TP_TMP3[i];
	}
	}
	

RP_U2[0] = TP_TMP4[0];
RP_U2[1] = TP_TMP4[1];
RP_U2[2] = TP_TMP4[2];


//////////////////////////////////////////////////////////////////////////////////////////////
cout << "Iteration = " << iter_count << endl << "\n";

U4_mod_min = _double(abs(Inf(TP_TMP4[0])));
U4_mod_max = _double(abs(Sup(TP_TMP4[0])));
U5_mod_min = _double(abs(Inf(TP_TMP4[1])));
U5_mod_max = _double(abs(Sup(TP_TMP4[1])));
U6_mod_min = _double(abs(Inf(TP_TMP4[2])));
U6_mod_max = _double(abs(Sup(TP_TMP4[2])));

cout << "U4 = " << TP_TMP4[0] << endl << "\n";
cout << "U5 = " << TP_TMP4[1] << endl << "\n";
cout << "U6 = " << TP_TMP4[2] << endl << "\n";


cout << "mod U4 = [ " << U4_mod_min << " ; " << U4_mod_max << " ]" << endl << "\n";
cout << "mod U5 = [ " << U5_mod_min << " ; " << U5_mod_max << " ]" << endl << "\n";
cout << "mod U6 = [ " << U6_mod_min << " ; " << U6_mod_max << " ]" << endl << "\n";
//////////////////////////////////////////////////////////////////////////////////////////////


}															// stop process
*/
}

//////////////////////////////////////////////////////////////////////////////////////////////
void getintvloss_triangle()
{

//  Расчет 1-ой схемы

ap::complex_2d_array Yb_1_shm, Yb_2_shm, M_shm, P_shm, Y_1_shm, Y_2_shm;
ap::complex One, Zero;
Yb_1_shm.setlength(18,18);
Yb_2_shm.setlength(18,18);
Y_1_shm.setlength(9,9);
Y_2_shm.setlength(9,9);
M_shm.setlength(9,18);
One.x = 1;
One.y = 0;
Zero.x = 0;
Zero.y = 0;
// i = число строк, j - число столбцов
// Yb
for (int i=0; i < 18; i++)
{
	for (int j=0; j < 18; j++)
	{
	Yb_1_shm(i,j) = Zero;
	Yb_2_shm(i,j) = Zero;
	}
}
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_1(i,j);
	Yb_2_shm(i,j) = Yc2_temp_1(i,j);
	}
}
for (int i=6; i < 12; i++)
{
	for (int j=6; j < 12; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_2(i-6,j-6);
	Yb_2_shm(i,j) = Yc2_temp_2(i-6,j-6);
	}
}
for (int i=12; i < 18; i++)
{
	for (int j=12; j < 18; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_3(i-12,j-12);
	Yb_2_shm(i,j) = Yc2_temp_3(i-12,j-12);
	}
}

P_shm.setlength(3,6);
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	P_shm(i,j) = Zero;
	}
}
P_shm(0,0) = One;
P_shm(1,1) = One;
P_shm(2,2) = One;

// i = число строк, j - число столбцов
// M
for (int i=0; i < 9; i++)
{
	for (int j=0; j < 18; j++)
	{
	M_shm(i,j) = Zero;
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = - P_shm(i,j);
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = - P_shm(i,j-6);
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = P_shm(i-3,j-6);
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = - P_shm(i-3,j-12);
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = P_shm(i-6,j);
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = P_shm(i-6,j-12);
	}
}

ap::complex_2d_array Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm_temp.setlength(9,18);
Y_2_shm_temp.setlength(9,18);
Y_1_shm_temp_reint.setlength(9,18);
Y_2_shm_temp_reint.setlength(9,18);

// M * Yb
for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 18; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		 Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}

for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}
// Transpose M

ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
M_shm_transposed.setlength(18,9);

for (int i=0; i < 9; i++)
	for (int j=0; j < 18 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}

// Y_temp * M^transpose
for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 9; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}

ap::complex_2d_array Yf1_1_shm, Yf1_2_shm, Yf1_1_shm_inv, Yf1_2_shm_inv, Yfb_1_shm, Yfb_2_shm;
Yf1_1_shm.setlength(6,6);
Yf1_2_shm.setlength(6,6);
Yf1_1_shm_inv.setlength(6,6);
Yf1_2_shm_inv.setlength(6,6);
Yfb_1_shm.setlength(6,3);
Yfb_2_shm.setlength(6,3);

for (int i = 0; i < 6; i++)
{
	for (int j = 0; j < 6; j++)
	{
	Yf1_1_shm(i,j) = Y_1_shm(i,j);
	Yf1_2_shm(i,j) = Y_2_shm(i,j);
	}

}

Yfb_1_shm(0,0) = Y_1_shm(0,6);
Yfb_1_shm(0,1) = Y_1_shm(0,7);
Yfb_1_shm(0,2) = Y_1_shm(0,8);
Yfb_1_shm(1,0) = Y_1_shm(1,6);
Yfb_1_shm(1,1) = Y_1_shm(1,7);
Yfb_1_shm(1,2) = Y_1_shm(1,8);
Yfb_1_shm(2,0) = Y_1_shm(2,6);
Yfb_1_shm(2,1) = Y_1_shm(2,7);
Yfb_1_shm(2,2) = Y_1_shm(2,8);
Yfb_1_shm(3,0) = Y_1_shm(3,6);
Yfb_1_shm(3,1) = Y_1_shm(3,7);
Yfb_1_shm(3,2) = Y_1_shm(3,8);
Yfb_1_shm(4,0) = Y_1_shm(4,6);
Yfb_1_shm(4,1) = Y_1_shm(4,7);
Yfb_1_shm(4,2) = Y_1_shm(4,8);
Yfb_1_shm(5,0) = Y_1_shm(5,6);
Yfb_1_shm(5,1) = Y_1_shm(5,7);
Yfb_1_shm(5,2) = Y_1_shm(5,8);

Yfb_2_shm(0,0) = Y_2_shm(0,6);
Yfb_2_shm(0,1) = Y_2_shm(0,7);
Yfb_2_shm(0,2) = Y_2_shm(0,8);
Yfb_2_shm(1,0) = Y_2_shm(1,6);
Yfb_2_shm(1,1) = Y_2_shm(1,7);
Yfb_2_shm(1,2) = Y_2_shm(1,8);
Yfb_2_shm(2,0) = Y_2_shm(2,6);
Yfb_2_shm(2,1) = Y_2_shm(2,7);
Yfb_2_shm(2,2) = Y_2_shm(2,8);
Yfb_2_shm(3,0) = Y_2_shm(3,6);
Yfb_2_shm(3,1) = Y_2_shm(3,7);
Yfb_2_shm(3,2) = Y_2_shm(3,8);
Yfb_2_shm(4,0) = Y_2_shm(4,6);
Yfb_2_shm(4,1) = Y_2_shm(4,7);
Yfb_2_shm(4,2) = Y_2_shm(4,8);
Yfb_2_shm(5,0) = Y_2_shm(5,6);
Yfb_2_shm(5,1) = Y_2_shm(5,7);
Yfb_2_shm(5,2) = Y_2_shm(5,8);

cmatrixinverse(Yf1_1_shm, 6);
Yf1_1_shm_inv = Yf1_1_shm;
cmatrixinverse(Yf1_2_shm, 6);
Yf1_2_shm_inv = Yf1_2_shm;

// Решение системы
ap::complex  Imax;
ap::complex_1d_array U_1_shm, U_2_shm, I_1_shm, I_2_shm, Ubal_1_shm, Ubal_2_shm;
U_1_shm.setlength(6);
U_2_shm.setlength(6);
I_1_shm.setlength(6);
I_2_shm.setlength(6);
Ubal_1_shm.setlength(3);
Ubal_2_shm.setlength(3);

ap::complex S1, S2, coef_S;
temp1_1.setlength(6);
temp2_1.setlength(6);
temp3_1.setlength(6);
temp1_2.setlength(6);
temp2_2.setlength(6);
temp3_2.setlength(6);

double IP_Ua_f = (IP_Ua_l * 1000.0) / sqrt(3.0);         ////////////////////////////////////////////////////////////////////
double IP_Ub_f = (IP_Ub_l * 1000.0) / sqrt(3.0); 
double IP_Uc_f = (IP_Uc_l * 1000.0) / sqrt(3.0); 

// Напряжение в начале ЛЭП
// Нижнее значение
// Ua1
Ubal_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Верхнее значение
// Ua1
Ubal_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);

// Начальное приближение
// Нижнее значение
// Ua1
U_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Верхнее значение
// Ua1
U_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);

// Нижнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1.x = IP_P1 * 1000000.0;
S1.y = IP_Q1 * 1000000.0;
S2.x = IP_P2 * 1000000.0;
S2.y = IP_Q2 * 1000000.0;

I_1_shm(0).x = - (S1 / U_1_shm(0)).x;
I_1_shm(0).y = (S1 / U_1_shm(0)).y;
I_1_shm(1).x = - (S1 / U_1_shm(1)).x;
I_1_shm(1).y = (S1 / U_1_shm(1)).y;
I_1_shm(2).x = - (S1 / U_1_shm(2)).x;
I_1_shm(2).y = (S1 / U_1_shm(2)).y;
I_1_shm(3).x = - (S2 / U_1_shm(3)).x;
I_1_shm(3).y = (S2 / U_1_shm(3)).y;
I_1_shm(4).x = - (S2 / U_1_shm(4)).x;
I_1_shm(4).y = (S2 / U_1_shm(4)).y;
I_1_shm(5).x = - (S2 / U_1_shm(5)).x;
I_1_shm(5).y = (S2 / U_1_shm(5)).y;

for(int j = 0 ; j < 6; j++)
{
temp1_1(j).x = 0;
temp1_1(j).y = 0;
}
for(int j = 0 ; j < 6; j++)
 {
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yfb_1_shm(j,i) * Ubal_1_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 6; i++)
	{
	temp2_1(i) = I_1_shm(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 6; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 6; i++)
		{
		temp3_1(j) += Yf1_1_shm_inv(j,i) * temp2_1(i);
		}
	}
for (int i = 0; i < 6; i++)
{
U_1_shm(i) = temp3_1(i);
}
};

// Верхнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{


S1.x = IP_P1 * 1000000.0;
S1.y = IP_Q1 * 1000000.0;
S2.x = IP_P2 * 1000000.0;
S2.y = IP_Q2 * 1000000.0;

I_2_shm(0).x = - (S1 / U_2_shm(0)).x;
I_2_shm(0).y = (S1 / U_2_shm(0)).y;
I_2_shm(1).x = - (S1 / U_2_shm(1)).x;
I_2_shm(1).y = (S1 / U_2_shm(1)).y;
I_2_shm(2).x = - (S1 / U_2_shm(2)).x;
I_2_shm(2).y = (S1 / U_2_shm(2)).y;
I_2_shm(3).x = - (S2 / U_2_shm(3)).x;
I_2_shm(3).y = (S2 / U_2_shm(3)).y;
I_2_shm(4).x = - (S2 / U_2_shm(4)).x;
I_2_shm(4).y = (S2 / U_2_shm(4)).y;
I_2_shm(5).x = - (S2 / U_2_shm(5)).x;
I_2_shm(5).y = (S2 / U_2_shm(5)).y;

for(int j = 0 ; j < 6; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yfb_2_shm(j,i) * Ubal_2_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 6; i++)
	{
	temp2_2(i) = I_2_shm(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 6; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 6; i++)
		{
		temp3_2(j) += Yf1_2_shm_inv(j,i) * temp2_2(i);
		}
	}

for (int i = 0; i < 6; i++)
{
U_2_shm(i) = temp3_2(i);
}

};

// Итоговый отчет
// Нижнее значение
ap::complex U_a1_1 = temp3_1(0) / 1000.0;
ap::complex U_b1_1 = temp3_1(1) / 1000.0;
ap::complex U_c1_1 = temp3_1(2) / 1000.0;
ap::complex U_a2_1 = temp3_1(3) / 1000.0;
ap::complex U_b2_1 = temp3_1(4) / 1000.0;
ap::complex U_c2_1 = temp3_1(5) / 1000.0;
U_a1_1_mod = sqrt(pow((temp3_1(0).x), 2.0)+pow((temp3_1(0).y), 2.0)) / 1000.0;
U_b1_1_mod = sqrt(pow((temp3_1(1).x), 2.0)+pow((temp3_1(1).y), 2.0)) / 1000.0;
U_c1_1_mod = sqrt(pow((temp3_1(2).x), 2.0)+pow((temp3_1(2).y), 2.0)) / 1000.0;
U_a2_1_mod = sqrt(pow((temp3_1(3).x), 2.0)+pow((temp3_1(3).y), 2.0)) / 1000.0;
U_b2_1_mod = sqrt(pow((temp3_1(4).x), 2.0)+pow((temp3_1(4).y), 2.0)) / 1000.0;
U_c2_1_mod = sqrt(pow((temp3_1(5).x), 2.0)+pow((temp3_1(5).y), 2.0)) / 1000.0;


// Выбор узла
if (triangle_ptype == 0)
{
U4f_1_mod = U_a1_1_mod;
U5f_1_mod = U_b1_1_mod;
U6f_1_mod = U_c1_1_mod;
}
if (triangle_ptype == 1)
{
U4f_1_mod = U_a2_1_mod;
U5f_1_mod = U_b2_1_mod;
U6f_1_mod = U_c2_1_mod;
}

// Верхнее значение
ap::complex U_a1_2 = temp3_2(0) / 1000.0;
ap::complex U_b1_2 = temp3_2(1) / 1000.0;
ap::complex U_c1_2 = temp3_2(2) / 1000.0;
ap::complex U_a2_2 = temp3_2(3) / 1000.0;
ap::complex U_b2_2 = temp3_2(4) / 1000.0;
ap::complex U_c2_2 = temp3_2(5) / 1000.0;
U_a1_2_mod = sqrt(pow((temp3_2(0).x), 2.0)+pow((temp3_2(0).y), 2.0)) / 1000.0;
U_b1_2_mod = sqrt(pow((temp3_2(1).x), 2.0)+pow((temp3_2(1).y), 2.0)) / 1000.0;
U_c1_2_mod = sqrt(pow((temp3_2(2).x), 2.0)+pow((temp3_2(2).y), 2.0)) / 1000.0;
U_a2_2_mod = sqrt(pow((temp3_2(3).x), 2.0)+pow((temp3_2(3).y), 2.0)) / 1000.0;
U_b2_2_mod = sqrt(pow((temp3_2(4).x), 2.0)+pow((temp3_2(4).y), 2.0)) / 1000.0;
U_c2_2_mod = sqrt(pow((temp3_2(5).x), 2.0)+pow((temp3_2(5).y), 2.0)) / 1000.0;

cout << "mod U_a1 = [ " << U_a1_1_mod << " , " << U_a1_2_mod << " ] кВ\n"; 
cout << "mod U_b1 = [ " << U_b1_1_mod << " , " << U_b1_2_mod << " ] кВ\n"; 
cout << "mod U_c1 = [ " << U_c1_1_mod << " , " << U_c1_2_mod << " ] кВ\n"; 
cout << "mod U_a2 = [ " << U_a2_1_mod << " , " << U_a2_2_mod << " ] кВ\n"; 
cout << "mod U_b2 = [ " << U_b2_1_mod << " , " << U_b2_2_mod << " ] кВ\n"; 
cout << "mod U_c2 = [ " << U_c2_1_mod << " , " << U_c2_2_mod << " ] кВ\n";

for (int i=0; i < 6; i++)
{
MF_U_low(i) = temp3_1(i);
MF_U_high(i) = temp3_2(i);
}

}


void setintvloss_triangle()
{
Yc1_temp_1.setlength(6,6);
Yc1_temp_2.setlength(6,6);
Yc1_temp_3.setlength(6,6);
Yc2_temp_1.setlength(6,6);
Yc2_temp_2.setlength(6,6);
Yc2_temp_3.setlength(6,6);

triangle_ptype = 0;
IP_P1 = 4.0;
IP_Q1 = 3.0;
IP_P2 = 12.0;
IP_Q2 = 9.0;
cout << "S1 = " << IP_P1 << " + j*" << IP_Q1 << endl <<"\n";
cout << "S2 = " << IP_P2 << " + j*" << IP_Q2 << endl <<"\n";
Length = 45.0;
cout << "L = " << Length << " ; [" << IP_L << "] км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_1(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_1(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_1(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_1(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 10.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	MF_ZM_low(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	MF_ZM_low(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	MF_ZM_high(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	MF_ZM_high(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	}
}

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_2(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_2(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_2(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_2(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 15.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_3(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_3(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_3(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_3(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

getintvloss_triangle();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void getintvloss_star()
{

//  Расчет 2-ой схемы

ap::complex_2d_array Yb_1_shm, Yb_2_shm, M_shm, P_shm, Y_1_shm, Y_2_shm;
ap::complex One, Zero;
Yb_1_shm.setlength(18,18);
Yb_2_shm.setlength(18,18);
Y_1_shm.setlength(12,12);
Y_2_shm.setlength(12,12);
M_shm.setlength(12,18);
One.x = 1;
One.y = 0;
Zero.x = 0;
Zero.y = 0;
// i = число строк, j - число столбцов
// Yb
for (int i=0; i < 18; i++)
{
	for (int j=0; j < 18; j++)
	{
	Yb_1_shm(i,j) = Zero;
	Yb_2_shm(i,j) = Zero;
	}
}
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_1(i,j);
	Yb_2_shm(i,j) = Yc2_temp_1(i,j);
	}
}
for (int i=6; i < 12; i++)
{
	for (int j=6; j < 12; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_2(i-6,j-6);
	Yb_2_shm(i,j) = Yc2_temp_2(i-6,j-6);
	}
}
for (int i=12; i < 18; i++)
{
	for (int j=12; j < 18; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_3(i-12,j-12);
	Yb_2_shm(i,j) = Yc2_temp_3(i-12,j-12);
	}
}

P_shm.setlength(3,6);
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	P_shm(i,j) = Zero;
	}
}
P_shm(0,0) = One;
P_shm(1,1) = One;
P_shm(2,2) = One;

// i = число строк, j - число столбцов
// M
for (int i=0; i < 12; i++)
{
	for (int j=0; j < 18; j++)
	{
	M_shm(i,j) = Zero;
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = - P_shm(i,j-6);
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = - P_shm(i-3,j-12);
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = - P_shm(i-6,j);
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = P_shm(i-6,j-6);
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = P_shm(i-6,j-12);
	}
}
for (int i=9; i < 12; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = P_shm(i-9,j);
	}
}

ap::complex_2d_array Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm_temp.setlength(12,18);
Y_2_shm_temp.setlength(12,18);
Y_1_shm_temp_reint.setlength(12,18);
Y_2_shm_temp_reint.setlength(12,18);

// M * Yb
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 18; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		 Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}

for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}

// Transpose M

ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
M_shm_transposed.setlength(18,12);

for (int i=0; i < 12; i++)
	for (int j=0; j < 18 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}

// Y_temp * M^transpose
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}

ap::complex_2d_array Yf1_1_shm, Yf1_2_shm, Yf1_1_shm_inv, Yf1_2_shm_inv, Yfb_1_shm, Yfb_2_shm;
Yf1_1_shm.setlength(9,9);
Yf1_2_shm.setlength(9,9);
Yf1_1_shm_inv.setlength(9,9);
Yf1_2_shm_inv.setlength(9,9);
Yfb_1_shm.setlength(9,3);
Yfb_2_shm.setlength(9,3);

for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 9; j++)
	{
	Yf1_1_shm(i,j) = Y_1_shm(i,j);
	Yf1_2_shm(i,j) = Y_2_shm(i,j);
	}

}

Yfb_1_shm(0,0) = Y_1_shm(0,9);
Yfb_1_shm(0,1) = Y_1_shm(0,10);
Yfb_1_shm(0,2) = Y_1_shm(0,11);
Yfb_1_shm(1,0) = Y_1_shm(1,9);
Yfb_1_shm(1,1) = Y_1_shm(1,10);
Yfb_1_shm(1,2) = Y_1_shm(1,11);
Yfb_1_shm(2,0) = Y_1_shm(2,9);
Yfb_1_shm(2,1) = Y_1_shm(2,10);
Yfb_1_shm(2,2) = Y_1_shm(2,11);
Yfb_1_shm(3,0) = Y_1_shm(3,9);
Yfb_1_shm(3,1) = Y_1_shm(3,10);
Yfb_1_shm(3,2) = Y_1_shm(3,11);
Yfb_1_shm(4,0) = Y_1_shm(4,9);
Yfb_1_shm(4,1) = Y_1_shm(4,10);
Yfb_1_shm(4,2) = Y_1_shm(4,11);
Yfb_1_shm(5,0) = Y_1_shm(5,9);
Yfb_1_shm(5,1) = Y_1_shm(5,10);
Yfb_1_shm(5,2) = Y_1_shm(5,11);
Yfb_1_shm(6,0) = Y_1_shm(6,9);
Yfb_1_shm(6,1) = Y_1_shm(6,10);
Yfb_1_shm(6,2) = Y_1_shm(6,11);
Yfb_1_shm(7,0) = Y_1_shm(7,9);
Yfb_1_shm(7,1) = Y_1_shm(7,10);
Yfb_1_shm(7,2) = Y_1_shm(7,11);
Yfb_1_shm(8,0) = Y_1_shm(8,9);
Yfb_1_shm(8,1) = Y_1_shm(8,10);
Yfb_1_shm(8,2) = Y_1_shm(8,11);

Yfb_2_shm(0,0) = Y_2_shm(0,9);
Yfb_2_shm(0,1) = Y_2_shm(0,10);
Yfb_2_shm(0,2) = Y_2_shm(0,11);
Yfb_2_shm(1,0) = Y_2_shm(1,9);
Yfb_2_shm(1,1) = Y_2_shm(1,10);
Yfb_2_shm(1,2) = Y_2_shm(1,11);
Yfb_2_shm(2,0) = Y_2_shm(2,9);
Yfb_2_shm(2,1) = Y_2_shm(2,10);
Yfb_2_shm(2,2) = Y_2_shm(2,11);
Yfb_2_shm(3,0) = Y_2_shm(3,9);
Yfb_2_shm(3,1) = Y_2_shm(3,10);
Yfb_2_shm(3,2) = Y_2_shm(3,11);
Yfb_2_shm(4,0) = Y_2_shm(4,9);
Yfb_2_shm(4,1) = Y_2_shm(4,10);
Yfb_2_shm(4,2) = Y_2_shm(4,11);
Yfb_2_shm(5,0) = Y_2_shm(5,9);
Yfb_2_shm(5,1) = Y_2_shm(5,10);
Yfb_2_shm(5,2) = Y_2_shm(5,11);
Yfb_2_shm(6,0) = Y_2_shm(6,9);
Yfb_2_shm(6,1) = Y_2_shm(6,10);
Yfb_2_shm(6,2) = Y_2_shm(6,11);
Yfb_2_shm(7,0) = Y_2_shm(7,9);
Yfb_2_shm(7,1) = Y_2_shm(7,10);
Yfb_2_shm(7,2) = Y_2_shm(7,11);
Yfb_2_shm(8,0) = Y_2_shm(8,9);
Yfb_2_shm(8,1) = Y_2_shm(8,10);
Yfb_2_shm(8,2) = Y_2_shm(8,11);

cmatrixinverse(Yf1_1_shm, 9);
Yf1_1_shm_inv = Yf1_1_shm;
cmatrixinverse(Yf1_2_shm, 9);
Yf1_2_shm_inv = Yf1_2_shm;

// Решение системы
ap::complex  Imax;
ap::complex_1d_array U_1_shm, U_2_shm, I_1_shm, I_2_shm, Ubal_1_shm, Ubal_2_shm;
U_1_shm.setlength(9);
U_2_shm.setlength(9);
I_1_shm.setlength(9);
I_2_shm.setlength(9);
Ubal_1_shm.setlength(3);
Ubal_2_shm.setlength(3);

ap::complex S1, S2, coef_S;
temp1_1.setlength(9);
temp2_1.setlength(9);
temp3_1.setlength(9);
temp1_2.setlength(9);
temp2_2.setlength(9);
temp3_2.setlength(9);

double IP_Ua_f = (IP_Ua_l * 1000.0) / sqrt(3.0);         ////////////////////////////////////////////////////////////////////
double IP_Ub_f = (IP_Ub_l * 1000.0) / sqrt(3.0); 
double IP_Uc_f = (IP_Uc_l * 1000.0) / sqrt(3.0); 

// Напряжение в начале ЛЭП
// Нижнее значение
// Ua1
Ubal_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Верхнее значение
// Ua1
Ubal_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);

// Начальное приближение
// Нижнее значение
// Ua1
U_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Верхнее значение
// Ua1
U_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Нижнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1.x = IP_P1 * 1000000.0;
S1.y = IP_Q1 * 1000000.0;
S2.x = IP_P2 * 1000000.0;
S2.y = IP_Q2 * 1000000.0;

I_1_shm(0).x = - (S1 / U_1_shm(0)).x;
I_1_shm(0).y = (S1 / U_1_shm(0)).y;
I_1_shm(1).x = - (S1 / U_1_shm(1)).x;
I_1_shm(1).y = (S1 / U_1_shm(1)).y;
I_1_shm(2).x = - (S1 / U_1_shm(2)).x;
I_1_shm(2).y = (S1 / U_1_shm(2)).y;
I_1_shm(3).x = - (S2 / U_1_shm(3)).x;
I_1_shm(3).y = (S2 / U_1_shm(3)).y;
I_1_shm(4).x = - (S2 / U_1_shm(4)).x;
I_1_shm(4).y = (S2 / U_1_shm(4)).y;
I_1_shm(5).x = - (S2 / U_1_shm(5)).x;
I_1_shm(5).y = (S2 / U_1_shm(5)).y;
I_1_shm(6).x = - (0 / U_1_shm(6)).x;
I_1_shm(6).y = (0 / U_1_shm(6)).y;
I_1_shm(7).x = - (0 / U_1_shm(7)).x;
I_1_shm(7).y = (0 / U_1_shm(7)).y;
I_1_shm(8).x = - (0 / U_1_shm(8)).x;
I_1_shm(8).y = (0 / U_1_shm(8)).y;

for(int j = 0 ; j < 9; j++)
{
temp1_1(j).x = 0;
temp1_1(j).y = 0;
}
for(int j = 0 ; j < 9; j++)
 {
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yfb_1_shm(j,i) * Ubal_1_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 9; i++)
	{
	temp2_1(i) = I_1_shm(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_1(j) += Yf1_1_shm_inv(j,i) * temp2_1(i);
		}
	}
for (int i = 0; i < 9; i++)
{
U_1_shm(i) = temp3_1(i);
}
};

// Верхнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{


S1.x = IP_P1 * 1000000.0;
S1.y = IP_Q1 * 1000000.0;
S2.x = IP_P2 * 1000000.0;
S2.y = IP_Q2 * 1000000.0;

I_2_shm(0).x = - (S1 / U_2_shm(0)).x;
I_2_shm(0).y = (S1 / U_2_shm(0)).y;
I_2_shm(1).x = - (S1 / U_2_shm(1)).x;
I_2_shm(1).y = (S1 / U_2_shm(1)).y;
I_2_shm(2).x = - (S1 / U_2_shm(2)).x;
I_2_shm(2).y = (S1 / U_2_shm(2)).y;
I_2_shm(3).x = - (S2 / U_2_shm(3)).x;
I_2_shm(3).y = (S2 / U_2_shm(3)).y;
I_2_shm(4).x = - (S2 / U_2_shm(4)).x;
I_2_shm(4).y = (S2 / U_2_shm(4)).y;
I_2_shm(5).x = - (S2 / U_2_shm(5)).x;
I_2_shm(5).y = (S2 / U_2_shm(5)).y;
I_2_shm(6).x = - (0 / U_2_shm(6)).x;
I_2_shm(6).y = (0 / U_2_shm(6)).y;
I_2_shm(7).x = - (0 / U_2_shm(7)).x;
I_2_shm(7).y = (0 / U_2_shm(7)).y;
I_2_shm(8).x = - (0 / U_2_shm(8)).x;
I_2_shm(8).y = (0 / U_2_shm(8)).y;

for(int j = 0 ; j < 9; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yfb_2_shm(j,i) * Ubal_2_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 9; i++)
	{
	temp2_2(i) = I_2_shm(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_2(j) += Yf1_2_shm_inv(j,i) * temp2_2(i);
		}
	}

for (int i = 0; i < 9; i++)
{
U_2_shm(i) = temp3_2(i);
}

};

// Итоговый отчет
// Нижнее значение
ap::complex U_a1_1 = temp3_1(0) / 1000.0;
ap::complex U_b1_1 = temp3_1(1) / 1000.0;
ap::complex U_c1_1 = temp3_1(2) / 1000.0;
ap::complex U_a2_1 = temp3_1(3) / 1000.0;
ap::complex U_b2_1 = temp3_1(4) / 1000.0;
ap::complex U_c2_1 = temp3_1(5) / 1000.0;
ap::complex U_a3_1 = temp3_1(6) / 1000.0;
ap::complex U_b3_1 = temp3_1(7) / 1000.0;
ap::complex U_c3_1 = temp3_1(8) / 1000.0;

U_a1_1_mod = sqrt(pow((temp3_1(0).x), 2.0)+pow((temp3_1(0).y), 2.0)) / 1000.0;
U_b1_1_mod = sqrt(pow((temp3_1(1).x), 2.0)+pow((temp3_1(1).y), 2.0)) / 1000.0;
U_c1_1_mod = sqrt(pow((temp3_1(2).x), 2.0)+pow((temp3_1(2).y), 2.0)) / 1000.0;
U_a2_1_mod = sqrt(pow((temp3_1(3).x), 2.0)+pow((temp3_1(3).y), 2.0)) / 1000.0;
U_b2_1_mod = sqrt(pow((temp3_1(4).x), 2.0)+pow((temp3_1(4).y), 2.0)) / 1000.0;
U_c2_1_mod = sqrt(pow((temp3_1(5).x), 2.0)+pow((temp3_1(5).y), 2.0)) / 1000.0;
U_a3_1_mod = sqrt(pow((temp3_1(6).x), 2.0)+pow((temp3_1(6).y), 2.0)) / 1000.0;
U_b3_1_mod = sqrt(pow((temp3_1(7).x), 2.0)+pow((temp3_1(7).y), 2.0)) / 1000.0;
U_c3_1_mod = sqrt(pow((temp3_1(8).x), 2.0)+pow((temp3_1(8).y), 2.0)) / 1000.0;

// Выбор узла
if (triangle_ptype == 0)
{
U4f_1_mod = U_a1_1_mod;
U5f_1_mod = U_b1_1_mod;
U6f_1_mod = U_c1_1_mod;
}
if (triangle_ptype == 1)
{
U4f_1_mod = U_a2_1_mod;
U5f_1_mod = U_b2_1_mod;
U6f_1_mod = U_c2_1_mod;
}

// Верхнее значение
ap::complex U_a1_2 = temp3_2(0) / 1000.0;
ap::complex U_b1_2 = temp3_2(1) / 1000.0;
ap::complex U_c1_2 = temp3_2(2) / 1000.0;
ap::complex U_a2_2 = temp3_2(3) / 1000.0;
ap::complex U_b2_2 = temp3_2(4) / 1000.0;
ap::complex U_c2_2 = temp3_2(5) / 1000.0;
ap::complex U_a3_2 = temp3_2(6) / 1000.0;
ap::complex U_b3_2 = temp3_2(7) / 1000.0;
ap::complex U_c3_2 = temp3_2(8) / 1000.0;

U_a1_2_mod = sqrt(pow((temp3_2(0).x), 2.0)+pow((temp3_2(0).y), 2.0)) / 1000.0;
U_b1_2_mod = sqrt(pow((temp3_2(1).x), 2.0)+pow((temp3_2(1).y), 2.0)) / 1000.0;
U_c1_2_mod = sqrt(pow((temp3_2(2).x), 2.0)+pow((temp3_2(2).y), 2.0)) / 1000.0;
U_a2_2_mod = sqrt(pow((temp3_2(3).x), 2.0)+pow((temp3_2(3).y), 2.0)) / 1000.0;
U_b2_2_mod = sqrt(pow((temp3_2(4).x), 2.0)+pow((temp3_2(4).y), 2.0)) / 1000.0;
U_c2_2_mod = sqrt(pow((temp3_2(5).x), 2.0)+pow((temp3_2(5).y), 2.0)) / 1000.0;
U_a3_2_mod = sqrt(pow((temp3_2(6).x), 2.0)+pow((temp3_2(6).y), 2.0)) / 1000.0;
U_b3_2_mod = sqrt(pow((temp3_2(7).x), 2.0)+pow((temp3_2(7).y), 2.0)) / 1000.0;
U_c3_2_mod = sqrt(pow((temp3_2(8).x), 2.0)+pow((temp3_2(8).y), 2.0)) / 1000.0;

cout << "mod U_a1 = [ " << U_a1_1_mod << " , " << U_a1_2_mod << " ] кВ\n"; 
cout << "mod U_b1 = [ " << U_b1_1_mod << " , " << U_b1_2_mod << " ] кВ\n"; 
cout << "mod U_c1 = [ " << U_c1_1_mod << " , " << U_c1_2_mod << " ] кВ\n"; 
cout << "mod U_a2 = [ " << U_a2_1_mod << " , " << U_a2_2_mod << " ] кВ\n"; 
cout << "mod U_b2 = [ " << U_b2_1_mod << " , " << U_b2_2_mod << " ] кВ\n"; 
cout << "mod U_c2 = [ " << U_c2_1_mod << " , " << U_c2_2_mod << " ] кВ\n"; 
cout << "mod U_a3 = [ " << U_a3_1_mod << " , " << U_a3_2_mod << " ] кВ\n"; 
cout << "mod U_b3 = [ " << U_b3_1_mod << " , " << U_b3_2_mod << " ] кВ\n"; 
cout << "mod U_c3 = [ " << U_c3_1_mod << " , " << U_c3_2_mod << " ] кВ\n"; 

for (int i=0; i < 9; i++)
{
MF_U_low(i) = temp3_1(i);
MF_U_high(i) = temp3_2(i);
}


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setintvloss_star()
{
Yc1_temp_1.setlength(6,6);
Yc1_temp_2.setlength(6,6);
Yc1_temp_3.setlength(6,6);
Yc2_temp_1.setlength(6,6);
Yc2_temp_2.setlength(6,6);
Yc2_temp_3.setlength(6,6);

triangle_ptype = 0;
IP_P1 = 12.0;
IP_Q1 = 9.0;
IP_P2 = 15.0;
IP_Q2 = 11.25;
//cout << "S1 = " << IP_P1 << " + j*" << IP_Q1 << endl <<"\n";
//cout << "S2 = " << IP_P2 << " + j*" << IP_Q2 << endl <<"\n";
Length = 30.0;
cout << "L = " << Length << " ; [" << IP_L << "] км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_1(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_1(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_1(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_1(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 25.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	MF_ZM_low(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	MF_ZM_low(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	MF_ZM_high(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	MF_ZM_high(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	}
}

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_2(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_2(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_2(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_2(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 20.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_3(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_3(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_3(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_3(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

getintvloss_star();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void getintvalue_00()
{
interval IP_GIV00_Y, IP_GIV00_A, IP_GIV00_B, IP_GIV00_C, IP_GIV00_X;
double IP_GPV00_X;
Inf(IP_GIV00_A) = 2.0;
Sup(IP_GIV00_A) = 3.0;
Inf(IP_GIV00_B) = 4.0;
Sup(IP_GIV00_B) = 5.0;
Inf(IP_GIV00_C) = 6.0;
Sup(IP_GIV00_C) = 8.0;
Inf(IP_GIV00_X) = 3.0;
Sup(IP_GIV00_X) = 3.0;
IP_GPV00_X = 0.0;
for (int i=0; i < 4; i++)
{
Inf(IP_GIV00_X) = IP_GPV00_X;
Sup(IP_GIV00_X) = IP_GPV00_X;
IP_GPV00_X += 1.0;
IP_GIV00_Y = (IP_GIV00_A * IP_GIV00_B * IP_GIV00_X) / IP_GIV00_C;
cout << "A = " << IP_GIV00_A << endl << "\n";
cout << "B = " << IP_GIV00_B << endl << "\n";
cout << "C = " << IP_GIV00_C << endl << "\n";
cout << "X = " << IP_GIV00_X << endl << "\n";
cout << "Y = (A * B * X) / C = " << IP_GIV00_Y << endl << "\n";
}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void getintvloss_3line()
{

//  Расчет 3-ей схемы

ap::complex_2d_array Yb_1_shm, Yb_2_shm, M_shm, P_shm, Y_1_shm, Y_2_shm;
ap::complex One, Zero;
Yb_1_shm.setlength(18,18);
Yb_2_shm.setlength(18,18);
Y_1_shm.setlength(12,12);
Y_2_shm.setlength(12,12);
M_shm.setlength(12,18);
One.x = 1;
One.y = 0;
Zero.x = 0;
Zero.y = 0;
// i = число строк, j - число столбцов
// Yb
for (int i=0; i < 18; i++)
{
	for (int j=0; j < 18; j++)
	{
	Yb_1_shm(i,j) = Zero;
	Yb_2_shm(i,j) = Zero;
	}
}
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_1(i,j);
	Yb_2_shm(i,j) = Yc2_temp_1(i,j);
	}
}
for (int i=6; i < 12; i++)
{
	for (int j=6; j < 12; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_2(i-6,j-6);
	Yb_2_shm(i,j) = Yc2_temp_2(i-6,j-6);
	}
}
for (int i=12; i < 18; i++)
{
	for (int j=12; j < 18; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_3(i-12,j-12);
	Yb_2_shm(i,j) = Yc2_temp_3(i-12,j-12);
	}
}

P_shm.setlength(3,6);
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	P_shm(i,j) = Zero;
	}
}
P_shm(0,0) = One;
P_shm(1,1) = One;
P_shm(2,2) = One;

// i = число строк, j - число столбцов
// M
for (int i=0; i < 12; i++)
{
	for (int j=0; j < 18; j++)
	{
	M_shm(i,j) = Zero;
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = - P_shm(i,j);  //*
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = P_shm(i,j-6);  //*
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = - P_shm(i-3,j-6); //*
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = P_shm(i-3,j-12); //*
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = - P_shm(i-6,j-12);
	}
}
for (int i=9; i < 12; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = P_shm(i-9,j);
	}
}

ap::complex_2d_array Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm_temp.setlength(12,18);
Y_2_shm_temp.setlength(12,18);
Y_1_shm_temp_reint.setlength(12,18);
Y_2_shm_temp_reint.setlength(12,18);

// M * Yb
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 18; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		 Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}

for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 18; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}

// Transpose M

ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
M_shm_transposed.setlength(18,12);

for (int i=0; i < 12; i++)
	for (int j=0; j < 18 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}

// Y_temp * M^transpose
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 18; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}
//cout << "Y1(0,0) = " << Y_1_shm(9,6).y << "\n";
ap::complex_2d_array Yf1_1_shm, Yf1_2_shm, Yf1_1_shm_inv, Yf1_2_shm_inv, Yfb_1_shm, Yfb_2_shm;
Yf1_1_shm.setlength(9,9);
Yf1_2_shm.setlength(9,9);
Yf1_1_shm_inv.setlength(9,9);
Yf1_2_shm_inv.setlength(9,9);
Yfb_1_shm.setlength(9,3);
Yfb_2_shm.setlength(9,3);

for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 9; j++)
	{
	Yf1_1_shm(i,j) = Y_1_shm(i,j);
	Yf1_2_shm(i,j) = Y_2_shm(i,j);
	}

}

Yfb_1_shm(0,0) = Y_1_shm(0,9);
Yfb_1_shm(0,1) = Y_1_shm(0,10);
Yfb_1_shm(0,2) = Y_1_shm(0,11);
Yfb_1_shm(1,0) = Y_1_shm(1,9);
Yfb_1_shm(1,1) = Y_1_shm(1,10);
Yfb_1_shm(1,2) = Y_1_shm(1,11);
Yfb_1_shm(2,0) = Y_1_shm(2,9);
Yfb_1_shm(2,1) = Y_1_shm(2,10);
Yfb_1_shm(2,2) = Y_1_shm(2,11);
Yfb_1_shm(3,0) = Y_1_shm(3,9);
Yfb_1_shm(3,1) = Y_1_shm(3,10);
Yfb_1_shm(3,2) = Y_1_shm(3,11);
Yfb_1_shm(4,0) = Y_1_shm(4,9);
Yfb_1_shm(4,1) = Y_1_shm(4,10);
Yfb_1_shm(4,2) = Y_1_shm(4,11);
Yfb_1_shm(5,0) = Y_1_shm(5,9);
Yfb_1_shm(5,1) = Y_1_shm(5,10);
Yfb_1_shm(5,2) = Y_1_shm(5,11);
Yfb_1_shm(6,0) = Y_1_shm(6,9);
Yfb_1_shm(6,1) = Y_1_shm(6,10);
Yfb_1_shm(6,2) = Y_1_shm(6,11);
Yfb_1_shm(7,0) = Y_1_shm(7,9);
Yfb_1_shm(7,1) = Y_1_shm(7,10);
Yfb_1_shm(7,2) = Y_1_shm(7,11);
Yfb_1_shm(8,0) = Y_1_shm(8,9);
Yfb_1_shm(8,1) = Y_1_shm(8,10);
Yfb_1_shm(8,2) = Y_1_shm(8,11);

Yfb_2_shm(0,0) = Y_2_shm(0,9);
Yfb_2_shm(0,1) = Y_2_shm(0,10);
Yfb_2_shm(0,2) = Y_2_shm(0,11);
Yfb_2_shm(1,0) = Y_2_shm(1,9);
Yfb_2_shm(1,1) = Y_2_shm(1,10);
Yfb_2_shm(1,2) = Y_2_shm(1,11);
Yfb_2_shm(2,0) = Y_2_shm(2,9);
Yfb_2_shm(2,1) = Y_2_shm(2,10);
Yfb_2_shm(2,2) = Y_2_shm(2,11);
Yfb_2_shm(3,0) = Y_2_shm(3,9);
Yfb_2_shm(3,1) = Y_2_shm(3,10);
Yfb_2_shm(3,2) = Y_2_shm(3,11);
Yfb_2_shm(4,0) = Y_2_shm(4,9);
Yfb_2_shm(4,1) = Y_2_shm(4,10);
Yfb_2_shm(4,2) = Y_2_shm(4,11);
Yfb_2_shm(5,0) = Y_2_shm(5,9);
Yfb_2_shm(5,1) = Y_2_shm(5,10);
Yfb_2_shm(5,2) = Y_2_shm(5,11);
Yfb_2_shm(6,0) = Y_2_shm(6,9);
Yfb_2_shm(6,1) = Y_2_shm(6,10);
Yfb_2_shm(6,2) = Y_2_shm(6,11);
Yfb_2_shm(7,0) = Y_2_shm(7,9);
Yfb_2_shm(7,1) = Y_2_shm(7,10);
Yfb_2_shm(7,2) = Y_2_shm(7,11);
Yfb_2_shm(8,0) = Y_2_shm(8,9);
Yfb_2_shm(8,1) = Y_2_shm(8,10);
Yfb_2_shm(8,2) = Y_2_shm(8,11);

cmatrixinverse(Yf1_1_shm, 9);
Yf1_1_shm_inv = Yf1_1_shm;
cmatrixinverse(Yf1_2_shm, 9);
Yf1_2_shm_inv = Yf1_2_shm;

// Решение системы
ap::complex  Imax;
ap::complex_1d_array U_1_shm, U_2_shm, I_1_shm, I_2_shm, Ubal_1_shm, Ubal_2_shm;
U_1_shm.setlength(9);
U_2_shm.setlength(9);
I_1_shm.setlength(9);
I_2_shm.setlength(9);
Ubal_1_shm.setlength(3);
Ubal_2_shm.setlength(3);

ap::complex S1_A, S1_B, S1_C, S2_A, S2_B, S2_C, S3_A, S3_B, S3_C, coef_S;
temp1_1.setlength(9);
temp2_1.setlength(9);
temp3_1.setlength(9);
temp1_2.setlength(9);
temp2_2.setlength(9);
temp3_2.setlength(9);

/*
double IP_Ua_f = (IP_Ua_l * 1000.0) / sqrt(3.0);       
double IP_Ub_f = (IP_Ub_l * 1000.0) / sqrt(3.0); 
double IP_Uc_f = (IP_Uc_l * 1000.0) / sqrt(3.0); 
*/
double IP_Ua_f = 133000.0;       
double IP_Ub_f = 133000.0;  
double IP_Uc_f = 133000.0;  

// Напряжение в начале ЛЭП
// Нижнее значение
// Ua1
Ubal_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Верхнее значение
// Ua1
Ubal_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);

// Начальное приближение
// Нижнее значение
// Ua1
U_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Верхнее значение
// Ua1
U_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Нижнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1_A.x = IP_P1_A_low * 1000000.0;
S1_A.y = IP_Q1_A_low * 1000000.0;
S1_B.x = IP_P1_B_low * 1000000.0;
S1_B.y = IP_Q1_B_low * 1000000.0;
S1_C.x = IP_P1_C_low * 1000000.0;
S1_C.y = IP_Q1_C_low * 1000000.0;

S2_A.x = IP_P2_A_low * 1000000.0;
S2_A.y = IP_Q2_A_low * 1000000.0;
S2_B.x = IP_P2_B_low * 1000000.0;
S2_B.y = IP_Q2_B_low * 1000000.0;
S2_C.x = IP_P2_C_low * 1000000.0;
S2_C.y = IP_Q2_C_low * 1000000.0;

S3_A.x = IP_P3_A_low * 1000000.0;
S3_A.y = IP_Q3_A_low * 1000000.0;
S3_B.x = IP_P3_B_low * 1000000.0;
S3_B.y = IP_Q3_B_low * 1000000.0;
S3_C.x = IP_P3_C_low * 1000000.0;
S3_C.y = IP_Q3_C_low * 1000000.0;

I_1_shm(0).x = - (S1_A / U_1_shm(0)).x;
I_1_shm(0).y = (S1_A / U_1_shm(0)).y;
I_1_shm(1).x = - (S1_B / U_1_shm(1)).x;
I_1_shm(1).y = (S1_B / U_1_shm(1)).y;
I_1_shm(2).x = - (S1_C / U_1_shm(2)).x;
I_1_shm(2).y = (S1_C / U_1_shm(2)).y;
I_1_shm(3).x = - (S2_A / U_1_shm(3)).x;
I_1_shm(3).y = (S2_A / U_1_shm(3)).y;
I_1_shm(4).x = - (S2_B / U_1_shm(4)).x;
I_1_shm(4).y = (S2_B / U_1_shm(4)).y;
I_1_shm(5).x = - (S2_C / U_1_shm(5)).x;
I_1_shm(5).y = (S2_C / U_1_shm(5)).y;
I_1_shm(6).x = - (S3_A / U_1_shm(6)).x;
I_1_shm(6).y = (S3_A / U_1_shm(6)).y;
I_1_shm(7).x = - (S3_B / U_1_shm(7)).x;
I_1_shm(7).y = (S3_B / U_1_shm(7)).y;
I_1_shm(8).x = - (S3_C / U_1_shm(8)).x;
I_1_shm(8).y = (S3_C / U_1_shm(8)).y;

for(int j = 0 ; j < 9; j++)
{
temp1_1(j).x = 0;
temp1_1(j).y = 0;
}
for(int j = 0 ; j < 9; j++)
 {
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yfb_1_shm(j,i) * Ubal_1_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 9; i++)
	{
	temp2_1(i) = I_1_shm(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_1(j) += Yf1_1_shm_inv(j,i) * temp2_1(i);
		}
	}
for (int i = 0; i < 9; i++)
{
U_1_shm(i) = temp3_1(i);
}
};

// Верхнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1_A.x = IP_P1_A_high * 1000000.0;
S1_A.y = IP_Q1_A_high * 1000000.0;
S1_B.x = IP_P1_B_high * 1000000.0;
S1_B.y = IP_Q1_B_high * 1000000.0;
S1_C.x = IP_P1_C_high * 1000000.0;
S1_C.y = IP_Q1_C_high * 1000000.0;

S2_A.x = IP_P2_A_high * 1000000.0;
S2_A.y = IP_Q2_A_high * 1000000.0;
S2_B.x = IP_P2_B_high * 1000000.0;
S2_B.y = IP_Q2_B_high * 1000000.0;
S2_C.x = IP_P2_C_high * 1000000.0;
S2_C.y = IP_Q2_C_high * 1000000.0;

S3_A.x = IP_P3_A_high * 1000000.0;
S3_A.y = IP_Q3_A_high * 1000000.0;
S3_B.x = IP_P3_B_high * 1000000.0;
S3_B.y = IP_Q3_B_high * 1000000.0;
S3_C.x = IP_P3_C_high * 1000000.0;
S3_C.y = IP_Q3_C_high * 1000000.0;

I_2_shm(0).x = - (S1_A / U_2_shm(0)).x;
I_2_shm(0).y = (S1_A / U_2_shm(0)).y;
I_2_shm(1).x = - (S1_B / U_2_shm(1)).x;
I_2_shm(1).y = (S1_B / U_2_shm(1)).y;
I_2_shm(2).x = - (S1_C / U_2_shm(2)).x;
I_2_shm(2).y = (S1_C / U_2_shm(2)).y;
I_2_shm(3).x = - (S2_A / U_2_shm(3)).x;
I_2_shm(3).y = (S2_A / U_2_shm(3)).y;
I_2_shm(4).x = - (S2_B / U_2_shm(4)).x;
I_2_shm(4).y = (S2_B / U_2_shm(4)).y;
I_2_shm(5).x = - (S2_C / U_2_shm(5)).x;
I_2_shm(5).y = (S2_C / U_2_shm(5)).y;
I_2_shm(6).x = - (S3_A / U_2_shm(6)).x;
I_2_shm(6).y = (S3_A / U_2_shm(6)).y;
I_2_shm(7).x = - (S3_B / U_2_shm(7)).x;
I_2_shm(7).y = (S3_B / U_2_shm(7)).y;
I_2_shm(8).x = - (S3_C / U_2_shm(8)).x;
I_2_shm(8).y = (S3_C / U_2_shm(8)).y;

for(int j = 0 ; j < 9; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yfb_2_shm(j,i) * Ubal_2_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 9; i++)
	{
	temp2_2(i) = I_2_shm(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_2(j) += Yf1_2_shm_inv(j,i) * temp2_2(i);
		}
	}

for (int i = 0; i < 9; i++)
{
U_2_shm(i) = temp3_2(i);
}

};




// Итоговый отчет
// Нижнее значение
ap::complex U_a1_1 = temp3_1(0) / 1000.0;
ap::complex U_b1_1 = temp3_1(1) / 1000.0;
ap::complex U_c1_1 = temp3_1(2) / 1000.0;
ap::complex U_a2_1 = temp3_1(3) / 1000.0;
ap::complex U_b2_1 = temp3_1(4) / 1000.0;
ap::complex U_c2_1 = temp3_1(5) / 1000.0;
ap::complex U_a3_1 = temp3_1(6) / 1000.0;
ap::complex U_b3_1 = temp3_1(7) / 1000.0;
ap::complex U_c3_1 = temp3_1(8) / 1000.0;

U_a1_1_mod = sqrt(pow((temp3_1(0).x), 2.0)+pow((temp3_1(0).y), 2.0)) / 1000.0;
U_b1_1_mod = sqrt(pow((temp3_1(1).x), 2.0)+pow((temp3_1(1).y), 2.0)) / 1000.0;
U_c1_1_mod = sqrt(pow((temp3_1(2).x), 2.0)+pow((temp3_1(2).y), 2.0)) / 1000.0;
U_a2_1_mod = sqrt(pow((temp3_1(3).x), 2.0)+pow((temp3_1(3).y), 2.0)) / 1000.0;
U_b2_1_mod = sqrt(pow((temp3_1(4).x), 2.0)+pow((temp3_1(4).y), 2.0)) / 1000.0;
U_c2_1_mod = sqrt(pow((temp3_1(5).x), 2.0)+pow((temp3_1(5).y), 2.0)) / 1000.0;
U_a3_1_mod = sqrt(pow((temp3_1(6).x), 2.0)+pow((temp3_1(6).y), 2.0)) / 1000.0;
U_b3_1_mod = sqrt(pow((temp3_1(7).x), 2.0)+pow((temp3_1(7).y), 2.0)) / 1000.0;
U_c3_1_mod = sqrt(pow((temp3_1(8).x), 2.0)+pow((temp3_1(8).y), 2.0)) / 1000.0;

// Выбор узла
if (triangle_ptype == 0)
{
U4f_1_mod = U_a1_1_mod;
U5f_1_mod = U_b1_1_mod;
U6f_1_mod = U_c1_1_mod;
}
if (triangle_ptype == 1)
{
U4f_1_mod = U_a2_1_mod;
U5f_1_mod = U_b2_1_mod;
U6f_1_mod = U_c2_1_mod;
}

// Верхнее значение
ap::complex U_a1_2 = temp3_2(0) / 1000.0;
ap::complex U_b1_2 = temp3_2(1) / 1000.0;
ap::complex U_c1_2 = temp3_2(2) / 1000.0;
ap::complex U_a2_2 = temp3_2(3) / 1000.0;
ap::complex U_b2_2 = temp3_2(4) / 1000.0;
ap::complex U_c2_2 = temp3_2(5) / 1000.0;
ap::complex U_a3_2 = temp3_2(6) / 1000.0;
ap::complex U_b3_2 = temp3_2(7) / 1000.0;
ap::complex U_c3_2 = temp3_2(8) / 1000.0;

U_a1_2_mod = sqrt(pow((temp3_2(0).x), 2.0)+pow((temp3_2(0).y), 2.0)) / 1000.0;
U_b1_2_mod = sqrt(pow((temp3_2(1).x), 2.0)+pow((temp3_2(1).y), 2.0)) / 1000.0;
U_c1_2_mod = sqrt(pow((temp3_2(2).x), 2.0)+pow((temp3_2(2).y), 2.0)) / 1000.0;
U_a2_2_mod = sqrt(pow((temp3_2(3).x), 2.0)+pow((temp3_2(3).y), 2.0)) / 1000.0;
U_b2_2_mod = sqrt(pow((temp3_2(4).x), 2.0)+pow((temp3_2(4).y), 2.0)) / 1000.0;
U_c2_2_mod = sqrt(pow((temp3_2(5).x), 2.0)+pow((temp3_2(5).y), 2.0)) / 1000.0;
U_a3_2_mod = sqrt(pow((temp3_2(6).x), 2.0)+pow((temp3_2(6).y), 2.0)) / 1000.0;
U_b3_2_mod = sqrt(pow((temp3_2(7).x), 2.0)+pow((temp3_2(7).y), 2.0)) / 1000.0;
U_c3_2_mod = sqrt(pow((temp3_2(8).x), 2.0)+pow((temp3_2(8).y), 2.0)) / 1000.0;


cout << "mod U_a1 = [ " << U_a1_1_mod << " , " << U_a1_2_mod << " ] кВ\n"; 
cout << "mod U_b1 = [ " << U_b1_1_mod << " , " << U_b1_2_mod << " ] кВ\n"; 
cout << "mod U_c1 = [ " << U_c1_1_mod << " , " << U_c1_2_mod << " ] кВ\n"; 
cout << "mod U_a2 = [ " << U_a2_1_mod << " , " << U_a2_2_mod << " ] кВ\n"; 
cout << "mod U_b2 = [ " << U_b2_1_mod << " , " << U_b2_2_mod << " ] кВ\n"; 
cout << "mod U_c2 = [ " << U_c2_1_mod << " , " << U_c2_2_mod << " ] кВ\n"; 
cout << "mod U_a3 = [ " << U_a3_1_mod << " , " << U_a3_2_mod << " ] кВ\n"; 
cout << "mod U_b3 = [ " << U_b3_1_mod << " , " << U_b3_2_mod << " ] кВ\n"; 
cout << "mod U_c3 = [ " << U_c3_1_mod << " , " << U_c3_2_mod << " ] кВ\n"; 

TEMP_U_1.setlength(3,1);
TEMP_U_2.setlength(3,1);
/*
for(int i = 0 ; i < 3; i++)
{
	for(int j = 0 ; j < 1; j++)
	{
	TEMP_U_1(i,j) = temp3_1(i);
	TEMP_U_2(i,j) = temp3_2(i);
	}
}
*/

for(int i = 3 ; i < 6; i++)
{
	for(int j = 0 ; j < 1; j++)
	{
	TEMP_U_1(i-3,j) = temp3_1(i);
	TEMP_U_2(i-3,j) = temp3_2(i);
	}
}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void setintvloss_3line()
{
cout << "Расчет участка сети.\n";
Yc1_temp_1.setlength(6,6);
Yc1_temp_2.setlength(6,6);
Yc1_temp_3.setlength(6,6);
Yc2_temp_1.setlength(6,6);
Yc2_temp_2.setlength(6,6);
Yc2_temp_3.setlength(6,6);

triangle_ptype = 0;
/*
IP_P1_A = 10.0;
IP_Q1_A = 7.5;
IP_P1_B = 10.0;
IP_Q1_B = 7.5;
IP_P1_C = 10.0;
IP_Q1_C = 7.5;

IP_P2_A = 15.0;
IP_Q2_A = 11.25;
IP_P2_B = 15.0;
IP_Q2_B = 11.25;
IP_P2_C = 15.0;
IP_Q2_C = 11.25;
*/
IP_P3_A = 0.0;
IP_Q3_A = 0.0;
IP_P3_B = 0.0;
IP_Q3_B = 0.0;
IP_P3_C = 0.0;
IP_Q3_C = 0.0;

/////////////////////////////
IP_P1_A_low = 9.0;
IP_P1_A_high = 11.0;
IP_Q1_A_low = 6.75;
IP_Q1_A_high = 8.25;

IP_P1_B_low = 7.2;
IP_P1_B_high = 8.8;
IP_Q1_B_low = 5.4;
IP_Q1_B_high = 6.6;

IP_P1_C_low = 10.8;
IP_P1_C_high = 13.2;
IP_Q1_C_low = 8.1;
IP_Q1_C_high = 9.9;

/////////////////////////////
IP_P2_A_low = 13.5;
IP_P2_A_high = 16.5;
IP_Q2_A_low = 10.125;
IP_Q2_A_high = 11.275;

IP_P2_B_low = 10.8;
IP_P2_B_high = 13.2;
IP_Q2_B_low = 8.1;
IP_Q2_B_high = 9.9;

IP_P2_C_low = 16.2;
IP_P2_C_high = 19.8;
IP_Q2_C_low = 12.15;
IP_Q2_C_high = 14.85;

/////////////////////////////
IP_P3_A_low = IP_P3_A;
IP_P3_A_high = IP_P3_A;
IP_Q3_A_low = IP_Q3_A;
IP_Q3_A_high = IP_Q3_A;

IP_P3_B_high = IP_P3_B;
IP_P3_B_low = IP_P3_B;
IP_Q3_B_high = IP_Q3_B;
IP_Q3_B_low = IP_Q3_B;

IP_P3_C_high = IP_P3_C;
IP_P3_C_low = IP_P3_C;
IP_Q3_C_high = IP_Q3_C;
IP_Q3_C_low = IP_Q3_C;
/////////////////////////////

get_PQ_5();

//cout << "S1 = " << IP_P1 << " + j*" << IP_Q1 << endl <<"\n";
//cout << "S2 = " << IP_P2 << " + j*" << IP_Q2 << endl <<"\n";
Length = 40.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_1(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_1(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_1(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_1(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 60.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_2(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_2(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_2(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_2(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 0.001;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_3(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_3(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_3(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_3(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

getintvloss_3line();
get_K2U();
cout << "Process is OK" << "\n";
}
//-----------------------------------------------------------------------------------------------------------
void setintvloss_3line_fault()
{
cout << "Расчет участка сети.\n";
Yc1_temp_1.setlength(6,6);
Yc1_temp_2.setlength(6,6);
Yc1_temp_3.setlength(6,6);
Yc2_temp_1.setlength(6,6);
Yc2_temp_2.setlength(6,6);
Yc2_temp_3.setlength(6,6);

triangle_ptype = 0;

IP_P1_A = 0.0;
IP_Q1_A = 0.0;
IP_P1_B = 0.0;
IP_Q1_B = 0.0;
IP_P1_C = 0.0;
IP_Q1_C = 0.0;

IP_P2_A = 0.0;
IP_Q2_A = 0.0;
IP_P2_B = 0.0;
IP_Q2_B = 0.0;
IP_P2_C = 0.0;
IP_Q2_C = 0.0;

IP_P3_A = 0.0;
IP_Q3_A = 0.0;
IP_P3_B = 36.0;
IP_Q3_B = 27.0;
IP_P3_C = 36.0;
IP_Q3_C = 27.0;

Length = 50.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_1(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_1(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_1(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_1(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 0.001;
cout << "L = " << Length << " км" << endl << "\n";

phase_fault = 1;
getintvloss();
phase_fault = 0;
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_2(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_2(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_2(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_2(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 45.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_3(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_3(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_3(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_3(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

getintvloss_3line();

cout << "Process is OK" << "\n";
}
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//====================================================================================================
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void getintvloss_triangle_4point()
{
// Process triangle fault

//  Расчет 3-ей схемы

ap::complex_2d_array Yb_1_shm, Yb_2_shm, M_shm, P_shm, Y_1_shm, Y_2_shm;
ap::complex One, Zero;
Yb_1_shm.setlength(24,24);
Yb_2_shm.setlength(24,24);
Y_1_shm.setlength(12,12);
Y_2_shm.setlength(12,12);
M_shm.setlength(12,24);
One.x = 1;
One.y = 0;
Zero.x = 0;
Zero.y = 0;
// i = число строк, j - число столбцов
// Yb
for (int i=0; i < 24; i++)
{
	for (int j=0; j < 24; j++)
	{
	Yb_1_shm(i,j) = Zero;
	Yb_2_shm(i,j) = Zero;
	}
}
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_1(i,j);
	Yb_2_shm(i,j) = Yc2_temp_1(i,j);
	}
}
for (int i=6; i < 12; i++)
{
	for (int j=6; j < 12; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_2(i-6,j-6);
	Yb_2_shm(i,j) = Yc2_temp_2(i-6,j-6);
	}
}
for (int i=12; i < 18; i++)
{
	for (int j=12; j < 18; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_3(i-12,j-12);
	Yb_2_shm(i,j) = Yc2_temp_3(i-12,j-12);
	}
}
for (int i=18; i < 24; i++)
{
	for (int j=18; j < 24; j++)
	{
	Yb_1_shm(i,j) = Yc1_temp_4(i-18,j-18);
	Yb_2_shm(i,j) = Yc2_temp_4(i-18,j-18);
	}
}
P_shm.setlength(3,6);
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	P_shm(i,j) = Zero;
	}
}
P_shm(0,0) = One;
P_shm(1,1) = One;
P_shm(2,2) = One;

// i = число строк, j - число столбцов
// M
for (int i=0; i < 12; i++)
{
	for (int j=0; j < 24; j++)
	{
	M_shm(i,j) = Zero;
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = - P_shm(i,j);  //ok
	}
}
for (int i=0; i < 3; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = - P_shm(i,j-6);  //ok
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=6; j < 12; j++)
	{
	M_shm(i,j) = P_shm(i-3,j-6); //ok
	}
}
for (int i=3; i < 6; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = - P_shm(i-3,j-12); //ok
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=12; j < 18; j++)
	{
	M_shm(i,j) = P_shm(i-6,j-12); //ok
	}
}
for (int i=6; i < 9; i++)
{
	for (int j=18; j < 24; j++)
	{
	M_shm(i,j) = - P_shm(i-6,j-18);//ok
	}
}
for (int i=9; i < 12; i++)
{
	for (int j=0; j < 6; j++)
	{
	M_shm(i,j) = P_shm(i-9,j);// ok
	}
}
for (int i=9; i < 12; i++)
{
	for (int j=18; j < 24; j++)
	{
	M_shm(i,j) = P_shm(i-9,j-18);// ok
	}
}

ap::complex_2d_array Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm_temp.setlength(12,24);
Y_2_shm_temp.setlength(12,24);
Y_1_shm_temp_reint.setlength(12,24);
Y_2_shm_temp_reint.setlength(12,24);

// M * Yb
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 24; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 24; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		 Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}

for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 24; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 12; i++)
{
	for (int j = 0; j < 24; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}

// Transpose M

ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
M_shm_transposed.setlength(24,12);

for (int i=0; i < 12; i++)
	for (int j=0; j < 24 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}

// Y_temp * M^transpose
for (int i = 0; i < 12; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 24; k++)  // Кол-во столбцов 1 матрицы
	  {
		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}
//cout << "Y1(0,0) = " << Y_1_shm(9,6).x << "\n";
//cout << "Y1(0,0) = " << Y_1_shm(9,6).y << "\n";
ap::complex_2d_array Yf1_1_shm, Yf1_2_shm, Yf1_1_shm_inv, Yf1_2_shm_inv, Yfb_1_shm, Yfb_2_shm;
Yf1_1_shm.setlength(9,9);
Yf1_2_shm.setlength(9,9);
Yf1_1_shm_inv.setlength(9,9);
Yf1_2_shm_inv.setlength(9,9);
Yfb_1_shm.setlength(9,3);
Yfb_2_shm.setlength(9,3);

for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 9; j++)
	{
	Yf1_1_shm(i,j) = Y_1_shm(i,j);
	Yf1_2_shm(i,j) = Y_2_shm(i,j);
	}

}

Yfb_1_shm(0,0) = Y_1_shm(0,9);
Yfb_1_shm(0,1) = Y_1_shm(0,10);
Yfb_1_shm(0,2) = Y_1_shm(0,11);
Yfb_1_shm(1,0) = Y_1_shm(1,9);
Yfb_1_shm(1,1) = Y_1_shm(1,10);
Yfb_1_shm(1,2) = Y_1_shm(1,11);
Yfb_1_shm(2,0) = Y_1_shm(2,9);
Yfb_1_shm(2,1) = Y_1_shm(2,10);
Yfb_1_shm(2,2) = Y_1_shm(2,11);
Yfb_1_shm(3,0) = Y_1_shm(3,9);
Yfb_1_shm(3,1) = Y_1_shm(3,10);
Yfb_1_shm(3,2) = Y_1_shm(3,11);
Yfb_1_shm(4,0) = Y_1_shm(4,9);
Yfb_1_shm(4,1) = Y_1_shm(4,10);
Yfb_1_shm(4,2) = Y_1_shm(4,11);
Yfb_1_shm(5,0) = Y_1_shm(5,9);
Yfb_1_shm(5,1) = Y_1_shm(5,10);
Yfb_1_shm(5,2) = Y_1_shm(5,11);
Yfb_1_shm(6,0) = Y_1_shm(6,9);
Yfb_1_shm(6,1) = Y_1_shm(6,10);
Yfb_1_shm(6,2) = Y_1_shm(6,11);
Yfb_1_shm(7,0) = Y_1_shm(7,9);
Yfb_1_shm(7,1) = Y_1_shm(7,10);
Yfb_1_shm(7,2) = Y_1_shm(7,11);
Yfb_1_shm(8,0) = Y_1_shm(8,9);
Yfb_1_shm(8,1) = Y_1_shm(8,10);
Yfb_1_shm(8,2) = Y_1_shm(8,11);

Yfb_2_shm(0,0) = Y_2_shm(0,9);
Yfb_2_shm(0,1) = Y_2_shm(0,10);
Yfb_2_shm(0,2) = Y_2_shm(0,11);
Yfb_2_shm(1,0) = Y_2_shm(1,9);
Yfb_2_shm(1,1) = Y_2_shm(1,10);
Yfb_2_shm(1,2) = Y_2_shm(1,11);
Yfb_2_shm(2,0) = Y_2_shm(2,9);
Yfb_2_shm(2,1) = Y_2_shm(2,10);
Yfb_2_shm(2,2) = Y_2_shm(2,11);
Yfb_2_shm(3,0) = Y_2_shm(3,9);
Yfb_2_shm(3,1) = Y_2_shm(3,10);
Yfb_2_shm(3,2) = Y_2_shm(3,11);
Yfb_2_shm(4,0) = Y_2_shm(4,9);
Yfb_2_shm(4,1) = Y_2_shm(4,10);
Yfb_2_shm(4,2) = Y_2_shm(4,11);
Yfb_2_shm(5,0) = Y_2_shm(5,9);
Yfb_2_shm(5,1) = Y_2_shm(5,10);
Yfb_2_shm(5,2) = Y_2_shm(5,11);
Yfb_2_shm(6,0) = Y_2_shm(6,9);
Yfb_2_shm(6,1) = Y_2_shm(6,10);
Yfb_2_shm(6,2) = Y_2_shm(6,11);
Yfb_2_shm(7,0) = Y_2_shm(7,9);
Yfb_2_shm(7,1) = Y_2_shm(7,10);
Yfb_2_shm(7,2) = Y_2_shm(7,11);
Yfb_2_shm(8,0) = Y_2_shm(8,9);
Yfb_2_shm(8,1) = Y_2_shm(8,10);
Yfb_2_shm(8,2) = Y_2_shm(8,11);

cmatrixinverse(Yf1_1_shm, 9);
Yf1_1_shm_inv = Yf1_1_shm;
cmatrixinverse(Yf1_2_shm, 9);
Yf1_2_shm_inv = Yf1_2_shm;

// Решение системы
ap::complex  Imax;
ap::complex_1d_array U_1_shm, U_2_shm, I_1_shm, I_2_shm, Ubal_1_shm, Ubal_2_shm;
U_1_shm.setlength(9);
U_2_shm.setlength(9);
I_1_shm.setlength(9);
I_2_shm.setlength(9);
Ubal_1_shm.setlength(3);
Ubal_2_shm.setlength(3);

ap::complex S1_A, S1_B, S1_C, S2_A, S2_B, S2_C, S3_A, S3_B, S3_C, S4_A, S4_B, S4_C, coef_S;
temp1_1.setlength(9);
temp2_1.setlength(9);
temp3_1.setlength(9);
temp1_2.setlength(9);
temp2_2.setlength(9);
temp3_2.setlength(9);


//double IP_Ua_f = (IP_Ua_l * 1000.0) / sqrt(3.0);       
//double IP_Ub_f = (IP_Ub_l * 1000.0) / sqrt(3.0); 
//double IP_Uc_f = (IP_Uc_l * 1000.0) / sqrt(3.0); 

double IP_Ua_f = 133000.0;       
double IP_Ub_f = 133000.0;  
double IP_Uc_f = 133000.0;  

// Напряжение в начале ЛЭП
// Нижнее значение
// Ua1
Ubal_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Верхнее значение
// Ua1
Ubal_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
Ubal_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
Ubal_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
Ubal_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
Ubal_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
Ubal_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);

// Начальное приближение
// Нижнее значение
// Ua1
U_1_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_1_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_1_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_1_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_1_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_1_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_1_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Верхнее значение
// Ua1
U_2_shm(0).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(0).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(1).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(1).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(2).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(2).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(3).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(3).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(4).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(4).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(5).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(5).y = IP_Uc_f * sin((120.0 / 180.0) * pi);
// Ua1
U_2_shm(6).x = IP_Ua_f * cos((0.0 / 180.0) * pi);
U_2_shm(6).y = IP_Ua_f * sin((0.0 / 180.0) * pi);
// Ub1
U_2_shm(7).x = IP_Ub_f * cos((-120.0 / 180.0) * pi);
U_2_shm(7).y = IP_Ub_f * sin((-120.0 / 180.0) * pi);
// Uc1
U_2_shm(8).x = IP_Uc_f * cos((120.0 / 180.0) * pi);
U_2_shm(8).y = IP_Uc_f * sin((120.0 / 180.0) * pi);


// Нижнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1_A.x = IP_P1_A * 1000000.0;
S1_A.y = IP_Q1_A * 1000000.0;
S1_B.x = IP_P1_B * 1000000.0;
S1_B.y = IP_Q1_B * 1000000.0;
S1_C.x = IP_P1_C * 1000000.0;
S1_C.y = IP_Q1_C * 1000000.0;

S2_A.x = IP_P2_A * 1000000.0;
S2_A.y = IP_Q2_A * 1000000.0;
S2_B.x = IP_P2_B * 1000000.0;
S2_B.y = IP_Q2_B * 1000000.0;
S2_C.x = IP_P2_C * 1000000.0;
S2_C.y = IP_Q2_C * 1000000.0;

S3_A.x = IP_P3_A * 1000000.0;
S3_A.y = IP_Q3_A * 1000000.0;
S3_B.x = IP_P3_B * 1000000.0;
S3_B.y = IP_Q3_B * 1000000.0;
S3_C.x = IP_P3_C * 1000000.0;
S3_C.y = IP_Q3_C * 1000000.0;

I_1_shm(0).x = - (S1_A / U_1_shm(0)).x;
I_1_shm(0).y = (S1_A / U_1_shm(0)).y;
I_1_shm(1).x = - (S1_B / U_1_shm(1)).x;
I_1_shm(1).y = (S1_B / U_1_shm(1)).y;
I_1_shm(2).x = - (S1_C / U_1_shm(2)).x;
I_1_shm(2).y = (S1_C / U_1_shm(2)).y;
I_1_shm(3).x = - (S2_A / U_1_shm(3)).x;
I_1_shm(3).y = (S2_A / U_1_shm(3)).y;
I_1_shm(4).x = - (S2_B / U_1_shm(4)).x;
I_1_shm(4).y = (S2_B / U_1_shm(4)).y;
I_1_shm(5).x = - (S2_C / U_1_shm(5)).x;
I_1_shm(5).y = (S2_C / U_1_shm(5)).y;
I_1_shm(6).x = - (S3_A / U_1_shm(6)).x;
I_1_shm(6).y = (S3_A / U_1_shm(6)).y;
I_1_shm(7).x = - (S3_B / U_1_shm(7)).x;
I_1_shm(7).y = (S3_B / U_1_shm(7)).y;
I_1_shm(8).x = - (S3_C / U_1_shm(8)).x;
I_1_shm(8).y = (S3_C / U_1_shm(8)).y;

for(int j = 0 ; j < 9; j++)
{
temp1_1(j).x = 0;
temp1_1(j).y = 0;
}
for(int j = 0 ; j < 9; j++)
 {
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yfb_1_shm(j,i) * Ubal_1_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 9; i++)
	{
	temp2_1(i) = I_1_shm(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_1(j) += Yf1_1_shm_inv(j,i) * temp2_1(i);
		}
	}
for (int i = 0; i < 9; i++)
{
U_1_shm(i) = temp3_1(i);
}
};

// Верхнее значение
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

S1_A.x = IP_P1_A * 1000000.0;
S1_A.y = IP_Q1_A * 1000000.0;
S1_B.x = IP_P1_B * 1000000.0;
S1_B.y = IP_Q1_B * 1000000.0;
S1_C.x = IP_P1_C * 1000000.0;
S1_C.y = IP_Q1_C * 1000000.0;

S2_A.x = IP_P2_A * 1000000.0;
S2_A.y = IP_Q2_A * 1000000.0;
S2_B.x = IP_P2_B * 1000000.0;
S2_B.y = IP_Q2_B * 1000000.0;
S2_C.x = IP_P2_C * 1000000.0;
S2_C.y = IP_Q2_C * 1000000.0;

S3_A.x = IP_P3_A * 1000000.0;
S3_A.y = IP_Q3_A * 1000000.0;
S3_B.x = IP_P3_B * 1000000.0;
S3_B.y = IP_Q3_B * 1000000.0;
S3_C.x = IP_P3_C * 1000000.0;
S3_C.y = IP_Q3_C * 1000000.0;

I_2_shm(0).x = - (S1_A / U_2_shm(0)).x;
I_2_shm(0).y = (S1_A / U_2_shm(0)).y;
I_2_shm(1).x = - (S1_B / U_2_shm(1)).x;
I_2_shm(1).y = (S1_B / U_2_shm(1)).y;
I_2_shm(2).x = - (S1_C / U_2_shm(2)).x;
I_2_shm(2).y = (S1_C / U_2_shm(2)).y;
I_2_shm(3).x = - (S2_A / U_2_shm(3)).x;
I_2_shm(3).y = (S2_A / U_2_shm(3)).y;
I_2_shm(4).x = - (S2_B / U_2_shm(4)).x;
I_2_shm(4).y = (S2_B / U_2_shm(4)).y;
I_2_shm(5).x = - (S2_C / U_2_shm(5)).x;
I_2_shm(5).y = (S2_C / U_2_shm(5)).y;
I_2_shm(6).x = - (S3_A / U_2_shm(6)).x;
I_2_shm(6).y = (S3_A / U_2_shm(6)).y;
I_2_shm(7).x = - (S3_B / U_2_shm(7)).x;
I_2_shm(7).y = (S3_B / U_2_shm(7)).y;
I_2_shm(8).x = - (S3_C / U_2_shm(8)).x;
I_2_shm(8).y = (S3_C / U_2_shm(8)).y;

for(int j = 0 ; j < 9; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yfb_2_shm(j,i) * Ubal_2_shm(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 9; i++)
	{
	temp2_2(i) = I_2_shm(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 9; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 9; i++)
		{
		temp3_2(j) += Yf1_2_shm_inv(j,i) * temp2_2(i);
		}
	}

for (int i = 0; i < 9; i++)
{
U_2_shm(i) = temp3_2(i);
}

};

// Итоговый отчет
// Нижнее значение
ap::complex U_a1_1 = temp3_1(0) / 1000.0;
ap::complex U_b1_1 = temp3_1(1) / 1000.0;
ap::complex U_c1_1 = temp3_1(2) / 1000.0;
ap::complex U_a2_1 = temp3_1(3) / 1000.0;
ap::complex U_b2_1 = temp3_1(4) / 1000.0;
ap::complex U_c2_1 = temp3_1(5) / 1000.0;
ap::complex U_a3_1 = temp3_1(6) / 1000.0;
ap::complex U_b3_1 = temp3_1(7) / 1000.0;
ap::complex U_c3_1 = temp3_1(8) / 1000.0;

U_a1_1_mod = sqrt(pow((temp3_1(0).x), 2.0)+pow((temp3_1(0).y), 2.0)) / 1000.0;
U_b1_1_mod = sqrt(pow((temp3_1(1).x), 2.0)+pow((temp3_1(1).y), 2.0)) / 1000.0;
U_c1_1_mod = sqrt(pow((temp3_1(2).x), 2.0)+pow((temp3_1(2).y), 2.0)) / 1000.0;
U_a2_1_mod = sqrt(pow((temp3_1(3).x), 2.0)+pow((temp3_1(3).y), 2.0)) / 1000.0;
U_b2_1_mod = sqrt(pow((temp3_1(4).x), 2.0)+pow((temp3_1(4).y), 2.0)) / 1000.0;
U_c2_1_mod = sqrt(pow((temp3_1(5).x), 2.0)+pow((temp3_1(5).y), 2.0)) / 1000.0;
U_a3_1_mod = sqrt(pow((temp3_1(6).x), 2.0)+pow((temp3_1(6).y), 2.0)) / 1000.0;
U_b3_1_mod = sqrt(pow((temp3_1(7).x), 2.0)+pow((temp3_1(7).y), 2.0)) / 1000.0;
U_c3_1_mod = sqrt(pow((temp3_1(8).x), 2.0)+pow((temp3_1(8).y), 2.0)) / 1000.0;

// Выбор узла
if (triangle_ptype == 0)
{
U4f_1_mod = U_a1_1_mod;
U5f_1_mod = U_b1_1_mod;
U6f_1_mod = U_c1_1_mod;
}
if (triangle_ptype == 1)
{
U4f_1_mod = U_a2_1_mod;
U5f_1_mod = U_b2_1_mod;
U6f_1_mod = U_c2_1_mod;
}

// Верхнее значение
ap::complex U_a1_2 = temp3_2(0) / 1000.0;
ap::complex U_b1_2 = temp3_2(1) / 1000.0;
ap::complex U_c1_2 = temp3_2(2) / 1000.0;
ap::complex U_a2_2 = temp3_2(3) / 1000.0;
ap::complex U_b2_2 = temp3_2(4) / 1000.0;
ap::complex U_c2_2 = temp3_2(5) / 1000.0;
ap::complex U_a3_2 = temp3_2(6) / 1000.0;
ap::complex U_b3_2 = temp3_2(7) / 1000.0;
ap::complex U_c3_2 = temp3_2(8) / 1000.0;

U_a1_2_mod = sqrt(pow((temp3_2(0).x), 2.0)+pow((temp3_2(0).y), 2.0)) / 1000.0;
U_b1_2_mod = sqrt(pow((temp3_2(1).x), 2.0)+pow((temp3_2(1).y), 2.0)) / 1000.0;
U_c1_2_mod = sqrt(pow((temp3_2(2).x), 2.0)+pow((temp3_2(2).y), 2.0)) / 1000.0;
U_a2_2_mod = sqrt(pow((temp3_2(3).x), 2.0)+pow((temp3_2(3).y), 2.0)) / 1000.0;
U_b2_2_mod = sqrt(pow((temp3_2(4).x), 2.0)+pow((temp3_2(4).y), 2.0)) / 1000.0;
U_c2_2_mod = sqrt(pow((temp3_2(5).x), 2.0)+pow((temp3_2(5).y), 2.0)) / 1000.0;
U_a3_2_mod = sqrt(pow((temp3_2(6).x), 2.0)+pow((temp3_2(6).y), 2.0)) / 1000.0;
U_b3_2_mod = sqrt(pow((temp3_2(7).x), 2.0)+pow((temp3_2(7).y), 2.0)) / 1000.0;
U_c3_2_mod = sqrt(pow((temp3_2(8).x), 2.0)+pow((temp3_2(8).y), 2.0)) / 1000.0;


cout << "mod U_a1 = [ " << U_a1_1_mod << " , " << U_a1_2_mod << " ] кВ\n"; 
cout << "mod U_b1 = [ " << U_b1_1_mod << " , " << U_b1_2_mod << " ] кВ\n"; 
cout << "mod U_c1 = [ " << U_c1_1_mod << " , " << U_c1_2_mod << " ] кВ\n"; 
cout << "mod U_a2 = [ " << U_a2_1_mod << " , " << U_a2_2_mod << " ] кВ\n"; 
cout << "mod U_b2 = [ " << U_b2_1_mod << " , " << U_b2_2_mod << " ] кВ\n"; 
cout << "mod U_c2 = [ " << U_c2_1_mod << " , " << U_c2_2_mod << " ] кВ\n"; 
cout << "mod U_a3 = [ " << U_a3_1_mod << " , " << U_a3_2_mod << " ] кВ\n"; 
cout << "mod U_b3 = [ " << U_b3_1_mod << " , " << U_b3_2_mod << " ] кВ\n"; 
cout << "mod U_c3 = [ " << U_c3_1_mod << " , " << U_c3_2_mod << " ] кВ\n"; 


}
//-----------------------------------------------------------------------------------------------------
void setintvloss_triangle_fault()
{
	//
cout << "Расчет участка сети. Треугольник.\n";
Yc1_temp_1.setlength(6,6);
Yc1_temp_2.setlength(6,6);
Yc1_temp_3.setlength(6,6);
Yc1_temp_4.setlength(6,6);

Yc2_temp_1.setlength(6,6);
Yc2_temp_2.setlength(6,6);
Yc2_temp_3.setlength(6,6);
Yc2_temp_4.setlength(6,6);

triangle_ptype = 0;

IP_P1_A = 36.0;
IP_Q1_A = 27.0;
IP_P1_B = 36.0;
IP_Q1_B = 27.0;
IP_P1_C = 36.0;
IP_Q1_C = 27.0;

IP_P2_A = 0.0;
IP_Q2_A = 0.0;
IP_P2_B = 0.0;
IP_Q2_B = 0.0;
IP_P2_C = 0.0;
IP_Q2_C = 0.0;

IP_P3_A = 0.0;
IP_Q3_A = 0.0;
IP_P3_B = 0.0;
IP_Q3_B = 0.0;
IP_P3_C = 0.0;
IP_Q3_C = 0.0;


Length = 45.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();
for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_1(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_1(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_1(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_1(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 25.0;
cout << "L = " << Length << " км" << endl << "\n";

getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_2(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_2(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_2(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_2(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 0.0001;
cout << "L = " << Length << " км" << endl << "\n";

phase_fault = 1;
getintvloss();
phase_fault = 0;

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_3(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_3(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_3(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_3(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

Length = 35.0;
cout << "L = " << Length << " км" << endl << "\n";
getintvloss();

for (int i=0; i < 6; i++)
{
	for (int j=0; j < 6; j++)
	{
	Yc1_temp_4(i,j).x = _double(Inf(Re(RP_Yc[i][j])));
	Yc1_temp_4(i,j).y = _double(Inf(Im(RP_Yc[i][j])));
	Yc2_temp_4(i,j).x = _double(Sup(Re(RP_Yc[i][j])));
	Yc2_temp_4(i,j).y = _double(Sup(Im(RP_Yc[i][j])));
	}
}

getintvloss_triangle_4point();

cout << "Process is OK" << "\n";
}

//=======================================================================================================
void get_shortcircuit_2f_sim()
{
// Расчет токов короткого замыкания. Двухфазное КЗ. Метод симметричных составляющих.

SC_UL = 10.5;
SC_L = 5.0;
SC_L_1 = 9.5;
SC_L_2 = 10.5;
SC_X0 = 0.4;
SC_X0_1 = 0.38;
SC_X0_2 = 0.41;

Inf(iSC_UL) = _real(SC_UL);
Sup(iSC_UL) = _real(SC_UL);
Inf(iSC_L) = _real(SC_L_1);
Sup(iSC_L) = _real(SC_L_2);
Inf(iSC_X0) = _real(SC_X0_1);
Sup(iSC_X0) = _real(SC_X0_2);
iSC_UF = iSC_UL / sqrt(3.0);
iSC_X1 = iSC_X0 * iSC_L;
iSC_X0 = 3.5 * iSC_X0;

Re(iSC_temp_01) = iSC_UF;
Im(iSC_temp_01) = 0.0;
Re(iSC_temp_02) = 0.0;
Im(iSC_temp_02) = 2.0 * iSC_X1;
iSC_IA1 = iSC_temp_01 / iSC_temp_02;

Re(iSC_temp_01) = 0.0;
Im(iSC_temp_01) = (2.0 / 3.0) * Pi();
iSC_a = exp(iSC_temp_01);

iSC_IA2 = - iSC_IA1;

Re(iSC_temp_01) = 0.0;
Im(iSC_temp_01) = iSC_X1;
iSC_UA1 = iSC_temp_01 * iSC_IA1;

iSC_UA2 = iSC_UA1;

iSC_IA = iSC_IA1 + iSC_IA2;
iSC_IB = power(iSC_a, 2.0) * iSC_IA1 + iSC_a * iSC_IA2;
iSC_IC = iSC_a * iSC_IA1 + power(iSC_a, 2.0) * iSC_IA2;
iSC_IS = iSC_IA + iSC_IB + iSC_IC;

iSC_UA = iSC_UA1 + iSC_UA2;
iSC_UB = power(iSC_a, 2.0) * iSC_UA1 + iSC_a * iSC_UA2;
iSC_UC = iSC_a * iSC_UA1 + power(iSC_a, 2.0) * iSC_UA2;

cout << "IA = " << iSC_IA << "\n";
cout << "IB = " << iSC_IB << "\n";
cout << "IC = " << iSC_IC << "\n";


}
//===============================================================================================================
//===============================================================================================================
void get_unsimmetric_fcoor()
{
complex cSC_Unom, cSC_U2, cSC_a, cSC_aa, cSC_2_3_of_Pi, cSC_30_180_of_Pi, cSC_fi2, cSC_B11, cSC_B12;
Re(cSC_Unom) = 6.1;
Im(cSC_Unom) = 0.0;
cSC_U2 = 0.06 * cSC_Unom;
Re(cSC_2_3_of_Pi) = 0.0;
Im(cSC_2_3_of_Pi) = double((2.0 / 3.0) * pi);
Re(cSC_30_180_of_Pi) = 0.0;
Im(cSC_30_180_of_Pi) = double((30.0 / 180.0) * pi);
Re(cSC_fi2) = 0.0;
Im(cSC_fi2) = double(0.0 * pi);
cSC_a = exp(cSC_2_3_of_Pi);
cSC_aa = pow(cSC_a, 2.0);
cSC_B11 = cSC_Unom * exp(cSC_30_180_of_Pi);
cSC_B12 = cSC_U2 * exp(cSC_fi2);

ap::complex_2d_array ISC_S, ISC_B1;
ISC_S.setlength(3,3);
ISC_B1.setlength(3,1);
ISC_U_unsim.setlength(3,1);
ISC_S(0,0).x = 1.0;
ISC_S(0,0).y = 0.0;
ISC_S(0,1).x = 1.0;
ISC_S(0,1).y = 0.0;
ISC_S(0,2).x = 1.0;
ISC_S(0,2).y = 0.0;
ISC_S(1,0).x = _double(Re(cSC_aa));
ISC_S(1,0).y = _double(Im(cSC_aa));
ISC_S(1,1).x = _double(Re(cSC_a));
ISC_S(1,1).y = _double(Im(cSC_a));
ISC_S(1,2).x = 1.0;
ISC_S(1,2).y = 0.0;
ISC_S(2,0).x = _double(Re(cSC_a));
ISC_S(2,0).y = _double(Im(cSC_a));
ISC_S(2,1).x = _double(Re(cSC_aa));
ISC_S(2,1).y = _double(Im(cSC_aa));
ISC_S(2,2).x = 1.0;
ISC_S(2,2).y = 0.0;

ISC_B1(0,0).x = _double(Re(cSC_B11));
ISC_B1(0,0).y = _double(Im(cSC_B11));
ISC_B1(1,0).x = _double(Re(cSC_B12));
ISC_B1(1,0).y = _double(Im(cSC_B12));
ISC_B1(2,0).x = 0.0;
ISC_B1(2,0).y = 0.0;

for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 ISC_U_unsim(i,j) += ISC_S(i,k) * ISC_B1(k,j);
	  }
   }
}

cinterval iISC_U_insim[3][1];
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_U_insim[i][j])) = real(ISC_U_unsim(i,j).x);
	Inf(Im(iISC_U_insim[i][j])) = real(ISC_U_unsim(i,j).y);
	Sup(Re(iISC_U_insim[i][j])) = real(ISC_U_unsim(i,j).x);
	Sup(Im(iISC_U_insim[i][j])) = real(ISC_U_unsim(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "U [" << i << "] = " << iISC_U_insim[i][j] << "\n";
	}
}

}
//===============================================================================================================
//===============================================================================================================
void get_shortcircuit_2f_fcoor()
{
// Расчет токов короткого замыкания. Двухфазное КЗ. Расчет в фазных координатах.
SC_UL = 10.5;
SC_UL_1 = 10.5;
SC_UL_2 = 10.5;
SC_L = 5.0;
SC_L_1 = 9.0;
SC_L_2 = 11.0;
SC_X0 = 0.4;
SC_X0_1 = 0.32;
SC_X0_2 = 0.48;
interval iSC_XL, iSC_XM, iSC_35;
cinterval iSC_ZL, iSC_ZM, iSC_C1;
cinterval iSC_Z1[3][3], iSC_Z2[3][3], iSC_Y1[3][3], iSC_Y2[3][3];
cinterval iSC_YC1[6][6], iSC_YC2[6][6], iSC_O[6][6], iSC_YB[12][12], iSC_YS1[9][9], iSC_M[9][12], iSC_M_trans[12][9];
interval iSC_P[3][6], iSC_O1[3][6];

ap::complex_2d_array apSC_Z1_inf;
ap::complex_2d_array apSC_Z1_sup;
apSC_Z1_inf.setlength(3,3);
apSC_Z1_sup.setlength(3,3);

ap::complex_2d_array apSC_Z2_inf;
ap::complex_2d_array apSC_Z2_sup;
apSC_Z2_inf.setlength(3,3);
apSC_Z2_sup.setlength(3,3);
//===============================================================================================================
Inf(iSC_UL) = _real(SC_UL_1);
Sup(iSC_UL) = _real(SC_UL_2);
Inf(iSC_L) = _real(SC_L_1);
Sup(iSC_L) = _real(SC_L_2);
cout << "L = " << iSC_L << "км \n";
Inf(iSC_X0) = _real(SC_X0_1);
Sup(iSC_X0) = _real(SC_X0_2);
Inf(iSC_35) = _real(3.5);
Sup(iSC_35) = _real(3.5);
iSC_UF = iSC_UL / sqrt(3.0);
iSC_X1 = iSC_X0 * iSC_L;
iSC_X0 = iSC_35 * iSC_X1;
iSC_XL = (iSC_X1 / 1.5) + (iSC_X0 / 3.0);
iSC_XM = iSC_XL - iSC_X1;
Re(iSC_ZL) = 0.0;
Im(iSC_ZL) = iSC_XL;
Re(iSC_ZM) = 0.0;
Im(iSC_ZM) = iSC_XM;
// Проверка
iSC_C1 = iSC_XL - iSC_XM;
//===============================================================================================================
iSC_Z1[0][0] = iSC_ZL;
iSC_Z1[0][1] = iSC_ZM;
iSC_Z1[0][2] = iSC_ZM;
iSC_Z1[1][0] = iSC_ZM;
iSC_Z1[1][1] = iSC_ZL;
iSC_Z1[1][2] = iSC_ZM;
iSC_Z1[2][0] = iSC_ZM;
iSC_Z1[2][1] = iSC_ZM;
iSC_Z1[2][2] = iSC_ZL;

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z1_inf(i,j).x = _double(Inf(Re(iSC_Z1[i][j])));
	apSC_Z1_inf(i,j).y = _double(Inf(Im(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).x = _double(Sup(Re(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).y = _double(Sup(Im(iSC_Z1[i][j])));
	}
}

cmatrixinverse(apSC_Z1_inf,3);
cmatrixinverse(apSC_Z1_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).x;
	Inf(Im(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).y;
	Sup(Re(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).x;
	Sup(Im(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).y;
	}
}
//===============================================================================================================
/*// Двухфазное КЗ (между фазами)
Re(iSC_Z2[0][0]) = 0.0;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 10000.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 10000.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 10000.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 0.0;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.001;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 10000.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.001;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 0.0;
Im(iSC_Z2[2][2]) = 0.0;
*/

// Однофазное КЗ на землю
Re(iSC_Z2[0][0]) = 100000.0;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 0.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 0.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 0.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 0.0;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.0;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 0.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.0;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 0.0;
Im(iSC_Z2[2][2]) = 0.0;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z2_inf(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	apSC_Z2_inf(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
	}
}

cmatrixinverse(apSC_Z2_inf,3);
cmatrixinverse(apSC_Z2_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).x;
	Inf(Im(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).y;
	Sup(Re(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).x;
	Sup(Im(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).y;
	}
}

//===============================================================================================================
// Матрица YC1
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i][j]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i][j]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i][j]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i][j-3]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i][j-3]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i][j-3]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i-3][j]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i-3][j]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i-3][j]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i-3][j-3]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i-3][j-3]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i-3][j-3]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i-3][j-3]));
	}
}
//===============================================================================================================
// Матрица YC2
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i][j]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i][j]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i][j]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i][j-3]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i][j-3]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i][j-3]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i-3][j]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i-3][j]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i-3][j]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i-3][j-3]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i-3][j-3]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i-3][j-3]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i-3][j-3]));
	}
}

//===============================================================================================================
// Матрица P
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_P[i][j] = 0.0;
	}
}
iSC_P[0][0] = 1.0;
iSC_P[1][1] = 1.0;
iSC_P[2][2] = 1.0;
//===============================================================================================================
// Матрица O1
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O1[i][j] = 0.0;
	}
}
//===============================================================================================================
// Матрица M
for (int i=0;i<9;i++)
{
	for (int j=0;j<12;j++)
	{
	Re(iSC_M[i][j]) = 0.0;
	Im(iSC_M[i][j]) = 0.0;
	}
}
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i][j];
	}
}
for (int i=0;i<3;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i][j-6];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-3][j];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i-3][j-6];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i-6][j];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-6][j-6];
	}
}

//===============================================================================================================
//Матрица O

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O[i][j] = 0.0;
	}
}
//===============================================================================================================
//Матрица YB
for (int i=0;i<12;i++)
{
	for (int j=0;j<12;j++)
	{
	iSC_YB[i][j] = 0.0;
	}
}
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_YC1[i][j];
	}
}
for (int i=0;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_O[i][j-6];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_O[i-6][j];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_YC2[i-6][j-6];
	}
}
//===============================================================================================================
// M * Yb = (9x12) * (12x12)
cinterval iSC_Ytemp[9][12];
ap::complex_2d_array Y_1_shm, Y_2_shm, Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm.setlength(9,9);
Y_2_shm.setlength(9,9);
Y_1_shm_temp.setlength(9,12);
Y_2_shm_temp.setlength(9,12);
Y_1_shm_temp_reint.setlength(9,12);
Y_2_shm_temp_reint.setlength(9,12);
ap::complex_2d_array Yb_1_shm, Yb_2_shm;
Yb_1_shm.setlength(12,12);
Yb_2_shm.setlength(12,12);
ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
ap::complex_2d_array M_shm;
M_shm_transposed.setlength(12,9);
M_shm.setlength(9,12);

for (int i = 0; i < 9; i++)     
{
   for (int j = 0; j < 12; j++) 
   {
	M_shm(i,j).x = _double(Inf(Re(iSC_M[i][j])));
	M_shm(i,j).y = _double(0.0);
   }
}
for (int i = 0; i < 12; i++)      
{
   for (int j = 0; j < 12; j++)  
   {
	Yb_1_shm(i,j).x = _double(Inf(Re(iSC_YB[i][j])));
	Yb_1_shm(i,j).y = _double(Inf(Im(iSC_YB[i][j])));
	Yb_2_shm(i,j).x = _double(Sup(Re(iSC_YB[i][j])));
	Yb_2_shm(i,j).y = _double(Sup(Im(iSC_YB[i][j])));
   }
}

for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {
		Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}


for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}


// Transpose M
for (int i=0; i < 9; i++)
	{
	for (int j=0; j < 12 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}
	}



// Y_temp * M^transpose = (9x12) * (12x9)
for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 9; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {

		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}
ap::complex_2d_array SC_Y1_1, SC_Y1_2, SC_Y2_1, SC_Y2_2, SC_Y12_1, SC_Y12_2, SC_Y21_1, SC_Y21_2;
SC_Y1_1.setlength(6,6);
SC_Y1_2.setlength(6,6);
SC_Y2_1.setlength(3,3);
SC_Y2_2.setlength(3,3);
SC_Y12_1.setlength(6,3);
SC_Y12_2.setlength(6,3);
SC_Y21_1.setlength(3,6);
SC_Y21_2.setlength(3,6);
for (int i = 0; i < 6; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y1_1(i,j) = Y_1_shm(i,j);
    SC_Y1_2(i,j) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y2_1(i-6,j-6) = Y_1_shm(i,j);
    SC_Y2_2(i-6,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 0; i < 6; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y12_1(i,j-6) = Y_1_shm(i,j);
    SC_Y12_2(i,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y21_1(i-6,j) = Y_1_shm(i,j);
    SC_Y21_2(i-6,j) = Y_2_shm(i,j);
    }
}
ap::complex_2d_array SC_e1, SC_e1_T, SC_Y12_1_T, SC_Y12_2_T, SC_C1_1, SC_C1_2, SC_C2_1, SC_C2_2, SC_C3_1, SC_C3_2;
ap::complex_2d_array SC_S1_1, SC_S1_2, SC_S2_1, SC_S2_2, SC_YE_1, SC_YE_2;
SC_S1_1.setlength(7,1);
SC_S1_2.setlength(7,1);
SC_S2_1.setlength(7,6);
SC_S2_2.setlength(7,6);
SC_YE_1.setlength(7,7);
SC_YE_2.setlength(7,7);
SC_e1.setlength(3,1);
SC_e1_T.setlength(1,3);
SC_C1_1.setlength(6,1);
SC_C1_2.setlength(6,1);
SC_C2_1.setlength(1,6);
SC_C2_2.setlength(1,6);
SC_C3_1.setlength(1,1);
SC_C3_2.setlength(1,1);
SC_Y12_1_T.setlength(3,6);
SC_Y12_2_T.setlength(3,6);
for (int i = 0; i < 3; i++)
{
	for (int j = 0; j < 1; j++)
	{
	SC_e1(i,j).x = 1.0;
    SC_e1(i,j).y = 0.0;
    }
}

// C1 = Y12 * e1
for (int i = 0; i < 6; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C1_1(i,j) += SC_Y12_1(i,k) * SC_e1(k,j);
		 SC_C1_2(i,j) += SC_Y12_2(i,k) * SC_e1(k,j);
	  }
   }
}

// Transpose e1
for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1 ; j++)
		{
		M_element = SC_e1(i,j);
		SC_e1_T(j,i) = M_element;
		}
	}

// Transpose Y12
ap::complex ELEM_1, ELEM_2;
for (int i=0; i < 6; i++)
	{
	for (int j=0; j < 3 ; j++)
		{
		ELEM_1 = SC_Y12_1(i,j);
		ELEM_2 = SC_Y12_2(i,j);
		SC_Y12_1_T(j,i) = ELEM_1;
		SC_Y12_2_T(j,i) = ELEM_2;
		}
	}

// C2 = e1^T * Y12^T
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 6; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C2_1(i,j) += SC_e1_T(i,k) * SC_Y12_1_T(k,j);
		 SC_C2_2(i,j) += SC_e1_T(i,k) * SC_Y12_2_T(k,j);
	  }
   }
}
// C3 = e1^T * Y2 * e1
ap::complex_2d_array SC_TEMP_1, SC_TEMP_2;
SC_TEMP_1.setlength(1,3);
SC_TEMP_2.setlength(1,3);
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 3; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_TEMP_1(i,j) += SC_e1_T(i,k) * SC_Y2_1(k,j);
		 SC_TEMP_2(i,j) += SC_e1_T(i,k) * SC_Y2_2(k,j);
	  }
   }
}
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C3_1(i,j) += SC_TEMP_1(i,k) * SC_e1(k,j);
		 SC_C3_2(i,j) += SC_TEMP_2(i,k) * SC_e1(k,j);
	  }
   }
}
//===============================================================================================================
// S1
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C1_1(i,j);
	SC_S1_2(i,j) = SC_C1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C3_1(i-6,j);
	SC_S1_2(i,j) = SC_C3_2(i-6,j);
   }
}
//===============================================================================================================
// S2 = STACK(Y1, C2) <=> 6x6 & 1x6
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_Y1_1(i,j);
	SC_S2_2(i,j) = SC_Y1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_C2_1(i-6,j);
	SC_S2_2(i,j) = SC_C2_2(i-6,j);
   }
}
//===============================================================================================================
// YE = augment(S2, S1) = 7x6 & 7x1
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_YE_1(i,j) = SC_S2_1(i,j);
	SC_YE_2(i,j) = SC_S2_2(i,j);
   }
}
for (int i = 0; i < 7; i++)    
{
   for (int j = 6; j < 7; j++)  
   {
	SC_YE_1(i,j) = SC_S1_1(i,j-6);
	SC_YE_2(i,j) = SC_S1_2(i,j-6);
   }
}
/*
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 7; j++)  
   {
	cout << "YE.x ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).x << "\n";
	cout << "YE.y ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).y << "\n";
   }
}
*/
//===============================================================================================================
// UBE
cinterval iSC_UBE[3][1];
ap::complex_2d_array SC_UBE_1, SC_UBE_2;
SC_UBE_1.setlength(3,1);
SC_UBE_2.setlength(3,1);
cinterval SC_2_3_of_Pi;
Re(SC_2_3_of_Pi) = 0.0;
Im(SC_2_3_of_Pi) = double((2.0 / 3.0) * pi);

if (param_unsim == 0)
{
// then 0
iSC_UBE[0][0] = 1.0;
iSC_UBE[1][0] = exp(- SC_2_3_of_Pi);
iSC_UBE[2][0] = exp(SC_2_3_of_Pi);
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	   iSC_UBE[i][j] = iSC_UBE[i][j] * iSC_UF;
   }
}
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_UBE_1(i,j).x = _double(Inf(Re(iSC_UBE[i][j])));
	SC_UBE_1(i,j).y = _double(Inf(Im(iSC_UBE[i][j])));
	SC_UBE_2(i,j).x = _double(Sup(Re(iSC_UBE[i][j])));
	SC_UBE_2(i,j).y = _double(Sup(Im(iSC_UBE[i][j])));
   }
}
//end 0
}
else if (param_unsim == 1)
{
// then 1
get_unsimmetric_fcoor();
	for (int i = 0; i < 3; i++)    
	{
		for (int j = 0; j < 1; j++)  
		{
		SC_UBE_1(i,j) = ISC_U_unsim(i,j);
		SC_UBE_2(i,j) = ISC_U_unsim(i,j);
		}
}

}

//cout << "y = " << iSC_UBE[1][0] << "\n";
//===============================================================================================================
// Y11
ap::complex_2d_array SC_Y11_1, SC_Y11_2, SC_Y11_1_T, SC_Y11_2_T, SC_Y22_1, SC_Y22_2, SC_V_1, SC_V_2; 
ap::complex SC_DYE_1, SC_DYE_2;
SC_Y11_1.setlength(4,4);
SC_Y11_2.setlength(4,4);
SC_Y11_1_T.setlength(4,4);
SC_Y11_2_T.setlength(4,4);
SC_Y22_1.setlength(4,3);
SC_Y22_2.setlength(4,3);
SC_V_1.setlength(4,1);
SC_V_2.setlength(4,1);
for (int i = 3; i < 7; i++)    
{
   for (int j = 3; j < 7; j++)  
   {
	   SC_Y11_1(i-3,j-3) = SC_YE_1(i,j);
	   SC_Y11_2(i-3,j-3) = SC_YE_2(i,j);
   }
}
for (int i = 3; i < 7; i++)    
{
   for (int j = 0; j < 3; j++)  
   {
	   SC_Y22_1(i-3,j) = SC_YE_1(i,j);
	   SC_Y22_2(i-3,j) = SC_YE_2(i,j);
   }
}

SC_DYE_1 = cmatrixdet(SC_Y11_1, 4);
SC_DYE_2 = cmatrixdet(SC_Y11_2, 4);
cout << "DYE_1.x = " << SC_DYE_1.x << "\n";
cout << "DYE_1.y = " << SC_DYE_1.y << "\n";
cout << "DYE_2.x = " << SC_DYE_2.x << "\n";
cout << "DYE_2.y = " << SC_DYE_2.y << "\n";
/*
ap::complex a;
ap::complex_2d_array b = ;
a = cmatrixdet(b, 2);
*/

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_V_1(i,j) += SC_Y22_1(i,k) * SC_UBE_1(k,j);
		 SC_V_2(i,j) += SC_Y22_2(i,k) * SC_UBE_2(k,j);
	  }
   }
}

for (int i = 0; i < 4; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "V.x [" << i << "] = " << SC_V_1(i,j).x << "\n";
	cout << "V.y [" << i << "] = " << SC_V_1(i,j).y << "\n";
	}
}

//===============================================================================================================
// US
ap::complex_2d_array SC_US_1, SC_US_2, SC_US1_1, SC_US1_2, SC_IK_1, SC_IK_2, SC_Z2_1, SC_Z2_2;
SC_US_1.setlength(4,1);
SC_US_2.setlength(4,1);
SC_US1_1.setlength(4,1);
SC_US1_2.setlength(4,1);
SC_IK_1.setlength(3,1);
SC_IK_2.setlength(3,1);
SC_Z2_1.setlength(3,3);
SC_Z2_2.setlength(3,3);

// Inverse Z2
cmatrixinverse(SC_Y11_1,4);
cmatrixinverse(SC_Y11_2,4);

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 4; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_US_1(i,j) += SC_Y11_1(i,k) * - SC_V_1(k,j);
		 SC_US_2(i,j) += SC_Y11_2(i,k) * - SC_V_2(k,j);
	  }
   }
}

for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1; j++)
		{
		SC_US1_1(i,j) = SC_US_1(i,j);
		SC_US1_2(i,j) = SC_US_2(i,j);
		}
	}

for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_Z2_1(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	SC_Z2_1(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	SC_Z2_2(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	SC_Z2_2(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
   }
}

// Inverse Z2
//cmatrixinverse(SC_Z2_1,3);
//cmatrixinverse(SC_Z2_2,3);

for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_IK_1(i,j) += apSC_Z2_inf(i,k) * SC_US1_1(k,j);
		 SC_IK_2(i,j) += apSC_Z2_sup(i,k) * SC_US1_2(k,j);
	  }
   }
}

//cout << "US1.x [0,0] = " << SC_US1_1(0,0).x << "\n";
//cout << "US1.y [0,0] = " << SC_US1_1(0,0).y << "\n";
//cout << "IK.x [0,0] = " << SC_IK_1(0,0).x << "\n";
//cout << "IK.y [0,0] = " << SC_IK_1(0,0).y << "\n";
//===============================================================================================================
// End
cinterval iSC_IK[3][1];

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iSC_IK[i][j])) = real(SC_IK_1(i,j).x);
	Inf(Im(iSC_IK[i][j])) = real(SC_IK_1(i,j).y);
	Sup(Re(iSC_IK[i][j])) = real(SC_IK_2(i,j).x);
	Sup(Im(iSC_IK[i][j])) = real(SC_IK_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "IK [" << i << "] = " << iSC_IK[i][j] << "\n";
	}
}

/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2 [" << i << "," << j << "] = " << iSC_Z2[i][j] << "\n";
	}
}
*/
/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2.x [" << i << "," << j << "] = " << SC_Z2_1(i,j).x << "\n";
	cout << "Z2.y [" << i << "," << j << "] = " << SC_Z2_1(i,j).y << "\n";
	}
}
*/
complex cSC_a, cSC_aa, cSC_2_3_of_Pi;
Re(cSC_2_3_of_Pi) = 0.0;
Im(cSC_2_3_of_Pi) = double((2.0 / 3.0) * pi);
cSC_a = exp(cSC_2_3_of_Pi);
cSC_aa = pow(cSC_a, 2.0);
ap::complex_2d_array ISC_S, ISC_I120_1, ISC_I120_2;
ISC_S.setlength(3,3);
ISC_I120_1.setlength(3,1);
ISC_I120_2.setlength(3,1);
ISC_S(0,0).x = 1.0;
ISC_S(0,0).y = 0.0;
ISC_S(0,1).x = 1.0;
ISC_S(0,1).y = 0.0;
ISC_S(0,2).x = 1.0;
ISC_S(0,2).y = 0.0;
ISC_S(1,0).x = _double(Re(cSC_aa));
ISC_S(1,0).y = _double(Im(cSC_aa));
ISC_S(1,1).x = _double(Re(cSC_a));
ISC_S(1,1).y = _double(Im(cSC_a));
ISC_S(1,2).x = 1.0;
ISC_S(1,2).y = 0.0;
ISC_S(2,0).x = _double(Re(cSC_a));
ISC_S(2,0).y = _double(Im(cSC_a));
ISC_S(2,1).x = _double(Re(cSC_aa));
ISC_S(2,1).y = _double(Im(cSC_aa));
ISC_S(2,2).x = 1.0;
ISC_S(2,2).y = 0.0;
cmatrixinverse(ISC_S, 3);
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 ISC_I120_1(i,j) += ISC_S(i,k) * SC_IK_1(k,j);
		 ISC_I120_2(i,j) += ISC_S(i,k) * SC_IK_2(k,j);
	  }
   }
}
cinterval iISC_I120[3][1];
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_I120[i][j])) = real(ISC_I120_1(i,j).x);
	Inf(Im(iISC_I120[i][j])) = real(ISC_I120_1(i,j).y);
	Sup(Re(iISC_I120[i][j])) = real(ISC_I120_2(i,j).x);
	Sup(Im(iISC_I120[i][j])) = real(ISC_I120_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "I120 [" << i << "] = " << iISC_I120[i][j] << "\n";
	}
}
cout << "End of code" << "\n";
}

/*
DWORD getFreqProcessor()

    SYSTEM_INFO systemInfo;
    HKEY hKey;
    char Buffer[_MAX_PATH];
    long lError;
    DWORD BufSize = _MAX_PATH;
    DWORD dwMHz = _MAX_PATH;

    GetSystemInfo(&systemInfo);

    lError = RegOpenKeyEx(HKEY_LOCAL_MACHINE,
                        "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
                        0,
                        KEY_READ,
                        &hKey);


    if(lError != ERROR_SUCCESS)
    {// if the key is not found, tell the user why:
            FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM,
                NULL,
                lError,
                0,
                Buffer,
                _MAX_PATH,
                0);
            return (LRESULT)"N/A";
    }

        // query the key:
        RegQueryValueEx(hKey, "~MHz", NULL, NULL, (LPBYTE) &dwMHz, &BufSize);
        return dwMHz;
}
*/

//-----------------------------------------------------------------------------------------------------
void get_magnetic_field()
{
// Интервальный анализ электромагнитных полей

cout << "Расчет электромагнитного поля\n";

double MF_X0, MF_X1, MF_X2, MF_X3, MF_Y0, MF_Y1, MF_Y2, MF_Y3, MF_R, MF_MIS;

ap::complex_2d_array	MF_I_1, MF_dU_1, MF_I_2, 
						MF_dU_2, MF_T_low, MF_T_high,
						MF_VU_1, MF_VU_2;

cinterval iMF_R;
MF_I_1.setlength(3,1);
MF_I_2.setlength(3,1);
MF_dU_1.setlength(3,1);
MF_dU_2.setlength(3,1);

MF_R = 0.012;
Re(iMF_R) = MF_R;
Im(iMF_R) = NULL;

MF_MIS = 5.0;			// Ошибка %

MF_X0 = 2.0;
MF_Y0 = 1.8;

MF_X1 = -5.0;
MF_X2 = 0.0;
MF_X3 = 5.0;

MF_Y1 = 19.0;
MF_Y2 = 19.0;
MF_Y3 = 19.0;

cinterval iMF_X[4], iMF_Y[4];
Inf(Re(iMF_X[0])) = _real(min((MF_X0 * (1.0 - MF_MIS / 100.0) ), MF_X0 * (1.0 + MF_MIS / 100.0)) );
Inf(Re(iMF_X[1])) = _real(MF_X1);
Inf(Re(iMF_X[2])) = _real(MF_X2);
Inf(Re(iMF_X[3])) = _real(MF_X3);
Sup(Re(iMF_X[0])) = _real(max((MF_X0 * (1.0 - MF_MIS / 100.0) ), MF_X0 * (1.0 + MF_MIS / 100.0)) );
Sup(Re(iMF_X[1])) = _real(MF_X1);
Sup(Re(iMF_X[2])) = _real(MF_X2);
Sup(Re(iMF_X[3])) = _real(MF_X3);

//cout << "X1 = " << iMF_X[1] << "\n";

Im(iMF_X[0]) = _real(0.0);
Im(iMF_X[1]) = _real(0.0);
Im(iMF_X[2]) = _real(0.0);
Im(iMF_X[3]) = _real(0.0);

Inf(Re(iMF_Y[0])) = _real(min((MF_Y0 * (1.0 - MF_MIS / 100.0) ), MF_Y0 * (1.0 + MF_MIS / 100.0)) );
Inf(Re(iMF_Y[1])) = _real(MF_Y1);
Inf(Re(iMF_Y[2])) = _real(MF_Y2);
Inf(Re(iMF_Y[3])) = _real(MF_Y3);
Sup(Re(iMF_Y[0])) = _real(max((MF_Y0 * (1.0 - MF_MIS / 100.0) ), MF_Y0 * (1.0 + MF_MIS / 100.0)) );
Sup(Re(iMF_Y[1])) = _real(MF_Y1);
Sup(Re(iMF_Y[2])) = _real(MF_Y2);
Sup(Re(iMF_Y[3])) = _real(MF_Y3);

Im(iMF_Y[0]) = _real(0.0);
Im(iMF_Y[1]) = _real(0.0);
Im(iMF_Y[2]) = _real(0.0);
Im(iMF_Y[3]) = _real(0.0);

setintvloss_star();
//setintvloss_triangle();

double MF_E0;
MF_E0 = 8.85418782*pow(10.0,-12.0);

cinterval iMF_E0;
Re(iMF_E0) = _real(MF_E0);
Im(iMF_E0) = 0.0;


//cout << "E0 = " << MF_E0 << "\n";

cinterval iMF_X0, iMF_Y0;
Inf(Re(iMF_X0)) = _real(MF_X0);
Inf(Im(iMF_X0)) = _real(0.0);
Sup(Re(iMF_X0)) = _real(MF_X0);
Sup(Im(iMF_X0)) = _real(0.0);
Inf(Re(iMF_Y0)) = _real(MF_Y0);
Inf(Im(iMF_Y0)) = _real(0.0);
Sup(Re(iMF_Y0)) = _real(MF_Y0);
Sup(Im(iMF_Y0)) = _real(0.0);

cinterval iMF_U_lev[3], iMF_U_prav[3], iMF_dU[3], iMF_ZM[3][3];

/* FOR TRIANGLE
for (int i = 0; i < 3; i++)    
{
Inf(Re(iMF_U_lev[i])) = _real(min(MF_U_low(i+3).x, MF_U_high(i+3).x));
Inf(Im(iMF_U_lev[i])) = _real(min(MF_U_low(i+3).y, MF_U_high(i+3).y));
Sup(Re(iMF_U_lev[i])) = _real(max(MF_U_low(i+3).x, MF_U_high(i+3).x));
Sup(Im(iMF_U_lev[i])) = _real(max(MF_U_low(i+3).y, MF_U_high(i+3).y));

Inf(Re(iMF_U_prav[i])) = _real(min(MF_U_low(i).x, MF_U_high(i).x));
Inf(Im(iMF_U_prav[i])) = _real(min(MF_U_low(i).y, MF_U_high(i).y));
Sup(Re(iMF_U_prav[i])) = _real(max(MF_U_low(i).x, MF_U_high(i).x));
Sup(Im(iMF_U_prav[i])) = _real(max(MF_U_low(i).y, MF_U_high(i).y));
}
*/
// FOR STAR
for (int i = 0; i < 3; i++)    
{
Inf(Re(iMF_U_lev[i])) = _real(min(MF_U_low(i+6).x, MF_U_high(i+6).x));
Inf(Im(iMF_U_lev[i])) = _real(min(MF_U_low(i+6).y, MF_U_high(i+6).y));
Sup(Re(iMF_U_lev[i])) = _real(max(MF_U_low(i+6).x, MF_U_high(i+6).x));
Sup(Im(iMF_U_lev[i])) = _real(max(MF_U_low(i+6).y, MF_U_high(i+6).y));

Inf(Re(iMF_U_prav[i])) = _real(min(MF_U_low(i).x, MF_U_high(i).x));
Inf(Im(iMF_U_prav[i])) = _real(min(MF_U_low(i).y, MF_U_high(i).y));
Sup(Re(iMF_U_prav[i])) = _real(max(MF_U_low(i).x, MF_U_high(i).x));
Sup(Im(iMF_U_prav[i])) = _real(max(MF_U_low(i).y, MF_U_high(i).y));
}


for (int i = 0; i < 3; i++)    
{
iMF_dU[i] = iMF_U_lev[i] - iMF_U_prav[i];
}

for (int i = 0; i < 3; i++)    
{
MF_dU_1(i,0).x = min(_double(Inf(Re(iMF_dU[i]))),_double(Sup(Re(iMF_dU[i]))));
MF_dU_1(i,0).y = min(_double(Inf(Im(iMF_dU[i]))),_double(Sup(Im(iMF_dU[i]))));
MF_dU_2(i,0).x = max(_double(Inf(Re(iMF_dU[i]))),_double(Sup(Re(iMF_dU[i]))));
MF_dU_2(i,0).y = max(_double(Inf(Im(iMF_dU[i]))),_double(Sup(Im(iMF_dU[i]))));
}

//cout << "dU1 = " << iMF_dU[0] << "\n";
//cout << "dU2 = " << iMF_dU[1] << "\n";
//cout << "dU2 = " << iMF_dU[2] << "\n";




cmatrixinverse(MF_ZM_low, 3);
cmatrixinverse(MF_ZM_high, 3);

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++) 
	{
	MF_ZM_low(i,j).x = min(MF_ZM_low(i,j).x, MF_ZM_high(i,j).x);
	MF_ZM_low(i,j).y = min(MF_ZM_low(i,j).y, MF_ZM_high(i,j).y);
	MF_ZM_high(i,j).x = max(MF_ZM_low(i,j).x, MF_ZM_high(i,j).x);
	MF_ZM_high(i,j).y = max(MF_ZM_low(i,j).y, MF_ZM_high(i,j).y);
	}
}

//cout << "Z(0,0) MIN = [" << MF_ZM_low(0,0).x << "]+j*[" << MF_ZM_low(0,0).y << "]\n";
//cout << "Z(0,0) MAX = [" << MF_ZM_high(0,0).x << "]+j*[" << MF_ZM_high(0,0).y << "]\n";

//MF_I[3][1] = MF_ZM[3][3] * MD_dU[3][1]
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 MF_I_1(i,j) += MF_ZM_low(i,k) * MF_dU_1(k,j);
		 MF_I_2(i,j) += MF_ZM_high(i,k) * MF_dU_2(k,j);
	  }
   }
}

/*
for (int i = 0; i < 3; i++)     
{
	cout << "MF_I(" << i+1 << ").X = " << MF_I_1(i,0).x << "\n";
	cout << "MF_I(" << i+1 << ").Y = " << MF_I_1(i,0).y << "\n";
}
*/

MF_T_low.setlength(3,1);
MF_T_high.setlength(3,1);
MF_VU_1.setlength(3,1);
MF_VU_2.setlength(3,1);

/*// FOR TRIANGLE
for (int i = 0; i < 3; i++)    
{
MF_VU_1(i,0).x = min(MF_U_low(i+3).x, MF_U_high(i+3).x);
MF_VU_1(i,0).y = min(MF_U_low(i+3).y, MF_U_high(i+3).y);
MF_VU_2(i,0).x = max(MF_U_low(i+3).x, MF_U_high(i+3).x);
MF_VU_2(i,0).y = max(MF_U_low(i+3).y, MF_U_high(i+3).y);
}
*/
// FOR STAR
for (int i = 0; i < 3; i++)    
{
MF_VU_1(i,0).x = min(MF_U_low(i+6).x, MF_U_high(i+6).x);
MF_VU_1(i,0).y = min(MF_U_low(i+6).y, MF_U_high(i+6).y);
MF_VU_2(i,0).x = max(MF_U_low(i+6).x, MF_U_high(i+6).x);
MF_VU_2(i,0).y = max(MF_U_low(i+6).y, MF_U_high(i+6).y);
}


ap::complex_2d_array MF_A, MF_A_MIN, MF_A_MAX, MF_A_MIN_INV, MF_A_MAX_INV;
MF_A.setlength(3,3);
MF_A_MIN.setlength(3,3);
MF_A_MAX.setlength(3,3);
MF_A_MIN_INV.setlength(3,3);
MF_A_MAX_INV.setlength(3,3);
MF_A(0,0).x = (1/(2*pi*MF_E0)) * _double(ln(_real((2*MF_Y1)/MF_R)));
MF_A(0,1).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X1 - MF_X2, 2.0) + pow(MF_Y1 + MF_Y2, 2.0)) )  /  sqrt( (pow(MF_X1 - MF_X2, 2.0) + pow(MF_Y1 - MF_Y2, 2.0)) ))));
MF_A(0,2).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X1 - MF_X3, 2.0) + pow(MF_Y1 + MF_Y3, 2.0)) )  /  sqrt( (pow(MF_X1 - MF_X3, 2.0) + pow(MF_Y1 - MF_Y3, 2.0)) ))));
MF_A(1,0).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X2 - MF_X1, 2.0) + pow(MF_Y2 + MF_Y1, 2.0)) )  /  sqrt( (pow(MF_X2 - MF_X1, 2.0) + pow(MF_Y2 - MF_Y1, 2.0)) ))));
MF_A(1,1).x = (1/(2*pi*MF_E0)) * _double(ln(_real((2*MF_Y2)/MF_R)));
MF_A(1,2).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X2 - MF_X3, 2.0) + pow(MF_Y2 + MF_Y3, 2.0)) )  /  sqrt( (pow(MF_X2 - MF_X3, 2.0) + pow(MF_Y2 - MF_Y3, 2.0)) ))));
MF_A(2,0).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X3 - MF_X1, 2.0) + pow(MF_Y3 + MF_Y1, 2.0)) )  /  sqrt( (pow(MF_X3 - MF_X1, 2.0) + pow(MF_Y3 - MF_Y1, 2.0)) ))));
MF_A(2,1).x = (1/(2*pi*MF_E0)) * _double(ln(_real(sqrt( (pow(MF_X3 - MF_X2, 2.0) + pow(MF_Y3 + MF_Y2, 2.0)) )  /  sqrt( (pow(MF_X3 - MF_X2, 2.0) + pow(MF_Y3 - MF_Y2, 2.0)) ))));
MF_A(2,2).x = (1/(2*pi*MF_E0)) * _double(ln(_real((2*MF_Y3)/MF_R)));

// INTERVAL CASE
cinterval iMF_A[3][3], iMF_PI;
Re(iMF_PI) = Pi();
Im(iMF_PI) = NULL;

iMF_A[0][0] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln((2.0*iMF_Y[1])/iMF_R);
iMF_A[0][1] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[1] - iMF_X[2])*(iMF_X[1] - iMF_X[2])  + (iMF_Y[1] + iMF_Y[2])*(iMF_Y[1] + iMF_Y[2]) ) )  /  sqrt( (iMF_X[1] - iMF_X[2])*(iMF_X[1] - iMF_X[2]) + (iMF_Y[1] - iMF_Y[2])*(iMF_Y[1] - iMF_Y[2])) );
iMF_A[0][2] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[1] - iMF_X[3])*(iMF_X[1] - iMF_X[3])  + (iMF_Y[1] + iMF_Y[3])*(iMF_Y[1] + iMF_Y[3]) ) )  /  sqrt( (iMF_X[1] - iMF_X[3])*(iMF_X[1] - iMF_X[3]) + (iMF_Y[1] - iMF_Y[3])*(iMF_Y[1] - iMF_Y[3])) );
iMF_A[1][0] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[2] - iMF_X[1])*(iMF_X[2] - iMF_X[1])  + (iMF_Y[2] + iMF_Y[1])*(iMF_Y[2] + iMF_Y[1]) ) )  /  sqrt( (iMF_X[2] - iMF_X[1])*(iMF_X[2] - iMF_X[1]) + (iMF_Y[2] - iMF_Y[1])*(iMF_Y[2] - iMF_Y[1])) );
iMF_A[1][1] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln((2.0*iMF_Y[2])/iMF_R);
iMF_A[1][2] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[2] - iMF_X[3])*(iMF_X[2] - iMF_X[3])  + (iMF_Y[2] + iMF_Y[3])*(iMF_Y[2] + iMF_Y[3]) ) )  /  sqrt( (iMF_X[2] - iMF_X[3])*(iMF_X[2] - iMF_X[3]) + (iMF_Y[2] - iMF_Y[3])*(iMF_Y[2] - iMF_Y[3])) );
iMF_A[2][0] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[3] - iMF_X[1])*(iMF_X[3] - iMF_X[1])  + (iMF_Y[3] + iMF_Y[1])*(iMF_Y[3] + iMF_Y[1]) ) )  /  sqrt( (iMF_X[3] - iMF_X[1])*(iMF_X[3] - iMF_X[1]) + (iMF_Y[3] - iMF_Y[1])*(iMF_Y[3] - iMF_Y[1])) );
iMF_A[2][1] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln(sqrt( ( (iMF_X[3] - iMF_X[2])*(iMF_X[3] - iMF_X[2])  + (iMF_Y[3] + iMF_Y[2])*(iMF_Y[3] + iMF_Y[2]) ) )  /  sqrt( (iMF_X[3] - iMF_X[2])*(iMF_X[3] - iMF_X[2]) + (iMF_Y[3] - iMF_Y[2])*(iMF_Y[3] - iMF_Y[2])) );
iMF_A[2][2] = (1.0/(2.0*iMF_PI*iMF_E0)) * ln((2.0*iMF_Y[3])/iMF_R);
//cout << iMF_A[0][1] << "\n";
// END OF INTERVAL CASE

/*
for (int i = 0; i < 3; i++)     
{
	for (int j = 0; j < 3; j++) 
	{
	cout << "MF_A(" << i+1 << "," << j+1 << ").X = " << MF_A(i,j).x << "\n";
	}
}
*/

for (int i = 0; i < 3; i++)     
{
	for (int j = 0; j < 3; j++) 
	{
	MF_A_MIN(i,j).x = _double(Inf(Re(iMF_A[i][j])));
	MF_A_MIN(i,j).y = _double(Inf(Im(iMF_A[i][j])));
	MF_A_MAX(i,j).x = _double(Sup(Re(iMF_A[i][j])));
	MF_A_MAX(i,j).y = _double(Sup(Im(iMF_A[i][j])));
	}
}

cmatrixinverse(MF_A,3);
cmatrixinverse(MF_A_MIN,3);
cmatrixinverse(MF_A_MAX,3);

for (int i = 0; i < 3; i++)     
{
	for (int j = 0; j < 3; j++) 
	{
	MF_A_MIN_INV(i,j).x = min(MF_A_MIN(i,j).x, MF_A_MAX(i,j).x);
	MF_A_MIN_INV(i,j).y = min(MF_A_MIN(i,j).y, MF_A_MAX(i,j).y);
	MF_A_MAX_INV(i,j).x = max(MF_A_MIN(i,j).x, MF_A_MAX(i,j).x);
	MF_A_MAX_INV(i,j).y = max(MF_A_MIN(i,j).y, MF_A_MAX(i,j).y);
	}
}

for (int i = 0; i < 3; i++)     
{
	for (int j = 0; j < 3; j++) 
	{
	MF_A_MIN(i,j).x = MF_A_MIN_INV(i,j).x;
	MF_A_MIN(i,j).y = MF_A_MIN_INV(i,j).y;
	MF_A_MAX(i,j).x = MF_A_MAX_INV(i,j).x;
	MF_A_MAX(i,j).y = MF_A_MAX_INV(i,j).y;
	}
}
//cout << "A(0,0) MIN = [" << MF_A_MIN(0,0).x << "]+j*[" << MF_A_MIN(0,0).y << "]\n";
//cout << "A(0,0) MAX = [" << MF_A_MAX(0,0).x << "]+j*[" << MF_A_MAX(0,0).y << "]\n";

cinterval iMF_Ex, iMF_Ey, iMF_Hx, iMF_Hy, iMF_E_ef, iMF_H_ef, iMF_T[3][1], iMF_I[3][1];

// T [3,1] = a_inv [3,3] * U [3,1]

for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 MF_T_low(i,j) += MF_A_MIN(i,k) * MF_VU_1(k,j);
		 MF_T_high(i,j) += MF_A_MAX(i,k) * MF_VU_2(k,j);
	  }
   }
}

for (int i = 0; i < 3; i++)     
{
Inf(Re(iMF_T[i][0])) = min(MF_T_low(i,0).x, MF_T_high(i,0).x);
Inf(Im(iMF_T[i][0])) = min(MF_T_low(i,0).y, MF_T_high(i,0).y);
Sup(Re(iMF_T[i][0])) = max(MF_T_low(i,0).x, MF_T_high(i,0).x);
Sup(Im(iMF_T[i][0])) = max(MF_T_low(i,0).y, MF_T_high(i,0).y);

Inf(Re(iMF_I[i][0])) = min(MF_I_1(i,0).x, MF_I_2(i,0).x);
Inf(Im(iMF_I[i][0])) = min(MF_I_1(i,0).y, MF_I_2(i,0).y);
Sup(Re(iMF_I[i][0])) = max(MF_I_1(i,0).x, MF_I_2(i,0).x);
Sup(Im(iMF_I[i][0])) = max(MF_I_1(i,0).y, MF_I_2(i,0).y);
}

/*
for (int i = 0; i < 3; i++)     
{
cout << "MF_T(" << i+1 << ") = [" << MF_T_low(i,0).x << "+j*[" << MF_T_low(i,0).y << "],[" << MF_T_high(i,0).x << "+j*[" << MF_T_high(i,0).y << "]" << "\n";
}
*/

iMF_Ex = ((2.0)/(Pi()*iMF_E0)) * ( iMF_T[0][0] * ( (iMF_X[0]-iMF_X[1])*iMF_Y[0]*iMF_Y[1]) / ( ((iMF_X[0]-iMF_X[1])*(iMF_X[0]-iMF_X[1]) + (iMF_Y[0]+iMF_Y[1])*(iMF_Y[0]+iMF_Y[1]) ) * ((iMF_X[0]-iMF_X[1])*(iMF_X[0]-iMF_X[1]) + (iMF_Y[0]-iMF_Y[1])*(iMF_Y[0]-iMF_Y[1])) )  
							   +   iMF_T[1][0] * ( (iMF_X[0]-iMF_X[2])*iMF_Y[0]*iMF_Y[2]) / ( ((iMF_X[0]-iMF_X[2])*(iMF_X[0]-iMF_X[2]) + (iMF_Y[0]+iMF_Y[2])*(iMF_Y[0]+iMF_Y[2]) ) * ((iMF_X[0]-iMF_X[2])*(iMF_X[0]-iMF_X[2]) + (iMF_Y[0]-iMF_Y[2])*(iMF_Y[0]-iMF_Y[2])) )  
		                       +   iMF_T[2][0] * ( (iMF_X[0]-iMF_X[3])*iMF_Y[0]*iMF_Y[3]) / ( ((iMF_X[0]-iMF_X[3])*(iMF_X[0]-iMF_X[3]) + (iMF_Y[0]+iMF_Y[3])*(iMF_Y[0]+iMF_Y[3]) ) * ((iMF_X[0]-iMF_X[3])*(iMF_X[0]-iMF_X[3]) + (iMF_Y[0]-iMF_Y[3])*(iMF_Y[0]-iMF_Y[3])) ) );

iMF_Ey = ((- 1.0)/(Pi()*iMF_E0)) * (   iMF_T[0][0] * ( (iMF_Y[1]*((iMF_X[0]-iMF_X[1])*(iMF_X[0]-iMF_X[1])-iMF_Y[0]*iMF_Y[0]+iMF_Y[1]*iMF_Y[1])) / ( ((iMF_X[0]-iMF_X[1])*(iMF_X[0]-iMF_X[1]) + (iMF_Y[0]+iMF_Y[1])*(iMF_Y[0]+iMF_Y[1]) ) * ((iMF_X[0]-iMF_X[1])*(iMF_X[0]-iMF_X[1]) + (iMF_Y[0]-iMF_Y[1])*(iMF_Y[0]-iMF_Y[1]))) ) 
		                             + iMF_T[1][0] * ( (iMF_Y[2]*((iMF_X[0]-iMF_X[2])*(iMF_X[0]-iMF_X[2])-iMF_Y[0]*iMF_Y[0]+iMF_Y[2]*iMF_Y[2])) / ( ((iMF_X[0]-iMF_X[2])*(iMF_X[0]-iMF_X[2]) + (iMF_Y[0]+iMF_Y[2])*(iMF_Y[0]+iMF_Y[2]) ) * ((iMF_X[0]-iMF_X[2])*(iMF_X[0]-iMF_X[2]) + (iMF_Y[0]-iMF_Y[2])*(iMF_Y[0]-iMF_Y[2]))) ) 
									 + iMF_T[2][0] * ( (iMF_Y[3]*((iMF_X[0]-iMF_X[3])*(iMF_X[0]-iMF_X[3])-iMF_Y[0]*iMF_Y[0]+iMF_Y[3]*iMF_Y[3])) / ( ((iMF_X[0]-iMF_X[3])*(iMF_X[0]-iMF_X[3]) + (iMF_Y[0]+iMF_Y[3])*(iMF_Y[0]+iMF_Y[3]) ) * ((iMF_X[0]-iMF_X[3])*(iMF_X[0]-iMF_X[3]) + (iMF_Y[0]-iMF_Y[3])*(iMF_Y[0]-iMF_Y[3]))) ) );


iMF_Hx = (1.0 / (2.0 * Pi())) *  (	iMF_I[0][0] * ( ( iMF_Y[0]-iMF_Y[1] ) / ( (iMF_X[1]-iMF_X[0])*(iMF_X[1]-iMF_X[0])+(iMF_Y[1]-iMF_Y[0])*(iMF_Y[1]-iMF_Y[0]) ) )  
							 +	iMF_I[1][0] * ( ( iMF_Y[0]-iMF_Y[2] ) / ( (iMF_X[2]-iMF_X[0])*(iMF_X[2]-iMF_X[0])+(iMF_Y[2]-iMF_Y[0])*(iMF_Y[2]-iMF_Y[0]) ) )  
							 +	iMF_I[2][0] * ( ( iMF_Y[0]-iMF_Y[3] ) / ( (iMF_X[3]-iMF_X[0])*(iMF_X[3]-iMF_X[0])+(iMF_Y[3]-iMF_Y[0])*(iMF_Y[3]-iMF_Y[0]) ) ) );

iMF_Hy = (-1.0/ (2.0 * Pi())) *  (	iMF_I[0][0] * ( ( iMF_X[0]-iMF_X[1] ) / ( (iMF_X[1]-iMF_X[0])*(iMF_X[1]-iMF_X[0])+(iMF_Y[1]-iMF_Y[0])*(iMF_Y[1]-iMF_Y[0]) ) )  
							 +	iMF_I[1][0] * ( ( iMF_X[0]-iMF_X[2] ) / ( (iMF_X[2]-iMF_X[0])*(iMF_X[2]-iMF_X[0])+(iMF_Y[2]-iMF_Y[0])*(iMF_Y[2]-iMF_Y[0]) ) )  
							 +	iMF_I[2][0] * ( ( iMF_X[0]-iMF_X[3] ) / ( (iMF_X[3]-iMF_X[0])*(iMF_X[3]-iMF_X[0])+(iMF_Y[3]-iMF_Y[0])*(iMF_Y[3]-iMF_Y[0]) ) ) );

/*
iMF_E_ef = sqrt (
			sqrt ( Max(Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex))) * Max(Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex))) + Max(Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex))) * Max(Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex))) ) *
			sqrt ( Max(Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex))) * Max(Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex))) + Max(Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex))) * Max(Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex))) ) +
			sqrt ( Max(Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey))) * Max(Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey))) + Max(Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey))) * Max(Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey))) ) *
			sqrt ( Max(Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey))) * Max(Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey))) + Max(Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey))) * Max(Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey))) ) 
			);

iMF_H_ef = sqrt (
			sqrt ( Max(Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx))) * Max(Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx))) + Max(Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx))) * Max(Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx))) ) *
			sqrt ( Max(Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx))) * Max(Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx))) + Max(Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx))) * Max(Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx))) ) +
			sqrt ( Max(Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy))) * Max(Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy))) + Max(Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy))) * Max(Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy))) ) *
			sqrt ( Max(Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy))) * Max(Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy))) + Max(Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy))) * Max(Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy))) ) 
			);
*/

cinterval iMF_Ex_2, iMF_Ey_2, iMF_Hx_2, iMF_Hy_2;

Inf(Re(iMF_Ex_2)) = min (Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex)));
Inf(Im(iMF_Ex_2)) = min (Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex)));
Sup(Re(iMF_Ex_2)) = max (Inf(Re(iMF_Ex)),Sup(Re(iMF_Ex)));
Sup(Im(iMF_Ex_2)) = max (Inf(Im(iMF_Ex)),Sup(Im(iMF_Ex)));

Inf(Re(iMF_Ey_2)) = min (Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey)));
Inf(Im(iMF_Ey_2)) = min (Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey)));
Sup(Re(iMF_Ey_2)) = max (Inf(Re(iMF_Ey)),Sup(Re(iMF_Ey)));
Sup(Im(iMF_Ey_2)) = max (Inf(Im(iMF_Ey)),Sup(Im(iMF_Ey)));

Inf(Re(iMF_Hx_2)) = min (Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx)));
Inf(Im(iMF_Hx_2)) = min (Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx)));
Sup(Re(iMF_Hx_2)) = max (Inf(Re(iMF_Hx)),Sup(Re(iMF_Hx)));
Sup(Im(iMF_Hx_2)) = max (Inf(Im(iMF_Hx)),Sup(Im(iMF_Hx)));

Inf(Re(iMF_Hy_2)) = min (Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy)));
Inf(Im(iMF_Hy_2)) = min (Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy)));
Sup(Re(iMF_Hy_2)) = max (Inf(Re(iMF_Hy)),Sup(Re(iMF_Hy)));
Sup(Im(iMF_Hy_2)) = max (Inf(Im(iMF_Hy)),Sup(Im(iMF_Hy)));


//iMF_E_ef = sqrt ( abs(iMF_Ex_2) * abs(iMF_Ex_2) + abs(iMF_Ey_2) * abs(iMF_Ey_2) );
//iMF_H_ef = sqrt ( abs(iMF_Hx_2) * abs(iMF_Hx_2) + abs(iMF_Hy_2) * abs(iMF_Hy_2) );

iMF_E_ef = sqr(iMF_Ex_2 * iMF_Ex_2  + iMF_Ey_2 * iMF_Ey_2);
//cout << "iMF_Ex_2 = " << iMF_Ex_2 << "\n";
//cout << "iMF_Ey_2 = " << iMF_Ey_2 << "\n";
//cout << "iMF_E_ef = " << iMF_E_ef << "\n";

iMF_H_ef = sqr(iMF_Hx_2 * iMF_Hx_2  + iMF_Hy_2 * iMF_Hy_2);
//cout << "iMF_H_ef = " << iMF_H_ef << "\n";


//cout << "Hx = " << iMF_Hx << "\n";
//cout << "Hy = " << iMF_Hy << "\n";
//cout << "Hef = " << iMF_H_ef << "\n";
//cout << "Hx2 = " << iMF_Hx_2 << "\n";
//cout << "Hy2 = " << iMF_Hy_2 << "\n";

double MF_Ex_low_mod, MF_Ex_high_mod, MF_Ey_low_mod, MF_Ey_high_mod, MF_Hx_low_mod, MF_Hx_high_mod, MF_Hy_low_mod, MF_Hy_high_mod;

MF_Ex_low_mod  = _double(sqrt(Inf(Re(iMF_Ex_2))*Inf(Re(iMF_Ex_2)) + Inf(Im(iMF_Ex_2))*Inf(Im(iMF_Ex_2))));
MF_Ex_high_mod = _double(sqrt(Sup(Re(iMF_Ex_2))*Sup(Re(iMF_Ex_2)) + Sup(Im(iMF_Ex_2))*Sup(Im(iMF_Ex_2))));
MF_Ey_low_mod  = _double(sqrt(Inf(Re(iMF_Ey_2))*Inf(Re(iMF_Ey_2)) + Inf(Im(iMF_Ey_2))*Inf(Im(iMF_Ey_2))));
MF_Ey_high_mod = _double(sqrt(Sup(Re(iMF_Ey_2))*Sup(Re(iMF_Ey_2)) + Sup(Im(iMF_Ey_2))*Sup(Im(iMF_Ey_2))));

MF_Hx_low_mod  = _double(sqrt(Inf(Re(iMF_Hx_2))*Inf(Re(iMF_Hx_2)) + Inf(Im(iMF_Hx_2))*Inf(Im(iMF_Hx_2))));
MF_Hx_high_mod = _double(sqrt(Sup(Re(iMF_Hx_2))*Sup(Re(iMF_Hx_2)) + Sup(Im(iMF_Hx_2))*Sup(Im(iMF_Hx_2))));
MF_Hy_low_mod  = _double(sqrt(Inf(Re(iMF_Hy_2))*Inf(Re(iMF_Hy_2)) + Inf(Im(iMF_Hy_2))*Inf(Im(iMF_Hy_2))));
MF_Hy_high_mod = _double(sqrt(Sup(Re(iMF_Hy_2))*Sup(Re(iMF_Hy_2)) + Sup(Im(iMF_Hy_2))*Sup(Im(iMF_Hy_2))));

/*
cout << "Ex = " << iMF_Ex << "\n";
cout << "Ey = " << iMF_Ey << "\n";
cout << "Hx = " << iMF_Hx << "\n";
cout << "Hy = " << iMF_Hy << "\n";
*/

iGIG_A = iMF_Ex_2;
get_interval_grad();
interval iMF_Ex_grad = iGIG_GRAD;

iGIG_A = iMF_Ey_2;
get_interval_grad();
interval iMF_Ey_grad = iGIG_GRAD;

iGIG_A = iMF_Hx_2;
get_interval_grad();
interval iMF_Hx_grad = iGIG_GRAD;

iGIG_A = iMF_Hy_2;
get_interval_grad();
interval iMF_Hy_grad = iGIG_GRAD;

cout << "MIN |Ex| = " << MF_Ex_low_mod << "\n";
cout << "MAX |Ex| = " << MF_Ex_high_mod << "\n";
cout << "Ex grad = " << iMF_Ex_grad << "\n";
cout << "MIN |Ey| = " << MF_Ey_low_mod << "\n";
cout << "MAX |Ey| = " << MF_Ey_high_mod << "\n";
cout << "Ey grad = " << iMF_Ey_grad << "\n";
cout << "E effective = " << abs(iMF_E_ef) << "\n";

cout << "MIN |Hx| = " << MF_Hx_low_mod << "\n";
cout << "MAX |Hx| = " << MF_Hx_high_mod << "\n";
cout << "Hx grad = " << iMF_Hx_grad << "\n";
cout << "MIN |Hy| = " << MF_Hy_low_mod << "\n";
cout << "MAX |Hy| = " << MF_Hy_high_mod << "\n";
cout << "Hy grad = " << iMF_Hy_grad << "\n";
cout << "H effective = " << abs(iMF_H_ef) << "\n";

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_dls_cramer()
{
//
cout << "Решение СЛАУ методом Крамера\n";
interval iDLS_Cramer_x1, iDLS_Cramer_x2;
interval iDLS_Cramer_a11, iDLS_Cramer_a12, iDLS_Cramer_a21, iDLS_Cramer_a22, iDLS_Cramer_b1, iDLS_Cramer_b2;
interval iDLS_Cramer_D, iDLS_Cramer_D1, iDLS_Cramer_D2;
Inf(iDLS_Cramer_a11) = 2.0;
Sup(iDLS_Cramer_a11) = 4.0;
Inf(iDLS_Cramer_a12) = - 2.0;
Sup(iDLS_Cramer_a12) = 1.0;
Inf(iDLS_Cramer_a21) = - 1.0;
Sup(iDLS_Cramer_a21) = 2.0;
Inf(iDLS_Cramer_a22) = 2.0;
Sup(iDLS_Cramer_a22) = 4.0;
Inf(iDLS_Cramer_b1) = - 2.0;
Sup(iDLS_Cramer_b1) = 2.0;
Inf(iDLS_Cramer_b2) = - 2.0;
Sup(iDLS_Cramer_b2) = 2.0;
iDLS_Cramer_D = iDLS_Cramer_a11 * iDLS_Cramer_a22 - iDLS_Cramer_a12 * iDLS_Cramer_a21;
iDLS_Cramer_D1 = iDLS_Cramer_b1 * iDLS_Cramer_a22 - iDLS_Cramer_b2 * iDLS_Cramer_a21;
iDLS_Cramer_D2 = iDLS_Cramer_a11 * iDLS_Cramer_b2 - iDLS_Cramer_a12 * iDLS_Cramer_b1;
iDLS_Cramer_x1 = iDLS_Cramer_D1 / iDLS_Cramer_D;
iDLS_Cramer_x2 = iDLS_Cramer_D2 / iDLS_Cramer_D;

cout << "X1 = " << iDLS_Cramer_x1 << "\n";
cout << "X2 = " << iDLS_Cramer_x2 << "\n";
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_shortcircuit_mis()
{
cinterval iSCM_Y1, iSCM_U, iSCM_YB1, iSCM_UB;



}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_shortcircuit_1f_fcoor_ground()
{

// Расчет токов короткого замыкания. Двухфазное КЗ. Расчет в фазных координатах.
SC_UL = 110.0;
SC_UL_1 = 110.0;
SC_UL_2 = 110.0;
SC_L = 100.0;
SC_L_1 = 97.5;
SC_L_2 = 102.5;
SC_X0 = 0.4;
SC_X0_1 = 0.38; //0.32;
SC_X0_2 = 0.41; //0.48;
interval iSC_XL, iSC_XM, iSC_35;
cinterval iSC_ZL, iSC_ZM, iSC_C1;
cinterval iSC_Z1[3][3], iSC_Z2[3][3], iSC_Y1[3][3], iSC_Y2[3][3];
cinterval iSC_YC1[6][6], iSC_YC2[6][6], iSC_O[6][6], iSC_YB[12][12], iSC_YS1[9][9], iSC_M[9][12], iSC_M_trans[12][9];
interval iSC_P[3][6], iSC_O1[3][6];

ap::complex_2d_array apSC_Z1_inf;
ap::complex_2d_array apSC_Z1_sup;
apSC_Z1_inf.setlength(3,3);
apSC_Z1_sup.setlength(3,3);

ap::complex_2d_array apSC_Z2_inf;
ap::complex_2d_array apSC_Z2_sup;
apSC_Z2_inf.setlength(3,3);
apSC_Z2_sup.setlength(3,3);
//===============================================================================================================
Inf(iSC_UL) = _real(SC_UL_1);
Sup(iSC_UL) = _real(SC_UL_2);
Inf(iSC_L) = _real(SC_L_1);
Sup(iSC_L) = _real(SC_L_2);
cout << "L = " << iSC_L << "км \n";
Inf(iSC_X0) = _real(SC_X0_1);
Sup(iSC_X0) = _real(SC_X0_2);
Inf(iSC_35) = _real(3.5);
Sup(iSC_35) = _real(3.5);
iSC_UF = iSC_UL / sqrt(3.0);
iSC_X1 = iSC_X0 * iSC_L;
iSC_X0 = iSC_35 * iSC_X1;
iSC_XL = (iSC_X1 / 1.5) + (iSC_X0 / 3.0);
iSC_XM = iSC_XL - iSC_X1;
Re(iSC_ZL) = 0.0;
Im(iSC_ZL) = iSC_XL;
Re(iSC_ZM) = 0.0;
Im(iSC_ZM) = iSC_XM;
// Проверка
iSC_C1 = iSC_XL - iSC_XM;
//===============================================================================================================
iSC_Z1[0][0] = iSC_ZL;
iSC_Z1[0][1] = iSC_ZM;
iSC_Z1[0][2] = iSC_ZM;
iSC_Z1[1][0] = iSC_ZM;
iSC_Z1[1][1] = iSC_ZL;
iSC_Z1[1][2] = iSC_ZM;
iSC_Z1[2][0] = iSC_ZM;
iSC_Z1[2][1] = iSC_ZM;
iSC_Z1[2][2] = iSC_ZL;

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z1_inf(i,j).x = _double(Inf(Re(iSC_Z1[i][j])));
	apSC_Z1_inf(i,j).y = _double(Inf(Im(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).x = _double(Sup(Re(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).y = _double(Sup(Im(iSC_Z1[i][j])));
	}
}

cmatrixinverse(apSC_Z1_inf,3);
cmatrixinverse(apSC_Z1_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).x;
	Inf(Im(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).y;
	Sup(Re(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).x;
	Sup(Im(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).y;
	}
}
//===============================================================================================================
/*// Двухфазное КЗ (между фазами)
Re(iSC_Z2[0][0]) = 0.0;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 10000.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 10000.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 10000.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 0.0;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.001;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 10000.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.001;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 0.0;
Im(iSC_Z2[2][2]) = 0.0;
*/

// Однофазное КЗ на землю
Re(iSC_Z2[0][0]) = 0.000001;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 0.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 0.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 0.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 10000000.0;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.0;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 0.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.0;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 10000000.0;
Im(iSC_Z2[2][2]) = 0.0;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z2_inf(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	apSC_Z2_inf(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
	}
}

cmatrixinverse(apSC_Z2_inf,3);
cmatrixinverse(apSC_Z2_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).x;
	Inf(Im(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).y;
	Sup(Re(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).x;
	Sup(Im(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).y;
	}
}

//===============================================================================================================
// Матрица YC1
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i][j]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i][j]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i][j]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i][j-3]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i][j-3]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i][j-3]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i-3][j]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i-3][j]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i-3][j]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i-3][j-3]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i-3][j-3]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i-3][j-3]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i-3][j-3]));
	}
}
//===============================================================================================================
// Матрица YC2
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i][j]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i][j]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i][j]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i][j-3]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i][j-3]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i][j-3]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i-3][j]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i-3][j]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i-3][j]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i-3][j-3]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i-3][j-3]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i-3][j-3]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i-3][j-3]));
	}
}

//===============================================================================================================
// Матрица P
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_P[i][j] = 0.0;
	}
}
iSC_P[0][0] = 1.0;
iSC_P[1][1] = 1.0;
iSC_P[2][2] = 1.0;
//===============================================================================================================
// Матрица O1
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O1[i][j] = 0.0;
	}
}
//===============================================================================================================
// Матрица M
for (int i=0;i<9;i++)
{
	for (int j=0;j<12;j++)
	{
	Re(iSC_M[i][j]) = 0.0;
	Im(iSC_M[i][j]) = 0.0;
	}
}
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i][j];
	}
}
for (int i=0;i<3;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i][j-6];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-3][j];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i-3][j-6];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i-6][j];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-6][j-6];
	}
}

//===============================================================================================================
//Матрица O

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O[i][j] = 0.0;
	}
}
//===============================================================================================================
//Матрица YB
for (int i=0;i<12;i++)
{
	for (int j=0;j<12;j++)
	{
	iSC_YB[i][j] = 0.0;
	}
}
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_YC1[i][j];
	}
}
for (int i=0;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_O[i][j-6];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_O[i-6][j];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_YC2[i-6][j-6];
	}
}
//===============================================================================================================
// M * Yb = (9x12) * (12x12)
cinterval iSC_Ytemp[9][12];
ap::complex_2d_array Y_1_shm, Y_2_shm, Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm.setlength(9,9);
Y_2_shm.setlength(9,9);
Y_1_shm_temp.setlength(9,12);
Y_2_shm_temp.setlength(9,12);
Y_1_shm_temp_reint.setlength(9,12);
Y_2_shm_temp_reint.setlength(9,12);
ap::complex_2d_array Yb_1_shm, Yb_2_shm;
Yb_1_shm.setlength(12,12);
Yb_2_shm.setlength(12,12);
ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
ap::complex_2d_array M_shm;
M_shm_transposed.setlength(12,9);
M_shm.setlength(9,12);

for (int i = 0; i < 9; i++)     
{
   for (int j = 0; j < 12; j++) 
   {
	M_shm(i,j).x = _double(Inf(Re(iSC_M[i][j])));
	M_shm(i,j).y = _double(0.0);
   }
}
for (int i = 0; i < 12; i++)      
{
   for (int j = 0; j < 12; j++)  
   {
	Yb_1_shm(i,j).x = _double(Inf(Re(iSC_YB[i][j])));
	Yb_1_shm(i,j).y = _double(Inf(Im(iSC_YB[i][j])));
	Yb_2_shm(i,j).x = _double(Sup(Re(iSC_YB[i][j])));
	Yb_2_shm(i,j).y = _double(Sup(Im(iSC_YB[i][j])));
   }
}

for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {
		Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}


for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}


// Transpose M
for (int i=0; i < 9; i++)
	{
	for (int j=0; j < 12 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}
	}



// Y_temp * M^transpose = (9x12) * (12x9)
for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 9; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {

		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}
ap::complex_2d_array SC_Y1_1, SC_Y1_2, SC_Y2_1, SC_Y2_2, SC_Y12_1, SC_Y12_2, SC_Y21_1, SC_Y21_2;
SC_Y1_1.setlength(6,6);
SC_Y1_2.setlength(6,6);
SC_Y2_1.setlength(3,3);
SC_Y2_2.setlength(3,3);
SC_Y12_1.setlength(6,3);
SC_Y12_2.setlength(6,3);
SC_Y21_1.setlength(3,6);
SC_Y21_2.setlength(3,6);
for (int i = 0; i < 6; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y1_1(i,j) = Y_1_shm(i,j);
    SC_Y1_2(i,j) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y2_1(i-6,j-6) = Y_1_shm(i,j);
    SC_Y2_2(i-6,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 0; i < 6; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y12_1(i,j-6) = Y_1_shm(i,j);
    SC_Y12_2(i,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y21_1(i-6,j) = Y_1_shm(i,j);
    SC_Y21_2(i-6,j) = Y_2_shm(i,j);
    }
}
ap::complex_2d_array SC_e1, SC_e1_T, SC_Y12_1_T, SC_Y12_2_T, SC_C1_1, SC_C1_2, SC_C2_1, SC_C2_2, SC_C3_1, SC_C3_2;
ap::complex_2d_array SC_S1_1, SC_S1_2, SC_S2_1, SC_S2_2, SC_YE_1, SC_YE_2;
SC_S1_1.setlength(7,1);
SC_S1_2.setlength(7,1);
SC_S2_1.setlength(7,6);
SC_S2_2.setlength(7,6);
SC_YE_1.setlength(7,7);
SC_YE_2.setlength(7,7);
SC_e1.setlength(3,1);
SC_e1_T.setlength(1,3);
SC_C1_1.setlength(6,1);
SC_C1_2.setlength(6,1);
SC_C2_1.setlength(1,6);
SC_C2_2.setlength(1,6);
SC_C3_1.setlength(1,1);
SC_C3_2.setlength(1,1);
SC_Y12_1_T.setlength(3,6);
SC_Y12_2_T.setlength(3,6);
for (int i = 0; i < 3; i++)
{
	for (int j = 0; j < 1; j++)
	{
	SC_e1(i,j).x = 1.0;
    SC_e1(i,j).y = 0.0;
    }
}

// C1 = Y12 * e1
for (int i = 0; i < 6; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C1_1(i,j) += SC_Y12_1(i,k) * SC_e1(k,j);
		 SC_C1_2(i,j) += SC_Y12_2(i,k) * SC_e1(k,j);
	  }
   }
}

// Transpose e1
for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1 ; j++)
		{
		M_element = SC_e1(i,j);
		SC_e1_T(j,i) = M_element;
		}
	}

// Transpose Y12
ap::complex ELEM_1, ELEM_2;
for (int i=0; i < 6; i++)
	{
	for (int j=0; j < 3 ; j++)
		{
		ELEM_1 = SC_Y12_1(i,j);
		ELEM_2 = SC_Y12_2(i,j);
		SC_Y12_1_T(j,i) = ELEM_1;
		SC_Y12_2_T(j,i) = ELEM_2;
		}
	}

// C2 = e1^T * Y12^T
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 6; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C2_1(i,j) += SC_e1_T(i,k) * SC_Y12_1_T(k,j);
		 SC_C2_2(i,j) += SC_e1_T(i,k) * SC_Y12_2_T(k,j);
	  }
   }
}
// C3 = e1^T * Y2 * e1
ap::complex_2d_array SC_TEMP_1, SC_TEMP_2;
SC_TEMP_1.setlength(1,3);
SC_TEMP_2.setlength(1,3);
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 3; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_TEMP_1(i,j) += SC_e1_T(i,k) * SC_Y2_1(k,j);
		 SC_TEMP_2(i,j) += SC_e1_T(i,k) * SC_Y2_2(k,j);
	  }
   }
}
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C3_1(i,j) += SC_TEMP_1(i,k) * SC_e1(k,j);
		 SC_C3_2(i,j) += SC_TEMP_2(i,k) * SC_e1(k,j);
	  }
   }
}
//===============================================================================================================
// S1
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C1_1(i,j);
	SC_S1_2(i,j) = SC_C1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C3_1(i-6,j);
	SC_S1_2(i,j) = SC_C3_2(i-6,j);
   }
}
//===============================================================================================================
// S2 = STACK(Y1, C2) <=> 6x6 & 1x6
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_Y1_1(i,j);
	SC_S2_2(i,j) = SC_Y1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_C2_1(i-6,j);
	SC_S2_2(i,j) = SC_C2_2(i-6,j);
   }
}
//===============================================================================================================
// YE = augment(S2, S1) = 7x6 & 7x1
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_YE_1(i,j) = SC_S2_1(i,j);
	SC_YE_2(i,j) = SC_S2_2(i,j);
   }
}
for (int i = 0; i < 7; i++)    
{
   for (int j = 6; j < 7; j++)  
   {
	SC_YE_1(i,j) = SC_S1_1(i,j-6);
	SC_YE_2(i,j) = SC_S1_2(i,j-6);
   }
}
/*
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 7; j++)  
   {
	cout << "YE.x ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).x << "\n";
	cout << "YE.y ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).y << "\n";
   }
}
*/
//===============================================================================================================
// UBE
cinterval iSC_UBE[3][1];
ap::complex_2d_array SC_UBE_1, SC_UBE_2;
SC_UBE_1.setlength(3,1);
SC_UBE_2.setlength(3,1);
cinterval SC_2_3_of_Pi;
Re(SC_2_3_of_Pi) = 0.0;
Im(SC_2_3_of_Pi) = double((2.0 / 3.0) * pi);

if (param_unsim == 0)
{
// then 0
iSC_UBE[0][0] = 1.0;
iSC_UBE[1][0] = exp(- SC_2_3_of_Pi);
iSC_UBE[2][0] = exp(SC_2_3_of_Pi);
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	   iSC_UBE[i][j] = iSC_UBE[i][j] * iSC_UF;
   }
}
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_UBE_1(i,j).x = _double(Inf(Re(iSC_UBE[i][j])));
	SC_UBE_1(i,j).y = _double(Inf(Im(iSC_UBE[i][j])));
	SC_UBE_2(i,j).x = _double(Sup(Re(iSC_UBE[i][j])));
	SC_UBE_2(i,j).y = _double(Sup(Im(iSC_UBE[i][j])));
   }
}
//end 0
}
else if (param_unsim == 1)
{
// then 1
get_unsimmetric_fcoor();
	for (int i = 0; i < 3; i++)    
	{
		for (int j = 0; j < 1; j++)  
		{
		SC_UBE_1(i,j) = ISC_U_unsim(i,j);
		SC_UBE_2(i,j) = ISC_U_unsim(i,j);
		}
}

}

//cout << "y = " << iSC_UBE[1][0] << "\n";
//===============================================================================================================
// Y11
ap::complex_2d_array SC_Y11_1, SC_Y11_2, SC_Y11_det_1, SC_Y11_det_2, SC_Y11_1_T, SC_Y11_2_T, SC_Y22_1, SC_Y22_2, SC_V_1, SC_V_2; 
ap::complex SC_DYE_1, SC_DYE_2;
SC_Y11_1.setlength(4,4);
SC_Y11_2.setlength(4,4);
SC_Y11_det_1.setlength(4,4);
SC_Y11_det_2.setlength(4,4);
SC_Y11_1_T.setlength(4,4);
SC_Y11_2_T.setlength(4,4);
SC_Y22_1.setlength(4,3);
SC_Y22_2.setlength(4,3);
SC_V_1.setlength(4,1);
SC_V_2.setlength(4,1);
for (int i = 3; i < 7; i++)    
{
   for (int j = 3; j < 7; j++)  
   {
	   SC_Y11_1(i-3,j-3) = SC_YE_1(i,j);
	   SC_Y11_2(i-3,j-3) = SC_YE_2(i,j);
	   SC_Y11_det_1(i-3,j-3) = SC_YE_1(i,j);
	   SC_Y11_det_2(i-3,j-3) = SC_YE_2(i,j);
   }
}
for (int i = 3; i < 7; i++)    
{
   for (int j = 0; j < 3; j++)  
   {
	   SC_Y22_1(i-3,j) = SC_YE_1(i,j);
	   SC_Y22_2(i-3,j) = SC_YE_2(i,j);
   }
}


//cout << "Y11 [0][0].X = " << SC_Y11_1(0,0).x << "\n"; 
//cout << "Y11 [0][0].Y = " << SC_Y11_2(0,0).y << "\n"; 

SC_DYE_1 = cmatrixdet(SC_Y11_det_1, 4);
SC_DYE_2 = cmatrixdet(SC_Y11_det_2, 4);

//cout << "DYE_1.x = " << SC_DYE_1.x << "\n";
//cout << "DYE_1.y = " << SC_DYE_1.y << "\n";
//cout << "DYE_2.x = " << SC_DYE_2.x << "\n";
//cout << "DYE_2.y = " << SC_DYE_2.y << "\n";
/*
ap::complex a;
ap::complex_2d_array b = ;
a = cmatrixdet(b, 2);
*/

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_V_1(i,j) += SC_Y22_1(i,k) * SC_UBE_1(k,j);
		 SC_V_2(i,j) += SC_Y22_2(i,k) * SC_UBE_2(k,j);
	  }
   }
}

for (int i = 0; i < 4; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "V.x [" << i << "] = " << SC_V_1(i,j).x << "\n";
	cout << "V.y [" << i << "] = " << SC_V_1(i,j).y << "\n";
	}
}

//===============================================================================================================
// US
ap::complex_2d_array SC_US_1, SC_US_2, SC_US1_1, SC_US1_2, SC_IK_1, SC_IK_2, SC_Z2_1, SC_Z2_2;
ap::complex_2d_array SC_Y11_inv_1, SC_Y11_inv_2;
SC_Y11_inv_1.setlength(4,4);
SC_Y11_inv_2.setlength(4,4);
SC_US_1.setlength(4,1);
SC_US_2.setlength(4,1);
SC_US1_1.setlength(4,1);
SC_US1_2.setlength(4,1);
SC_IK_1.setlength(3,1);
SC_IK_2.setlength(3,1);
SC_Z2_1.setlength(3,3);
SC_Z2_2.setlength(3,3);

for (int i=0; i < 4; i++)
	{
	for (int j=0; j < 4; j++)
		{
		SC_Y11_inv_1(i,j) = SC_Y11_1(i,j);
		SC_Y11_inv_2(i,j) = SC_Y11_2(i,j);
		}
	}


// Inverse Z2
cmatrixinverse(SC_Y11_inv_1,4);
cmatrixinverse(SC_Y11_inv_2,4);

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 4; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_US_1(i,j) += SC_Y11_inv_1(i,k) * - SC_V_1(k,j);
		 SC_US_2(i,j) += SC_Y11_inv_2(i,k) * - SC_V_2(k,j);
	  }
   }
}

for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1; j++)
		{
		SC_US1_1(i,j) = SC_US_1(i,j);
		SC_US1_2(i,j) = SC_US_2(i,j);
		}
	}

for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_Z2_1(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	SC_Z2_1(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	SC_Z2_2(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	SC_Z2_2(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
   }
}


// Inverse Z2
//cmatrixinverse(SC_Z2_1,3);
//cmatrixinverse(SC_Z2_2,3);


ap::complex_2d_array SC_Y11E_1, SC_Y11E_2, SC_Y11E2_1, SC_Y11E2_2, SC_Y11E_INV_1, SC_Y11E_INV_2, SC_V1_1, SC_V1_2, SC_UUU_1, SC_UUU_2;
SC_Y11E_1.setlength(3,3);
SC_Y11E_2.setlength(3,3);
SC_Y11E2_1.setlength(3,3);
SC_Y11E2_2.setlength(3,3);
SC_Y11E_INV_1.setlength(3,3);
SC_Y11E_INV_2.setlength(3,3);
SC_V1_1.setlength(3,1);
SC_V1_2.setlength(3,1);
SC_UUU_1.setlength(3,1);
SC_UUU_2.setlength(3,1);



for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 3; j++)  
   {
	SC_Y11E_1(i,j) = SC_Y11_1(i,j);
	SC_Y11E_2(i,j) = SC_Y11_2(i,j);
	SC_Y11E2_1(i,j) = SC_Y11_1(i,j);
	SC_Y11E2_2(i,j) = SC_Y11_2(i,j);
	}
}


for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_V1_1(i,j) = SC_V_1(i,j);
	SC_V1_2(i,j) = SC_V_2(i,j);
	}
}

// Y11E
cmatrixinverse(SC_Y11E2_1, 3);
cmatrixinverse(SC_Y11E2_2, 3);


// CUT
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	SC_Y11E_INV_1(i,j).x = min(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x);
	SC_Y11E_INV_1(i,j).y = min(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y);
	SC_Y11E_INV_2(i,j).x = max(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x);
	SC_Y11E_INV_2(i,j).y = max(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y);
	}
}
// END OF CUT



// -1 * Y11E
cinterval iSC_Y11E_INV[3][3];


// CUT
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	Inf(Re(iSC_Y11E_INV[i][j])) = _real(min(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x));
	Inf(Im(iSC_Y11E_INV[i][j])) = _real(min(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y));
	Sup(Re(iSC_Y11E_INV[i][j])) = _real(max(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x));
	Sup(Im(iSC_Y11E_INV[i][j])) = _real(max(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y));
	}
}
// END OF CUT


for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	iSC_Y11E_INV[i][j] = -1.0 * iSC_Y11E_INV[i][j];
	}
}


for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	SC_Y11E_1(i,j).x = _double(Inf(Re(iSC_Y11E_INV[i][j])));
	SC_Y11E_1(i,j).y = _double(Inf(Im(iSC_Y11E_INV[i][j])));
	SC_Y11E_2(i,j).x = _double(Sup(Re(iSC_Y11E_INV[i][j])));
	SC_Y11E_2(i,j).y = _double(Sup(Im(iSC_Y11E_INV[i][j])));
	}
}

// UUU = - Y11E ^ (-1) * V1
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_UUU_1(i,j) += SC_Y11E_1(i,k) * SC_V1_1(k,j);
		 SC_UUU_2(i,j) += SC_Y11E_2(i,k) * SC_V1_2(k,j);
	  }
   }
}


// IK = Z2 ^ (-1) * UUU
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_IK_1(i,j) += apSC_Z2_inf(i,k) * SC_UUU_1(k,j);
		 SC_IK_2(i,j) += apSC_Z2_sup(i,k) * SC_UUU_2(k,j);
	  }
   }
}
 
//cout << "US1.x [0,0] = " << SC_US1_1(0,0).x << "\n";
//cout << "US1.y [0,0] = " << SC_US1_1(0,0).y << "\n";
//cout << "IK.x [0,0] = " << SC_IK_1(0,0).x << "\n";
//cout << "IK.y [0,0] = " << SC_IK_1(0,0).y << "\n";
//===============================================================================================================
// End
cinterval iSC_IK[3][1];

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iSC_IK[i][j])) = real(SC_IK_1(i,j).x);
	Inf(Im(iSC_IK[i][j])) = real(SC_IK_1(i,j).y);
	Sup(Re(iSC_IK[i][j])) = real(SC_IK_2(i,j).x);
	Sup(Im(iSC_IK[i][j])) = real(SC_IK_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "IK [" << i << "] = " << iSC_IK[i][j] << "\n";
	}
}

/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2 [" << i << "," << j << "] = " << iSC_Z2[i][j] << "\n";
	}
}
*/
/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2.x [" << i << "," << j << "] = " << SC_Z2_1(i,j).x << "\n";
	cout << "Z2.y [" << i << "," << j << "] = " << SC_Z2_1(i,j).y << "\n";
	}
}
*/
complex cSC_a, cSC_aa, cSC_2_3_of_Pi;
Re(cSC_2_3_of_Pi) = 0.0;
Im(cSC_2_3_of_Pi) = double((2.0 / 3.0) * pi);
cSC_a = exp(cSC_2_3_of_Pi);
cSC_aa = pow(cSC_a, 2.0);
ap::complex_2d_array ISC_S, ISC_I120_1, ISC_I120_2;
ISC_S.setlength(3,3);
ISC_I120_1.setlength(3,1);
ISC_I120_2.setlength(3,1);
ISC_S(0,0).x = 1.0;
ISC_S(0,0).y = 0.0;
ISC_S(0,1).x = 1.0;
ISC_S(0,1).y = 0.0;
ISC_S(0,2).x = 1.0;
ISC_S(0,2).y = 0.0;
ISC_S(1,0).x = _double(Re(cSC_aa));
ISC_S(1,0).y = _double(Im(cSC_aa));
ISC_S(1,1).x = _double(Re(cSC_a));
ISC_S(1,1).y = _double(Im(cSC_a));
ISC_S(1,2).x = 1.0;
ISC_S(1,2).y = 0.0;
ISC_S(2,0).x = _double(Re(cSC_a));
ISC_S(2,0).y = _double(Im(cSC_a));
ISC_S(2,1).x = _double(Re(cSC_aa));
ISC_S(2,1).y = _double(Im(cSC_aa));
ISC_S(2,2).x = 1.0;
ISC_S(2,2).y = 0.0;
cmatrixinverse(ISC_S, 3);
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 ISC_I120_1(i,j) += ISC_S(i,k) * SC_IK_1(k,j);
		 ISC_I120_2(i,j) += ISC_S(i,k) * SC_IK_2(k,j);
	  }
   }
}
cinterval iISC_I120[3][1];
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_I120[i][j])) = real(ISC_I120_1(i,j).x);
	Inf(Im(iISC_I120[i][j])) = real(ISC_I120_1(i,j).y);
	Sup(Re(iISC_I120[i][j])) = real(ISC_I120_2(i,j).x);
	Sup(Im(iISC_I120[i][j])) = real(ISC_I120_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "I120 [" << i << "] = " << iISC_I120[i][j] << "\n";
	}
}
cout << "End of code" << "\n";

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_shortcircuit_2f_fcoor_ground()
{

// Расчет токов короткого замыкания. Двухфазное КЗ. Расчет в фазных координатах.
SC_UL = 110.0;
SC_UL_1 = 110.0;
SC_UL_2 = 110.0;
SC_L = 100.0;
SC_L_1 = 95.0;
SC_L_2 = 105.0;
SC_X0 = 0.4;
SC_X0_1 = 0.32; //0.32;
SC_X0_2 = 0.48; //0.48;
interval iSC_XL, iSC_XM, iSC_35;
cinterval iSC_ZL, iSC_ZM, iSC_C1;
cinterval iSC_Z1[3][3], iSC_Z2[3][3], iSC_Y1[3][3], iSC_Y2[3][3];
cinterval iSC_YC1[6][6], iSC_YC2[6][6], iSC_O[6][6], iSC_YB[12][12], iSC_YS1[9][9], iSC_M[9][12], iSC_M_trans[12][9];
interval iSC_P[3][6], iSC_O1[3][6];

ap::complex_2d_array apSC_Z1_inf;
ap::complex_2d_array apSC_Z1_sup;
apSC_Z1_inf.setlength(3,3);
apSC_Z1_sup.setlength(3,3);

ap::complex_2d_array apSC_Z2_inf;
ap::complex_2d_array apSC_Z2_sup;
apSC_Z2_inf.setlength(3,3);
apSC_Z2_sup.setlength(3,3);
//===============================================================================================================
Inf(iSC_UL) = _real(SC_UL_1);
Sup(iSC_UL) = _real(SC_UL_2);
Inf(iSC_L) = _real(SC_L_1);
Sup(iSC_L) = _real(SC_L_2);
cout << "L = " << iSC_L << "км \n";
Inf(iSC_X0) = _real(SC_X0_1);
Sup(iSC_X0) = _real(SC_X0_2);
Inf(iSC_35) = _real(3.5);
Sup(iSC_35) = _real(3.5);
iSC_UF = iSC_UL / sqrt(3.0);
iSC_X1 = iSC_X0 * iSC_L;
iSC_X0 = iSC_35 * iSC_X1;
iSC_XL = (iSC_X1 / 1.5) + (iSC_X0 / 3.0);
iSC_XM = iSC_XL - iSC_X1;
Re(iSC_ZL) = 0.0;
Im(iSC_ZL) = iSC_XL;
Re(iSC_ZM) = 0.0;
Im(iSC_ZM) = iSC_XM;
// Проверка
iSC_C1 = iSC_XL - iSC_XM;
//===============================================================================================================
iSC_Z1[0][0] = iSC_ZL;
iSC_Z1[0][1] = iSC_ZM;
iSC_Z1[0][2] = iSC_ZM;
iSC_Z1[1][0] = iSC_ZM;
iSC_Z1[1][1] = iSC_ZL;
iSC_Z1[1][2] = iSC_ZM;
iSC_Z1[2][0] = iSC_ZM;
iSC_Z1[2][1] = iSC_ZM;
iSC_Z1[2][2] = iSC_ZL;

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z1_inf(i,j).x = _double(Inf(Re(iSC_Z1[i][j])));
	apSC_Z1_inf(i,j).y = _double(Inf(Im(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).x = _double(Sup(Re(iSC_Z1[i][j])));
	apSC_Z1_sup(i,j).y = _double(Sup(Im(iSC_Z1[i][j])));
	}
}

cmatrixinverse(apSC_Z1_inf,3);
cmatrixinverse(apSC_Z1_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).x;
	Inf(Im(iSC_Y1[i][j])) = apSC_Z1_inf(i,j).y;
	Sup(Re(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).x;
	Sup(Im(iSC_Y1[i][j])) = apSC_Z1_sup(i,j).y;
	}
}
//===============================================================================================================
/*// Двухфазное КЗ (между фазами)
Re(iSC_Z2[0][0]) = 0.0;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 10000.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 10000.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 10000.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 0.0;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.001;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 10000.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.001;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 0.0;
Im(iSC_Z2[2][2]) = 0.0;
*/

// Однофазное КЗ на землю
Re(iSC_Z2[0][0]) = 10000000.0;
Im(iSC_Z2[0][0]) = 0.0;
Re(iSC_Z2[0][1]) = 0.0;
Im(iSC_Z2[0][1]) = 0.0;
Re(iSC_Z2[0][2]) = 0.0;
Im(iSC_Z2[0][2]) = 0.0;
Re(iSC_Z2[1][0]) = 0.0;
Im(iSC_Z2[1][0]) = 0.0;
Re(iSC_Z2[1][1]) = 0.0000001;
Im(iSC_Z2[1][1]) = 0.0;
Re(iSC_Z2[1][2]) = 0.0;
Im(iSC_Z2[1][2]) = 0.0;
Re(iSC_Z2[2][0]) = 0.0;
Im(iSC_Z2[2][0]) = 0.0;
Re(iSC_Z2[2][1]) = 0.0;
Im(iSC_Z2[2][1]) = 0.0;
Re(iSC_Z2[2][2]) = 0.0000001;
Im(iSC_Z2[2][2]) = 0.0;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	apSC_Z2_inf(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	apSC_Z2_inf(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	apSC_Z2_sup(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
	}
}

cmatrixinverse(apSC_Z2_inf,3);
cmatrixinverse(apSC_Z2_sup,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).x;
	Inf(Im(iSC_Y2[i][j])) = apSC_Z2_inf(i,j).y;
	Sup(Re(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).x;
	Sup(Im(iSC_Y2[i][j])) = apSC_Z2_sup(i,j).y;
	}
}

//===============================================================================================================
// Матрица YC1
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i][j]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i][j]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i][j]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i][j-3]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i][j-3]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i][j-3]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = - Inf(Re(iSC_Y1[i-3][j]));
	Inf(Im(iSC_YC1[i][j])) = - Inf(Im(iSC_Y1[i-3][j]));
	Sup(Re(iSC_YC1[i][j])) = - Sup(Re(iSC_Y1[i-3][j]));
	Sup(Im(iSC_YC1[i][j])) = - Sup(Im(iSC_Y1[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC1[i][j])) = Inf(Re(iSC_Y1[i-3][j-3]));
	Inf(Im(iSC_YC1[i][j])) = Inf(Im(iSC_Y1[i-3][j-3]));
	Sup(Re(iSC_YC1[i][j])) = Sup(Re(iSC_Y1[i-3][j-3]));
	Sup(Im(iSC_YC1[i][j])) = Sup(Im(iSC_Y1[i-3][j-3]));
	}
}
//===============================================================================================================
// Матрица YC2
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i][j]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i][j]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i][j]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i][j]));
	}
}
for (int i=0;i<3;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i][j-3]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i][j-3]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i][j-3]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i][j-3]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = - Inf(Re(iSC_Y2[i-3][j]));
	Inf(Im(iSC_YC2[i][j])) = - Inf(Im(iSC_Y2[i-3][j]));
	Sup(Re(iSC_YC2[i][j])) = - Sup(Re(iSC_Y2[i-3][j]));
	Sup(Im(iSC_YC2[i][j])) = - Sup(Im(iSC_Y2[i-3][j]));
	}
}
for (int i=3;i<6;i++)
{
	for (int j=3;j<6;j++)
	{
	Inf(Re(iSC_YC2[i][j])) = Inf(Re(iSC_Y2[i-3][j-3]));
	Inf(Im(iSC_YC2[i][j])) = Inf(Im(iSC_Y2[i-3][j-3]));
	Sup(Re(iSC_YC2[i][j])) = Sup(Re(iSC_Y2[i-3][j-3]));
	Sup(Im(iSC_YC2[i][j])) = Sup(Im(iSC_Y2[i-3][j-3]));
	}
}

//===============================================================================================================
// Матрица P
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_P[i][j] = 0.0;
	}
}
iSC_P[0][0] = 1.0;
iSC_P[1][1] = 1.0;
iSC_P[2][2] = 1.0;
//===============================================================================================================
// Матрица O1
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O1[i][j] = 0.0;
	}
}
//===============================================================================================================
// Матрица M
for (int i=0;i<9;i++)
{
	for (int j=0;j<12;j++)
	{
	Re(iSC_M[i][j]) = 0.0;
	Im(iSC_M[i][j]) = 0.0;
	}
}
for (int i=0;i<3;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i][j];
	}
}
for (int i=0;i<3;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i][j-6];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-3][j];
	}
}
for (int i=3;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = iSC_P[i-3][j-6];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(iSC_M[i][j]) = iSC_O1[i-6][j];
	}
}
for (int i=6;i<9;i++)
{
	for (int j=6;j<12;j++)
	{
	Re(iSC_M[i][j]) = - iSC_P[i-6][j-6];
	}
}

//===============================================================================================================
//Матрица O

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_O[i][j] = 0.0;
	}
}
//===============================================================================================================
//Матрица YB
for (int i=0;i<12;i++)
{
	for (int j=0;j<12;j++)
	{
	iSC_YB[i][j] = 0.0;
	}
}
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_YC1[i][j];
	}
}
for (int i=0;i<6;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_O[i][j-6];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=0;j<6;j++)
	{
	iSC_YB[i][j] = iSC_O[i-6][j];
	}
}
for (int i=6;i<12;i++)
{
	for (int j=6;j<12;j++)
	{
	iSC_YB[i][j] = iSC_YC2[i-6][j-6];
	}
}
//===============================================================================================================
// M * Yb = (9x12) * (12x12)
cinterval iSC_Ytemp[9][12];
ap::complex_2d_array Y_1_shm, Y_2_shm, Y_1_shm_temp, Y_2_shm_temp, Y_1_shm_temp_reint, Y_2_shm_temp_reint;
Y_1_shm.setlength(9,9);
Y_2_shm.setlength(9,9);
Y_1_shm_temp.setlength(9,12);
Y_2_shm_temp.setlength(9,12);
Y_1_shm_temp_reint.setlength(9,12);
Y_2_shm_temp_reint.setlength(9,12);
ap::complex_2d_array Yb_1_shm, Yb_2_shm;
Yb_1_shm.setlength(12,12);
Yb_2_shm.setlength(12,12);
ap::complex M_element;
ap::complex_2d_array M_shm_transposed;
ap::complex_2d_array M_shm;
M_shm_transposed.setlength(12,9);
M_shm.setlength(9,12);

for (int i = 0; i < 9; i++)     
{
   for (int j = 0; j < 12; j++) 
   {
	M_shm(i,j).x = _double(Inf(Re(iSC_M[i][j])));
	M_shm(i,j).y = _double(0.0);
   }
}
for (int i = 0; i < 12; i++)      
{
   for (int j = 0; j < 12; j++)  
   {
	Yb_1_shm(i,j).x = _double(Inf(Re(iSC_YB[i][j])));
	Yb_1_shm(i,j).y = _double(Inf(Im(iSC_YB[i][j])));
	Yb_2_shm(i,j).x = _double(Sup(Re(iSC_YB[i][j])));
	Yb_2_shm(i,j).y = _double(Sup(Im(iSC_YB[i][j])));
   }
}

for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 12; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {
		Y_1_shm_temp(i,j) += M_shm(i,k) * Yb_1_shm(k,j);
		Y_2_shm_temp(i,j) += M_shm(i,k) * Yb_2_shm(k,j);
	  }
   }
}


for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp_reint(i,j).x = min(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_1_shm_temp_reint(i,j).y = min(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	Y_2_shm_temp_reint(i,j).x = max(Y_1_shm_temp(i,j).x,Y_2_shm_temp(i,j).x);
	Y_2_shm_temp_reint(i,j).y = max(Y_1_shm_temp(i,j).y,Y_2_shm_temp(i,j).y);
	}
}
for (int i = 0; i < 9; i++)
{
	for (int j = 0; j < 12; j++)
	{
	Y_1_shm_temp(i,j) = Y_1_shm_temp_reint(i,j);
    Y_2_shm_temp(i,j) = Y_2_shm_temp_reint(i,j);
    }
}


// Transpose M
for (int i=0; i < 9; i++)
	{
	for (int j=0; j < 12 ; j++)
		{
			M_element = M_shm(i,j);
			M_shm_transposed(j,i) = M_element;
		}
	}



// Y_temp * M^transpose = (9x12) * (12x9)
for (int i = 0; i < 9; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 9; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 12; k++)  // Кол-во столбцов 1 матрицы
	  {

		 Y_1_shm(i,j) += Y_1_shm_temp(i,k) * M_shm_transposed(k,j);
		 Y_2_shm(i,j) += Y_2_shm_temp(i,k) * M_shm_transposed(k,j);
	  }
   }
}
ap::complex_2d_array SC_Y1_1, SC_Y1_2, SC_Y2_1, SC_Y2_2, SC_Y12_1, SC_Y12_2, SC_Y21_1, SC_Y21_2;
SC_Y1_1.setlength(6,6);
SC_Y1_2.setlength(6,6);
SC_Y2_1.setlength(3,3);
SC_Y2_2.setlength(3,3);
SC_Y12_1.setlength(6,3);
SC_Y12_2.setlength(6,3);
SC_Y21_1.setlength(3,6);
SC_Y21_2.setlength(3,6);
for (int i = 0; i < 6; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y1_1(i,j) = Y_1_shm(i,j);
    SC_Y1_2(i,j) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y2_1(i-6,j-6) = Y_1_shm(i,j);
    SC_Y2_2(i-6,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 0; i < 6; i++)
{
	for (int j = 6; j < 9; j++)
	{
	SC_Y12_1(i,j-6) = Y_1_shm(i,j);
    SC_Y12_2(i,j-6) = Y_2_shm(i,j);
    }
}
for (int i = 6; i < 9; i++)
{
	for (int j = 0; j < 6; j++)
	{
	SC_Y21_1(i-6,j) = Y_1_shm(i,j);
    SC_Y21_2(i-6,j) = Y_2_shm(i,j);
    }
}
ap::complex_2d_array SC_e1, SC_e1_T, SC_Y12_1_T, SC_Y12_2_T, SC_C1_1, SC_C1_2, SC_C2_1, SC_C2_2, SC_C3_1, SC_C3_2;
ap::complex_2d_array SC_S1_1, SC_S1_2, SC_S2_1, SC_S2_2, SC_YE_1, SC_YE_2;
SC_S1_1.setlength(7,1);
SC_S1_2.setlength(7,1);
SC_S2_1.setlength(7,6);
SC_S2_2.setlength(7,6);
SC_YE_1.setlength(7,7);
SC_YE_2.setlength(7,7);
SC_e1.setlength(3,1);
SC_e1_T.setlength(1,3);
SC_C1_1.setlength(6,1);
SC_C1_2.setlength(6,1);
SC_C2_1.setlength(1,6);
SC_C2_2.setlength(1,6);
SC_C3_1.setlength(1,1);
SC_C3_2.setlength(1,1);
SC_Y12_1_T.setlength(3,6);
SC_Y12_2_T.setlength(3,6);
for (int i = 0; i < 3; i++)
{
	for (int j = 0; j < 1; j++)
	{
	SC_e1(i,j).x = 1.0;
    SC_e1(i,j).y = 0.0;
    }
}

// C1 = Y12 * e1
for (int i = 0; i < 6; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C1_1(i,j) += SC_Y12_1(i,k) * SC_e1(k,j);
		 SC_C1_2(i,j) += SC_Y12_2(i,k) * SC_e1(k,j);
	  }
   }
}

// Transpose e1
for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1 ; j++)
		{
		M_element = SC_e1(i,j);
		SC_e1_T(j,i) = M_element;
		}
	}

// Transpose Y12
ap::complex ELEM_1, ELEM_2;
for (int i=0; i < 6; i++)
	{
	for (int j=0; j < 3 ; j++)
		{
		ELEM_1 = SC_Y12_1(i,j);
		ELEM_2 = SC_Y12_2(i,j);
		SC_Y12_1_T(j,i) = ELEM_1;
		SC_Y12_2_T(j,i) = ELEM_2;
		}
	}

// C2 = e1^T * Y12^T
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 6; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C2_1(i,j) += SC_e1_T(i,k) * SC_Y12_1_T(k,j);
		 SC_C2_2(i,j) += SC_e1_T(i,k) * SC_Y12_2_T(k,j);
	  }
   }
}
// C3 = e1^T * Y2 * e1
ap::complex_2d_array SC_TEMP_1, SC_TEMP_2;
SC_TEMP_1.setlength(1,3);
SC_TEMP_2.setlength(1,3);
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 3; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_TEMP_1(i,j) += SC_e1_T(i,k) * SC_Y2_1(k,j);
		 SC_TEMP_2(i,j) += SC_e1_T(i,k) * SC_Y2_2(k,j);
	  }
   }
}
for (int i = 0; i < 1; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_C3_1(i,j) += SC_TEMP_1(i,k) * SC_e1(k,j);
		 SC_C3_2(i,j) += SC_TEMP_2(i,k) * SC_e1(k,j);
	  }
   }
}
//===============================================================================================================
// S1
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C1_1(i,j);
	SC_S1_2(i,j) = SC_C1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_S1_1(i,j) = SC_C3_1(i-6,j);
	SC_S1_2(i,j) = SC_C3_2(i-6,j);
   }
}
//===============================================================================================================
// S2 = STACK(Y1, C2) <=> 6x6 & 1x6
for (int i = 0; i < 6; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_Y1_1(i,j);
	SC_S2_2(i,j) = SC_Y1_2(i,j);
   }
}
for (int i = 6; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_S2_1(i,j) = SC_C2_1(i-6,j);
	SC_S2_2(i,j) = SC_C2_2(i-6,j);
   }
}
//===============================================================================================================
// YE = augment(S2, S1) = 7x6 & 7x1
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 6; j++)  
   {
	SC_YE_1(i,j) = SC_S2_1(i,j);
	SC_YE_2(i,j) = SC_S2_2(i,j);
   }
}
for (int i = 0; i < 7; i++)    
{
   for (int j = 6; j < 7; j++)  
   {
	SC_YE_1(i,j) = SC_S1_1(i,j-6);
	SC_YE_2(i,j) = SC_S1_2(i,j-6);
   }
}
/*
for (int i = 0; i < 7; i++)    
{
   for (int j = 0; j < 7; j++)  
   {
	cout << "YE.x ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).x << "\n";
	cout << "YE.y ["<< i << ", "<< j << "] = " << SC_YE_1(i,j).y << "\n";
   }
}
*/
//===============================================================================================================
// UBE
cinterval iSC_UBE[3][1];
ap::complex_2d_array SC_UBE_1, SC_UBE_2;
SC_UBE_1.setlength(3,1);
SC_UBE_2.setlength(3,1);
cinterval SC_2_3_of_Pi;
Re(SC_2_3_of_Pi) = 0.0;
Im(SC_2_3_of_Pi) = double((2.0 / 3.0) * pi);

if (param_unsim == 0)
{
// then 0
iSC_UBE[0][0] = 1.0;
iSC_UBE[1][0] = exp(- SC_2_3_of_Pi);
iSC_UBE[2][0] = exp(SC_2_3_of_Pi);
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	   iSC_UBE[i][j] = iSC_UBE[i][j] * iSC_UF;
   }
}
for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_UBE_1(i,j).x = _double(Inf(Re(iSC_UBE[i][j])));
	SC_UBE_1(i,j).y = _double(Inf(Im(iSC_UBE[i][j])));
	SC_UBE_2(i,j).x = _double(Sup(Re(iSC_UBE[i][j])));
	SC_UBE_2(i,j).y = _double(Sup(Im(iSC_UBE[i][j])));
   }
}
//end 0
}
else if (param_unsim == 1)
{
// then 1
get_unsimmetric_fcoor();
	for (int i = 0; i < 3; i++)    
	{
		for (int j = 0; j < 1; j++)  
		{
		SC_UBE_1(i,j) = ISC_U_unsim(i,j);
		SC_UBE_2(i,j) = ISC_U_unsim(i,j);
		}
}

}

//cout << "y = " << iSC_UBE[1][0] << "\n";
//===============================================================================================================
// Y11
ap::complex_2d_array SC_Y11_1, SC_Y11_2, SC_Y11_det_1, SC_Y11_det_2, SC_Y11_1_T, SC_Y11_2_T, SC_Y22_1, SC_Y22_2, SC_V_1, SC_V_2; 
ap::complex SC_DYE_1, SC_DYE_2;
SC_Y11_1.setlength(4,4);
SC_Y11_2.setlength(4,4);
SC_Y11_det_1.setlength(4,4);
SC_Y11_det_2.setlength(4,4);
SC_Y11_1_T.setlength(4,4);
SC_Y11_2_T.setlength(4,4);
SC_Y22_1.setlength(4,3);
SC_Y22_2.setlength(4,3);
SC_V_1.setlength(4,1);
SC_V_2.setlength(4,1);
for (int i = 3; i < 7; i++)    
{
   for (int j = 3; j < 7; j++)  
   {
	   SC_Y11_1(i-3,j-3) = SC_YE_1(i,j);
	   SC_Y11_2(i-3,j-3) = SC_YE_2(i,j);
	   SC_Y11_det_1(i-3,j-3) = SC_YE_1(i,j);
	   SC_Y11_det_2(i-3,j-3) = SC_YE_2(i,j);
   }
}
for (int i = 3; i < 7; i++)    
{
   for (int j = 0; j < 3; j++)  
   {
	   SC_Y22_1(i-3,j) = SC_YE_1(i,j);
	   SC_Y22_2(i-3,j) = SC_YE_2(i,j);
   }
}


//cout << "Y11 [0][0].X = " << SC_Y11_1(0,0).x << "\n"; 
//cout << "Y11 [0][0].Y = " << SC_Y11_2(0,0).y << "\n"; 

SC_DYE_1 = cmatrixdet(SC_Y11_det_1, 4);
SC_DYE_2 = cmatrixdet(SC_Y11_det_2, 4);

//cout << "DYE_1.x = " << SC_DYE_1.x << "\n";
//cout << "DYE_1.y = " << SC_DYE_1.y << "\n";
//cout << "DYE_2.x = " << SC_DYE_2.x << "\n";
//cout << "DYE_2.y = " << SC_DYE_2.y << "\n";
/*
ap::complex a;
ap::complex_2d_array b = ;
a = cmatrixdet(b, 2);
*/

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_V_1(i,j) += SC_Y22_1(i,k) * SC_UBE_1(k,j);
		 SC_V_2(i,j) += SC_Y22_2(i,k) * SC_UBE_2(k,j);
	  }
   }
}

for (int i = 0; i < 4; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "V.x [" << i << "] = " << SC_V_1(i,j).x << "\n";
	cout << "V.y [" << i << "] = " << SC_V_1(i,j).y << "\n";
	}
}

//===============================================================================================================
// US
ap::complex_2d_array SC_US_1, SC_US_2, SC_US1_1, SC_US1_2, SC_IK_1, SC_IK_2, SC_Z2_1, SC_Z2_2;
ap::complex_2d_array SC_Y11_inv_1, SC_Y11_inv_2;
SC_Y11_inv_1.setlength(4,4);
SC_Y11_inv_2.setlength(4,4);
SC_US_1.setlength(4,1);
SC_US_2.setlength(4,1);
SC_US1_1.setlength(4,1);
SC_US1_2.setlength(4,1);
SC_IK_1.setlength(3,1);
SC_IK_2.setlength(3,1);
SC_Z2_1.setlength(3,3);
SC_Z2_2.setlength(3,3);

for (int i=0; i < 4; i++)
	{
	for (int j=0; j < 4; j++)
		{
		SC_Y11_inv_1(i,j) = SC_Y11_1(i,j);
		SC_Y11_inv_2(i,j) = SC_Y11_2(i,j);
		}
	}


// Inverse Z2
cmatrixinverse(SC_Y11_inv_1,4);
cmatrixinverse(SC_Y11_inv_2,4);

for (int i = 0; i < 4; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 4; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_US_1(i,j) += SC_Y11_inv_1(i,k) * - SC_V_1(k,j);
		 SC_US_2(i,j) += SC_Y11_inv_2(i,k) * - SC_V_2(k,j);
	  }
   }
}

for (int i=0; i < 3; i++)
	{
	for (int j=0; j < 1; j++)
		{
		SC_US1_1(i,j) = SC_US_1(i,j);
		SC_US1_2(i,j) = SC_US_2(i,j);
		}
	}

for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_Z2_1(i,j).x = _double(Inf(Re(iSC_Z2[i][j])));
	SC_Z2_1(i,j).y = _double(Inf(Im(iSC_Z2[i][j])));
	SC_Z2_2(i,j).x = _double(Sup(Re(iSC_Z2[i][j])));
	SC_Z2_2(i,j).y = _double(Sup(Im(iSC_Z2[i][j])));
   }
}


// Inverse Z2
//cmatrixinverse(SC_Z2_1,3);
//cmatrixinverse(SC_Z2_2,3);


ap::complex_2d_array SC_Y11E_1, SC_Y11E_2, SC_Y11E2_1, SC_Y11E2_2, SC_Y11E_INV_1, SC_Y11E_INV_2, SC_V1_1, SC_V1_2, SC_UUU_1, SC_UUU_2;
SC_Y11E_1.setlength(3,3);
SC_Y11E_2.setlength(3,3);
SC_Y11E2_1.setlength(3,3);
SC_Y11E2_2.setlength(3,3);
SC_Y11E_INV_1.setlength(3,3);
SC_Y11E_INV_2.setlength(3,3);
SC_V1_1.setlength(3,1);
SC_V1_2.setlength(3,1);
SC_UUU_1.setlength(3,1);
SC_UUU_2.setlength(3,1);



for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 3; j++)  
   {
	SC_Y11E_1(i,j) = SC_Y11_1(i,j);
	SC_Y11E_2(i,j) = SC_Y11_2(i,j);
	SC_Y11E2_1(i,j) = SC_Y11_1(i,j);
	SC_Y11E2_2(i,j) = SC_Y11_2(i,j);
	}
}


for (int i = 0; i < 3; i++)    
{
   for (int j = 0; j < 1; j++)  
   {
	SC_V1_1(i,j) = SC_V_1(i,j);
	SC_V1_2(i,j) = SC_V_2(i,j);
	}
}

// Y11E
cmatrixinverse(SC_Y11E2_1, 3);
cmatrixinverse(SC_Y11E2_2, 3);

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	SC_Y11E_INV_1(i,j).x = min(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x);
	SC_Y11E_INV_1(i,j).y = min(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y);
	SC_Y11E_INV_2(i,j).x = max(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x);
	SC_Y11E_INV_2(i,j).y = max(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y);
	}
}

// -1 * Y11E
cinterval iSC_Y11E_INV[3][3];


for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	Inf(Re(iSC_Y11E_INV[i][j])) = _real(min(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x));
	Inf(Im(iSC_Y11E_INV[i][j])) = _real(min(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y));
	Sup(Re(iSC_Y11E_INV[i][j])) = _real(max(SC_Y11E2_1(i,j).x, SC_Y11E2_2(i,j).x));
	Sup(Im(iSC_Y11E_INV[i][j])) = _real(max(SC_Y11E2_1(i,j).y, SC_Y11E2_2(i,j).y));
	}
}


for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	iSC_Y11E_INV[i][j] = -1.0 * iSC_Y11E_INV[i][j];
	}
}


for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	SC_Y11E_1(i,j).x = _double(Inf(Re(iSC_Y11E_INV[i][j])));
	SC_Y11E_1(i,j).y = _double(Inf(Im(iSC_Y11E_INV[i][j])));
	SC_Y11E_2(i,j).x = _double(Sup(Re(iSC_Y11E_INV[i][j])));
	SC_Y11E_2(i,j).y = _double(Sup(Im(iSC_Y11E_INV[i][j])));
	}
}

// UUU = - Y11E ^ (-1) * V1
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_UUU_1(i,j) += SC_Y11E_1(i,k) * SC_V1_1(k,j);
		 SC_UUU_2(i,j) += SC_Y11E_2(i,k) * SC_V1_2(k,j);
	  }
   }
}


// IK = Z2 ^ (-1) * UUU
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 SC_IK_1(i,j) += apSC_Z2_inf(i,k) * SC_UUU_1(k,j);
		 SC_IK_2(i,j) += apSC_Z2_sup(i,k) * SC_UUU_2(k,j);
	  }
   }
}
 
//cout << "US1.x [0,0] = " << SC_US1_1(0,0).x << "\n";
//cout << "US1.y [0,0] = " << SC_US1_1(0,0).y << "\n";
//cout << "IK.x [0,0] = " << SC_IK_1(0,0).x << "\n";
//cout << "IK.y [0,0] = " << SC_IK_1(0,0).y << "\n";
//===============================================================================================================
// End
cinterval iSC_IK[3][1];

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iSC_IK[i][j])) = real(SC_IK_1(i,j).x);
	Inf(Im(iSC_IK[i][j])) = real(SC_IK_1(i,j).y);
	Sup(Re(iSC_IK[i][j])) = real(SC_IK_2(i,j).x);
	Sup(Im(iSC_IK[i][j])) = real(SC_IK_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "IK [" << i << "] = " << iSC_IK[i][j] << "\n";
	}
}

/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2 [" << i << "," << j << "] = " << iSC_Z2[i][j] << "\n";
	}
}
*/
/*
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 3; j++)  
	{
	cout << "Z2.x [" << i << "," << j << "] = " << SC_Z2_1(i,j).x << "\n";
	cout << "Z2.y [" << i << "," << j << "] = " << SC_Z2_1(i,j).y << "\n";
	}
}
*/
complex cSC_a, cSC_aa, cSC_2_3_of_Pi;
Re(cSC_2_3_of_Pi) = 0.0;
Im(cSC_2_3_of_Pi) = double((2.0 / 3.0) * pi);
cSC_a = exp(cSC_2_3_of_Pi);
cSC_aa = pow(cSC_a, 2.0);
ap::complex_2d_array ISC_S, ISC_I120_1, ISC_I120_2;
ISC_S.setlength(3,3);
ISC_I120_1.setlength(3,1);
ISC_I120_2.setlength(3,1);
ISC_S(0,0).x = 1.0;
ISC_S(0,0).y = 0.0;
ISC_S(0,1).x = 1.0;
ISC_S(0,1).y = 0.0;
ISC_S(0,2).x = 1.0;
ISC_S(0,2).y = 0.0;
ISC_S(1,0).x = _double(Re(cSC_aa));
ISC_S(1,0).y = _double(Im(cSC_aa));
ISC_S(1,1).x = _double(Re(cSC_a));
ISC_S(1,1).y = _double(Im(cSC_a));
ISC_S(1,2).x = 1.0;
ISC_S(1,2).y = 0.0;
ISC_S(2,0).x = _double(Re(cSC_a));
ISC_S(2,0).y = _double(Im(cSC_a));
ISC_S(2,1).x = _double(Re(cSC_aa));
ISC_S(2,1).y = _double(Im(cSC_aa));
ISC_S(2,2).x = 1.0;
ISC_S(2,2).y = 0.0;
cmatrixinverse(ISC_S, 3);
for (int i = 0; i < 3; i++)      // Кол-во строк 1 матрицы
{
   for (int j = 0; j < 1; j++)  // Кол-во столбцов 2 матрицы
   {
	  for (int k = 0; k < 3; k++)  // Кол-во столбцов 1 матрицы
	  {
		 ISC_I120_1(i,j) += ISC_S(i,k) * SC_IK_1(k,j);
		 ISC_I120_2(i,j) += ISC_S(i,k) * SC_IK_2(k,j);
	  }
   }
}
cinterval iISC_I120[3][1];
for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	Inf(Re(iISC_I120[i][j])) = real(ISC_I120_1(i,j).x);
	Inf(Im(iISC_I120[i][j])) = real(ISC_I120_1(i,j).y);
	Sup(Re(iISC_I120[i][j])) = real(ISC_I120_2(i,j).x);
	Sup(Im(iISC_I120[i][j])) = real(ISC_I120_2(i,j).y);
	}
}

for (int i = 0; i < 3; i++)    
{
	for (int j = 0; j < 1; j++)  
	{
	cout << "I120 [" << i << "] = " << iISC_I120[i][j] << "\n";
	}
}
cout << "End of code" << "\n";

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_ils()
{
// Interval linear system
cout << "Решение интервальной системы линейных уравнений\n";
cinterval ILS_Y, ILS_U, ILS_I, ILS_H, ILS_d, ILS_D;
/*
Inf(ILS_D) = 
Sup(ILS_D) = 

*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_bar()
{
// БРУСЫ
cout << "Перебор возможных решений ИСЛАУ\n";
//double BAR_Length[6], BAR_Yground[6], BAR_X1[11], BAR_X2[11], BAR_X3[11], BAR_Height[11], BAR_temperature[17];

double BAR_Length[6] = {	45.0, 
							47.0, 
							49.0, 
							51.0, 
							53.0, 
							55.0 };

double BAR_Yground[6] = {	0.05, 
							0.1, 
							0.15, 
							0.2, 
							0.25, 
							0.3 };

double BAR_X1[6] = {		-4.375,
							-4.625,
							-4.875,
							-5.125,
							-5.375,
							-5.625 };

double BAR_X2[6] = {		0.0,
							0.0,
							0.0,
							0.0,
							0.0,
							0.0 };

double BAR_X3[6] = {		4.375,
							4.625,
							4.875,
							5.125,
							5.375,
							5.625 };

double BAR_Height[6] = {	16.625,
							17.575,
							18.525,
							19.475,
							20.425,
							21.375 };

double BAR_temperature[9] = {
							-40.0,
							-30.0,
							-20.0,
							-10.0,
							0.0,
							10.0,
							20.0,
							30.0,
							40.0 };







int BAR_thread = 6*6*6*6*6*6*9;
int BAR_count = 0;

std::ofstream file;// can be merged to std::ofstream file("file.txt");
file.open("INTCALC_data.csv");

for (int i0 = 0; i0 < 6; i0++)	// Length
{
	for (int i1 = 0; i1 < 6; i1++)	// Yground
	{
		for (int i2 = 0; i2 < 6; i2++)	// X1    
		{
			for (int i3 = 0; i3 < 6; i3++)	// X2  
			{
				for (int i4 = 0; i4 < 6; i4++)	// x3    
				{
					for (int i5 = 0; i5 < 6; i5++)	// Height  
					{
						for (int i6 = 0; i6 < 9; i6++)    
						{
							Length = BAR_Length[i0];			// SET TO 6
							Yground = BAR_Yground[i1];			// SET TO 6
							X1 = BAR_X1[i2];					// SET TO 6
							X2 = BAR_X2[i3];					// SET TO 6
							X3 = BAR_X3[i4];					// SET TO 6
							Height = BAR_Height[i5];			// SET TO 6
							temperature = BAR_temperature[i6];	// SET TO 9
							getintvloss();
							BAR_UA[BAR_count] = RP_Umod_max[0];
							BAR_UB[BAR_count] = RP_Umod_max[1];
							BAR_UC[BAR_count] = RP_Umod_max[2];
							file << BAR_UA[BAR_count] << ";" << BAR_UB[BAR_count] << ";" << BAR_UC[BAR_count] << "\n";
							BAR_count = BAR_count++;
							cout << "THREAD: " << BAR_count << " / " << BAR_thread << "\n";
						}
					}
				}
			}
		}
	}
}



file.close();// is not necessary because the destructor closes the open file by default





}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void process_nonsinusoidal()
{

// http://lvvas.co.cc/index.php?option=com_content&view=article&id=135&Itemid=157
// initial parameters

double MagConst;

interval IP_MagConst, IP_WC_temperature, IP_Yground, IP_freq, IP_cfreq, IP_WC_R0, IP_WC_S, IP_WC_R,
		 IP_WC_R_Equal, IP_X1, IP_X2, IP_X3, IP_Y1, IP_Y2, IP_Y3;

interval RP_d12, RP_d13, RP_d23, RP_D12, RP_D13, RP_D23;

interval IP_U_A_l, IP_U_B_l, IP_U_C_l, IP_U_A_f, IP_U_B_f, IP_U_C_f;

cinterval IP_S;

Re(IP_S) = 10.0 * 1000000.0;
Im(IP_S) = 10.0 * 1000000.0;

MagConst = 0.00000125663706;

IP_U_A_l = IP_Ua_l * 1000.0;
IP_U_B_l = IP_Ub_l * 1000.0;
IP_U_C_l = IP_Uc_l * 1000.0;
//cout << "U" << IP_U_A_l;
IP_U_A_f = IP_U_A_l / sqrt(3.0);
IP_U_B_f = IP_U_B_l / sqrt(3.0);
IP_U_C_f = IP_U_C_l / sqrt(3.0);


if (param_interval == 0)
{

IP_MagConst = MagConst;	// магнитная постоянная
IP_Yground = Yground;//(0.05,0.3); // проводимость земли
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
IP_L = Length;
IP_WC_R0 = 0.1197*(1.0 + 0.004*(temperature-20.0));
IP_WC_S = 240.0;
IP_WC_R = 0.012;
IP_X1 = X1;
IP_X2 = X2;
IP_X3 = X3;
IP_Y1 = Height;
IP_Y2 = Height;
IP_Y3 = Height;
}

if (param_interval == 1)
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_WC_temperature) = -40.0;
Sup(IP_WC_temperature) = 40.0;
Inf(IP_Yground) = 0.05; // проводимость земли
Sup(IP_Yground) = 0.3;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.1197*(1.0 + 0.004*(IP_WC_temperature-20.0));
IP_WC_S = 240.0;
IP_WC_R = 0.012;
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}

if (param_interval == 2)	// fazanord min & max parameters
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_Yground) = 0.01; // проводимость земли
Sup(IP_Yground) = 0.01;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.122;
IP_WC_S = 240.0;
IP_WC_R = 0.012;
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}
//cout << IP_X1 << "\n";
//cout << IP_X2 << "\n";

//cout << "L = " << IP_L << " км" << "\n";
//cout << "L = " << IP_L << " км" << "\n";
//cout << IP_X2 << "\n";
interval IP_X12;
Inf(IP_X12) = Inf(IP_X1) - Inf(IP_X2);
Sup(IP_X12) = Sup(IP_X1) - Sup(IP_X2);
Inf(IP_X12) = Inf(IP_X12) * Inf(IP_X12);
Sup(IP_X12) = Sup(IP_X12) * Sup(IP_X12);

RP_d12 = sqrt(IP_X12 + Power(IP_Y1 - IP_Y2, 2.0));
RP_d13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 - IP_Y3), 2.0));
RP_d23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 - IP_Y3), 2.0));

RP_D12 = sqrt(IP_X12 + Power(IP_Y1 + IP_Y2, 2.0));
RP_D13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 + IP_Y3), 2.0));
RP_D23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 + IP_Y3), 2.0));
//RP_D12 = RP_D23; 
//cout << RP_D12 << "\n";
cinterval RP_Z11, RP_Z12, RP_Z13, RP_Z21, RP_Z22, RP_Z23, RP_Z31, RP_Z32, RP_Z33, RP_Z11out, RP_Z11inn;
interval pow_0_755, pow_0_83;
pow_0_755 = 0.755;
pow_0_83 = 0.83;

if (phase_fault == 0)
{

Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;
}

if (phase_fault == 1)
{
Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;

// -----------------------------------------------------------------------------------------
Re(RP_Z11) = 1000000.0;
Im(RP_Z11) = 0.0;
}

cinterval RP_d_p[3][3], RP_d_n[3][3];
RP_Z[0][0] = RP_Z11;
RP_Z[0][1] = RP_Z12;
RP_Z[0][2] = RP_Z13;
RP_Z[1][0] = RP_Z21;
RP_Z[1][1] = RP_Z22;
RP_Z[1][2] = RP_Z23;
RP_Z[2][0] = RP_Z31;
RP_Z[2][1] = RP_Z32;
RP_Z[2][2] = RP_Z33;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	RP_Z[i][j] = IP_L * RP_Z[i][j];

	}
}

//cout << RP_Z[0][0] << "\n";

ap::complex_2d_array Z_inf;
ap::complex_2d_array Z_sup;
ap::complex_2d_array D_inf;
ap::complex_2d_array D_sup;
Z_inf.setlength(3,3);
Z_sup.setlength(3,3);
D_inf.setlength(3,3);
D_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Z_inf(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	Z_inf(i,j).y = _double(Inf(Im(RP_Z[i][j])));
	Z_sup(i,j).x = _double(Sup(Re(RP_Z[i][j])));
	Z_sup(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	}
}
//cout << Z_inf(0,0).x << "\n";
cmatrixinverse(Z_inf,3);
cmatrixinverse(Z_sup,3);
//cout << Z_inf(0,0).x << "\n";
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_d_p[i][j])) = Z_inf(i,j).x;
	Inf(Im(RP_d_p[i][j])) = Z_inf(i,j).y;
	Sup(Re(RP_d_p[i][j])) = Z_sup(i,j).x;
	Sup(Im(RP_d_p[i][j])) = Z_sup(i,j).y;
	Inf(Re(RP_d_n[i][j])) = - Z_inf(i,j).x;// -Z
	Inf(Im(RP_d_n[i][j])) = - Z_inf(i,j).y;
	Sup(Re(RP_d_n[i][j])) = - Z_sup(i,j).x;
	Sup(Im(RP_d_n[i][j])) = - Z_sup(i,j).y;
	}
}

/*
//inv
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_Z[i][j]) = Inf(RP_Z[i][j]) / 100.0;
	Sup(RP_Z[i][j]) = Sup(RP_Z[i][j]) / 100.0;
	}
}

cinterval M_A[3][3], M_B[3][3], M_X[3][3], M_J[3][3];

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
		M_A[i][j] = RP_Z[i][j]; 
		Re(M_J[i][j]) = 0.0;
		Im(M_J[i][j]) = 0.0;
	}
}

Re(M_J[0][0]) = 1.0;
Im(M_J[0][0]) = 0.0;
Re(M_J[1][1]) = 1.0;
Im(M_J[1][1]) = 0.0;
Re(M_J[2][2]) = 1.0;
Im(M_J[2][2]) = 0.0;

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	M_B[i][j] = M_J[i][j] - M_A[i][j]; 
	}
}

cout << abs(M_B[0][0]) << "\n";


// норма комплексной матрицы

double M_B_norm_row, M_B_norm_col;

interval M_B_row_1, M_B_row_2, M_B_row_3;

//M_B_row_1;
for (int j=0; j<3;j++)
{
M_B_row_1 += abs(M_B[0][j]);
M_B_row_2 += abs(M_B[1][j]);
M_B_row_3 += abs(M_B[2][j]);
}

 
//cout << M_B_row_2 << "\n";

const int len = 3;
const int start = 0;
double arr[len];
arr[0] = _double(Sup(M_B_row_1));
arr[1] = _double(Sup(M_B_row_2));
arr[2] = _double(Sup(M_B_row_3));
double* first = arr;
double* last = arr + len;
double* asd;
asd = max_element(first,last);


cout << "Matrix norm" << *asd << "\n";

M_B_norm_col = 1.0 / (1.0 - *asd);


Inf(M_X[0][0]) = - M_B_norm_col;
Sup(M_X[0][0]) = M_B_norm_col + 2.0;
Inf(M_X[0][1]) = - M_B_norm_col;
Sup(M_X[0][1]) = M_B_norm_col;
Inf(M_X[0][2]) = - M_B_norm_col;
Sup(M_X[0][2]) = M_B_norm_col;
Inf(M_X[1][0]) = - M_B_norm_col;
Sup(M_X[1][0]) = M_B_norm_col;
Inf(M_X[1][1]) = - M_B_norm_col;
Sup(M_X[1][1]) = M_B_norm_col + 2.0;
Inf(M_X[1][2]) = - M_B_norm_col;
Sup(M_X[1][2]) = M_B_norm_col;
Inf(M_X[2][0]) = - M_B_norm_col;
Sup(M_X[2][0]) = M_B_norm_col;
Inf(M_X[2][1]) = - M_B_norm_col;
Sup(M_X[2][1]) = M_B_norm_col;
Inf(M_X[2][2]) = - M_B_norm_col;
Sup(M_X[2][2]) = M_B_norm_col + 2.0;
cout << "M_X" << M_X[0][0] << "\n";

cimatrix iM_X(3,3), iM_J(3,3), iM_A(3,3);

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	iM_X(i+1,j+1) = M_X[i][j];
	}
}
cout << "M_X" << iM_X(1,1) << "\n";


	for (int i=0; i<10; i++)
	{
	M_X[i][j] = mid(M_X) + M_X * (M_J - M_A * mid(M_X));
	}
	
//cout << Re(M_X(1,1)) << "\n";
*/
// Yrc

cinterval RP_Yrc[6][6];
RP_Yrc[0][0] = RP_d_p[0][0];
RP_Yrc[0][1] = RP_d_p[0][1];
RP_Yrc[0][2] = RP_d_p[0][2];
RP_Yrc[0][3] = RP_d_n[0][3];
RP_Yrc[0][4] = RP_d_n[0][1];
RP_Yrc[0][5] = RP_d_n[0][2];
RP_Yrc[1][0] = RP_d_p[1][0];
RP_Yrc[1][1] = RP_d_p[1][1];
RP_Yrc[1][2] = RP_d_p[1][2];
RP_Yrc[1][3] = RP_d_n[1][0];
RP_Yrc[1][4] = RP_d_n[1][1];
RP_Yrc[1][5] = RP_d_n[1][2];
RP_Yrc[2][0] = RP_d_p[2][0];
RP_Yrc[2][1] = RP_d_p[2][1];
RP_Yrc[2][2] = RP_d_p[2][2];
RP_Yrc[2][3] = RP_d_n[2][0];
RP_Yrc[2][4] = RP_d_n[2][1];
RP_Yrc[2][5] = RP_d_n[2][2];
RP_Yrc[3][0] = RP_d_n[0][0];
RP_Yrc[3][1] = RP_d_n[0][1];
RP_Yrc[3][2] = RP_d_n[0][2];
RP_Yrc[3][3] = RP_d_p[0][0];
RP_Yrc[3][4] = RP_d_p[0][1];
RP_Yrc[3][5] = RP_d_p[0][2];
RP_Yrc[4][0] = RP_d_n[1][0];
RP_Yrc[4][1] = RP_d_n[1][1];
RP_Yrc[4][2] = RP_d_n[1][2];
RP_Yrc[4][3] = RP_d_p[1][0];
RP_Yrc[4][4] = RP_d_p[1][1];			
RP_Yrc[4][5] = RP_d_p[1][2];
RP_Yrc[5][0] = RP_d_n[2][0];
RP_Yrc[5][1] = RP_d_n[2][1];
RP_Yrc[5][2] = RP_d_n[2][2];
RP_Yrc[5][3] = RP_d_p[2][0];
RP_Yrc[5][4] = RP_d_p[2][1];
RP_Yrc[5][5] = RP_d_p[2][2];

const int n_ABC=3;
interval RP_A[n_ABC][n_ABC], RP_B[n_ABC][n_ABC], RP_C[n_ABC][n_ABC];

RP_A[0][0] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y1) / (IP_WC_R * 100.0));
RP_A[0][1] = 1.8 * power(10.0, 7.0) * ln(RP_D12 / RP_d12);
RP_A[0][2] = 1.8 * power(10.0, 7.0) * ln(RP_D13 / RP_d13);
RP_A[1][0] = RP_A[0][1];
RP_A[1][1] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y2) / (IP_WC_R * 100.0));
RP_A[1][2] = 1.8 * power(10.0, 7.0) * ln(RP_D23 / RP_d23);
RP_A[2][0] = RP_A[0][2];
RP_A[2][1] = RP_A[1][2];
RP_A[2][2] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y3) / (IP_WC_R * 100.0));

ap::real_2d_array a_inf;
ap::real_2d_array a_sup;
a_inf.setlength(3,3);
a_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	a_inf(i,j) = _double(Inf(RP_A[i][j]));
	a_sup(i,j) = _double(Sup(RP_A[i][j]));
	}
}
//cout << a_inf(0,0) << "\n";
rmatrixinverse(a_inf,3);
rmatrixinverse(a_sup,3);
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_B[i][j]) = a_inf(i,j);
	Sup(RP_B[i][j]) = a_sup(i,j);
	}
}

/*
RP_C[0][0] = RP_B[0][0] + RP_B[0][1] + RP_B[0][2];
RP_C[0][1] = - RP_B[0][1];
RP_C[0][2] = - RP_B[0][2];
RP_C[1][0] = - RP_B[1][0];
RP_C[1][1] = RP_B[1][0] + RP_B[1][1] + RP_B[1][2];
RP_C[1][2] = - RP_B[1][2];
RP_C[2][0] = - RP_B[2][0];
RP_C[2][1] = - RP_B[2][1];
RP_C[2][2] = RP_B[2][0] + RP_B[2][1] + RP_B[2][2];
*/

//const int n_CY=6;
interval RP_Cy[6][6];
cinterval RP_Yc_inv[6][6]; 
interval i_zero(0.0);

RP_Cy[0][0] = RP_B[0][0];
RP_Cy[0][1] = RP_B[0][1];
RP_Cy[0][2] = RP_B[0][2];
RP_Cy[0][3] = i_zero;
RP_Cy[0][4] = i_zero;
RP_Cy[0][5] = i_zero;
RP_Cy[1][0] = RP_B[1][0];
RP_Cy[1][1] = RP_B[1][1];
RP_Cy[1][2] = RP_B[1][2];
RP_Cy[1][3] = i_zero;
RP_Cy[1][4] = i_zero;
RP_Cy[1][5] = i_zero;
RP_Cy[2][0] = RP_B[2][0];
RP_Cy[2][1] = RP_B[2][1];
RP_Cy[2][2] = RP_B[2][2];
RP_Cy[2][3] = i_zero;
RP_Cy[2][4] = i_zero; 
RP_Cy[2][5] = i_zero;
RP_Cy[3][0] = i_zero;
RP_Cy[3][1] = i_zero;
RP_Cy[3][2] = i_zero;
RP_Cy[3][3] = RP_B[0][0];
RP_Cy[3][4] = RP_B[0][1];
RP_Cy[3][5] = RP_B[0][2];
RP_Cy[4][0] = i_zero;
RP_Cy[4][1] = i_zero;
RP_Cy[4][2] = i_zero;
RP_Cy[4][3] = RP_B[1][0];
RP_Cy[4][4] = RP_B[1][1];
RP_Cy[4][5] = RP_B[1][2];
RP_Cy[5][0] = i_zero;
RP_Cy[5][1] = i_zero;
RP_Cy[5][2] = i_zero;
RP_Cy[5][3] = RP_B[2][0];
RP_Cy[5][4] = RP_B[2][1];
RP_Cy[5][5] = RP_B[2][2];

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Cy[i][j] = RP_Cy[i][j] * 0.5 * IP_L * IP_cfreq;
	//cout << "Cy" << RP_Cy[i][j] << "\n";
	}
}

cinterval RP_m_temp[6][6];
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(RP_m_temp[i][j]) = i_zero;
	Im(RP_m_temp[i][j]) = i_zero;
	Re(RP_Yc[i][j]) = i_zero;
	Im(RP_Yc[i][j]) = i_zero;
	}
}


for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Im(RP_m_temp[i][j]) = RP_Cy[i][j];
	}
}

double RP_Yc_I_Re, RP_Yc_I_Im, RP_Yc_S_Re, RP_Yc_S_Im; 
double RP_Yrc_I_Re, RP_Yrc_I_Im, RP_Yrc_S_Re, RP_Yrc_S_Im; 
double RP_m_temp_I_Re, RP_m_temp_I_Im, RP_m_temp_S_Re, RP_m_temp_S_Im;

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Yrc_I_Re = _double(Inf(Re(RP_Yrc[i][j])));
	RP_Yrc_I_Im = _double(Inf(Im(RP_Yrc[i][j])));
	RP_Yrc_S_Re = _double(Sup(Re(RP_Yrc[i][j])));
	RP_Yrc_S_Im = _double(Sup(Im(RP_Yrc[i][j])));
	
	RP_m_temp_I_Re = _double(Inf(Re(RP_m_temp[i][j])));
	RP_m_temp_I_Im = _double(Inf(Im(RP_m_temp[i][j])));
	RP_m_temp_S_Re = _double(Sup(Re(RP_m_temp[i][j])));
	RP_m_temp_S_Im = _double(Sup(Im(RP_m_temp[i][j])));

	RP_Yc_I_Re = RP_Yrc_I_Re + RP_m_temp_I_Re;
	RP_Yc_I_Im = RP_Yrc_I_Im + RP_m_temp_I_Im;
	RP_Yc_S_Re = RP_Yrc_S_Re + RP_m_temp_S_Re;
	RP_Yc_S_Im = RP_Yrc_S_Im + RP_m_temp_S_Im;

	Inf(Re(RP_Yc[i][j])) = _real(RP_Yc_I_Re);
	Inf(Im(RP_Yc[i][j])) = _real(RP_Yc_I_Im);
	Sup(Re(RP_Yc[i][j])) = _real(RP_Yc_S_Re);
	Sup(Im(RP_Yc[i][j])) = _real(RP_Yc_S_Im);
	}
}

//cout << "Yrc = " << RP_Yrc[0][0] << "\n";
//cout << "w * Cy = " << RP_m_temp[0][0] << "\n";
//cout << "Yc = " << RP_Yc[0][0] << "\n";

cinterval RP_Yc22[3][3], RP_Yc12[3][3], RP_Yc22_inv[3][3];

RP_Yc22[0][0] = RP_Yc[3][3];
RP_Yc22[0][1] = RP_Yc[3][4];
RP_Yc22[0][2] = RP_Yc[3][5];
RP_Yc22[1][0] = RP_Yc[4][3];
RP_Yc22[1][1] = RP_Yc[4][4];
RP_Yc22[1][2] = RP_Yc[4][5];
RP_Yc22[2][0] = RP_Yc[5][3];
RP_Yc22[2][1] = RP_Yc[5][4];
RP_Yc22[2][2] = RP_Yc[5][5];

RP_Yc12[0][0] = RP_Yc[3][0];
RP_Yc12[0][1] = RP_Yc[3][1];
RP_Yc12[0][2] = RP_Yc[3][2];
RP_Yc12[1][0] = RP_Yc[4][0];
RP_Yc12[1][1] = RP_Yc[4][1];
RP_Yc12[1][2] = RP_Yc[4][2];
RP_Yc12[2][0] = RP_Yc[5][0];
RP_Yc12[2][1] = RP_Yc[5][1];
RP_Yc12[2][2] = RP_Yc[5][2];
//cout << "Yc22 00"<< RP_Yc[3][3];// <-----------[-0.01]

ap::complex_2d_array Yc12_1, Yc12_2, Yc22_1, Yc22_2;
Yc12_1.setlength(3,3);
Yc12_2.setlength(3,3);
Yc22_1.setlength(3,3);
Yc22_2.setlength(3,3);
//cout << "noninv1" << RP_Yc22[0][0] << endl;//<----[-0.01]
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Yc22_1(i,j).x = _double(Inf(Re(RP_Yc22[i][j])));
	Yc22_1(i,j).y = _double(Inf(Im(RP_Yc22[i][j])));
	Yc22_2(i,j).x = _double(Sup(Re(RP_Yc22[i][j])));
	Yc22_2(i,j).y = _double(Sup(Im(RP_Yc22[i][j])));
	}
}
//cout << "noninv" << Yc22_1(0,0).y << endl;//<----[-0.01]
cmatrixinverse(Yc22_1, 3);
cmatrixinverse(Yc22_2, 3);
//!!!!!!!!!!
//cout << "inv" << Yc22_1(0,0).x << endl;
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Inf(Im(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).y,Yc22_2(i,j).y));
	Sup(Re(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Sup(Im(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).y,Yc22_2(i,j).y));
	}
}

/*
cout << "00" << RP_Yc22[0][0] << endl;
cout << "01" << RP_Yc22[0][1] << endl;
cout << "02" << RP_Yc22[0][2] << endl;
cout << "10" << RP_Yc22[1][0] << endl;
cout << "11" << RP_Yc22[1][1] << endl;
cout << "12" << RP_Yc22[1][2] << endl;
cout << "20" << RP_Yc22[2][0] << endl;
cout << "21" << RP_Yc22[2][1] << endl;
cout << "22" << RP_Yc22[2][2] << endl;
*/

cinterval RP_U1[3], RP_U2[3];
interval IP_arg_A, IP_arg_B, IP_arg_C, IP_180;
interval RP1_cos_A, RP1_cos_B, RP1_cos_C, RP1_sin_A, RP1_sin_B, RP1_sin_C;
interval RP2_cos_A, RP2_cos_B, RP2_cos_C, RP2_sin_A, RP2_sin_B, RP2_sin_C;

IP_arg_A = 0.0;
IP_arg_B = - 120.0;
IP_arg_C = 120.0;
//cout << "120" << IP_arg_B << "\n" << endl;
Inf(IP_180) = 180.0;
Sup(IP_180) = 180.0;
RP1_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP1_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP1_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP1_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP1_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP1_sin_C = sin((IP_arg_C / IP_180) * Pi());

RP2_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP2_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP2_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP2_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP2_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP2_sin_C = sin((IP_arg_C / IP_180) * Pi());
/*
interval fd;
Inf(fd) = 0.5;
Sup(fd) = 0.5;
interval fdd = cos(fd);
*/
//cout << "cos" << fdd << "\n" << endl;

//cout << "grad" << RP2_cos_B * Pi() << "\n" << endl;
//cout << "PI" << Pi() << "\n" << endl;
Re(RP_U1[0]) = IP_U_A_f * RP1_cos_A;
Im(RP_U1[0]) = IP_U_A_f * RP1_sin_A;
Re(RP_U1[1]) = IP_U_B_f * RP1_cos_B;
Im(RP_U1[1]) = IP_U_B_f * RP1_sin_B;
Re(RP_U1[2]) = IP_U_C_f * RP1_cos_C;
Im(RP_U1[2]) = IP_U_C_f * RP1_sin_C;

Re(RP_U2[0]) = IP_U_A_f * RP2_cos_A;
Im(RP_U2[0]) = IP_U_A_f * RP2_sin_A;
Re(RP_U2[1]) = IP_U_B_f * RP2_cos_B;
Im(RP_U2[1]) = IP_U_B_f * RP2_sin_B;
Re(RP_U2[2]) = IP_U_C_f * RP2_cos_C;
Im(RP_U2[2]) = IP_U_C_f * RP2_sin_C;

cinterval RP_I2[3], TP_SU[3];
cinterval TP_TMP1[3], TP_TMP2[3], TP_TMP3[3], TP_TMP4[3], TP_TMP5[3];

//--------------------------------------------------------------------------------------
ap::complex_1d_array temp1_1, temp2_1, temp3_1;
ap::complex_1d_array temp1_2, temp2_2, temp3_2;
ap::complex_1d_array U1_1, U1_2, U2_1, U2_2, I2_1, I2_2;
U1_1.setlength(3);
temp1_1.setlength(3);
temp2_1.setlength(3);
temp3_1.setlength(3);
I2_1.setlength(3);

U1_2.setlength(3);
temp1_2.setlength(3);
temp2_2.setlength(3);
temp3_2.setlength(3);
I2_2.setlength(3);

ap::complex_1d_array U456_1, U456_2;
U456_1.setlength(3);
U456_2.setlength(3);
ap::complex I4_1, I5_1, I6_1, U1_1_c, U2_1_c, U3_1_c, U4_1_c, U5_1_c, U6_1_c, Imax;
ap::complex I4_2, I5_2, I6_2, U1_2_c, U2_2_c, U3_2_c, U4_2_c, U5_2_c, U6_2_c;
ap::complex S;
S.x = _double(Inf(Re(IP_S)));
S.y = _double(Inf(Im(IP_S)));

ap::complex_2d_array Yc22_1_inv, Yc22_2_inv;
Yc22_1_inv.setlength(3,3);
Yc22_2_inv.setlength(3,3);

for(int i = 0 ; i < 3; i++)
{
	for(int j = 0 ; j < 3; j++)
	{
	Yc22_1_inv(i,j).x = _double(Inf(Re(RP_Yc22_inv[i][j])));
	Yc22_1_inv(i,j).y = _double(Inf(Im(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).x = _double(Sup(Re(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).y = _double(Sup(Im(RP_Yc22_inv[i][j])));

	Yc12_1(i,j).x = _double(Inf(Re(RP_Yc12[i][j])));
	Yc12_1(i,j).y = _double(Inf(Im(RP_Yc12[i][j])));
	Yc12_2(i,j).x = _double(Sup(Re(RP_Yc12[i][j])));
	Yc12_2(i,j).y = _double(Sup(Im(RP_Yc12[i][j])));
	}
}

U1_1_c.x = _double(Inf(Re(RP_U1[0]))); 
U1_1_c.y = _double(Inf(Im(RP_U1[0]))); 
U2_1_c.x = _double(Inf(Re(RP_U1[1]))); 
U2_1_c.y = _double(Inf(Im(RP_U1[1])));
U3_1_c.x = _double(Inf(Re(RP_U1[2]))); 
U3_1_c.y = _double(Inf(Im(RP_U1[2])));
U1_2_c.x = _double(Sup(Re(RP_U1[0]))); 
U1_2_c.y = _double(Sup(Im(RP_U1[0]))); 
U2_2_c.x = _double(Sup(Re(RP_U1[1]))); 
U2_2_c.y = _double(Sup(Im(RP_U1[1])));
U3_2_c.x = _double(Sup(Re(RP_U1[2]))); 
U3_2_c.y = _double(Sup(Im(RP_U1[2])));

U4_1_c.x = _double(Inf(Re(RP_U2[0]))); 
U4_1_c.y = _double(Inf(Im(RP_U2[0]))); 
U5_1_c.x = _double(Inf(Re(RP_U2[1]))); 
U5_1_c.y = _double(Inf(Im(RP_U2[1])));
U6_1_c.x = _double(Inf(Re(RP_U2[2]))); 
U6_1_c.y = _double(Inf(Im(RP_U2[2])));
U4_2_c.x = _double(Sup(Re(RP_U2[0]))); 
U4_2_c.y = _double(Sup(Im(RP_U2[0]))); 
U5_2_c.x = _double(Sup(Re(RP_U2[1]))); 
U5_2_c.y = _double(Sup(Im(RP_U2[1])));
U6_2_c.x = _double(Sup(Re(RP_U2[2]))); 
U6_2_c.y = _double(Sup(Im(RP_U2[2])));


//cout << "U4 = " << U4_2_c.x << endl << "\n";
// Нижнее значение <-------------------------------------------------------- Run process 1
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

I4_1.x = - (S / U4_1_c).x;
I4_1.y = (S / U4_1_c).y;
I5_1.x = - (S / U5_1_c).x;
I5_1.y = (S / U5_1_c).y;
I6_1.x = - (S / U6_1_c).x;
I6_1.y = (S / U6_1_c).y;

I2_1(0) = I4_1;
I2_1(1) = I5_1;
I2_1(2) = I6_1;

U1_1(0) = U1_1_c;
U1_1(1) = U2_1_c;
U1_1(2) = U3_1_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_1(j).x = 0;
  temp1_1(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yc12_1(j,i) * U1_1(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_1(i) = I2_1(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_1(j) += Yc22_1_inv(j,i) * temp2_1(i);
		}
	}

U4_1_c = temp3_1(0);
U5_1_c = temp3_1(1);
U6_1_c = temp3_1(2);

// Верхнее значение <-------------------------------------------------------- Run process 2

I4_2.x = - (S / U4_2_c).x;
I4_2.y = (S / U4_2_c).y;
I5_2.x = - (S / U5_2_c).x;
I5_2.y = (S / U5_2_c).y;
I6_2.x = - (S / U6_2_c).x;
I6_2.y = (S / U6_2_c).y;

I2_2(0) = I4_2;
I2_2(1) = I5_2;
I2_2(2) = I6_2;

U1_2(0) = U1_2_c;
U1_2(1) = U2_2_c;
U1_2(2) = U3_2_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yc12_2(j,i) * U1_2(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_2(i) = I2_2(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_2(j) += Yc22_2_inv(j,i) * temp2_2(i);
		}
	}

U4_2_c = temp3_2(0);
U5_2_c = temp3_2(1);
U6_2_c = temp3_2(2);



//cout << "Iteration = " << iter_count << endl << "\n";

for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}
/*
for(int i = 0 ; i < 3; i++)
{
cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}
*/
};



for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}


for(int i = 0 ; i < 3; i++)
{
//cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ap::complex_1d_array apRP_U2_low, apRP_U2_high, apRP_SU_low, apRP_SU_high;
apRP_U2_low.setlength(3);
apRP_U2_high.setlength(3);
apRP_SU_low.setlength(3);
apRP_SU_high.setlength(3);

for (int iter_count = 0; iter_count < iter_max; iter_count++)							// run process
{

for (int i = 0; i < 3; i++)
{
apRP_U2_low(i).x = _double(Inf(Re(RP_U2[i])));
apRP_U2_low(i).y = _double(Inf(Im(RP_U2[i])));
apRP_U2_high(i).x = _double(Sup(Re(RP_U2[i])));
apRP_U2_high(i).y = _double(Sup(Im(RP_U2[i])));
}

apRP_SU_low(0) = S / apRP_U2_low(0);
apRP_SU_low(1) = S / apRP_U2_low(1);
apRP_SU_low(2) = S / apRP_U2_low(2);
apRP_SU_high(0) = S / apRP_U2_high(0);
apRP_SU_high(1) = S / apRP_U2_high(1);
apRP_SU_high(2) = S / apRP_U2_high(2);

for (int i = 0; i < 3; i++)
{
Inf(Re(TP_SU[i])) = _real(min(apRP_SU_low(i).x,apRP_SU_high(i).x));
Inf(Im(TP_SU[i])) = _real(min(apRP_SU_low(i).y,apRP_SU_high(i).y));
Sup(Re(TP_SU[i])) = _real(max(apRP_SU_low(i).x,apRP_SU_high(i).x));
Sup(Im(TP_SU[i])) = _real(max(apRP_SU_low(i).y,apRP_SU_high(i).y));
}

	
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_U2_low [" << i << "] = [" << apRP_U2_low(i).x << "]+j*[" << apRP_U2_low(i).y << "]" << endl << "\n";
//cout << "apRP_U2_high [" << i << "] = [" << apRP_U2_high(i).x << "]+j*[" << apRP_U2_high(i).y << "]" << endl << "\n";
//}
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_SU_low [" << i << "] = [" << apRP_SU_low(i).x << "]+j*[" << apRP_SU_low(i).y << "]" << endl << "\n";
//cout << "apRP_SU_high [" << i << "] = [" << apRP_SU_high(i).x << "]+j*[" << apRP_SU_high(i).y << "]" << endl << "\n";
//}



Re(RP_I2[0]) = -1.0 * Re(TP_SU[0]); // ERR
Im(RP_I2[0]) = Im(TP_SU[0]); // OK
Re(RP_I2[1]) = -1.0 * Re(TP_SU[1]); // ERR
Im(RP_I2[1]) = Im(TP_SU[1]);  // OK
Re(RP_I2[2]) = -1.0 * Re(TP_SU[2]);  // OK
Im(RP_I2[2]) = Im(TP_SU[2]); // OK


double TP_TMP1_S_Re, TP_TMP1_S_Im, TP_TMP1_I_Re, TP_TMP1_I_Im;
double TP_TMP2_S_Re, TP_TMP2_S_Im, TP_TMP2_I_Re, TP_TMP2_I_Im;
double TP_TMP3_S_Re, TP_TMP3_S_Im, TP_TMP3_I_Re, TP_TMP3_I_Im;
double TP_TMP_S_Re, TP_TMP_S_Im, TP_TMP_I_Re, TP_TMP_I_Im;

for(int j = 0 ; j < 3; j++)
 {
 Re(TP_TMP2[j]) = 0.0;
 Im(TP_TMP2[j]) = 0.0;
 TP_TMP1[0] = RP_Yc12[j][0] * RP_U1[0];
 TP_TMP1[1] = RP_Yc12[j][1] * RP_U1[1];
 TP_TMP1[2] = RP_Yc12[j][2] * RP_U1[2];
 TP_TMP1_I_Re = _double(Inf(Re(TP_TMP1[0])));
 TP_TMP1_I_Im = _double(Inf(Im(TP_TMP1[0])));
 TP_TMP1_S_Re = _double(Sup(Re(TP_TMP1[0])));
 TP_TMP1_S_Im = _double(Sup(Im(TP_TMP1[0])));

 TP_TMP2_I_Re = _double(Inf(Re(TP_TMP1[1])));
 TP_TMP2_I_Im = _double(Inf(Im(TP_TMP1[1])));
 TP_TMP2_S_Re = _double(Sup(Re(TP_TMP1[1])));
 TP_TMP2_S_Im = _double(Sup(Im(TP_TMP1[1])));

 TP_TMP3_I_Re = _double(Inf(Re(TP_TMP1[2])));
 TP_TMP3_I_Im = _double(Inf(Im(TP_TMP1[2])));
 TP_TMP3_S_Re = _double(Sup(Re(TP_TMP1[2])));
 TP_TMP3_S_Im = _double(Sup(Im(TP_TMP1[2])));

 TP_TMP_I_Re = TP_TMP1_I_Re + TP_TMP2_I_Re + TP_TMP3_I_Re;
 TP_TMP_I_Im = TP_TMP1_I_Im + TP_TMP2_I_Im + TP_TMP3_I_Im;
 TP_TMP_S_Re = TP_TMP1_S_Re + TP_TMP2_S_Re + TP_TMP3_S_Re;
 TP_TMP_S_Im = TP_TMP1_S_Im + TP_TMP2_S_Im + TP_TMP3_S_Im;

 Inf(Re(TP_TMP2[j])) = _real(min(TP_TMP_I_Re,TP_TMP_S_Re));
 Inf(Im(TP_TMP2[j])) = _real(min(TP_TMP_I_Im,TP_TMP_S_Im));
 Sup(Re(TP_TMP2[j])) = _real(max(TP_TMP_I_Re,TP_TMP_S_Re));
 Sup(Im(TP_TMP2[j])) = _real(max(TP_TMP_I_Im,TP_TMP_S_Im));
 }


// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 3; i++)
	{
	TP_TMP3[i] = RP_I2[i] - TP_TMP2[i];
	}	
	
	for(int j = 0 ; j < 3; j++)
	{
	Inf(Re(TP_TMP4[j])) = real(0.0);
	Inf(Im(TP_TMP4[j])) = real(0.0);
	Sup(Re(TP_TMP4[j])) = real(0.0);
	Sup(Im(TP_TMP4[j])) = real(0.0);
	for(int i = 0 ; i < 3; i++)
	{
	TP_TMP4[j] += RP_Yc22_inv[j][i] * TP_TMP3[i];
	}
	}
	

RP_U2[0] = TP_TMP4[0];
RP_U2[1] = TP_TMP4[1];
RP_U2[2] = TP_TMP4[2];


//////////////////////////////////////////////////////////////////////////////////////////////
//cout << "Iteration = " << iter_count << endl << "\n";
double U4_mod_min, U4_mod_max, U5_mod_min, U5_mod_max, U6_mod_min, U6_mod_max;
U4_mod_min = _double(abs(Inf(TP_TMP4[0])));
U4_mod_max = _double(abs(Sup(TP_TMP4[0])));
U5_mod_min = _double(abs(Inf(TP_TMP4[1])));
U5_mod_max = _double(abs(Sup(TP_TMP4[1])));
U6_mod_min = _double(abs(Inf(TP_TMP4[2])));
U6_mod_max = _double(abs(Sup(TP_TMP4[2])));

//cout << "U4 = " << TP_TMP4[0] << endl << "\n";
//cout << "U5 = " << TP_TMP4[1] << endl << "\n";
//cout << "U6 = " << TP_TMP4[2] << endl << "\n";


cout << "mod U4 = [ " << U4_mod_min << " ; " << U4_mod_max << " ]" << endl << "\n";
cout << "mod U5 = [ " << U5_mod_min << " ; " << U5_mod_max << " ]" << endl << "\n";
cout << "mod U6 = [ " << U6_mod_min << " ; " << U6_mod_max << " ]" << endl << "\n";
//////////////////////////////////////////////////////////////////////////////////////////////
}															// stop process

}
//////////////////////////////////////////////////////////////////////////////////////////////
void get_nonsinusoidal()
{
cout << "Расчет несинусоидальных режимов.\n";

char* NS_Harmonic[14];
NS_Harmonic[0] = "1";
NS_Harmonic[1] = "5";
NS_Harmonic[2] = "7";
NS_Harmonic[3] = "11";
NS_Harmonic[4] = "13";
NS_Harmonic[5] = "17";
NS_Harmonic[6] = "19";
NS_Harmonic[7] = "23";
NS_Harmonic[8] = "25";
NS_Harmonic[9] = "29";
NS_Harmonic[10] = "31";
NS_Harmonic[11] = "35";
NS_Harmonic[12] = "37";
NS_Harmonic[13] = "41";

for (int i = 0; i < 14; i++)    
	{
	cout << "Расчет гармоники: " << NS_Harmonic[i] << ".\n";
	//process_nonsinusoidal();

};


cout << "OK.\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void getintvloss_fuzzy()
{

double MagConst;

interval IP_MagConst, IP_WC_temperature, IP_Yground, IP_freq, IP_cfreq, IP_WC_R0, IP_WC_S, IP_WC_R,
		 IP_WC_R_Equal, IP_X1, IP_X2, IP_X3, IP_Y1, IP_Y2, IP_Y3;

interval RP_d12, RP_d13, RP_d23, RP_D12, RP_D13, RP_D23;

interval IP_U_A_l, IP_U_B_l, IP_U_C_l, IP_U_A_f, IP_U_B_f, IP_U_C_f;

cinterval IP_S;

Re(IP_S) = 15.0 * 1000000.0;
Im(IP_S) = 7.5 * 1000000.0;

MagConst = 0.00000125663706;

IP_U_A_l = IP_Ua_l * 1000.0;
IP_U_B_l = IP_Ub_l * 1000.0;
IP_U_C_l = IP_Uc_l * 1000.0;
//cout << "U" << IP_U_A_l;
IP_U_A_f = IP_U_A_l / sqrt(3.0);
IP_U_B_f = IP_U_B_l / sqrt(3.0);
IP_U_C_f = IP_U_C_l / sqrt(3.0);


if (param_interval == 0)
{

IP_MagConst = MagConst;	// магнитная постоянная
IP_Yground = Yground;//(0.05,0.3); // проводимость земли
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
IP_L = Length;
IP_WC_R0 = 0.1197*(1.0 + 0.004*(temperature-20.0));
IP_WC_R0 = 0.2;		// сопротивление 1 км провода, См/м
IP_WC_S = 150.0;	// 240.0
IP_WC_R = 0.0085;	// 0.012, радиус провода
IP_X1 = X1;
IP_X2 = X2;
IP_X3 = X3;
IP_Y1 = Y1;
IP_Y2 = Y2;
IP_Y3 = Y3;
}

if (param_interval == 1)
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_WC_temperature) = -40.0;	// -40
Sup(IP_WC_temperature) = 40.0;	// 40
Inf(IP_Yground) = 0.05; // проводимость земли
Sup(IP_Yground) = 0.3;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.1197*(1.0 + 0.004*(IP_WC_temperature-20.0));
IP_WC_S = 150.0;	// 240.0
IP_WC_R = 0.0085;	// 0.012
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}

if (param_interval == 2)	// fazanord min & max parameters
{
IP_MagConst = 0.00000125663706;	// магнитная постоянная
Inf(IP_Yground) = 0.01; // проводимость земли
Sup(IP_Yground) = 0.01;
IP_freq = 50.0; // частота тока
IP_cfreq = 2.0 * Pi() * IP_freq; // круговая частота
Inf(IP_L) = ((1.0 - IP_L_mistake) * Length);
Sup(IP_L) = ((1.0 + IP_L_mistake) * Length);
IP_WC_R0 = 0.122;
IP_WC_S = 240.0;
IP_WC_R = 0.012;
Inf(IP_X1) = _real(0.975 * X1);
Sup(IP_X1) = _real(1.025 * X1);
Inf(IP_X2) = (0.975 * X2);
Sup(IP_X2) = (1.025 * X2);
Inf(IP_X3) = (0.975 * X3);
Sup(IP_X3) = (1.025 * X3);
Inf(IP_Y1) = (0.875 * Height);
Sup(IP_Y1) = (1.125 * Height);
Inf(IP_Y2) = (0.875 * Height);
Sup(IP_Y2) = (1.125 * Height);
Inf(IP_Y3) = (0.875 * Height);
Sup(IP_Y3) = (1.125 * Height);
}
//cout << IP_X1 << "\n";
//cout << IP_X2 << "\n";

//cout << "L = " << IP_L << " км" << "\n";
//cout << "L = " << IP_L << " км" << "\n";
//cout << IP_X2 << "\n";
interval IP_X12;
Inf(IP_X12) = Inf(IP_X1) - Inf(IP_X2);
Sup(IP_X12) = Sup(IP_X1) - Sup(IP_X2);
Inf(IP_X12) = Inf(IP_X12) * Inf(IP_X12);
Sup(IP_X12) = Sup(IP_X12) * Sup(IP_X12);

RP_d12 = sqrt(IP_X12 + Power(IP_Y1 - IP_Y2, 2.0));
RP_d13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 - IP_Y3), 2.0));
RP_d23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 - IP_Y3), 2.0));

RP_D12 = sqrt(IP_X12 + Power(IP_Y1 + IP_Y2, 2.0));
RP_D13 = sqrt(Power((IP_X1 - IP_X3), 2.0) + Power((IP_Y1 + IP_Y3), 2.0));
RP_D23 = sqrt(Power((IP_X2 - IP_X3), 2.0) + Power((IP_Y2 + IP_Y3), 2.0));
//RP_D12 = RP_D23; 
//cout << RP_D12 << "\n";
cinterval RP_Z11, RP_Z12, RP_Z13, RP_Z21, RP_Z22, RP_Z23, RP_Z31, RP_Z32, RP_Z33, RP_Z11out, RP_Z11inn;
interval pow_0_755, pow_0_83;
pow_0_755 = 0.755;
pow_0_83 = 0.83;

if (phase_fault == 0)
{

Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;
}

if (phase_fault == 1)
{
Re(RP_Z11out) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z11out) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (IP_WC_R * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
Re(RP_Z11inn) = IP_WC_R0 * (0.9 + 0.0063 * pow(IP_freq,pow_0_755));
Im(RP_Z11inn) = 0.001 * ((0.033 - 0.00107 * pow(IP_freq,pow_0_83)) * IP_WC_S + ( 1.07 * pow(IP_freq,pow_0_83) - 13.5 ));
RP_Z11 = RP_Z11out * 1000.0 + RP_Z11inn;


Re(RP_Z12) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z12) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d12 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z12 = RP_Z12 * 1000.0;

Re(RP_Z13) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z13) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d13 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z13 = RP_Z13 * 1000.0;

RP_Z21 = RP_Z12;

RP_Z22 = RP_Z11;

Re(RP_Z23) = IP_cfreq * IP_MagConst / 8.0;  
Im(RP_Z23) = (IP_cfreq * IP_MagConst) / ( 2.0 * Pi() ) * ln( 1.85 / (RP_d23 * sqrt(IP_Yground * IP_cfreq * IP_MagConst)));
RP_Z23 = RP_Z23 * 1000.0;

RP_Z31 = RP_Z13;

RP_Z32 = RP_Z23;

RP_Z33 = RP_Z11;

// -----------------------------------------------------------------------------------------
Re(RP_Z11) = 1000000.0;
Im(RP_Z11) = 0.0;
}

cinterval RP_d_p[3][3], RP_d_n[3][3];
RP_Z[0][0] = RP_Z11;
RP_Z[0][1] = RP_Z12;
RP_Z[0][2] = RP_Z13;
RP_Z[1][0] = RP_Z21;
RP_Z[1][1] = RP_Z22;
RP_Z[1][2] = RP_Z23;
RP_Z[2][0] = RP_Z31;
RP_Z[2][1] = RP_Z32;
RP_Z[2][2] = RP_Z33;


for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	RP_Z[i][j] = IP_L * RP_Z[i][j];

	}
}

//cout << RP_Z[0][0] << "\n";

ap::complex_2d_array Z_inf;
ap::complex_2d_array Z_sup;
ap::complex_2d_array D_inf;
ap::complex_2d_array D_sup;
Z_inf.setlength(3,3);
Z_sup.setlength(3,3);
D_inf.setlength(3,3);
D_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Z_inf(i,j).x = _double(Inf(Re(RP_Z[i][j])));
	Z_inf(i,j).y = _double(Inf(Im(RP_Z[i][j])));
	Z_sup(i,j).x = _double(Sup(Re(RP_Z[i][j])));
	Z_sup(i,j).y = _double(Sup(Im(RP_Z[i][j])));
	}
}
//cout << Z_inf(0,0).x << "\n";
cmatrixinverse(Z_inf,3);
cmatrixinverse(Z_sup,3);
//cout << Z_inf(0,0).x << "\n";
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_d_p[i][j])) = Z_inf(i,j).x;
	Inf(Im(RP_d_p[i][j])) = Z_inf(i,j).y;
	Sup(Re(RP_d_p[i][j])) = Z_sup(i,j).x;
	Sup(Im(RP_d_p[i][j])) = Z_sup(i,j).y;
	Inf(Re(RP_d_n[i][j])) = - Z_inf(i,j).x;// -Z
	Inf(Im(RP_d_n[i][j])) = - Z_inf(i,j).y;
	Sup(Re(RP_d_n[i][j])) = - Z_sup(i,j).x;
	Sup(Im(RP_d_n[i][j])) = - Z_sup(i,j).y;
	}
}

/*
//inv
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_Z[i][j]) = Inf(RP_Z[i][j]) / 100.0;
	Sup(RP_Z[i][j]) = Sup(RP_Z[i][j]) / 100.0;
	}
}

cinterval M_A[3][3], M_B[3][3], M_X[3][3], M_J[3][3];

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
		M_A[i][j] = RP_Z[i][j]; 
		Re(M_J[i][j]) = 0.0;
		Im(M_J[i][j]) = 0.0;
	}
}

Re(M_J[0][0]) = 1.0;
Im(M_J[0][0]) = 0.0;
Re(M_J[1][1]) = 1.0;
Im(M_J[1][1]) = 0.0;
Re(M_J[2][2]) = 1.0;
Im(M_J[2][2]) = 0.0;

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	M_B[i][j] = M_J[i][j] - M_A[i][j]; 
	}
}

cout << abs(M_B[0][0]) << "\n";


// норма комплексной матрицы

double M_B_norm_row, M_B_norm_col;

interval M_B_row_1, M_B_row_2, M_B_row_3;

//M_B_row_1;
for (int j=0; j<3;j++)
{
M_B_row_1 += abs(M_B[0][j]);
M_B_row_2 += abs(M_B[1][j]);
M_B_row_3 += abs(M_B[2][j]);
}

 
//cout << M_B_row_2 << "\n";

const int len = 3;
const int start = 0;
double arr[len];
arr[0] = _double(Sup(M_B_row_1));
arr[1] = _double(Sup(M_B_row_2));
arr[2] = _double(Sup(M_B_row_3));
double* first = arr;
double* last = arr + len;
double* asd;
asd = max_element(first,last);


cout << "Matrix norm" << *asd << "\n";

M_B_norm_col = 1.0 / (1.0 - *asd);


Inf(M_X[0][0]) = - M_B_norm_col;
Sup(M_X[0][0]) = M_B_norm_col + 2.0;
Inf(M_X[0][1]) = - M_B_norm_col;
Sup(M_X[0][1]) = M_B_norm_col;
Inf(M_X[0][2]) = - M_B_norm_col;
Sup(M_X[0][2]) = M_B_norm_col;
Inf(M_X[1][0]) = - M_B_norm_col;
Sup(M_X[1][0]) = M_B_norm_col;
Inf(M_X[1][1]) = - M_B_norm_col;
Sup(M_X[1][1]) = M_B_norm_col + 2.0;
Inf(M_X[1][2]) = - M_B_norm_col;
Sup(M_X[1][2]) = M_B_norm_col;
Inf(M_X[2][0]) = - M_B_norm_col;
Sup(M_X[2][0]) = M_B_norm_col;
Inf(M_X[2][1]) = - M_B_norm_col;
Sup(M_X[2][1]) = M_B_norm_col;
Inf(M_X[2][2]) = - M_B_norm_col;
Sup(M_X[2][2]) = M_B_norm_col + 2.0;
cout << "M_X" << M_X[0][0] << "\n";

cimatrix iM_X(3,3), iM_J(3,3), iM_A(3,3);

for (int i=0; i<3; i++)
{
	for (int j=0; j<3; j++)
	{
	iM_X(i+1,j+1) = M_X[i][j];
	}
}
cout << "M_X" << iM_X(1,1) << "\n";


	for (int i=0; i<10; i++)
	{
	M_X[i][j] = mid(M_X) + M_X * (M_J - M_A * mid(M_X));
	}
	
//cout << Re(M_X(1,1)) << "\n";
*/
// Yrc

cinterval RP_Yrc[6][6];
RP_Yrc[0][0] = RP_d_p[0][0];
RP_Yrc[0][1] = RP_d_p[0][1];
RP_Yrc[0][2] = RP_d_p[0][2];
RP_Yrc[0][3] = RP_d_n[0][3];
RP_Yrc[0][4] = RP_d_n[0][1];
RP_Yrc[0][5] = RP_d_n[0][2];
RP_Yrc[1][0] = RP_d_p[1][0];
RP_Yrc[1][1] = RP_d_p[1][1];
RP_Yrc[1][2] = RP_d_p[1][2];
RP_Yrc[1][3] = RP_d_n[1][0];
RP_Yrc[1][4] = RP_d_n[1][1];
RP_Yrc[1][5] = RP_d_n[1][2];
RP_Yrc[2][0] = RP_d_p[2][0];
RP_Yrc[2][1] = RP_d_p[2][1];
RP_Yrc[2][2] = RP_d_p[2][2];
RP_Yrc[2][3] = RP_d_n[2][0];
RP_Yrc[2][4] = RP_d_n[2][1];
RP_Yrc[2][5] = RP_d_n[2][2];
RP_Yrc[3][0] = RP_d_n[0][0];
RP_Yrc[3][1] = RP_d_n[0][1];
RP_Yrc[3][2] = RP_d_n[0][2];
RP_Yrc[3][3] = RP_d_p[0][0];
RP_Yrc[3][4] = RP_d_p[0][1];
RP_Yrc[3][5] = RP_d_p[0][2];
RP_Yrc[4][0] = RP_d_n[1][0];
RP_Yrc[4][1] = RP_d_n[1][1];
RP_Yrc[4][2] = RP_d_n[1][2];
RP_Yrc[4][3] = RP_d_p[1][0];
RP_Yrc[4][4] = RP_d_p[1][1];			
RP_Yrc[4][5] = RP_d_p[1][2];
RP_Yrc[5][0] = RP_d_n[2][0];
RP_Yrc[5][1] = RP_d_n[2][1];
RP_Yrc[5][2] = RP_d_n[2][2];
RP_Yrc[5][3] = RP_d_p[2][0];
RP_Yrc[5][4] = RP_d_p[2][1];
RP_Yrc[5][5] = RP_d_p[2][2];

const int n_ABC=3;
interval RP_A[n_ABC][n_ABC], RP_B[n_ABC][n_ABC], RP_C[n_ABC][n_ABC];

RP_A[0][0] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y1) / (IP_WC_R * 100.0));
RP_A[0][1] = 1.8 * power(10.0, 7.0) * ln(RP_D12 / RP_d12);
RP_A[0][2] = 1.8 * power(10.0, 7.0) * ln(RP_D13 / RP_d13);
RP_A[1][0] = RP_A[0][1];
RP_A[1][1] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y2) / (IP_WC_R * 100.0));
RP_A[1][2] = 1.8 * power(10.0, 7.0) * ln(RP_D23 / RP_d23);
RP_A[2][0] = RP_A[0][2];
RP_A[2][1] = RP_A[1][2];
RP_A[2][2] = 1.8 * power(10.0, 7.0) * ln((200.0 * IP_Y3) / (IP_WC_R * 100.0));

ap::real_2d_array a_inf;
ap::real_2d_array a_sup;
a_inf.setlength(3,3);
a_sup.setlength(3,3);

for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	a_inf(i,j) = _double(Inf(RP_A[i][j]));
	a_sup(i,j) = _double(Sup(RP_A[i][j]));
	}
}
//cout << a_inf(0,0) << "\n";
rmatrixinverse(a_inf,3);
rmatrixinverse(a_sup,3);
for (int i=0;i<3;i++)
{
	for (int j=0;j<3;j++)
	{
	Inf(RP_B[i][j]) = a_inf(i,j);
	Sup(RP_B[i][j]) = a_sup(i,j);
	}
}

/*
RP_C[0][0] = RP_B[0][0] + RP_B[0][1] + RP_B[0][2];
RP_C[0][1] = - RP_B[0][1];
RP_C[0][2] = - RP_B[0][2];
RP_C[1][0] = - RP_B[1][0];
RP_C[1][1] = RP_B[1][0] + RP_B[1][1] + RP_B[1][2];
RP_C[1][2] = - RP_B[1][2];
RP_C[2][0] = - RP_B[2][0];
RP_C[2][1] = - RP_B[2][1];
RP_C[2][2] = RP_B[2][0] + RP_B[2][1] + RP_B[2][2];
*/

//const int n_CY=6;
interval RP_Cy[6][6];
cinterval RP_Yc_inv[6][6]; 
interval i_zero(0.0);

RP_Cy[0][0] = RP_B[0][0];
RP_Cy[0][1] = RP_B[0][1];
RP_Cy[0][2] = RP_B[0][2];
RP_Cy[0][3] = i_zero;
RP_Cy[0][4] = i_zero;
RP_Cy[0][5] = i_zero;
RP_Cy[1][0] = RP_B[1][0];
RP_Cy[1][1] = RP_B[1][1];
RP_Cy[1][2] = RP_B[1][2];
RP_Cy[1][3] = i_zero;
RP_Cy[1][4] = i_zero;
RP_Cy[1][5] = i_zero;
RP_Cy[2][0] = RP_B[2][0];
RP_Cy[2][1] = RP_B[2][1];
RP_Cy[2][2] = RP_B[2][2];
RP_Cy[2][3] = i_zero;
RP_Cy[2][4] = i_zero; 
RP_Cy[2][5] = i_zero;
RP_Cy[3][0] = i_zero;
RP_Cy[3][1] = i_zero;
RP_Cy[3][2] = i_zero;
RP_Cy[3][3] = RP_B[0][0];
RP_Cy[3][4] = RP_B[0][1];
RP_Cy[3][5] = RP_B[0][2];
RP_Cy[4][0] = i_zero;
RP_Cy[4][1] = i_zero;
RP_Cy[4][2] = i_zero;
RP_Cy[4][3] = RP_B[1][0];
RP_Cy[4][4] = RP_B[1][1];
RP_Cy[4][5] = RP_B[1][2];
RP_Cy[5][0] = i_zero;
RP_Cy[5][1] = i_zero;
RP_Cy[5][2] = i_zero;
RP_Cy[5][3] = RP_B[2][0];
RP_Cy[5][4] = RP_B[2][1];
RP_Cy[5][5] = RP_B[2][2];

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Cy[i][j] = RP_Cy[i][j] * 0.5 * IP_L * IP_cfreq;
	//cout << "Cy" << RP_Cy[i][j] << "\n";
	}
}

cinterval RP_m_temp[6][6];
for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Re(RP_m_temp[i][j]) = i_zero;
	Im(RP_m_temp[i][j]) = i_zero;
	Re(RP_Yc[i][j]) = i_zero;
	Im(RP_Yc[i][j]) = i_zero;
	}
}


for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	Im(RP_m_temp[i][j]) = RP_Cy[i][j];
	}
}

double RP_Yc_I_Re, RP_Yc_I_Im, RP_Yc_S_Re, RP_Yc_S_Im; 
double RP_Yrc_I_Re, RP_Yrc_I_Im, RP_Yrc_S_Re, RP_Yrc_S_Im; 
double RP_m_temp_I_Re, RP_m_temp_I_Im, RP_m_temp_S_Re, RP_m_temp_S_Im;

for (int i=0;i<6;i++)
{
	for (int j=0;j<6;j++)
	{
	RP_Yrc_I_Re = _double(Inf(Re(RP_Yrc[i][j])));
	RP_Yrc_I_Im = _double(Inf(Im(RP_Yrc[i][j])));
	RP_Yrc_S_Re = _double(Sup(Re(RP_Yrc[i][j])));
	RP_Yrc_S_Im = _double(Sup(Im(RP_Yrc[i][j])));
	
	RP_m_temp_I_Re = _double(Inf(Re(RP_m_temp[i][j])));
	RP_m_temp_I_Im = _double(Inf(Im(RP_m_temp[i][j])));
	RP_m_temp_S_Re = _double(Sup(Re(RP_m_temp[i][j])));
	RP_m_temp_S_Im = _double(Sup(Im(RP_m_temp[i][j])));

	RP_Yc_I_Re = RP_Yrc_I_Re + RP_m_temp_I_Re;
	RP_Yc_I_Im = RP_Yrc_I_Im + RP_m_temp_I_Im;
	RP_Yc_S_Re = RP_Yrc_S_Re + RP_m_temp_S_Re;
	RP_Yc_S_Im = RP_Yrc_S_Im + RP_m_temp_S_Im;

	Inf(Re(RP_Yc[i][j])) = _real(RP_Yc_I_Re);
	Inf(Im(RP_Yc[i][j])) = _real(RP_Yc_I_Im);
	Sup(Re(RP_Yc[i][j])) = _real(RP_Yc_S_Re);
	Sup(Im(RP_Yc[i][j])) = _real(RP_Yc_S_Im);
	}
}

//cout << "Yrc = " << RP_Yrc[0][0] << "\n";
//cout << "w * Cy = " << RP_m_temp[0][0] << "\n";
//cout << "Yc = " << RP_Yc[0][0] << "\n";

cinterval RP_Yc22[3][3], RP_Yc12[3][3], RP_Yc22_inv[3][3];

RP_Yc22[0][0] = RP_Yc[3][3];
RP_Yc22[0][1] = RP_Yc[3][4];
RP_Yc22[0][2] = RP_Yc[3][5];
RP_Yc22[1][0] = RP_Yc[4][3];
RP_Yc22[1][1] = RP_Yc[4][4];
RP_Yc22[1][2] = RP_Yc[4][5];
RP_Yc22[2][0] = RP_Yc[5][3];
RP_Yc22[2][1] = RP_Yc[5][4];
RP_Yc22[2][2] = RP_Yc[5][5];

RP_Yc12[0][0] = RP_Yc[3][0];
RP_Yc12[0][1] = RP_Yc[3][1];
RP_Yc12[0][2] = RP_Yc[3][2];
RP_Yc12[1][0] = RP_Yc[4][0];
RP_Yc12[1][1] = RP_Yc[4][1];
RP_Yc12[1][2] = RP_Yc[4][2];
RP_Yc12[2][0] = RP_Yc[5][0];
RP_Yc12[2][1] = RP_Yc[5][1];
RP_Yc12[2][2] = RP_Yc[5][2];
//cout << "Yc22 00"<< RP_Yc[3][3];// <-----------[-0.01]

ap::complex_2d_array Yc12_1, Yc12_2, Yc22_1, Yc22_2;
Yc12_1.setlength(3,3);
Yc12_2.setlength(3,3);
Yc22_1.setlength(3,3);
Yc22_2.setlength(3,3);
//cout << "noninv1" << RP_Yc22[0][0] << endl;//<----[-0.01]
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Yc22_1(i,j).x = _double(Inf(Re(RP_Yc22[i][j])));
	Yc22_1(i,j).y = _double(Inf(Im(RP_Yc22[i][j])));
	Yc22_2(i,j).x = _double(Sup(Re(RP_Yc22[i][j])));
	Yc22_2(i,j).y = _double(Sup(Im(RP_Yc22[i][j])));
	}
}
//cout << "noninv" << Yc22_1(0,0).y << endl;//<----[-0.01]
cmatrixinverse(Yc22_1, 3);
cmatrixinverse(Yc22_2, 3);
//!!!!!!!!!!
//cout << "inv" << Yc22_1(0,0).x << endl;
for (int i=0;i <3;i++)
{   
	for (int j=0;j<3;j++)
	{
	Inf(Re(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Inf(Im(RP_Yc22_inv[i][j])) = _real(min(Yc22_1(i,j).y,Yc22_2(i,j).y));
	Sup(Re(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).x,Yc22_2(i,j).x));
	Sup(Im(RP_Yc22_inv[i][j])) = _real(max(Yc22_1(i,j).y,Yc22_2(i,j).y));
	}
}

/*
cout << "00" << RP_Yc22[0][0] << endl;
cout << "01" << RP_Yc22[0][1] << endl;
cout << "02" << RP_Yc22[0][2] << endl;
cout << "10" << RP_Yc22[1][0] << endl;
cout << "11" << RP_Yc22[1][1] << endl;
cout << "12" << RP_Yc22[1][2] << endl;
cout << "20" << RP_Yc22[2][0] << endl;
cout << "21" << RP_Yc22[2][1] << endl;
cout << "22" << RP_Yc22[2][2] << endl;
*/

cinterval RP_U1[3], RP_U2[3];
interval IP_arg_A, IP_arg_B, IP_arg_C, IP_180;
interval RP1_cos_A, RP1_cos_B, RP1_cos_C, RP1_sin_A, RP1_sin_B, RP1_sin_C;
interval RP2_cos_A, RP2_cos_B, RP2_cos_C, RP2_sin_A, RP2_sin_B, RP2_sin_C;

IP_arg_A = 0.0;
IP_arg_B = - 120.0;
IP_arg_C = 120.0;
//cout << "120" << IP_arg_B << "\n" << endl;
Inf(IP_180) = 180.0;
Sup(IP_180) = 180.0;
RP1_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP1_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP1_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP1_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP1_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP1_sin_C = sin((IP_arg_C / IP_180) * Pi());

RP2_cos_A = cos((IP_arg_A / IP_180) * Pi());
RP2_cos_B = cos((IP_arg_B / IP_180) * Pi());
RP2_cos_C = cos((IP_arg_C / IP_180) * Pi());
RP2_sin_A = sin((IP_arg_A / IP_180) * Pi());
RP2_sin_B = sin((IP_arg_B / IP_180) * Pi());
RP2_sin_C = sin((IP_arg_C / IP_180) * Pi());
/*
interval fd;
Inf(fd) = 0.5;
Sup(fd) = 0.5;
interval fdd = cos(fd);
*/
//cout << "cos" << fdd << "\n" << endl;

//cout << "grad" << RP2_cos_B * Pi() << "\n" << endl;
//cout << "PI" << Pi() << "\n" << endl;
Re(RP_U1[0]) = IP_U_A_f * RP1_cos_A;
Im(RP_U1[0]) = IP_U_A_f * RP1_sin_A;
Re(RP_U1[1]) = IP_U_B_f * RP1_cos_B;
Im(RP_U1[1]) = IP_U_B_f * RP1_sin_B;
Re(RP_U1[2]) = IP_U_C_f * RP1_cos_C;
Im(RP_U1[2]) = IP_U_C_f * RP1_sin_C;

Re(RP_U2[0]) = IP_U_A_f * RP2_cos_A;
Im(RP_U2[0]) = IP_U_A_f * RP2_sin_A;
Re(RP_U2[1]) = IP_U_B_f * RP2_cos_B;
Im(RP_U2[1]) = IP_U_B_f * RP2_sin_B;
Re(RP_U2[2]) = IP_U_C_f * RP2_cos_C;
Im(RP_U2[2]) = IP_U_C_f * RP2_sin_C;

cinterval RP_I2[3], TP_SU[3];
cinterval TP_TMP1[3], TP_TMP2[3], TP_TMP3[3], TP_TMP4[3], TP_TMP5[3];

//--------------------------------------------------------------------------------------
ap::complex_1d_array temp1_1, temp2_1, temp3_1;
ap::complex_1d_array temp1_2, temp2_2, temp3_2;
ap::complex_1d_array U1_1, U1_2, U2_1, U2_2, I2_1, I2_2;
U1_1.setlength(3);
temp1_1.setlength(3);
temp2_1.setlength(3);
temp3_1.setlength(3);
I2_1.setlength(3);

U1_2.setlength(3);
temp1_2.setlength(3);
temp2_2.setlength(3);
temp3_2.setlength(3);
I2_2.setlength(3);

ap::complex_1d_array U456_1, U456_2;
U456_1.setlength(3);
U456_2.setlength(3);
ap::complex I4_1, I5_1, I6_1, U1_1_c, U2_1_c, U3_1_c, U4_1_c, U5_1_c, U6_1_c, Imax;
ap::complex I4_2, I5_2, I6_2, U1_2_c, U2_2_c, U3_2_c, U4_2_c, U5_2_c, U6_2_c;
ap::complex S;
S.x = _double(Inf(Re(IP_S)));
S.y = _double(Inf(Im(IP_S)));

ap::complex_2d_array Yc22_1_inv, Yc22_2_inv;
Yc22_1_inv.setlength(3,3);
Yc22_2_inv.setlength(3,3);

for(int i = 0 ; i < 3; i++)
{
	for(int j = 0 ; j < 3; j++)
	{
	Yc22_1_inv(i,j).x = _double(Inf(Re(RP_Yc22_inv[i][j])));
	Yc22_1_inv(i,j).y = _double(Inf(Im(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).x = _double(Sup(Re(RP_Yc22_inv[i][j])));
	Yc22_2_inv(i,j).y = _double(Sup(Im(RP_Yc22_inv[i][j])));

	Yc12_1(i,j).x = _double(Inf(Re(RP_Yc12[i][j])));
	Yc12_1(i,j).y = _double(Inf(Im(RP_Yc12[i][j])));
	Yc12_2(i,j).x = _double(Sup(Re(RP_Yc12[i][j])));
	Yc12_2(i,j).y = _double(Sup(Im(RP_Yc12[i][j])));
	}
}

U1_1_c.x = _double(Inf(Re(RP_U1[0]))); 
U1_1_c.y = _double(Inf(Im(RP_U1[0]))); 
U2_1_c.x = _double(Inf(Re(RP_U1[1]))); 
U2_1_c.y = _double(Inf(Im(RP_U1[1])));
U3_1_c.x = _double(Inf(Re(RP_U1[2]))); 
U3_1_c.y = _double(Inf(Im(RP_U1[2])));
U1_2_c.x = _double(Sup(Re(RP_U1[0]))); 
U1_2_c.y = _double(Sup(Im(RP_U1[0]))); 
U2_2_c.x = _double(Sup(Re(RP_U1[1]))); 
U2_2_c.y = _double(Sup(Im(RP_U1[1])));
U3_2_c.x = _double(Sup(Re(RP_U1[2]))); 
U3_2_c.y = _double(Sup(Im(RP_U1[2])));

U4_1_c.x = _double(Inf(Re(RP_U2[0]))); 
U4_1_c.y = _double(Inf(Im(RP_U2[0]))); 
U5_1_c.x = _double(Inf(Re(RP_U2[1]))); 
U5_1_c.y = _double(Inf(Im(RP_U2[1])));
U6_1_c.x = _double(Inf(Re(RP_U2[2]))); 
U6_1_c.y = _double(Inf(Im(RP_U2[2])));
U4_2_c.x = _double(Sup(Re(RP_U2[0]))); 
U4_2_c.y = _double(Sup(Im(RP_U2[0]))); 
U5_2_c.x = _double(Sup(Re(RP_U2[1]))); 
U5_2_c.y = _double(Sup(Im(RP_U2[1])));
U6_2_c.x = _double(Sup(Re(RP_U2[2]))); 
U6_2_c.y = _double(Sup(Im(RP_U2[2])));


//cout << "U4 = " << U4_2_c.x << endl << "\n";
// Нижнее значение <-------------------------------------------------------- Run process 1
for (int iter_count = 0; iter_count < iter_max; iter_count++)
{

I4_1.x = - (S / U4_1_c).x;
I4_1.y = (S / U4_1_c).y;
I5_1.x = - (S / U5_1_c).x;
I5_1.y = (S / U5_1_c).y;
I6_1.x = - (S / U6_1_c).x;
I6_1.y = (S / U6_1_c).y;

I2_1(0) = I4_1;
I2_1(1) = I5_1;
I2_1(2) = I6_1;

U1_1(0) = U1_1_c;
U1_1(1) = U2_1_c;
U1_1(2) = U3_1_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_1(j).x = 0;
  temp1_1(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_1(j) += Yc12_1(j,i) * U1_1(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_1(i) = I2_1(i) - temp1_1(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_1(j).x = 0;
	temp3_1(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_1(j) += Yc22_1_inv(j,i) * temp2_1(i);
		}
	}

U4_1_c = temp3_1(0);
U5_1_c = temp3_1(1);
U6_1_c = temp3_1(2);

// Верхнее значение <-------------------------------------------------------- Run process 2

I4_2.x = - (S / U4_2_c).x;
I4_2.y = (S / U4_2_c).y;
I5_2.x = - (S / U5_2_c).x;
I5_2.y = (S / U5_2_c).y;
I6_2.x = - (S / U6_2_c).x;
I6_2.y = (S / U6_2_c).y;

I2_2(0) = I4_2;
I2_2(1) = I5_2;
I2_2(2) = I6_2;

U1_2(0) = U1_2_c;
U1_2(1) = U2_2_c;
U1_2(2) = U3_2_c;

for(int j = 0 ; j < 3; j++)
 {
  temp1_2(j).x = 0;
  temp1_2(j).y = 0;
  for(int i = 0; i < 3; i++)
   temp1_2(j) += Yc12_2(j,i) * U1_2(i);
 }

// Вычитание I2 - [ Yc(12) * U(1) ]

	for (int i=0; i < 3; i++)
	{
	temp2_2(i) = I2_2(i) - temp1_2(i);
	}

// Перемножение Yc22_inv * temp2

for(int j = 0 ; j < 3; j++)
	{
	temp3_2(j).x = 0;
	temp3_2(j).y = 0;
		for(int i = 0; i < 3; i++)
		{
		temp3_2(j) += Yc22_2_inv(j,i) * temp2_2(i);
		}
	}

U4_2_c = temp3_2(0);
U5_2_c = temp3_2(1);
U6_2_c = temp3_2(2);



//cout << "Iteration = " << iter_count << endl << "\n";

for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}
/*
for(int i = 0 ; i < 3; i++)
{
cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}
*/
};



for(int i = 0 ; i < 3; i++)
{
RP_Umod_min[i] = min(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
RP_Umod_max[i] = max(sqrt(pow((temp3_1(i).x), 2.0) + pow((temp3_1(i).y), 2.0)), sqrt(pow((temp3_2(i).x), 2.0) + pow((temp3_2(i).y), 2.0)));
}


for(int i = 0 ; i < 3; i++)
{
cout << "U [" << i << "] = [" << RP_Umod_min[i] << ", " << RP_Umod_max[i] << "]" << endl << "\n";
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
ap::complex_1d_array apRP_U2_low, apRP_U2_high, apRP_SU_low, apRP_SU_high;
apRP_U2_low.setlength(3);
apRP_U2_high.setlength(3);
apRP_SU_low.setlength(3);
apRP_SU_high.setlength(3);

for (int iter_count = 0; iter_count < iter_max; iter_count++)							// run process
{

for (int i = 0; i < 3; i++)
{
apRP_U2_low(i).x = _double(Inf(Re(RP_U2[i])));
apRP_U2_low(i).y = _double(Inf(Im(RP_U2[i])));
apRP_U2_high(i).x = _double(Sup(Re(RP_U2[i])));
apRP_U2_high(i).y = _double(Sup(Im(RP_U2[i])));
}

apRP_SU_low(0) = S / apRP_U2_low(0);
apRP_SU_low(1) = S / apRP_U2_low(1);
apRP_SU_low(2) = S / apRP_U2_low(2);
apRP_SU_high(0) = S / apRP_U2_high(0);
apRP_SU_high(1) = S / apRP_U2_high(1);
apRP_SU_high(2) = S / apRP_U2_high(2);

for (int i = 0; i < 3; i++)
{
Inf(Re(TP_SU[i])) = _real(min(apRP_SU_low(i).x,apRP_SU_high(i).x));
Inf(Im(TP_SU[i])) = _real(min(apRP_SU_low(i).y,apRP_SU_high(i).y));
Sup(Re(TP_SU[i])) = _real(max(apRP_SU_low(i).x,apRP_SU_high(i).x));
Sup(Im(TP_SU[i])) = _real(max(apRP_SU_low(i).y,apRP_SU_high(i).y));
}

	
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_U2_low [" << i << "] = [" << apRP_U2_low(i).x << "]+j*[" << apRP_U2_low(i).y << "]" << endl << "\n";
//cout << "apRP_U2_high [" << i << "] = [" << apRP_U2_high(i).x << "]+j*[" << apRP_U2_high(i).y << "]" << endl << "\n";
//}
//for(int i = 0 ; i < 3; i++)
//{
//cout << "apRP_SU_low [" << i << "] = [" << apRP_SU_low(i).x << "]+j*[" << apRP_SU_low(i).y << "]" << endl << "\n";
//cout << "apRP_SU_high [" << i << "] = [" << apRP_SU_high(i).x << "]+j*[" << apRP_SU_high(i).y << "]" << endl << "\n";
//}



Re(RP_I2[0]) = -1.0 * Re(TP_SU[0]); // ERR
Im(RP_I2[0]) = Im(TP_SU[0]); // OK
Re(RP_I2[1]) = -1.0 * Re(TP_SU[1]); // ERR
Im(RP_I2[1]) = Im(TP_SU[1]);  // OK
Re(RP_I2[2]) = -1.0 * Re(TP_SU[2]);  // OK
Im(RP_I2[2]) = Im(TP_SU[2]); // OK


double TP_TMP1_S_Re, TP_TMP1_S_Im, TP_TMP1_I_Re, TP_TMP1_I_Im;
double TP_TMP2_S_Re, TP_TMP2_S_Im, TP_TMP2_I_Re, TP_TMP2_I_Im;
double TP_TMP3_S_Re, TP_TMP3_S_Im, TP_TMP3_I_Re, TP_TMP3_I_Im;
double TP_TMP_S_Re, TP_TMP_S_Im, TP_TMP_I_Re, TP_TMP_I_Im;

for(int j = 0 ; j < 3; j++)
 {
 Re(TP_TMP2[j]) = 0.0;
 Im(TP_TMP2[j]) = 0.0;
 TP_TMP1[0] = RP_Yc12[j][0] * RP_U1[0];
 TP_TMP1[1] = RP_Yc12[j][1] * RP_U1[1];
 TP_TMP1[2] = RP_Yc12[j][2] * RP_U1[2];
 TP_TMP1_I_Re = _double(Inf(Re(TP_TMP1[0])));
 TP_TMP1_I_Im = _double(Inf(Im(TP_TMP1[0])));
 TP_TMP1_S_Re = _double(Sup(Re(TP_TMP1[0])));
 TP_TMP1_S_Im = _double(Sup(Im(TP_TMP1[0])));

 TP_TMP2_I_Re = _double(Inf(Re(TP_TMP1[1])));
 TP_TMP2_I_Im = _double(Inf(Im(TP_TMP1[1])));
 TP_TMP2_S_Re = _double(Sup(Re(TP_TMP1[1])));
 TP_TMP2_S_Im = _double(Sup(Im(TP_TMP1[1])));

 TP_TMP3_I_Re = _double(Inf(Re(TP_TMP1[2])));
 TP_TMP3_I_Im = _double(Inf(Im(TP_TMP1[2])));
 TP_TMP3_S_Re = _double(Sup(Re(TP_TMP1[2])));
 TP_TMP3_S_Im = _double(Sup(Im(TP_TMP1[2])));

 TP_TMP_I_Re = TP_TMP1_I_Re + TP_TMP2_I_Re + TP_TMP3_I_Re;
 TP_TMP_I_Im = TP_TMP1_I_Im + TP_TMP2_I_Im + TP_TMP3_I_Im;
 TP_TMP_S_Re = TP_TMP1_S_Re + TP_TMP2_S_Re + TP_TMP3_S_Re;
 TP_TMP_S_Im = TP_TMP1_S_Im + TP_TMP2_S_Im + TP_TMP3_S_Im;

 Inf(Re(TP_TMP2[j])) = _real(min(TP_TMP_I_Re,TP_TMP_S_Re));
 Inf(Im(TP_TMP2[j])) = _real(min(TP_TMP_I_Im,TP_TMP_S_Im));
 Sup(Re(TP_TMP2[j])) = _real(max(TP_TMP_I_Re,TP_TMP_S_Re));
 Sup(Im(TP_TMP2[j])) = _real(max(TP_TMP_I_Im,TP_TMP_S_Im));
 }


// Вычитание I2 - [ Yc(12) * U(1) ]
	for (int i=0; i < 3; i++)
	{
	TP_TMP3[i] = RP_I2[i] - TP_TMP2[i];
	}	
	
	for(int j = 0 ; j < 3; j++)
	{
	Inf(Re(TP_TMP4[j])) = real(0.0);
	Inf(Im(TP_TMP4[j])) = real(0.0);
	Sup(Re(TP_TMP4[j])) = real(0.0);
	Sup(Im(TP_TMP4[j])) = real(0.0);
	for(int i = 0 ; i < 3; i++)
	{
	TP_TMP4[j] += RP_Yc22_inv[j][i] * TP_TMP3[i];
	}
	}
	

RP_U2[0] = TP_TMP4[0];
RP_U2[1] = TP_TMP4[1];
RP_U2[2] = TP_TMP4[2];


//////////////////////////////////////////////////////////////////////////////////////////////
cout << "Iteration = " << iter_count << endl << "\n";

U4_mod_min = _double(abs(Inf(TP_TMP4[0])));
U4_mod_max = _double(abs(Sup(TP_TMP4[0])));
U5_mod_min = _double(abs(Inf(TP_TMP4[1])));
U5_mod_max = _double(abs(Sup(TP_TMP4[1])));
U6_mod_min = _double(abs(Inf(TP_TMP4[2])));
U6_mod_max = _double(abs(Sup(TP_TMP4[2])));

cout << "U4 = " << TP_TMP4[0] << endl << "\n";
cout << "U5 = " << TP_TMP4[1] << endl << "\n";
cout << "U6 = " << TP_TMP4[2] << endl << "\n";


cout << "mod U4 = [ " << U4_mod_min << " ; " << U4_mod_max << " ]" << endl << "\n";
cout << "mod U5 = [ " << U5_mod_min << " ; " << U5_mod_max << " ]" << endl << "\n";
cout << "mod U6 = [ " << U6_mod_min << " ; " << U6_mod_max << " ]" << endl << "\n";
//////////////////////////////////////////////////////////////////////////////////////////////


}															// stop process
*/
}

//////////////////////////////////////////////////////////////////////////////////////////////
void get_fuzzy_interval()
{
cout << "Интервально-нечеткое моделирование.\n";
getintvloss_fuzzy();

cout << "OK.\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void init_new()
{
cout << "Выберите расчет:\n";
cout << "1. Расчет участка схемы.\n";
cout << "2. Расчет схемы.\n";
cout << "3. Расчет короткого замыкания.\n";
cout << "4. Расчет неполнофазных режимов.\n";
cout << "5. Расчет несимметричных режимов.\n";
cout << "6. Моделирование электромагнитного поля.\n";
cout << "7. Дополнительные материалы.\n";
cout << "8. Выход.\n";
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(void)
{
MF_U_low.setlength(9);
MF_U_high.setlength(9);
MF_ZM_low.setlength(3,3);
MF_ZM_high.setlength(3,3);

setlocale (LC_ALL,"Russian");
cout << "Interval Calculation 2.0 Console\n";
cout << "Uses C-XSC, Alglib\n";
cout << "Code: Alexander Litvintsev, IrGUPS\n";
cout << "Uses ALGLIB Open Source\n";
cout << "Uses CXSC Interval library\n";
cout << "Build: 25.01.2015\n";
cout << "\n";
param_interval=1;	// 0 - non-interval, 1 - interval
param_unsim=0;		// 0 - normal, 1 - non-simmetric

	//getintvloss();
	//setintvloss_triangle();
	//setintvloss_star();
	//getintvalue_00();
	//setintvloss_3line();
	//get_K2U();
	//get_shortcircuit_2f_sim();
	//get_unsimmetric_fcoor();
	//get_shortcircuit_2f_fcoor();
	//get_shortcircuit_1f_fcoor_ground();
	//get_shortcircuit_2f_fcoor_ground();
	//setintvloss_3line_fault();
	//setintvloss_triangle_fault();
	//get_magnetic_field();
	//get_dls_cramer();
	//get_shortcircuit_mis();
	//get_ils();	// Решение ИСЛАУ
	//get_bar();	// Брусы
	//process_nonsinusoidal();
	//get_nonsinusoidal();	// Расчет несинусоидальных режимов
	  get_fuzzy_interval();	// Интервально-нечеткое моделирование
	  
	//init_new();	// Меню

system("pause");

setlocale(LC_ALL, ".OCP");
return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

