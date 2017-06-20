/************************************************
���� ������ ������������ ������������ AlgoPascal.
************************************************/

#include "ap.h"

void ludecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots);
void ludecompositionunpacked(ap::real_2d_array a,
     int m,
     int n,
     ap::real_2d_array& l,
     ap::real_2d_array& u,
     ap::integer_1d_array& pivots);

/*************************************************************************
LU-���������� ������� ������ ���� ������� M x N

������������ ��������� LU-���������� ������������� ������� ������  ����  �
��������� ������� �������� �������� (� �������������� �����).

������� ���������:
    A       -   ������� A. ��������� ���������: [1..M, 1..N]
    M       -   ����� ����� � ������� A
    N       -   ����� �������� � ������� A

�������� ���������:
    A       -   ������� L � U � ���������� ����� (��. ����).
                ��������� ���������: [1..M, 1..N]
    Pivots  -   ������� ������������ � ���������� ����� (��. ����).
                ��������� ���������: [1..Min(M,N)]
                
������� A ��������������, ��� A = P * L * U, ��� P - ������� ������������,
������� L - ���������������� (��� ��������������������, ���� M>N) �������,
U - ����������������� (��� ���������������������, ���� M<N) �������.

���������� ���������� ����� �������� �� ������� ��� M=4, N=3:

                   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
                   ( L31 L32  1  )    (         U33  )
                   ( L41 L42 L43 )
                   
����� ������� L  �����  ������  M  x  Min(M,N),  �������  U  �����  ������
Min(M,N) x N, �������  P(i)  ����������  �����  ������������  �  ���������
������� �������� M x M ����� � �������� I � Pivots[I]

����������� ������ ��������� �������� ������ Pivots  �  ��������� �������,
����������  �������  A,  �  �����������  � ���������� ����� ������� L � U
(������ �������� ��� M=4, N=3):

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

��� �����, ��������� ��������� ������� L  ��  �����������.
���� N>M, �� �������������� �������� ������� ������ � ������������
���������.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
void ludecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    int i;
    int j;
    int jp;
    ap::real_1d_array t1;
    double s;

    pivots.setbounds(1, ap::minint(m, n));
    t1.setbounds(1, ap::maxint(m, n));
    ap::ap_error::make_assertion(m>=0&&n>=0);
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 1; j <= ap::minint(m, n); j++)
    {
        jp = j;
        for(i = j+1; i <= m; i++)
        {
            if( fabs(a(i,j))>fabs(a(jp,j)) )
            {
                jp = i;
            }
        }
        pivots(j) = jp;
        if( a(jp,j)!=0 )
        {
            if( jp!=j )
            {
                ap::vmove(t1.getvector(1, n), a.getrow(j, 1, n));
                ap::vmove(a.getrow(j, 1, n), a.getrow(jp, 1, n));
                ap::vmove(a.getrow(jp, 1, n), t1.getvector(1, n));
            }
            if( j<m )
            {
                jp = j+1;
                s = double(1)/double(a(j,j));
                ap::vmul(a.getcolumn(j, jp, m), s);
            }
        }
        if( j<ap::minint(m, n) )
        {
            jp = j+1;
            for(i = j+1; i <= m; i++)
            {
                s = a(i,j);
                ap::vsub(a.getrow(i, jp, n), a.getrow(j, jp, n), s);
            }
        }
    }
}


/*************************************************************************
LU-���������� ������� ������ ���� ������� M x N

����������  LUDecomposition.   ��  ����������������  ����������  ���,  ���
�������  �������  L  �  U �� � ���������� �����, � � ���� ��������� ������
������ ����, ����������� � ��������������� ������ �������� ����������.

������������ ��������� ������������� ��� ������������ ����, ���
"���������������" ��������� ������ ������������ LUDecomposition

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void ludecompositionunpacked(ap::real_2d_array a,
     int m,
     int n,
     ap::real_2d_array& l,
     ap::real_2d_array& u,
     ap::integer_1d_array& pivots)
{
    int i;
    int j;
    int minmn;

    if( m==0||n==0 )
    {
        return;
    }
    minmn = ap::minint(m, n);
    l.setbounds(1, m, 1, minmn);
    u.setbounds(1, minmn, 1, n);
    ludecomposition(a, m, n, pivots);
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= minmn; j++)
        {
            if( j>i )
            {
                l(i,j) = 0;
            }
            if( j==i )
            {
                l(i,j) = 1;
            }
            if( j<i )
            {
                l(i,j) = a(i,j);
            }
        }
    }
    for(i = 1; i <= minmn; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if( j<i )
            {
                u(i,j) = 0;
            }
            if( j>=i )
            {
                u(i,j) = a(i,j);
            }
        }
    }
}


