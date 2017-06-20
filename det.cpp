/************************************************
Этот модуль сгенерирован транслятором AlgoPascal.
************************************************/

#include "ap.h"
#include "det.h"

/*************************************************************************
Вычисление определителя матрицы, заданной LU-разложением.

Входные параметры:
    A       -   LU-разложение  матрицы   (результат   работы  подпрограммы
                LUDecomposition).
    Pivots  -   таблица перестановок,  произведенных в ходе LU-разложения.
                (результат работы подпрограммы LUDecomposition).
    N       -   размерность матрицы
    
Результат: определитель матрицы

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
double determinantlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n)
{
    double result;
    int i;
    int s;

    result = 1;
    s = 1;
    for(i = 1; i <= n; i++)
    {
//        result = result*a(i,i);
        result = result*(a(i,i)/100.);
        if( pivots(i)!=i )
        {
            s = -s;
        }
    }
    result = result*s;
    return result;
}


/*************************************************************************
Вычисление определителя матрицы общего вида

Входные параметры:
    A   -   матрица. массив с нумерацией элементов [1..N, 1..N]
    N   -   размерность матрицы A

Результат: определитель матрицы A

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
double determinant(ap::real_2d_array a, int n)
{
    double result;
    ap::integer_1d_array pivots;

    ludecomposition(a, n, n, pivots);
    result = determinantlu(a, pivots, n);
    return result;
}






//lu.cpp


/*************************************************************************
LU-разложение матрицы общего вида размера M x N

Подпрограмма вычисляет LU-разложение прямоугольной матрицы общего  вида  с
частичным выбором ведущего элемента (с перестановками строк).

Входные параметры:
    A       -   матрица A. Нумерация элементов: [1..M, 1..N]
    M       -   число строк в матрице A
    N       -   число столбцов в матрице A

Выходные параметры:
    A       -   матрицы L и U в компактной форме (см. ниже).
                Нумерация элементов: [1..M, 1..N]
    Pivots  -   матрица перестановок в компактной форме (см. ниже).
                Нумерация элементов: [1..Min(M,N)]
                
Матрица A представляется, как A = P * L * U, где P - матрица перестановок,
матрица L - нижнетреугольная (или нижнетрапецоидальная, если M>N) матрица,
U - верхнетреугольная (или верхнетрапецоидальная, если M<N) матрица.

Рассмотрим разложение более подробно на примере при M=4, N=3:

                   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
                   ( L31 L32  1  )    (         U33  )
                   ( L41 L42 L43 )

Здесь матрица L  имеет  размер  M  x  Min(M,N),  матрица  U  имеет  размер
Min(M,N) x N, матрица  P(i)  получается  путем  перестановки  в  единичной
матрице размером M x M строк с номерами I и Pivots[I]

Результатом работы алгоритма являются массив Pivots  и  следующая матрица,
замещающая  матрицу  A,  и  сохраняющая  в компактной форме матрицы L и U
(пример приведен для M=4, N=3):

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

Как видно, единичная диагональ матрицы L  не  сохраняется.
Если N>M, то соответственно меняются размеры матриц и расположение
элементов.

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
LU-разложение матрицы общего вида размера M x N

Использует  LUDecomposition.   По  функциональности  отличается  тем,  что
выводит  матрицы  L  и  U не в компактной форме, а в виде отдельных матриц
общего вида, заполненных в соответствующих местах нулевыми элементами.

Подпрограмма приведена исключительно для демонстрации того, как
"распаковывается" результат работы подпрограммы LUDecomposition

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


void safesolvetriangular(const ap::real_2d_array& a,
     int n,
     ap::real_1d_array& x,
     double& s,
     bool isupper,
     bool istrans,
     bool isunit,
     bool normin,
     ap::real_1d_array& cnorm)
{
    int i;
    int imax;
    int j;
    int jfirst;
    int jinc;
    int jlast;
    int jm1;
    int jp1;
    int ip1;
    int im1;
    int k;
    int flg;
    double v;
    double vd;
    double bignum;
    double grow;
    double rec;
    double smlnum;
    double sumj;
    double tjj;
    double tjjs;
    double tmax;
    double tscal;
    double uscal;
    double xbnd;
    double xj;
    double xmax;
    bool notran;
    bool upper;
    bool nounit;

    upper = isupper;
    notran = !istrans;
    nounit = !isunit;
    if( n==0 )
    {
        return;
    }
    smlnum = ap::minrealnumber/(ap::machineepsilon*2);
    bignum = double(1)/double(smlnum);
    s = 1;
    if( !normin )
    {
        cnorm.setbounds(1, n);
        if( upper )
        {
            for(j = 1; j <= n; j++)
            {
                v = 0;
                for(k = 1; k <= j-1; k++)
                {
                    v = v+fabs(a(k,j));
                }
                cnorm(j) = v;
            }
        }
        else
        {
            for(j = 1; j <= n-1; j++)
            {
                v = 0;
                for(k = j+1; k <= n; k++)
                {
                    v = v+fabs(a(k,j));
                }
                cnorm(j) = v;
            }
            cnorm(n) = 0;
        }
    }
    imax = 1;
    for(k = 2; k <= n; k++)
    {
        if( cnorm(k)>cnorm(imax) )
        {
            imax = k;
        }
    }
    tmax = cnorm(imax);
    if( tmax<=bignum )
    {
        tscal = 1;
    }
    else
    {
        tscal = double(1)/double(smlnum*tmax);
        ap::vmul(cnorm.getvector(1, n), tscal);
    }
    j = 1;
    for(k = 2; k <= n; k++)
    {
        if( fabs(x(k))>fabs(x(j)) )
        {
            j = k;
        }
    }
    xmax = fabs(x(j));
    xbnd = xmax;
    if( notran )
    {
        if( upper )
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        else
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        if( tscal!=1 )
        {
            grow = 0;
        }
        else
        {
            if( nounit )
            {
                grow = double(1)/double(ap::maxreal(xbnd, smlnum));
                xbnd = grow;
                j = jfirst;
                while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
                {
                    if( grow<=smlnum )
                    {
                        break;
                    }
                    tjj = fabs(a(j,j));
                    xbnd = ap::minreal(xbnd, ap::minreal(1, tjj)*grow);
                    if( tjj+cnorm(j)>=smlnum )
                    {
                        grow = grow*(tjj/(tjj+cnorm(j)));
                    }
                    else
                    {
                        grow = 0;
                    }
                    if( j==jlast )
                    {
                        grow = xbnd;
                    }
                    j = j+jinc;
                }
            }
            else
            {
                grow = ap::minreal(1, double(1)/double(ap::maxreal(xbnd, smlnum)));
                j = jfirst;
                while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
                {
                    if( grow<=smlnum )
                    {
                        break;
                    }
                    grow = grow*(double(1)/double(1+cnorm(j)));
                    j = j+jinc;
                }
            }
        }
    }
    else
    {
        if( upper )
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        else
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        if( tscal!=1 )
        {
            grow = 0;
        }
        else
        {
            if( nounit )
            {
                grow = double(1)/double(ap::maxreal(xbnd, smlnum));
                xbnd = grow;
                j = jfirst;
                while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
                {
                    if( grow<=smlnum )
                    {
                        break;
                    }
                    xj = 1+cnorm(j);
                    grow = ap::minreal(grow, xbnd/xj);
                    tjj = fabs(a(j,j));
                    if( xj>tjj )
                    {
                        xbnd = xbnd*(tjj/xj);
                    }
                    if( j==jlast )
                    {
                        grow = ap::minreal(grow, xbnd);
                    }
                    j = j+jinc;
                }
            }
            else
            {
                grow = ap::minreal(1, double(1)/double(ap::maxreal(xbnd, smlnum)));
                j = jfirst;
                while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
                {
                    if( grow<=smlnum )
                    {
                        break;
                    }
                    xj = 1+cnorm(j);
                    grow = grow/xj;
                    j = j+jinc;
                }
            }
        }
    }
    if( grow*tscal>smlnum )
    {
        if( upper&&notran||!upper&&!notran )
        {
            if( nounit )
            {
                vd = a(n,n);
            }
            else
            {
                vd = 1;
            }
            x(n) = x(n)/vd;
            for(i = n-1; i >= 1; i--)
            {
                ip1 = i+1;
                if( upper )
                {
                    v = ap::vdotproduct(a.getrow(i, ip1, n), x.getvector(ip1, n));
                }
                else
                {
                    v = ap::vdotproduct(a.getcolumn(i, ip1, n), x.getvector(ip1, n));
                }
                if( nounit )
                {
                    vd = a(i,i);
                }
                else
                {
                    vd = 1;
                }
                x(i) = (x(i)-v)/vd;
            }
        }
        else
        {
            if( nounit )
            {
                vd = a(1,1);
            }
            else
            {
                vd = 1;
            }
            x(1) = x(1)/vd;
            for(i = 2; i <= n; i++)
            {
                im1 = i-1;
                if( upper )
                {
                    v = ap::vdotproduct(a.getcolumn(i, 1, im1), x.getvector(1, im1));
                }
                else
                {
                    v = ap::vdotproduct(a.getrow(i, 1, im1), x.getvector(1, im1));
                }
                if( nounit )
                {
                    vd = a(i,i);
                }
                else
                {
                    vd = 1;
                }
                x(i) = (x(i)-v)/vd;
            }
        }
    }
    else
    {
        if( xmax>bignum )
        {
            s = bignum/xmax;
            ap::vmul(x.getvector(1, n), s);
            xmax = bignum;
        }
        if( notran )
        {
            j = jfirst;
            while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
            {
                xj = fabs(x(j));
                flg = 0;
                if( nounit )
                {
                    tjjs = a(j,j)*tscal;
                }
                else
                {
                    tjjs = tscal;
                    if( tscal==1 )
                    {
                        flg = 100;
                    }
                }
                if( flg!=100 )
                {
                    tjj = fabs(tjjs);
                    if( tjj>smlnum )
                    {
                        if( tjj<1 )
                        {
                            if( xj>tjj*bignum )
                            {
                                rec = double(1)/double(xj);
                                ap::vmul(x.getvector(1, n), rec);
                                s = s*rec;
                                xmax = xmax*rec;
                            }
                        }
                        x(j) = x(j)/tjjs;
                        xj = fabs(x(j));
                    }
                    else
                    {
                        if( tjj>0 )
                        {
                            if( xj>tjj*bignum )
                            {
                                rec = tjj*bignum/xj;
                                if( cnorm(j)>1 )
                                {
                                    rec = rec/cnorm(j);
                                }
                                ap::vmul(x.getvector(1, n), rec);
                                s = s*rec;
                                xmax = xmax*rec;
                            }
                            x(j) = x(j)/tjjs;
                            xj = fabs(x(j));
                        }
                        else
                        {
                            for(i = 1; i <= n; i++)
                            {
                                x(i) = 0;
                            }
                            x(j) = 1;
                            xj = 1;
                            s = 0;
                            xmax = 0;
                        }
                    }
                }
                if( xj>1 )
                {
                    rec = double(1)/double(xj);
                    if( cnorm(j)>(bignum-xmax)*rec )
                    {
                        rec = rec*0.5;
                        ap::vmul(x.getvector(1, n), rec);
                        s = s*rec;
                    }
                }
                else
                {
                    if( xj*cnorm(j)>bignum-xmax )
                    {
                        ap::vmul(x.getvector(1, n), 0.5);
                        s = s*0.5;
                    }
                }
                if( upper )
                {
                    if( j>1 )
                    {
                        v = x(j)*tscal;
                        jm1 = j-1;
                        ap::vsub(x.getvector(1, jm1), a.getcolumn(j, 1, jm1), v);
                        i = 1;
                        for(k = 2; k <= j-1; k++)
                        {
                            if( fabs(x(k))>fabs(x(i)) )
                            {
                                i = k;
                            }
                        }
                        xmax = fabs(x(i));
                    }
                }
                else
                {
                    if( j<n )
                    {
                        jp1 = j+1;
                        v = x(j)*tscal;
                        ap::vsub(x.getvector(jp1, n), a.getcolumn(j, jp1, n), v);
                        i = j+1;
                        for(k = j+2; k <= n; k++)
                        {
                            if( fabs(x(k))>fabs(x(i)) )
                            {
                                i = k;
                            }
                        }
                        xmax = fabs(x(i));
                    }
                }
                j = j+jinc;
            }
        }
        else
        {
            j = jfirst;
            while(jinc>0&&j<=jlast||jinc<0&&j>=jlast)
            {
                xj = fabs(x(j));
                uscal = tscal;
                rec = double(1)/double(ap::maxreal(xmax, 1));
                if( cnorm(j)>(bignum-xj)*rec )
                {
                    rec = rec*0.5;
                    if( nounit )
                    {
                        tjjs = a(j,j)*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                    }
                    tjj = fabs(tjjs);
                    if( tjj>1 )
                    {
                        rec = ap::minreal(1, rec*tjj);
                        uscal = uscal/tjjs;
                    }
                    if( rec<1 )
                    {
                        ap::vmul(x.getvector(1, n), rec);
                        s = s*rec;
                        xmax = xmax*rec;
                    }
                }
                sumj = 0;
                if( uscal==1 )
                {
                    if( upper )
                    {
                        if( j>1 )
                        {
                            jm1 = j-1;
                            sumj = ap::vdotproduct(a.getcolumn(j, 1, jm1), x.getvector(1, jm1));
                        }
                        else
                        {
                            sumj = 0;
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            jp1 = j+1;
                            sumj = ap::vdotproduct(a.getcolumn(j, jp1, n), x.getvector(jp1, n));
                        }
                    }
                }
                else
                {
                    if( upper )
                    {
                        for(i = 1; i <= j-1; i++)
                        {
                            v = a(i,j)*uscal;
                            sumj = sumj+v*x(i);
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            for(i = j+1; i <= n; i++)
                            {
                                v = a(i,j)*uscal;
                                sumj = sumj+v*x(i);
                            }
                        }
                    }
                }
                if( uscal==tscal )
                {
                    x(j) = x(j)-sumj;
                    xj = fabs(x(j));
                    flg = 0;
                    if( nounit )
                    {
                        tjjs = a(j,j)*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                        if( tscal==1 )
                        {
                            flg = 150;
                        }
                    }
                    if( flg!=150 )
                    {
                        tjj = fabs(tjjs);
                        if( tjj>smlnum )
                        {
                            if( tjj<1 )
                            {
                                if( xj>tjj*bignum )
                                {
                                    rec = double(1)/double(xj);
                                    ap::vmul(x.getvector(1, n), rec);
                                    s = s*rec;
                                    xmax = xmax*rec;
                                }
                            }
                            x(j) = x(j)/tjjs;
                        }
                        else
                        {
                            if( tjj>0 )
                            {
                                if( xj>tjj*bignum )
                                {
                                    rec = tjj*bignum/xj;
                                    ap::vmul(x.getvector(1, n), rec);
                                    s = s*rec;
                                    xmax = xmax*rec;
                                }
                                x(j) = x(j)/tjjs;
                            }
                            else
                            {
                                for(i = 1; i <= n; i++)
                                {
                                    x(i) = 0;
                                }
                                x(j) = 1;
                                s = 0;
                                xmax = 0;
                            }
                        }
                    }
                }
                else
                {
                    x(j) = x(j)/tjjs-sumj;
                }
                xmax = ap::maxreal(xmax, fabs(x(j)));
                j = j+jinc;
            }
        }
        s = s/tscal;
    }
    if( tscal!=1 )
    {
        v = double(1)/double(tscal);
        ap::vmul(cnorm.getvector(1, n), v);
    }
}

/*************************************************************************
Оценка числа обусловленности матрицы (1-норма).

Алгоритм вычисляет нижнюю границу числа обусловленности матрицы. При  этом
возвращается не  сама  нижняя  граница числа  обусловленности,  а  единица,
деленная на оценку числа  обусловленности  (чтобы  избежать  переполнения,
если матрица вырождена).

Входные параметры:
    A   -   матрица. Массив с нумерацией элементов [1..N, 1..N]
    N   -   размерность матрицы.
    
Результат: 1/LowerBound(cond(A))
*************************************************************************/
double rcond1(ap::real_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;

    nrm = 0;
    for(j = 1; j <= n; j++)
    {
        v = 0;
        for(i = 1; i <= n; i++)
        {
            v = v+fabs(a(i,j));
        }
        nrm = ap::maxreal(nrm, v);
    }
    ludecomposition(a, n, n, pivots);
    internalestimatercondlu(a, n, true, true, nrm, v);
    result = v;
    return result;
}


/*************************************************************************
Оценка числа обусловленности матрицы, заданной LU-разложением (1-норма).

Алгоритм вычисляет нижнюю границу числа обусловленности матрицы. При  этом
возвращается не  сама  нижняя  граница числа  обусловленности,  а  единица,
деленная на оценку числа  обусловленности  (чтобы  избежать  переполнения,
если матрица вырождена).

Входные параметры:
    LU  -   LU-разложение матрицы в упакованной форме. Результат работы
            подпрограммы LUDecomposition
    N   -   размерность матрицы.

Результат: 1/LowerBound(cond(A))
*************************************************************************/
double rcond1lu(const ap::real_2d_array& lu, int n)
{
    double result;
    double nrm;
    double v;

    internalestimatercondlu(lu, n, true, false, nrm, v);
    result = v;
    return result;
}


/*************************************************************************
Оценка числа обусловленности матрицы (бесконечная норма).

Алгоритм вычисляет нижнюю границу числа обусловленности матрицы. При  этом
возвращается не  сама  нижняя  граница числа  обусловленности,  а  единица,
деленная на оценку числа  обусловленности  (чтобы  избежать  переполнения,
если матрица вырождена).

Входные параметры:
    A   -   матрица. Массив с нумерацией элементов [1..N, 1..N]
    N   -   размерность матрицы.

Результат: 1/LowerBound(cond(A))
*************************************************************************/
double rcondinf(ap::real_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;

    nrm = 0;
    for(i = 1; i <= n; i++)
    {
        v = 0;
        for(j = 1; j <= n; j++)
        {
            v = v+fabs(a(i,j));
        }
        nrm = ap::maxreal(nrm, v);
    }
    ludecomposition(a, n, n, pivots);
    internalestimatercondlu(a, n, false, true, nrm, v);
    result = v;
    return result;
}


/*************************************************************************
Оценка числа обусловленности матрицы, заданной LU-разложением
(бесконечная норма).

Алгоритм вычисляет нижнюю границу числа обусловленности матрицы. При  этом
возвращается не  сама  нижняя  граница числа  обусловленности,  а  единица,
деленная на оценку числа  обусловленности  (чтобы  избежать  переполнения,
если матрица вырождена).

Входные параметры:
    LU  -   LU-разложение матрицы в упакованной форме. Результат работы
            подпрограммы LUDecomposition
    N   -   размерность матрицы.

Результат: 1/LowerBound(cond(A))
*************************************************************************/
double rcondinflu(const ap::real_2d_array& lu, int n)
{
    double result;
    double nrm;
    double v;

    internalestimatercondlu(lu, n, false, false, nrm, v);
    result = v;
    return result;
}


/*************************************************************************
Оценка числа обусловленности матрицы.

Служебная подпрограмма. Вызывается функциями RCond1,  RCond1LU,  RCondInf,
RCondInfLU.

Алгоритм вычисляет нижнюю границу числа обусловленности матрицы,  заданной
LU-разожением. В зависимости от переданных параметров вычисляется  1-норма
или бесконечная норма. При этом возвращается не само число обусловленности
а единица, деленная на это число (чтобы избежать переполнения, если матрица
вырождена).

Входные параметры:
    LU  -       LU-разложение матрицы, норма которой вычисляется. Результат
                работы подпрограммы LUDecomposition.
    N   -       размерность матрицы.
    OneNorm-    флаг, задающий тип нормы (1-норма или бесконечна норма).
    IsANormProvided-флаг,  сообщающий, располагает ли алгоритм информацией
                о точной норме матрицы A. Если параметр равен True,  то  в
                параметре ANORM находится точное значение нормы матрицы A,
                LU-разложение которой нам передается. В  противном  случае
                алгоритм самостоятельно оценивает норму исходной матрицы.
    ANORM-      норма матрицы A, если IsANormProvided = True.
    
Выходные параметры:
    RCOND-      1/(LowerBound(cond(A)))

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
void internalestimatercondlu(const ap::real_2d_array& lu,
     int n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rcond)
{
    ap::real_1d_array work0;
    ap::real_1d_array work1;
    ap::real_1d_array work2;
    ap::real_1d_array work3;
    ap::integer_1d_array iwork;
    double v;
    bool normin;
    int i;
    int j;
    int im1;
    int ip1;
    int ix;
    int kase;
    int kase1;
    double ainvnm;
    double ascale;
    double sl;
    double smlnum;
    double su;
    bool mupper;
    bool mtrans;
    bool munit;

    if( n==0 )
    {
        rcond = 1;
        return;
    }
    if( onenorm )
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    mupper = true;
    mtrans = true;
    munit = true;
    work0.setbounds(1, n);
    work1.setbounds(1, n);
    work2.setbounds(1, n);
    work3.setbounds(1, n);
    iwork.setbounds(1, n);
    if( !isanormprovided )
    {
        kase = 0;
        anorm = 0;
        while(true)
        {
            internalestimatenorm(n, work1, work0, iwork, anorm, kase);
            if( kase==0 )
            {
                break;
            }
            if( kase==kase1 )
            {
                for(i = 1; i <= n; i++)
                {
                    v = ap::vdotproduct(lu.getrow(i, i, n), work0.getvector(i, n));
                    work0(i) = v;
                }
                for(i = n; i >= 1; i--)
                {
                    im1 = i-1;
                    if( i>1 )
                    {
                        v = ap::vdotproduct(lu.getrow(i, 1, im1), work0.getvector(1, im1));
                    }
                    else
                    {
                        v = 0;
                    }
                    work0(i) = work0(i)+v;
                }
            }
            else
            {
                for(i = 1; i <= n; i++)
                {
                    ip1 = i+1;
                    v = ap::vdotproduct(lu.getcolumn(i, ip1, n), work0.getvector(ip1, n));
                    work0(i) = work0(i)+v;
                }
                for(i = n; i >= 1; i--)
                {
                    v = ap::vdotproduct(lu.getcolumn(i, 1, i), work0.getvector(1, i));
                    work0(i) = v;
                }
            }
        }
    }
    rcond = 0;
    if( anorm==0 )
    {
        return;
    }
    smlnum = ap::minrealnumber;
    ainvnm = 0;
    normin = false;
    kase = 0;
    while(true)
    {
        internalestimatenorm(n, work1, work0, iwork, ainvnm, kase);
        if( kase==0 )
        {
            break;
        }
        if( kase==kase1 )
        {
            safesolvetriangular(lu, n, work0, sl, !mupper, !mtrans, munit, normin, work2);
            safesolvetriangular(lu, n, work0, su, mupper, !mtrans, !munit, normin, work3);
        }
        else
        {
            safesolvetriangular(lu, n, work0, su, mupper, mtrans, !munit, normin, work3);
            safesolvetriangular(lu, n, work0, sl, !mupper, mtrans, munit, normin, work2);
        }
        ascale = sl*su;
        normin = true;
        if( ascale!=1 )
        {
            ix = 1;
            for(i = 2; i <= n; i++)
            {
                if( fabs(work0(i))>fabs(work0(ix)) )
                {
                    ix = i;
                }
            }
            if( ascale<fabs(work0(ix))*smlnum||ascale==0 )
            {
                return;
            }
            for(i = 1; i <= n; i++)
            {
                work0(i) = work0(i)/ascale;
            }
        }
    }
    if( ainvnm!=0 )
    {
        rcond = double(1)/double(ainvnm);
        rcond = rcond/anorm;
    }
}


/*************************************************************************

Служебная подпрограмма, проводящая оценку нормы матрицы.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
void internalestimatenorm(int n,
     ap::real_1d_array& v,
     ap::real_1d_array& x,
     ap::integer_1d_array& isgn,
     double& est,
     int& kase)
{
    int itmax;
    int i;
    double t;
    bool flg;
    int positer;
    int posj;
    int posjlast;
    int posjump;
    int posaltsgn;
    int posestold;
    int postemp;

    itmax = 5;
    posaltsgn = n+1;
    posestold = n+2;
    postemp = n+3;
    positer = n+1;
    posj = n+2;
    posjlast = n+3;
    posjump = n+4;
    if( kase==0 )
    {
        v.setbounds(1, n+3);
        x.setbounds(1, n);
        isgn.setbounds(1, n+4);
        t = double(1)/double(n);
        for(i = 1; i <= n; i++)
        {
            x(i) = t;
        }
        kase = 1;
        isgn(posjump) = 1;
        return;
    }
    if( isgn(posjump)==1 )
    {
        if( n==1 )
        {
            v(1) = x(1);
            est = fabs(v(1));
            kase = 0;
            return;
        }
        est = 0;
        for(i = 1; i <= n; i++)
        {
            est = est+fabs(x(i));
        }
        for(i = 1; i <= n; i++)
        {
            if( x(i)>=0 )
            {
                x(i) = 1;
            }
            else
            {
                x(i) = -1;
            }
            isgn(i) = ap::sign(x(i));
        }
        kase = 2;
        isgn(posjump) = 2;
        return;
    }
    if( isgn(posjump)==2 )
    {
        isgn(posj) = 1;
        for(i = 2; i <= n; i++)
        {
            if( fabs(x(i))>fabs(x(isgn(posj))) )
            {
                isgn(posj) = i;
            }
        }
        isgn(positer) = 2;
        for(i = 1; i <= n; i++)
        {
            x(i) = 0;
        }
        x(isgn(posj)) = 1;
        kase = 1;
        isgn(posjump) = 3;
        return;
    }
    if( isgn(posjump)==3 )
    {
        ap::vmove(v.getvector(1, n), x.getvector(1, n));
        v(posestold) = est;
        est = 0;
        for(i = 1; i <= n; i++)
        {
            est = est+fabs(v(i));
        }
        flg = false;
        for(i = 1; i <= n; i++)
        {
            if( x(i)>=0&&isgn(i)<0||x(i)<0&&isgn(i)>=0 )
            {
                flg = true;
            }
        }
        if( !flg||est<=v(posestold) )
        {
            v(posaltsgn) = 1;
            for(i = 1; i <= n; i++)
            {
                x(i) = v(posaltsgn)*(1+double(i-1)/double(n-1));
                v(posaltsgn) = -v(posaltsgn);
            }
            kase = 1;
            isgn(posjump) = 5;
            return;
        }
        for(i = 1; i <= n; i++)
        {
            if( x(i)>=0 )
            {
                x(i) = 1;
                isgn(i) = 1;
            }
            else
            {
                x(i) = -1;
                isgn(i) = -1;
            }
        }
        kase = 2;
        isgn(posjump) = 4;
        return;
    }
    if( isgn(posjump)==4 )
    {
        isgn(posjlast) = isgn(posj);
        isgn(posj) = 1;
        for(i = 2; i <= n; i++)
        {
            if( fabs(x(i))>fabs(x(isgn(posj))) )
            {
                isgn(posj) = i;
            }
        }
        if( x(isgn(posjlast))!=fabs(x(isgn(posj)))&&isgn(positer)<itmax )
        {
            isgn(positer) = isgn(positer)+1;
            for(i = 1; i <= n; i++)
            {
                x(i) = 0;
            }
            x(isgn(posj)) = 1;
            kase = 1;
            isgn(posjump) = 3;
            return;
        }
        v(posaltsgn) = 1;
        for(i = 1; i <= n; i++)
        {
            x(i) = v(posaltsgn)*(1+double(i-1)/double(n-1));
            v(posaltsgn) = -v(posaltsgn);
        }
        kase = 1;
        isgn(posjump) = 5;
        return;
    }
    if( isgn(posjump)==5 )
    {
        v(postemp) = 0;
        for(i = 1; i <= n; i++)
        {
            v(postemp) = v(postemp)+fabs(x(i));
        }
        v(postemp) = double(2*v(postemp))/double(3*n);
        if( v(postemp)>est )
        {
            ap::vmove(v.getvector(1, n), x.getvector(1, n));
            est = v(postemp);
        }
        kase = 0;
        return;
    }
}


bool solvesystemlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x)
{
    bool result;
    ap::real_1d_array y;
    int i;
    int j;
    double v;
    int ip1;
    int im1;

    y.setbounds(1, n);
    x.setbounds(1, n);
    result = true;
    for(i = 1; i <= n; i++)
    {
        if( a(i,i)==0 )
        {
            result = false;
            return result;
        }
    }
    for(i = 1; i <= n; i++)
    {
        if( pivots(i)!=i )
        {
            v = b(i);
            b(i) = b(pivots(i));
            b(pivots(i)) = v;
        }
    }
    y(1) = b(1);
    for(i = 2; i <= n; i++)
    {
        im1 = i-1;
        v = ap::vdotproduct(a.getrow(i, 1, im1), y.getvector(1, im1));
        y(i) = b(i)-v;
    }
    x(n) = y(n)/a(n,n);
    for(i = n-1; i >= 1; i--)
    {
        ip1 = i+1;
        v = ap::vdotproduct(a.getrow(i, ip1, n), x.getvector(ip1, n));
        x(i) = (y(i)-v)/a(i,i);
    }
    return result;
}


/*************************************************************************
Решение системы  линейных  уравнений

Алгоритм решает систему линейных уравнений с использованием LU-разложения.
Алгоритм решает только системы уравнений с квадратной матрицей.

Входные параметры:
    A   -   Матрица системы.
            Массив с нумерацией элементов [1..N, 1..N].
    B   -   Правая часть.
            Массив с нумерацией элементов [1..N, 1..N].
    N   -   размерность системы.

Выходные параметры:
    X   -   решение системы. Массив с нумерацией элементов [1..N]

Результат:
    True, если система не вырождена (но, возможно, близка к вырожденной).
    False, если система вырождена. В таком случае X не содержит решение.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
bool solvesystem(ap::real_2d_array a,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x)
{
    bool result;
    ap::integer_1d_array pivots;
    int i;

    ludecomposition(a, n, n, pivots);
    result = solvesystemlu(a, pivots, b, n, x);
    return result;
}


