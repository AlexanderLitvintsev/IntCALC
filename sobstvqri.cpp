/************************************************
Этот модуль сгенерирован транслятором AlgoPascal.
************************************************/

#include "ap.h"
#include "sobstvqri.h"


/*************************************************************************
Получение собственных пар произвольной вещественной матрицы
приведением к матрице Хессенберга и использованием QR алгоритма
с неявными сдвигами.

Входные данные:
    a - массив с нумерацией элементов [1..n, 1..n]. Содержит элементы
        матрицы.
    n - размерность задачи

Выходные данные:
    EVr - Массив реальных частей собственных значений.
          Нумерация элементов [1..n].
    EVi - Массив мнимых   частей собственных значений.
          Нумерация элементов [1..n].
    V   -     Матрица, столбцы которой являются собственными векторами.
          Нумерация элементов [1..n, 1..n].
              Если J-ый элемент массива EVi равен 0, то мы имеем
          вещественное собственное число EVr[J], которому соответствует
          собственный вектор - J-ый столбец матрицы V.
              Если J-ый элемент массива EVi больше 0, то комплексному
          собственному числу соответствует вектор V(j) + i*V(j+1), где
          V(j), V(j+1) - столбцы матрицы V.
              Если J-ый элемент массива EVi меньше 0, то комплексному
          собственному числу соответствует вектор V(j-1) - i*V(j), где
          V(j), V(j-1) - столбцы матрицы V.

Результат:
    True,  если алгоритм корректно завершил работу.
    False, если в течение 30 итераций не была достигнута сходимость.

*************************************************************************/
bool hessenbergqrieigenvaluesandvectors(ap::real_2d_array a,
     int n,
     ap::real_1d_array& evr,
     ap::real_1d_array& evi,
     ap::real_2d_array& v)
{
    bool result;
    double r1;
    double r2;
    double r3;
    double q1r;
    double q2r;
    double q3r;
    double q1i;
    double q2i;
    double q3i;
    ap::integer_1d_array irab1;
    int ierr;
    int b;
    int c;
    int d;
    int e;
    int f;
    int g;
    int h;
    int i1;
    int i2;
    int i3;
    int k;
    int l;
    int m;
    int o;
    int p;
    int q;
    int r;
    int s;
    int t;
    int u;
    int w;
    int x;
    int y;
    int z;
    int ba;
    int bb;
    int bc;
    int bd;
    int be;
    int bf;
    int bg;
    int ir1;
    int ir2;
    int flagjump19;
    bool bw;
    double i;
    double j;
    double bh;
    double bi;
    double bj;
    double bk;
    double bl;
    double bm;
    double bn;
    double bo;
    double bp;
    double bq;
    double br;
    double bs;
    double bt;
    double bu;

    evr.setbounds(1, n);
    evi.setbounds(1, n);
    v.setbounds(1, n, 1, n);
    irab1.setbounds(1, n);
    ir1 = 1;
    ir2 = n;
    e = ir2-1;
    f = ir1+1;
    if( e>=f )
    {
        i1 = e;
        for(d = f; d <= i1; d++)
        {
            g = d-1;
            i = 0;
            b = d;
            i2 = ir2;
            for(c = d; c <= i2; c++)
            {
                if( fabs(a(c,g))<=fabs(i) )
                {
                    continue;
                }
                i = a(c,g);
                b = c;
            }
            irab1(d) = b;
            if( b!=d )
            {
                i2 = n;
                for(c = g; c <= i2; c++)
                {
                    j = a(b,c);
                    a(b,c) = a(d,c);
                    a(d,c) = j;
                }
                i2 = ir2;
                for(c = 1; c <= i2; c++)
                {
                    j = a(c,b);
                    a(c,b) = a(c,d);
                    a(c,d) = j;
                }
            }
            if( i==0 )
            {
                continue;
            }
            h = d+1;
            i2 = ir2;
            for(b = h; b <= i2; b++)
            {
                j = a(b,g);
                if( j==0 )
                {
                    continue;
                }
                j = j/i;
                a(b,g) = j;
                i3 = n;
                for(c = d; c <= i3; c++)
                {
                    a(b,c) = a(b,c)-j*a(d,c);
                }
                i3 = ir2;
                for(c = 1; c <= i3; c++)
                {
                    a(c,d) = a(c,d)+j*a(c,b);
                }
            }
        }
    }
    i1 = n;
    for(k = 1; k <= i1; k++)
    {
        i2 = n;
        for(l = 1; l <= i2; l++)
        {
            v(k,l) = 0;
        }
        v(k,k) = 1;
    }
    m = ir2-ir1-1;
    if( m>=1 )
    {
        i1 = m;
        for(o = 1; o <= i1; o++)
        {
            p = ir2-o;
            q = p+1;
            i2 = ir2;
            for(k = q; k <= i2; k++)
            {
                v(k,p) = a(k,p-1);
            }
            k = irab1(p);
            if( k==p )
            {
                continue;
            }
            i2 = ir2;
            for(l = p; l <= i2; l++)
            {
                v(p,l) = v(k,l);
                v(k,l) = 0;
            }
            v(k,p) = 1;
        }
    }
    ierr = 0;
    bu = 0;
    t = 1;
    i1 = n;
    for(r = 1; r <= i1; r++)
    {
        i2 = n;
        for(s = t; s <= i2; s++)
        {
            bu = bu+fabs(a(r,s));
        }
        t = r;
        if( r>=ir1&&r<=ir2 )
        {
            continue;
        }
        evr(r) = a(r,r);
        evi(r) = 0;
    }
    x = ir2;
    bl = 0;
    flagjump19 = 0;
    while(true)
    {
        if( flagjump19==0 )
        {
            if( x<ir1 )
            {
                if( bu==0 )
                {
                    result = ierr==0;
                    return result;
                }
                i1 = n;
                for(bd = 1; bd <= i1; bd++)
                {
                    x = n+1-bd;
                    bh = evr(x);
                    bi = evi(x);
                    bc = x-1;
                    if( bi>=0 )
                    {
                        if( bi!=0 )
                        {
                            continue;
                        }
                        w = x;
                        a(x,x) = 1;
                        if( bc==0 )
                        {
                            continue;
                        }
                        i2 = bc;
                        for(y = 1; y <= i2; y++)
                        {
                            r = x-y;
                            bm = a(r,r)-bh;
                            bj = a(r,x);
                            if( w<=bc )
                            {
                                i3 = bc;
                                for(s = w; s <= i3; s++)
                                {
                                    bj = bj+a(r,s)*a(s,x);
                                }
                            }
                            if( evi(r)<0 )
                            {
                                bt = bm;
                                bk = bj;
                                continue;
                            }
                            w = r;
                            if( evi(r)==0 )
                            {
                                bl = bm;
                                if( bm==0 )
                                {
                                    bl = 100*ap::machineepsilon*bu;
                                }
                                a(r,x) = -bj/bl;
                                continue;
                            }
                            bn = a(r,r+1);
                            bo = a(r+1,r);
                            bi = (evr(r)-bh)*(evr(r)-bh)+evi(r)*evi(r);
                            bl = (bn*bk-bt*bj)/bi;
                            a(r,x) = bl;
                            if( fabs(bn)<=fabs(bt) )
                            {
                                a(r+1,x) = (-bk-bo*bl)/bt;
                            }
                            else
                            {
                                a(r+1,x) = (-bj-bm*bl)/bn;
                            }
                        }
                        continue;
                    }
                    w = bc;
                    if( fabs(a(x,bc))<=fabs(a(bc,x)) )
                    {
                        r1 = -a(bc,x);
                        q2r = 0;
                        q2i = r1;
                        r2 = a(bc,bc)-bh;
                        q3r = r2;
                        q3i = bi;
                        divcomplex(q1r, q1i, q2r, q2i, q3r, q3i);
                        a(bc,bc) = q1r;
                        a(bc,x) = q1i;
                    }
                    else
                    {
                        a(bc,bc) = bi/a(x,bc);
                        a(bc,x) = -(a(x,x)-bh)/a(x,bc);
                    }
                    a(x,bc) = 0;
                    a(x,x) = 1;
                    bg = bc-1;
                    if( bg==0 )
                    {
                        continue;
                    }
                    i2 = bg;
                    for(y = 1; y <= i2; y++)
                    {
                        r = bc-y;
                        bm = a(r,r)-bh;
                        bp = 0;
                        bq = a(r,x);
                        i3 = bc;
                        for(s = w; s <= i3; s++)
                        {
                            bp = bp+a(r,s)*a(s,bc);
                            bq = bq+a(r,s)*a(s,x);
                        }
                        if( evi(r)<0 )
                        {
                            bt = bm;
                            bj = bp;
                            bk = bq;
                            continue;
                        }
                        w = r;
                        if( evi(r)==0 )
                        {
                            r1 = -bp;
                            r2 = -bq;
                            q2r = r1;
                            q2i = r2;
                            q3r = bm;
                            q3i = bi;
                            divcomplex(q1r, q1i, q2r, q2i, q3r, q3i);
                            a(r,bc) = q1r;
                            a(r,x) = q1i;
                            continue;
                        }
                        bn = a(r,r+1);
                        bo = a(r+1,r);
                        bs = (evr(r)-bh)*(evr(r)-bh)+evi(r)*evi(r)-bi*bi;
                        br = (evr(r)-bh)*2*bi;
                        if( bs==0&&br==0 )
                        {
                            bs = 100*ap::machineepsilon*bu*(fabs(bm)+fabs(bi)+fabs(bn)+fabs(bo)+fabs(bt));
                        }
                        r1 = bn*bj-bt*bp+bi*bq;
                        r2 = bn*bk-bt*bq-bi*bp;
                        q2r = r1;
                        q2i = r2;
                        q3r = bs;
                        q3i = br;
                        divcomplex(q1r, q1i, q2r, q2i, q3r, q3i);
                        a(r,bc) = q1r;
                        a(r,x) = q1i;
                        if( fabs(bn)<=fabs(bt)+fabs(bi) )
                        {
                            r1 = -bj-bo*a(r,bc);
                            r2 = -bk-bo*a(r,x);
                            q2r = r1;
                            q2i = r2;
                            q3r = bt;
                            q3i = bi;
                            divcomplex(q1r, q1i, q2r, q2i, q3r, q3i);
                            a(r+1,bc) = q1r;
                            a(r+1,x) = q1i;
                        }
                        else
                        {
                            a(r+1,bc) = (-bp-bm*a(r,bc)+bi*a(r,x))/bn;
                            a(r+1,x) = (-bq-bm*a(r,x)-bi*a(r,bc))/bn;
                        }
                    }
                }
                i1 = n;
                for(r = 1; r <= i1; r++)
                {
                    if( r>=ir1&&r<=ir2 )
                    {
                        continue;
                    }
                    i2 = n;
                    for(s = r; s <= i2; s++)
                    {
                        v(r,s) = a(r,s);
                    }
                }
                i1 = n;
                for(z = ir1; z <= i1; z++)
                {
                    s = n+ir1-z;
                    if( s<ir2 )
                    {
                        w = s;
                    }
                    else
                    {
                        w = ir2;
                    }
                    i2 = ir2;
                    for(r = ir1; r <= i2; r++)
                    {
                        bt = 0;
                        i3 = w;
                        for(t = ir1; t <= i3; t++)
                        {
                            bt = bt+v(r,t)*a(t,s);
                        }
                        v(r,s) = bt;
                    }
                }
                result = ierr==0;
                return result;
            }
            be = 0;
            bc = x-1;
            bg = bc-1;
        }
        flagjump19 = 0;
        i1 = x;
        for(ba = ir1; ba <= i1; ba++)
        {
            u = x+ir1-ba;
            if( u==ir1 )
            {
                break;
            }
            bk = fabs(a(u-1,u-1))+fabs(a(u,u));
            if( bk==0 )
            {
                bk = bu;
            }
            if( fabs(a(u,u-1))<=100*ap::machineepsilon )
            {
                break;
            }
        }
        bn = a(x,x);
        if( u==x )
        {
            a(x,x) = bn+bl;
            evr(x) = a(x,x);
            evi(x) = 0;
            x = bc;
            continue;
        }
        bo = a(bc,bc);
        bm = a(x,bc)*a(bc,x);
        if( u==bc )
        {
            bh = double(bo-bn)/double(2);
            bi = bh*bh+bm;
            bt = sqrt(fabs(bi));
            a(x,x) = bn+bl;
            bn = a(x,x);
            a(bc,bc) = bo+bl;
            if( bi<0 )
            {
                evr(bc) = bn+bh;
                evr(x) = bn+bh;
                evi(bc) = bt;
                evi(x) = -bt;
            }
            else
            {
                if( bh>=0 )
                {
                    bt = bh+fabs(bt);
                }
                else
                {
                    bt = bh-fabs(bt);
                }
                evr(bc) = bn+bt;
                evr(x) = evr(bc);
                if( bt!=0 )
                {
                    evr(x) = bn-bm/bt;
                }
                evi(bc) = 0;
                evi(x) = 0;
                bn = a(x,bc);
                bk = fabs(bn)+fabs(bt);
                bh = bn/bk;
                bi = bt/bk;
                bj = sqrt(bh*bh+bi*bi);
                bh = bh/bj;
                bi = bi/bj;
                i1 = n;
                for(s = bc; s <= i1; s++)
                {
                    bt = a(bc,s);
                    a(bc,s) = bi*bt+bh*a(x,s);
                    a(x,s) = bi*a(x,s)-bh*bt;
                }
                i1 = x;
                for(r = 1; r <= i1; r++)
                {
                    bt = a(r,bc);
                    a(r,bc) = bi*bt+bh*a(r,x);
                    a(r,x) = bi*a(r,x)-bh*bt;
                }
                i1 = ir2;
                for(r = ir1; r <= i1; r++)
                {
                    bt = v(r,bc);
                    v(r,bc) = bi*bt+bh*v(r,x);
                    v(r,x) = bi*v(r,x)-bh*bt;
                }
            }
            x = bg;
            continue;
        }
        if( be==30 )
        {
            ierr = x;
            result = ierr==0;
            return result;
        }
        if( !(be!=10&&be!=20) )
        {
            bl = bl+bn;
            i1 = x;
            for(r = ir1; r <= i1; r++)
            {
                a(r,r) = a(r,r)-bn;
            }
            bk = fabs(a(x,bc))+fabs(a(bc,bg));
            bn = bk*0.75;
            bo = bn;
            bm = -bk*0.4375*bk;
        }
        be = be+1;
        i1 = bg;
        for(bb = u; bb <= i1; bb++)
        {
            w = bg+u-bb;
            bt = a(w,w);
            bj = bn-bt;
            bk = bo-bt;
            bh = (bj*bk-bm)/a(w+1,w)+a(w,w+1);
            bi = a(w+1,w+1)-bt-bj-bk;
            bj = a(w+2,w+1);
            bk = fabs(bh)+fabs(bi)+fabs(bj);
            bh = bh/bk;
            bi = bi/bk;
            bj = bj/bk;
            if( w==u )
            {
                break;
            }
            if( fabs(a(w,w-1))*(fabs(bi)+fabs(bj))<=100*ap::machineepsilon*fabs(bh)*(fabs(a(w-1,w-1))+fabs(bt)+fabs(a(w+1,w+1))) )
            {
                break;
            }
        }
        bf = w+2;
        i1 = x;
        for(r = bf; r <= i1; r++)
        {
            a(r,r-2) = 0;
            if( r==bf )
            {
                continue;
            }
            a(r,r-3) = 0;
        }
        i1 = bc;
        for(t = w; t <= i1; t++)
        {
            bw = t!=bc;
            if( t!=w )
            {
                bh = a(t,t-1);
                bi = a(t+1,t-1);
                bj = 0;
                if( bw )
                {
                    bj = a(t+2,t-1);
                }
                bn = fabs(bh)+fabs(bi)+fabs(bj);
                if( bn==0 )
                {
                    continue;
                }
                bh = bh/bn;
                bi = bi/bn;
                bj = bj/bn;
            }
            r1 = sqrt(bh*bh+bi*bi+bj*bj);
            if( bh>=0 )
            {
                bk = r1;
            }
            else
            {
                bk = -r1;
            }
            if( t==w )
            {
                if( u!=w )
                {
                    a(t,t-1) = -a(t,t-1);
                }
            }
            else
            {
                a(t,t-1) = -bk*bn;
            }
            bh = bh+bk;
            bn = bh/bk;
            bo = bi/bk;
            bt = bj/bk;
            bi = bi/bh;
            bj = bj/bh;
            i2 = n;
            for(s = t; s <= i2; s++)
            {
                bh = a(t,s)+bi*a(t+1,s);
                if( bw )
                {
                    bh = bh+bj*a(t+2,s);
                    a(t+2,s) = a(t+2,s)-bh*bt;
                }
                a(t+1,s) = a(t+1,s)-bh*bo;
                a(t,s) = a(t,s)-bh*bn;
            }
            i2 = x;
            i3 = t+3;
            if( i2<i3 )
            {
                s = i2;
            }
            else
            {
                s = i3;
            }
            i2 = s;
            for(r = 1; r <= i2; r++)
            {
                bh = bn*a(r,t)+bo*a(r,t+1);
                if( bw )
                {
                    bh = bh+bt*a(r,t+2);
                    a(r,t+2) = a(r,t+2)-bh*bj;
                }
                a(r,t+1) = a(r,t+1)-bh*bi;
                a(r,t) = a(r,t)-bh;
            }
            i2 = ir2;
            for(r = ir1; r <= i2; r++)
            {
                bh = bn*v(r,t)+bo*v(r,t+1);
                if( bw )
                {
                    bh = bh+bt*v(r,t+2);
                    v(r,t+2) = v(r,t+2)-bh*bj;
                }
                v(r,t+1) = v(r,t+1)-bh*bi;
                v(r,t) = v(r,t)-bh;
            }
        }
        flagjump19 = 1;
        continue;
    }
    return result;
}


/*************************************************************************
Служебная подпрограмма деления комплексных чисел
*************************************************************************/
void divcomplex(double& e, double& f, double a, double b, double c, double d)
{
    double r;
    double d1;

    if( fabs(c)<fabs(d) )
    {
        r = c/d;
        d1 = d+r*c;
        e = (a*r+b)/d1;
        f = (b*r-a)/d1;
    }
    else
    {
        r = d/c;
        d1 = c+r*d;
        e = (a+b*r)/d1;
        f = (b-a*r)/d1;
    }
}


