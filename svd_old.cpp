/************************************************
���� ������ ������������ ������������ AlgoPascal.
************************************************/

#include "ap.h"
#include "svd_old.h"

/*************************************************************************
�������� ������������ ���������� ������� �������� MxN.

�������� ��������� ������� A � �������� � � ����: A = U*W*Transpone(V).
������� U ����� ������ MxN, ������� W - ������������,  �������  V  ����-
�������������.

������� ������:
    *������� �, ������ � 1 �� M, ������� � 1 �� N.
    *M, N - ������ �������.
��������� ������ ����������:
    *������� U �������� ������� A, ������ � 1 �� M, ������� � 1 �� N.
    *��������� ������� W �������� � ���������� W (�������� � 1 �� N).
    *������� V (�� �����������������) �������� � ���������� V. ������
     � 1 �� N, ������� � 1 �� N.

���� �������� �������� ������ ���������, �� ������� ���������� True. ����
�� �������� ��������� �� ������� � ���������� (��������� ������), �� ���
��������� False;
*************************************************************************/
bool svddecomposition_old(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& w,
     ap::real_2d_array& v)
{
    bool result;
    int nm;
    int minmn;
    int l;
    int k;
    int j;
    int jj;
    int its;
    int i;
    double z;
    double y;
    double x;
    double vscale;
    double s;
    double h;
    double g;
    double f;
    double c;
    double anorm;
    ap::real_1d_array rv1;
    bool flag;

    rv1.setbounds(1, n);
    w.setbounds(1, n);
    v.setbounds(1, n, 1, n);
    result = true;
    if( m<n )
    {
        minmn = m;
    }
    else
    {
        minmn = n;
    }
    g = 0.0;
    vscale = 0.0;
    anorm = 0.0;
    for(i = 1; i <= n; i++)
    {
        l = i+1;
        rv1(i) = vscale*g;
        g = 0;
        s = 0;
        vscale = 0;
        if( i<=m )
        {
            for(k = i; k <= m; k++)
            {
                vscale = vscale+fabs(a(k,i));
            }
            if( vscale!=0.0 )
            {
                for(k = i; k <= m; k++)
                {
                    a(k,i) = a(k,i)/vscale;
                    s = s+a(k,i)*a(k,i);
                }
                f = a(i,i);
                g = -extsign(sqrt(s), f);
                h = f*g-s;
                a(i,i) = f-g;
                if( i!=n )
                {
                    for(j = l; j <= n; j++)
                    {
                        s = 0.0;
                        for(k = i; k <= m; k++)
                        {
                            s = s+a(k,i)*a(k,j);
                        }
                        f = s/h;
                        for(k = i; k <= m; k++)
                        {
                            a(k,j) = a(k,j)+f*a(k,i);
                        }
                    }
                }
                for(k = i; k <= m; k++)
                {
                    a(k,i) = vscale*a(k,i);
                }
            }
        }
        w(i) = vscale*g;
        g = 0.0;
        s = 0.0;
        vscale = 0.0;
        if( i<=m&&i!=n )
        {
            for(k = l; k <= n; k++)
            {
                vscale = vscale+fabs(a(i,k));
            }
            if( vscale!=0.0 )
            {
                for(k = l; k <= n; k++)
                {
                    a(i,k) = a(i,k)/vscale;
                    s = s+a(i,k)*a(i,k);
                }
                f = a(i,l);
                g = -extsign(sqrt(s), f);
                h = f*g-s;
                a(i,l) = f-g;
                for(k = l; k <= n; k++)
                {
                    rv1(k) = a(i,k)/h;
                }
                if( i!=m )
                {
                    for(j = l; j <= m; j++)
                    {
                        s = 0.0;
                        for(k = l; k <= n; k++)
                        {
                            s = s+a(j,k)*a(i,k);
                        }
                        for(k = l; k <= n; k++)
                        {
                            a(j,k) = a(j,k)+s*rv1(k);
                        }
                    }
                }
                for(k = l; k <= n; k++)
                {
                    a(i,k) = vscale*a(i,k);
                }
            }
        }
        anorm = mymax(anorm, fabs(w(i))+fabs(rv1(i)));
    }
    for(i = n; i >= 1; i--)
    {
        if( i<n )
        {
            if( g!=0.0 )
            {
                for(j = l; j <= n; j++)
                {
                    v(j,i) = a(i,j)/a(i,l)/g;
                }
                for(j = l; j <= n; j++)
                {
                    s = 0.0;
                    for(k = l; k <= n; k++)
                    {
                        s = s+a(i,k)*v(k,j);
                    }
                    for(k = l; k <= n; k++)
                    {
                        v(k,j) = v(k,j)+s*v(k,i);
                    }
                }
            }
            for(j = l; j <= n; j++)
            {
                v(i,j) = 0.0;
                v(j,i) = 0.0;
            }
        }
        v(i,i) = 1.0;
        g = rv1(i);
        l = i;
    }
    for(i = minmn; i >= 1; i--)
    {
        l = i+1;
        g = w(i);
        if( i<n )
        {
            for(j = l; j <= n; j++)
            {
                a(i,j) = 0.0;
            }
        }
        if( g!=0.0 )
        {
            g = 1.0/g;
            if( i!=n )
            {
                for(j = l; j <= n; j++)
                {
                    s = 0.0;
                    for(k = l; k <= m; k++)
                    {
                        s = s+a(k,i)*a(k,j);
                    }
                    f = s/a(i,i)*g;
                    for(k = i; k <= m; k++)
                    {
                        a(k,j) = a(k,j)+f*a(k,i);
                    }
                }
            }
            for(j = i; j <= m; j++)
            {
                a(j,i) = a(j,i)*g;
            }
        }
        else
        {
            for(j = i; j <= m; j++)
            {
                a(j,i) = 0.0;
            }
        }
        a(i,i) = a(i,i)+1.0;
    }
    for(k = n; k >= 1; k--)
    {
        for(its = 1; its <= maxsvditerations; its++)
        {
            flag = true;
            for(l = k; l >= 1; l--)
            {
                nm = l-1;
                if( fabs(rv1(l))+anorm==anorm )
                {
                    flag = false;
                    break;
                }
                if( fabs(w(nm))+anorm==anorm )
                {
                    break;
                }
            }
            if( flag )
            {
                c = 0.0;
                s = 1.0;
                for(i = l; i <= k; i++)
                {
                    f = s*rv1(i);
                    if( fabs(f)+anorm!=anorm )
                    {
                        g = w(i);
                        h = pythag(f, g);
                        w(i) = h;
                        h = 1.0/h;
                        c = g*h;
                        s = -f*h;
                        for(j = 1; j <= m; j++)
                        {
                            y = a(j,nm);
                            z = a(j,i);
                            a(j,nm) = y*c+z*s;
                            a(j,i) = -y*s+z*c;
                        }
                    }
                }
            }
            z = w(k);
            if( l==k )
            {
                if( z<0.0 )
                {
                    w(k) = -z;
                    for(j = 1; j <= n; j++)
                    {
                        v(j,k) = -v(j,k);
                    }
                }
                break;
            }
            if( its==maxsvditerations )
            {
                result = false;
                return result;
            }
            x = w(l);
            nm = k-1;
            y = w(nm);
            g = rv1(nm);
            h = rv1(k);
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f, 1);
            f = ((x-z)*(x+z)+h*(y/(f+extsign(g, f))-h))/x;
            c = 1.0;
            s = 1.0;
            for(j = l; j <= nm; j++)
            {
                i = j+1;
                g = rv1(i);
                y = w(i);
                h = s*g;
                g = c*g;
                z = pythag(f, h);
                rv1(j) = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = -x*s+g*c;
                h = y*s;
                y = y*c;
                for(jj = 1; jj <= n; jj++)
                {
                    x = v(jj,j);
                    z = v(jj,i);
                    v(jj,j) = x*c+z*s;
                    v(jj,i) = -x*s+z*c;
                }
                z = pythag(f, h);
                w(j) = z;
                if( z!=0.0 )
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g+s*y;
                x = -s*g+c*y;
                for(jj = 1; jj <= m; jj++)
                {
                    y = a(jj,j);
                    z = a(jj,i);
                    a(jj,j) = y*c+z*s;
                    a(jj,i) = -y*s+z*c;
                }
            }
            rv1(l) = 0.0;
            rv1(k) = f;
            w(k) = x;
        }
    }
    return result;
}


/*************************************************************************
��������� ������������
*************************************************************************/
double extsign(double a, double b)
{
    double result;

    if( b>=0 )
    {
        result = fabs(a);
    }
    else
    {
        result = -fabs(a);
    }
    return result;
}


double mymax(double a, double b)
{
    double result;

    if( a>b )
    {
        result = a;
    }
    else
    {
        result = b;
    }
    return result;
}


double pythag(double a, double b)
{
    double result;

    if( fabs(a)<fabs(b) )
    {
        result = fabs(b)*sqrt(1+ap::sqr(a/b));
    }
    else
    {
        result = fabs(a)*sqrt(1+ap::sqr(b/a));
    }
    return result;
}

