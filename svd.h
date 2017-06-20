bool svddecomposition(ap::real_2d_array a,
     int m,
     int n,
     int uneeded,
     int vtneeded,
     ap::real_1d_array& w,
     ap::real_2d_array& u,
     ap::real_2d_array& vt);

void applyrotationsfromtheleft(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_2d_array& a,
     ap::real_1d_array& work);
void applyrotationsfromtheright(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_2d_array& a,
     ap::real_1d_array& work);
void generaterotation(double f, double g, double& cs, double& sn, double& r);
void generatereflection(ap::real_1d_array& x, int n, double& tau);
void applyreflectionfromtheleft(ap::real_2d_array& c,
     double tau,
     const ap::real_1d_array& v,
     int m1,
     int m2,
     int n1,
     int n2,
     ap::real_1d_array& work);
void applyreflectionfromtheright(ap::real_2d_array& c,
     double tau,
     const ap::real_1d_array& v,
     int m1,
     int m2,
     int n1,
     int n2,
     ap::real_1d_array& work);
void qrdecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& tau);
void unpackqfromqr(const ap::real_2d_array& qr,
     int m,
     int n,
     const ap::real_1d_array& tau,
     int qcolumns,
     ap::real_2d_array& q);
void qrdecompositionunpacked(ap::real_2d_array a,
     int m,
     int n,
     ap::real_2d_array& q,
     ap::real_2d_array& r);
void lqdecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& tau);
void unpackqfromlq(const ap::real_2d_array& lq,
     int m,
     int n,
     const ap::real_1d_array& tau,
     int qrows,
     ap::real_2d_array& q);
void lqdecompositionunpacked(ap::real_2d_array a,
     int m,
     int n,
     ap::real_2d_array& l,
     ap::real_2d_array& q);
bool bidiagonalsvddecomposition(ap::real_1d_array& d,
     ap::real_1d_array e,
     int n,
     bool isupper,
     bool isfractionalaccuracyrequired,
     ap::real_2d_array& u,
     int nru,
     ap::real_2d_array& vt,
     int ncvt);
double extsignbdsqr(double a, double b);
void svd2x2(double f, double g, double h, double& ssmin, double& ssmax);
void svdv2x2(double f,
     double g,
     double h,
     double& ssmin,
     double& ssmax,
     double& snr,
     double& csr,
     double& snl,
     double& csl);
void tobidiagonal(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& tauq,
     ap::real_1d_array& taup);
void unpackqfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& tauq,
     int qcolumns,
     ap::real_2d_array& q);
void unpackptfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& taup,
     int ptrows,
     ap::real_2d_array& pt);
void unpackdiagonalsfrombidiagonal(const ap::real_2d_array& b,
     int m,
     int n,
     bool& isupper,
     ap::real_1d_array& d,
     ap::real_1d_array& e);

void balancematrix(ap::real_2d_array& a, int n, ap::real_1d_array& d);

static const int radix = 2;
