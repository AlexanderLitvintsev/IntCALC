/*-----------------------------------------------
��� ������������ ������ ���������� �����������:

void ludecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots);
-----------------------------------------------*/

double determinantlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n);
double determinant(ap::real_2d_array a, int n);
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

void safesolvetriangular(const ap::real_2d_array& a,
     int n,
     ap::real_1d_array& x,
     double& s,
     bool isupper,
     bool istrans,
     bool isunit,
     bool normin,
     ap::real_1d_array& cnorm);
double rcond1(ap::real_2d_array a, int n);
double rcond1lu(const ap::real_2d_array& lu, int n);
double rcondinf(ap::real_2d_array a, int n);
double rcondinflu(const ap::real_2d_array& lu, int n);
void internalestimatercondlu(const ap::real_2d_array& lu,
     int n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rcond);
void internalestimatenorm(int n,
     ap::real_1d_array& v,
     ap::real_1d_array& x,
     ap::integer_1d_array& isgn,
     double& est,
     int& kase);

//*************************************************************************

bool solvesystemlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x);
bool solvesystem(ap::real_2d_array a,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x); 