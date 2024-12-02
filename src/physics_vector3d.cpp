#include "physics_vector3d.h"

/*!
@brief inner product
*/
double Dot(const CVector3D &arg1, const CVector3D &arg2)
{
    return arg1.x*arg2.x + arg1.y*arg2.y + arg1.z*arg2.z;
}

/*!
@brief cross product
*/
CVector3D Cross(const CVector3D& arg1, const CVector3D& arg2)
{
    CVector3D temp;
    temp.x = arg1.y*arg2.z - arg1.z*arg2.y;
    temp.y = arg1.z*arg2.x - arg1.x*arg2.z;
    temp.z = arg1.x*arg2.y - arg1.y*arg2.x;
    return temp;
}


//! add
CVector3D operator+ (const CVector3D& lhs, const CVector3D& rhs){
    CVector3D temp = lhs;
    temp += rhs;
    return temp;
}

//! subtract
CVector3D operator- (const CVector3D& lhs, const CVector3D& rhs){
    CVector3D temp = lhs;
    temp -= rhs;
    return temp;
}

//! scale
CVector3D operator* (double d, const CVector3D& rhs){
    CVector3D temp = rhs;
    temp *= d;
    return temp;
}

//! scale
CVector3D operator* (const CVector3D& vec, double d){
  CVector3D temp = vec;
  temp *= d;
  return temp;
}

//! mult
double operator* (const CVector3D& lhs, const CVector3D& rhs){
    return Dot(lhs,rhs);
}


//! divide by real number
CVector3D operator/ (const CVector3D& vec, double d){
    CVector3D temp = vec;
    temp /= d;
    return temp;
}

//! mult
CVector3D operator^ (const CVector3D& lhs, const CVector3D& rhs){
    return Cross(lhs,rhs);
}

std::ostream &operator<<(std::ostream &output, const CVector3D& v)
{
  output.setf(std::ios::scientific);
  output << v.x << " " << v.y << " " << v.z;
  return output;
}


double ScalarTripleProduct(const CVector3D& a, const CVector3D& b, const CVector3D& c){
  return a.x*(b.y*c.z - b.z*c.y) + a.y*(b.z*c.x - b.x*c.z) + a.z*(b.x*c.y - b.y*c.x);
}


bool operator== (const CVector3D& lhs, const CVector3D& rhs){
    if( fabs(lhs.x - rhs.x) < NEARLY_ZERO
            && fabs(lhs.y - rhs.y) < NEARLY_ZERO
            && fabs(lhs.z - rhs.z) < NEARLY_ZERO )
        return true;
    else return false;
}

bool operator!= (const CVector3D& lhs, const CVector3D& rhs){
    if( lhs == rhs )	return false;
    else return true;
}



void CVector3D::SetNormalizedVector()
{
    double invmag = 1.0/Length();
    x *= invmag;
    y *= invmag;
    z *= invmag;
}

void CVector3D::SetZero()
{
    x = 0.0;
    y = 0.0;
    z = 0.0;
}

double ScalarTripleProduct3D(const double a[], const double b[], const double c[]){
    return a[0]*(b[1]*c[2] - b[2]*c[1])
        +a[1]*(b[2]*c[0] - b[0]*c[2])
        +a[2]*(b[0]*c[1] - b[1]*c[0]);
}

double Dot3D(const double a[], const double b[]){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void Cross3D(double r[3], const double v1[3], const double v2[3]){
  r[0] = v1[1]*v2[2] - v2[1]*v1[2];
  r[1] = v1[2]*v2[0] - v2[2]*v1[0];
  r[2] = v1[0]*v2[1] - v2[0]*v1[1];
}

double Length3D(const double v[3]){
    return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

void Normalize3D(double v[3]){
    double len = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
}

double SquareLength3D(const double v[3]){
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

double SquareDistance3D(const double p0[3], const double p1[3]){
  return (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
}

double Distance3D(const double p0[3], const double p1[3]){
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

double TriArea3D(const double v1[3], const double v2[3], const double v3[3]){
  double x, y, z;
  x = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  y = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return 0.5*sqrt( x*x + y*y + z*z );
}

void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
    n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
    n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
    n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
    a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
    const double invlen = 0.5/a;
    n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

void NormalTri3D(double n[3], const double v1[3], const double v2[3], const double v3[3]){
  n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
  n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
  n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
}

double TetVolume3D(const double v1[3],
                                 const double v2[3],
                                 const double v3[3],
                                 const double v4[3] )
{
  return
  (   ( v2[0] - v1[0] )*( ( v3[1] - v1[1] )*( v4[2] - v1[2] ) - ( v4[1] - v1[1] )*( v3[2] - v1[2] ) )
   -	( v2[1] - v1[1] )*( ( v3[0] - v1[0] )*( v4[2] - v1[2] ) - ( v4[0] - v1[0] )*( v3[2] - v1[2] ) )
   +	( v2[2] - v1[2] )*( ( v3[0] - v1[0] )*( v4[1] - v1[1] ) - ( v4[0] - v1[0] )*( v3[1] - v1[1] ) )
   ) * 0.16666666666666666666666666666667;
}



void GetVertical2Vector3D
(const double vec_n[3],
 double vec_x[3], double vec_y[3])
{
  const double vec_s[3] = {0,1,0};
  Cross3D(vec_x,vec_s,vec_n);
  const double len = Length3D(vec_x);
  if( len < 1.0e-10 ){
    const double vec_t[3] = {1,0,0};
    Cross3D(vec_x,vec_t,vec_n);  // z????
    Cross3D(vec_y,vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x[0] *= invlen;
    vec_x[1] *= invlen;
    vec_x[2] *= invlen;
    Cross3D(vec_y,vec_n,vec_x);
  }
}

//! Hight of a tetrahedra
double Height(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3, const CVector3D& v4){
    // get normal vector
    double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
    double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
    double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);

    // normalize normal vector
    const double dtmp1 = 1.0 / sqrt( dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z );
    dtmp_x *= dtmp1;
    dtmp_y *= dtmp1;
    dtmp_z *= dtmp1;

    return (v4.x-v1.x)*dtmp_x+(v4.y-v1.y)*dtmp_y+(v4.z-v1.z)*dtmp_z;
}


////////////////////////////////////////////////////////////////

void GetVertical2Vector
(const CVector3D& vec_n,
 CVector3D& vec_x, CVector3D& vec_y)
{
  vec_x = ::Cross(CVector3D(0,1,0),vec_n);
  const double len = vec_x.Length();
  if( len < 1.0e-10 ){
    vec_x = ::Cross(CVector3D(1,0,0),vec_n);  // z????
    vec_x.SetNormalizedVector();
    vec_y = ::Cross(vec_n,vec_x);  // x????
  }
  else{
    const double invlen = 1.0/len;
    vec_x *= invlen;
    vec_y = ::Cross(vec_n,vec_x);
  }
}

CVector3D GetMinDist_LinePoint(const CVector3D& p, // point
                                           const CVector3D& s, // source
                                           const CVector3D& d) // direction
{
  assert( Dot(d,d) > 1.0e-20 );
  const CVector3D ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  return s+t*d;
}

CVector3D GetMinDist_LineSegPoint(const CVector3D& p, // point
                                  const CVector3D& s, // source
                                  const CVector3D& e) // direction
{
  CVector3D d = e-s;
  assert( Dot(d,d) > 1.0e-20 );
  const CVector3D ps = s-p;
  double a = Dot(d,d);
  double b = Dot(d,s-p);
  double t = -b/a;
  if( t < 0 ) t = 0;
  if( t > 1 ) t = 1;
  return s+t*d;
}



//! Volume of a tetrahedra
double TetVolume(const CVector3D& v1,
                        const CVector3D& v2,
                        const CVector3D& v3,
                        const CVector3D& v4 )
{
    return
        (   ( v2.x - v1.x )*( ( v3.y - v1.y )*( v4.z - v1.z ) - ( v4.y - v1.y )*( v3.z - v1.z ) )
        -	( v2.y - v1.y )*( ( v3.x - v1.x )*( v4.z - v1.z ) - ( v4.x - v1.x )*( v3.z - v1.z ) )
        +	( v2.z - v1.z )*( ( v3.x - v1.x )*( v4.y - v1.y ) - ( v4.x - v1.x )*( v3.y - v1.y ) )
        ) * 0.16666666666666666666666666666667;
}

//! élñ ëÃÇÃëÃêœ
double TetVolume( int iv1, int iv2, int iv3, int iv4,
                        const std::vector<CVector3D>& node)
{
    return TetVolume(node[iv1],node[iv2],node[iv3],node[iv4]);
}

////////////////////////////////////////////////

//! äOê⁄ÉxÉNÉgÉã
void Cross( CVector3D& lhs, const CVector3D& v1, const CVector3D& v2 ){
    lhs.x = v1.y*v2.z - v2.y*v1.z;
    lhs.y = v1.z*v2.x - v2.z*v1.x;
    lhs.z = v1.x*v2.y - v2.x*v1.y;
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœ
double TriArea(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3)
{
    double x, y, z;
    x = ( v2.y - v1.y )*( v3.z - v1.z ) - ( v3.y - v1.y )*( v2.z - v1.z );
    y = ( v2.z - v1.z )*( v3.x - v1.x ) - ( v3.z - v1.z )*( v2.x - v1.x );
    z = ( v2.x - v1.x )*( v3.y - v1.y ) - ( v3.x - v1.x )*( v2.y - v1.y );
    return 0.5*sqrt( x*x + y*y + z*z );
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœ
double TriArea(const int iv1, const int iv2, const int iv3,
                      const std::vector<CVector3D>& node )
{
    return TriArea(node[iv1],node[iv2],node[iv3]);
}

//! ÇRéüå≥ÇRäpå`ÇÃñ êœÇÃÇQèÊ
double SquareTriArea(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3)
{
    double dtmp_x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
    double dtmp_y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
    double dtmp_z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
    return (dtmp_x*dtmp_x + dtmp_y*dtmp_y + dtmp_z*dtmp_z)*0.25;
}

////////////////////////////////////////////////

//! í∑Ç≥ÇÃÇQèÊ
double SquareDistance(const CVector3D& ipo0, const CVector3D& ipo1)
{
    return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) + ( ipo1.z - ipo0.z )*( ipo1.z - ipo0.z );
}

//! í∑Ç≥ÇÃÇQèÊ
double SquareLength(const CVector3D& point)
{
    return	point.x*point.x + point.y*point.y + point.z*point.z;
}

////////////////////////////////////////////////

//! length of vector
double Length(const CVector3D& point)
{
    return	sqrt( point.x*point.x + point.y*point.y + point.z*point.z );
}

//! distance between two points
double Distance(const CVector3D& ipo0, const CVector3D& ipo1)
{
    return	sqrt( SquareDistance(ipo0,ipo1) );
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡í∑Ç¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
double SqareLongestEdgeLength(
        const CVector3D& ipo0,
        const CVector3D& ipo1,
        const CVector3D& ipo2,
        const CVector3D& ipo3 )
{
    double edge1, edge2;
    edge1 = SquareDistance( ipo0, ipo1 );
    edge2 = SquareDistance( ipo0, ipo2 );
    if( edge2 > edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo0, ipo3 );
    if( edge2 > edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo1, ipo2 );
    if( edge2 > edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo1, ipo3 );
    if( edge2 > edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo2, ipo3 );
    if( edge2 > edge1 ) edge1 = edge2;
    return edge1;
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡í∑Ç¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
double LongestEdgeLength(
        const CVector3D& ipo0,
        const CVector3D& ipo1,
        const CVector3D& ipo2,
        const CVector3D& ipo3 )
{
    return sqrt( SqareLongestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡íZÇ¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
double SqareShortestEdgeLength(const CVector3D& ipo0,
                          const CVector3D& ipo1,
                          const CVector3D& ipo2,
                          const CVector3D& ipo3 )
{
    double edge1, edge2;
    edge1 = SquareDistance( ipo0, ipo1 );
    edge2 = SquareDistance( ipo0, ipo2 );
    if( edge2 < edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo0, ipo3 );
    if( edge2 < edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo1, ipo2 );
    if( edge2 < edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo1, ipo3 );
    if( edge2 < edge1 ) edge1 = edge2;
    edge2 = SquareDistance( ipo2, ipo3 );
    if( edge2 < edge1 ) edge1 = edge2;
    return edge1;
}

////////////////////////////////////////////////


//! ÇSì_Çå›Ç¢Ç…åãÇ‘ÇUÇ¬ÇÃê¸ï™ÇÃÇ§Çøç≈Ç‡íZÇ¢Ç‡ÇÃÇÃí∑Ç≥Åiélñ ëÃÇÃéøï]âøÇ≈ópÇ¢ÇÈÅj
double ShortestEdgeLength(
        const CVector3D& ipo0,
        const CVector3D& ipo1,
        const CVector3D& ipo2,
        const CVector3D& ipo3 )
{
    return sqrt( SqareShortestEdgeLength(ipo0,ipo1,ipo2,ipo3) );
}

////////////////////////////////////////////////

//! ñ@ê¸ÉxÉNÉgÉã
void Normal(
        CVector3D& vnorm,
        const CVector3D& v1,
        const CVector3D& v2,
        const CVector3D& v3)
{
    vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
    vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
    vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
}

////////////////////////////////////////////////

//! íPà ñ@ê¸ÉxÉNÉgÉã
void UnitNormal(
        CVector3D& vnorm,
        const CVector3D& v1,
        const CVector3D& v2,
        const CVector3D& v3)
{
    vnorm.x = (v2.y-v1.y)*(v3.z-v1.z)-(v2.z-v1.z)*(v3.y-v1.y);
    vnorm.y = (v2.z-v1.z)*(v3.x-v1.x)-(v2.x-v1.x)*(v3.z-v1.z);
    vnorm.z = (v2.x-v1.x)*(v3.y-v1.y)-(v2.y-v1.y)*(v3.x-v1.x);
    const double dtmp1 = 1.0 / Length(vnorm);
    vnorm.x *= dtmp1;
    vnorm.y *= dtmp1;
    vnorm.z *= dtmp1;
}

////////////////////////////////////////////////

/*!
äOê⁄ãÖÇÃîºåa
*/
double SquareCircumradius(
        const CVector3D& ipo0,
        const CVector3D& ipo1,
        const CVector3D& ipo2,
        const CVector3D& ipo3)
{
    double base[3][3] = {
        { ipo1.x-ipo0.x, ipo1.y-ipo0.y, ipo1.z-ipo0.z },
        { ipo2.x-ipo0.x, ipo2.y-ipo0.y, ipo2.z-ipo0.z },
        { ipo3.x-ipo0.x, ipo3.y-ipo0.y, ipo3.z-ipo0.z }
    };
    double s[6] = {
        base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
        base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
        base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
        base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
        base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
        base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
    };
    const double vol = TetVolume(ipo0,ipo1,ipo2,ipo3)*6.0;
    if( vol < 1.0e-20 ){ assert(0); }
    const double inv_det = 1.0 / (vol*vol);
    double t[6] = {
        (s[1]*s[2]-s[3]*s[3])*0.5*inv_det,
        (s[2]*s[0]-s[4]*s[4])*0.5*inv_det,
        (s[0]*s[1]-s[5]*s[5])*0.5*inv_det,
        (s[4]*s[5]-s[0]*s[3])*0.5*inv_det,
        (s[5]*s[3]-s[1]*s[4])*0.5*inv_det,
        (s[3]*s[4]-s[2]*s[5])*0.5*inv_det,
    };
    double u[3] = {
        t[0]*s[0]+t[5]*s[1]+t[4]*s[2],
        t[5]*s[0]+t[1]*s[1]+t[3]*s[2],
        t[4]*s[0]+t[3]*s[1]+t[2]*s[2],
    };
    return  0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
    /*
    const double square_radius = 0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
    CVector3D vec1;
    vec1.x = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.x;
    vec1.y = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.y;
    vec1.z = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.z;
    std::cout << square_radius << " ";
    std::cout << SquareLength(vec1,ipo0) << " ";
    std::cout << SquareLength(vec1,ipo1) << " ";
    std::cout << SquareLength(vec1,ipo2) << " ";
    std::cout << SquareLength(vec1,ipo3) << std::endl;;
    return square_radius;
    */
}

CVector3D CircumCenter(
        const CVector3D& ipo0,
        const CVector3D& ipo1,
        const CVector3D& ipo2,
        const CVector3D& ipo3)
{

    double base[3][3] = {
        { ipo1.x-ipo0.x, ipo1.y-ipo0.y, ipo1.z-ipo0.z },
        { ipo2.x-ipo0.x, ipo2.y-ipo0.y, ipo2.z-ipo0.z },
        { ipo3.x-ipo0.x, ipo3.y-ipo0.y, ipo3.z-ipo0.z }
    };
    double s[6] = {
        base[0][0]*base[0][0]+base[0][1]*base[0][1]+base[0][2]*base[0][2],
        base[1][0]*base[1][0]+base[1][1]*base[1][1]+base[1][2]*base[1][2],
        base[2][0]*base[2][0]+base[2][1]*base[2][1]+base[2][2]*base[2][2],
        base[1][0]*base[2][0]+base[1][1]*base[2][1]+base[1][2]*base[2][2],
        base[2][0]*base[0][0]+base[2][1]*base[0][1]+base[2][2]*base[0][2],
        base[0][0]*base[1][0]+base[0][1]*base[1][1]+base[0][2]*base[1][2],
    };
    const double vol = TetVolume(ipo0,ipo1,ipo2,ipo3)*6.0;
    if( vol < 1.0e-20 ){ assert(0); }
    const double inv_det = 1.0 / (vol*vol);
    double t[6] = {
        (s[1]*s[2]-s[3]*s[3])*0.5*inv_det,
        (s[2]*s[0]-s[4]*s[4])*0.5*inv_det,
        (s[0]*s[1]-s[5]*s[5])*0.5*inv_det,
        (s[4]*s[5]-s[0]*s[3])*0.5*inv_det,
        (s[5]*s[3]-s[1]*s[4])*0.5*inv_det,
        (s[3]*s[4]-s[2]*s[5])*0.5*inv_det,
    };
    double u[3] = {
        t[0]*s[0]+t[5]*s[1]+t[4]*s[2],
        t[5]*s[0]+t[1]*s[1]+t[3]*s[2],
        t[4]*s[0]+t[3]*s[1]+t[2]*s[2],
    };
//    const double square_radius = 0.5*(u[0]*s[0]+u[1]*s[1]+u[2]*s[2]);
    CVector3D vec1;
    vec1.x = base[0][0]*u[0]+base[1][0]*u[1]+base[2][0]*u[2] + ipo0.x;
    vec1.y = base[0][1]*u[0]+base[1][1]*u[1]+base[2][1]*u[2] + ipo0.y;
    vec1.z = base[0][2]*u[0]+base[1][2]*u[1]+base[2][2]*u[2] + ipo0.z;
  return vec1;
}

////////////////////////////////////////////////

/*!
curcumradius of a tetrahedra
*/
double Circumradius(const CVector3D& ipo0,
                           const CVector3D& ipo1,
                           const CVector3D& ipo2,
                           const CVector3D& ipo3){
    return sqrt( SquareCircumradius(ipo0,ipo1,ipo2,ipo3) );
}


CVector3D RotateVector(const CVector3D& vec0, const CVector3D& rot )
{
    const double theta = rot.Length();
    if( theta < 1.0e-30 ){
        return vec0;
    }
    CVector3D e0 = rot;
    e0.SetNormalizedVector();
    CVector3D e2 = ::Cross(e0,vec0);
    if( e2.Length() < 1.0e-30 ){
        return vec0;
    }
    e2.SetNormalizedVector();
    CVector3D e1 = ::Cross(e2,e0);
    assert( fabs( e1.Length() - 1 ) < 1.0e-10 );
//	assert( e2.x*vec_0.x + e2.y*vec_0.y + e2.z*vec_0.z < 1.0e-10 );
    const double dot00 = Dot(vec0,e0);
    const double dot01 = Dot(vec0,e1);
    const double cost = cos(theta);
    const double sint = sin(theta);
    CVector3D vec1;
    vec1.x = dot00*e0.x + dot01*cost*e1.x + dot01*sint*e2.x;
    vec1.y = dot00*e0.y + dot01*cost*e1.y + dot01*sint*e2.y;
    vec1.z = dot00*e0.z + dot01*cost*e1.z + dot01*sint*e2.z;
    return vec1;
}
