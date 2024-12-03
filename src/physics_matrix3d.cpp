#include "physics_matrix3d.h"

CMatrix3 operator* (double d, const CMatrix3& rhs){
  CMatrix3 temp = rhs;
  temp *= d;
  return temp;
}

CMatrix3 operator* (const CMatrix3& m, double d){
  CMatrix3 temp = m;
  temp *= d;
  return temp;
}

CVector3D operator* (const CMatrix3& m, const CVector3D& v){
  CVector3D r = m.MatVec(v);
  return r;
}

CMatrix3 operator+ (const CMatrix3& lhs, const CMatrix3& rhs){
  CMatrix3 temp = lhs;
  temp += rhs;
  return temp;
}

CMatrix3 operator* (const CMatrix3& lhs, const CMatrix3& rhs){
  return lhs.MatMat(rhs);
}

CMatrix3 operator- (const CMatrix3& lhs, const CMatrix3& rhs){
  CMatrix3 temp = lhs;
  temp -= rhs;
  return temp;
}
