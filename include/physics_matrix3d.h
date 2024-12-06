#ifndef PHYSICS_MATRIX_3D_H
#define PHYSICS_MATRIX_3D_H

#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>

#include "physics_vector3d.h"

#define NEARLY_ZERO 1.e-16

class CMatrix3;  // this pre-definition is needed for following functions
CMatrix3 operator*(double d, const CMatrix3& rhs);
CMatrix3 operator*(const CMatrix3& m, double d);
CVector3D operator*(const CMatrix3& m, const CVector3D& v);
CMatrix3 operator+(const CMatrix3& lhs, const CMatrix3& rhs);
CMatrix3 operator*(const CMatrix3& lhs, const CMatrix3& rhs);
CMatrix3 operator-(const CMatrix3& lhs, const CMatrix3& rhs);

// CMatrix3 operator*(double, const CMatrix3&);
// CMatrix3 operator*(const CMatrix3&, double);
// CMatrix3 operator*(const CMatrix3&, const CMatrix3&);

//!< 3x3 matrix
class CMatrix3
{
public:
    CMatrix3(const double s)
    {
        mat[0 * 3 + 0] = s;
        mat[0 * 3 + 1] = 0;
        mat[0 * 3 + 2] = 0;
        mat[1 * 3 + 0] = 0;
        mat[1 * 3 + 1] = s;
        mat[1 * 3 + 2] = 0;
        mat[2 * 3 + 0] = 0;
        mat[2 * 3 + 1] = 0;
        mat[2 * 3 + 2] = s;
    }
    CMatrix3(const CVector3D& vec0) { this->SetSpinTensor(vec0); }
    CMatrix3(const CVector3D& vec0, const CVector3D& vec1)
    {
        this->SetOuterProduct(vec0, vec1);
    }
    CMatrix3(const CVector3D& vec0, const CVector3D& vec1,
             const CVector3D& vec2)
    {
        mat[0 * 3 + 0] = vec0.x;
        mat[0 * 3 + 1] = vec1.x;
        mat[0 * 3 + 2] = vec2.x;
        mat[1 * 3 + 0] = vec0.y;
        mat[1 * 3 + 1] = vec1.y;
        mat[1 * 3 + 2] = vec2.y;
        mat[2 * 3 + 0] = vec0.z;
        mat[2 * 3 + 1] = vec1.z;
        mat[2 * 3 + 2] = vec2.z;
    }
    CMatrix3(double x, double y, double z)
    {
        mat[0 * 3 + 0] = x;
        mat[0 * 3 + 1] = 0;
        mat[0 * 3 + 2] = 0;
        mat[1 * 3 + 0] = 0;
        mat[1 * 3 + 1] = y;
        mat[1 * 3 + 2] = 0;
        mat[2 * 3 + 0] = 0;
        mat[2 * 3 + 1] = 0;
        mat[2 * 3 + 2] = z;
    }
    CMatrix3()
    {
        for (unsigned int i = 0; i < 9; i++) {
            mat[i] = 0;
        }
    }
    CMatrix3(const double m[9])
    {
        for (unsigned int i = 0; i < 9; i++) {
            mat[i] = m[i];
        }
    }
    ////
    void SetDiag(const CVector3D& d)
    {
        mat[0 * 3 + 0] = d.x;
        mat[1 * 3 + 1] = d.y;
        mat[2 * 3 + 2] = d.z;
    }
    void GetElements(double m[9]) const
    {
        for (unsigned int i = 0; i < 9; i++) {
            m[i] = mat[i];
        }
    }
    void GetAffineTransMatElements(double m[16]) const
    {
        m[0 * 4 + 0] = mat[0];
        m[1 * 4 + 0] = mat[1];
        m[2 * 4 + 0] = mat[2];
        m[3 * 4 + 0] = 0;

        m[0 * 4 + 1] = mat[3];
        m[1 * 4 + 1] = mat[4];
        m[2 * 4 + 1] = mat[5];
        m[3 * 4 + 1] = 0;

        m[0 * 4 + 2] = mat[6];
        m[1 * 4 + 2] = mat[7];
        m[2 * 4 + 2] = mat[8];
        m[3 * 4 + 2] = 0;

        m[0 * 4 + 3] = 0;
        m[1 * 4 + 3] = 0;
        m[2 * 4 + 3] = 0;
        m[3 * 4 + 3] = 1;
    }
    CVector3D MatVec(const CVector3D& vec0) const
    {
        CVector3D vec1;
        vec1.x = mat[0] * vec0.x + mat[1] * vec0.y + mat[2] * vec0.z;
        vec1.y = mat[3] * vec0.x + mat[4] * vec0.y + mat[5] * vec0.z;
        vec1.z = mat[6] * vec0.x + mat[7] * vec0.y + mat[8] * vec0.z;
        return vec1;
    }
    void MatVec(const double vec0[], double vec1[]) const
    {
        vec1[0] = mat[0] * vec0[0] + mat[1] * vec0[1] + mat[2] * vec0[2];
        vec1[1] = mat[3] * vec0[0] + mat[4] * vec0[1] + mat[5] * vec0[2];
        vec1[2] = mat[6] * vec0[0] + mat[7] * vec0[1] + mat[8] * vec0[2];
    }
    void MatVecTrans(const double vec0[], double vec1[]) const
    {
        vec1[0] = mat[0] * vec0[0] + mat[3] * vec0[1] + mat[6] * vec0[2];
        vec1[1] = mat[1] * vec0[0] + mat[4] * vec0[1] + mat[7] * vec0[2];
        vec1[2] = mat[2] * vec0[0] + mat[5] * vec0[1] + mat[8] * vec0[2];
    }
    CVector3D MatVecTrans(const CVector3D& vec0) const
    {
        CVector3D vec1;
        vec1.x = mat[0] * vec0.x + mat[3] * vec0.y + mat[6] * vec0.z;
        vec1.y = mat[1] * vec0.x + mat[4] * vec0.y + mat[7] * vec0.z;
        vec1.z = mat[2] * vec0.x + mat[5] * vec0.y + mat[8] * vec0.z;
        return vec1;
    }
    CMatrix3 MatMat(const CMatrix3& mat0) const
    {
        CMatrix3 m;
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                m.mat[i * 3 + j] = mat[i * 3 + 0] * mat0.mat[0 * 3 + j] +
                                   mat[i * 3 + 1] * mat0.mat[1 * 3 + j] +
                                   mat[i * 3 + 2] * mat0.mat[2 * 3 + j];
            }
        }
        return m;
    }
    CMatrix3 MatMatTrans(const CMatrix3& mat0) const
    {
        CMatrix3 m;
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                m.mat[i * 3 + j] = +mat[0 * 3 + i] * mat0.mat[0 * 3 + j] +
                                   mat[1 * 3 + i] * mat0.mat[1 * 3 + j] +
                                   mat[2 * 3 + i] * mat0.mat[2 * 3 + j];
            }
        }
        return m;
    }
    inline const CMatrix3 operator-() const { return (*this) * (-1.0); }
    inline const CMatrix3 operator+() const { return (*this); }
    inline CMatrix3& operator+=(const CMatrix3& rhs)
    {
        for (unsigned int i = 0; i < 9; i++) {
            mat[i] += rhs.mat[i];
        }
        return *this;
    }
    inline CMatrix3& operator-=(const CMatrix3& rhs)
    {
        for (unsigned int i = 0; i < 9; i++) {
            mat[i] -= rhs.mat[i];
        }
        return *this;
    }
    inline CMatrix3& operator*=(double d)
    {
        for (unsigned int i = 0; i < 9; i++) {
            mat[i] *= d;
        }
        return *this;
    }
    void SetRotMatrix_Cartesian(const CVector3D& v)
    {
        const double vec[3] = {v.x, v.y, v.z};
        this->SetRotMatrix_Cartesian(vec);
    }
    void SetRotMatrix_Cartesian(const double vec[])
    {
        double sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        if (sqt < 1.0e-20) {  // infinitesmal rotation approximation
            mat[0] = 1;
            mat[1] = -vec[2];
            mat[2] = +vec[1];
            mat[3] = +vec[2];
            mat[4] = 1;
            mat[5] = -vec[0];
            mat[6] = -vec[1];
            mat[7] = +vec[0];
            mat[8] = 1;
            return;
        }
        double t = sqrt(sqt);
        double invt = 1.0 / t;
        double n[3] = {vec[0] * invt, vec[1] * invt, vec[2] * invt};
        const double c0 = cos(t);
        const double s0 = sin(t);
        mat[0 * 3 + 0] = c0 + (1 - c0) * n[0] * n[0];
        mat[0 * 3 + 1] = -n[2] * s0 + (1 - c0) * n[0] * n[1];
        mat[0 * 3 + 2] = +n[1] * s0 + (1 - c0) * n[0] * n[2];
        mat[1 * 3 + 0] = +n[2] * s0 + (1 - c0) * n[1] * n[0];
        mat[1 * 3 + 1] = c0 + (1 - c0) * n[1] * n[1];
        mat[1 * 3 + 2] = -n[0] * s0 + (1 - c0) * n[1] * n[2];
        mat[2 * 3 + 0] = -n[1] * s0 + (1 - c0) * n[2] * n[0];
        mat[2 * 3 + 1] = +n[0] * s0 + (1 - c0) * n[2] * n[1];
        mat[2 * 3 + 2] = c0 + (1 - c0) * n[2] * n[2];
    }
    void SetRotMatrix_Rodrigues(const double vec[])
    {
        const double sqlen =
            vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        const double tmp1 = 1.0 / (1 + 0.25 * sqlen);
        mat[0] = 1 + tmp1 * (+0.5 * vec[0] * vec[0] - 0.5 * sqlen);
        mat[1] = +tmp1 * (-vec[2] + 0.5 * vec[0] * vec[1]);
        mat[2] = +tmp1 * (+vec[1] + 0.5 * vec[0] * vec[2]);
        mat[3] = +tmp1 * (+vec[2] + 0.5 * vec[1] * vec[0]);
        mat[4] = 1 + tmp1 * (+0.5 * vec[1] * vec[1] - 0.5 * sqlen);
        mat[5] = +tmp1 * (-vec[0] + 0.5 * vec[1] * vec[2]);
        mat[6] = +tmp1 * (-vec[1] + 0.5 * vec[2] * vec[0]);
        mat[7] = +tmp1 * (+vec[0] + 0.5 * vec[2] * vec[1]);
        mat[8] = 1 + tmp1 * (+0.5 * vec[2] * vec[2] - 0.5 * sqlen);
    }
    void SetRotMatrix_CRV(const double crv[])
    {
        const double c0 = 0.125 * (16.0 - crv[0] * crv[0] - crv[1] * crv[1] -
                                   crv[2] * crv[2]);
        const double tmp = 1.0 / ((4.0 - c0) * (4.0 - c0));
        mat[0 * 3 + 0] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[0] * crv[0]);
        mat[0 * 3 + 1] = tmp * (2 * crv[0] * crv[1] - 2 * c0 * crv[2]);
        mat[0 * 3 + 2] = tmp * (2 * crv[0] * crv[2] + 2 * c0 * crv[1]);
        mat[1 * 3 + 0] = tmp * (2 * crv[1] * crv[0] + 2 * c0 * crv[2]);
        mat[1 * 3 + 1] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[1] * crv[1]);
        mat[1 * 3 + 2] = tmp * (2 * crv[1] * crv[2] - 2 * c0 * crv[0]);
        mat[2 * 3 + 0] = tmp * (2 * crv[2] * crv[0] - 2 * c0 * crv[1]);
        mat[2 * 3 + 1] = tmp * (2 * crv[2] * crv[1] + 2 * c0 * crv[0]);
        mat[2 * 3 + 2] = tmp * ((c0 * c0 + 8 * c0 - 16) + 2 * crv[2] * crv[2]);
    }
    void SetRotMatrix_BryantAngle(double rx, double ry, double rz)
    {
        CMatrix3 mx;
        double rvx[3] = {rx, 0, 0};
        mx.SetRotMatrix_Cartesian(rvx);
        CMatrix3 my;
        double rvy[3] = {0, ry, 0};
        my.SetRotMatrix_Cartesian(rvy);
        CMatrix3 mz;
        double rvz[3] = {0, 0, rz};
        mz.SetRotMatrix_Cartesian(rvz);
        CMatrix3 m = mz;
        m = m.MatMat(my);
        m = m.MatMat(mx);
        *this = m;
    }
    void GetCRV_RotMatrix(double crv[]) const
    {
        const double smat[16] = {
            1 + mat[0 * 3 + 0] + mat[1 * 3 + 1] + mat[2 * 3 + 2],
            mat[2 * 3 + 1] - mat[1 * 3 + 2],
            mat[0 * 3 + 2] - mat[2 * 3 + 0],
            mat[1 * 3 + 0] - mat[0 * 3 + 1],
            mat[2 * 3 + 1] - mat[1 * 3 + 2],
            1 + mat[0 * 3 + 0] - mat[1 * 3 + 1] - mat[2 * 3 + 2],
            mat[0 * 3 + 1] + mat[1 * 3 + 0],
            mat[0 * 3 + 2] + mat[2 * 3 + 0],
            mat[0 * 3 + 2] - mat[2 * 3 + 0],
            mat[1 * 3 + 0] + mat[0 * 3 + 1],
            1 - mat[0 * 3 + 0] + mat[1 * 3 + 1] - mat[2 * 3 + 2],
            mat[1 * 3 + 2] + mat[2 * 3 + 1],
            mat[1 * 3 + 0] - mat[0 * 3 + 1],
            mat[0 * 3 + 2] + mat[2 * 3 + 0],
            mat[1 * 3 + 2] + mat[2 * 3 + 1],
            1 - mat[0 * 3 + 0] - mat[1 * 3 + 1] + mat[2 * 3 + 2],
        };

        unsigned int imax;
        imax = (smat[0 * 4 + 0] > smat[1 * 4 + 1]) ? 0 : 1;
        imax = (smat[imax * 4 + imax] > smat[2 * 4 + 2]) ? imax : 2;
        imax = (smat[imax * 4 + imax] > smat[3 * 4 + 3]) ? imax : 3;

        double eparam2[4];  // eular param
        eparam2[imax] = 0.5 * sqrt(smat[imax * 4 + imax]);
        for (unsigned int k = 0; k < 4; k++) {
            if (k == imax) continue;
            eparam2[k] = smat[imax * 4 + k] * 0.25 / eparam2[imax];
        }
        crv[0] = 4 * eparam2[1] / (1 + eparam2[0]);
        crv[1] = 4 * eparam2[2] / (1 + eparam2[0]);
        crv[2] = 4 * eparam2[3] / (1 + eparam2[0]);
    }
    void SetSpinTensor(const CVector3D& vec0)
    {
        mat[0] = 0;
        mat[1] = -vec0.z;
        mat[2] = +vec0.y;
        mat[3] = +vec0.z;
        mat[4] = 0;
        mat[5] = -vec0.x;
        mat[6] = -vec0.y;
        mat[7] = +vec0.x;
        mat[8] = 0;
    }
    void SetOuterProduct(const CVector3D& vec0, const CVector3D& vec1)
    {
        mat[0] = vec0.x * vec1.x;
        mat[1] = vec0.x * vec1.y;
        mat[2] = vec0.x * vec1.z;
        mat[3] = vec0.y * vec1.x;
        mat[4] = vec0.y * vec1.y;
        mat[5] = vec0.y * vec1.z;
        mat[6] = vec0.z * vec1.x;
        mat[7] = vec0.z * vec1.y;
        mat[8] = vec0.z * vec1.z;
    }
    void SetIdentity(double scale = 1)
    {
        mat[0] = scale;
        mat[1] = 0;
        mat[2] = 0;
        mat[3] = 0;
        mat[4] = scale;
        mat[5] = 0;
        mat[6] = 0;
        mat[7] = 0;
        mat[8] = scale;
    }
    CMatrix3 Trans() const
    {
        CMatrix3 m;
        m.mat[0] = mat[0];
        m.mat[1] = mat[3];
        m.mat[2] = mat[6];
        m.mat[3] = mat[1];
        m.mat[4] = mat[4];
        m.mat[5] = mat[7];
        m.mat[6] = mat[2];
        m.mat[7] = mat[5];
        m.mat[8] = mat[8];
        return m;
    }

    static CMatrix3 OuterProduct(const CVector3D& vec0, const CVector3D& vec1)
    {
        CMatrix3 m;
        m.SetOuterProduct(vec0, vec1);
        return m;
    }

    static CMatrix3 Spin(const CVector3D& vec0)
    {
        CMatrix3 m;
        m.SetSpinTensor(vec0);
        return m;
    }
    static CMatrix3 Identity(double scale = 1)
    {
        CMatrix3 m;
        m.SetIdentity(scale);
        return m;
    }

public:
    double mat[9];
};

#endif  // PHYSICS_MATRIX_3D_H
