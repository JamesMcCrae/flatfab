/*! @file
@brief 3 dimentional vector clas (CVector3D)
@author Nobuyuki Umetani
*/

#ifndef PHYSICS_VECTOR_3D_H
#define PHYSICS_VECTOR_3D_H

#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>

#define NEARLY_ZERO 1.e-16

// rule about naming, the method starts "Set" change it self (not const)
class CVector3D;

double Dot(const CVector3D& arg1, const CVector3D& arg2);
CVector3D Cross(const CVector3D& arg1, const CVector3D& arg2);
CVector3D operator+(const CVector3D& lhs, const CVector3D& rhs);
CVector3D operator-(const CVector3D& lhs, const CVector3D& rhs);
CVector3D operator*(double d, const CVector3D& rhs);
CVector3D operator*(const CVector3D& vec, double d);
double operator*(const CVector3D& lhs, const CVector3D& rhs);
CVector3D operator/(const CVector3D& vec, double d);
CVector3D operator^(const CVector3D& lhs, const CVector3D& rhs);
std::ostream& operator<<(std::ostream& output, const CVector3D& v);
double ScalarTripleProduct(const CVector3D& a, const CVector3D& b,
                           const CVector3D& c);
bool operator==(const CVector3D& lhs, const CVector3D& rhs);
bool operator!=(const CVector3D& lhs, const CVector3D& rhs);
double ScalarTripleProduct3D(const double a[], const double b[],
                             const double c[]);
double Dot3D(const double a[], const double b[]);
void Cross3D(double r[3], const double v1[3], const double v2[3]);
double Length3D(const double v[3]);
void Normalize3D(double v[3]);
double SquareLength3D(const double v[3]);
double SquareDistance3D(const double p0[3], const double p1[3]);
double Distance3D(const double p0[3], const double p1[3]);
double TriArea3D(const double v1[3], const double v2[3], const double v3[3]);
void UnitNormalAreaTri3D(double n[3], double& a, const double v1[3],
                         const double v2[3], const double v3[3]);
void NormalTri3D(double n[3], const double v1[3], const double v2[3],
                 const double v3[3]);
double TetVolume3D(const double v1[3], const double v2[3], const double v3[3],
                   const double v4[3]);
void GetVertical2Vector3D(const double vec_n[3], double vec_x[3],
                          double vec_y[3]);
double Height(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3,
              const CVector3D& v4);
void GetVertical2Vector(const CVector3D& vec_n, CVector3D& vec_x,
                        CVector3D& vec_y);
CVector3D GetMinDist_LinePoint(const CVector3D& p,  // point
                               const CVector3D& s,  // source
                               const CVector3D& d);
CVector3D GetMinDist_LineSegPoint(const CVector3D& p,  // point
                                  const CVector3D& s,  // source
                                  const CVector3D& e);
double TetVolume(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3,
                 const CVector3D& v4);
double TetVolume(int iv1, int iv2, int iv3, int iv4,
                 const std::vector<CVector3D>& node);
void Cross(CVector3D& lhs, const CVector3D& v1, const CVector3D& v2);
double TriArea(const CVector3D& v1, const CVector3D& v2, const CVector3D& v3);
double TriArea(const int iv1, const int iv2, const int iv3,
               const std::vector<CVector3D>& node);
double SquareTriArea(const CVector3D& v1, const CVector3D& v2,
                     const CVector3D& v3);
double SquareDistance(const CVector3D& ipo0, const CVector3D& ipo1);
double SquareLength(const CVector3D& point);
double Length(const CVector3D& point);
double Distance(const CVector3D& ipo0, const CVector3D& ipo1);
double SquareLongestEdgeLength(const CVector3D& ipo0, const CVector3D& ipo1,
                               const CVector3D& ipo2, const CVector3D& ipo3);
double LongestEdgeLength(const CVector3D& ipo0, const CVector3D& ipo1,
                         const CVector3D& ipo2, const CVector3D& ipo3);
double SquareShortestEdgeLength(const CVector3D& ipo0, const CVector3D& ipo1,
                                const CVector3D& ipo2, const CVector3D& ipo3);
double ShortestEdgeLength(const CVector3D& ipo0, const CVector3D& ipo1,
                          const CVector3D& ipo2, const CVector3D& ipo3);
void Normal(CVector3D& vnorm, const CVector3D& v1, const CVector3D& v2,
            const CVector3D& v3);
void UnitNormal(CVector3D& vnorm, const CVector3D& v1, const CVector3D& v2,
                const CVector3D& v3);
double SquareCircumradius(const CVector3D& ipo0, const CVector3D& ipo1,
                          const CVector3D& ipo2, const CVector3D& ipo3);
CVector3D CircumCenter(const CVector3D& ipo0, const CVector3D& ipo1,
                       const CVector3D& ipo2, const CVector3D& ipo3);
double Circumradius(const CVector3D& ipo0, const CVector3D& ipo1,
                    const CVector3D& ipo2, const CVector3D& ipo3);
CVector3D RotateVector(const CVector3D& vec0, const CVector3D& rot);

// 3D vector class
class CVector3D
{
public:
    CVector3D(double vx, double vy, double vz) : x(vx), y(vy), z(vz) {}
    CVector3D() : x(0.0), y(0.0), z(0.0) {}
    CVector3D(const CVector3D& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
    }
    virtual ~CVector3D() {}

    void SetVector(double vx, double vy, double vz)
    {
        x = vx;
        y = vy;
        z = vz;
    }

    inline const CVector3D operator-() const { return -1.0 * (*this); }
    inline const CVector3D operator+() const { return *this; }
    inline CVector3D& operator=(const CVector3D& rhs)
    {
        if (this != &rhs) {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
        }
        return *this;
    }
    inline CVector3D& operator+=(const CVector3D& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
    inline CVector3D& operator-=(const CVector3D& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
    inline CVector3D& operator*=(double d)
    {
        x *= d;
        y *= d;
        z *= d;
        return *this;
    }
    inline CVector3D& operator/=(double d)
    {
        if (fabs(d) < NEARLY_ZERO) {
            return *this;
        }
        x /= d;
        y /= d;
        z /= d;
        return *this;
    }
    inline double operator[](int i) const
    {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        return 0;
    }
    inline double& operator[](int i)
    {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        assert(0);
        return x;
    }
    inline CVector3D operator+() { return *this; }
    inline CVector3D operator-() { return CVector3D(-x, -y, -z); }

    friend bool operator==(const CVector3D&, const CVector3D&);
    friend bool operator!=(const CVector3D&, const CVector3D&);

    friend CVector3D Cross(const CVector3D&, const CVector3D&);
    friend double Dot(const CVector3D&, const CVector3D&);

    inline double Length() const { return sqrt(x * x + y * y + z * z); }
    inline double DLength() const { return x * x + y * y + z * z; }
    void SetNormalizedVector();
    void SetZero();

public:
    double x;  //!< x axis coordinate
    double y;  //!< y axis coordinate
    double z;  //!< z axis coordinate
};

//! 3D bounding box class
class CBoundingBox3D
{
public:
    CBoundingBox3D()
    {
        x_min = 0;
        x_max = 0;
        y_min = 0;
        y_max = 0;
        z_min = 0;
        z_max = 0;
        isnt_empty = false;
    }
    CBoundingBox3D(double x_min0, double x_max0, double y_min0, double y_max0,
                   double z_min0, double z_max0)
        : x_min(x_min0),
          x_max(x_max0),
          y_min(y_min0),
          y_max(y_max0),
          z_min(z_min0),
          z_max(z_max0)
    {
        assert(x_min <= x_max);
        assert(y_min <= y_max);
        assert(z_min <= z_max);
        isnt_empty = true;
    }
    CBoundingBox3D(const CBoundingBox3D& bb)
        : x_min(bb.x_min),
          x_max(bb.x_max),
          y_min(bb.y_min),
          y_max(bb.y_max),
          z_min(bb.z_min),
          z_max(bb.z_max),
          isnt_empty(bb.isnt_empty)
    {
    }
    void SetCenterWidth(double cx, double cy, double cz, double wx, double wy,
                        double wz)
    {
        x_min = cx - wx * 0.5;
        x_max = cx + wx * 0.5;
        y_min = cy - wy * 0.5;
        y_max = cy + wy * 0.5;
        z_min = cz - wz * 0.5;
        z_max = cz + wz * 0.5;
    }
    void GetCenterWidth(double& cx, double& cy, double& cz, double& wx,
                        double& wy, double& wz)
    {
        cx = (x_max + x_min) * 0.5;
        cy = (y_max + y_min) * 0.5;
        cz = (z_max + z_min) * 0.5;
        wx = (x_max - x_min);
        wy = (y_max - y_min);
        wz = (z_max - z_min);
    }
    CBoundingBox3D& operator+=(const CBoundingBox3D& bb)
    {
        if (!bb.isnt_empty) return *this;
        if (!isnt_empty) {
            x_max = bb.x_max;
            x_min = bb.x_min;
            y_max = bb.y_max;
            y_min = bb.y_min;
            z_max = bb.z_max;
            z_min = bb.z_min;
            this->isnt_empty = bb.isnt_empty;
            return *this;
        }
        x_max = (x_max > bb.x_max) ? x_max : bb.x_max;
        x_min = (x_min < bb.x_min) ? x_min : bb.x_min;
        y_max = (y_max > bb.y_max) ? y_max : bb.y_max;
        y_min = (y_min < bb.y_min) ? y_min : bb.y_min;
        z_max = (z_max > bb.z_max) ? z_max : bb.z_max;
        z_min = (z_min < bb.z_min) ? z_min : bb.z_min;
        return *this;
    }
    bool IsIntersect(const CBoundingBox3D& bb) const
    {
        if (!isnt_empty) return false;
        if (!bb.isnt_empty) return false;
        if (x_max < bb.x_min) return false;
        if (x_min > bb.x_max) return false;
        if (y_max < bb.y_min) return false;
        if (y_min > bb.y_max) return false;
        if (z_max < bb.z_min) return false;
        if (z_min > bb.z_max) return false;
        return true;
    }
    bool IsIntersectRay(const CVector3D& org, const CVector3D& dir) const
    {
        {
            double min_r, max_r;
            this->ProjectOnLine(min_r, max_r, org, dir);
            if (max_r < 0) return false;
        }
        if (dir.x > 0 && org.x > x_max) {
            return false;
        }
        if (dir.x < 0 && org.x < x_min) {
            return false;
        }
        if (dir.y > 0 && org.y > y_max) {
            return false;
        }
        if (dir.y < 0 && org.y < y_min) {
            return false;
        }
        if (dir.z > 0 && org.z > z_max) {
            return false;
        }
        if (dir.z < 0 && org.z < z_min) {
            return false;
        }
        const CVector3D& nx = ::Cross(dir, CVector3D(1, 0, 0));
        const CVector3D& ny = ::Cross(dir, CVector3D(0, 1, 0));
        const CVector3D& nz = ::Cross(dir, CVector3D(0, 0, 1));
        {
            double min_r, max_r;
            this->ProjectOnLine(min_r, max_r, org, nx);
            if (max_r * min_r > 0) return false;
        }
        {
            double min_r, max_r;
            this->ProjectOnLine(min_r, max_r, org, ny);
            if (max_r * min_r > 0) return false;
        }
        {
            double min_r, max_r;
            this->ProjectOnLine(min_r, max_r, org, nz);
            if (max_r * min_r > 0) return false;
        }
        return true;
    }
    CBoundingBox3D& operator+=(const CVector3D& v)
    {
        if (!isnt_empty) {
            x_max = v.x;
            x_min = v.x;
            y_max = v.y;
            y_min = v.y;
            z_max = v.z;
            z_min = v.z;
            this->isnt_empty = true;
            return *this;
        }
        x_max = (x_max > v.x) ? x_max : v.x;
        x_min = (x_min < v.x) ? x_min : v.x;
        y_max = (y_max > v.y) ? y_max : v.y;
        y_min = (y_min < v.y) ? y_min : v.y;
        z_max = (z_max > v.z) ? z_max : v.z;
        z_min = (z_min < v.z) ? z_min : v.z;
        return *this;
    }
    bool IsInside(const CVector3D& vec) const
    {
        if (!isnt_empty) return false;
        if (vec.x >= x_min && vec.x <= x_max && vec.y >= y_min &&
            vec.y <= y_max && vec.z >= z_min && vec.z <= z_max)
            return true;
        return false;
    }
    bool IsPossibilityIntersectSphere(const CVector3D& vec,
                                      const double radius) const
    {
        if (!isnt_empty) return false;
        if (vec.x < x_min - radius || vec.x > x_max + radius ||
            vec.y < y_min - radius || vec.y > y_max + radius ||
            vec.z < z_min - radius || vec.z > z_max + radius)
            return false;
        return true;
    }
    bool AddPoint(const CVector3D& vec, double eps)
    {
        if (eps <= 0) {
            return false;
        }
        if (isnt_empty) {
            x_min = (x_min < vec.x - eps) ? x_min : vec.x - eps;
            y_min = (y_min < vec.y - eps) ? y_min : vec.y - eps;
            z_min = (z_min < vec.z - eps) ? z_min : vec.z - eps;
            x_max = (x_max > vec.x + eps) ? x_max : vec.x + eps;
            y_max = (y_max > vec.y + eps) ? y_max : vec.y + eps;
            z_max = (z_max > vec.z + eps) ? z_max : vec.z + eps;
        } else {
            isnt_empty = true;
            x_min = vec.x - eps;
            y_min = vec.y - eps;
            z_min = vec.z - eps;
            x_max = vec.x + eps;
            y_max = vec.y + eps;
            z_max = vec.z + eps;
        }
        return true;
    }
    void SetValueToArray(double bb[8]) const
    {
        bb[0] = x_min;
        bb[2] = y_min;
        bb[4] = z_min;
        bb[1] = x_max;
        bb[3] = y_max;
        bb[5] = z_max;
    }
    void ProjectOnLine(double& min_r, double& max_r, const CVector3D& org,
                       const CVector3D& dir) const
    {
        const double d[8] = {Dot(dir, CVector3D(x_min, y_min, z_min) - org),
                             Dot(dir, CVector3D(x_max, y_min, z_min) - org),
                             Dot(dir, CVector3D(x_min, y_max, z_min) - org),
                             Dot(dir, CVector3D(x_max, y_max, z_min) - org),
                             Dot(dir, CVector3D(x_min, y_min, z_max) - org),
                             Dot(dir, CVector3D(x_max, y_min, z_max) - org),
                             Dot(dir, CVector3D(x_min, y_max, z_max) - org),
                             Dot(dir, CVector3D(x_max, y_max, z_max) - org)};
        min_r = max_r = d[0];
        for (unsigned int i = 1; i < 8; i++) {
            min_r = (d[i] < min_r) ? d[i] : min_r;
            max_r = (d[i] > max_r) ? d[i] : max_r;
        }
    }
    double MinimumDistance(const CVector3D& vec) const
    {
        double x0, y0, z0;
        if (vec.x < x_min) {
            x0 = x_min;
        } else if (vec.x < x_max) {
            x0 = vec.x;
        } else {
            x0 = x_max;
        }
        if (vec.y < y_min) {
            y0 = y_min;
        } else if (vec.y < y_max) {
            y0 = vec.y;
        } else {
            y0 = y_max;
        }
        if (vec.z < z_min) {
            z0 = z_min;
        } else if (vec.z < z_max) {
            z0 = vec.z;
        } else {
            z0 = z_max;
        }
        return (CVector3D(x0, y0, z0) - vec).Length();
    }
    double MaximumDistance(const CVector3D& vec) const
    {
        const double d[8] = {(CVector3D(x_min, y_min, z_min) - vec).Length(),
                             (CVector3D(x_max, y_min, z_min) - vec).Length(),
                             (CVector3D(x_min, y_max, z_min) - vec).Length(),
                             (CVector3D(x_max, y_max, z_min) - vec).Length(),
                             (CVector3D(x_min, y_min, z_max) - vec).Length(),
                             (CVector3D(x_max, y_min, z_max) - vec).Length(),
                             (CVector3D(x_min, y_max, z_max) - vec).Length(),
                             (CVector3D(x_max, y_max, z_max) - vec).Length()};
        double max_d = d[0];
        for (unsigned int i = 1; i < 8; i++) {
            max_d = (d[i] > max_d) ? d[i] : max_d;
        }
        return max_d;
    }

public:
    double x_min, x_max, y_min, y_max, z_min, z_max;
    bool isnt_empty;  //!< false if there is nothing inside
};

#endif  // PHYSICS_VECTOR_3D_H
