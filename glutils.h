#ifndef GLUTILS_H
#define GLUTILS_H

#include <QtOpenGL>
#include <QtGui>
#include <cfloat>

#ifdef __APPLE__
    #include <OpenGL/glu.h>
#else
    #include <GL/glu.h>
#endif

//#include "algebra3.h"
//#include "tnbframe.h"

class GLutils
{

public:

    GLutils();

    //drawing related
    static void DrawPointGL(const float x, const float y, const float rad);
    //static void DrawTNBFrame(const vec3 & o, const TNBFrame & tnb, const float scale);
    static void DrawArrow(const QVector3D & p1, const QVector3D & p2);
    static void DrawArrowFixedLength(const QVector3D & p1, const QVector3D & p2, const float len);
    static void DrawCylinder(const QVector3D & p1, const QVector3D & p2, const float radius);
    static void DrawCylinderFixedLength(const QVector3D & p1, const QVector3D & p2, const float radius, const float len);
    static void DrawRing(const QVector3D & p1, const QVector3D & p2, const float inner_rad, const float outer_rad);
    static void DrawSphere(const QVector3D & p, const float radius);


    //drawing (setting colour)
    static void SetColorByHue(const float val, const float range);
    static QVector3D GetColorByHue(const float val, const float range);
    static void ColorByIndex(const int ind);
    static QString ColorByIndexStr(const int ind);
    static void ColorByNearestMajorAxis(const QVector3D & v);

    //drawing (blending)
    static void EnableBlending();
    static void DisableBlending();

    //projection related
    //static vec2 ProjectPoint(const vec3 & p);
    static QVector3D ProjectPoint(const QVector3D & p);
    //static QVector <vec2> ProjectPoints(const QVector <vec3> & pts);
    //static vec2 ProjectVector(const vec3 & p1, const vec3 & p2);
    //static vec2 ProjectUnitVector(const vec3 & p1, const vec3 & p2);
    //static void ProjectCurve(const QList <vec4> & c, QList <vec4> & proj_c);
    //static void UnProjectPoint(const vec2 & p, const float depth_value, vec3 & unproj_pt);
    static void UnProjectPoint(const QVector2D & v, const float depth_value, QVector3D & unproj_pt);
    //static void UnProjectCurve(const QList <vec2> & c, const QVector <float> depth_values, QList <vec3> & unproj_curve);

    static void ReadPixelColor_FrontBuffer(const int x, const int y, unsigned char & r, unsigned char & g, unsigned char & b);
    static void ReadPixelColor_BackBuffer(const int x, const int y, unsigned char & r, unsigned char & g, unsigned char & b);
    static void PixelColorToIndex(const unsigned char r, const unsigned char g, const unsigned char b, int & index);
    static void IndexToPixelColor(const int index, unsigned char & r, unsigned char & g, unsigned char & b);
    static void SetPickColor(const int index);

    static int GetWindowWidth();
    static int GetWindowHeight();
    static void GetViewport(int v[4]);
    static void SetViewport(const int v[4]);
    static void SetViewport(const int x, const int y, const int w, const int h);

    //vector and geometry useful methods
    static QVector3D RotateVector(const QVector3D & v, const QVector3D & axis, const float angle_radians);
    static QVector3D MirrorVector(const QVector3D & v, const QVector3D & mirror_dir);
    static void GetBoundingBox(const QVector <QVector3D> & vs, QVector3D & bbox_min, QVector3D & bbox_max);
    static float AngleBetweenDeg(const QVector3D & v1, const QVector3D & v2);
    static float AngleBetweenRad(const QVector3D & v1, const QVector3D & v2);
    static float SignedAngleBetweenRad(const QVector3D & v1, const QVector3D & v2, const QVector3D & axis);
    static bool LineLineIntersection(const QVector2D & p1, const QVector2D & p2, const QVector2D & p3, const QVector2D & p4, QVector2D & intersect);
    static bool LineRayIntersection(const QVector2D & p1, const QVector2D & p2, const QVector2D & lp, const QVector2D & ld, QVector2D & intersect);
    static bool LineBoxIntersection(const QVector2D & p1, const QVector2D & p2, const QVector2D & lp, const QVector2D & ld);
    static bool LinePlaneIntersection(const QVector3D & p0, const QVector3D & n, const QVector3D & l0, const QVector3D & l1, QVector3D & intersect);
    static bool LineSegmentPlaneIntersection(const QVector3D & p0, const QVector3D & n, const QVector3D & l0, const QVector3D & l1, QVector3D & intersect);
    static float PointLineSignedDistance(const QVector2D & l1, const QVector2D & l2, const QVector2D & p);
    static float PointLineSegmentDistance(const QVector3D & p0, const QVector3D & p1, const QVector3D & p2);
    static float PointLineSegmentDistance(const QVector3D & p0, const QVector3D & p1, const QVector3D & p2, QVector3D & best_pt);
    static bool PlanePlaneIntersection(const QVector3D & n1, const QVector3D & p1, const QVector3D & n2, const QVector3D & p2, QVector3D & lp, QVector3D & ld);
    static void ConvexHull_GiftWrapping(const QList <QVector3D> & pts, QList <int> & hull);
    static bool ConvexHull_PointInside(const QList <QVector3D> & hull, const QVector3D & p);    
    static int GetClosestPoint(const QList <QVector3D> & pts, const QVector3D & p);
    static void GetClosestPairOfPoints(const QList <QVector3D> & pts1, const QList <QVector3D> & pts2, int & ind1, int & ind2);
    static QVector3D GetVectorNewBasis(const QVector3D & old_x, const QVector3D & old_y, const QVector3D & old_z, const QVector3D & old_p,
                                       const QVector3D & new_x, const QVector3D & new_y, const QVector3D & new_z, const QVector3D & new_p,
                                       const QVector3D & p, const float scale_x, const float scale_y, const float scale_z);
    static void SortPointsAlongDirection2D(const QVector2D & dir, QList <QVector2D> & isecs);
    static void SortPointsAlongDirection3D(const QVector3D & dir, QList <QVector3D> & pts);
    static void SortPointsAlongDirection3DExtra(const QVector3D dir, QList <QVector3D> & isecs, QList <bool> & isecs_which);

    static void GetSortedIntersectionPoints(const QList <QVector2D> & pts, const QVector2D ray_dir, const QVector2D ray_p, QList <QVector2D> & sorted_isecs);

    static QVector3D GetOrthoVec(const QVector3D & v);
    static void GetModelviewMatrixFromVec(const QVector3D & z, const QVector3D & p, float modelview[16]);


    static void DrawDisc(const QVector3D & p1, const QVector3D & p2, const float inner_rad, const float outer_rad);


private:

    static void ReadPixelColor(const int x, const int y, unsigned char & r, unsigned char & g, unsigned char & b);
    static int point_disp_list;
    static int draw_slices;

    static GLUquadric *quadric;

};

#endif // GLUTILS_H
