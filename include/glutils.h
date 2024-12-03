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
    static QVector3D ProjectPoint(const QVector3D & p);
    static void UnProjectPoint(const QVector2D & v, const float depth_value, QVector3D & unproj_pt);

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
    static void DrawArc(float cx, float cy, float cz, float innerR, float outerR, float start_angle, float arc_angle, int num_segments);
    static void DrawSemiRing(const QVector3D & p1, const QVector3D & p2, float start_angle, float arc_angle, const float inner_rad, const float outer_rad, int segNumber);
    static void DrawArcLineStrip(float cx, float cy, float cz, float radius, float start_angle, float arc_angle, int num_segments, const float thickness);
    static void DrawArcLineStrip(QVector3D c, QVector3D dir, QVector3D axis, float arc_angle, int num_segments, float thickness, bool dashed = false);
    static void DrawSemiRingLineStrip(const QVector3D & p1, const QVector3D & p2, float start_angle, float arc_angle, const float radius, const float thickness, int segNumber);

    static void CurvesFromLineSegments(QList <QVector2D> & segments, QList <QList <QVector2D> > & curves);
    static int GetClosestPoint(const QList <QVector2D> & pts, const QVector2D & p);

    inline static void glColor(QVector3D c) { glColor3f(c.x(), c.y(), c.z()); }
    inline static void glColor(QVector4D c) { glColor4f(c.x(), c.y(), c.z(), c.w()); }

private:

    static void ReadPixelColor(const int x, const int y, unsigned char & r, unsigned char & g, unsigned char & b);
    static int point_disp_list;
    static int draw_slices;

    static GLUquadric *quadric;

};

#endif // GLUTILS_H
