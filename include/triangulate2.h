#ifndef TRIANGULATE2_H
#define TRIANGULATE2_H

#include <QtGui>
#include <cfloat>

class Triangulate2
{

public:

    Triangulate2();

    static bool Process(const QList <QVector2D> & contour, QList <QVector2D> & result, QList <QVector2D> & poly);
    static bool Process(const QList <QList <QVector2D> > & contours, QList <QVector2D> & result, QList <QVector2D> & poly);

private:

    static bool GetClockwise(const QList <QVector2D> & contour);
    static void Reverse(QList <QVector2D> & list);

    static float Cross(const QVector2D & a, const QVector2D & b);
    static bool IsReflexAngle(const QVector2D & a, const QVector2D & b, const QVector2D & c);
    static bool IsReflexIndex(const QList <QVector2D> & coords, const QList <int> polygon, const int pcurr);
    static float IntersectSegmentX(const QVector2D & p0, const QVector2D & p1, const float y);
    static void GetNeighbours(const QList <int> & polygon, const int pcurr, int & pprev, int & pnext);
    static bool PointInTri(const QVector2D tri[3], const QVector2D p);
    static void GetSlice(const QList <QVector2D> & coords, const QList <int> & polygon, const QList <QVector2D> & hole, int slice[2]);
    static bool CheckEar(const QList <QVector2D> & coords, const QList <int> & polygon, const QList <bool> & reflex, const int reflexCount, const int pcurr);

    static bool Process(QList <QVector2D> & coords, const QList <QList <QVector2D> > & holes, QList <QVector2D> & triangles, QList <QVector2D> & polygon);

};

#endif // TRIANGULATE2_H
