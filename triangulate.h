#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <QtGui>

#include "glutils.h"

extern "C" {
#include "triangle.h"
}

class Triangulate
{
public:

  // triangulate a contour/polygon, places results in STL vector
  // as series of triangles.
  static bool Process(const QList <QVector2D> & contour, QList <QVector2D> & result);
  static bool Process(const QList <QList <QVector2D> > & contours, QList <QVector2D> & result);

  // compute area of a contour/polygon
  static float Area(const QList <QVector2D> & contour);

  // decide if point Px/Py is inside triangle defined by
  // (Ax,Ay) (Bx,By) (Cx,Cy)
  static bool InsideTriangle(float Ax, float Ay,
                      float Bx, float By,
                      float Cx, float Cy,
                      float Px, float Py);

  static QVector2D Centroid(const QList <QVector2D> & contour);


private:
  static bool Snip(const QList <QVector2D> & contour, int u, int v ,int w, int n,  const QVector <int> & V);

};

#endif // TRIANGULATE_H
