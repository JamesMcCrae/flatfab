#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include "algebra3.h"

class Triangulate
{
public:
    // triangulate a contour/polygon, places results in STL vector
    // as series of triangles.
    static bool Process(const QList<vec2> & contour, QList<vec2> & result);

    // compute area of a contour/polygon
    static float Area(const QList<vec2> & contour);

    // decide if point Px/Py is inside triangle defined by
    // (Ax,Ay) (Bx,By) (Cx,Cy)
    static bool InsideTriangle(float Ax, float Ay, float Bx, float By, float Cx,
                               float Cy, float Px, float Py);

private:
    static bool Snip(const QList<vec2> & contour, int u, int v, int w, int n,
                     const QVector<int> & V);
};

#endif  // TRIANGULATE_H
