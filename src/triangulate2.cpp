#include "triangulate2.h"

Triangulate2::Triangulate2() {}

bool Triangulate2::GetClockwise(const QList<QVector2D> & contour)
{
    float total = 0.0f;

    for (int i0 = 0; i0 < contour.size(); ++i0) {
        const int i1 = (i0 + 1) % contour.size();
        const int i2 = (i0 + 2) % contour.size();

        const QVector2D v0 = contour[i1] - contour[i0];
        const QVector2D v1 = contour[i2] - contour[i1];

        total += (v0.x() * v1.y()) - (v1.x() * v0.y());
    }

    return (total > 0.0f);
}

float Triangulate2::Cross(const QVector2D & a, const QVector2D & b)
{
    return a.x() * b.y() - a.y() * b.x();
}

bool Triangulate2::IsReflexAngle(const QVector2D & a, const QVector2D & b,
                                 const QVector2D & c)
{
    QVector2D ac = c - a;
    QVector2D ab = b - a;
    return (0 > Cross(ac, ab));
}

bool Triangulate2::IsReflexIndex(const QList<QVector2D> & coords,
                                 const QList<int> polygon, const int pcurr)
{
    int pprev;
    int pnext;
    GetNeighbours(polygon, pcurr, pprev, pnext);
    QVector2D a = coords[polygon[pprev]];
    QVector2D b = coords[polygon[pcurr]];
    QVector2D c = coords[polygon[pnext]];
    return IsReflexAngle(a, b, c);
}

float Triangulate2::IntersectSegmentX(const QVector2D & p0,
                                      const QVector2D & p1, const float y)
{
    if (p0.y() == p1.y()) {
        return p0.x();
    }

    if (p0.y() < p1.y()) {
        const float t = (y - p0.y()) / (p1.y() - p0.y());
        return p0.x() + t * (p1.x() - p0.x());
    } else {
        const float t = (y - p1.y()) / (p0.y() - p1.y());
        return p1.x() + t * (p0.x() - p1.x());
    }
}

void Triangulate2::GetNeighbours(const QList<int> & polygon, const int pcurr,
                                 int & pprev, int & pnext)
{
    pprev = (pcurr + polygon.size() - 1) % polygon.size();
    pnext = (pcurr + 1) % polygon.size();
}

bool Triangulate2::PointInTri(const QVector2D tri[3], const QVector2D p)
{
    QVector2D ab = tri[1] - tri[0];
    QVector2D bc = tri[2] - tri[1];
    QVector2D ca = tri[0] - tri[2];
    QVector2D ap = p - tri[0];
    QVector2D bp = p - tri[1];
    QVector2D cp = p - tri[2];
    float a = Cross(ab, ap);
    float b = Cross(bc, bp);
    float c = Cross(ca, cp);
    return (a < 0 && b < 0 && c < 0) || (a > 0 && b > 0 && c > 0);
}

void Triangulate2::GetSlice(const QList<QVector2D> & coords,
                            const QList<int> & polygon,
                            const QList<QVector2D> & hole, int slice[2])
{
    float min_dist = FLT_MAX;
    slice[0] = -1;
    slice[1] = -1;

    for (int i = 0; i < polygon.size(); ++i) {
        for (int j = 0; j < hole.size(); ++j) {
            const float len_sq = (coords[polygon[i]] - hole[j]).lengthSquared();
            if (len_sq < min_dist) {
                min_dist = len_sq;
                slice[0] = i;
                slice[1] = j;
            }
        }
    }
}

bool Triangulate2::CheckEar(const QList<QVector2D> & coords,
                            const QList<int> & polygon,
                            const QList<bool> & reflex, const int reflexCount,
                            const int pcurr)
{
    if (reflexCount == 0) {
        return true;
    }

    int pprev, pnext;
    GetNeighbours(polygon, pcurr, pprev, pnext);

    int ptriangle[3];
    ptriangle[0] = pprev;
    ptriangle[1] = pcurr;
    ptriangle[2] = pnext;

    int ntriangle[3];
    ntriangle[0] = polygon[ptriangle[0]];
    ntriangle[1] = polygon[ptriangle[1]];
    ntriangle[2] = polygon[ptriangle[2]];

    QVector2D tricoords[3];
    tricoords[0] = coords[ntriangle[0]];
    tricoords[1] = coords[ntriangle[1]];
    tricoords[2] = coords[ntriangle[2]];

    bool isEar = true;

    for (int p = 0; p < polygon.size(); ++p) {
        int n = polygon[p];

        if (ntriangle[0] == n || ntriangle[1] == n || ntriangle[2] == n) {
            continue;
        }

        if (p >= reflex.size() || !reflex[p]) {
            continue;
        }

        if (PointInTri(tricoords, coords[n])) {
            isEar = false;
            break;
        }
    }

    return isEar;
}

void Triangulate2::Reverse(QList<QVector2D> & list)
{
    for (int i = 0; i < list.size() / 2; ++i) {
        qSwap(list[i], list[list.size() - i - 1]);
    }
}

bool Triangulate2::Process(const QList<QVector2D> & contour,
                           QList<QVector2D> & result, QList<QVector2D> & poly)
{
    QList<QVector2D> coords = contour;
    QList<QList<QVector2D> > holes;
    if (GetClockwise(coords)) {
        Reverse(coords);
    }

    return Process(coords, holes, result, poly);
}

bool Triangulate2::Process(const QList<QList<QVector2D> > & contours,
                           QList<QVector2D> & result, QList<QVector2D> & poly)
{
    if (contours.empty()) {
        return false;
    }

    QList<QVector2D> coords = contours.first();
    QList<QList<QVector2D> > holes;

    if (GetClockwise(coords)) {
        Reverse(coords);
    }

    for (int i = 1; i < contours.size(); ++i) {
        QList<QVector2D> hole = contours[i];
        if (!GetClockwise(hole)) {
            // qDebug() << "hole" << i << "clockwise false";
            Reverse(hole);
        } else {
            // qDebug() << "hole" << i << "clockwise true";
        }
        holes.push_back(hole);
    }

    return Process(coords, holes, result, poly);
}

bool Triangulate2::Process(QList<QVector2D> & coords,
                           const QList<QList<QVector2D> > & holes,
                           QList<QVector2D> & triangles,
                           QList<QVector2D> & poly)
{
    // QVector2D ab, bc, ca, ap, bp, cp, ac;
    triangles.clear();

    QList<bool> reflex;

    // polygon = [0...coords.length] // the "..." means exclusive not inclusive
    // range
    QList<int> polygon;
    for (int i = 0; i < coords.size(); ++i) {
        polygon.push_back(i);
    }
    int reflexCount = 0;

    // 1.  handle trivial cases
    if (coords.size() < 3) {
        return false;
    }

    if (coords.size() == 3 && holes.empty()) {
        triangles.push_back(coords[0]);
        triangles.push_back(coords[1]);
        triangles.push_back(coords[2]);
        return true;
    }

    // Repair the polygon if it has holes by duplicating verts along
    // a new edge called "slice".  The slice verts must be
    // visible to each other.  (ie, no edge intersections)
    for (int c = 0; c < holes.size(); ++c) {
        if (holes[c].size() < 3) {
            continue;
        }

        reflex.clear();
        for (int p = 0; p < polygon.size(); ++p) {
            // int n = polygon[p];
            reflex.push_back(IsReflexIndex(coords, polygon, p));
        }

        const QList<QVector2D> & hole =
            holes[c];  // NOTE: this needs to generalize to more than 1 hole

        // Find any two mutually visible vertices, the first
        // from the outer contour, the second from the hole.
        int slice[2];
        GetSlice(coords, polygon, hole, slice);
        // qDebug() << "Slicing with index" << slice[0] << slice[1];

        // Clone the outer contour and append the hole verts.
        // QList <QVector2D> coords2 = coords;
        int holeStart = coords.size();
        coords += hole;

        // Perform a rotational shift of the indices in 'polygon'
        // such that the slice vertex winds up at the end of the list.
        QList<int> newPolygon;
        int i = (slice[0] + 1) % polygon.size();
        for (int p = 0; p < polygon.size(); ++p) {
            newPolygon.push_back(polygon[i]);
            i = (i + 1) % polygon.size();
        }

        // Similarly shift the indices of 'hole' and append
        // them to the new polygon.
        i = slice[1];
        for (int j = 0; j < hole.size(); ++j) {
            newPolygon.push_back(holeStart + i);
            i = (i + 1) % hole.size();
        }

        // Insert the two duplicated verts that occur along
        // the new "slice" edge.
        newPolygon.push_back(newPolygon[polygon.length()]);
        newPolygon.push_back(newPolygon[polygon.length() - 1]);
        // qDebug() << polygon << newPolygon;
        polygon = newPolygon;
    }

    // We're now ready for the ear clipping algorithm.  First,
    // find all reflex verts.
    QList<int> convex;
    reflex.clear();
    reflexCount = 0;

    for (int p = 0; p < polygon.size(); ++p) {
        // int n = polygon[p];
        if (IsReflexIndex(coords, polygon, p)) {
            reflex.push_back(true);
            ++reflexCount;
        } else {
            reflex.push_back(false);
            convex.push_back(p);
        }
    }

    // Next find all the initial ears, which are verts that form triangles that
    // don't contain any other verts.  This is a n-squared operation.
    QList<int> ears;
    for (int j = 0; j < convex.size(); ++j) {
        int p = convex[j];
        if (CheckEar(coords, polygon, reflex, reflexCount, p)) {
            ears.push_back(p);
        }
    }

    // Diagnostic output.
    // verbose = false
    // if verbose
    // console.info ""
    // console.info "ears    #{ears}"
    // console.info "reflex  #{reflex}"
    // console.info "convex  #{convex}"
    // qDebug() << "Ears:" << ears.size() << "Reflex:" << reflex.size() <<
    // reflexCount << "Convex:" << convex.size() << "Coords:" << coords.size();

    poly.clear();
    for (int i = 0; i < polygon.size(); ++i) {
        poly.push_back(coords[polygon[i]]);
    }

    // Remove ears, one by one.
    // int watchdog = 100000;
    triangles.clear();
    while (polygon.size() > 0) {
        // qDebug() << "Polygon:" << polygon.size() << "Ears:" << ears.size();

        // Remove the index from the ear list.
        if (ears.empty()) {
            // qDebug() << "Triangulate2::Process() - Warning, ears list empty,
            // halting";
            break;
        }
        int pcurr = ears.last();
        ears.pop_back();
        //--watchdog;

        // Insert the ear into the triangle list that we're building.
        int pprev, pnext;
        GetNeighbours(polygon, pcurr, pprev, pnext);
        int ptriangle[3];
        ptriangle[0] = pprev;
        ptriangle[1] = pcurr;
        ptriangle[2] = pnext;

        // qDebug() << "PolyIndexes:" << pprev << pcurr << pnext;

        int ntriangle[3];
        ntriangle[0] = polygon[ptriangle[0]];
        ntriangle[1] = polygon[ptriangle[1]];
        ntriangle[2] = polygon[ptriangle[2]];

        triangles.push_back(coords[ntriangle[2]]);
        triangles.push_back(coords[ntriangle[1]]);
        triangles.push_back(coords[ntriangle[0]]);

        // Remove the ear vertex from the clipped polygon.
        polygon.removeAt(pcurr);
        reflex.removeAt(pcurr);
        for (int i = 0; i < ears.size(); ++i) {
            if (ears[i] > pcurr) {
                ears[i] = ears[i] - 1;
            }
        }
        if (pnext > pcurr) {
            --pnext;
        }
        if (pprev > pcurr) {
            --pprev;
        }

        // if (polygon.empty() || reflex.empty()) {
        // qDebug() << "LEAVING FUNCTION NOW!";
        // break;
        //}

        // qDebug() << "PolyIndexes after removal:" << pprev << pcurr << pnext;

        // Removing an ear changes the configuration as follows:
        //  - If the neighbor is reflex, it might become convex and possibly an
        //  ear.
        //  - If the neighbor is convex, it remains convex and might become an
        //  ear.
        //  - If the neighbor is an ear, it might not stay an ear.
        for (int i = 0; i < 2; ++i) {
            int neighbour = ((i == 0) ? pprev : pnext);

            // qDebug() << "reflex test" << reflex.size() << neighbour;
            if (neighbour >= 0 && neighbour < reflex.size() &&
                reflex[neighbour] &&
                !IsReflexIndex(coords, polygon, neighbour)) {
                reflex[neighbour] = false;
                --reflexCount;
            }

            if (neighbour < 0 || neighbour >= reflex.size() ||
                !reflex[neighbour]) {
                bool isEar =
                    CheckEar(coords, polygon, reflex, reflexCount, neighbour);
                int earIndex = ears.indexOf(neighbour);
                bool wasEar = (earIndex != -1);
                if (isEar && !wasEar) {
                    ears.push_back(neighbour);
                } else if (!isEar && wasEar) {
                    ears.removeAt(earIndex);
                }
            }
        }
    }

    // qDebug() << "POLY:" << poly;

    return true;
}
