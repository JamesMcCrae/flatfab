#include "triangulate.h"

//static const float EPSILON=0.000001f;
static const float EPSILON=0.00000001f;

float Triangulate::Area(const QList <QVector2D> & contour)
{

    int n = contour.size();

    float A = 0.0f;

    for(int p=n-1,q=0; q<n; p=q++) {

        A += contour[p].x() * contour[q].y() - contour[q].x() * contour[p].y();

    }

    return A * 0.5f;

}

/*
     InsideTriangle decides if a point P is Inside of the triangle
     defined by A, B, C.
   */
bool Triangulate::InsideTriangle(float Ax, float Ay,
                                 float Bx, float By,
                                 float Cx, float Cy,
                                 float Px, float Py)

{

    float ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
    float cCROSSap, bCROSScp, aCROSSbp;

    ax = Cx - Bx;  ay = Cy - By;
    bx = Ax - Cx;  by = Ay - Cy;
    cx = Bx - Ax;  cy = By - Ay;
    apx= Px - Ax;  apy= Py - Ay;
    bpx= Px - Bx;  bpy= Py - By;
    cpx= Px - Cx;  cpy= Py - Cy;

    aCROSSbp = ax*bpy - ay*bpx;
    cCROSSap = cx*apy - cy*apx;
    bCROSScp = bx*cpy - by*cpx;

    return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
};

bool Triangulate::Snip(const QList <QVector2D> & contour, int u, int v, int w, int n, const QVector <int> & V)
{
    int p;
    float Ax, Ay, Bx, By, Cx, Cy, Px, Py;

    Ax = contour[V[u]].x();
    Ay = contour[V[u]].y();

    Bx = contour[V[v]].x();
    By = contour[V[v]].y();

    Cx = contour[V[w]].x();
    Cy = contour[V[w]].y();

    if ( EPSILON > (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) {
        return false;
    }

    for (p=0; p<n; p++) {

        if( (p == u) || (p == v) || (p == w) ) {
            continue;
        }

        Px = contour[V[p]].x();
        Py = contour[V[p]].y();

        if (InsideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) {
            return false;
        }
    }

    return true;
}

bool Triangulate::Process(const QList <QVector2D> & contour, QList <QVector2D> & result)
{

    result.clear();

    /* allocate and initialize list of Vertices in polygon */
    int n = contour.size();
    if ( n < 3 ) return false;

    QVector <int> V(n);

    /* we want a counter-clockwise polygon in V */

    if ( 0.0f < Area(contour) ) {
        for (int v=0; v<n; v++) {
            V[v] = v;
        }
    }
    else {
        for(int v=0; v<n; v++) {
            V[v] = (n-1)-v;
        }
    }

    int nv = n;

    /*  remove nv-2 Vertices, creating 1 triangle every time */
    int count = 2*nv;   /* error detection */

    for(int m=0, v=nv-1; nv>2; )
    {
        /* if we loop, it is probably a non-simple polygon */
        if (0 >= (count--))
        {
            //** Triangulate: ERROR - probably bad polygon!
            //qDebug() << "Triangulate::Process() - ERROR - probably bad polygon!";
            return false;
        }

        /* three consecutive vertices in current polygon, <u,v,w> */
        int u = v;
        if (nv <= u) {
            u = 0;     /* previous */
        }

        v = u+1;
        if (nv <= v) {
            v = 0;     /* new v    */
        }

        int w = v+1;
        if (nv <= w) {
            w = 0;     /* next     */
        }

        if (Snip(contour, u, v, w, nv, V)) {
            int a,b,c,s,t;

            /* true names of the vertices */
            a = V[u]; b = V[v]; c = V[w];

            /* output Triangle */
            result.push_back( contour[a] );
            result.push_back( contour[b] );
            result.push_back( contour[c] );

            m++;

            /* remove v from remaining polygon */
            for(s=v,t=v+1; t<nv; s++,t++) {
                V[s] = V[t];
            }

            nv--;

            /* reset error detection counter */
            count = 2*nv;

        }
    }

    return true;
}

bool Triangulate::Process(const QList <QList <QVector2D> > & contours, QList <QVector2D> & result)
{

    /*
    float smallest = FLT_MAX;
    for (int i=0; i<contours.size(); ++i) {
        for (int j=0; j<contours[i].size(); ++j) {

            float each = (contours[i][j] - contours[i][(j+1)%contours[i].size()]).length();

            if (each < smallest) {
                smallest = each;
            }

        }
    }

    qDebug() << "Triangulate::Process smallest length" << smallest;
    */

    result.clear();

    //since there's multiple curves (ones after the first are holes)
    //we use the triangle library rather than ear method, which can handle PSLGs with holes

    //
    // Do 2D PSLG triangulation on curves formed by intersecting plane with mesh
    //
    triangulateio in, out;

    in.numberofpoints = 0;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    int cIndex=0;

    for (int c=0; c<contours.size(); ++c) {

        in.numberofsegments += contours[c].size();

        if (c > 0) {
            ++in.numberofholes;
        }
    }

    // Define input points.
    in.numberofpoints = in.numberofsegments;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;

    for (int c=0; c<contours.size(); ++c) {

        for (int i=0; i<contours[c].size(); ++i) {

            in.pointlist[cIndex*2]=contours[c][i].x() * 1000.0f;
            in.pointlist[cIndex*2+1]=contours[c][i].y() * 1000.0f;

            //qDebug() << "Chain" << c << "projected point" << i << "is" << projchains[c][i].x << projchains[c][i].y;
            cIndex++;

        }
    }


    //define PSLG line segments
    in.segmentlist=(int *)malloc(in.numberofsegments*2*sizeof(int));
    in.segmentmarkerlist=NULL;

    unsigned int origIndex;

    cIndex=0;

    for (int c=0; c<contours.size(); ++c) {

        origIndex=cIndex;

        for (int i=0; i<contours[c].size(); ++i) {

            in.segmentlist[cIndex*2]=cIndex;

            if ( i+1 < contours[c].size()) {
                in.segmentlist[cIndex*2+1]=cIndex+1;
            }
            else {
                in.segmentlist[cIndex*2+1]=origIndex;
            }

            ++cIndex;
        }
    }

    //define midpoints of those contours with holes, as being holes
    in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));

    cIndex=0;
    for (int c=1; c<contours.size(); ++c) {

        //const QVector2D each_centroid = Centroid(contours[c]);

        const QList <QVector2D> & pts = contours[c];

        //rather than use centroid, which doesnt work all the time, use a ray trace and use midpoint of first sorted segments
        //QVector2D ray_dir = pts[pts.size()/2] - pts.first();
        QVector2D ray_dir(1, 0);
        QVector2D ray_p = (pts[1] + pts[0]) * 0.5f;

        QList <QVector2D> sorted_isecs;
        GLutils::GetSortedIntersectionPoints(contours[c], ray_dir, ray_p, sorted_isecs);

        //qDebug() << "sorted_isecs" << sorted_isecs;

        if (sorted_isecs.size() >= 2) {
            QVector2D hole_pt = (sorted_isecs[1] + sorted_isecs[0]) * 0.5f;

            in.holelist[cIndex*2]=hole_pt.x() * 1000.0f;
            in.holelist[cIndex*2+1]=hole_pt.y() * 1000.0f;
            ++cIndex;
        }

    }

    in.numberofregions = 0;

    // Make necessary initializations so that Triangle can return a
    //   triangulation in `mid' and a voronoi diagram in `vorout'.

    out.pointlist = (REAL *) NULL;            // Not needed if -N switch used.
    // Not needed if -N switch used or number of point attributes is zero:
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; // Not needed if -N or -B switch used.
    out.trianglelist = (int *) NULL;          // Not needed if -E switch used.
    //Not needed if -E switch used or number of triangle attributes is zero:
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;         // Needed only if -n switch used.
    //Needed only if segments are output (-p or -c) and -P not used:
    out.segmentlist = (int *) NULL;
    // Needed only if segments are output (-p or -c) and -P and -B not used:
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;             // Needed only if -e switch used.
    out.edgemarkerlist = (int *) NULL;   // Needed if -e used and -B not used.

    if (in.numberofsegments>=3)
        triangulate((char *)"pzQY", &in, &out, NULL); // not delaunay, fewer tris
    else
        out.numberoftriangles=0;

    //
    // Map every 2D planar point back into 3-space
    //

    for (int i=0; i<out.numberoftriangles; ++i) {

        for (int j=0; j<3; ++j) {

            float x = out.pointlist[out.trianglelist[i*3+j]*2] / 1000.0f;
            float y = out.pointlist[out.trianglelist[i*3+j]*2+1] / 1000.0f;

            result.push_back(QVector2D(x, y));

        }

    }

    free(in.pointlist);
    free(in.segmentlist);
    free(in.holelist);

    if (out.numberoftriangles>0) {
        free(out.pointlist);
        free(out.trianglelist);
    }

    return true;

}

QVector2D Triangulate::Centroid(const QList <QVector2D> & contour)
{

    QVector2D centroid(0, 0);

    for (int i=0; i<contour.size(); ++i) {
        centroid += contour[i];
    }

    centroid /= float(contour.size());

    return centroid;

}
