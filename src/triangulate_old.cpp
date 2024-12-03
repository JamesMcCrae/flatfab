#include "triangulate.h"

static const float EPSILON=0.000001f;

float Triangulate::Area(const QList <vec2> & contour)
{

    int n = contour.size();

    float A = 0.0f;

    for(int p=n-1,q=0; q<n; p=q++) {

        A+= contour[p][0]*contour[q][1] - contour[q][0]*contour[p][1];

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

bool Triangulate::Snip(const QList <vec2> & contour, int u, int v, int w, int n, const QVector <int> & V)
{
    int p;
    float Ax, Ay, Bx, By, Cx, Cy, Px, Py;

    Ax = contour[V[u]][0];
    Ay = contour[V[u]][1];

    Bx = contour[V[v]][0];
    By = contour[V[v]][1];

    Cx = contour[V[w]][0];
    Cy = contour[V[w]][1];

    if ( EPSILON > (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) return false;

    for (p=0;p<n;p++)
    {
        if( (p == u) || (p == v) || (p == w) ) continue;
        Px = contour[V[p]][0];
        Py = contour[V[p]][1];
        if (InsideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
    }

    return true;
}

bool Triangulate::Process(const QList <vec2> & contour, QList <vec2> & result)
{
    /* allocate and initialize list of Vertices in polygon */

    int n = contour.size();
    if ( n < 3 ) return false;

    QVector <int> V(n);

    /* we want a counter-clockwise polygon in V */

    if ( 0.0f < Area(contour) )
        for (int v=0; v<n; v++) V[v] = v;
    else
        for(int v=0; v<n; v++) V[v] = (n-1)-v;

    int nv = n;

    /*  remove nv-2 Vertices, creating 1 triangle every time */
    int count = 2*nv;   /* error detection */

    for(int m=0, v=nv-1; nv>2; )
    {
        /* if we loop, it is probably a non-simple polygon */
        if (0 >= (count--))
        {
            //** Triangulate: ERROR - probably bad polygon!
            return false;
        }

        /* three consecutive vertices in current polygon, <u,v,w> */
        int u = v  ; if (nv <= u) u = 0;     /* previous */
        v = u+1; if (nv <= v) v = 0;     /* new v    */
        int w = v+1; if (nv <= w) w = 0;     /* next     */

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
            for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

            /* reset error detection counter */
            count = 2*nv;
        }
    }

    return true;
}
