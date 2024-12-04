#include "bezierfit.h"

BezierFit::BezierFit() {}

/*
An Algorithm for Automatically Fitting Digitized Curves
by Philip J. Schneider
from "Graphics Gems", Academic Press, 1990
*/

/*
 *  main:
 *	Example of how to use the curve-fitting code.  Given an array
 *   of points and a tolerance (squared error between points and
 *	fitted curve), the algorithm will generate a piecewise
 *	cubic Bezier representation that approximates the points.
 *	When a cubic is generated, the routine "DrawBezierCurve"
 *	is called, which outputs the Bezier curve just created
 *	(arguments are the degree and the control points, respectively).
 *	Users will have to implement this function themselves
 *   ascii output, etc.
 *
 */
/*
main()
{
    static Point2 d[7] = {	//  Digitized points
    { 0.0, 0.0 },
    { 0.0, 0.5 },
    { 1.1, 1.4 },
    { 2.1, 1.6 },
    { 3.2, 1.1 },
    { 4.0, 0.2 },
    { 4.0, 0.0 },
    };
    double	error = 4.0;		//  Squared error
    FitCurve(d, 7, error);		//  Fit the Bezier curves
}
    */

/*
 *  FitCurve :
 *  	Fit a Bezier curve to a set of digitized points
 */
/*  Array of digitized points	*/
/*  Number of digitized points	*/
/*  User-defined error squared	*/
void BezierFit::FitCurve(const QList<QVector2D> & d, const double error,
                         QList<QVector2D> & fit_curve)
{
    // fit_curve.clear();

    // const QVector2D tHat1 = ComputeLeftTangent(d, 0); /*  Unit tangent
    // vectors at endpoints */ const QVector2D tHat2 = ComputeRightTangent(d,
    // d.size() - 1); FitCubic(d, 0, d.size() - 1, tHat1, tHat2, error,
    // fit_curve);

    // compute curve length
    // float len = 0.0f;
    // for (int i=1; i<d.size(); ++i) {
    //     len += (d[i] - d[i-1]).length();
    // }
    // qDebug() << d.size() << len;

    /*
    const int n_split = 40;
    for (int i=0; i<d.size()-1; i+=n_split) {

        const int i2 = qMin(d.size()-1, i+n_split);
        //qDebug() << i << i2 << d.size();

        const QVector2D tHat1 = ComputeLeftTangent(d, i); // Unit tangent
    vectors at endpoints const QVector2D tHat2 = ComputeRightTangent(d, i2);
        FitCubic(d, i, i2, tHat1, tHat2, error, fit_curve);

    }
    */

    const float smoothThreshold = 0.6f;
    const float curvatureThreshold = 0.6f;
    const int iterations = 3;

    // 1. Compute curvatures for all points on curve
    QVector<float> ks;
    ComputeCurvatures(d, ks);

    // 2.  Iterative passes to smooth the curve via local averaging
    QList<QVector2D> d_smooth = d;

    for (int j = 0; j < iterations; j++) {
        for (int i = 1; i < d.size() - 1; i++) {
            if (ks[i] < smoothThreshold) {
                d_smooth[i] =
                    (d_smooth[i - 1] + d_smooth[i] + d_smooth[i + 1]) / 3.0f;
            }
        }
    }

    // 3.  Recompute the curvature for the smoothed curve
    ComputeCurvatures(d_smooth, ks);

    // 4.  Compute list of indexes to points of high curvature
    //(these are used to segment the sketched curve at corners/discontinuities)
    QList<int> segment_indexes;
    segment_indexes.push_back(0);  // always include the first point
    for (int i = 2; i < ks.size() - 2; ++i) {
        if (ks[i] > curvatureThreshold) {
            segment_indexes.push_back(i);
            ++i;  // skip a point
        }
    }
    segment_indexes.push_back(d.size() - 1);  // always include the last point

    // qDebug() << "segment_indexes" << segment_indexes;

    // 5.  Fit cubic bezier splines to each segmented curve
    for (int i = 0; i < segment_indexes.size() - 1; ++i) {
        const int i1 = segment_indexes[i];
        const int i2 = segment_indexes[i + 1];

        const QVector2D tHat1 = ComputeLeftTangent(
            d_smooth, i1);  // Unit tangent vectors at endpoints
        const QVector2D tHat2 = ComputeRightTangent(d_smooth, i2);
        FitCubic(d_smooth, i1, i2, tHat1, tHat2, error, fit_curve);
    }
}

/*
 *  FitCubic :
 *  	Fit a Bezier curve to a (sub)set of digitized points
 */
// Point2	*d;			  Array of digitized points
// int		first, last;	 Indices of first and last pts in region
// Vector2	tHat1, tHat2;	 Unit tangent vectors at endpoints
// double	error;		  User-defined error squared
void BezierFit::FitCubic(const QList<QVector2D> & d, const int first,
                         const int last, const QVector2D & tHat1,
                         const QVector2D & tHat2, const double error,
                         QList<QVector2D> & fit_curve)

{
    QVector<QVector2D> bezCurve;  // Control points of fitted Bezier curve
    QVector<double> u;            //  Parameter values for point
    QVector<double> uPrime;       //  Improved parameter values
    double maxError;              //  Maximum fitting error
    int splitPoint;               //  Point to split point set at
    int nPts;                     //  Number of points in subset
    double iterationError;        // Error below which you try iterating
    int maxIterations = 4;        //  Max times to try iterating
    QVector2D tHatCenter;         // Unit tangent vector at splitPoint

    iterationError = error * error;
    nPts = last - first + 1;

    /*  Use heuristic if region only has two points in it */
    if (nPts == 2) {
        double dist =
            (d[last] - d[first])
                .length();  // double dist = V2DistanceBetween2Points(&d[last],
                            // &d[first]) / 3.0;

        bezCurve.resize(4);  // bezCurve = (Point2 *)malloc(4 * sizeof(Point2));
        bezCurve[0] = d[first];
        bezCurve[3] = d[last];
        bezCurve[1] = bezCurve[0] +
                      (tHat1 * dist);  // V2Add(&bezCurve[0], V2Scale(&tHat1,
                                       // dist), &bezCurve[1]);
        bezCurve[2] = bezCurve[3] +
                      (tHat2 * dist);  // V2Add(&bezCurve[3], V2Scale(&tHat2,
                                       // dist), &bezCurve[2]);
        // DrawBezierCurve(3, bezCurve);
        // qDebug() << "BezierFit::FitCubic - line 99" << first << last;
        fit_curve.push_back(bezCurve[0]);
        fit_curve.push_back(bezCurve[1]);
        fit_curve.push_back(bezCurve[2]);
        fit_curve.push_back(bezCurve[3]);
        // qDebug() << bezCurve;
        // free((void *)bezCurve);
        return;
    }

    /*  Parameterize points, and attempt to fit curve */
    ChordLengthParameterize(d, first, last, u);
    GenerateBezier(d, first, last, u, tHat1, tHat2, bezCurve);

    /*  Find max deviation of points to fitted curve */
    maxError = ComputeMaxError(d, first, last, bezCurve, u, splitPoint);
    if (maxError < error) {
        // TODO: do something here to provide curves
        // qDebug() << "BezierFit::FitCubic - line 113" << first << last;
        fit_curve.push_back(bezCurve[0]);
        fit_curve.push_back(bezCurve[1]);
        fit_curve.push_back(bezCurve[2]);
        fit_curve.push_back(bezCurve[3]);
        // qDebug() << bezCurve;
        // DrawBezierCurve(3, bezCurve);
        // free((void *)u);
        // free((void *)bezCurve);
        return;
    }

    /*  If error not too large, try some reparameterization  */
    /*  and iteration */
    if (maxError < iterationError) {
        for (int i = 0; i < maxIterations; i++) {
            Reparameterize(d, first, last, u, bezCurve, uPrime);
            // free((void *)bezCurve);

            GenerateBezier(d, first, last, uPrime, tHat1, tHat2, bezCurve);
            maxError =
                ComputeMaxError(d, first, last, bezCurve, uPrime, splitPoint);

            if (maxError < error) {
                // TODO: do something with this curve
                // qDebug() << "BezierFit::FitCubic - line 134" << first <<
                // last;
                fit_curve.push_back(bezCurve[0]);
                fit_curve.push_back(bezCurve[1]);
                fit_curve.push_back(bezCurve[2]);
                fit_curve.push_back(bezCurve[3]);
                // qDebug() << bezCurve;
                // DrawBezierCurve(3, bezCurve);
                // free((void *)u);
                // free((void *)bezCurve);
                // free((void *)uPrime);
                return;
            }

            // free((void *)u);
            u = uPrime;
        }
    }

    /* Fitting failed -- split at max error point and fit recursively */
    // free((void *)u);
    // free((void *)bezCurve);
    tHatCenter = ComputeCenterTangent(d, splitPoint);
    FitCubic(d, first, splitPoint, tHat1, tHatCenter, error, fit_curve);
    tHatCenter = -tHatCenter;  // V2Negate(&tHatCenter);
    FitCubic(d, splitPoint, last, tHatCenter, tHat2, error, fit_curve);
}

/*
 *  GenerateBezier :
 *  Use least-squares method to find Bezier control points for region.
 *
 */
void BezierFit::GenerateBezier(const QList<QVector2D> & d, const int first,
                               const int last, const QVector<double> & uPrime,
                               const QVector2D tHat1, const QVector2D tHat2,
                               QVector<QVector2D> & bezCurve)
//    Point2	*d;			//  Array of digitized points
//    int		first, last;		//  Indices defining region
//    double	*uPrime;		//  Parameter values for region
//    Vector2	tHat1, tHat2;	//  Unit tangents at endpoints
{
    // int 	i;

    int nPts;         /* Number of pts in sub-curve */
    double C[2][2];   /* Matrix C		*/
    double X[2];      /* Matrix X			*/
    double det_C0_C1, /* Determinants of matrices	*/
        det_C0_X, det_X_C1;
    double alpha_l, /* Alpha values, left and right	*/
        alpha_r;
    // QVector2D 	tmp;			/* Utility variable
    // */ BezierCurve	bezCurve;	/* RETURN bezier curve ctl pts	*/

    bezCurve.resize(4);  // bezCurve = (Point2 *)malloc(4 * sizeof(Point2));
    nPts = last - first + 1;

    // QVector2D A[nPts][2];	/* Precomputed rhs for eqn	*/
    QVector<QVector<QVector2D> > A(nPts, QVector<QVector2D>(2));

    /* Compute the A's	*/
    for (int i = 0; i < nPts; i++) {
        QVector2D v1, v2;
        /*
        v1 = tHat1;
        v2 = tHat2;
        V2Scale(&v1, B1(uPrime[i]));
        V2Scale(&v2, B2(uPrime[i]));
        */
        v1 = tHat1 * B1(uPrime[i]);
        v2 = tHat2 * B2(uPrime[i]);
        A[i][0] = v1;
        A[i][1] = v2;
    }

    /* Create the C and X matrices	*/
    C[0][0] = 0.0;
    C[0][1] = 0.0;
    C[1][0] = 0.0;
    C[1][1] = 0.0;
    X[0] = 0.0;
    X[1] = 0.0;

    for (int i = 0; i < nPts; i++) {
        C[0][0] += QVector2D::dotProduct(
            A[i][0], A[i][0]);  // V2Dot(&A[i][0], &A[i][0]);
        C[0][1] += QVector2D::dotProduct(
            A[i][0], A[i][1]);  // V2Dot(&A[i][0], &A[i][1]);
        /*					C[1][0] += V2Dot(&A[i][0],
         * &A[i][1]);*/
        C[1][0] = C[0][1];
        C[1][1] += QVector2D::dotProduct(
            A[i][1], A[i][1]);  // V2Dot(&A[i][1], &A[i][1]);

        /*
        tmp = V2SubII(d[first + i],
                V2AddII(
                    V2ScaleIII(d[first], B0(uPrime[i])),
                    V2AddII(
                        V2ScaleIII(d[first], B1(uPrime[i])),
                        V2AddII(
                            V2ScaleIII(d[last], B2(uPrime[i])),
                            V2ScaleIII(d[last], B3(uPrime[i]))))));
                            */

        const QVector2D p0 = d[first] * B0(uPrime[i]);
        const QVector2D p1 = d[first] * B1(uPrime[i]);
        const QVector2D p2 = d[last] * B2(uPrime[i]);
        const QVector2D p3 = d[last] * B3(uPrime[i]);

        const QVector2D tmp = d[first + i] - (p0 + p1 + p2 + p3);

        X[0] += QVector2D::dotProduct(A[i][0], tmp);  // V2Dot(&A[i][0], &tmp);
        X[1] += QVector2D::dotProduct(A[i][1], tmp);  // V2Dot(&A[i][1], &tmp);
    }

    /* Compute the determinants of C and X	*/
    det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
    det_C0_X = C[0][0] * X[1] - C[1][0] * X[0];
    det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];

    /* Finally, derive alpha values	*/
    alpha_l = (det_C0_C1 == 0) ? 0.0 : det_X_C1 / det_C0_C1;
    alpha_r = (det_C0_C1 == 0) ? 0.0 : det_C0_X / det_C0_C1;

    /* If alpha negative, use the Wu/Barsky heuristic (see text) */
    /* (if alpha is 0, you get coincident control points that lead to
     * divide by zero in any subsequent NewtonRaphsonRootFind() call. */
    double segLength =
        (d[last] - d[first])
            .length();  // V2DistanceBetween2Points(&d[last], &d[first]);
    double epsilon = 1.0e-6 * segLength;
    if (alpha_l < epsilon || alpha_r < epsilon) {
        /* fall back on standard (probably inaccurate) formula, and subdivide
         * further if needed. */
        double dist = segLength / 3.0;
        bezCurve[0] = d[first];
        bezCurve[3] = d[last];
        bezCurve[1] = bezCurve[0] +
                      (tHat1 * dist);  // V2Add(&bezCurve[0], V2Scale(&tHat1,
                                       // dist), &bezCurve[1]);
        bezCurve[2] = bezCurve[3] +
                      (tHat2 * dist);  // V2Add(&bezCurve[3], V2Scale(&tHat2,
                                       // dist), &bezCurve[2]);
        return;                        // return (bezCurve);
    }

    /*  First and last control points of the Bezier curve are */
    /*  positioned exactly at the first and last data points */
    /*  Control points 1 and 2 are positioned an alpha distance out */
    /*  on the tangent vectors, left and right, respectively */
    bezCurve[0] = d[first];
    bezCurve[3] = d[last];
    bezCurve[1] =
        bezCurve[0] + (tHat1 * alpha_l);  // V2Add(&bezCurve[0], V2Scale(&tHat1,
                                          // alpha_l), &bezCurve[1]);
    bezCurve[2] =
        bezCurve[3] + (tHat2 * alpha_r);  // V2Add(&bezCurve[3], V2Scale(&tHat2,
                                          // alpha_r), &bezCurve[2]);
    return;                               // return (bezCurve);
}

/*
 *  Reparameterize:
 *	Given set of points and their parameterization, try to find
 *   a better parameterization.
 *
 */
void BezierFit::Reparameterize(const QList<QVector2D> & d, const int first,
                               const int last, const QVector<double> & u,
                               const QVector<QVector2D> & bezCurve,
                               QVector<double> & uPrime)
//    Point2	*d;			/*  Array of digitized points	*/
//    int		first, last;		/*  Indices defining region
//    */ double	*u;			/*  Current parameter values	*/
//    BezierCurve	bezCurve;	/*  Current fitted curve	*/
{
    int nPts = last - first + 1;
    // double	*uPrime;		/*  New parameter values	*/

    // uPrime = (double *)malloc(nPts * sizeof(double));
    uPrime.resize(nPts);
    for (int i = first; i <= last; i++) {
        uPrime[i - first] = NewtonRaphsonRootFind(bezCurve, d[i], u[i - first]);
    }
    // return (uPrime);
}

/*
 *  NewtonRaphsonRootFind :
 *	Use Newton-Raphson iteration to find better root.
 */
double BezierFit::NewtonRaphsonRootFind(const QVector<QVector2D> & Q,
                                        const QVector2D & P, const double u)
//    BezierCurve	Q;			/*  Current fitted curve
//    */ Point2 		P;		/*  Digitized point
//    */ double 		u;		/*  Parameter value for "P"
//    */
{
    double numerator, denominator;
    // QVector2D 		Q1[3], Q2[2];	/*  Q' and Q'' */
    QVector<QVector2D> Q1;
    Q1.resize(3);
    QVector<QVector2D> Q2; /*  Q' and Q''			*/
    Q2.resize(2);
    QVector2D Q_u, Q1_u, Q2_u; /*u evaluated at Q, Q', & Q''	*/
    double uPrime;             /*  Improved u			*/
    // int 		i;

    /* Compute Q(u)	*/
    Q_u = BezierII(3, Q, u);

    /* Generate control vertices for Q'	*/
    for (int i = 0; i <= 2; i++) {
        // Q1[i].x = (Q[i+1].x - Q[i].x) * 3.0;
        // Q1[i].y = (Q[i+1].y - Q[i].y) * 3.0;
        Q1[i] = QVector2D((Q[i + 1].x() - Q[i].x()) * 3.0,
                          (Q[i + 1].y() - Q[i].y()) * 3.0);
    }

    /* Generate control vertices for Q'' */
    for (int i = 0; i <= 1; i++) {
        // Q2[i].x = (Q1[i+1].x - Q1[i].x) * 2.0;
        // Q2[i].y = (Q1[i+1].y - Q1[i].y) * 2.0;
        Q2[i] = QVector2D((Q1[i + 1].x() - Q1[i].x()) * 2.0,
                          (Q1[i + 1].y() - Q1[i].y()) * 2.0);
    }

    /* Compute Q'(u) and Q''(u)	*/
    Q1_u = BezierII(2, Q1, u);
    Q2_u = BezierII(1, Q2, u);

    /* Compute f(u)/f'(u) */
    numerator = (Q_u.x() - P.x()) * (Q1_u.x()) + (Q_u.y() - P.y()) * (Q1_u.y());
    denominator = (Q1_u.x()) * (Q1_u.x()) + (Q1_u.y()) * (Q1_u.y()) +
                  (Q_u.x() - P.x()) * (Q2_u.x()) +
                  (Q_u.y() - P.y()) * (Q2_u.y());

    if (denominator == 0.0f) {
        return u;
    }

    /* u = u - f(u)/f'(u) */
    uPrime = u - (numerator / denominator);
    return (uPrime);
}

/*
 *  Bezier :
 *  	Evaluate a Bezier curve at a particular parameter value
 *
 */
QVector2D BezierFit::BezierII(const int degree, const QVector<QVector2D> & V,
                              const double t)
//    int		degree;		/* The degree of the bezier curve
//    */ Point2 	*V;		/* Array of control points
//    */ double 	t;		/* Parametric value to find point for
//    */
{
    // int 	i, j;
    QVector2D Q;              /* Point on curve at parameter t	*/
    QVector<QVector2D> Vtemp; /* Local copy of control points		*/

    /* Copy array	*/
    Vtemp.resize(degree + 1);  // Vtemp = (Point2 *)malloc((unsigned)((degree+1)
                               // * sizeof (Point2)));

    for (int i = 0; i <= degree; i++) {
        Vtemp[i] = V[i];
    }

    /* Triangle computation	*/
    for (int i = 1; i <= degree; i++) {
        for (int j = 0; j <= degree - i; j++) {
            // Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j+1].x;
            // Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j+1].y;
            Vtemp[j] =
                QVector2D((1.0 - t) * Vtemp[j].x() + t * Vtemp[j + 1].x(),
                          (1.0 - t) * Vtemp[j].y() + t * Vtemp[j + 1].y());
        }
    }

    Q = Vtemp[0];
    // free((void *)Vtemp);
    return Q;
}

/*
 *  B0, B1, B2, B3 :
 *	Bezier multipliers
 */
double BezierFit::B0(const double u)
{
    double tmp = 1.0 - u;
    return (tmp * tmp * tmp);
}

double BezierFit::B1(const double u)
{
    double tmp = 1.0 - u;
    return (3 * u * (tmp * tmp));
}

double BezierFit::B2(const double u)
{
    double tmp = 1.0 - u;
    return (3 * u * u * tmp);
}

double BezierFit::B3(const double u) { return (u * u * u); }

/*
 * ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent :
 *Approximate unit tangents at endpoints and "center" of digitized curve
 */
QVector2D BezierFit::ComputeLeftTangent(const QList<QVector2D> & d,
                                        const int end)
//    Point2	*d;			//  Digitized points
//    int		end;		// Index to "left" end of region
{
    QVector2D tHat1;
    tHat1 = d[end + 1] - d[end];  // V2SubII(d[end+1], d[end]);
    tHat1.normalize();            // tHat1 = *V2Normalize(&tHat1);
    return tHat1;
}

QVector2D BezierFit::ComputeRightTangent(const QList<QVector2D> & d,
                                         const int end)
//    Point2	*d;			//  Digitized points
//    int		end;		// Index to "left" end of region
{
    QVector2D tHat2;
    tHat2 = d[end - 1] - d[end];  // V2SubII(d[end-1], d[end]);
    tHat2.normalize();            // tHat2 = *V2Normalize(&tHat2);
    return tHat2;
}

QVector2D BezierFit::ComputeCenterTangent(const QList<QVector2D> & d,
                                          const int center)
//    Point2	*d;			//  Digitized points
//    int		center;		//  Index to point inside region
{
    QVector2D V1, V2, tHatCenter;

    V1 = d[center - 1] - d[center];  // V2SubII(d[center-1], d[center]);
    V2 = d[center] - d[center + 1];  // V2SubII(d[center], d[center+1]);
    // tHatCenter.x = (V1.x + V2.x)/2.0;
    // tHatCenter.y = (V1.y + V2.y)/2.0;
    tHatCenter = QVector2D((V1.x() + V2.x()) / 2.0, (V1.y() + V2.y()) / 2.0);
    tHatCenter.normalize();  // tHatCenter = *V2Normalize(&tHatCenter);
    return tHatCenter;
}

/*
 *  ChordLengthParameterize :
 *	Assign parameter values to digitized points
 *	using relative distances between points.
 */
void BezierFit::ChordLengthParameterize(const QList<QVector2D> & d,
                                        const int first, const int last,
                                        QVector<double> & u)
//    Point2	*d;			// Array of digitized points
//    int		first, last;		//  Indices defining region
{
    int i;
    // double	*u;			/*  Parameterization		*/

    // u = (double *)malloc((unsigned)(last-first+1) * sizeof(double));
    u.resize(last - first + 1);

    u[0] = 0.0;
    for (i = first + 1; i <= last; i++) {
        u[i - first] =
            u[i - first - 1] +
            (d[i] - d[i - 1])
                .length();  // V2DistanceBetween2Points(&d[i], &d[i-1]);
    }

    for (i = first + 1; i <= last; i++) {
        u[i - first] = u[i - first] / u[last - first];
    }

    // return(u);
}

/*
 *  ComputeMaxError :
 *	Find the maximum squared distance of digitized points
 *	to fitted curve.
 */
double BezierFit::ComputeMaxError(const QList<QVector2D> & d, const int first,
                                  const int last,
                                  const QVector<QVector2D> & bezCurve,
                                  const QVector<double> & u, int & splitPoint)
//    Point2	*d;			//  Array of digitized points
//    int		first, last;		//  Indices defining region
//    BezierCurve	bezCurve;		//  Fitted Bezier curve
//    double	*u;			//  Parameterization of points
//    int		*splitPoint;		//  Point of maximum error
{
    double maxDist; /*  Maximum error		*/
    double dist;    /*  Current error		*/
    QVector2D P;    /*  Point on curve		*/
    QVector2D v;    /*  Vector from point to curve	*/

    splitPoint = (last - first + 1) / 2;
    maxDist = 0.0;
    for (int i = first + 1; i < last; i++) {
        P = BezierII(3, bezCurve, u[i - first]);
        v = P - d[i];              // V2SubII(P, d[i]);
        dist = v.lengthSquared();  // V2SquaredLength(&v);
        if (dist >= maxDist) {
            maxDist = dist;
            splitPoint = i;
        }
    }
    return (maxDist);
}
/*
static Vector2 V2AddII(a, b)
Vector2 a, b;
{
Vector2	c;
c.x = a.x + b.x;  c.y = a.y + b.y;
return (c);
}
static Vector2 V2ScaleIII(v, s)
Vector2	v;
double	s;
{
Vector2 result;
result.x = v.x * s; result.y = v.y * s;
return (result);
}

static Vector2 V2SubII(a, b)
Vector2	a, b;
{
Vector2	c;
c.x = a.x - b.x; c.y = a.y - b.y;
return (c);
}
*/

void BezierFit::ComputeCurvatures(const QList<QVector2D> & d,
                                  QVector<float> & ks)
{
    ks.clear();
    ks.reserve(d.size());
    for (int i = 0; i < d.size(); i++) {
        ks.push_back(CurvatureAtPoint(d, i));
    }
}

float BezierFit::CurvatureAtPoint(const QList<QVector2D> & d, const int index)
{
    if (index <= 0 || index >= d.size() - 1) {
        return 0;
    }

    float a = 1.0f;
    float b = 1.0f;
    float c = ((d[index + 1] - d[index]).normalized() -
               (d[index - 1] - d[index]).normalized())
                  .length();

    if (c > 2.0f - 0.00001f) {
        return 0;
    }

    float result = 1.0f / ((a * b * c) / sqrtf((a + b + c) * (b + c - a) *
                                               (c + a - b) * (a + b - c)));

    return result;
}
