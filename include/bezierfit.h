#ifndef BEZIERFIT_H
#define BEZIERFIT_H

#include <QtGui>
#include <QVector>

class BezierFit
{

public:

    BezierFit();

    static void FitCurve(const QList <QVector2D> & d, const double error, QList <QVector2D> & fit_curve);
    static void FitCubic(const QList <QVector2D> & d, const int first, const int last, const QVector2D & tHat1, const QVector2D & tHat2, const double error, QList <QVector2D> & fit_curve);
    static void Reparameterize(const QList <QVector2D> & d, const int first, const int last, const QVector <double> & u, const QVector <QVector2D> & bezCurve, QVector <double> & uPrime);
    static double NewtonRaphsonRootFind(const QVector <QVector2D> & Q, const QVector2D & P, const double u);
    static QVector2D BezierII(const int degree, const QVector <QVector2D> & V, const double t);
    static QVector2D ComputeLeftTangent(const QList <QVector2D> & d, const int end);
    static QVector2D ComputeRightTangent(const QList <QVector2D> & d, const int end);
    static QVector2D ComputeCenterTangent(const QList <QVector2D> & d, const int center);
    static double ComputeMaxError(const QList <QVector2D> & d, const int first, const int last, const QVector <QVector2D> & bezCurve, const QVector <double> & u, int & splitPoint);
    static void ChordLengthParameterize(const QList <QVector2D> & d, const int first, const int last, QVector <double> & u);
    static void GenerateBezier(const QList <QVector2D> & d, const int first, const int last, const QVector <double> & uPrime, const QVector2D tHat1, const QVector2D tHat2, QVector <QVector2D> & bezCurve);

    static void ComputeCurvatures(const QList <QVector2D> & d, QVector <float> & curvs);
    static float CurvatureAtPoint(const QList <QVector2D> & d, const int index);

    static double B0(const double u);
    static double B1(const double u);
    static double B2(const double u);
    static double B3(const double u);

};

#endif // BEZIERFIT_H
