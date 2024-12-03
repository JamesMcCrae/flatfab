#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include <QtGui>
#include <QtOpenGL>
#include <cfloat>

#include "glutils.h"

struct BezierCurvePoint {
    QVector2D point;
    int segment; //index to the POINT, so this will be some multile of 3
    float t; //t =[0,1]

    bool operator< (const BezierCurvePoint & p) const {
        return segment < p.segment || (segment == p.segment && t < p.t);
    }

    bool operator== (const BezierCurvePoint & p) const {
        return (segment == p.segment) && (t == p.t);
    }
};

class BezierCurve
{

public:

    BezierCurve();

    void AddPoint(const QVector2D & v);
    void SetPoints(const QList <QVector2D> & ps);    
    void SetPointsFromPolyline(const QList <QVector2D> & polyline);
    void SetPoint(const int index, const QVector2D & v);
    void TranslateControlPoints(const QVector2D & v);
    void RotateControlPoints(const float rotate_radians);
    void ScaleControlPoints(const float x, const float y);
    int GetNumControlPoints() const;
    QVector2D GetControlPointCentroid() const;
    const QList <QVector2D> & Points() const;
    const QVector2D Point(const int i) const;
    void ClearPoints();

    void LoadFromSVGData(const QString & s);

    int GetPreviousSegmentIndex(const int i) const;
    int GetNextSegmentIndex(const int i) const;

    void GetPointsTangentsAlongCurve(const int num, QVector <QVector2D> & points, QVector <QVector2D> & tangents) const;

    void SetClosed(bool b);
    bool IsClosed();
    void Draw2D();

    void SetSamplesPerSegment(const int i);
    void UpdateSamples();
    const QList <QVector2D> & Samples() const;    

    void SelectPoint(const QVector2D & p, const float select_size);
    void SetSelectedPoint(const int index);
    int SelectedPoint() const;
    void MoveSelectedPoint(const QVector2D & p, const bool keep_g1, const bool equal_lengths = false);
    void DeletePoint(const int segment_index);
    void DeleteSelectedPoint();
    void InsertPoint(const int segment_index);
    void InsertSelectedPoint();
    void UnselectPoint();

    void GetLineIntersections(const QVector2D & line_p, const QVector2D & line_d, QList <BezierCurvePoint> & intersects);
    void GetLineSegmentIntersections(const QVector2D & line_p1, const QVector2D & line_p2, QList <BezierCurvePoint> & intersects);
    bool IsPointInside(const QVector2D & p);
    void SplitAlongLine(const QVector2D & split_p, const QVector2D & split_d, QList <BezierCurve> & curves);

    void SubdivideLongestSegment();    

    //used for curve interpolation
    void GetCurvePointCorrespondence(const BezierCurve & other, int & corresp_offset, bool & corresp_forward) const;
    void InterpolateCurvePointCorrespondence(const BezierCurve & other, const int corresp_offset, const bool corresp_forward, const float t, BezierCurve & interp_curve) const;

    float Length() const;
    void Reverse();
    void GetSubCurve(const BezierCurvePoint & b1, const BezierCurvePoint & b2, BezierCurve & sub_curve);


private:

    void SamplePointsForSegment(const int i, QVector <QVector2D> & sample_segment) const;

    void Subdivide(const int i, const float t);
    static void Subdivide(const BezierCurve & c, const int i, const float t, QList <QVector2D> & subdiv_pts);
    static void SamplePoint(const QVector2D & p1, const QVector2D & t1, const QVector2D & t2, const QVector2D & p2, const float t, QVector2D & point);
    static void SampleTangent(const QVector2D & p1, const QVector2D & t1, const QVector2D & t2, const QVector2D & p2, const float t, QVector2D & tangent);

    QList <QVector2D> pts; //points
    QList <QVector2D> samples; //points (with segment intermediate samples)

    bool closed;

    //what do segments have
    //initial point
    //   "    tangent
    //end     point
    //   "    tangent

    int samples_per_segment;
    int tests_sample_rate;
    int selected;

};

#endif // BEZIERCURVE_H
