#include "beziercurve.h"

BezierCurve::BezierCurve() :
    closed(true),
    samples_per_segment(15),
    tests_sample_rate(200),
    selected(-1)
{
}

void BezierCurve::AddPoint(const QVector2D & v)
{
    pts.push_back(v);
}

void BezierCurve::SetPoints(const QList <QVector2D> & ps)
{
    pts = ps;
}

void BezierCurve::SetPointsFromPolyline(const QList <QVector2D> & polyline)
{

    pts.clear();

    for (int i=0; i<polyline.size()-1; ++i) {

        if (i == 0) {
            pts.push_back(polyline[i]);
        }
        pts.push_back(polyline[i] * 0.6666f + polyline[i+1] * 0.3333f);
        pts.push_back(polyline[i] * 0.3333f + polyline[i+1] * 0.6666f);
        pts.push_back(polyline[i+1]);

    }

}

void BezierCurve::SetPoint(const int index, const QVector2D & v)
{
    if (index < 0 || index >= pts.size()) {
        qDebug() << "BezierCurve::SetPoint - Warning, index"
                 << index << "out of range";
        return;
    }
    pts[index] = v;
}

void BezierCurve::TranslateControlPoints(const QVector2D & v)
{

    for (int i=0; i<pts.size(); ++i) {
        pts[i] += v;
    }

}

void BezierCurve::RotateControlPoints(const float rotate_radians)
{

    const float cos_rad = cosf(rotate_radians);
    const float sin_rad = sinf(rotate_radians);

    for (int i=0; i<pts.size(); ++i) {

        const float pt_x = cos_rad*pts[i].x() - sin_rad*pts[i].y();
        const float pt_y = sin_rad*pts[i].x() + cos_rad*pts[i].y();

        pts[i] = QVector2D(pt_x, pt_y);
    }
}

void BezierCurve::ScaleControlPoints(const float x, const float y)
{    
    for (int i=0; i<pts.size(); ++i) {
        pts[i].setX( pts[i].x() * x);
        pts[i].setY( pts[i].y() * y);
    }
}

int BezierCurve::GetNumControlPoints() const
{
    return pts.size();
}

QVector2D BezierCurve::GetControlPointCentroid() const
{
    QVector2D centroid(0, 0);
    for (int i=0; i<pts.size(); ++i) {
        centroid += pts[i];
    }
    centroid *= 1.0f / float(pts.size());
    return centroid;
}

void BezierCurve::ClearPoints()
{
    pts.clear();
    samples.clear();
}

void BezierCurve::LoadFromSVGData(const QString & s)
{

    pts.clear();
    samples.clear();

    QString s_filtered = s;
    s_filtered.replace(QString(","), QString(" "));

    QStringList path_data_tokens = s_filtered.split(" ");
    QVector2D current_pos(0, 0);
    QString last_command("");

    while (path_data_tokens.size() > 0) {

        //examine next token
        QString token = path_data_tokens.first();
        //qDebug() << "Parsing token" << path_data_tokens.size() << ":" << token;

        //check if token is a command token
        if (QString::compare(token, "M") == 0 || QString::compare(token, "m") == 0 || //M - "moveto" (absolute), m (relative)
            QString::compare(token, "Z") == 0 || QString::compare(token, "z") == 0 || //Z - close the curve
            QString::compare(token, "C") == 0 || QString::compare(token, "c") == 0 || //C - cubic Bezier "curveto" (absolute), c (relative)
            QString::compare(token, "S") == 0 || QString::compare(token, "s") == 0 || //S - cubic Bezier "shorthand curveto" (absolute), s (relative)
            QString::compare(token, "L") == 0 || QString::compare(token, "l") == 0) { //L - line (absolute), l (relative)

            last_command = token;
            path_data_tokens.pop_front(); //discard the command token

        }

        //parse the coordinates that come after a command token
        if (QString::compare(last_command, "M") == 0 || QString::compare(last_command, "m") == 0) { //M - "moveto" (absolute), m (relative)

            QString tok1 = path_data_tokens.first(); path_data_tokens.pop_front();
            QString tok2 = path_data_tokens.first(); path_data_tokens.pop_front();

            if (QString::compare(last_command, "M") == 0) {
                current_pos = QVector2D(tok1.toFloat(), tok2.toFloat());
            }
            else {
                current_pos = current_pos + QVector2D(tok1.toFloat(), tok2.toFloat());
            }

        }
        else if (QString::compare(last_command, "Z") == 0 || QString::compare(last_command, "z") == 0) { //Z/z - close the curve
            ;
        }
        else if (QString::compare(last_command, "C") == 0 || QString::compare(last_command, "c") == 0 || //C - cubic Bezier "curveto" (absolute), c (relative)
                 QString::compare(last_command, "S") == 0 || QString::compare(last_command, "s") == 0) { //S - cubic Bezier "shorthand curveto" (absolute), s (relative)

            QVector2D tan1;
            QVector2D tan2;
            QVector2D pos;

            if (QString::compare(last_command, "C") == 0 || QString::compare(last_command, "c") == 0) {

                //expect 3 params total
                QString tok1_tan1 = path_data_tokens.first(); path_data_tokens.pop_front();
                QString tok2_tan1 = path_data_tokens.first(); path_data_tokens.pop_front();

                QString tok1_tan2 = path_data_tokens.first(); path_data_tokens.pop_front();
                QString tok2_tan2 = path_data_tokens.first(); path_data_tokens.pop_front();

                QString tok1 = path_data_tokens.first(); path_data_tokens.pop_front();
                QString tok2 = path_data_tokens.first(); path_data_tokens.pop_front();

                tan1 = QVector2D(tok1_tan1.toFloat(), tok2_tan1.toFloat());
                tan2 = QVector2D(tok1_tan2.toFloat(), tok2_tan2.toFloat());
                pos = QVector2D(tok1.toFloat(), tok2.toFloat());

            }
            else {

                //expect only 2 params
                QString tok1_tan2 = path_data_tokens.first(); path_data_tokens.pop_front();
                QString tok2_tan2 = path_data_tokens.first(); path_data_tokens.pop_front();

                QString tok1 = path_data_tokens.first(); path_data_tokens.pop_front();
                QString tok2 = path_data_tokens.first(); path_data_tokens.pop_front();

                //http://www.w3.org/TR/SVG/paths.html#PathDataCubicBezierCommands
                //"The first control point is assumed to be the reflection of the second
                //control point on the previous command relative to the current point. "
                if (QString::compare(token, "S") == 0) {
                    tan1 = current_pos + (current_pos - pts[pts.size()-2]);
                }
                else {
                    tan1 = current_pos - pts[pts.size()-2];
                }
                tan2 = QVector2D(tok1_tan2.toFloat(), tok2_tan2.toFloat());
                pos = QVector2D(tok1.toFloat(), tok2.toFloat());

            }

            if (QString::compare(last_command, "C") == 0 || QString::compare(last_command, "S") == 0) { //absolute

                //add a bezier segment
                qDebug() << "Adding absolute Bezier segment" << current_pos << tan1 << tan2 << pos;
                if (pts.empty()) {
                    pts.push_back(current_pos);
                }
                pts.push_back(tan1);
                pts.push_back(tan2);
                pts.push_back(pos);

                current_pos = pos;

            }
            else { //relative

                //add a bezier segment
                qDebug() << "Adding relative Bezier segment" << current_pos << current_pos+tan1 << current_pos+tan2 << current_pos+pos;
                if (pts.empty()) {
                    pts.push_back(current_pos);
                }
                pts.push_back(current_pos + tan1);
                pts.push_back(current_pos + tan2);
                pts.push_back(current_pos + pos);

                current_pos = current_pos + pos;

            }

        }
        else if (QString::compare(last_command, "L") == 0 || QString::compare(last_command, "l") == 0) { //L - line (absolute), l (relative)

            QString tok1 = path_data_tokens.first(); path_data_tokens.pop_front();
            QString tok2 = path_data_tokens.first(); path_data_tokens.pop_front();

            QVector2D last_pos = current_pos;

            if (QString::compare(last_command, "L") == 0) {
                current_pos = QVector2D(tok1.toFloat(), tok2.toFloat());
            }
            else {
                current_pos = current_pos + QVector2D(tok1.toFloat(), tok2.toFloat());
            }

            QVector2D tan1 = last_pos + (current_pos - last_pos) * 0.25;
            QVector2D tan2 = last_pos + (current_pos - last_pos) * 0.75;

            //add a bezier segment which closes this path
            pts.push_back(last_pos);
            pts.push_back(tan1);
            pts.push_back(tan2);
            pts.push_back(current_pos);

        }

    };

    //next we normalize this
    QVector2D bmin(FLT_MAX, FLT_MAX);
    QVector2D bmax(-FLT_MAX, -FLT_MAX);

    for (int i=0; i<pts.size(); ++i) {
        bmin.setX(qMin(bmin.x(), pts[i].x()));
        bmin.setY(qMin(bmin.y(), pts[i].y()));
        bmax.setX(qMax(bmax.x(), pts[i].x()));
        bmax.setY(qMax(bmax.y(), pts[i].y()));
    }

    const float scale_fac = 4.0f;
    const float scale_val = qMax(bmax.x() - bmin.x(), bmax.y() - bmin.y()); //scale both X and Y with same val to preserve aspect ratio
    for (int i=0; i<pts.size(); ++i) {
        pts[i] = (pts[i] - bmin);
        pts[i] /= scale_val;
        pts[i].setX(pts[i].x() * scale_fac - scale_fac/2.0f);
        pts[i].setY(scale_fac - (pts[i].y() * scale_fac));
    }
}

int BezierCurve::GetPreviousSegmentIndex(const int i) const
{
    int next = i - 3;

    if (next < 0) {
        next = pts.size() - 4;
    }

    return next;
}

int BezierCurve::GetNextSegmentIndex(const int i) const
{
    int next = i + 3;

    if (next >= pts.size()-1) {
        next = 0;
    }

    return next;
}

void BezierCurve::SetClosed(bool b)
{
    closed = b;
}

void BezierCurve::Draw2D()
{
    glColor3f(0.75f, 0.75f, 0.75f);
    glBegin(GL_LINES);
    for (int i=0; i+3<pts.size(); i+=3) {

        glVertex2f(pts[i].x(), pts[i].y());
        glVertex2f(pts[i+1].x(), pts[i+1].y());

        glVertex2f(pts[i+2].x(), pts[i+2].y());
        glVertex2f(pts[i+3].x(), pts[i+3].y());

    }
    glEnd();

    glColor3f(0.5f, 0.5f, 0.5f);
    closed ? glBegin(GL_LINE_LOOP) : glBegin(GL_LINE_STRIP);

    for (int i=0; i<samples.size(); ++i) {
        glVertex2f(samples[i].x(), samples[i].y());
    }

    glEnd();

    glColor3f(0.25f, 1.0f, 0.25f);
    glPointSize(6.0f);
    glBegin(GL_POINTS);
    for (int i=0; i<pts.size(); ++i) {
        glVertex2f(pts[i].x(), pts[i].y());
    }

    glColor3f(0.25f, 0.25f, 1.0f);
    if (selected >= 0) {
        glVertex2f(pts[selected].x(), pts[selected].y());
    }

    glEnd();
}


bool BezierCurve::IsClosed()
{
    return closed;
}

void BezierCurve::SetSamplesPerSegment(const int i)
{
    samples_per_segment = i;
}

void BezierCurve::UpdateSamples()
{
    samples.clear();

    for (int i=0; i+3<pts.size(); i+=3) {

        for (int j=0; j<samples_per_segment; ++j) {

            QVector2D p;
            const float interp = float(j) / float(samples_per_segment);
            SamplePoint(pts[i], pts[i+1], pts[i+2], pts[i+3], interp, p);
            samples.push_back(p);
        }
    }
}

const QList <QVector2D> & BezierCurve::Samples() const
{
    return samples;
}

const QList <QVector2D> & BezierCurve::Points() const
{
    return pts;
}

const QVector2D BezierCurve::Point(const int i) const
{
    return pts[i];
}

void BezierCurve::SelectPoint(const QVector2D & p, const float select_size)
{
    selected = -1;
    float min_dist = select_size;

    for (int i=0; i<pts.size(); ++i) {

        float each_dist = (pts[i] - p).lengthSquared();

        if (each_dist < min_dist) {
            min_dist = each_dist;
            selected = i;
        }
    }
}

void BezierCurve::SetSelectedPoint(const int index)
{
    selected = index;
}

int BezierCurve::SelectedPoint() const
{
    return selected;
}

void BezierCurve::MoveSelectedPoint(const QVector2D & p, const bool keep_g1, const bool equal_lengths)
{
    if (selected == -1) {
        qDebug() << "BezierCurve::MoveSelectedPoint - Warning, trying to move "
                    "no selected point";
        return;
    }
    else if (selected >= pts.size()) {
        qDebug() << "BezierCurve::MoveSelectedPoint - Warning, selected point "
                    "index >= size of points";
    }

    if (selected % 3 != 0) { //if tangent point we may need to keep collinear
        int other_index = -1;
        int ctrl_pt = -1;

        if (selected % 3 == 1) {

            ctrl_pt = selected-1;

            if (selected == 1) { //it's the first segment
                other_index = closed ? pts.size()-2 : -1;
            }
            else {
                other_index = selected-2;
            }

        }
        else if (selected % 3 == 2) {

            ctrl_pt = selected+1;

            if (selected == pts.size()-2) { //it's the last segment
                other_index = closed ? 1 : -1;
            }
            else {
                other_index = selected+2;
            }
        }

        //make other_index point collinear
        if (keep_g1 && other_index >= 0) {

            float dist = (pts[other_index]-pts[ctrl_pt]).length();
            if(equal_lengths)
                pts[other_index] = pts[ctrl_pt] - (p-pts[ctrl_pt]);
            else
                pts[other_index] = pts[ctrl_pt] + (p-pts[ctrl_pt]).normalized() * (-dist);

        }

        pts[selected] = p;
    }
    else if (selected % 3 == 0) {
        QVector2D movement = p - pts[selected];

        if (closed && (selected == 0 || selected == pts.size()-1)) {
            pts[0] += movement;
            pts[pts.size()-1] += movement;
            if(keep_g1)
            {
                pts[1] += movement;
                pts[pts.size()-2] += movement;
            }
        }
        else {
            if(keep_g1)
            {
                if (selected > 0) {
                    pts[selected-1] += movement;
                }
                if (selected < pts.size()-1) {
                    pts[selected+1] += movement;
                }
            }
            pts[selected] += movement;
        }
    }
}

void BezierCurve::DeletePoint(const int segment_index)
{

    if (closed && pts.size() <= 7) {
        qDebug() << "BezierCurve::DeleteSelectedPoint() - Warning, not "
                    "deleting since curve closed and only"
                 << pts.size() << "pts.";
        return;
    }
    if (!closed && pts.size() <= 4) {
        qDebug() << "BezierCurve::DeleteSelectedPoint() - Warning, not "
                    "deleting since curve not closed and only"
                 << pts.size() << "pts.";
        return;
    }

    // move to nearest endpoint
    int index = segment_index;
    if (index % 3 == 1) {
        --index;
    }
    if (index % 3 == 2) {
        ++index;
    }

    // now we remove the endpoint, plus the tangent before and after
    // TODO: this could be more optimal in that the shape is better
    // preserved following removal
    if (index > 0 && index < pts.size()-1) {
        pts.removeAt(index-1);
        pts.removeAt(index-1);
        pts.removeAt(index-1);
    }
    else {
        if (!closed && index == 0) {
            pts.removeFirst();
            pts.removeFirst();
            pts.removeFirst();
        }
        else if (!closed && index == pts.size()-1) {
            pts.removeLast();
            pts.removeLast();
            pts.removeLast();
        }
        else { //closed and at the endpoint
            QVector2D t2 = pts[2];

            pts.removeLast();
            pts.removeLast();

            pts.removeFirst();
            pts.removeFirst();
            pts.removeFirst();

            pts.push_back(t2);
            pts.push_back(pts.first());
        }
    }
}

void BezierCurve::DeleteSelectedPoint()
{
    DeletePoint(selected);
    selected = -1;
}

void BezierCurve::InsertPoint(const int segment_index)
{
    if (segment_index < 0 || segment_index >= pts.size()) {
        qDebug() << "BezierCurve::InsertPoint() - Warning, selected is"
                 << segment_index;
        return;
    }

    int index = segment_index;
    index -= (index % 3);

    QList <QVector2D> subdiv_pts;
    Subdivide(index, 0.5f);
}

void BezierCurve::InsertSelectedPoint()
{
    InsertPoint(selected);
}

void BezierCurve::UnselectPoint()
{
    selected = -1;
}

void BezierCurve::SamplePoint(const QVector2D & p1, const QVector2D & t1, const QVector2D & t2, const QVector2D & p2, const float t, QVector2D & point)
{
    const float omt = 1.0f - t;
    const float omt_2 = omt * omt;
    const float t_2 = t * t;

    const float a0 = omt_2 * omt;
    const float a1 = 3 * omt_2 * t;
    const float a2 = 3 * omt * t_2;
    const float a3 = t_2 * t;

    point = (p1 * a0) + (t1 * a1) + (t2 * a2) + (p2 * a3);
}

void BezierCurve::SampleTangent(const QVector2D & p1, const QVector2D & t1, const QVector2D & t2, const QVector2D & p2, const float t, QVector2D & tangent)
{
    const float omt = 1.0f - t;

    const QVector2D p10 = p1 * omt + t1 * t;
    const QVector2D p11 = t1 * omt  + t2 * t;
    const QVector2D p12 = t2 * omt  + p2 * t;

    const QVector2D p20 = p10 * omt  + p11 * t;
    const QVector2D p21 = p11 * omt  + p12 * t;

    tangent = (p21 - p20).normalized();
}

bool BezierCurve::IsPointInside(const QVector2D & p)
{

    //use ray casting algorithm, count intersection points
    //with x value less than p.x
    int intersects_left_of_p = 0;

    for (int i=0; i<samples.size(); ++i) {
        int i1 = i;
        int i2 = (i+1)%samples.size();
        QVector2D intersect;

        if (samples[i1].x() < p.x() || samples[i2].x() < p.x()) {
            if (GLutils::LineRayIntersection(samples[i1], samples[i2], p, QVector2D(1, 0), intersect)) {
                if (intersect.x() < p.x()) {
                    ++intersects_left_of_p;
                }
            }
        }
    }

    return (intersects_left_of_p % 2) == 1;

}

void BezierCurve::Subdivide(const int i, const float t)
{
    const float omt = 1.0f - t;

    const QVector2D p00 = pts[i];
    const QVector2D p01 = pts[i+1];
    const QVector2D p02 = pts[i+2];
    const QVector2D p03 = pts[i+3];

    const QVector2D p10 = p00 * omt + p01 * t;
    const QVector2D p11 = p01 * omt  + p02 * t;
    const QVector2D p12 = p02 * omt  + p03 * t;

    const QVector2D p20 = p10 * omt  + p11 * t;
    const QVector2D p21 = p11 * omt  + p12 * t;

    const QVector2D p30 = p20 * omt  + p21 * t;

    //Use these:
    //first ctrl point - p00
    //start tangent - p10
    //end tangent - p20
    //new endpoint - p30
    //start tangent - p21
    //end tangent - p12
    //end ctrl point - p03

    //add in the middle ctrl points
    pts.removeAt(i+1);
    pts.removeAt(i+1);

    pts.insert(i+1, p10);
    pts.insert(i+2, p20);
    pts.insert(i+3, p30);
    pts.insert(i+4, p21);
    pts.insert(i+5, p12);
}

void BezierCurve::Subdivide(const BezierCurve & c, const int i, const float t, QList <QVector2D> & subdiv_pts)
{
    const float omt = 1.0f - t;

    const QVector2D p00 = c.pts[i];
    const QVector2D p01 = c.pts[i+1];
    const QVector2D p02 = c.pts[i+2];
    const QVector2D p03 = c.pts[i+3];

    const QVector2D p10 = p00 * omt + p01 * t;
    const QVector2D p11 = p01 * omt  + p02 * t;
    const QVector2D p12 = p02 * omt  + p03 * t;

    const QVector2D p20 = p10 * omt  + p11 * t;
    const QVector2D p21 = p11 * omt  + p12 * t;

    const QVector2D p30 = p20 * omt  + p21 * t;

    //Use these:
    //first ctrl point - p00
    //start tangent - p10
    //end tangent - p20
    //new endpoint - p30
    //start tangent - p21
    //end tangent - p12
    //end ctrl point - p03
    subdiv_pts.push_back(p00);
    subdiv_pts.push_back(p10);
    subdiv_pts.push_back(p20);
    subdiv_pts.push_back(p30);
    subdiv_pts.push_back(p21);
    subdiv_pts.push_back(p12);
    subdiv_pts.push_back(p03);

}

void BezierCurve::GetLineIntersections(const QVector2D & line_p, const QVector2D & line_d, QList <BezierCurvePoint> & intersects)
{

    intersects.clear();

    QVector <float> segment_line_dists = QVector <float> (tests_sample_rate+1);

    QVector2D line_n(-line_d.y(), line_d.x());

    for (int i=0; i+3<pts.size(); i+=3) {

        QVector <QVector2D> segment_samples;
        SamplePointsForSegment(i, segment_samples);

        for (int t=0; t<=tests_sample_rate; ++t) {

            //compute distance from line along direction split_n (ortho to split_d)
            const float each_dist = QVector2D::dotProduct(segment_samples[t] - line_p, line_n);
            segment_line_dists[t] = each_dist;

        }

        //now we iterate through and observe changes in signed distance
        for (int t=0; t<tests_sample_rate; ++t) {
            if (segment_line_dists[t] * segment_line_dists[t+1] > 0.0f) {
                continue;
            }

            BezierCurvePoint bcp;
            bcp.segment = i;
            bcp.point = (segment_samples[t] + segment_samples[t+1]) * 0.5f;
            bcp.t = (float(t) + 0.5f)/float(tests_sample_rate);
            intersects.push_back(bcp);
        }
    }

    //sort intersections along the ray direction
    for (int i=0; i<intersects.size(); ++i) {
        for (int j=i+1; j<intersects.size(); ++j) {

            if (QVector2D::dotProduct(intersects[i].point, line_d) > QVector2D::dotProduct(intersects[j].point, line_d)) {
                intersects.swapItemsAt(i, j);
            }
        }
    }
}

void BezierCurve::GetLineSegmentIntersections(const QVector2D & line_p1, const QVector2D & line_p2, QList <BezierCurvePoint> & intersects)
{
     for (int i=0; i+3<pts.size(); i+=3) {
        QVector <QVector2D> segment_samples;
        SamplePointsForSegment(i, segment_samples);

        for (int t=0; t<tests_sample_rate; ++t) {

            QVector2D intersect;

            if (GLutils::LineLineIntersection(line_p1, line_p2, segment_samples[t], segment_samples[t+1], intersect)) {
                BezierCurvePoint bcp;
                bcp.segment = i;
                bcp.point = (segment_samples[t] + segment_samples[t+1]) * 0.5f;
                bcp.t = (float(t) + 0.5f)/float(tests_sample_rate);
                intersects.push_back(bcp);
            }
        }
    }
}

void BezierCurve::SplitAlongLine(const QVector2D & split_p, const QVector2D & split_d, QList <BezierCurve> & curves)
{  
    //1.  compute the intersections with the bezier curve
    QList <BezierCurvePoint> split;
    GetLineIntersections(split_p, split_d, split);


    //2.  go through making loops starting at each split segment (in both directions)    
    //if the line didn't split this curve, we just return
    if (split.empty()) {
        return;
    }

    //if there's an odd number of intersections, we return
    if (split.size() % 2 == 1) {
        qDebug() << "BezierCurve::SplitAlongLine() - Problem, grid line "
                    "made odd number of intersections";
        return;
    }

    //subdivide if two adjacent split point share the same segment
    //TODO!!!: special case where split_segment[start] and
    //split_segment[start+1] are the same
    //(i.e. no position ctrl points on one side, but interpolated curve crosses)
    for (int i=0; i+1<split.size(); i+=2) {
        if (split[i].segment == split[i+1].segment) {
            Subdivide(split[i].segment, (split[i].t + split[i+1].t) * 0.5f);

            //reprocess intersections (since curve parameterization changed)
            GetLineIntersections(split_p, split_d, split);
        }
    }

    //create a new curve which follows this split
    QVector <bool> visited_forward = QVector <bool> (split.size(), false);
    QVector <bool> visited_backward = QVector <bool> (split.size(), false);

    //visit all split lines twice (in both directions)
    while (true) {

        int start = -1;
        bool go_forward = false;

        //figure out split to start at and direction
        for (int i=0; i<visited_backward.size(); i+=2) {
            if (!visited_forward[i]) {
                start = i;
                go_forward = true;
                break;
            }
            if (!visited_backward[i]) {
                start = i;
                go_forward = false;
                break;
            }
        }

        if (start == -1) {
            break;
        }

        //ok, now we go through making a new bezier curve
        int i;
        int i0;

        BezierCurve new_curve;
        new_curve.SetClosed(true);       

        QList <QVector2D> subdiv_pts;
        Subdivide((*this), split[start].segment, split[start].t, subdiv_pts);
        QList <QVector2D> next_subdiv_pts;
        Subdivide((*this), split[start+1].segment, split[start+1].t, next_subdiv_pts);

        if (go_forward) {
            new_curve.AddPoint(subdiv_pts[0]);
            new_curve.AddPoint(subdiv_pts[1]);
            new_curve.AddPoint(subdiv_pts[2]);
            new_curve.AddPoint(subdiv_pts[3]);
            new_curve.AddPoint(subdiv_pts[3] * 0.75 + next_subdiv_pts[3] * 0.25);
            new_curve.AddPoint(subdiv_pts[3] * 0.25 + next_subdiv_pts[3] * 0.75);
            new_curve.AddPoint(next_subdiv_pts[3]);
            new_curve.AddPoint(next_subdiv_pts[4]);
            new_curve.AddPoint(next_subdiv_pts[5]);
            new_curve.AddPoint(next_subdiv_pts[6]);

            visited_forward[start] = true;
            visited_forward[start+1] = true;

            i0 = split[start].segment;
            i = split[start+1].segment;
        }
        else {           
            new_curve.AddPoint(next_subdiv_pts[0]);
            new_curve.AddPoint(next_subdiv_pts[1]);
            new_curve.AddPoint(next_subdiv_pts[2]);
            new_curve.AddPoint(next_subdiv_pts[3]);
            new_curve.AddPoint(subdiv_pts[3] * 0.25 + next_subdiv_pts[3] * 0.75);
            new_curve.AddPoint(subdiv_pts[3] * 0.75 + next_subdiv_pts[3] * 0.25);
            new_curve.AddPoint(subdiv_pts[3]);
            new_curve.AddPoint(subdiv_pts[4]);
            new_curve.AddPoint(subdiv_pts[5]);
            new_curve.AddPoint(subdiv_pts[6]);

            visited_backward[start] = true;
            visited_backward[start+1] = true;

            i0 = split[start+1].segment;
            i = split[start].segment;
        }

        i = GetNextSegmentIndex(i);

        //halt when the curve is closed
        while (i0 != i) {
            //check if this segment is in the split list
            int ind = -1;
            for (int j=0; j<split.size(); ++j) {
                if (i == split[j].segment) {
                    ind = j;
                    break;
                }
            }

            if (ind == -1) { //the segment is NOT in the split_list
                //no subdivision
                new_curve.AddPoint(pts[i+1]);
                new_curve.AddPoint(pts[i+2]);
                new_curve.AddPoint(pts[i+3]);
            }
            else { //the segment IS in the split_list

                if (go_forward && visited_forward[ind]) {
                    break;
                }
                else if (!go_forward && visited_backward[ind]) {
                    break;
                }

                //int next_ind = go_forward ? ind+1 : ind-1;
                int next_ind = (ind % 2) == 0 ? ind+1 : ind-1;

                QList <QVector2D> subdiv_pts;
                Subdivide((*this), split[ind].segment, split[ind].t, subdiv_pts);
                QList <QVector2D> next_subdiv_pts;
                Subdivide((*this), split[next_ind].segment, split[next_ind].t, next_subdiv_pts);

                //segment to line
                new_curve.AddPoint(subdiv_pts[1]);
                new_curve.AddPoint(subdiv_pts[2]);
                new_curve.AddPoint(subdiv_pts[3]);

                //intermediate along line
                new_curve.AddPoint(subdiv_pts[3] * 0.75 + next_subdiv_pts[3] * 0.25);
                new_curve.AddPoint(subdiv_pts[3] * 0.25 + next_subdiv_pts[3] * 0.75);
                new_curve.AddPoint(next_subdiv_pts[3]);

                //next segment
                new_curve.AddPoint(next_subdiv_pts[4]);
                new_curve.AddPoint(next_subdiv_pts[5]);
                new_curve.AddPoint(next_subdiv_pts[6]);

                //mark as visited
                if (go_forward) {
                    visited_forward[ind] = true;
                    visited_forward[next_ind] = true;
                }
                else {
                    visited_backward[ind] = true;
                    visited_backward[next_ind] = true;
                }

                i = split[next_ind].segment;
            }

            //increment i (move along the loop)
            i = GetNextSegmentIndex(i);
        }

        //add this curve to the list
        new_curve.UpdateSamples();
        curves.push_back(new_curve);
    }
}

void BezierCurve::SubdivideLongestSegment()
{
    //1.  figure out the longest segment
    int max_segment_index = -1;
    float max_segment_length = -FLT_MAX;

    for (int i=0; i<(pts.size()-1)/3; ++i) {

        float each_segment_len = 0.0f;

        for (int j=0; j<samples_per_segment-1; ++j) {

            QVector2D & sample1 = samples[i*samples_per_segment+j+0];
            QVector2D & sample2 = samples[i*samples_per_segment+j+1];

            each_segment_len += (sample2 - sample1).length();

        }

        if (each_segment_len > max_segment_length) {
            max_segment_length = each_segment_len;
            max_segment_index = i*3;
        }

    }

    //2.  split the longest segment
    InsertPoint(max_segment_index);

    //3.  update samples
    UpdateSamples();

}

void BezierCurve::GetCurvePointCorrespondence(const BezierCurve & other, int & corresp_offset, bool & corresp_forward) const
{

    //1.  ensure the two curves have the same # of control points
    if (pts.size() != other.pts.size()) {
        qDebug() << "BezierCurve::GetCurvePointCorrespondence() - Warning, "
                    "two curves do not have same # of control points";
        corresp_offset = 0;
        return;
    }

    //2.  iterate around the ring trying each offset,
    //    evaluate distances to control points
    int min_offset_index = -1;
    float min_offset_dist = FLT_MAX;
    bool min_offset_forward = true;

    //remove the duplicated last point for curves c1, c2
    BezierCurve c1 = (*this);
    BezierCurve c2 = other;
    c1.pts.removeLast();
    c2.pts.removeLast();

    //translate c1, c2 control points by their centroid
    c1.TranslateControlPoints(-c1.GetControlPointCentroid());
    c2.TranslateControlPoints(-c2.GetControlPointCentroid());

    //march backward and forward looking for best correspondence
    //k = 0 -> go forward
    //k = 1 -> go backward
    for (int k=0; k<2; ++k) {

        for (int i=0; i<c1.pts.size(); i+=3) {
            //note: (i+=3) only doing correspondences matching position
            //ctrl point <--> position ctrl point

            //sum up lengths between all control points using offset i
            //(even the tangent-related ones)
            float each_offset_dist = 0.0f;
            for (int j=0; j<c1.pts.size(); ++j) {

                int ind1 = j;
                int ind2;

                if (k == 0) {
                    ind2 = (j + i) % c1.pts.size();
                }
                else {
                    ind2 = (-j + i + c1.pts.size()) % c1.pts.size();
                }

                each_offset_dist += (c1.pts[ind1] - c2.pts[ind2]).length();
            }

            if (each_offset_dist < min_offset_dist) {
                min_offset_dist = each_offset_dist;
                min_offset_index = i;
                min_offset_forward = (k == 0);
            }
        }
    }

    //3.  use the offset where the sum of distances to control points is minimized
    corresp_offset = min_offset_index;
    corresp_forward = min_offset_forward;
}

void BezierCurve::InterpolateCurvePointCorrespondence(const BezierCurve & other, const int corresp_offset, const bool corresp_forward, const float t, BezierCurve & interp_curve) const
{
    //1.  ensure the two curves have the same # of control points
    if (pts.size() != other.pts.size()) {
        qDebug() << "BezierCurve::InterpolateCurvePointCorrespondnece() - "
                    "Warning, two curves do not have same # of control points";
        return;
    }

    //2.  do linear interpolation of control point positions for new curve
    interp_curve.ClearPoints();

    //remove the duplicated last point for curves c1, c2
    BezierCurve c1 = (*this);
    BezierCurve c2 = other;
    c1.pts.removeLast();
    c2.pts.removeLast();

    for (int i=0; i<c1.pts.size(); ++i) {

        int ind1 = i;
        int ind2;
        if (corresp_forward) {
            ind2 = (i + corresp_offset) % c1.pts.size();
        }
        else {
            ind2 = (-i + corresp_offset + c1.pts.size()) % c1.pts.size();
        }

        QVector2D interp_pt = c1.pts[ind1] * (1.0f - t) + c2.pts[ind2] * t;

        interp_curve.AddPoint(interp_pt);
    }

    //duplicate last pt for interpolated curve
    interp_curve.pts.push_back(interp_curve.pts.first());

    //3.  update samples of interpolated curve
    interp_curve.UpdateSamples();

}

float BezierCurve::Length() const
{
    float len = 0.0f;

    for (int i=0; i<samples.size()-1; ++i) {        
        len += (samples[i+1] - samples[i]).length();
    }

    return len;
}

void BezierCurve::Reverse()
{
    for(int i = 0; i < (pts.size()/2); i++) {
        pts.swapItemsAt(i, pts.size()-(1+i));
    }

    for(int i = 0; i < (samples.size()/2); i++) {
        samples.swapItemsAt(i, samples.size()-(1+i));
    }
}

void BezierCurve::GetSubCurve(const BezierCurvePoint & b1, const BezierCurvePoint & b2, BezierCurve & sub_curve)
{

    QList <QVector2D> b1_subdiv;
    QList <QVector2D> b2_subdiv;

    //2.  construct the curve
    //  2a) the 2 points are in the same segment (final curve has just 4 pts then)
    if (b1.segment == b2.segment && b1.t != b2.t) {

        BezierCurve segment_curve;

        //we now need to adjust p0's t value since we took (1-p1.t) away from it
        //p0's new t value is p0.t / p1.t (note that p1.t > p0.t)
        if (b1.t < b2.t) {
            Subdivide((*this), b2.segment, b2.t, b1_subdiv);
            segment_curve.SetPoints(b1_subdiv);
            Subdivide(segment_curve, 0, b1.t / b2.t, b2_subdiv);
        }
        else {
            Subdivide((*this), b1.segment, b1.t, b1_subdiv);
            segment_curve.SetPoints(b1_subdiv);
            Subdivide(segment_curve, 0, b2.t / b1.t, b2_subdiv);
        }

        sub_curve.ClearPoints();
        sub_curve.AddPoint(b2_subdiv[3]);
        sub_curve.AddPoint(b2_subdiv[4]);
        sub_curve.AddPoint(b2_subdiv[5]);
        sub_curve.AddPoint(b2_subdiv[6]);

        if (b1.t > b2.t) {
            sub_curve.Reverse();
        }       

    }
    else {

        //2b) the two points are in different bezier curve segments (the "easy" case)
        Subdivide((*this), b1.segment, b1.t, b1_subdiv);
        Subdivide((*this), b2.segment, b2.t, b2_subdiv);

        sub_curve.ClearPoints();

        //i) add last part of first segment's subdivision
        sub_curve.AddPoint(b1_subdiv[3]);
        sub_curve.AddPoint(b1_subdiv[4]);
        sub_curve.AddPoint(b1_subdiv[5]);
        sub_curve.AddPoint(b1_subdiv[6]);

        //ii) add the intermediate segments
        int i = GetNextSegmentIndex(b1.segment);

        while (i != b2.segment) {
            sub_curve.AddPoint(pts[i+1]);
            sub_curve.AddPoint(pts[i+2]);
            sub_curve.AddPoint(pts[i+3]);

            i = GetNextSegmentIndex(i);

        }

        //iii) add the first part of the last segment's subdivision
        sub_curve.AddPoint(b2_subdiv[1]);
        sub_curve.AddPoint(b2_subdiv[2]);
        sub_curve.AddPoint(b2_subdiv[3]);

    }

    //make samples dense so we get a close subcurve
    sub_curve.SetSamplesPerSegment(tests_sample_rate);
    sub_curve.UpdateSamples();

}

void BezierCurve::GetPointsTangentsAlongCurve(const int num, QVector <QVector2D> & points, QVector <QVector2D> & tangents) const
{

    if (num < 2) {
        qDebug() << "BezierCurve::GetPointsTangentsAlongCurve() - "
                    "Warning, parameter less than 2 passed in";
        return;
    }

    if (pts.size() < 4) {
        qDebug() << "BezierCurve::GetPointsTangentsAlongCurve() - "
                    "Warning, less than 4 control points";
        return;
    }

    points.resize(num+1);
    tangents.resize(num+1);

    //assign first point/tangent
    points[0] = pts.first();
    SampleTangent(pts[0], pts[1], pts[2], pts[3], 0.0f, tangents[0]);

    //now deal with the middle ones by sampling along arc-length of curve
    int n_index = 1;
    float cur_len = 0.0f;
    const float total_len = Length();    

    //go through each segment
    for (int i=0; i+3<pts.size(); i+=3) {

        QVector <QVector2D> segment_samples;
        SamplePointsForSegment(i, segment_samples);       

        //sampling densely along it and counting the length we've travelled
        for (int j=0; j<tests_sample_rate; ++j) {

            cur_len += (segment_samples[j+1]-segment_samples[j]).length();

            const float len_needed = (float(n_index) / float(num)) * total_len;

            //if we've travelled the length of the spacing between samples,
            //we sample here
            if (cur_len >= len_needed) {

                const float t = float(j + 0.5f) / float(tests_sample_rate+1);

                SamplePoint(pts[i], pts[i+1], pts[i+2], pts[i+3], t, points[n_index]);
                SampleTangent(pts[i], pts[i+1], pts[i+2], pts[i+3], t, tangents[n_index]);

                ++n_index;
                if (n_index >= num) {
                    break;
                }
            }
        }      
    }

    //add the final endpoint as an extra point
    int i = pts.size()-4;
    SamplePoint(pts[i], pts[i+1], pts[i+2], pts[i+3], 1.0f, points[num]);
    SampleTangent(pts[i], pts[i+1], pts[i+2], pts[i+3], 1.0f, tangents[num]);
}

void BezierCurve::SamplePointsForSegment(const int i, QVector <QVector2D> & sample_segment) const
{

    sample_segment = QVector <QVector2D> (tests_sample_rate+1);

    for (int t=0; t<=tests_sample_rate; ++t) {
        //sample at each point along t
        SamplePoint(pts[i], pts[i+1], pts[i+2], pts[i+3], float(t)/float(tests_sample_rate), sample_segment[t]);
    }

}
