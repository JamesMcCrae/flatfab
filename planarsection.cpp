#include "planarsection.h"

QList <QVector3D> PlanarSection::convex_hull;

PlanarSection::PlanarSection()
{

    p = QVector3D(0.0f, 0.0f, 0.0f);

    t = QVector3D(1.0f, 0.0f, 0.0f);
    b = QVector3D(0.0f, 1.0f, 0.0f);
    n = QVector3D(0.0f, 0.0f, 1.0f);

    editing_sketch = false;
    slab_thickness = 0.06f;

    centroid = QVector2D(0, 0);
    area = 0.0f;

    part_of_cycle = false;
    connected = true;

    selected_weight = -1;

    quality_samples = 15;

    //update_slab_disp_list = false;
    //update_slabcurves_disp_list = false;
    //slab_disp_list = 0;
    //slabcurves_disp_list = 0;

    AddNewCurve();

    num_radial_sectors = 0;
    num_radial_points_per_sector = 0;

    radial = false;
    local_symmetry_mode = false;

}

BezierCurve & PlanarSection::GetCurve(const int i)
{
    return bez_curve[i];
}

int PlanarSection::GetNumCurves() const
{
    return bez_curve.size();
}

void PlanarSection::SetCurve(const int i, const BezierCurve & b)
{
    bez_curve[i] = b;
}

QList <BezierCurve> & PlanarSection::GetCurves()
{
    return bez_curve;
}

void PlanarSection::SetCurves(const QList <BezierCurve> & curves)
{
    bez_curve = curves;
}

void PlanarSection::AddNewCurve()
{
    BezierCurve b;
    QList <QVector2D> pts;

    sketch_pts.push_back(pts);
    bez_curve.push_back(b);
}

void PlanarSection::RemoveCurve(const int i)
{
    sketch_pts.removeAt(i);
    bez_curve.removeAt(i);
}

void PlanarSection::RemoveHoles()
{
    while (sketch_pts.size() > 1) {
        sketch_pts.pop_back();
    }
    while (bez_curve.size() > 1) {
        bez_curve.pop_back();
    }
}


void PlanarSection::SetNewBasis(const QVector3D & new_t, const QVector3D & new_n, const QVector3D & new_b, const QVector3D & new_p)
{   

    for (int j=0; j<bez_curve.size(); ++j) {

        const QList <QVector2D> & pts = bez_curve[j].Points();

        for (int i=0; i<pts.size(); ++i) {

            QVector2D each_new_pt = GetPoint2D(GLutils::GetVectorNewBasis(t, n, b, p, new_t, new_n, new_b, new_p, GetPoint3D(pts[i]), 1.0f, 1.0f, 1.0f));
            bez_curve[j].SetPoint(i, each_new_pt);

        }
    }

    t = new_t;
    n = new_n;
    b = new_b;
    p = new_p;

}

void PlanarSection::SetP(const QVector3D & v)
{
    p = v;
}

void PlanarSection::CreateSquare(const float size)
{

    QVector2D v(size, size);
    CreateRectangle(-v, v);

}

void PlanarSection::CreateRectangle(QVector2D min_v, QVector2D max_v)
{

    const float xsize = max_v.x() - min_v.x();
    const float ysize = max_v.y() - min_v.y();

    QVector2D p1 = min_v;
    QVector2D p2 = min_v + QVector2D(xsize, 0);
    QVector2D p3 = max_v;
    QVector2D p4 = min_v + QVector2D(0, ysize);   

    bez_curve[0].ClearPoints();

    bez_curve[0].AddPoint(p1 * 1.0f + p2 * 0.0f);
    bez_curve[0].AddPoint(p1 * 0.75f + p2 * 0.25f);
    bez_curve[0].AddPoint(p1 * 0.25f + p2 * 0.75f);
    bez_curve[0].AddPoint(p1 * 0.0f + p2 * 1.0f);
    bez_curve[0].AddPoint(p2 * 0.75f + p3 * 0.25f);
    bez_curve[0].AddPoint(p2 * 0.25f + p3 * 0.75f);
    bez_curve[0].AddPoint(p2 * 0.0f + p3 * 1.0f);
    bez_curve[0].AddPoint(p3 * 0.75f + p4 * 0.25f);
    bez_curve[0].AddPoint(p3 * 0.25f + p4 * 0.75f);
    bez_curve[0].AddPoint(p3 * 0.0f + p4 * 1.0f);
    bez_curve[0].AddPoint(p4 * 0.75f + p1 * 0.25f);
    bez_curve[0].AddPoint(p4 * 0.25f + p1 * 0.75f);
    bez_curve[0].AddPoint(p4 * 0.0f + p1 * 1.0f);

    UpdateCurveTrisSlab();

}

void PlanarSection::CreateCircle(QVector2D centre, const float radius, int num_sections)
{

    QList <QVector2D> pts;

    const int circ_samples = num_sections;
    const float pi2 = (2.0f * 3.14159f);

    for (int i=0; i<circ_samples; ++i) {

        const float t1 = float(i) / float(circ_samples) * pi2;
        const float t2 = (float(i) + 0.333f) / float(circ_samples) * pi2;
        const float t3 = (float(i) + 0.666f) / float(circ_samples) * pi2;

        pts.push_back(QVector2D(cosf(t1), sinf(t1)));
        pts.push_back(QVector2D(cosf(t2), sinf(t2)) * 1.01f);
        pts.push_back(QVector2D(cosf(t3), sinf(t3)) * 1.01f);

    }

    bez_curve[0].ClearPoints();

    for (int i=0; i<=pts.size(); ++i) {
        bez_curve[0].AddPoint(pts[i%pts.size()] * radius + centre);
    }

    UpdateCurveTrisSlab();

}

void PlanarSection::CreateRadial(QVector2D centre, const float base_rad, const int num_sectors, const float radii[9])
{

    QList <QVector2D> pts;

    const float pi2 = (2.0f * 3.14159f);

    for (int i=0; i<num_sectors; ++i) {

        for (int j=0; j<9; ++j) {

            const float t = (float(i) + float(j) / 9.0f) / float(num_sectors) * pi2;
            pts.push_back(QVector2D(cosf(t), sinf(t)) * radii[j]);

        }

    }

    bez_curve[0].ClearPoints();

    for (int i=0; i<=pts.size(); ++i) {
        bez_curve[0].AddPoint(pts[i%pts.size()] * base_rad + centre);
    }

    UpdateCurveTrisSlab();

}




// This rotates and copies each sector - assumes the
void PlanarSection::CopySectorRadially(QVector2D centre, int num_sectors, int points_per_sector)
{

    if (!radial) {
        qDebug()<< "Error - Section is not radial";
        return;
    }

    QList <QLineF> lines;
    QList <QVector2D> pts;

    if (bez_curve[0].GetNumControlPoints() < points_per_sector) {
        return;
    }

    // Convert QVector2Ds to QLineFs
    for(int i = 0; i < points_per_sector; i++) {
        lines.push_back(QLineF(centre.toPointF(), bez_curve[0].Point(i).toPointF()));
        pts.push_back(bez_curve[0].Point(i));
    }

    const float angularLength = 360.0f/num_sectors;

//    // This ensures that keep_g1 is satisfied
//    if(bez_curve[0].SelectedPoint() == 0)
//    {
//        lines.back() = QLineF(centre.toPointF(), bez_curve[0].Point(bez_curve[0].Points().size() - 2).toPointF());
//        lines.back().setAngle(lines.back().angle() - angularLength);
//        pts.back() = QVector2D(lines.back().x2(), lines.back().y2());
//    }

    for (int i=1; i<num_sectors; ++i) {

        for (int j=0; j<points_per_sector; ++j) {

            lines[j].setAngle(lines[j].angle() - angularLength);
            pts.push_back(QVector2D(lines[j].x2(), lines[j].y2()));

        }

    }
    pts.push_back(pts.front());

    bez_curve[0].ClearPoints();

    for (int i=0; i<pts.size(); ++i) {
        bez_curve[0].AddPoint(pts[i]);
    }

    UpdateCurveTrisSlab();
}


void PlanarSection::CreateRadialHoles(QVector2D centre, const float base_rad, const int num_sectors, const float radii[9])
{

    QVector <QList <QVector2D> > pts;
    pts.resize(num_sectors);

    const float pi2 = (2.0f * 3.14159f);

    if (num_sectors > 1) {
        for (int i=0; i<num_sectors; ++i) {

            const float t1 = float(i) / float(num_sectors) * pi2;
            const float t2 = float(i + 0.25f) / float(num_sectors) * pi2;
            const float t3 = float(i + 0.50f) / float(num_sectors) * pi2;
            const float t4 = float(i + 0.75f) / float(num_sectors) * pi2;

            QVector2D c1(cosf(t1), sinf(t1));
            QVector2D c2(cosf(t2), sinf(t2));
            QVector2D c3(cosf(t3), sinf(t3));
            QVector2D c4(cosf(t4), sinf(t4));

            pts[i].push_back(c1 * radii[0]);
            pts[i].push_back(c1 * (radii[0] * 0.75f + radii[3] * 0.25f));
            pts[i].push_back(c1 * (radii[0] * 0.25f + radii[3] * 0.75f));
            pts[i].push_back(c1 * radii[3]);
            pts[i].push_back(c2 * radii[2]);
            pts[i].push_back(c3 * radii[2]);
            pts[i].push_back(c4 * radii[3]);
            pts[i].push_back(c4 * (radii[3] * 0.75f + radii[0] * 0.25f));
            pts[i].push_back(c4 * (radii[3] * 0.25f + radii[0] * 0.75f));
            pts[i].push_back(c4 * radii[0]);
            pts[i].push_back(c3 * radii[1]);
            pts[i].push_back(c2 * radii[1]);
            pts[i].push_back(c1 * radii[0]);

        }
    }
    else {

    }

    for (int i=0; i<num_sectors; ++i) {

        AddNewCurve();

        for (int j=0; j<pts[i].size(); ++j) {
            bez_curve.last().AddPoint(pts[i][j] * base_rad + centre);
        }

    }

    UpdateCurveTrisSlab();

}

void PlanarSection::MoveP(const QVector3D & v)
{

    const QVector2D translate = GetPoint2D(p) - GetPoint2D(v);
    for (int c=0; c<bez_curve.size(); ++c) {
        bez_curve[c].TranslateControlPoints(translate);
    }
    p = v;

}

void PlanarSection::AddWeight(const QVector2D & pos_2d, const float mass)
{

    ExternalWeight weight;

    weight.mass_kg = mass;
    weight.p = pos_2d;
    //weight.rad = powf(mass, 1.0f / 3.0f) / 3.14159f;
    weight.rad = powf(mass, 1.0f / 3.0f) / 3.14159f * 5.0f; //exaggerate size for some figures
    weight.height = weight.rad;

    weights.push_back(weight);

}

void PlanarSection::AddWeightAtMousePos(const QVector2D & mouse_pos, const float mass)
{

    QVector3D intersect;
    MouseRayIntersect(mouse_pos, intersect);

    AddWeight(GetPoint2D(intersect), mass);

}

void PlanarSection::RemoveWeights()
{
    selected_weight = -1;
    weights.clear();
}

int PlanarSection::GetNumWeights()
{
    return weights.size();
}

QVector3D PlanarSection::GetWeightPosition(const int i)
{
    return GetPoint3D(weights[i].p);
}

QVector3D PlanarSection::GetWeightForce(const int i)
{
    return QVector3D(0, -9.8, 0) * weights[i].mass_kg;
}

float PlanarSection::GetWeightMass(const int i)
{
    return weights[i].mass_kg;
}

void PlanarSection::Scale(const float x, const float y)
{
    for (int i=0; i<bez_curve.size(); ++i) {
        bez_curve[i].ScaleControlPoints(x, y);
    }
}

void PlanarSection::Rotate(const float angle_rad)
{
    for (int i=0; i<bez_curve.size(); ++i) {
        bez_curve[i].RotateControlPoints(angle_rad);
    }
}

void PlanarSection::FlipN()
{
    t = -t;
    n = -n;
    for (int i=0; i<bez_curve.size(); ++i) {
        bez_curve[i].ScaleControlPoints(-1.0f, 1.0f);
    }
}

void PlanarSection::FlipTB()
{
    t = -t;
    b = -b;
    for (int i=0; i<bez_curve.size(); ++i) {
        bez_curve[i].ScaleControlPoints(-1.0f, -1.0f);
    }
}

void PlanarSection::SetT(const QVector3D & v)
{
    t = v;
}

void PlanarSection::SetN(const QVector3D & v)
{
    n = v;
}

void PlanarSection::SetB(const QVector3D & v)
{
    b = v;
}

QVector3D PlanarSection::P() const
{
    return p;
}

QVector3D PlanarSection::T() const
{
    return t;
}

QVector3D PlanarSection::N() const
{
    return n;
}

QVector3D PlanarSection::B() const
{
    return b;
}

void PlanarSection::DrawInputPolyline()
{

    glPointSize(3.0f);
    for (int c=0; c<sketch_pts.size(); ++c) {

        const QList <QVector2D> & samples = sketch_pts[c];

        glBegin(GL_POINTS);
        for (int i=0; i<samples.size(); ++i) {

            QVector3D v = p + (t * samples[i].x()) + (b * samples[i].y());
            glVertex3f(v.x(), v.y(), v.z());

        }
        glEnd();

    }
}

void PlanarSection::DrawCurve()
{

    //shows bridge edges for "repairing holes"
//    glBegin(GL_LINE_LOOP);
//    for (int i=0; i<poly.size(); ++i) {
//        QVector3D v = p + (t * poly[i].x()) + (b * poly[i].y());
//        glVertex3f(v.x(), v.y(), v.z());
//    }
//    glEnd();

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & samples = bez_curve[c].Samples();

        glBegin(GL_LINE_LOOP);
        for (int i=0; i<samples.size(); ++i) {
            QVector3D v = p + (t * samples[i].x()) + (b * samples[i].y());
            glVertex3f(v.x(), v.y(), v.z());
        }

        glEnd();

    }

}

void PlanarSection::DrawCurveControlPolygon()
{

    for (int i=0; i<bez_curve.size(); ++i) {

        const QList <QVector2D> & pts = bez_curve[i].Points();

        glBegin(GL_LINE_STRIP);
        for (int i=0; i<pts.size(); ++i) {
            const QVector3D v = p + (t * pts[i].x()) + (b * pts[i].y());
            glVertex3f(v.x(), v.y(), v.z());
        }
        glEnd();

    }

}


void PlanarSection::DrawCurveControlPoints(const float cam_width)
{

    for (int c=0; c<bez_curve.size(); ++c) {
        const QList <QVector2D> & pts = bez_curve[c].Points();

        glBegin(GL_LINES);
        for (int i=0; i<pts.size(); ++i) {

            const float s = (i == bez_curve[c].SelectedPoint()) ? 0.015f * cam_width : 0.0075f * cam_width;

            QVector3D v1, v2, v3, v4;

            if (i%3 == 0) {

                v1 = p + (t * (pts[i].x() - s)) + (b * (pts[i].y() - s));
                v2 = p + (t * (pts[i].x() + s)) + (b * (pts[i].y() + s));
                v3 = p + (t * (pts[i].x() - s)) + (b * (pts[i].y() + s));
                v4 = p + (t * (pts[i].x() + s)) + (b * (pts[i].y() - s));

                //X
                glVertex3f(v1.x(), v1.y(), v1.z());
                glVertex3f(v2.x(), v2.y(), v2.z());
                glVertex3f(v3.x(), v3.y(), v3.z());
                glVertex3f(v4.x(), v4.y(), v4.z());


            }
            else {

                v1 = p + (t * (pts[i].x() - s)) + (b * (pts[i].y() - s));
                v2 = p + (t * (pts[i].x() + s)) + (b * (pts[i].y() - s));
                v3 = p + (t * (pts[i].x() + s)) + (b * (pts[i].y() + s));
                v4 = p + (t * (pts[i].x() - s)) + (b * (pts[i].y() + s));

                //square
                glVertex3f(v1.x(), v1.y(), v1.z());
                glVertex3f(v2.x(), v2.y(), v2.z());
                glVertex3f(v2.x(), v2.y(), v2.z());
                glVertex3f(v3.x(), v3.y(), v3.z());
                glVertex3f(v3.x(), v3.y(), v3.z());
                glVertex3f(v4.x(), v4.y(), v4.z());
                glVertex3f(v4.x(), v4.y(), v4.z());
                glVertex3f(v1.x(), v1.y(), v1.z());
            }

        }
        glEnd();

    }

}




void PlanarSection::DrawCurveControlPointsHandleStyle(const float cam_width, const QVector3D cam_pos)
{
    QVector3D colour1(1.0, 0.0, .5);
    QVector3D colour2(1.0, 0.5, 0.8);;
    QVector3D colour_rad(0.2, 0.2, 0.6);

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & pts = bez_curve[c].Points();

        if (pts.empty()) {
            continue;
        }

        float s;
        //int numSegs;

        QVector3D lastPoint = p + (t * (pts[0].x())) + (b * (pts[0].y()));
        QVector3D point = p + (t * (pts[1].x())) + (b * (pts[1].y()));
        QVector3D nextPoint;

        glLineWidth(1.5);
        GLutils::glColor(colour1);

        int num_points = pts.size() - 1;

        glDisable(GL_DEPTH_TEST);
        glBegin(GL_LINES);
        for (int i=1; i<num_points; ++i) {
            if (radial && i == num_radial_points_per_sector) {
                GLutils::glColor(colour_rad);
            }

            nextPoint = p + (t * (pts[i+1].x())) + (b * (pts[i+1].y()));

            if (i%3 == 1) {
                glVertex3f(point.x(),point.y(),point.z());
                glVertex3f(lastPoint.x(),lastPoint.y(),lastPoint.z());
            }
            else if(i%3 == 2) {
                glVertex3f(point.x(),point.y(),point.z());
                glVertex3f(nextPoint.x(),nextPoint.y(),nextPoint.z());
            }
            lastPoint = point;
            point = nextPoint;

        }
        glEnd();
        glEnable(GL_DEPTH_TEST);
        glLineWidth(1);

        for (int i=0; i<num_points; ++i) {

            s = (i == bez_curve[c].SelectedPoint()) ? 0.005f * cam_width : 0.0035f * cam_width;

            //numSegs = (i == bez_curve[c].SelectedPoint()) ? 20 : 10;

            point = p + (t * (pts[i].x())) + (b * (pts[i].y()));

            if (i%3 == 0) {

                if (radial && i >= num_radial_points_per_sector) {
                    GLutils::glColor(colour_rad);
                    GLutils::DrawDisc(point, point + (cam_pos - point).normalized()*.2, 0.0f, 2*s);
                }
                else {
                    GLutils::glColor(colour2);

                    // for debug purposes
                    //glColor3f(0, float(i)/num_points, 1 - float(i)/num_points);

                    GLutils::DrawDisc(point, point + (cam_pos - point).normalized()*.2, 0.0f, 1.5*s);

                    GLutils::glColor(colour1);
                    GLutils::DrawDisc(point, point + (cam_pos - point).normalized()*.2, 1.5*s, 2*s);
                }


            }
            else {

                if (radial && i >= num_radial_points_per_sector) {
                    GLutils::glColor(colour_rad);
                }
                else {
                    GLutils::glColor(colour1);
                }

                // for debug purposes
                //glColor3f(0, float(i)/num_points, 1 - float(i)/num_points);

                GLutils::DrawDisc(point, point + (cam_pos - point).normalized()*.2, 0, s);
            }

        }

    }

}

void PlanarSection::DrawTNBFrame()
{

    glColor3f(1, 0, 0);
    GLutils::DrawArrowFixedLength(p, p+t, 0.5f);
    glColor3f(0, 1, 0);
    GLutils::DrawArrowFixedLength(p, p+n, 0.5f);
    glColor3f(0, 0, 1);
    GLutils::DrawArrowFixedLength(p, p+b, 0.5f);

    //const QVector3D c = GetPoint3D(centroid);
    /*
    const QVector3D c = p;

    glPushMatrix();

    glTranslatef(c.x(), c.y(), c.z());

    //draw tnb frame
    glBegin(GL_LINES);    

    //this draws at centroid
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(t.x(), t.y(), t.z());
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(n.x(), n.y(), n.z());
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(b.x(), b.y(), b.z());
    glEnd();

    glPopMatrix();
    */

}

void PlanarSection::DrawShadow()
{

    glColor3f(.75f, .75f, .75f);
    glBegin(GL_TRIANGLES);
    glNormal3f(0, 1, 0);
    for (int i=0; i<slab_vert.size(); ++i) {

        const QVector3D & v = slab_vert[i];

        glVertex3f(v.x() - qMax(float(v.y()), 0.0f), 0.0f, v.z() + qMax(float(v.y()), 0.0f)/2.0f);

    }
    glEnd();

}

void PlanarSection::DrawXZPlaneLine()
{

    //draw line where this plane and xz plane meet (if they do meet, i.e. plane is not along xz)
    const QVector3D ld = QVector3D::crossProduct(n, QVector3D(0, 1, 0)).normalized();

    //any point on this plane with y=0
    QVector3D lp = QVector3D(0, 0, 0);

    if (fabsf(n.x()) > fabsf(n.z())) {
        lp.setX(QVector3D::dotProduct(n, p) / n.x());
    }
    else {
        lp.setZ(QVector3D::dotProduct(n, p) / n.z());
    }

    glBegin(GL_LINES);
    glVertex3f(lp.x() - ld.x() * 100.0f, lp.y() - ld.y()* 100.0f, lp.z() - ld.z()* 100.0f);
    glVertex3f(lp.x() + ld.x() * 100.0f, lp.y() + ld.y()* 100.0f, lp.z() + ld.z()* 100.0f);
    glEnd();
    //qDebug() << lp << ld;

}

void PlanarSection::DrawDeadzone(const QVector3D & slot_start, const QVector3D & slot_end, const float deadzone_radius)
{

    const QVector2D p1 = GetPoint2D(slot_start);
    const QVector2D p2 = GetPoint2D(slot_end);

    for (int j=0; j<2; ++j) {

        if (j == 0) {
            glLineWidth(1.0f);
            glBegin(GL_LINE_LOOP);
            glColor3f(0.75f, 0.75f, 0.75f);
        }
        else {

            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glBegin(GL_TRIANGLE_FAN);
            glColor4f(0.75f, 0.75f, 0.75f, 0.5f);
        }

        for (int i=0; i<360; i+=5) {

            const float f = float(i) / 180.0f * 3.14159f;

            QVector2D p = QVector2D(cosf(f), sinf(f)) * deadzone_radius;

            if (QVector2D::dotProduct(p2-p1, p) >= 0.0f) {
                p += p2;
            }
            else {
                p += p1;
            }

            const QVector3D p3 = GetPoint3D(p);
            glVertex3f(p3.x(), p3.y(), p3.z());

        }
        glEnd();

    }

    glLineWidth(1.0f);
    glDisable(GL_BLEND);

}

void PlanarSection::DrawConvexHull(const QList <QVector3D> & hull, const bool draw_green)
{

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glLineWidth(3.0f);
    if (draw_green) {
        glColor4f(0.0f, 0.25f, 0.0f, 0.75f);
    }
    else {
        glColor4f(0.25f, 0.0f, 0.0f, 0.75f);
    }
    glBegin(GL_LINE_LOOP);
    for (int i=0; i<hull.size(); ++i) {
        glVertex3f(hull[i].x(), hull[i].y(), hull[i].z());
    }
    glEnd();

    glLineWidth(1.0f);
    if (draw_green) {
        glColor4f(0.5f, 0.75f, 0.5f, 0.5f);
    }
    else {
        glColor4f(0.75f, 0.5f, 0.5f, 0.5f);
    }
    glBegin(GL_TRIANGLE_FAN);
    for (int i=0; i<hull.size(); ++i) {
        glVertex3f(hull[i].x(), hull[i].y(), hull[i].z());
    }
    glEnd();

    glDisable(GL_BLEND);

}

void PlanarSection::Draw()
{

    DrawSketch();
    DrawCurve();
    DrawTris();
    DrawTNBFrame();

}

void PlanarSection::DrawForPicking(const int pickval)
{

    GLutils::SetPickColor(pickval);
    //DrawTris();
    DrawSlab();

}

void PlanarSection::Update(const int c, const double error_tolerance)
{

    //const double error_tolerance = 0.04;
    //const double error_tolerance = 0.02;
    //const double error_tolerance = 0.01;    
    if (sketch_pts[c].size() < 2) {
        return;
    }

    QList <QVector2D> fit_curve;
    bez_fit.FitCurve(sketch_pts[c], error_tolerance, fit_curve);

    bez_curve[c].ClearPoints();

    for (int i=0; i+3<fit_curve.size(); i+=4) {
        bez_curve[c].AddPoint(fit_curve[i]);
        bez_curve[c].AddPoint(fit_curve[i+1]);
        bez_curve[c].AddPoint(fit_curve[i+2]);
    }

    //we do a hermite-closing
    //QVector2D tan_end = (fit_curve[fit_curve.size()-1] - fit_curve[fit_curve.size()-2]).normalized();
    QVector2D tan_end;
    if (sketch_pts[c].size() > 5) {
        tan_end = (sketch_pts[c][sketch_pts[c].size()-1] - sketch_pts[c][sketch_pts[c].size()-5]).normalized();
    }
    else {
        tan_end = (sketch_pts[c][sketch_pts[c].size()-1] - sketch_pts[c][sketch_pts[c].size()-2]).normalized();
    }
    QVector2D tan_start = (fit_curve[0] - fit_curve[1]).normalized();

    float dist = (fit_curve.last() - fit_curve.first()).length();

    //tan_end *= dist * 0.25;
    //tan_start *= dist * 0.25;
    tan_end *= dist * 0.35f;
    tan_start *= dist * 0.35f;

    bez_curve[c].AddPoint(fit_curve.last());
    bez_curve[c].AddPoint(fit_curve.last() + tan_end);
    bez_curve[c].AddPoint(fit_curve.first() + tan_start);
    bez_curve[c].AddPoint(fit_curve.first());

    bez_curve[c].SetClosed(true);

    UpdateCurveTrisSlab();

}

void PlanarSection::UpdateCurveTrisSlab()
{      

    for (int c=0; c<bez_curve.size(); ++c) {
        bez_curve[c].UpdateSamples();
    }

    if (bez_curve.size() == 1) {

        //we triangulate only 1 curve using ear method
        Triangulate2::Process(bez_curve[0].Samples(), tris, poly);

    }
    else {

        //we triangulate multiple curves (indicating holes) using triangle lib
        /*
        QList <QList <QVector2D> > samples;
        for (int c=0; c<bez_curve.size(); ++c) {
            if (bez_curve[c].Samples().size() > 5) {
                samples.push_back(bez_curve[c].Samples());
            }
        }
        Triangulate::Process(samples, tris);
        */
        QList <QList <QVector2D> > contours;
        for (int c=0; c<bez_curve.size(); ++c) {
            //if (bez_curve[c].Samples().size() > 5) {
            contours.push_back(bez_curve[c].Samples());
            //}
        }
        Triangulate2::Process(contours, tris, poly);

    }

    ComputeCentroidAndArea();

    if (!editing_sketch) {
        ComputeSlab();
    }   

    /*
    if (slab_disp_list > 0) {
        glDeleteLists(1, slab_disp_list);
        slab_disp_list = 0;
    }

    if (slab_curves_disp_list > 0) {
        glDeleteLists(1, slab_curves_disp_list);
        slab_curves_disp_list = 0;
    }

    update_slab_disp_list = true;
    update_slabcurves_disp_list = true;
        */

}

void PlanarSection::SaveHeaderToFile(const QList <PlanarSection> & sections, QTextStream & ofs)
{
    ofs << "HeaderParameters 1\n";
    ofs << "FlatFab " << sections.size() << "\n";
    //ofs << "RotateDuration " << rot_dur << "\n";
    //ofs << "SlabThickness " << slab_thick << "\n";
}

int PlanarSection::LoadHeaderFromFile(QTextStream & ifs)
{

    QStringList line;

    line = ifs.readLine().split(" ");
    if (line.size() != 2) {
        qDebug() << "PlanarSection::LoadHeaderFromFile - Invalid format";
        return 0;
    }
    if (QString::compare(line.first(), "HeaderParameters") != 0) {
        qDebug() << "PlanarSection::LoadHeaderFromFile - Invalid format";
        return 0;
    }

    line = ifs.readLine().split(" ");
    if (QString::compare(line.first(), "PlaneSketch") == 0 || QString::compare(line.first(), "FlatFab") == 0) {
        return line.last().toInt();
    }
    else {
        qDebug() << "PlanarSection::LoadHeaderFromFile - Missing important ";
        return 0;
    }

}

void PlanarSection::SaveToFile(QTextStream & ofs)
{

    //P and TNB frame
    ofs << "P " << p.x() << " " << p.y() << " " << p.z() << "\n";
    ofs << "T " << t.x() << " " << t.y() << " " << t.z() << "\n";
    ofs << "N " << n.x() << " " << n.y() << " " << n.z() << "\n";
    ofs << "B " << b.x() << " " << b.y() << " " << b.z() << "\n";

    //sketch points    
    ofs << "BoundaryCurves " << sketch_pts.size() << "\n";
    for (int j=0; j<sketch_pts.size(); ++j) {

        ofs << "SketchPoints " << sketch_pts[j].size() << "\n";
        for (int i=0; i<sketch_pts[j].size(); ++i) {
            ofs << sketch_pts[j][i].x() << " " << sketch_pts[j][i].y() << "\n";
        }

        //bezier curve control points
        const QList <QVector2D> & bez_pts = bez_curve[j].Points();
        ofs << "BezierPoints " << bez_pts.size() << "\n";
        for (int i=0; i<bez_pts.size(); ++i) {
            ofs << bez_pts[i].x() << " " << bez_pts[i].y() << "\n";
        }

    }

}

void PlanarSection::LoadFromFile(QTextStream & ifs)
{

    QStringList line;

    line = ifs.readLine().split(" ");
    p = QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat());

    line = ifs.readLine().split(" ");
    t = QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat());

    line = ifs.readLine().split(" ");
    n = QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat());

    line = ifs.readLine().split(" ");
    b = QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat());

    //determine what to do here... if we see "boundarycurves", holes are supported
    //if not, there are no holes, just a single boundary contour
    QStringList curves_line = ifs.readLine().split(" ");

    if (QString::compare(curves_line.first(), "BoundaryCurves") == 0) {

        const int nCurves = curves_line.last().toInt();

        for (int j=0; j<nCurves; ++j) {

            if (j > 0) {
                AddNewCurve();
            }

            //sketch points
            const int nSketchPts = ifs.readLine().split(" ").last().toInt();
            for (int i=0; i<nSketchPts; ++i) {
                line = ifs.readLine().split(" ");
                sketch_pts[j].push_back(QVector2D(line[0].toFloat(), line[1].toFloat()));
            }

            //bezier curve control points
            const int nBezPts = ifs.readLine().split(" ").last().toInt();
            for (int i=0; i<nBezPts; ++i) {
                line = ifs.readLine().split(" ");
                bez_curve[j].AddPoint(QVector2D(line[0].toFloat(), line[1].toFloat()));
            }
        }

    }
    else {

        //old way, before holes were supported
        //const int nSketchPts = ifs.readLine().split(" ").last().toInt();
        const int nSketchPts = curves_line.last().toInt();
        for (int i=0; i<nSketchPts; ++i) {
            line = ifs.readLine().split(" ");
            sketch_pts[0].push_back(QVector2D(line[0].toFloat(), line[1].toFloat()));
        }

        //bezier curve control points
        const int nBezPts = ifs.readLine().split(" ").last().toInt();
        for (int i=0; i<nBezPts; ++i) {
            line = ifs.readLine().split(" ");
            bez_curve[0].AddPoint(QVector2D(line[0].toFloat(), line[1].toFloat()));
        }

    }

}

void PlanarSection::ComputePacking(const QList <PlanarSection> & sections, const double width_height_ratio, QVector <QVector2D> & pack_pts, QVector2D & bbox)
{

    qDebug() << "Packing" << sections.size() << "sections into rectangle of ratio" << width_height_ratio << "...";

    //First Fit Decreasing (FFD) strategy, operates by first sorting the items to be inserted in decreasing order by their sizes,
    //and then inserting each item into the first bin in the list with sufficient remaining space.

    //1.  get each rectangle's dimensions and compute the total area
    QList <QPair <int, float> > pieces;

    float total_area = 0.0f;
    QVector <QVector2D> bb = QVector <QVector2D> (sections.size());

    for (int i=0; i<sections.size(); ++i) {

        QVector2D min_v, max_v;
        sections[i].GetBoundingBox2D(min_v, max_v);
        bb[i] = QVector2D(max_v.x()-min_v.x(), max_v.y()-min_v.y());

        total_area += bb[i].x() * bb[i].y();

        QPair <int, float> new_piece = QPair <int, float> (i, qMax(bb[i].x(), bb[i].y()));
        pieces.push_back(new_piece);

    }

    //1.5 sort by size
    for (int i=0; i<pieces.size(); ++i) {
        for (int j=i+1; j<pieces.size(); ++j) {
            if (pieces[i].second < pieces[j].second) {
                qSwap(pieces[i], pieces[j]);
            }
        }
    }

    //2.  make a rectangle (double the area of this square) with boolean occupancies
    pack_pts = QVector <QVector2D> (sections.size(), QVector2D(0, 0));

    /*
    const int rw = sqrtf(total_area) * 2.0f;
    const int rh = sqrtf(total_area);

    QVector <QVector <bool> > occ = QVector <QVector <bool> > (rw);
    for (int i=0; i<rw; ++i) {
        occ[i] = QVector <bool> (rh, false);
    }
    */

    bbox = QVector2D(0, 0);

    //2.  start packing (using the axis-aligned bounding box)
    for (int i=0; i<sections.size(); ++i) {

        const int p1 = pieces[i].first;

        //2a.  find spot within existing bounding area
        bool found = false;

        for (int x=0; x+bb[p1].x()<bbox.x(); ++x) {
            for (int y=0; y+bb[p1].y()<bbox.y(); ++y) {

                pack_pts[p1] = QVector2D(x, y);

                bool collide = false;
                for (int j=0; j<i; ++j) { //check all pieces already placed for collision

                    const int p2 = pieces[j].first;

                    if (pack_pts[p1].x() < pack_pts[p2].x() + bb[p2].x() &&
                            pack_pts[p1].x() + bb[p1].x() > pack_pts[p2].x() &&
                            pack_pts[p1].y() < pack_pts[p2].y() + bb[p2].y() &&
                            pack_pts[p1].y() + bb[p1].y() > pack_pts[p2].y()) {
                        collide = true;
                        break;
                    }
                }

                if (!collide) {
                    found = true;
                    break;
                }

            }

            if (found) {
                break;
            }

        }

        //2b.  put spot in if we can
        if (!found) {

            //concatenate this one, expand the box
            //try either more to the right, or more down
            float new_right = bbox.x() + bb[p1].x();
            float new_up = bbox.y() + bb[p1].y();

            if (new_right / width_height_ratio < new_up) {
                pack_pts[p1] = QVector2D(bbox.x(), 0);
            }
            else {
                pack_pts[p1] = QVector2D(0, bbox.y());
            }

        }

        //qDebug() << "Expanding square" << bb[p1].x() << bb[p1].y() << "at spot" << pack_pts[p1].x() << pack_pts[p1].y() << "of size" << bb[p1].x() << bb[p1].y();

        bbox.setX(qMax(bbox.x(), pack_pts[p1].x() + bb[p1].x()));
        bbox.setY(qMax(bbox.y(), pack_pts[p1].y() + bb[p1].y()));

    }

    qDebug() << "Completed packing inside rectangle of dimensions" << bbox.x() << bbox.y();

}

void PlanarSection::SaveToDXF(const QList <PlanarSection> & sections, QTextStream & ofs, const double metres_per_unit, const int dpi, const double width_height_ratio, const bool use_numeric_labels, const QList <QString> & numeric_labels)
{

    //so say we have this 2x2 box... we want each unit to be 1 inch
    //the scaling factor is...
    //multiply by dpi.... then multiply by factor with inches
    //const double scale_inches = metres_per_unit / 0.0254;
    //const double scale = double(dpi) * scale_inches;
    //const double scale = scale_inches;
    const double scale = metres_per_unit * 1000.0 * 1000.0; //format calls for mm I believe

    //1. TODO pack each section's shape into the rectangular box
    QVector <QVector2D> pack_pts;
    QVector2D bbox;
    ComputePacking(sections, width_height_ratio, pack_pts, bbox);

    //if we are writing out the planar sections with labels, we require a connectivity matrix
    QVector <QVector <QString> > graph_labels;
    if (use_numeric_labels) {
        ComputeSlitNumericLabels(sections, numeric_labels, graph_labels);
    }

    //2. Write out the DXF data
    //2a.  COMMENTS
    ofs << "999\n";
    ofs << "DXF created using FlatFab\n";
    //2b.  HEADER
    ofs << "0\n";
    ofs << "SECTION\n";
    ofs << "2\n";
    ofs << "HEADER\n";
    ofs << "9\n";
    ofs << "$ACADVER\n";
    ofs << "1\n";
    ofs << "AC1006\n";
    ofs << "9\n";
    //ofs << "$MEASUREMENT\n";
    //ofs << "70\n";
    //ofs << "1\n";
    //ofs << "9\n";
    //ofs << "$INSUNITS\n";
    //ofs << "70\n";
    //ofs << "6\n"; //NOTE SHOULD BE 4 for mm, 6 is metres
    //ofs << "9\n";
    //ofs << "$LUNITS\n";
    //ofs << "70\n";
    //ofs << "6\n";
    //ofs << "9\n";
    ofs << "$INSBASE\n";
    ofs << "10\n";
    ofs << "0.0\n";
    ofs << "20\n";
    ofs << "0.0\n";
    ofs << "30\n";
    ofs << "0.0\n";
    ofs << "9\n";
    ofs << "$EXTMIN\n";
    ofs << "10\n";
    ofs << "0.0\n";
    ofs << "20\n";
    ofs << "0.0\n";
    ofs << "9\n";
    ofs << "$EXTMAX\n";
    ofs << "10\n";
    ofs << "1000.0\n";
    //ofs << bbox.x() * scale << "\n";
    ofs << "20\n";
    ofs << "1000.0\n";
    //ofs << bbox.y() * scale  << "\n";
    ofs << "0\n";
    ofs << "ENDSEC\n";
    //2c. TABLES
    ofs << "0\n";
    ofs << "SECTION\n";
    ofs << "2\n";
    ofs << "TABLES\n";
    ofs << "0\n";
    ofs << "TABLE\n";
    ofs << "2\n";
    ofs << "LTYPE\n";
    ofs << "70\n";
    ofs << "1\n";
    ofs << "0\n";
    ofs << "LTYPE\n";
    ofs << "2\n";
    ofs << "CONTINUOUS\n";
    ofs << "70\n";
    ofs << "64\n";
    ofs << "3\n";
    ofs << "Solid line\n";
    ofs << "72\n";
    ofs << "65\n";
    ofs << "73\n";
    ofs << "0\n";
    ofs << "40\n";
    ofs << "0.000000\n";
    ofs << "0\n";
    ofs << "ENDTAB\n";
    ofs << "0\n";
    ofs << "TABLE\n";
    ofs << "2\n";
    ofs << "LAYER\n";
    ofs << "70\n";
    ofs << "6\n";
    ofs << "0\n";
    ofs << "LAYER\n";
    ofs << "2\n";
    ofs << "1\n";
    ofs << "70\n";
    ofs << "64\n";
    ofs << "62\n";
    ofs << "7\n";
    ofs << "6\n";
    ofs << "CONTINUOUS\n";
    ofs << "0\n";
    ofs << "LAYER\n";
    ofs << "2\n";
    ofs << "2\n";
    ofs << "70\n";
    ofs << "64\n";
    ofs << "62\n";
    ofs << "7\n";
    ofs << "6\n";
    ofs << "CONTINUOUS\n";
    ofs << "0\n";
    ofs << "ENDTAB\n";
    ofs << "0\n";
    ofs << "TABLE\n";
    ofs << "2\n";
    ofs << "STYLE\n";
    ofs << "70\n";
    ofs << "0\n";
    ofs << "0\n";
    ofs << "ENDTAB\n";
    ofs << "0\n";
    ofs << "ENDSEC\n";
    //2d. BLOCKS
    ofs << "0\n";
    ofs << "SECTION\n";
    ofs << "2\n";
    ofs << "BLOCKS\n";
    ofs << "0\n";
    ofs << "ENDSEC\n";
    //2e. ENTITIES (where geometry is specified)
    ofs << "0\n";
    ofs << "SECTION\n";
    ofs << "2\n";
    ofs << "ENTITIES\n";

    for (int i=0; i<sections.size(); ++i) {

        QVector2D min_v, max_v;
        sections[i].GetBoundingBox2D(min_v, max_v);

        for (int c=0; c<sections[i].bez_curve.size(); ++c) {

            const QList <QVector2D> & samples = sections[i].bez_curve[c].Samples();

            if (samples.empty()) {
                continue;
            }

            for (int j=0; j<samples.size(); ++j) {
                const int i1 = j;
                const int i2 = (j+1) % samples.size();

                QVector3D p1 = (samples[i1] + pack_pts[i] - min_v) * scale;
                QVector3D p2 = (samples[i2] + pack_pts[i] - min_v) * scale;

                DXFWriteLine(ofs, p1, p2, 1, 1);

            }

        }

        //draw the slotlines
        for (int j=0; j<sections.size(); ++j) {

            if (i == j) {
                continue;
            }

            QList <QVector2D> slots_i;
            QList <QList <QVector2D> > slots_i_rect;
            sections[i].GetSlots(sections[j], slots_i, slots_i_rect);

            for (int k=0; k<slots_i_rect.size(); k+=5) {


                for (int l=1; l<slots_i_rect[k].size(); ++l) {
                    const int i1 = l-1;
                    const int i2 = l;

                    QVector3D p1 = (slots_i_rect[k][i1] + pack_pts[i] - min_v) * scale;
                    QVector3D p2 = (slots_i_rect[k][i2] + pack_pts[i] - min_v) * scale;

                    DXFWriteLine(ofs, p1, p2, 0, 0);

                }

                /*
                if (use_numeric_labels) {
                    //ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#ff0000;stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
                }
                else {
                    //ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#" << GLutils::ColorByIndexStr(j) << ";stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
                }
                */

                //draw numeric labels
                if (use_numeric_labels && !slots_i_rect[k].empty()) {

                    QVector2D dir = (slots_i_rect[k][2] - slots_i_rect[k][1]).normalized();
                    QVector3D p1 = ((slots_i_rect[k][2] + slots_i_rect[k][3]) * 0.5f + pack_pts[i] - min_v) * scale + dir * 20.0f;
                    DXFWriteText(ofs, p1, graph_labels[i][j], 1, 1, 100.0f);
                    //ofs << "<text x=\"" << p.x() << "\" y=\"" << p.y() << "\" text-anchor=\"middle\" ";
                    //ofs << "style=\"font-size:16px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill:#0000ff;fill-opacity:1;stroke:none;font-family:Sans\">\n";
                    //ofs << graph_labels[i][j] << "\n";
                    //ofs << "</text>\n";

                }

            }

        }

    }

    ofs << "0\n";
    ofs << "ENDSEC\n";
    //2f. EOF
    ofs << "0\n";
    ofs << "EOF\n";
}

void PlanarSection::DXFWriteLine(QTextStream & ofs, const QVector3D & p1, const QVector3D & p2, const int layer_index, const int colour_index)
{

    ofs << "0\n";
    ofs << "LINE\n";
    ofs << "8\n";
    ofs << layer_index << "\n";
    ofs << "62\n";
    ofs << colour_index << "\n";
    ofs << "10\n";
    ofs << p1.x() << "\n";
    ofs << "20\n";
    ofs << p1.y() << "\n";
    ofs << "30\n";
    ofs << p1.z() << "\n";
    ofs << "11\n";
    ofs << p2.x() << "\n";
    ofs << "21\n";
    ofs << p2.y() << "\n";
    ofs << "31\n";
    ofs << p2.z() << "\n";

}

void PlanarSection::DXFWriteText(QTextStream & ofs, const QVector3D & p1, const QString & text, const int layer_index, const int colour_index, const float height)
{

    ofs << "0\n";
    ofs << "TEXT\n";
    ofs << "8\n";
    ofs << layer_index << "\n";
    ofs << "62\n";
    ofs << colour_index << "\n";
    ofs << "10\n";
    ofs << p1.x() << "\n";
    ofs << "20\n";
    ofs << p1.y() << "\n";
    ofs << "30\n";
    ofs << p1.z() << "\n";
    ofs << "40\n";
    ofs << height << "\n";
    ofs << "1\n";
    ofs << text << "\n";

}

void PlanarSection::ComputeSlitNumericLabels(const QList <PlanarSection> & sections, const QList <QString> & numeric_labels, QVector <QVector <QString> > & graph_labels)
{

    //if we are writing out the planar sections with labels, we require a connectivity matrix
    int current_numeric_label = 1;
    QVector <QVector <bool> > graph;

    ComputeIntersectionGraph(sections, graph);

    graph_labels.resize(graph.size());
    for (int i=0; i<graph.size(); ++i) {
        graph_labels[i].resize(graph[i].size());
    }

    for (int i=0; i<graph.size(); ++i) {
        for (int j=i+1; j<graph[i].size(); ++j) {
            if (graph[i][j]) {
                if (j < numeric_labels.size()) {
                    graph_labels[i][j] = numeric_labels[j];
                    graph_labels[j][i] = graph_labels[i][j];
                }
                else {
                    graph_labels[i][j] = QString::number(current_numeric_label);
                    graph_labels[j][i] = graph_labels[i][j];
                    ++current_numeric_label;
                }
            }
        }
    }

}

void PlanarSection::SaveToSVG(const QList <PlanarSection> & sections, QTextStream & ofs, const double metres_per_unit, const int dpi, const double width_height_ratio, const bool use_numeric_labels, const QList <QString> & numeric_labels)
{

    //so say we have this 2x2 box... we want each unit to be 1 inch
    //the scaling factor is...
    //multiply by dpi.... then multiply by factor with inches
    const double scale_inches = metres_per_unit / 0.0254;
    const double scale = double(dpi) * scale_inches;

    //1. TODO pack each section's shape into the rectangular box
    QVector <QVector2D> pack_pts;
    QVector2D bbox;
    ComputePacking(sections, width_height_ratio, pack_pts, bbox);

    //if we are writing out the planar sections with labels, we require a connectivity matrix
    QVector <QVector <QString> > graph_labels;
    if (use_numeric_labels) {
        ComputeSlitNumericLabels(sections, numeric_labels, graph_labels);
    }

    //2. Write out the SVG data
    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    ofs << "<svg\n";
    ofs << "xmlns:svg=\"http://www.w3.org/2000/svg\"\n";
    ofs << "xmlns=\"http://www.w3.org/2000/svg\"\n";
    ofs << "version=\"1.0\"\n";
    ofs << "width=\"" << bbox.x() * scale_inches << "in\"\n";
    ofs << "height=\"" << bbox.y() * scale_inches << "in\"\n";
    ofs << "viewbox=\"0 0 " << bbox.x() * scale_inches << " " << bbox.y() * scale_inches << "\"\n";
    ofs << "id=\"svg2\">\n";
    ofs << "<defs\n";
    ofs << "id=\"defs4\" />\n";   

    //just draw the contours (pretty easy)
    for (int i=0; i<sections.size(); ++i) {

        QVector2D min_v, max_v;
        sections[i].GetBoundingBox2D(min_v, max_v);

        ofs << "<g id=\"layer" << i << "\">\n";

        for (int c=0; c<sections[i].bez_curve.size(); ++c) {

            const QList <QVector2D> & samples = sections[i].bez_curve[c].Samples();

            if (samples.empty()) {
                continue;
            }


            //ofs << "<g id=\"layer" << i << "curve" << c << "\">\n";
            ofs << "<path d=\"";

            QVector2D p = samples[0] + pack_pts[i] - min_v;
            p *= scale;
            ofs << "M " << p.x() << " " << p.y();

            for (int j=1; j<samples.size(); ++j) {
                p = samples[j] + pack_pts[i] - min_v;
                p *= scale;
                ofs << " L " << p.x() << " " << p.y();
            }

            p = samples[0] + pack_pts[i] - min_v;
            p *= scale;
            ofs << " L " << p.x() << " " << p.y();

            QColor color;
            color.setRgb(i);

            ofs << "\" id=\"layer" << i << "curve" << c << "\"\n";
            if (use_numeric_labels) {
                ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#ff0000;stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
            }
            else {
                ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#" << GLutils::ColorByIndexStr(i) << ";stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
            }

        }

        //draw the slotlines
        for (int j=0; j<sections.size(); ++j) {

            if (i == j) {
                continue;
            }

            QList <QVector2D> slots_i;
            QList <QList <QVector2D> > slots_i_rect;
            sections[i].GetSlots(sections[j], slots_i, slots_i_rect);

            //if (!slots_i.empty()) {
            //    qDebug() << "Section" << i << "intersects" << j << "?" << !slots_i.empty() << slots_i.size();
            //}

            for (int k=0; k<slots_i_rect.size(); k+=5) {
                //ofs << "<g id=\"layerslot" << i << "to" << j << "at" << k << "\">\n";
                ofs << "<path d=\"";

                QVector2D p = (slots_i_rect[k][0] + pack_pts[i] - min_v) * scale;
                ofs << "M " << p.x() << " " << p.y();

                for (int l=1; l<slots_i_rect[k].size(); ++l) {
                    p = (slots_i_rect[k][l] + pack_pts[i] - min_v) * scale;
                    ofs << " L " << p.x() << " " << p.y();
                }

                ofs << "\" id=\"layer" << i << "slot" << j << "\"\n";
                if (use_numeric_labels) {
                    ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#ff0000;stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
                }
                else {
                    ofs << "style=\"fill:none;fill-rule:evenodd;stroke:#" << GLutils::ColorByIndexStr(j) << ";stroke-width:0.02in;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" />\n";
                }
                //ofs << "</g>\n";

                //draw numeric labels
                if (use_numeric_labels && !slots_i_rect[k].empty()) {

                    QVector2D dir = (slots_i_rect[k][2] - slots_i_rect[k][1]).normalized();
                    p = ((slots_i_rect[k][2] + slots_i_rect[k][3]) * 0.5f + pack_pts[i] - min_v) * scale + dir * 20.0f;
                    ofs << "<text x=\"" << p.x() << "\" y=\"" << p.y() << "\" text-anchor=\"middle\" ";
                    ofs << "style=\"font-size:16px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill:#0000ff;fill-opacity:1;stroke:none;font-family:Sans\">\n";                  
                    ofs << graph_labels[i][j] << "\n";
                    ofs << "</text>\n";

                }

            }

        }

        ofs << "</g>\n";

    }

    ofs << "</svg>\n";

}

void PlanarSection::GetSlots(const PlanarSection & other, QList <QVector2D> & my_slots, QList <QList <QVector2D> > & my_slots_rect) const
{

    //qDebug() << "PlanarSection::GetSlots() - calling with " << index << other_index;
    my_slots.clear();
    my_slots_rect.clear();

    //1.  get the slotline and ensure it has nonzero length (a ray is properly defined)
    QVector2D p0_slot_start, p0_slot_end;
    QVector2D p1_slot_start, p1_slot_end;
    if (!GetIntersectionSlotLine(other, p0_slot_start, p0_slot_end)) {
        return;
    }   

    p1_slot_start = other.GetPoint2D(GetPoint3D(p0_slot_start));
    p1_slot_end = other.GetPoint2D(GetPoint3D(p0_slot_end));

    //2.  determine the contour intersection points for THIS section and the OTHER section
    //qDebug() << "PlanarSection::GetSlots()" << p0_slot_start << p0_slot_end;
    QList <QVector2D> p0_isecs; //note: 2D coordinates are in that plane's local coordinate system
    QList <QVector2D> p1_isecs;

    QVector2D p0_slot_dir = (p0_slot_end - p0_slot_start).normalized();
    QVector2D p0_slot_ortho = QVector2D(-p0_slot_dir.y(), p0_slot_dir.x());
    QVector2D p1_slot_dir = (p1_slot_end - p1_slot_start).normalized();

    ComputeContourLineIntersects(p0_slot_start, p0_slot_dir, p0_isecs);
    other.ComputeContourLineIntersects(p1_slot_start, p1_slot_dir, p1_isecs);

    //2a.  make sure there's intervals for both p0 and p1 (the THIS and OTHER plane)
    //     if there's not, silently return
    if (p0_isecs.empty() || p1_isecs.empty()) {
        //qDebug() << "PlanarSection::GetSlots() - was empty, no big deal" << p0_isecs << p1_isecs;
        return;
    }

    //2b.  if there's an odd number of intersection points, we've got an issue (or we just tangentially grazed a point)
    //     (for now, we abort when this happens)
    /*
    if ((p0_isecs.size() % 2) == 1 || (p1_isecs.size() % 2) == 1) {
        qDebug() << "PlanarSection::GetSlots() - odd number of intersection points along ray encountered, aborting" << p0_isecs.size() << p1_isecs.size();
        qDebug() << p0_isecs;
        qDebug() << p1_isecs;

        return;
    }
    */

    //3.  sort these points, using the common ray in 3D
    QVector3D slot_start_3d = GetPoint3D(p0_slot_start);
    QVector3D slot_end_3d = GetPoint3D(p0_slot_end);
    QVector3D slot_dir_3d = (slot_end_3d - slot_start_3d).normalized();

    QList <QVector3D> p0_isecs_3d;
    QList <QVector3D> p1_isecs_3d;

    for (int i=0; i<p0_isecs.size(); ++i) {
        p0_isecs_3d.push_back(GetPoint3D(p0_isecs[i]));
    }
    for (int i=0; i<p1_isecs.size(); ++i) {
        p1_isecs_3d.push_back(other.GetPoint3D(p1_isecs[i]));
    }

    GLutils::SortPointsAlongDirection3D(slot_dir_3d, p0_isecs_3d);
    GLutils::SortPointsAlongDirection3D(slot_dir_3d, p1_isecs_3d);

    //DEBUGGGGGGGGGGGGGGGGGGGGGG//DEBUGGGGGGGGGGGGGGGGGGGGGG//DEBUGGGGGGGGGGGGGGGGGGGGGG
    /*
    qDebug() << "PlanarSection::GetSlots() - HELLO!!!" << p0_isecs << p1_isecs;
    QList <float> p0l, p1l;
    for (int i=0; i<p0_isecs_3d.size(); i+=2) { //iterate over THIS's intervals

        const float p0_start = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i]);
        const float p0_end = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i+1]);
        p0l.push_back(p0_start);
        p0l.push_back(p0_end);
    }

    for (int j=0; j<p1_isecs_3d.size(); j+=2) {//iterate over OTHER's intervals, see if there's some overlap

            const float p1_start = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j]);
            const float p1_end = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j+1]);
            p1l.push_back(p1_start);
            p1l.push_back(p1_end);
    }
    qDebug() << "p0l" << p0l;
    qDebug() << "p1l" << p1l;
    */
    //DEBUGGGGGGGGGGGGGGGGGGGGGG//DEBUGGGGGGGGGGGGGGGGGGGGGG//DEBUGGGGGGGGGGGGGGGGGGGGGG

    //4.  with the sorted points, and since these are closed bounary intersection points, each sequential pair defines
    //    an interior region.  so we want to know the intersections of those interior regions/intervals, because at each
    //    one we will want to create a slot rectangle.
    //    also, we only care about doing this for the THIS plane, not the OTHER plane

    //qDebug() << "PlanarSection::GetSlots() - made it here" << p0_isecs_3d.size() << p1_isecs_3d.size();

    //bool found_anything = false;

    for (int i=0; i+1<p0_isecs_3d.size(); i+=2) { //iterate over THIS's intervals

        const float p0_start = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i]);
        const float p0_end = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i+1]);

        for (int j=0; j+1<p1_isecs_3d.size(); j+=2) {//iterate over OTHER's intervals, see if there's some overlap

            const float p1_start = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j]);
            const float p1_end = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j+1]);

            //qDebug() << "THE VALUES" << p0_start << p0_end << p1_start << p1_end;
            //so there's probably some set of conditions we're looking for here...
            //a) p1's startpoint is within p0's start and end points
            //b) p1's endpoint is within p0's start and end points

            bool do_add_slot = false;
            QVector3D slot_mid;
            QVector2D slot_start_2d;
            QVector2D slot_mid_2d;

            if (p0_start < p1_start && p0_end > p1_start) { //condition a)

                do_add_slot = true;

                //add the slot
                slot_mid = (p1_isecs_3d[j] + p0_isecs_3d[i+1]) * 0.5f;
                slot_start_2d = GetPoint2D(p0_isecs_3d[i+1]);
                slot_mid_2d = GetPoint2D(slot_mid);

            }
            else if (p0_start < p1_end && p0_end > p1_end) { //condition b)

                do_add_slot = true;

                //add the slot
                slot_mid = (p0_isecs_3d[i] + p1_isecs_3d[j+1]) * 0.5f;
                slot_start_2d = GetPoint2D(p0_isecs_3d[i]);
                slot_mid_2d = GetPoint2D(slot_mid);
            }

            if (do_add_slot) { //one of the conditions was satisfied, and the values were set

                my_slots.push_back(slot_start_2d);
                my_slots.push_back(slot_mid_2d);

                //add a slot rectangle
                QList <QVector2D> slot_rect1;
                slot_rect1.push_back(slot_start_2d + p0_slot_ortho * other.SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_start_2d - p0_slot_ortho * other.SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_mid_2d - p0_slot_ortho * other.SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_mid_2d + p0_slot_ortho * other.SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_start_2d + p0_slot_ortho * other.SlabThickness() * 0.5f);
                my_slots_rect.push_back(slot_rect1);

                //found_anything = true;

            }

        }

    }

    /*
    if (!found_anything) {
        qDebug() << "PlanarSection::GetSlots() - WARNING! Seemed good, but didn't find anything" << p0_isecs.size() << p1_isecs.size();
        QList <float> p0l, p1l;
        for (int i=0; i<p0_isecs_3d.size(); i+=2) { //iterate over THIS's intervals

            const float p0_start = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i]);
            const float p0_end = QVector3D::dotProduct(slot_dir_3d, p0_isecs_3d[i+1]);
            p0l.push_back(p0_start);
            p0l.push_back(p0_end);
        }

        for (int j=0; j<p1_isecs_3d.size(); j+=2) {//iterate over OTHER's intervals, see if there's some overlap

                const float p1_start = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j]);
                const float p1_end = QVector3D::dotProduct(slot_dir_3d, p1_isecs_3d[j+1]);
                p1l.push_back(p1_start);
                p1l.push_back(p1_end);
        }
        qDebug() << "p0l" << p0l;
        qDebug() << "p1l" << p1l;
    }
    */

    //old method
    /*
    const QVector2D slot_dir = (slot_end - slot_start).normalized() * 1.005;
    const QVector2D slot_ortho = QVector2D(-slot_dir.y(), slot_dir.x());

    QList <QVector3D> isecs;
    QList <bool> isecs_which;

    //order is consistent for i-j intersection this way
    if (index < other_index) {
        GetContourIntersections(other, isecs, isecs_which);
    }
    else {
        other.GetContourIntersections((*this), isecs, isecs_which);
    }

    //we now have an ordered set of contour-plane intersections between the two planes
    bool inside_i = false;
    bool inside_j = false;

    for (int k=0; k<isecs.size()-3; ++k) {

        const float dist1 = (isecs[k+1]-isecs[k]).length();
        const float dist2 = (isecs[k+3]-isecs[k+2]).length();

        bool flip = false;

        //as we move across axis, we go in and out of the planar contours
        if (isecs_which[k]) {
            inside_i = !inside_i;
        }
        else {
            inside_j = !inside_j;
        }

        //we want to make a slot when we are inside one, and we encounter another one
        if (isecs_which[k] != isecs_which[k+1] && ((inside_i && !inside_j) || (inside_j && !inside_i))) {

            if (dist1 > dist2) {

                if (index < other_index) {
                    if (isecs_which[k]) {
                        flip = true;
                        //qDebug() << "S1";
                    }
                    else {
                        flip = false;
                        //qDebug() << "S2";
                    }
                }
                else {
                    if (isecs_which[k]) {
                        flip = false;
                        //qDebug() << "S3";
                    }
                    else {
                        flip = true;
                        //qDebug() << "S4";
                    }
                }


            }
            else {

                if (index < other_index) {
                    if (isecs_which[k]) {
                        flip = false;
                        //qDebug() << "S5";
                    }
                    else {
                        flip = false;
                        //qDebug() << "S6";
                    }
                }
                else {
                    if (isecs_which[k]) {
                        flip = true;
                        //qDebug() << "S7";
                    }
                    else {
                        flip = true;
                        //qDebug() << "S8";
                    }
                }

            }

            QVector2D slot_start1, slot_start2, slot_mid;
            QList <QVector2D> slot_rect1, slot_rect2;

            slot_start1 = flip ? GetPoint2D(isecs[k+2]) : GetPoint2D(isecs[k+1]);
            slot_start2 = flip ? GetPoint2D(isecs[k+1]) : GetPoint2D(isecs[k+2]);
            slot_mid = GetPoint2D((isecs[k+1] + isecs[k+2]) * 0.5);

            float slot_len = (slot_start1 - slot_mid).length();
            slot_start1 = slot_mid + (slot_start1 - slot_mid).normalized() * slot_len * 1.01f;
            slot_start2 = slot_mid + (slot_start2 - slot_mid).normalized() * slot_len * 1.01f;

            slot_rect1.push_back(slot_start1 + slot_ortho * SlabThickness() * 0.5f);
            slot_rect1.push_back(slot_start1 - slot_ortho * SlabThickness() * 0.5f);
            slot_rect1.push_back(slot_mid - slot_ortho * SlabThickness() * 0.5f);
            slot_rect1.push_back(slot_mid + slot_ortho * SlabThickness() * 0.5f);
            slot_rect1.push_back(slot_start1 + slot_ortho * SlabThickness() * 0.5f);

            slot_rect2.push_back(slot_start2 + slot_ortho * SlabThickness() * 0.5f);
            slot_rect2.push_back(slot_start2 - slot_ortho * SlabThickness() * 0.5f);
            slot_rect2.push_back(slot_mid - slot_ortho * SlabThickness() * 0.5f);
            slot_rect2.push_back(slot_mid + slot_ortho * SlabThickness() * 0.5f);
            slot_rect2.push_back(slot_start2 + slot_ortho * SlabThickness() * 0.5f);

            if (LinesIntersectContour(slot_rect1) && !LinesIntersectContour(slot_rect2)) {

                my_slots.push_back(slot_start1);
                my_slots.push_back(slot_mid);
                my_slots_rect.push_back(slot_rect1);

            }
            else if (!LinesIntersectContour(slot_rect1) && LinesIntersectContour(slot_rect2)) {

                my_slots.push_back(slot_start2);
                my_slots.push_back(slot_mid);
                my_slots_rect.push_back(slot_rect2);

            }
            else {

                //qDebug() << "PlanarSection::GetSlots() - Problem! slot_rect1 and slot_rect2 both inadequate for boundary intersection";

                float slot_len = (slot_start1 - slot_mid).length();
                slot_start1 = slot_mid + (slot_start1 - slot_mid).normalized() * slot_len * 1.5f;
                slot_start2 = slot_mid + (slot_start2 - slot_mid).normalized() * slot_len * 1.5f;

                slot_rect1.clear();
                slot_rect2.clear();

                slot_rect1.push_back(slot_start1 + slot_ortho * SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_start1 - slot_ortho * SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_mid - slot_ortho * SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_mid + slot_ortho * SlabThickness() * 0.5f);
                slot_rect1.push_back(slot_start1 + slot_ortho * SlabThickness() * 0.5f);

                slot_rect2.push_back(slot_start2 + slot_ortho * SlabThickness() * 0.5f);
                slot_rect2.push_back(slot_start2 - slot_ortho * SlabThickness() * 0.5f);
                slot_rect2.push_back(slot_mid - slot_ortho * SlabThickness() * 0.5f);
                slot_rect2.push_back(slot_mid + slot_ortho * SlabThickness() * 0.5f);
                slot_rect2.push_back(slot_start2 + slot_ortho * SlabThickness() * 0.5f);

                if (LinesIntersectContour(slot_rect1) && !LinesIntersectContour(slot_rect2)) {

                    my_slots.push_back(slot_start1);
                    my_slots.push_back(slot_mid);
                    my_slots_rect.push_back(slot_rect1);

                }
                else if (!LinesIntersectContour(slot_rect1) && LinesIntersectContour(slot_rect2)) {

                    my_slots.push_back(slot_start2);
                    my_slots.push_back(slot_mid);
                    my_slots_rect.push_back(slot_rect2);

                }
                else {
                    //we need to expand out the rectangle some now
                    qDebug() << "slotlen" << slot_len << slot_start1 << slot_mid << slot_start2 << isecs[k+1] << isecs[k+2] << isecs_which[k] << isecs_which[k+1] << isecs_which[k+2] << isecs_which[k+3] << isecs_which;
                    qDebug() << "PlanarSection::GetSlots() - Problem! slot_rect1 and slot_rect2 both inadequate for boundary intersection";
                }


            }

            //markers.push_back(GetPoint3D(slot_start1));
            //markers_col.push_back(QVector3D(1, 0, 0));
            //markers.push_back(GetPoint3D(slot_mid));
            //markers_col.push_back(QVector3D(1, 1, 0));
            //markers.push_back(GetPoint3D(slot_start2));
            //markers_col.push_back(QVector3D(0, 1, 0));

            k += 2;

        }

    }
    */

}

bool PlanarSection::GetIntersectionSlotLine(const PlanarSection & other, QVector2D & slot_start, QVector2D & slot_end) const
{

    QVector3D lp, ld; //line point, line dir
    if (!GLutils::PlanePlaneIntersection(n, p, other.n, other.p, lp, ld)) {
        return false;
    }

    slot_start = GetPoint2D(lp);
    slot_end = GetPoint2D(lp + ld);

    return true;

}

void PlanarSection::GetContourIntersections(const PlanarSection & other, QList <QVector3D> & isecs, QList <bool> & isecs_which) const
{

    isecs.clear();
    isecs_which.clear();

    QVector3D intersect;

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & samples = bez_curve[c].Samples();
        for (int i=0; i<samples.size(); ++i) {
            const int j = (i+1) % samples.size();
            if (GLutils::LineSegmentPlaneIntersection(other.p, other.n, GetPoint3D(samples[i]), GetPoint3D(samples[j]), intersect)) {
                isecs.push_back(intersect);
                isecs_which.push_back(true);
            }
        }

    }

    for (int c=0; c<other.bez_curve.size(); ++c) {

        const QList <QVector2D> & other_samples = other.bez_curve[c].Samples();
        for (int i=0; i<other_samples.size(); ++i) {
            const int j = (i+1) % other_samples.size();
            if (GLutils::LineSegmentPlaneIntersection(p, n, other.GetPoint3D(other_samples[i]), other.GetPoint3D(other_samples[j]), intersect)) {
                isecs.push_back(intersect);
                isecs_which.push_back(false);
            }
        }

    }

    //now sort them    
    if (isecs.size() >= 2) {

        const QVector3D d = isecs[1] - isecs[0];
        //qDebug() << "PlanarSection::GetContourIntersections() - d length" << d.length();
        GLutils::SortPointsAlongDirection3DExtra(d, isecs, isecs_which);

    }

    //now remove isecs where there is no alternation    
    for (int i=0; i+1<isecs.size(); i+=2) {
    //for (int i=0; i+1<isecs.size(); ++i) {
        if (isecs_which[i] == isecs_which[i+1]) {
            isecs.removeAt(i);
            isecs.removeAt(i);
            isecs_which.removeAt(i);
            isecs_which.removeAt(i);
            i -= 2;
        }
    }

}

bool PlanarSection::LinesIntersectContour(const QList <QVector2D> & lines) const
{

    for (int c=0; c<bez_curve.size(); ++c) {
        const QList <QVector2D> & samples = bez_curve[c].Samples();
        for (int i=0; i<lines.size(); ++i) {
            const int l = (i+1) % lines.size();
            for (int j=0; j<samples.size(); ++j) {
                const int k = (j+1) % samples.size();

                QVector2D intersect;
                if (GLutils::LineLineIntersection(lines[i], lines[l], samples[j], samples[k], intersect)) {
                    return true;
                }

            }

        }
    }

    return false;

}

void PlanarSection::GetBoundingBox2D(QVector2D & min_v, QVector2D & max_v) const
{

    min_v = QVector2D(FLT_MAX, FLT_MAX);
    max_v = QVector2D(-FLT_MAX, -FLT_MAX);

    const QList <QVector2D> & samples = bez_curve[0].Samples();

    for (int i=0; i<samples.size(); ++i) {

        min_v.setX(qMin(min_v.x(), samples[i].x()));
        min_v.setY(qMin(min_v.y(), samples[i].y()));

        max_v.setX(qMax(max_v.x(), samples[i].x()));
        max_v.setY(qMax(max_v.y(), samples[i].y()));

    }

}

void PlanarSection::GetBoundingInterval(const QVector2D & direction, float & min_v, float & max_v) const
{

    min_v = FLT_MAX;
    max_v = -FLT_MAX;

    const QList <QVector2D> & samples = bez_curve[0].Samples();

    for (int i=0; i<samples.size(); ++i) {

        const float proj_length = QVector2D::dotProduct(direction, samples[i]);

        min_v = qMin(min_v, proj_length);
        max_v = qMax(max_v, proj_length);

    }

}

void PlanarSection::ComputeContourLineIntersects(const QVector2D & lp, const QVector2D & ld, QList <QVector2D> & intersects) const
{

    intersects.clear();

    const QList <QVector2D> & samples = bez_curve[0].Samples();
    for (int i=0; i<samples.size(); ++i) {

        const int ind1 = i;
        const int ind2 = (i+1) % samples.size();

        QVector2D intersect;
        if (GLutils::LineRayIntersection(samples[ind1], samples[ind2], lp, ld, intersect)) {
            //note: we do not add duplicates... e.g. the endpoint of a line segment was hit
            //if (!intersects.contains(intersect)) {
                intersects.push_back(intersect);
            //}
        }

    }

    //eliminate any "almost" duplicates
    for (int i=0; i<intersects.size(); ++i) {
        for (int j=i+1; j<intersects.size(); ++j) {
            if ((intersects[i]-intersects[j]).length() < 0.0001f) {
                intersects.removeAt(j);
                --j;
            }
        }
    }

}

void PlanarSection::ComputeIntersectors(const QList <PlanarSection> & sections, const int which_section, QVector <bool> & intersectors)
{

    intersectors = QVector <bool> (sections.size(), false);

    for (int j=0; j<sections.size(); ++j) {

        if (which_section == j) {
            continue;
        }

        intersectors[j] = sections[which_section].IsIntersectingSection(sections[j]);

    }

}

void PlanarSection::ComputeIntersectionGraph(const QList <PlanarSection> & sections, QVector <QVector <bool> > & graph)
{

    //init graph
    graph.clear();
    graph.resize(sections.size());    

    //compute connectivity
    for (int i=0; i<sections.size(); ++i) {
        ComputeIntersectors(sections, i, graph[i]);

    }

}

void PlanarSection::TestConnectedness(QList <PlanarSection> & sections, QList <QList <int> > & all_paths)
{

    if (sections.empty()) {
        return;
    }

    sections[0].SetConnected(true);
    for (int i=1; i<sections.size(); ++i) {

        sections[i].SetConnected(false);
        for (int j=0; j<all_paths.size(); ++j) {

            if (all_paths[j].contains(0) && all_paths[j].contains(i)) {
                sections[i].SetConnected(true);
                break;
            }

        }

        //qDebug() << "section" << i << "connected" << sections[i].Connected();

    }

}

void PlanarSection::SetPartOfCycle(const bool b)
{
    part_of_cycle = b;
}

bool PlanarSection::PartOfCycle() const
{
    return part_of_cycle;
}

bool PlanarSection::CycleAssemblable(const QList <PlanarSection> & sections, const QList <int> & cycle)
{

    QVector <QVector3D> int_dir;
    int_dir.resize(cycle.size());

    //get the intersection lines in the cycle
    for (int i=0; i<cycle.size(); ++i) {
        const QVector3D n1 = sections[cycle[i]].N();
        const QVector3D n2 = sections[cycle[(i+1)%cycle.size()]].N();
        int_dir[i] = QVector3D::crossProduct(n1, n2).normalized();
    }

    //do pairwise comparison of intersection lines, seeing if two are roughly parallel
    for (int i=0; i<cycle.size(); ++i) {
        for (int j=i+1; j<cycle.size(); ++j) {

            const float dot_prod = QVector3D::dotProduct(int_dir[i], int_dir[j]);
            if (fabsf(dot_prod) >= 0.999f) {
                return true;
            }

        }
    }

    //didnt find two about-parallel lines of intersection, cycle cannot be assembled
    return false;

}

void PlanarSection::SetConnected(const bool b)
{
    connected = b;
}

bool PlanarSection::Connected() const
{
    return connected;
}

void PlanarSection::SetQuality(const int i)
{

    quality_samples = i;

    for (int i=0; i<bez_curve.size(); ++i) {
        bez_curve[i].SetSamplesPerSegment(quality_samples);
    }

}

void PlanarSection::SketchClear(const int c)
{
    sketch_pts[c].clear();
}

void PlanarSection::SketchAdd(const int c, const QVector2D & v)
{
    sketch_pts[c].push_back(v);
}

void PlanarSection::SketchSetCurve(const int c, const QList <QVector2D> & v)
{
    sketch_pts[c] = v;
}

void PlanarSection::SketchSetEditing(const bool b)
{
    editing_sketch = b;
}

int PlanarSection::SketchNumPoints() const
{    
    return sketch_pts[0].size();
}

bool PlanarSection::SketchEditing()
{
    return editing_sketch;
}

bool PlanarSection::CurveClockwise(const int c)
{

    //const QList <QVector2D> & pts = bez_curve[c].Points();
    const QList <QVector2D> & pts = bez_curve[c].Samples();

    if (pts.size() < 2) {
        return false;
    }

    float total = 0.0f;
    for (int i=0; i<pts.size()-1; ++i) {
        total += (pts[i+1].x() - pts[i].x()) * (pts[i+1].y()+pts[i].y());
    }
    total += (pts.first().x() - pts.last().x()) * (pts.first().y()+pts.last().y());

    return total >= 0.0f;

}

void PlanarSection::SketchSymmetryTest()
{

    const float magic_ratio = 1.0f;
    // original value was 6.0f
    //const float magic_ratio = 6.0f;

    //get region1/region2 distances
    float dist_r1 = 0.0f;
    float dist_r2 = 0.0f;

    for (int i=0; i<sketch_pts[0].size(); ++i) {

        //wrt local coordinate system, the region sidedness is based on the x-coordinate of each point
        //(a line where x=0 divides the two regions in local 2D parameter space)
        //qDebug() << i << sketch_pts[i].x();
        dist_r1 = qMin(dist_r1, float(sketch_pts[0][i].x()));
        dist_r2 = qMax(dist_r2, float(sketch_pts[0][i].x()));

    }

    //qDebug() << "dist_r1" << dist_r1 << "dist_r2" << dist_r2;

    float sign = 1.0f; //this stores whether points to preserve are on positive side (sign=1) or negative side (sign=-1)

    if ((dist_r1 > 0.0f && dist_r2 > 0.0f) || (fabsf(dist_r2) / fabsf(dist_r1) > magic_ratio)) { //mirror to region1 (negative half)
        sign = 1.0f;
    }
    else if ((dist_r1 < 0.0f && dist_r2 < 0.0f) || (fabsf(dist_r1) / fabsf(dist_r2) > magic_ratio)) { //mirror to region2 (positive half)
        sign = -1.0f;
    }
    else {
        return;
    }

    //first remove points on negative side
    for (int i=0; i<sketch_pts[0].size(); ++i) {
        if (sketch_pts[0][i].x() * sign < 0.0f) {
            sketch_pts[0].removeAt(i);
            --i;
        }
    }

    //add some intermediate points
    QVector2D pt1 = sketch_pts[0].last();
    QVector2D pt2(0, pt1.y());
    const int pts_to_add = 10;
    for (int i=1; i<pts_to_add; ++i) {
        const float f = float(i) / float(pts_to_add);
        QVector2D new_pt = pt1 * (1.0f - f) + pt2 * f;
        sketch_pts[0].push_back(new_pt);
    }

    const int nPts = sketch_pts[0].size();

    //reverse and append it
    for (int i=nPts-1; i>=0; --i) {
        sketch_pts[0].push_back(QVector2D(-sketch_pts[0][i].x(), sketch_pts[0][i].y()));
    }

}

void PlanarSection::MirrorControlPoints()
{

    if (bez_curve.empty()) {
        return;
    }

    if (bez_curve.first().GetNumControlPoints() < 2) {
        return;
    }

    BezierCurve & curve = bez_curve.first();

    QList <QVector2D> points = curve.Points();
    QList <QVector2D> new_points = points;

    new_points.pop_back(); //remove last 2 ctrl points
    new_points.pop_back();

    for (int i=points.size()-3; i>=0; --i) { //copy the middle bit
        new_points.push_back(points[i]);
        new_points.last().setX(-new_points.last().x()); //flip it
    }

    new_points.push_back(points[points.size()-2]); //copy the tangent going into the "start/end" point
    new_points.last().setX(-new_points.last().x());

    new_points.push_back(points[points.size()-2]); //copy the two final points on the original side to cap it off
    new_points.push_back(points[points.size()-1]);

    curve.SetPoints(new_points);
}

void PlanarSection::CreateLocalSymmetry()
{

    if (bez_curve.empty()) {
        return;
    }

    if (bez_curve.first().GetNumControlPoints() < 2) {
        return;
    }

    QList <QVector2D> new_points;

    int sign = 1;
    if (bez_curve[0].Point(0).x() < 0) {
        sign = -1;
    }

    // first half of curve
    new_points.push_back(bez_curve[0].Point(0));
    new_points.push_back(bez_curve[0].Point(1));

    // forward walk for the first section
    int index = 3;
    int smallest_x = 0;

    while (sign*bez_curve[0].Point(index).x() > 0.0f && index<bez_curve[0].GetNumControlPoints()-2) {

        new_points.push_back(bez_curve[0].Point(index-1));
        new_points.push_back(bez_curve[0].Point(index));
        new_points.push_back(bez_curve[0].Point(index+1));

        if (bez_curve[0].Point(index).x() < bez_curve[0].Point(smallest_x).x()) {
            smallest_x = index;
        }

        index+=3;

    }

    // reverse the first section
    for (int i=new_points.size()-1; i>=0; i--)
    {
       new_points.push_back(QVector2D(-new_points[i].x(), new_points[i].y()));
    }

    int lowerSectionLastPointIndex = new_points.size()-1;

    // reverse second last point
    new_points.push_back(QVector2D(-bez_curve[0].Point(bez_curve[0].GetNumControlPoints()-2).x(), bez_curve[0].Point(bez_curve[0].GetNumControlPoints()-2).y()));

    // backward walk to get reversed last part
    index = bez_curve[0].GetNumControlPoints()-4;
    while(sign*bez_curve[0].Point(index).x() > 0.0f && index>2)
    {
        new_points.push_back(bez_curve[0].Point(index+1));
        new_points.push_back(bez_curve[0].Point(index));
        new_points.push_back(bez_curve[0].Point(index-1));

        index-=3;
    }

    // reverse the last section
    for (int i=new_points.size()-1; i>=lowerSectionLastPointIndex; i--)
    {
       new_points.push_back(QVector2D(-new_points[i].x(), new_points[i].y()));
    }

    bez_curve[0].SetPoints(new_points);


}



void PlanarSection::DrawSketch()
{    

    for (int c=0; c<sketch_pts.size(); ++c) {
        glBegin(GL_LINE_STRIP);
        for (int i=0; i<sketch_pts[c].size(); ++i) {
            QVector3D v = p + (t * sketch_pts[c][i].x()) + (b * sketch_pts[c][i].y());
            glVertex3f(v.x(), v.y(), v.z());
        }
        glEnd();
    }

}

bool PlanarSection::AddMouseRayIntersect(const int c, const QVector2D & v)
{   

    //const float min_dist_threshold = 0.025f;
    const float min_dist_threshold = 0.05f;

    QVector3D intersect;
    bool result = MouseRayIntersect(v, intersect);

    if (result) {

        //enforce minimum distance
        const QVector2D new_pt = GetPoint2D(intersect);
        if (sketch_pts[c].empty() || (new_pt - sketch_pts[c].last()).length() > min_dist_threshold) {
            sketch_pts[c].push_back(GetPoint2D(intersect));
        }

    }

    //qDebug() << "PlanarSection::AddMouseRayIntersect - " << intersect;
    return result;

}


bool PlanarSection::AddCtrlPointPenPress(const int c, const QVector2D & v)
{
    QVector3D intersect;
    bool result = MouseRayIntersect(v, intersect);

    QList <QVector2D> bez_points = bez_curve[c].Points();

    QVector2D last_tangent;
    QVector2D intersect2D = GetPoint2D(intersect);

    if (result) {

        // Delete the end added from the previous times
        if(bez_points.size() == 0)
        {
            bez_points.push_back(intersect2D);
            bez_points.push_back(intersect2D);
            last_tangent = intersect2D;
        }
        else
        {

            // removing the points added to the back
            bez_points.pop_back();
            last_tangent = bez_points.back();
            bez_points.pop_back();

            // adding the new points
            bez_points.push_back(intersect2D);
            bez_points.push_back(intersect2D);
            bez_points.push_back(intersect2D);
        }


        // Add the last the last 2 points to connect the curve
        bez_points.push_back(last_tangent);
        //bez_curve[c].AddPoint( bez_curve[c].Point(0) - (bez_curve[c].Point(1) - bez_curve[c].Point(0) ) );
        bez_points.push_back(bez_points.front());

        bez_curve[c].SetPoints(bez_points);

        // Set closed if large enough
        if(bez_curve[c].Points().size() > 3)
            bez_curve[c].SetClosed(true);


        UpdateCurveTrisSlab();

    }

    //qDebug() << "PlanarSection::AddMouseRayIntersect - " << intersect;
    return result;
}



QVector2D PlanarSection::GetPoint2D(const QVector3D & v) const
{

    const QVector3D p_0 = v - p;
    return QVector2D(QVector3D::dotProduct(p_0, t),
                     QVector3D::dotProduct(p_0, b));

}

QVector3D PlanarSection::GetPoint3D(const QVector2D & v) const
{

    return (t * v.x()) + (b * v.y()) + p;

}

QVector3D PlanarSection::GetTangent3D(const QVector2D & v) const
{

    return (t * v.x()) + (b * v.y());

}

QVector3D PlanarSection::Centroid(const QList <PlanarSection> & sections)
{

    QVector3D cent = QVector3D(0, 0, 0);
    float total_area = 0.0f;

    for (int i=0; i<sections.size(); ++i) {

        cent += sections[i].GetCentroid3D() * fabsf(sections[i].SignedArea());
        total_area += fabsf(sections[i].SignedArea());

    }

    cent /= total_area;

    if (total_area <= 0.0f) {
        cent = QVector3D(FLT_MAX, 0, FLT_MAX);
    }

    return cent;

}

QVector2D PlanarSection::Centroid() const
{
    return centroid;
}

float PlanarSection::SignedArea() const
{
    return area;
}

QVector3D PlanarSection::GetCentroid3D() const
{
    return GetPoint3D(centroid);
}

void PlanarSection::ContourVertices3D(QList <QVector3D> & contour_verts) const
{

    contour_verts.clear();

    QList <QVector2D> samples = bez_curve[0].Samples();
    for (int i=0; i<samples.size(); ++i) {
        contour_verts.push_back(GetPoint3D(samples[i]));
    }

}

const QList <QVector2D> & PlanarSection::SliceTriangles() const
{
    return tris;
}

const QList <QVector3D> & PlanarSection::SlabVertices() const
{
    return slab_vert;
}

const QList <QVector3D> & PlanarSection::SlabNormals() const
{
    return slab_norm;
}

float PlanarSection::SlabThickness() const
{
    return slab_thickness;
}

void PlanarSection::SetSlabThickness(const float f)
{
    slab_thickness = f;
}

bool PlanarSection::RayIntersect(const QVector3D & p0, const QVector3D & dir, QVector3D & intersect)
{

    float denom = QVector3D::dotProduct(n, dir);

    if (fabs(denom) < 0.0001f) {
        return false;
    }

    float u = QVector3D::dotProduct(n, (p - p0)) / denom;
    intersect = p0 + (dir * u);

    return true;

}

bool PlanarSection::MouseOutsideDeadzone(const QVector2D & v, const QVector3D & slot_start, const QVector3D & slot_end, const float deadzone_radius)
{

    QVector3D intersect;
    MouseRayIntersect(v, intersect);

    return GLutils::PointLineSegmentDistance(intersect, slot_start, slot_end) >  deadzone_radius;

}

void PlanarSection::SelectMouseRayIntersect(const QVector2D & v, const float cam_width)
{

    QVector3D intersect;
    MouseRayIntersect(v, intersect);

    QVector2D int_2d = GetPoint2D(intersect);

    for (int c=0; c<bez_curve.size(); ++c) {
        bez_curve[c].UnselectPoint();
    }

    selected_weight = -1;

    for (int i=0; i<weights.size(); ++i) {

        const float each_len = (weights[i].p - int_2d).length();
        if (each_len < weights[i].rad * 1.25f) {
            selected_weight = i;
            break;
        }

    }

    int bez_curve_size = bez_curve.size();

    if(radial)
        bez_curve_size = 1;

    bool point_found;
    if (selected_weight < 0) {

        for (int c=0; c<bez_curve_size; ++c) {
            point_found = false;
            bez_curve[c].SelectPoint(int_2d, 0.003f * cam_width);
            if (bez_curve[c].SelectedPoint() >= 0) {

                if(!radial || bez_curve[c].SelectedPoint() <  num_radial_points_per_sector) // ensure that it's from the first sector when radial
                {
                    point_found = true;
                    break;
                }
            }
            if(!point_found)
                bez_curve[c].UnselectPoint();
        }



    }

}

bool PlanarSection::IsCtrlPointSelected()
{
    for (int c=0; c<bez_curve.size(); ++c) {
        if (bez_curve[c].SelectedPoint() >= 0) {
            return true;
        }
    }

    return false;
}

int PlanarSection::SelectedCtrlPoint()
{
    for (int c=0; c<bez_curve.size(); ++c) {
        if (bez_curve[c].SelectedPoint() >= 0) {
            return bez_curve[c].SelectedPoint();
        }
    }

    return -1;
}

void PlanarSection::InsertCtrlPoint()
{

    //qDebug() << "PlanarSection::InsertCtrlPoint() " << bez_curve.SelectedPoint();
    for (int c=0; c<bez_curve.size(); ++c) {
        if (bez_curve[c].SelectedPoint() >= 0) {
            if(radial && bez_curve[c].SelectedPoint() <  num_radial_points_per_sector)
            {
                bez_curve[c].InsertSelectedPoint();
                num_radial_points_per_sector+=3;
                UpdateRadial();
            }
            else
                bez_curve[c].InsertSelectedPoint();


        }
    }

}

void PlanarSection::DeleteSelectedCtrlPoint()
{

    //qDebug() << "PlanarSection::DeleteCtrlPoint()" << bez_curve.SelectedPoint();
    for (int c=0; c<bez_curve.size(); ++c) {
        if (bez_curve[c].SelectedPoint() >= 0) {

            if (c > 0 && bez_curve[c].GetNumControlPoints() <= 7) {
                RemoveCurve(c);
                --c;
            }
            else {
                if(radial && bez_curve[c].SelectedPoint() <  num_radial_points_per_sector)
                {
                    if(num_radial_points_per_sector >= 6)
                    {
                        bez_curve[c].DeleteSelectedPoint();
                        num_radial_points_per_sector-=3;
                        UpdateRadial();
                    }
                }
                else
                    bez_curve[c].DeleteSelectedPoint();
            }

        }
    }

}


void PlanarSection::DeleteCtrlPoint(int bez_curve_index, int ctrl_point_index)
{
    if (bez_curve_index < 0)
    {
        qDebug() << "Error - Curve index less than 0";
        return;
    }

//    if (ctrl_point_index < 0 || ctrl_point_index)
//    {
//        qDebug() << "Error - Control point index less than 0";
//        return;
//    }

    if ( bez_curve[bez_curve_index].GetNumControlPoints() <= 7) {
        if(bez_curve_index > 0)
        {
            RemoveCurve(bez_curve_index);
            --bez_curve_index;
        }
        else
        {
            if(bez_curve[bez_curve_index].GetNumControlPoints() <= 4)
                bez_curve[bez_curve_index] = BezierCurve();
            else
            {
                bool is_closed = bez_curve[bez_curve_index].IsClosed();
                // Set closed allows DeletePoint to work for smaller than 7 points
                bez_curve[bez_curve_index].SetClosed(false);
                bez_curve[bez_curve_index].DeletePoint(ctrl_point_index);
                bez_curve[bez_curve_index].SetClosed(is_closed);
            }
        }
    }
    else {
        bez_curve[bez_curve_index].DeletePoint(ctrl_point_index);
    }


}


bool PlanarSection::IsWeightSelected()
{
    return selected_weight >= 0;
}

void PlanarSection::InsertWeight(const QVector2D & mouse_pos, const float mass)
{

    QVector3D intersect;
    MouseRayIntersect(mouse_pos, intersect);

    QVector2D int_2d = GetPoint2D(intersect);

    selected_weight = weights.size();
    AddWeight(int_2d, mass);

}

void PlanarSection::DeleteWeight()
{
    if (selected_weight >= 0) {
        weights.removeAt(selected_weight);
    }
    selected_weight = -1;
}

void PlanarSection::MoveWeightMouseRayIntersect(const QVector2D & mouse_pos)
{

    QVector3D intersect;
    MouseRayIntersect(mouse_pos, intersect);

    QVector2D int_2d = GetPoint2D(intersect);

    if (selected_weight >= 0) {
        weights[selected_weight].p = int_2d;
    }


}

void PlanarSection::MoveCtrlPointMouseRayIntersect(const QVector2D & mouse_pos, const bool keep_g1, const bool equal_lengths)
{

    QVector3D intersect;
    MouseRayIntersect(mouse_pos, intersect);

    QVector2D int_2d = GetPoint2D(intersect);

    for (int c=0; c<bez_curve.size(); ++c) {
        if (bez_curve[c].SelectedPoint() >= 0) {

            if (!radial || bez_curve[c].SelectedPoint() <  num_radial_points_per_sector) {
                if (radial && bez_curve[c].SelectedPoint() == 0) {
                    bez_curve[c].MoveSelectedPoint(int_2d, false, equal_lengths);
                }
                else {
                    bez_curve[c].MoveSelectedPoint(int_2d, keep_g1, equal_lengths);
                }
            }
        }
    }

    if (radial) {
        UpdateRadial();
    }

}


void PlanarSection::SelectCtrlPoint(const int ctrl_point_index, const int bez_curve_index)
{
    bez_curve[bez_curve_index].SetSelectedPoint(ctrl_point_index);
}

void PlanarSection::UnselectCtrlPoint()
{    
    for (int c=0; c<bez_curve.size(); ++c) {
        bez_curve[c].UnselectPoint();
    }
}

void PlanarSection::SetCtrlPointsAboveXZPlane()
{

    const float just_touch_delta = 0.005f;

    //figure out the line parameters in the 2D parameter space
    QVector3D ld_3d = QVector3D::crossProduct(n, QVector3D(0, 1, 0)).normalized();
    QVector3D lp_3d = QVector3D(0, 0, 0);
    (fabsf(n.x()) > fabsf(n.z())) ? lp_3d.setX(QVector3D::dotProduct(n, p) / n.x()) : lp_3d.setZ(QVector3D::dotProduct(n, p) / n.z());

    //line parameters in 2d space defined by two points, lp1 and lp2
    QVector2D lp1 = GetPoint2D(lp_3d);
    QVector2D lp2 = GetPoint2D(lp_3d + ld_3d);

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & points = bez_curve[c].Points();

        for (int i=0; i<points.size(); ++i) {

            const float dist = GLutils::PointLineSignedDistance(lp1, lp2, points[i]);

            if (dist > 0.0f) { //this means point is below ground plane in the local parameter space

                //set new point location
                QVector2D norm_2d = (lp2 - lp1).normalized();
                norm_2d = QVector2D(-norm_2d.y(), norm_2d.x());

                bez_curve[c].SetPoint(i, points[i] - norm_2d * (dist - just_touch_delta));

                //qDebug() << points[i] << points[i] - norm_2d * (dist - just_touch_delta);

            }

        }

    }

}

void PlanarSection::MirrorPlanarSectionX(PlanarSection & new_section) const
{

    new_section = *this;

    //since we mirror TNB frame, we need to swap T and B to keep handedness of frame the same
    new_section.p = QVector3D(-p.x(), p.y(), p.z());
    new_section.t = QVector3D(-b.x(), b.y(), b.z());
    new_section.n = QVector3D(-n.x(), n.y(), n.z());
    new_section.b = QVector3D(-t.x(), t.y(), t.z());

    for (int c=0; c<bez_curve.size(); ++c) {
        //we take the sketchpts in 3D, mirror them, then reproject onto plane
        for (int i=0; i<new_section.sketch_pts[c].size(); ++i) {
            QVector3D v = GetPoint3D(sketch_pts[c][i]);
            v.setX(-v.x());
            new_section.sketch_pts[c][i] = new_section.GetPoint2D(v);
        }

        const QList <QVector2D> & bez_pts = new_section.bez_curve[c].Points();
        for (int i=0; i<bez_pts.size(); ++i) {
            QVector3D v = GetPoint3D(bez_curve[c].Point(i));
            v.setX(-v.x());
            new_section.bez_curve[c].SetPoint(i, new_section.GetPoint2D(v));
        }
    }

}

void PlanarSection::MirrorPlanarSectionZ(PlanarSection & new_section) const
{

    new_section = *this;

    //since we mirror TNB frame, we need to swap T and B to keep handedness of frame the same
    new_section.p = QVector3D(p.x(), p.y(), -p.z());
    new_section.t = QVector3D(b.x(), b.y(), -b.z());
    new_section.n = QVector3D(n.x(), n.y(), -n.z());
    new_section.b = QVector3D(t.x(), t.y(), -t.z());

    for (int c=0; c<bez_curve.size(); ++c) {

        //we take the sketchpts in 3D, mirror them, then reproject onto plane
        for (int i=0; i<new_section.sketch_pts[c].size(); ++i) {
            QVector3D v = GetPoint3D(sketch_pts[c][i]);
            v.setZ(-v.z());
            new_section.sketch_pts[c][i] = new_section.GetPoint2D(v);
        }

        const QList <QVector2D> & bez_pts = new_section.bez_curve[c].Points();
        for (int i=0; i<bez_pts.size(); ++i) {
            QVector3D v = GetPoint3D(bez_curve[c].Point(i));
            v.setZ(-v.z());
            new_section.bez_curve[c].SetPoint(i, new_section.GetPoint2D(v));
        }

    }

}

void PlanarSection::RotatePlanarSectionY(PlanarSection & new_section) const
{

    new_section = *this;

    //since we mirror TNB frame, we need to swap T and B to keep handedness of frame the same
    new_section.p = QVector3D(p.z(), p.y(), -p.x());
    new_section.t = QVector3D(t.z(), t.y(), -t.x());
    new_section.n = QVector3D(n.z(), n.y(), -n.x());
    new_section.b = QVector3D(b.z(), b.y(), -b.x());

    //we take the sketchpts in 3D, mirror them, then reproject onto plane
    /*
    for (int i=0; i<new_section.sketch_pts.size(); ++i) {
        QVector3D v = GetPoint3D(sketch_pts[i]);
        v.setZ(-v.z());
        new_section.sketch_pts[i] = new_section.GetPoint2D(v);
    }

    const QList <QVector2D> & bez_pts = new_section.bez_curve.Points();
    for (int i=0; i<bez_pts.size(); ++i) {
        QVector3D v = GetPoint3D(bez_curve.Point(i));
        v.setZ(-v.z());
        new_section.bez_curve.SetPoint(i, new_section.GetPoint2D(v));
    }
    */

}

bool PlanarSection::MouseRayIntersect(const QVector2D & v, QVector3D & intersect)
{

    QVector3D pt1, pt2;
    GLutils::UnProjectPoint(v, 0.0, pt1);
    GLutils::UnProjectPoint(v, 0.1, pt2);

    return RayIntersect(pt1, pt2-pt1, intersect);

}

void PlanarSection::DrawTris()
{

    glBegin(GL_TRIANGLES);
    glNormal3f(n.x(), n.y(), n.z());
    for (int i=0; i<tris.size(); i+=3) {
        for (int j=i; j<i+3; ++j) {
            QVector3D v = p + (t * tris[j].x()) + (b * tris[j].y());
            glVertex3f(v.x(), v.y(), v.z());
        }
    }
    glEnd();

}

void PlanarSection::DrawSlab()
{

    /*
    if (update_slab_disp_list) {

        if (slab_disp_list > 0) {
            glDeleteLists(slab_disp_list, 1);
        }

        slab_disp_list = 0;
        update_slab_disp_list = false;

    }

    if (slab_disp_list > 0) {
        glCallList(slab_disp_list);
    }
    else {

        slab_disp_list = glGenLists(1);
        glNewList(slab_disp_list, GL_COMPILE_AND_EXECUTE);
        */

        glBegin(GL_TRIANGLES);
        for (int i=0; i<slab_norm.size(); ++i) {
            for (int j=0; j<3; ++j) {
                glNormal3f(slab_norm[i].x(), slab_norm[i].y(), slab_norm[i].z());
                glVertex3f(slab_vert[i*3+j].x(), slab_vert[i*3+j].y(), slab_vert[i*3+j].z());
            }
        }
        glEnd();

        /*
        glEndList();

    }
    */

}

void PlanarSection::DrawSlabCurves()
{

    /*
    if (update_slabcurves_disp_list) {

        if (slabcurves_disp_list > 0) {
            glDeleteLists(update_slabcurves_disp_list, 1);
        }

        slabcurves_disp_list = 0;
        update_slabcurves_disp_list = false;

    }

    if (slabcurves_disp_list > 0) {
        glCallList(slabcurves_disp_list);
    }
    else {

        slabcurves_disp_list = glGenLists(1);
        glNewList(slabcurves_disp_list, GL_COMPILE_AND_EXECUTE);
        */

    const QVector3D n_offset = n * (slab_thickness / 2.0f);

    glBegin(GL_LINES);
    for (int c=0; c<bez_curve.size(); ++c) {
        const QList <QVector2D> & samples = bez_curve[c].Samples();

        for (int j=0; j<2; ++j) {

            for (int i=0; i<samples.size(); ++i) {

                //const float colour = float(i) / float(samples.size());
                //glColor3f(colour, 1.0f - colour, 0.0f);
                const int i2 = (i+1) % samples.size();
                QVector3D v1 = p + (t * samples[i].x()) + (b * samples[i].y());
                QVector3D v2 = p + (t * samples[i2].x()) + (b * samples[i2].y());

                if (j == 0) {
                    v1 += n_offset;
                    v2 += n_offset;
                }
                else {
                    v1 -= n_offset;
                    v2 -= n_offset;
                }

                glVertex3f(v1.x(), v1.y(), v1.z());
                glVertex3f(v2.x(), v2.y(), v2.z());

            }

        }
    }
    glEnd();

    /*
        glEndList();

    }
    */

}

void PlanarSection::ComputeCentroidAndArea()
{

    area = 0.0f;
    centroid = QVector2D(0, 0);

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & samples = bez_curve[c].Samples();

        //http://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
        //http://paulbourke.net/geometry/polygonmesh/

        //compute signed area
        float each_area = 0.0f;
        QVector2D each_centroid = QVector2D(0, 0);
        for (int i=0; i<samples.size(); ++i) {

            const int j = (i+1) % samples.size();
            const float cross_prod = samples[i].x() * samples[j].y() - samples[j].x() * samples[i].y();

            each_area += cross_prod;

            each_centroid.setX(each_centroid.x() + (samples[i].x() + samples[j].x()) * cross_prod);
            each_centroid.setY(each_centroid.y() + (samples[i].y() + samples[j].y()) * cross_prod);

        }

        each_area *= 0.5f;
        each_centroid /= 6.0f * each_area;

        //when c > 0, the curves represent holes
        if (c > 0) {
            each_area = -fabsf(each_area);
        }

        //total area is a sum of areas (negative area for holes)
        area += each_area;

        //total centroid is a weighted average of the c points
        centroid += each_centroid * each_area;

    }

    centroid /= area;

}

void PlanarSection::ComputeSlab()
{

    slab_vert.clear();
    slab_norm.clear();   

    const QVector3D n_offset = n * (slab_thickness / 2.0f);

    //compute slabs
    //two sides
    for (int i=0; i<tris.size(); ++i) {

        QVector3D each = GetPoint3D(tris[i]) + n_offset;
        slab_vert.push_back(each);

        if (i % 3 == 0) {
            slab_norm.push_back(n);
        }

    }
    for (int i=tris.size()-1; i>=0; --i) {

        QVector3D each = GetPoint3D(tris[i]) - n_offset;
        slab_vert.push_back(each);

        if (i % 3 == 0) {
            slab_norm.push_back(-n);
        }

    }

    //now the ring

    for (int c=0; c<bez_curve.size(); ++c) {

        const QList <QVector2D> & samples = bez_curve[c].Samples();

        bool clockwise = CurveClockwise(c);
        if (c >= 1) {
            clockwise = !clockwise; //for holes, this relationship is inverted
        }

        for (int i=0; i<samples.size(); ++i) {

            QVector3D v1, v2, v3, v4;

            if (i == samples.size()-1) {

                v1 = GetPoint3D(samples[i]) + n_offset;
                v2 = GetPoint3D(samples[0]) + n_offset;
                v3 = GetPoint3D(samples[i]) - n_offset;
                v4 = GetPoint3D(samples[0]) - n_offset;

            }
            else {

                v1 = GetPoint3D(samples[i]) + n_offset;
                v2 = GetPoint3D(samples[i+1]) + n_offset;
                v3 = GetPoint3D(samples[i]) - n_offset;
                v4 = GetPoint3D(samples[i+1]) - n_offset;

            }

            QVector3D cross1 = (v2-v1).normalized();
            QVector3D cross2 = (v3-v1).normalized();
            QVector3D each_n = QVector3D::crossProduct(cross2, cross1).normalized();

            //normal direction on ring based on curve being clockwise
            if (clockwise) {
                slab_vert.push_back(v1);
                slab_vert.push_back(v2);
                slab_vert.push_back(v3);

                slab_vert.push_back(v2);
                slab_vert.push_back(v4);
                slab_vert.push_back(v3);

                slab_norm.push_back(-each_n);
                slab_norm.push_back(-each_n);
            }
            else {
                slab_vert.push_back(v3);
                slab_vert.push_back(v2);
                slab_vert.push_back(v1);

                slab_vert.push_back(v3);
                slab_vert.push_back(v4);
                slab_vert.push_back(v2);

                slab_norm.push_back(each_n);
                slab_norm.push_back(each_n);
            }

        }

    }

}

void PlanarSection::GetXZPlanePoints(QList <QVector3D> & plane_pts) const
{

    plane_pts.clear();

    const QList <QVector2D> & tris = SliceTriangles();
    const QVector3D intersect_line = QVector3D::crossProduct(n, QVector3D(0, 1, 0));

    float min_val = FLT_MAX;
    QVector3D min_pt;
    float max_val = -FLT_MAX;
    QVector3D max_pt;

    //for each plane, compute all the XZ plane intersections
    for (int j=0; j<tris.size(); j+= 3) {

        for (int k=0; k<3; ++k) { //for each line of triangle

            const QVector3D v1 = GetPoint3D(tris[j+k]);
            const QVector3D v2 = GetPoint3D(tris[j+((k+1)%3)]);

            QVector3D i1;

            if (GLutils::LineSegmentPlaneIntersection(QVector3D(0, 0, 0), QVector3D(0, 1, 0), v1, v2, i1)) {

                i1.setY(0.0f); //just to be sure

                const float each_val = QVector3D::dotProduct(intersect_line, i1);

                if (each_val < min_val) {
                    min_val = each_val;
                    min_pt = i1;
                }
                if (each_val > max_val) {
                    max_val = each_val;
                    max_pt = i1;
                }

            }

        }

    }

    //we are only interested in the two extremal points on this line (line formed by intersection of each plane with ground plane)
    if (min_val < FLT_MAX && max_val > -FLT_MAX) {
        plane_pts.push_back(min_pt);
        plane_pts.push_back(max_pt);
    }

}

void PlanarSection::GetXZPlanePoints(const QList <PlanarSection> & sections, QList <QVector3D> & plane_pts)
{

    plane_pts.clear();

    for (int i=0; i<sections.size(); ++i) {

        QList <QVector3D> each_plane_pts;
        sections[i].GetXZPlanePoints(each_plane_pts);

        plane_pts.append(each_plane_pts);

    }

}

void PlanarSection::GetXZConvexHull(const QList <PlanarSection> & sections, QList <QVector3D> & convex_hull)
{

    QList <QVector3D> plane_pts;
    QList <int> hull;

    PlanarSection::GetXZPlanePoints(sections, plane_pts);
    GLutils::ConvexHull_GiftWrapping(plane_pts, hull);

    convex_hull.clear();
    for (int i=0; i<hull.size(); ++i) {
        convex_hull.push_back(plane_pts[hull[i]]);
    }

}

bool PlanarSection::GetXZPointInsideConvexHull(const QList <QVector3D> & hull, const QVector3D & p)
{

    const QVector3D p2 = QVector3D(p.x(), 0, p.z());
    return GLutils::ConvexHull_PointInside(hull, p2);

}

void PlanarSection::SnapSketchToTemplateCut(const QList <QVector3D> & cut_segments, const float snap_dist_3d)
{

    if (sketch_pts[0].empty()) {
        return;
    }

    const int i = sketch_pts[0].size()-1;

    //bring sketch pt to 3D
    const QVector3D pt3d = GetPoint3D(sketch_pts[0][i]);

    QVector3D last_pt3d;
    if (i > 0) {
        last_pt3d = GetPoint3D(sketch_pts[0][i-1]);
    }
    else {
        last_pt3d = pt3d;
    }

    //find the closest line segment from the cut
    float best_dist = FLT_MAX;
    QVector3D best_pt;

    for (int j=0; j<cut_segments.size(); j+=2) {

        QVector3D each_pt;
        QVector3D last_each_pt;
        float each_dist = GLutils::PointLineSegmentDistance(pt3d, cut_segments[j], cut_segments[j+1], each_pt) +
                GLutils::PointLineSegmentDistance(last_pt3d, cut_segments[j], cut_segments[j+1], last_each_pt);
        if (each_dist < best_dist) {
            best_dist = each_dist;
            best_pt = each_pt;
        }

    }

    //snap if within a threshold
    if (best_dist < snap_dist_3d) {
        sketch_pts[0][i] = GetPoint2D(best_pt);
    }

}

void PlanarSection::SaveConnectionsToOBJ(const QList <PlanarSection> & sections, QTextStream & ofs)
{

    QList <int> conn1, conn2;
    QList <QVector3D> conn_start, conn_end;

    for (int i=0; i<sections.size(); ++i) {
        for (int j=0; j<sections.size(); ++j) {

            if (i == j) {
                continue;
            }

            QList <QVector2D> my_slots;
            QList <QList <QVector2D> > my_slots_rect;
            sections[i].GetSlots(sections[j], my_slots, my_slots_rect);

            for (int k=0; k<my_slots.size(); k+=2) {
                conn1.push_back(i);
                conn2.push_back(j);
                conn_start.push_back(sections[i].GetPoint3D(my_slots[k]));
                conn_end.push_back(sections[i].GetPoint3D(my_slots[k+1]));
            }

        }
    }

    ofs << "# PlanarSections " << sections.size() << "\n";
    for (int i=0; i<sections.size(); ++i) {
        ofs << "# " << sections[i].N().x() << " " << sections[i].N().y() << " " << sections[i].N().z() << " "
            << sections[i].P().x() << " " << sections[i].P().y() << " " << sections[i].P().z() << "\n";
    }

    ofs << "# Connections " << conn1.size() << "\n";
    for (int i=0; i<conn1.size(); ++i) {
        ofs << "# " << conn1[i] << " " << conn2[i] << " "
            << conn_start[i].x() << " " << conn_start[i].y() << " " << conn_start[i].z() << " "
            << conn_end[i].x() << " " << conn_end[i].y() << " " << conn_end[i].z() << "\n";
    }

}

void PlanarSection::SplitAlongLine(const QVector2D & split_p, const QVector2D & split_d, QList <PlanarSection> & split_sections)
{

    //test if bounding box crosses the line
    QVector2D min_v, max_v;
    GetBoundingBox2D(min_v, max_v);

    if (!GLutils::LineBoxIntersection(min_v, max_v, split_p, split_d)) {
        return;
    }

    QList <QVector2D> int_pt;
    QList <QVector2D> int_index;

    //algorithm to split (roughly):
    //1. march clockwise along curve (starting anywhere). during marching test to see if split line is crossed (if never crossed, no split)
    //2. if the splitting line is crossed, march along splitting line (in clockwise direction) until another part of the curve is encountered.  eventually we get back to 1st point
    //we also want to continue along that point of intersection along the curve, as this defines another piece on the other side

    QList <BezierCurve> curves;
    bez_curve[0].SplitAlongLine(split_p, split_d, curves);

    for (int i=0; i<curves.size(); ++i) {
        PlanarSection new_section = (*this);
        new_section.bez_curve[0] = curves[i];
        new_section.UpdateCurveTrisSlab();
        split_sections.push_back(new_section);
    }

}

bool PlanarSection::IsIntersectingSection(const PlanarSection & other) const
{

    QList <QVector2D> my_slots;
    QList <QList <QVector2D> > my_slots_rect;
    this->GetSlots(other, my_slots, my_slots_rect);

    return !(my_slots.empty());

    //NOTE: Below code is BUGGY
    //test plane parallelism
    /*
    QVector2D slot_start, slot_end;
    if (GetIntersectionSlotLine(other, slot_start, slot_end)) {

        //planes not parallel
        //test if there is a planar section intersection
        QList <QVector3D> isecs;
        QList <bool> isecs_which;
        GetContourIntersections(other, isecs, isecs_which);

        bool inside_i = false;
        bool inside_j = false;
        for (int k=0; k<isecs.size()-3; ++k) {

            //as we move across axis, we go in and out of the planar contours
            if (isecs_which[k]) {
                inside_i = !inside_i;
            }
            else {
                inside_j = !inside_j;
            }

            //we want to make a slot when we are inside one, and we encounter another one
            if (isecs_which[k] != isecs_which[k+1] && ((inside_i && !inside_j) || (inside_j && !inside_i))) {

                return true;

            }

        }        

    }

    return false;
    */

}

void PlanarSection::WidenSlot(const float factor, QList <QVector2D> & slot)
{

    if (slot.size() != 2) {
        return;
    }

    QVector2D slot_mid = (slot[0] + slot[1]) * 0.5f;
    float slot_len = (slot[0] - slot_mid).length();
    slot[0] = slot_mid + (slot[0] - slot_mid).normalized() * slot_len * factor;
    slot[1] = slot_mid + (slot[1] - slot_mid).normalized() * slot_len * factor;

}

void PlanarSection::GetCurveBetweenSections(const PlanarSection & section1, const PlanarSection & section2, BezierCurve & curve)
{

    //1.  Find the segments and t values on the 2 curves that these sections intersect.
    //  1a) First we need the slot lines for each

    QList <QVector2D> slots1, slots2;
    QList <QList <QVector2D> > slots1_rect, slots2_rect;
    GetSlots(section1, slots1, slots1_rect);
    GetSlots(section2, slots2, slots2_rect);

    if (slots1.empty() || slots2.empty()) {
        return;
    }

    WidenSlot(1.5f, slots1);
    WidenSlot(1.5f, slots2);

    //  1b) Then we test the intersections between Bezier curve and slot lines
    QList <BezierCurvePoint> sec1_curve_ints;
    QList <BezierCurvePoint> sec2_curve_ints;

    bez_curve[0].GetLineSegmentIntersections(slots1.first(), slots1.last(), sec1_curve_ints);
    bez_curve[0].GetLineSegmentIntersections(slots2.first(), slots2.last(), sec2_curve_ints);

    if (sec1_curve_ints.empty() || sec2_curve_ints.empty()) {
        qDebug() << "PlanarSection::GetCurveBetweenSections() - Warning, GetLineSegmentIntersections returning empty";
        return;
    }    

    //2.  Construct a new curve going from section1's intersection to section2's intersection
    //    (note1: we try marching both forward and backward for the subcurve)
    //    (note2: we want the curve in the direction that is shorter)
    BezierCurve sub_curve1;
    BezierCurve sub_curve2;

    bez_curve[0].GetSubCurve(sec1_curve_ints.first(), sec2_curve_ints.first(), sub_curve1); //march forward from 1->2    
    bez_curve[0].GetSubCurve(sec2_curve_ints.first(), sec1_curve_ints.first(), sub_curve2);
    sub_curve2.Reverse(); //want to preserve order of going from 1 -> 2, so reversing this is like marching backward from 1->2

    curve = (sub_curve1.Length() < sub_curve2.Length()) ? sub_curve1 : sub_curve2;

}

void PlanarSection::GetCurveAroundSection(const PlanarSection & section, BezierCurve & curve)
{

    //1.  Find the segments and t values on the 2 curves that these sections intersect.
    //  1a) First we need the slot lines for each
    QVector2D slot_start;
    QVector2D slot_end;
    QList <QVector3D> isecs;
    QList <bool> isecs_which;

    if (!GetIntersectionSlotLine(section, slot_start, slot_end)) {
        qDebug() << "PlanarSection::GetCurveAroundSection() - Warning, plane intersection with section1 and section2 not found";
        return;
    }

    //qDebug() << "PlanarSection::GetCurveAroundSection() - slot_start/end" << slot_start << slot_end;

    GetContourIntersections(section, isecs, isecs_which);

    if (isecs.size() < 2) {
        qDebug() << "PlanarSection::GetCurveAroundSection() - Warning, GetContourIntersections returning less than 2";
        return;
    }

    //qDebug() << "PlanarSection::GetCurveAroundSection - number of isecs" << isecs.size() << isecs_which;

    //  1b) Then we test the intersections between Bezier curve and slot lines
    QList <BezierCurvePoint> sec_curve_ints;

    bez_curve[0].GetLineSegmentIntersections(GetPoint2D(isecs.first()), GetPoint2D(isecs.last()), sec_curve_ints);

    if (sec_curve_ints.empty()) {
        qDebug() << "PlanarSection::GetCurveAroundSection() - Warning, GetLineIntersections returning empty";
        return;
    }

    //qDebug() << "PlanarSection::GetCurveAroundSection - sec_curve_ints" << sec_curve_ints.size() << sec_curve_ints.first().point << sec_curve_ints.last().point
    //         << GetPoint3D(sec_curve_ints.first().point) << GetPoint3D(sec_curve_ints.last().point) << section.p;

    //qDebug() << "intersections" << sec_curve_ints.first().segment << sec_curve_ints.first().t;

    //2.  Construct a new curve going from section1's intersection to section2's intersection
    //    (note1: we try marching both forward and backward for the subcurve)
    //    (note2: we want the curve in the direction that is shorter)
    BezierCurve sub_curve;

    //qDebug() << "Curve" << sec_curve_ints.first().segment << sec_curve_ints.first().segment;
    bez_curve[0].GetSubCurve(sec_curve_ints.first(), sec_curve_ints.first(), curve); //march forward from 1->2
    curve.UpdateSamples();

    //qDebug() << "Curve length" << curve.Length();

}

void PlanarSection::GetSectionsAlongCurve(const PlanarSection & start_section, const BezierCurve & curve, const int num_frames, QList <PlanarSection> & sections)
{

    //1.  get the points and tangents on the bezier curve
    QVector <QVector2D> points;
    QVector <QVector2D> tangents;

    curve.GetPointsTangentsAlongCurve(num_frames, points, tangents);

    //2.  create TNB frames from these samples
    for (int i=0; i<num_frames+1; ++i) {

        PlanarSection new_section;

        const QVector3D each_p = GetPoint3D(points[i]);
        const QVector3D each_t = start_section.T();
        const QVector3D each_n = GetTangent3D(tangents[i]).normalized();
        const QVector3D each_b = QVector3D::crossProduct(each_n, each_t).normalized();

        new_section.SetP(each_p);
        new_section.SetT(each_t);
        new_section.SetN(each_n);
        new_section.SetB(each_b);

        sections.push_back(new_section);

    }

}

void PlanarSection::DrawWeights()
{

    for (int i=0; i<weights.size(); ++i) {

        if (selected_weight == i) {
            glColor3f(0.55f, 0.55f, 0.55f);
        }
        else {
            glColor3f(0.35f, 0.35f, 0.35f);
        }

        const float handle_thick = weights[i].rad / 5.0f;
        const float handle_inner = weights[i].rad * 4.0f / 8.0f;
        const float handle_outer = weights[i].rad * 5.0f / 8.0f;

        QVector3D p1 = GetPoint3D(weights[i].p);
        QVector3D p2 = p1 + QVector3D(0, weights[i].height, 0);

        GLutils::DrawCylinder(p1, p2, weights[i].rad);
        GLutils::DrawRing(p2 - QVector3D(handle_thick, 0, 0), p2 + QVector3D(handle_thick, 0, 0), handle_inner, handle_outer);
    }

}

float PlanarSection::GetPointDistance(const QVector3D & v)
{
    return QVector3D::dotProduct(n, (v - p));
}


void PlanarSection::SetRadial(const bool b)
{
    radial = b;
}

bool PlanarSection::IsRadial()
{
    return radial;
}


void PlanarSection::SetNumRadialSectors(int num_sectors)
{
    num_radial_sectors = num_sectors;
}

int PlanarSection::GetNumRadialSectors()
{
    return num_radial_sectors;
}

void PlanarSection::MakeRadial(int num_sectors, int points_per_sector)
{
    QVector2D min_v;
    QVector2D max_v;
    GetBoundingBox2D(min_v, max_v);

    const QVector2D centre = (min_v + max_v) * 0.5f;
    const float rad = (max_v.x() - min_v.x() + max_v.y() - min_v.y()) * 0.25f;
    CreateCircle(centre, rad, num_sectors);

    radial = true;
    num_radial_sectors = num_sectors;
    num_radial_points_per_sector = points_per_sector;
    radial_centre = centre;
}

void PlanarSection::UpdateRadial()
{
    CopySectorRadially(radial_centre, num_radial_sectors, num_radial_points_per_sector);
}

void PlanarSection::SetLocalSymmetryMode(const bool b)
{
    local_symmetry_mode = b;
}
