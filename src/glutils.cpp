#include "glutils.h"

int GLutils::point_disp_list = 0;
int GLutils::draw_slices = 16;
GLUquadric *GLutils::quadric = gluNewQuadric();

GLutils::GLutils() {}

void GLutils::DrawPointGL(const float x, const float y, const float rad)
{
    glPushMatrix();

    glTranslatef(x, y, 0);
    glScalef(rad, rad, 1);

    if (point_disp_list > 0) {
        glCallList(point_disp_list);
    } else {
        point_disp_list = glGenLists(1);
        glNewList(point_disp_list, GL_COMPILE_AND_EXECUTE);

        glBegin(GL_TRIANGLE_FAN);

        glVertex2f(0, 0);
        for (int i = 0; i <= 50; ++i) {
            float radval = 3.14159f * 2.0f * float(i) / float(50);

            glVertex2f(sinf(radval), cosf(radval));
        }

        glEnd();

        glEndList();
    }

    glPopMatrix();
}

void GLutils::SetColorByHue(const float val, const float range)
{
    float hue = val / range * 255.0f;
    QColor col = QColor::fromHsv(hue, 255, 255);
    glColor3f(col.redF(), col.greenF(), col.blueF());
}

QVector3D GLutils::GetColorByHue(const float val, const float range)
{
    float hue = val / range * 255.0f;
    QColor col = QColor::fromHsv(hue, 255, 255);
    return QVector3D(col.redF(), col.greenF(), col.blueF());
}

void GLutils::ColorByNearestMajorAxis(const QVector3D &v)
{
    QVector3D axes[6];
    axes[0] = QVector3D(1, 0, 0);
    axes[1] = QVector3D(0, 1, 0);
    axes[2] = QVector3D(0, 0, 1);
    axes[3] = QVector3D(-1, 0, 0);
    axes[4] = QVector3D(0, -1, 0);
    axes[5] = QVector3D(0, 0, -1);

    ColorByIndex(6);

    for (int i = 0; i < 6; ++i) {
        if (v == axes[i]) {
            ColorByIndex(i);
        }
    }
}

QString GLutils::ColorByIndexStr(const int ind)
{
    const int i = ind % 11;

    switch (i) {
        case 0:
            return QString("FF6600");
            break;
        case 1:
            return QString("CC0033");
            break;
        case 2:
            return QString("FF0099");
            break;
        case 3:
            return QString("9900CC");
            break;
        case 4:
            return QString("6600FF");
            break;
        case 5:
            return QString("0033CC");
            break;
        case 6:
            return QString("0099FF");
            break;
        case 7:
            return QString("00CC99");
            break;
        case 8:
            return QString("00FF66");
            break;
        case 9:
            return QString("CCFF00");
            break;
        case 10:
        default:
            return QString("CC9900");
            break;
    }
}

void GLutils::ColorByIndex(const int ind)
{
    int numcol = 7;

    if (ind % numcol == 0) {
        glColor3ub(61, 83, 140);
    } else if (ind % numcol == 1) {
        glColor3ub(249, 117, 114);
    } else if (ind % numcol == 2) {
        glColor3ub(159, 21, 141);
    } else if (ind % numcol == 3) {
        glColor3ub(255, 135, 8);
    } else if (ind % numcol == 4) {
        glColor3ub(19, 198, 127);
    } else if (ind % numcol == 5) {
        glColor3ub(222, 221, 120);
    } else if (ind % numcol == 6) {
        glColor3ub(242, 4, 8);
    }
}

/*
void GLutils::DrawTNBFrame(const vec3 & o, const TNBFrame & tnb, const float
scale)
{

    vec3 t = o + tnb.T() * scale;
    vec3 n = o + tnb.N() * scale;
    vec3 b = o + tnb.B() * scale;

    glBegin(GL_LINES);

    glColor3f(1,0,0);
    glVertex3d(o[0], o[1], o[2]);
    glVertex3d(t[0], t[1], t[2]);

    glColor3f(0,1,0);
    glVertex3d(o[0], o[1], o[2]);
    glVertex3d(n[0], n[1], n[2]);

    glColor3f(0,0,1);
    glVertex3d(o[0], o[1], o[2]);
    glVertex3d(b[0], b[1], b[2]);

    glEnd();

}
*/

void GLutils::DrawArrowFixedLength(const QVector3D &p1, const QVector3D &p2,
                                   const float len)
{
    const float cur_len = (p2 - p1).length();
    DrawArrow(p1, p1 + (p2 - p1) / cur_len * len);
}

void GLutils::DrawArrow(const QVector3D &p1, const QVector3D &p2)
{
    const float line_len = (p2 - p1).length();
    const float base_len = line_len * 2.0f / 3.0f;
    const float tip_len = line_len / 3.0f;
    const float base_thick = line_len / 20.0f;
    const float tip_thick = line_len / 10.0f;

    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    gluDisk(quadric, 0.0, base_thick, draw_slices, 1);
    gluCylinder(quadric, base_thick, base_thick, base_len, draw_slices, 1);
    glTranslatef(0, 0, base_len);
    gluDisk(quadric, base_thick, tip_thick, draw_slices, 1);
    gluCylinder(quadric, tip_thick, 0.0, tip_len, draw_slices, 1);

    glPopMatrix();
}

void GLutils::DrawSphere(const QVector3D &p, const float radius)
{
    glPushMatrix();
    glTranslated(p.x(), p.y(), p.z());
    gluSphere(quadric, radius, 16, 16);
    glPopMatrix();
}

void GLutils::DrawCylinder(const QVector3D &p1, const QVector3D &p2,
                           const float radius)
{
    const float line_len = (p2 - p1).length();
    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    // qDebug() << "DrawCylinder theta_deg" << theta_deg << "axis" << axis << p1
    // << p2;

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    gluDisk(quadric, 0.0, radius, draw_slices, 1);
    gluCylinder(quadric, radius, radius, line_len, draw_slices, 1);
    glTranslatef(0, 0, line_len);
    gluDisk(quadric, 0.0, radius, draw_slices, 1);

    glPopMatrix();
}

void GLutils::DrawRing(const QVector3D &p1, const QVector3D &p2,
                       const float inner_rad, const float outer_rad)
{
    const float line_len = (p2 - p1).length();
    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    gluDisk(quadric, inner_rad, outer_rad, draw_slices, 1);
    gluCylinder(quadric, inner_rad, inner_rad, line_len, draw_slices, 1);
    gluCylinder(quadric, outer_rad, outer_rad, line_len, draw_slices, 1);
    glTranslatef(0, 0, line_len);
    gluDisk(quadric, inner_rad, outer_rad, draw_slices, 1);

    glPopMatrix();
}

void GLutils::DrawCylinderFixedLength(const QVector3D &p1, const QVector3D &p2,
                                      const float radius, const float len)
{
    const float cur_len = (p2 - p1).length();
    DrawCylinder(p1, p1 + (p2 - p1) / cur_len * len, radius);
}

/*
vec2 GLutils::ProjectPoint(const vec3 & p)
{

    GLdouble model_view[16];
    GLdouble projection[16];
    GLint viewport[4];
    GLdouble winx, winy, winz;

    glGetDoublev(GL_MODELVIEW_MATRIX, model_view);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    gluProject(p[0], p[1], p[2], model_view, projection, viewport, &winx, &winy,
&winz);

    return vec2(winx, winy);

}
*/

QVector3D GLutils::ProjectPoint(const QVector3D &p)
{
    GLdouble model_view[16];
    GLdouble projection[16];
    GLint viewport[4];
    GLdouble winx, winy, winz;

    glGetDoublev(GL_MODELVIEW_MATRIX, model_view);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    gluProject(p.x(), p.y(), p.z(), model_view, projection, viewport, &winx,
               &winy, &winz);

    return QVector3D(winx, winy, winz);
}

void GLutils::UnProjectPoint(const QVector2D &v, const float depth_value,
                             QVector3D &unproj_pt)
{
    GLdouble model_view[16];
    GLdouble projection[16];
    GLint viewport[4];
    GLdouble objx, objy, objz;

    glGetDoublev(GL_MODELVIEW_MATRIX, model_view);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    gluUnProject(v.x(), v.y(), depth_value, model_view, projection, viewport,
                 &objx, &objy, &objz);

    unproj_pt = QVector3D(objx, objy, objz);
}

void GLutils::EnableBlending()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void GLutils::DisableBlending() { glDisable(GL_BLEND); }

void GLutils::ReadPixelColor(const int x, const int y, unsigned char &r,
                             unsigned char &g, unsigned char &b)
{
    glReadPixels(x, y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &r);
    glReadPixels(x, y, 1, 1, GL_GREEN, GL_UNSIGNED_BYTE, &g);
    glReadPixels(x, y, 1, 1, GL_BLUE, GL_UNSIGNED_BYTE, &b);
}

void GLutils::ReadPixelColor_FrontBuffer(const int x, const int y,
                                         unsigned char &r, unsigned char &g,
                                         unsigned char &b)
{
    glReadBuffer(GL_FRONT);
    ReadPixelColor(x, y, r, g, b);
}

void GLutils::ReadPixelColor_BackBuffer(const int x, const int y,
                                        unsigned char &r, unsigned char &g,
                                        unsigned char &b)
{
    glReadBuffer(GL_BACK);
    ReadPixelColor(x, y, r, g, b);
}

void GLutils::PixelColorToIndex(const unsigned char r, const unsigned char g,
                                const unsigned char b, int &index)
{
    index = int(b) + (int(g) << 8) + (int(r) << 16);
}

void GLutils::IndexToPixelColor(const int index, unsigned char &r,
                                unsigned char &g, unsigned char &b)
{
    r = (unsigned char)(index >> 16);
    g = (unsigned char)(index >> 8);
    b = (unsigned char)(index);
}

void GLutils::SetPickColor(const int index)
{
    unsigned char r = (unsigned char)(index >> 16);
    unsigned char g = (unsigned char)(index >> 8);
    unsigned char b = (unsigned char)(index);

    glColor3ub(r, g, b);
}

int GLutils::GetWindowWidth()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    return viewport[2];
}

int GLutils::GetWindowHeight()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    return viewport[3];
}

void GLutils::GetViewport(int v[4]) { glGetIntegerv(GL_VIEWPORT, v); }

void GLutils::SetViewport(const int v[4])
{
    glViewport(v[0], v[1], v[2], v[3]);
}

void GLutils::SetViewport(const int x, const int y, const int w, const int h)
{
    glViewport(x, y, w, h);
}

QVector3D GLutils::RotateVector(const QVector3D &v, const QVector3D &axis,
                                const float angle_radians)
{
    if (angle_radians == 0.0f) {
        return v;
    }

    const float sinAngle = sinf(angle_radians);
    const float cosAngle = cosf(angle_radians);
    const float oneSubCos = 1.0f - cosAngle;

    const float x = axis.x();
    const float y = axis.y();
    const float z = axis.z();

    QVector3D R1 = QVector3D(x * x + cosAngle * (1.0f - x * x),
                             x * y * oneSubCos - sinAngle * z,
                             x * z * oneSubCos + sinAngle * y);

    QVector3D R2 = QVector3D(x * y * oneSubCos + sinAngle * z,
                             y * y + cosAngle * (1.0f - y * y),
                             y * z * oneSubCos - sinAngle * x);

    QVector3D R3 = QVector3D(x * z * oneSubCos - sinAngle * y,
                             y * z * oneSubCos + sinAngle * x,
                             z * z + cosAngle * (1.0f - z * z));

    return QVector3D(QVector3D::dotProduct(v, R1), QVector3D::dotProduct(v, R2),
                     QVector3D::dotProduct(v, R3));
}

QVector3D GLutils::MirrorVector(const QVector3D &v, const QVector3D &mirror_dir)
{
    const float dot_prod = QVector3D::dotProduct(v, mirror_dir);
    return v - (mirror_dir * dot_prod * 2.0f);
}

void GLutils::GetBoundingBox(const QVector<QVector3D> &vs, QVector3D &bbox_min,
                             QVector3D &bbox_max)
{
    bbox_min = QVector3D(FLT_MAX, FLT_MAX, FLT_MAX);
    bbox_max = QVector3D(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (int i = 0; i < vs.size(); ++i) {
        bbox_min.setX(qMin(bbox_min.x(), vs[i].x()));
        bbox_min.setY(qMin(bbox_min.y(), vs[i].y()));
        bbox_min.setZ(qMin(bbox_min.z(), vs[i].z()));

        bbox_max.setX(qMax(bbox_max.x(), vs[i].x()));
        bbox_max.setY(qMax(bbox_max.y(), vs[i].y()));
        bbox_max.setZ(qMax(bbox_max.z(), vs[i].z()));
    }
}

float GLutils::AngleBetweenDeg(const QVector3D &v1, const QVector3D &v2)
{
    const float dotprod =
        QVector3D::dotProduct(v1.normalized(), v2.normalized());

    if (dotprod < 1.0f) {
        return acosf(dotprod) * 180.0f / 3.14159f;
    } else {
        return 0.0f;
    }
}

float GLutils::AngleBetweenRad(const QVector3D &v1, const QVector3D &v2)
{
    const float dotprod =
        QVector3D::dotProduct(v1.normalized(), v2.normalized());

    if (dotprod < 1.0f) {
        return acosf(dotprod);
    } else {
        return 0.0f;
    }
}

float GLutils::SignedAngleBetweenRad(const QVector3D &v1, const QVector3D &v2,
                                     const QVector3D &axis)
{
    const float dot = QVector3D::dotProduct(v1, v2);
    // det = x1*y2*zn + x2*yn*z1 + xn*y1*z2 - z1*y2*xn - z2*yn*x1 - zn*y1*x2
    const float det = v1.x() * v2.y() * axis.z() + v2.x() * axis.y() * v1.z() +
                      axis.x() * v1.y() * v2.z() - v1.z() * v2.y() * axis.x() -
                      v2.z() * axis.y() * v1.x() - axis.z() * v1.y() * v2.x();

    return atan2f(det, dot);
}

bool GLutils::LineLineIntersection(const QVector2D &p1, const QVector2D &p2,
                                   const QVector2D &p3, const QVector2D &p4,
                                   QVector2D &intersect)
{
    // http://paulbourke.net/geometry/pointlineplane/

    const float u_denom = (p4.y() - p3.y()) * (p2.x() - p1.x()) -
                          (p4.x() - p3.x()) * (p2.y() - p1.y());

    if (fabs(u_denom) < 0.0001f) {
        return false;
    }

    const float ua_num = (p4.x() - p3.x()) * (p1.y() - p3.y()) -
                         (p4.y() - p3.y()) * (p1.x() - p3.x());
    const float ub_num = (p2.x() - p1.x()) * (p1.y() - p3.y()) -
                         (p2.y() - p1.y()) * (p1.x() - p3.x());

    const float u1 = ua_num / u_denom;
    const float u2 = ub_num / u_denom;

    if (u1 >= 0.0f && u1 <= 1.0f && u2 >= 0.0f && u2 <= 1.0f) {
        intersect = p1 + (p2 - p1) * u1;
        return true;
    } else {
        return false;
    }
}

bool GLutils::LineRayIntersection(const QVector2D &p1, const QVector2D &p2,
                                  const QVector2D &lp, const QVector2D &ld,
                                  QVector2D &intersect)
{
    const float u_denom =
        ld.y() * (p2.x() - p1.x()) - ld.x() * (p2.y() - p1.y());

    if (fabs(u_denom) < 0.0001f) {
        return false;
    }

    const float ua_num =
        ld.x() * (p1.y() - lp.y()) - ld.y() * (p1.x() - lp.x());

    const float u1 = ua_num / u_denom;

    if (u1 >= 0.0f && u1 <= 1.0f) {
        intersect = p1 + (p2 - p1) * u1;
        return true;
    } else {
        return false;
    }
}

bool GLutils::LineBoxIntersection(const QVector2D &p1, const QVector2D &p2,
                                  const QVector2D &lp, const QVector2D &ld)
{
    // p3---- p2
    //|      |
    // p1-----p4

    QVector2D p3 = QVector2D(0, p2.y() - p1.y());
    QVector2D p4 = QVector2D(0, p2.x() - p1.x());

    QVector2D intersect;

    return LineRayIntersection(p1, p4, lp, ld, intersect) ||
           LineRayIntersection(p4, p2, lp, ld, intersect) ||
           LineRayIntersection(p2, p3, lp, ld, intersect) ||
           LineRayIntersection(p3, p1, lp, ld, intersect);
}

bool GLutils::LinePlaneIntersection(const QVector3D &p0, const QVector3D &n,
                                    const QVector3D &l0, const QVector3D &l1,
                                    QVector3D &intersect)
{
    // http://en.wikipedia.org/wiki/Line-plane_intersection
    if (QVector3D::dotProduct(l1 - l0, n) < 0.0001f) {
        return false;
    }

    const float p0n = QVector3D::dotProduct(p0, n);
    const float l0n = QVector3D::dotProduct(l0, n);
    const float l1n = QVector3D::dotProduct(l1, n);

    const float interp = (p0n - l0n) / (l1n - l0n);
    intersect = l0 + (l1 - l0) * interp;

    return true;
}

bool GLutils::LineSegmentPlaneIntersection(const QVector3D &p0,
                                           const QVector3D &n,
                                           const QVector3D &p1,
                                           const QVector3D &p2,
                                           QVector3D &intersect)
{
    // http://en.wikipedia.org/wiki/Line-plane_intersection
    const float p0n = QVector3D::dotProduct(p0, n);
    const float p1n = QVector3D::dotProduct(p1, n);
    const float p2n = QVector3D::dotProduct(p2, n);

    const float interp = (p0n - p1n) / (p2n - p1n);

    intersect = p1 + (p2 - p1) * interp;

    if (qMin(p1n, p2n) <= p0n && qMax(p1n, p2n) >= p0n) {
        return true;
    }

    return false;
}

float GLutils::PointLineSignedDistance(const QVector2D &l1, const QVector2D &l2,
                                       const QVector2D &p)
{
    QVector2D p0 = p - l1;
    QVector2D ldir = (l2 - l1).normalized();

    p0 = p0 - ldir * QVector2D::dotProduct(ldir, p0);

    return QVector2D::dotProduct(p0, QVector2D(-ldir.y(), ldir.x()));
}

// finds distance between point p0 and line segment defined by points p1, p2
float GLutils::PointLineSegmentDistance(const QVector3D &p0,
                                        const QVector3D &p1,
                                        const QVector3D &p2)
{
    // project p0 onto line p1, p2
    const QVector3D p1p0 = p0 - p1;
    const QVector3D p1p2 = p2 - p1;
    const float p1p2len = p1p2.length();

    // get projection's parameter along line
    const float t = QVector3D::dotProduct(p1p0, p1p2.normalized()) / p1p2len;

    if (t >= 1.0f) {
        return (p0 - p2).length();
    } else if (t <= 0.0f) {
        return p1p0.length();
    } else {
        // subtract the projection from p1p0, use that length (its the
        // orthogonal length from p0 to line)
        return (p1p0 - (p1p2 * t)).length();
    }
}

float GLutils::PointLineSegmentDistance(const QVector3D &p0,
                                        const QVector3D &p1,
                                        const QVector3D &p2, QVector3D &best_pt)
{
    // project p0 onto line p1, p2
    const QVector3D p1p0 = p0 - p1;
    const QVector3D p1p2 = p2 - p1;
    const float p1p2len = p1p2.length();

    // get projection's parameter along line
    const float t = QVector3D::dotProduct(p1p0, p1p2.normalized()) / p1p2len;

    if (t >= 1.0f) {
        best_pt = p2;
        return (p0 - p2).length();
    } else if (t <= 0.0f) {
        best_pt = p1;
        return p1p0.length();
    } else {
        // subtract the projection from p1p0, use that length (its the
        // orthogonal length from p0 to line)
        best_pt = p1 + p1p2 * t;
        return (p1p0 - (p1p2 * t)).length();
    }
}

bool GLutils::PlanePlaneIntersection(const QVector3D &n1, const QVector3D &p1,
                                     const QVector3D &n2, const QVector3D &p2,
                                     QVector3D &lp, QVector3D &ld)
{
    // http://paulbourke.net/geometry/pointlineplane/

    ld = QVector3D::crossProduct(n1, n2);

    if (ld.length() < 0.0001f) {
        return false;
    }

    float N11 = QVector3D::dotProduct(n1, n1);
    float N22 = QVector3D::dotProduct(n2, n2);
    float N12 = QVector3D::dotProduct(n1, n2);

    float d1 = QVector3D::dotProduct(n1, p1);
    float d2 = QVector3D::dotProduct(n2, p2);

    float det = N11 * N22 - (N12 * N12);
    float c1 = (d1 * N22 - d2 * N12) / det;
    float c2 = (d2 * N11 - d1 * N12) / det;

    lp = n1 * c1 + n2 * c2;

    return true;
}

void GLutils::ConvexHull_GiftWrapping(const QList<QVector3D> &pts,
                                      QList<int> &hull)
{
    hull.clear();

    if (pts.empty()) {
        return;
    } else if (pts.size() == 1) {
        hull.push_back(0);
        return;
    } else if (pts.size() == 2) {
        hull.push_back(0);
        hull.push_back(1);
        return;
    }
    // http://en.wikipedia.org/wiki/Gift_wrapping_algorithm

    // get left-most point (smallest X)

    float min_x = FLT_MAX;
    int point_on_hull = 0;

    for (int i = 0; i < pts.size(); ++i) {
        if (pts[i].x() < min_x) {
            min_x = pts[i].x();
            point_on_hull = i;
        }
    }

    hull.push_back(point_on_hull);

    // modified to guarantee to halt
    // there should be no more than pts.size() iterations
    for (int i = 0; i < pts.size(); ++i) {
        // we now select point pi such that all other points are "to the right"
        // of the line between p0 and p*
        int endpoint = 0;

        for (int j = 0; j < pts.size(); ++j) {
            if (endpoint ==
                hull.last()) {  // we need two distinct points to form a line
                endpoint = j;
            } else {  // endpoint and point_on_hull form a line
                // if line they form is such that pts[j] is "left" of it
                // then pts[j] is our new endpoint and candidate for next hull
                // point
                const QVector3D v1 = pts[endpoint] - pts[hull.last()];
                const QVector3D v2 = pts[j] - pts[hull.last()];

                if (QVector3D::crossProduct(v1, v2).y() > 0.0f) {
                    endpoint = j;
                }
            }
        }

        if (endpoint == hull.first()) {
            // if we are back where we started, halt
            break;
        }

        // otherwise, add this point to hull, keep going
        hull.push_back(endpoint);
    }
}

bool GLutils::ConvexHull_PointInside(const QList<QVector3D> &hull,
                                     const QVector3D &p)
{
    if (hull.size() == 2) {
        const float cross_prod =
            QVector3D::crossProduct(hull[1] - p, hull[0] - p).y();
        const QVector3D hull_vec = (hull[1] - hull[0]).normalized();
        const float hull_length = (hull[1] - hull[0]).length();
        const float p_dot_prod = QVector3D::dotProduct(p - hull[0], hull_vec);

        if (fabsf(cross_prod) < 0.0001f) {
            // points collinear
            // p falls within interval of hull[0] and hull[1]
            if (p_dot_prod >= 0.0f && p_dot_prod <= hull_length) {
                return true;
            }
        }

        // qDebug() << cross_prod << p_dot_prod << hull_length;
        return false;

    } else if (hull.size() >= 3) {
        for (int i = 0; i < hull.size(); ++i) {
            const int j = (i + 1) % hull.size();

            const QVector3D v1 = hull[j] - hull[i];
            const QVector3D v2 = p - hull[i];

            // ensure that this point is on the right side of this line
            if (QVector3D::crossProduct(v1, v2).y() > 0.0f) {
                return false;
            }
        }

        return true;
    }

    return false;
}

QVector3D GLutils::GetVectorNewBasis(
    const QVector3D &old_x, const QVector3D &old_y, const QVector3D &old_z,
    const QVector3D &old_p, const QVector3D &new_x, const QVector3D &new_y,
    const QVector3D &new_z, const QVector3D &new_p, const QVector3D &p,
    const float scale_x, const float scale_y, const float scale_z)
{
    const QVector3D p2 = p - old_p;

    const float x = QVector3D::dotProduct(old_x, p2) * scale_x;
    const float y = QVector3D::dotProduct(old_y, p2) * scale_y;
    const float z = QVector3D::dotProduct(old_z, p2) * scale_z;

    return new_x * x + new_y * y + new_z * z + new_p;
}

void GLutils::SortPointsAlongDirection2D(const QVector2D &dir,
                                         QList<QVector2D> &isecs)
{
    if (isecs.size() < 2) {
        return;
    }

    // bubble sort
    for (int i = 0; i < isecs.size(); ++i) {
        for (int j = i + 1; j < isecs.size(); ++j) {
            const float fi = QVector2D::dotProduct(dir, isecs[i]);
            const float fj = QVector2D::dotProduct(dir, isecs[j]);

            // i bigger than j, swap these two
            if (fi > fj) {
                // do the swap
                QVector2D temp = isecs[i];
                isecs[i] = isecs[j];
                isecs[j] = temp;
            }
        }
    }
}

void GLutils::SortPointsAlongDirection3D(const QVector3D &dir,
                                         QList<QVector3D> &pts)
{
    if (pts.size() < 2) {
        return;
    }

    // bubble sort
    for (int i = 0; i < pts.size(); ++i) {
        for (int j = i + 1; j < pts.size(); ++j) {
            const float fi = QVector3D::dotProduct(dir, pts[i]);
            const float fj = QVector3D::dotProduct(dir, pts[j]);

            // i bigger than j, swap these two
            if (fi > fj) {
                // do the swap
                QVector3D temp = pts[i];
                pts[i] = pts[j];
                pts[j] = temp;
            }
        }
    }
}

void GLutils::SortPointsAlongDirection3DExtra(const QVector3D dir,
                                              QList<QVector3D> &isecs,
                                              QList<bool> &isecs_which)
{
    if (isecs.size() < 2) {
        return;
    }

    // bubble sort
    for (int i = 0; i < isecs.size(); ++i) {
        for (int j = i + 1; j < isecs.size(); ++j) {
            const float fi = QVector3D::dotProduct(dir, isecs[i]);
            const float fj = QVector3D::dotProduct(dir, isecs[j]);

            // i bigger than j, swap these two
            if (fi > fj) {
                // do the swap
                QVector3D temp = isecs[i];
                isecs[i] = isecs[j];
                isecs[j] = temp;

                bool tempb = isecs_which[i];
                isecs_which[i] = isecs_which[j];
                isecs_which[j] = tempb;
            }
        }
    }
}

int GLutils::GetClosestPoint(const QList<QVector3D> &pts, const QVector3D &p)
{
    if (pts.empty()) {
        qDebug() << "GLutils::GetClosestPoint - Warning, no points";
        return -1;
    }

    float min_dist = FLT_MAX;
    int min_index = 0;

    for (int i = 0; i < pts.size(); ++i) {
        float each_dist = (pts[i] - p).lengthSquared();

        if (each_dist < min_dist) {
            min_dist = each_dist;
            min_index = i;
        }
    }

    return min_index;
}

void GLutils::GetClosestPairOfPoints(const QList<QVector3D> &pts1,
                                     const QList<QVector3D> &pts2, int &ind1,
                                     int &ind2)
{
    if (pts1.empty() || pts2.empty()) {
        qDebug() << "GLutils::GetClosestPoint - Warning, no points";
        ind1 = -1;
        ind2 = -1;
        return;
    }

    float min_dist = FLT_MAX;

    for (int i = 0; i < pts1.size(); ++i) {
        for (int j = 0; j < pts2.size(); ++j) {
            float each_dist = (pts1[i] - pts2[j]).lengthSquared();

            if (each_dist < min_dist) {
                min_dist = each_dist;
                ind1 = i;
                ind2 = j;
            }
        }
    }
}

QVector3D GLutils::GetOrthoVec(const QVector3D &v)
{
    QVector3D v2 = v.normalized();

    // preference to Y-up
    QVector3D cross_vec;

    if (fabsf(QVector3D::dotProduct(v2, QVector3D(0, 1, 0))) > 0.9f) {
        cross_vec = QVector3D::crossProduct(QVector3D(0, 0, 1), v2);
    } else {
        cross_vec = QVector3D::crossProduct(QVector3D(0, 1, 0), v2);
    }

    cross_vec.normalize();

    return cross_vec;
}

void GLutils::GetModelviewMatrixFromVec(const QVector3D &v, const QVector3D &p,
                                        float modelview[16])
{
    const QVector3D z = v.normalized();
    const QVector3D y = GetOrthoVec(z);
    const QVector3D x = QVector3D::crossProduct(y, z);

    modelview[0] = x.x();
    modelview[1] = x.y();
    modelview[2] = x.z();
    modelview[3] = 0.0f;
    modelview[4] = y.x();
    modelview[5] = y.y();
    modelview[6] = y.z();
    modelview[7] = 0.0f;
    modelview[8] = z.x();
    modelview[9] = z.y();
    modelview[10] = z.z();
    modelview[11] = 0.0f;
    modelview[12] = p.x();
    modelview[13] = p.y();
    modelview[14] = p.z();
    modelview[15] = 1.0f;
}

void GLutils::GetSortedIntersectionPoints(const QList<QVector2D> &pts,
                                          const QVector2D ray_dir,
                                          const QVector2D ray_p,
                                          QList<QVector2D> &sorted_isecs)
{
    sorted_isecs.clear();

    // 1.  get the intersections
    for (int i = 0; i < pts.size(); ++i) {
        QVector2D p1 = pts[i];
        QVector2D p2 = pts[(i + 1) % pts.size()];

        QVector2D isec;

        if (LineRayIntersection(p1, p2, ray_p, ray_dir, isec)) {
            sorted_isecs.push_back(isec);
        }
    }

    // 2.  sort them
    if (sorted_isecs.size() >= 2) {
        SortPointsAlongDirection2D(sorted_isecs.last() - sorted_isecs.first(),
                                   sorted_isecs);
    }
}

void GLutils::DrawDisc(const QVector3D &p1, const QVector3D &p2,
                       const float inner_rad, const float outer_rad)
{
    const float line_len = (p2 - p1).length();
    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    gluDisk(quadric, inner_rad, outer_rad, draw_slices, 1);

    glPopMatrix();
}

void GLutils::DrawSemiRing(const QVector3D &p1, const QVector3D &p2,
                           float start_angle, float arc_angle,
                           const float inner_rad, const float outer_rad,
                           int segNumber)
{
    const float line_len = (p2 - p1).length();
    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    DrawArc(0, 0, 0, inner_rad, outer_rad, start_angle, arc_angle, segNumber);

    glPopMatrix();
}

void GLutils::DrawSemiRingLineStrip(const QVector3D &p1, const QVector3D &p2,
                                    float start_angle, float arc_angle,
                                    float radius, float thickness,
                                    int segNumber)
{
    const float line_len = (p2 - p1).length();
    const float theta_deg =
        GLutils::AngleBetweenDeg(QVector3D(0, 0, 1), (p2 - p1) / line_len);
    QVector3D axis =
        QVector3D::crossProduct(QVector3D(0, 0, 1), (p2 - p1) / line_len);

    glPushMatrix();

    glTranslatef(p1.x(), p1.y(), p1.z());
    if (theta_deg > 0.0f && theta_deg < 180.0f) {
        glRotatef(theta_deg, axis.x(), axis.y(), axis.z());
    } else {
        glRotatef(theta_deg, 1.0f, 0.0f, 0.0f);
    }

    DrawArcLineStrip(0, 0, 0, radius, start_angle, arc_angle, segNumber,
                     thickness);

    glPopMatrix();
}

// An efficient algorith for drawing an arc avoiding trig functions
// This is from the website: http://slabode.exofire.net/circle_draw.shtml
// The author has released it to the public

void GLutils::DrawArc(float cx, float cy, float cz, float innerR, float outerR,
                      float start_angle, float arc_angle, int num_segments)
{
    float theta =
        arc_angle /
        float(num_segments -
              1);  // theta is now calculated from the arc angle instead, the -
                   // 1 bit comes from the fact that the arc is open

    float tangetial_factor = tanf(theta);

    float radial_factor = cosf(theta);

    float x = cosf(start_angle);  // we now start at the start angle
    float y = sinf(start_angle);

    glBegin(GL_TRIANGLE_STRIP);  // since the arc is not a closed curve, this is
                                 // a strip now
    for (int ii = 0; ii < num_segments; ii++) {
        glVertex3f(innerR * x + cx, innerR * y + cy, cz);

        glVertex3f(outerR * x + cx, outerR * y + cy, cz);

        float tx = -y;
        float ty = x;

        x += tx * tangetial_factor;
        y += ty * tangetial_factor;

        x *= radial_factor;
        y *= radial_factor;
    }
    glEnd();
}

void GLutils::DrawArcLineStrip(float cx, float cy, float cz, float radius,
                               float start_angle, float arc_angle,
                               int num_segments, float thickness)
{
    glPushAttrib(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    float theta =
        arc_angle /
        float(num_segments -
              1);  // theta is now calculated from the arc angle instead, the -
                   // 1 bit comes from the fact that the arc is open

    float tangetial_factor = tanf(theta);

    float radial_factor = cosf(theta);

    float x = cosf(start_angle);  // we now start at the start angle
    float y = sinf(start_angle);

    glLineWidth(thickness);
    glBegin(GL_LINE_STRIP);  // since the arc is not a closed curve, this is a
                             // strip now
    for (int ii = 0; ii < num_segments; ii++) {
        glVertex3f(radius * x + cx, radius * y + cy, cz);

        float tx = -y;
        float ty = x;

        x += tx * tangetial_factor;
        y += ty * tangetial_factor;

        x *= radial_factor;
        y *= radial_factor;
    }
    glEnd();
    glPopAttrib();
    glPopAttrib();
}

void GLutils::DrawArcLineStrip(QVector3D c, QVector3D dir, QVector3D axis,
                               float arc_angle, int num_segments,
                               float thickness, bool dashed)
{
    glPushAttrib(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    float theta =
        arc_angle /
        float(num_segments -
              1);  // theta is now calculated from the arc angle instead, the -
                   // 1 bit comes from the fact that the arc is open

    glPushMatrix();
    glTranslatef(c.x(), c.y(), c.z());
    QVector3D rotVec = RotateVector(dir, axis, -0.5 * arc_angle);

    glPushAttrib(GL_LINE_BIT);
    glLineWidth(thickness);
    glBegin(
        GL_LINES);  // since the arc is not a closed curve, this is a strip now
    for (int ii = 0; ii < num_segments; ii++) {
        glVertex3f(rotVec.x(), rotVec.y(), rotVec.z());
        rotVec = RotateVector(rotVec, axis, theta);
        if (!dashed) glVertex3f(rotVec.x(), rotVec.y(), rotVec.z());
    }
    glEnd();
    glPopAttrib();
    glPopMatrix();
    glPopAttrib();
}

// Note: segements is assumed to be a list of pairs of vertices defining
// segments, an even number of verts is expected
void GLutils::CurvesFromLineSegments(QList<QVector2D> &segments,
                                     QList<QList<QVector2D> > &curves)
{
    // qDebug() << "GLutils::CurvesFromLineSegments - segments size" <<
    // segments.size();
    curves.clear();

    const float max_dist = 0.0005f;

    while (!segments.empty()) {
        QList<QVector2D> curve;
        curve.push_back(segments[0]);
        curve.push_back(segments[1]);

        segments.pop_front();
        segments.pop_front();

        // find the closest vertex to the end of curve
        bool find_next = true;

        do {
            int closest_index1 = GetClosestPoint(segments, curve.last());
            int closest_index2 = GetClosestPoint(segments, curve.first());

            float dist1 = (segments[closest_index1] - curve.last()).length();
            float dist2 = (segments[closest_index2] - curve.first()).length();

            if (dist1 <= dist2 && dist1 < max_dist) {
                find_next = true;

                // add the OTHER point to the current curve we are growing
                if (closest_index1 % 2 == 0) {
                    curve.push_back(segments[closest_index1 + 1]);

                    // remove the segment
                    segments.removeAt(closest_index1);
                    segments.removeAt(closest_index1);
                } else {
                    curve.push_back(segments[closest_index1 - 1]);

                    // remove the segment
                    segments.removeAt(closest_index1 - 1);
                    segments.removeAt(closest_index1 - 1);
                }
            } else if (dist2 <= dist1 && dist2 < max_dist) {
                find_next = true;

                // add the OTHER point to the current curve we are growing
                if (closest_index2 % 2 == 0) {
                    curve.push_front(segments[closest_index2 + 1]);

                    // remove the segment
                    segments.removeAt(closest_index2);
                    segments.removeAt(closest_index2);
                } else {
                    curve.push_front(segments[closest_index2 - 1]);

                    // remove the segment
                    segments.removeAt(closest_index2 - 1);
                    segments.removeAt(closest_index2 - 1);
                }
            } else {
                find_next = false;
            }

            // might have prematurely closed the loop
            if ((curve.first() - curve.last()).length() < max_dist) {
                find_next = false;
            }

        } while (find_next);

        curves.push_back(curve);
        // NOTE: curve closed test presently disabled!
        // push the curve back if it's closed (the endpoints are about the same)
        // if ((curve.first() - curve.last()).length() < max_dist) {
        //
        // }
    }

    // qDebug() << "GLutils::CurvesFromLineSegments - curves created " <<
    // curves.size(); if (curves.size() > 0) {
    //     qDebug() << curves[0];
    // }
}

int GLutils::GetClosestPoint(const QList<QVector2D> &pts, const QVector2D &p)
{
    float best_dist = FLT_MAX;
    int best_index = -1;

    for (int i = 0; i < pts.size(); ++i) {
        const float each_dist = (pts[i] - p).lengthSquared();
        if (each_dist < best_dist) {
            best_dist = each_dist;
            best_index = i;
        }
    }

    return best_index;
}
