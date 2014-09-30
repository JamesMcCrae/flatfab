#include "transformwidget.h"

TransformWidget::TransformWidget() :
    visible(true),
    state(NONE),
    p(0,0,0),
    x(1,0,0),
    y(0,1,0),
    z(0,0,1)
{
}

void TransformWidget::DrawGL(const float cam_width)
{
    //qDebug() << "TransformWidget::DrawGL()";
    if (!visible) {
        return;
    }

    const float s = 0.2f * cam_width;

    //draw various arrows and things
    //translate x
    glColor3f(1,0,0);
    GLutils::DrawArrowFixedLength(p, p+x, 1.0f * s);

    //translate y
    glColor3f(0,1,0);
    GLutils::DrawArrowFixedLength(p, p+y, 1.0f * s);

    //translate z
    glColor3f(0,0,1);
    GLutils::DrawArrowFixedLength(p, p+z, 1.0f * s);

    //rotate x
//    glColor3f(1,0,0);
//    GLutils::DrawRing(p - x*0.05f, p + x*0.05f, 0.45f, 0.55f);
    glColor3f(0,0,1);
    GLutils::DrawArcLineStrip(p, z * 1.5f * s, x, .5*M_PI, 30, 7.0f,true);

    //rotate y
//    glColor3f(0,1,0);
//    GLutils::DrawRing(p - y*0.05f, p + y*0.05f, 0.45f, 0.55f);
    glColor3f(1,0,0);
    GLutils::DrawArcLineStrip(p, x * 1.5f * s, y, .5*M_PI, 30, 7.0f,true);

    //rotate z
//    glColor3f(0,0,1);
//    GLutils::DrawRing(p - z*0.05f, p + z*0.05f, 0.45f, 0.55f);
    glColor3f(0,1,0);
    GLutils::DrawArcLineStrip(p, y * 1.5f * s, z, .5*M_PI, 30, 7.0f,true);

    //scale x
    glColor3f(1,0,0);
    GLutils::DrawSphere(p + x* 1.25f * s, 0.125f * s);

    //scale z
    glColor3f(0,0,1);
    GLutils::DrawSphere(p + z* 1.25f * s, 0.125f  * s);

}

void TransformWidget::DrawSelectionGL(const float cam_width)
{
    if (!visible) {
        return;
    }

    const float s = 0.2f * cam_width;

//    //qDebug() << "TransformWidget::DrawSelectionGL()";
    GLutils::SetPickColor(TRANS_X);
    GLutils::DrawArrowFixedLength(p, p+x, 1.0f * s);

    GLutils::SetPickColor(TRANS_Y);
    GLutils::DrawArrowFixedLength(p, p+y, 1.0f * s);

    GLutils::SetPickColor(TRANS_Z);
    GLutils::DrawArrowFixedLength(p, p+z, 1.0f * s);

    GLutils::SetPickColor(ROT_X);
//    GLutils::DrawRing(p - x*0.05f, p + x*0.05f, 0.45f, 0.55f);
    GLutils::DrawArcLineStrip(p, z * 1.5f * s, x, .5*M_PI, 30, 7.0f);

    GLutils::SetPickColor(ROT_Y);
//    GLutils::DrawRing(p - y*0.05f, p + y*0.05f, 0.45f, 0.55f);
    GLutils::DrawArcLineStrip(p, x * 1.5f * s, y, .5*M_PI, 30, 7.0f);

    GLutils::SetPickColor(ROT_Z);
//    GLutils::DrawRing(p - z*0.05f, p + z*0.05f, 0.45f, 0.55f);
    GLutils::DrawArcLineStrip(p, y * 1.5f * s, z, .5*M_PI, 30, 7.0f);

    GLutils::SetPickColor(SCALE_X);
    GLutils::DrawSphere(p + x * 1.25f * s, 0.125f  * s);

    GLutils::SetPickColor(SCALE_Z);
    GLutils::DrawSphere(p + z * 1.25f * s, 0.125f  * s);
}

void TransformWidget::mousePressEvent(const int dx, const int dy)
{

}

void TransformWidget::mouseMoveEvent(const int dx, const int dy)
{

}

void TransformWidget::mouseReleaseEvent(const int dx, const int dy)
{

}

void TransformWidget::SetVisible(const bool b)
{
    visible = b;
}

bool TransformWidget::GetVisible() const
{
    return visible;
}

void TransformWidget::SetP(const QVector3D & v)
{
    p = v;
}

QVector3D TransformWidget::GetP() const
{
    return p;
}

void TransformWidget::SetX(const QVector3D & v)
{
    x = v;
}

QVector3D TransformWidget::GetX() const
{
    return x;
}

void TransformWidget::SetY(const QVector3D & v)
{
    y = v;
}

QVector3D TransformWidget::GetY() const
{
    return y;
}

void TransformWidget::SetZ(const QVector3D & v)
{
    z = v;
}

QVector3D TransformWidget::GetZ() const
{
    return z;
}

void TransformWidget::SetState(const TransformWidgetState s)
{
    state = s;
}

TransformWidgetState TransformWidget::GetState() const
{
    return state;
}
