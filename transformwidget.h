#ifndef TRANSFORMWIDGET_H
#define TRANSFORMWIDGET_H

#include <QtOpenGL>

#include "glutils.h"

enum TransformWidgetState
{
    NONE,
    TRANS_X,
    TRANS_Y,
    TRANS_Z,
    ROT_X,
    ROT_Y,
    ROT_Z,
    SCALE_X,
    SCALE_Y,
    SCALE_Z,
    SCALE_ALL,
    NUM_STATES
};

class TransformWidget
{

public:

    TransformWidget();

    void DrawGL();
    void DrawSelectionGL();

    void SetVisible(const bool b);
    bool GetVisible() const;

    void SetP(const QVector3D & p);
    QVector3D GetP() const;

    void SetX(const QVector3D & p);
    QVector3D GetX() const;

    void SetY(const QVector3D & p);
    QVector3D GetY() const;

    void SetZ(const QVector3D & p);
    QVector3D GetZ() const;

    void SetState(const TransformWidgetState s);
    TransformWidgetState GetState() const;

    void mousePressEvent(const int dx, const int dy);
    void mouseMoveEvent(const int dx, const int dy);
    void mouseReleaseEvent(const int dx, const int dy);

private:

    bool visible;
    TransformWidgetState state;

    QVector3D p; //origin
    QVector3D x; //x direction
    QVector3D y; //y direction
    QVector3D z; //z direction

};

#endif // TRANSFORMWIDGET_H
