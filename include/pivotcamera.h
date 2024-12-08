#ifndef PIVOTCAMERA_H
#define PIVOTCAMERA_H

#include <QtOpenGL>

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "glutils.h"

class PivotCamera
{
public:
    PivotCamera();

    void DrawGL() const;
    void DrawGL_Ortho() const;
    void DrawGL_3DOrtho() const;

    QVector3D Up() const;
    QVector3D Eye() const;
    QVector3D LookAt() const;
    QVector3D GetRightVector() const;
    QVector3D ViewDir() const;

    void SetEye(const QVector3D & e);
    void SetLookAt(const QVector3D & p);

    void SetWinWidth(int i);
    void SetWinHeight(int i);
    void SetCamWidth(float f);
    float CamWidth() const;

    float InterpTime() const;
    bool InterpActive() const;
    void SetInterpActive(const bool b);

    void StartPivot(const QVector3D & p, const QVector3D & v,
                    const QVector3D & new_viewdir, const float duration);
    void Update(const QVector2D & mouse_pos);
    void Orbit(const QVector2D & delta);
    void Dolly(const QVector2D & delta);
    void Zoom(const QVector2D & delta);

private:
    void SetRightVector(const QVector3D & v);
    void SetUpRightAuto();
    QVector3D Slerp(const QVector3D & v1, const QVector3D & v2, const float t);

    // void UpdateUpVector();
    QVector3D GetPivotOffset(const QVector2D & mouse_pos) const;

    double neardist;
    double fardist;

    int win_width;
    int win_height;
    float cam_width;

    QVector3D eye;
    QVector3D right_vec;
    QVector3D up;
    QVector3D lookat;

    bool interp_active;
    QElapsedTimer interp_time;
    float interp_duration;

    // what is needed for interpolation
    QVector3D pivot_p;
    QVector3D pivot_v;
    QVector3D pivot_e;

    QVector3D pivot_up;
    QVector3D pivot_right;
    QVector3D pivot_final_up;
    QVector3D pivot_final_right;

    QVector3D pivot_rot_axis;
    float pivot_rot_theta;
    QVector3D pivot_offset;
};

#endif  // PIVOTCAMERA_H
