#include "pivotcamera.h"

PivotCamera::PivotCamera()
    : neardist(0.1),
      fardist(200.0),
      cam_width(2.0f),
      lookat(0, 0, 0),
      interp_active(false)
{
}

QVector3D PivotCamera::Eye() const { return eye; }

void PivotCamera::SetEye(const QVector3D & e)
{
    eye = e;
    SetUpRightAuto();
}

void PivotCamera::SetRightVector(const QVector3D & v) { right_vec = v; }

void PivotCamera::SetUpRightAuto()
{
    if (fabs(QVector3D::dotProduct((lookat - eye).normalized(),
                                   QVector3D(0, 1, 0))) > 0.99) {
        right_vec = QVector3D::crossProduct((lookat - eye).normalized(),
                                            QVector3D(1, 0, 0))
                        .normalized();
    } else {
        right_vec = QVector3D::crossProduct((lookat - eye).normalized(),
                                            QVector3D(0, 1, 0))
                        .normalized();
    }

    up = QVector3D::crossProduct(GetRightVector(), (lookat - eye).normalized())
             .normalized();
}

QVector3D PivotCamera::LookAt() const { return lookat; }

void PivotCamera::SetLookAt(const QVector3D & p)
{
    lookat = p;
    SetUpRightAuto();
}

QVector3D PivotCamera::ViewDir() const { return (lookat - Eye()).normalized(); }

void PivotCamera::SetWinWidth(int i) { win_width = i; }

void PivotCamera::SetWinHeight(int i) { win_height = i; }

void PivotCamera::SetCamWidth(float f) { cam_width = f; }

float PivotCamera::CamWidth() const { return cam_width; }

void PivotCamera::DrawGL() const
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(75.0, float(win_width) / float(win_height), neardist,
                   fardist);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const QVector3D up = Up();
    gluLookAt(eye.x(), eye.y(), eye.z(), lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());
}

void PivotCamera::DrawGL_Ortho() const
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, win_width, 0.0f, win_height, -1.0f, 1.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void PivotCamera::DrawGL_3DOrtho() const
{
    float aspectr = float(win_width) / float(win_height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-cam_width / 2.0 * aspectr, cam_width / 2.0 * aspectr,
            -cam_width / 2.0, cam_width / 2.0, 0.0, fardist);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const QVector3D up = Up();
    gluLookAt(eye.x(), eye.y(), eye.z(), lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());
}

float PivotCamera::InterpTime() const
{
    float interp = float(interp_time.elapsed()) / interp_duration;

    interp = qMin(1.0f, interp);
    interp = qMax(0.0f, interp);

    return interp;
}

bool PivotCamera::InterpActive() const { return interp_active; }

void PivotCamera::SetInterpActive(const bool b) { interp_active = b; }

// parameters:
//   p - the pivot point
//   v - a central point in the viewport, acts as an offset for the pivot
//   new_viewdir - the outward looking direction of the view at the end of the
//   rotation duration - time for the animation (in milliseconds)
void PivotCamera::StartPivot(const QVector3D & p, const QVector3D & v,
                             const QVector3D & new_viewdir,
                             const float duration)
{
    const QVector3D new_viewdir_norm = new_viewdir.normalized();

    pivot_p = p;
    pivot_v = v - p;
    pivot_e = eye - p;

    const QVector3D v_to_e_dir = (pivot_v - pivot_e).normalized();
    pivot_rot_axis =
        QVector3D::crossProduct(v_to_e_dir, new_viewdir_norm).normalized();
    pivot_rot_theta = GLutils::AngleBetweenRad(v_to_e_dir, new_viewdir_norm);

    pivot_up = up;
    pivot_right = right_vec;
    pivot_offset = QVector3D(0, 0, 0);

    interp_time.start();
    interp_duration = duration;
    interp_active = true;

    // do awesome stuff
    // 1.  find intersection line between 2 planes
    QVector3D lp, ld;
    GLutils::PlanePlaneIntersection(QVector3D(0, 1, 0), QVector3D(0, 0, 0),
                                    new_viewdir, QVector3D(0, 0, 0), lp, ld);

    QVector3D best_upvec;
    best_upvec = GLutils::RotateVector(ld, new_viewdir, 3.14159f / 2.0f);
    if (best_upvec.y() < 0.0f) {
        best_upvec = -best_upvec;
    }

    pivot_final_up = best_upvec;
    pivot_final_right =
        QVector3D::crossProduct(new_viewdir, best_upvec).normalized();
}

void PivotCamera::Update(const QVector2D & mouse_pos)
{
    // strategy:
    // project the pivot point p
    // compare with the mouse

    // unproject both out, at p's depth
    // whatever the difference is, make that the camera translation offset

    // so we rotate some vectors some amount, based on time... if time is
    // greater than duration we stop
    if (interp_active) {
        // 1.  zero out pivot offset
        // pivot_offset = QVector3D(0, 0, 0);

        // 2.  update eye/look/up
        float interp_time = InterpTime();

        // rotate pivot bases about the pivot point
        QVector3D pivoted_v = GLutils::RotateVector(
            pivot_v, pivot_rot_axis, pivot_rot_theta * interp_time);
        QVector3D pivoted_e = GLutils::RotateVector(
            pivot_e, pivot_rot_axis, pivot_rot_theta * interp_time);

        lookat = pivot_p + pivoted_v;
        eye = pivot_p + pivoted_e;

        // slerp the up and right vecs
        up = Slerp(pivot_up, pivot_final_up, interp_time);
        right_vec = Slerp(pivot_right, pivot_final_right, interp_time);

        if (InterpTime() >= 1.0f) {
            interp_active = false;
        }

        // 3.  test offset
        pivot_offset = GetPivotOffset(mouse_pos);

        // 4.  update eye and look
        eye += pivot_offset;
        lookat += pivot_offset;
    }
}

void PivotCamera::Orbit(const QVector2D & delta)
{
    const float factor = 0.005f;

    // rotate existing eye and up by delta
    /*
    eye = lookat + GLutils::RotateVector(eye - lookat, QVector3D(0, 1, 0),
    delta.x() * factor); const QVector3D eye_last = eye; eye = lookat +
    GLutils::RotateVector(eye - lookat, GetRightVector(), -delta.y() * factor);
    if (Up().y() < 0.0f) { //this rotation about up vector made up vector upside
    down, we don't allow it eye = eye_last;
    }
    */

    // rotate about y-axis (x movement)
    eye = lookat + GLutils::RotateVector(eye - lookat, QVector3D(0, 1, 0),
                                         delta.x() * factor);
    right_vec = GLutils::RotateVector(right_vec, QVector3D(0, 1, 0),
                                      delta.x() * factor);
    up = GLutils::RotateVector(up, QVector3D(0, 1, 0), delta.x() * factor);

    // rotate about rightvec-axis (y movement)
    eye = lookat +
          GLutils::RotateVector(eye - lookat, right_vec, -delta.y() * factor);
    up = GLutils::RotateVector(up, right_vec, -delta.y() * factor);
}

void PivotCamera::Dolly(const QVector2D & delta)
{
    const float dolly_factor = 0.00125f * cam_width;
    const QVector3D translate = GetRightVector() * delta.x() * dolly_factor +
                                Up() * delta.y() * dolly_factor;

    eye += translate;
    lookat += translate;
}

void PivotCamera::Zoom(const QVector2D & delta)
{
    const float zoom_factor = 0.005f;

    cam_width *= (1.0f - (delta.y() * zoom_factor));
    cam_width = qMax(cam_width, 1.0f);
    cam_width = qMin(cam_width, 80.0f);
}

QVector3D PivotCamera::Up() const { return up; }

QVector3D PivotCamera::GetRightVector() const { return right_vec; }

QVector3D PivotCamera::GetPivotOffset(const QVector2D & mouse_pos) const
{
    // finally, we need to correct for how we will be off
    this->DrawGL_3DOrtho();

    QVector3D pivot_p_project = GLutils::ProjectPoint(pivot_p);

    QVector3D mouse_pos_unproject;
    GLutils::UnProjectPoint(mouse_pos, pivot_p_project.z(),
                            mouse_pos_unproject);

    QVector3D offset = pivot_p - mouse_pos_unproject;

    // qDebug() << "mouse" << mouse_pos << "proj_p" << pivot_p_project <<
    // "pivot_offset" << pivot_offset; qDebug() << "offset_length" <<
    // offset.length();

    return offset;
}

QVector3D PivotCamera::Slerp(const QVector3D & v1, const QVector3D & v2,
                             const float t)
{
    const float angle = GLutils::AngleBetweenRad(v1, v2);
    const QVector3D axis = QVector3D::crossProduct(v1, v2).normalized();
    return GLutils::RotateVector(v1, axis, angle * t);
}
