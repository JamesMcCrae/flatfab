//
//  rigidbody.cpp
//  inter_plane
//
//  Created by Nobuyuki Umetani on 12/15/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#include <map>

#include "physics.h"

CPhysics::CPhysics()
{
    nitr = 30;
    damping_ratio = 0.01;
    ////
    n = CVector3D(0, 1,
                  0);  // normal direction of floor (should be an unit vector)
    gravity = CVector3D(0, -10, 0);  // gravity
    ////
    cont_stiff = 1.0e+4;   // contact_stiffness (insensitive)
    trans_stiff = 1.0e+4;  // joint_translation_stiffness (insensitive)
    rot_stiff = 1.0e+9;    // joint_rotation_stiffness (insensitive)

    scale_force = 0.05;   // 0.5;
    scale_torque = 0.05;  // 0.5;
    irb_draw = -1;

    is_draw_deformed = false;        // false
    is_draw_skeleton = false;        // false;
    is_draw_force = false;           // false;
    is_draw_section = true;          // true;
    is_draw_grid = false;            // true;
    is_draw_section_moment = false;  // false;

    is_slit = true;
    aDir2D.clear();
    aDir2D.push_back(1);
    aDir2D.push_back(0);
    aDir2D.push_back(0);
    aDir2D.push_back(1);
    aDir2D.push_back(0.70710678118);
    aDir2D.push_back(+0.70710678118);
    aDir2D.push_back(0.70710678118);
    aDir2D.push_back(-0.70710678118);
    nHeight = 20;
    // max_stress = 50000;
    max_stress =
        60000000;  // this is the value for acrylic, 60MPa (60 * 1,000,000)
}

// main class interface methods to set up problem

void CPhysics::AddPlate(const double p[3], const double t[3], const double n[3],
                        const double b[3], const double thickness,
                        const double rho, const std::vector<double>& pts)
{
    CPlate plt;

    plt.p = CVector3D(p[0], p[1], p[2]);
    plt.t = CVector3D(t[0], t[1], t[2]);
    plt.n = CVector3D(n[0], n[1], n[2]);
    plt.b = CVector3D(b[0], b[1], b[2]);
    plt.thickness = thickness;
    plt.rho = rho;
    plt.aXY = pts;

    aPlate.push_back(plt);
}

void CPhysics::AddRigidBody(const double centre_of_mass[3], const double mass,
                            const std::vector<double>& contact_points)
{
    CRigidBody rb;
    rb.cg = CVector3D(centre_of_mass[0], centre_of_mass[1], centre_of_mass[2]);
    rb.m = mass;
    const int ncp = (int)contact_points.size() / 3;
    for (int icp = 0; icp < ncp; icp++) {
        CVector3D pc(contact_points[icp * 3 + 0], contact_points[icp * 3 + 1],
                     contact_points[icp * 3 + 2]);
        rb.aCP.push_back(pc);
    }
    aRigidBody.push_back(rb);
}

void CPhysics::AddJoint(const double position[3], const double jp0[4],
                        const double jp1[4], const int body_index1,
                        const int body_index2)
{
    CJoint j;
    j.p = CVector3D(position[0], position[1], position[2]);
    j.irb0 = body_index1;
    j.irb1 = body_index2;

    for (int i = 0; i < 4; ++i) {
        j.jp0[i] = jp0[i];
        j.jp1[i] = jp1[i];
    }

    aJoint.push_back(j);
}

void CPhysics::AddWeightToRigidBody(const double mass_pos[3],
                                    const double mass_magnitude[3])
{
    if (aRigidBody.empty()) {
        return;
    }

    CRigidBody& c = aRigidBody[aRigidBody.size() - 1];

    const int exf_size = c.aExForce.size();
    c.aExForce.resize(exf_size + 1);

    std::pair<CVector3D, CVector3D> newWeight;
    newWeight.first = CVector3D(mass_pos[0], mass_pos[1], mass_pos[2]);
    newWeight.second =
        CVector3D(mass_magnitude[0], mass_magnitude[1], mass_magnitude[2]);
    c.aExForce[exf_size] = newWeight;
}

void CPhysics::ClearProblem()
{
    aRigidBody.clear();
    aJoint.clear();
    aPlate.clear();
}

void CPhysics::ReadFile(const char* fname)
{
    const double rho = 1;
    FILE* fp = fopen(fname, "r");
    char buff[256];
    int nr;
    fscanf(fp, "%s %d", buff, &nr);
    std::cout << "nrigid: " << nr << std::endl;
    aRigidBody.clear();
    aPlate.clear();
    aJoint.clear();
    aRigidBody.resize(nr);
    aPlate.resize(nr);
    for (int ir = 0; ir < nr; ir++) {
        CRigidBody& rb = aRigidBody[ir];
        CPlate& plt = aPlate[ir];
        int tmp1;
        double d1, d2, d3;
        fscanf(fp, "%s %d", buff, &tmp1);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // PlaneP
        plt.p = CVector3D(d1, d2, d3);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // PlaneT
        plt.t = CVector3D(d1, d2, d3);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // PlaneN
        plt.n = CVector3D(d1, d2, d3);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // PlaneB
        plt.b = CVector3D(d1, d2, d3);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // CenterofMass3D
        rb.cg = CVector3D(d1, d2, d3);
        double area, thickness;
        fscanf(fp, "%s %lf", buff, &area);       // PlaneArea
        fscanf(fp, "%s %lf", buff, &thickness);  // PlaneThickness
        plt.thickness = thickness;
        plt.rho = rho;
        rb.m = area * thickness * rho;
        int ncnt;
        fscanf(fp, "%s %d", buff, &ncnt);  // PlaneContactPoints
        rb.aCP.resize(ncnt);
        for (int icnt = 0; icnt < ncnt; icnt++) {
            fscanf(fp, "%lf %lf %lf", &d1, &d2, &d3);
            rb.aCP[icnt] = CVector3D(d1, d2, d3);
        }
        int np;
        fscanf(fp, "%s %d", buff, &np);  // PlaneBoundary2D
        plt.aXY.resize(np * 2);
        std::cout << "npoint section2d: " << np << std::endl;
        for (int ip = 0; ip < np; ip++) {
            fscanf(fp, "%lf %lf", &d1, &d2);
            plt.aXY[ip * 2 + 0] = d1;
            plt.aXY[ip * 2 + 1] = d2;
        }
        rb.R.SetIdentity();
        rb.u = CVector3D(0, 0, 0);
        rb.aCForce.resize(rb.aCP.size());
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    int nj;
    fscanf(fp, "%s %d", buff, &nj);  // Joints
    aJoint.resize(nj);
    for (int ij = 0; ij < nj; ij++) {
        CJoint& j = aJoint[ij];
        int tmp1, tmp2;
        double d1, d2, d3, d4;
        fscanf(fp, "%s %d %d", buff, &tmp1, &tmp2);  // JointPlaneIndexes
        j.irb0 = tmp1;
        j.irb1 = tmp2;
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // JointCentre3D
        j.p = CVector3D(d1, d2, d3);
        fscanf(fp, "%s %lf %lf %lf", buff, &d1, &d2, &d3);  // JointAxis3D
        fscanf(fp, "%s %lf %lf %lf %lf", buff, &d1, &d2, &d3,
               &d4);  // JointPlane1SlitEndpoints2D
        j.jp0[0] = d1;
        j.jp0[1] = d2;
        j.jp0[2] = d3;
        j.jp0[3] = d4;
        std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << std::endl;
        fscanf(fp, "%s %lf %lf %lf %lf", buff, &d1, &d2, &d3,
               &d4);  // JointPlane1SlitEndpoints2D
        j.jp1[0] = d1;
        j.jp1[1] = d2;
        j.jp1[2] = d3;
        j.jp1[3] = d4;
        std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    int nexf;
    fscanf(fp, "%s %d", buff, &nexf);  // ExternalForces
    for (int iexf = 0; iexf < nexf; iexf++) {
        int irb;
        double x0, y0, z0;
        double x1, y1, z1;
        fscanf(fp, "%d  %lf%lf%lf  %lf%lf%lf", &irb, &x0, &y0, &z0, &x1, &y1,
               &z1);
        assert(irb < int(aRigidBody.size()));
        int irb_exf = (int)aRigidBody[irb].aExForce.size();
        aRigidBody[irb].aExForce.resize(irb_exf + 1);
        aRigidBody[irb].aExForce[irb_exf].first = CVector3D(x0, y0, z0);
        aRigidBody[irb].aExForce[irb_exf].second = CVector3D(x1, y1, z1);
    }

    fclose(fp);
}

void CPhysics::Solve_InsidePlane()
{
    if (aPlate.size() != aRigidBody.size()) return;
    for (int irb = 0; irb < int(aRigidBody.size()); irb++) {
        {
            CPlate& plt = aPlate[irb];
            if (is_slit) {
                GetPlateLocalForce(aPlate[irb].alForce, irb, aRigidBody, aJoint,
                                   plt.p, plt.t, plt.b, plt.n, plt.aXY_slit);
            } else {
                GetPlateLocalForce(aPlate[irb].alForce, irb, aRigidBody, aJoint,
                                   plt.p, plt.t, plt.b, plt.n, plt.aXY);
            }
        }
        aPlate[irb].aPS.clear();
        if (is_slit) {
            AddWeakSection_Inside(aPlate[irb], aPlate[irb].aXY_slit,
                                  aRigidBody[irb], gravity, aDir2D, nHeight,
                                  max_stress);
            // Sampling Around Slit
            AddWeakSection_Slit(aPlate[irb], aPlate[irb].aXY_slit,
                                aRigidBody[irb], gravity, aDir2D, nHeight,
                                max_stress);
        } else {
            AddWeakSection_Inside(aPlate[irb], aPlate[irb].aXY, aRigidBody[irb],
                                  gravity, aDir2D, nHeight, max_stress);
        }
    }
}

void CPhysics::EdEd_Potential(double& energy, CVector3D& dEdu, CVector3D& dEdw,
                              const CVector3D&,  // cg
                              const CVector3D& u, const double mass,
                              const CVector3D& g  // floor normal
)
{
    energy = mass * (u * g);
    dEdu = mass * g;
    dEdw = CVector3D(0, 0, 0);
}

void CPhysics::EdEddE_Exforce(double& energy, CVector3D& dEdu, CVector3D& dEdw,
                              CMatrix3& ddEddu, CMatrix3& ddEddw,
                              CMatrix3& ddEdudw,
                              const CVector3D& cg,   // the center of gravity
                              const CVector3D& pex,  // external force position
                              const CVector3D& fex,  // external force
                              const CVector3D& u,    // displacement
                              const CMatrix3& R      // rigid rotation
)
{
    CVector3D Rv = R * (pex - cg);
    CVector3D qex = Rv + cg + u;  // current external force position
    energy = +qex * fex;
    dEdu = +fex;
    dEdw = +Rv ^ fex;
    ddEddu = CMatrix3(0.0);
    ddEddw = +CMatrix3::Spin(fex) * CMatrix3::Spin(Rv);
    ddEdudw = CMatrix3(0.0);
}

void CPhysics::EdEddE_Contact(double& energy, CVector3D& dEdu, CVector3D& dEdw,
                              CMatrix3& ddEddu, CMatrix3& ddEddw,
                              CMatrix3& ddEdudw,
                              const CVector3D& cg,  // the center of gravity
                              const CVector3D& cp,  // contact position
                              const CVector3D& u,   // displacement
                              const CMatrix3& R,    // rigid rotation
                              const double cont_stiff,
                              const CVector3D& n  // floor normal
)
{
    CVector3D Rv = R * (cp - cg);
    CVector3D cq = Rv + cg + u;
    energy = 0.5 * (cq * n) * (cq * n) * cont_stiff;
    dEdu = ((cq * n) * cont_stiff) * n;
    dEdw = ((cq * n) * cont_stiff) * (Rv ^ n);
    ddEddu = cont_stiff * CMatrix3::OuterProduct(n, n);
    ddEddw = cont_stiff * CMatrix3::OuterProduct(Rv ^ n, Rv ^ n) +
             ((cq * n) * cont_stiff) * CMatrix3::Spin(n) * CMatrix3::Spin(Rv);
    ddEdudw = cont_stiff * CMatrix3::OuterProduct(Rv ^ n, n);
}

void CPhysics::EdEddE_ContactFriction(
    double& energy, CVector3D& dEdu, CVector3D& dEdw, CMatrix3& ddEddu,
    CMatrix3& ddEddw, CMatrix3& ddEdudw,
    const CVector3D& cg,  // the center of gravity
    const CVector3D& cp,  // contact position
    const CVector3D& u,   // displacement
    const CMatrix3& R,    // rigid rotation
    const double cont_stiff)
{
    CVector3D Rv = R * (cp - cg);
    CVector3D cq = Rv + cg + u;
    energy = 0.5 * (cq - cp) * (cq - cp) * cont_stiff;
    dEdu = cont_stiff * (cq - cp);
    dEdw = cont_stiff * (Rv ^ (cq - cp));
    ddEddu = cont_stiff * CMatrix3::Identity();
    //    ddEddw  = -cont_stiff*CMatrix3::OuterProduct(Rv,Rv) +
    //    cont_stiff*CMatrix3::Spin(Rv)*CMatrix3::Spin(cq-cp);
    ddEddw = -cont_stiff * CMatrix3::Spin(Rv) * CMatrix3::Spin(Rv) +
             cont_stiff * CMatrix3::Spin(cq - cp) * CMatrix3::Spin(Rv);
    ddEdudw = cont_stiff * CMatrix3::Spin(Rv);
}

void CPhysics::EdEddE_Joint(double& energy, CVector3D& dEdu0, CVector3D& dEdw0,
                            CVector3D& dEdu1, CVector3D& dEdw1,
                            ////
                            CMatrix3& ddEdu0du0, CMatrix3& ddEdu0dw0,
                            CMatrix3& ddEdu0du1, CMatrix3& ddEdu0dw1,
                            CMatrix3& ddEdw0dw0, CMatrix3& ddEdw0du1,
                            CMatrix3& ddEdw0dw1, CMatrix3& ddEdu1du1,
                            CMatrix3& ddEdu1dw1, CMatrix3& ddEdw1dw1,
                            ////
                            const double trans_stiff, const double rot_stiff,
                            const CVector3D& pj,
                            ////
                            const CVector3D& cg0, const CVector3D& u0,
                            const CMatrix3& R0, const CVector3D& cg1,
                            const CVector3D& u1, const CMatrix3& R1)
{
    CVector3D Rv0 = R0 * (pj - cg0);
    CVector3D qj0 =
        Rv0 + cg0 + u0;  // after deformation joint pos relative to rigid body 0

    CVector3D Rv1 = R1 * (pj - cg1);
    CVector3D qj1 =
        Rv1 + cg1 + u1;  // after deformation joint pos relative to rigid body 1

    energy = 0.5 * trans_stiff * (qj0 - qj1).DLength();

    dEdu0 = trans_stiff * (qj0 - qj1);
    dEdw0 = trans_stiff * (Rv0 ^ (qj0 - qj1));

    dEdu1 = trans_stiff * (qj1 - qj0);
    dEdw1 = trans_stiff * (Rv1 ^ (qj1 - qj0));

    ddEdu0du0 = trans_stiff * CMatrix3::Identity();
    ddEdu0du1 = -trans_stiff * CMatrix3::Identity();
    ddEdu0dw0 = +trans_stiff * CMatrix3::Spin(Rv0);
    ddEdu0dw1 = -trans_stiff * CMatrix3::Spin(Rv1);

    ddEdw0dw0 = -trans_stiff * CMatrix3::Spin(Rv0) * CMatrix3::Spin(Rv0) +
                trans_stiff * CMatrix3::Spin(qj0 - qj1) * CMatrix3::Spin(Rv0);
    ddEdw0du1 = +trans_stiff * CMatrix3::Spin(Rv0);
    ddEdw0dw1 = +trans_stiff * CMatrix3::Spin(Rv1) * CMatrix3::Spin(Rv0);

    ddEdu1du1 = trans_stiff * CMatrix3::Identity();
    ddEdu1dw1 = +trans_stiff * CMatrix3::Spin(Rv1);

    ddEdw1dw1 = -trans_stiff * CMatrix3::Spin(Rv1) * CMatrix3::Spin(Rv1) +
                trans_stiff * CMatrix3::Spin(qj1 - qj0) * CMatrix3::Spin(Rv1);

    CVector3D av(0, 0, 0);
    CMatrix3 davdw0, davdw1;
    {
        for (unsigned int i = 0; i < 3; i++) {
            CVector3D r0(R0.mat[0 * 3 + i], R0.mat[1 * 3 + i],
                         R0.mat[2 * 3 + i]);
            CVector3D r1(R1.mat[0 * 3 + i], R1.mat[1 * 3 + i],
                         R1.mat[2 * 3 + i]);
            av += (r0 ^ r1);
            davdw0 += CMatrix3(r1) * CMatrix3(r0);
            davdw1 -= CMatrix3(r0) * CMatrix3(r1);
        }
        av *= 0.5;
        davdw0 *= 0.5;
        davdw1 *= 0.5;
    }
    energy += 0.5 * rot_stiff * av.DLength();

    CMatrix3 m0, m1, m2;
    for (unsigned int i = 0; i < 3; i++) {
        CVector3D r0(R0.mat[0 * 3 + i], R0.mat[1 * 3 + i], R0.mat[2 * 3 + i]);
        CVector3D r1(R1.mat[0 * 3 + i], R1.mat[1 * 3 + i], R1.mat[2 * 3 + i]);
        dEdw0 += 0.5 * rot_stiff * r0 ^ (r1 ^ av);
        dEdw1 -= 0.5 * rot_stiff * r1 ^ (r0 ^ av);
        ddEdw0dw0 += 0.5 * rot_stiff *
                     (CMatrix3::Spin(r1 ^ av) * CMatrix3::Spin(r0) +
                      CMatrix3::Spin(r0) * CMatrix3::Spin(r1) * davdw0);
        ddEdw1dw1 -= 0.5 * rot_stiff *
                     (CMatrix3::Spin(r0 ^ av) * CMatrix3::Spin(r1) +
                      CMatrix3::Spin(r1) * CMatrix3::Spin(r0) * davdw1);
        ddEdw0dw1 -=
            0.5 * rot_stiff *
            (CMatrix3::Spin(r1) * CMatrix3::Spin(av) * CMatrix3::Spin(r0) +
             CMatrix3::Spin(r1) * CMatrix3::Spin(r0) * davdw0);
    }
}

CVector3D CPhysics::rand_vec(double s)
{
    CVector3D v;
    v.x = s * rand() / (RAND_MAX + 1.0);
    v.y = s * rand() / (RAND_MAX + 1.0);
    v.z = s * rand() / (RAND_MAX + 1.0);
    return v;
}

CMatrix3 CPhysics::rand_rot()
{
    double s = 3.5;
    CVector3D v;
    v.x = s * rand() / (RAND_MAX + 1.0);
    v.y = s * rand() / (RAND_MAX + 1.0);
    v.z = s * rand() / (RAND_MAX + 1.0);
    CMatrix3 R;
    R.SetRotMatrix_Cartesian(v);
    return R;
}

void CPhysics::CheckDiff_Contact()
{
    double energy = 0.0;
    CVector3D dEdu, dEdw;
    CMatrix3 ddEddu, ddEddw, ddEdudw;
    double epsilon = 1.0e-5;
    const double cont_stiff = 1.0e+3;
    CVector3D cg = rand_vec(1.0);
    CVector3D cp = rand_vec(1.0);
    CVector3D u = rand_vec(1.0);
    CMatrix3 R = rand_rot();
    CVector3D n(0, 0, 1);
    EdEddE_Contact(energy, dEdu, dEdw, ddEddu, ddEddw, ddEdudw, cg, cp, u, R,
                   cont_stiff, n);
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D u_ = u;
        u_[idim] += epsilon;
        EdEddE_Contact(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_, ddEdudw_, cg,
                       cp, u_, R, cont_stiff, n);
        std::cout << "dEdu " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdu[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddu  " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEddu.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEdudw.mat[jdim * 3 + idim] << std::endl;
        }
    }
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D w(0, 0, 0);
        w[idim] = epsilon;
        CMatrix3 dR;
        dR.SetRotMatrix_Cartesian(w);
        CMatrix3 R_ = dR * R;
        EdEddE_Contact(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_, ddEdudw_, cg,
                       cp, u, R_, cont_stiff, n);
        std::cout << "dEdw " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdw[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEdudw.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddw  " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEddw.mat[jdim * 3 + idim] << std::endl;
        }
    }
}

void CPhysics::CheckDiff_ContactFriction()
{
    double energy = 0.0;
    CVector3D dEdu, dEdw;
    CMatrix3 ddEddu, ddEddw, ddEdudw;
    double epsilon = 1.0e-5;
    const double cont_stiff = 1.0e+3;
    CVector3D cg = rand_vec(1.0);
    CVector3D cp = rand_vec(1.0);
    CVector3D u = rand_vec(1.0);
    CMatrix3 R = rand_rot();
    EdEddE_ContactFriction(energy, dEdu, dEdw, ddEddu, ddEddw, ddEdudw, cg, cp,
                           u, R, cont_stiff);
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D u_ = u;
        u_[idim] += epsilon;
        EdEddE_ContactFriction(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_,
                               ddEdudw_, cg, cp, u_, R, cont_stiff);
        std::cout << "dEdu " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdu[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddu  " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEddu.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEdudw.mat[jdim * 3 + idim] << std::endl;
        }
    }
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D w(0, 0, 0);
        w[idim] = epsilon;
        CMatrix3 dR;
        dR.SetRotMatrix_Cartesian(w);
        CMatrix3 R_ = dR * R;
        EdEddE_ContactFriction(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_,
                               ddEdudw_, cg, cp, u, R_, cont_stiff);
        std::cout << "dEdw " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdw[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEdudw.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddw  " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEddw.mat[jdim * 3 + idim] << std::endl;
        }
    }
}

void CPhysics::CheckDiff_Exforce()
{
    double energy = 0.0;
    CVector3D dEdu, dEdw;
    CMatrix3 ddEddu, ddEddw, ddEdudw;
    double epsilon = 1.0e-5;
    CVector3D cg = rand_vec(1.0);
    CVector3D u = rand_vec(1.0);
    CVector3D fex = rand_vec(1.0);
    CVector3D pex = rand_vec(1.0);
    CMatrix3 R = rand_rot();
    CVector3D n(0, 0, 1);
    EdEddE_Exforce(energy, dEdu, dEdw, ddEddu, ddEddw, ddEdudw, cg, pex, fex, u,
                   R);
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D u_ = u;
        u_[idim] += epsilon;
        EdEddE_Exforce(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_, ddEdudw_, cg,
                       pex, fex, u_, R);
        std::cout << "dEdu " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdu[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddu  " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEddu.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEdudw.mat[jdim * 3 + idim] << std::endl;
        }
    }
    for (unsigned int idim = 0; idim < 3; idim++) {
        CVector3D dEdu_, dEdw_;
        CMatrix3 ddEddu_, ddEddw_, ddEdudw_;
        double energy_ = 0.0;
        CVector3D w(0, 0, 0);
        w[idim] = epsilon;
        CMatrix3 dR;
        dR.SetRotMatrix_Cartesian(w);
        CMatrix3 R_ = dR * R;
        EdEddE_Exforce(energy_, dEdu_, dEdw_, ddEddu_, ddEddw_, ddEdudw_, cg,
                       pex, fex, u, R_);
        std::cout << "dEdw " << idim << " --> " << (energy_ - energy) / epsilon
                  << " " << dEdw[idim] << std::endl;
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEdudw " << idim << " " << jdim << " --> "
                      << (dEdu_[jdim] - dEdu[jdim]) / epsilon << " "
                      << ddEdudw.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int jdim = 0; jdim < 3; jdim++) {
            std::cout << " ddEddw  " << idim << " " << jdim << " --> "
                      << (dEdw_[jdim] - dEdw[jdim]) / epsilon << " "
                      << ddEddw.mat[jdim * 3 + idim] << std::endl;
        }
    }
}

void CPhysics::CheckDiff_Joint()
{
    double energy = 0.0;
    CVector3D dEdu0(0, 0, 0);
    CVector3D dEdw0(0, 0, 0);
    CVector3D dEdu1(0, 0, 0);
    CVector3D dEdw1(0, 0, 0);
    CMatrix3 ddEdu0du0;
    CMatrix3 ddEdu0dw0;
    CMatrix3 ddEdu0du1;
    CMatrix3 ddEdu0dw1;
    CMatrix3 ddEdw0dw0;
    CMatrix3 ddEdw0du1;
    CMatrix3 ddEdw0dw1;
    CMatrix3 ddEdu1du1;
    CMatrix3 ddEdu1dw1;
    CMatrix3 ddEdw1dw1;
    double epsilon = 1.0e-5;
    const double trans_stiff = 1.0e+3;
    const double rot_stiff = 1.0e+3;
    CVector3D pj = rand_vec(1.0);
    ////
    CVector3D cg0 = rand_vec(1.0);
    CVector3D cp0 = rand_vec(1.0);
    CVector3D u0 = rand_vec(1.0);
    CMatrix3 R0 = rand_rot();
    ////
    CVector3D cg1 = rand_vec(1.0);
    CVector3D cp1 = rand_vec(1.0);
    CVector3D u1 = rand_vec(1.0);
    CMatrix3 R1 = rand_rot();
    ////
    EdEddE_Joint(energy, dEdu0, dEdw0, dEdu1, dEdw1, ddEdu0du0, ddEdu0dw0,
                 ddEdu0du1, ddEdu0dw1, ddEdw0dw0, ddEdw0du1, ddEdw0dw1,
                 ddEdu1du1, ddEdu1dw1, ddEdw1dw1, trans_stiff, rot_stiff, pj,
                 cg0, u0, R0, cg1, u1, R1);
    for (unsigned int kdim = 0; kdim < 3; kdim++) {
        CVector3D dEdu0_, dEdw0_, dEdu1_, dEdw1_;
        CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_,
            ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
        double energy_ = 0.0;
        CVector3D u0_ = u0;
        u0_[kdim] += epsilon;
        EdEddE_Joint(energy_, dEdu0_, dEdw0_, dEdu1_, dEdw1_, ddEdu0du0_,
                     ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_,
                     ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                     trans_stiff, rot_stiff, pj, cg0, u0_, R0, cg1, u1, R1);
        std::cout << "dEdu0: " << kdim << " -->  "
                  << (energy_ - energy) / epsilon << " " << dEdu0[kdim]
                  << std::endl;
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu0du0 " << kdim << " " << idim << " --> "
                      << (dEdu0_[idim] - dEdu0[idim]) / epsilon << " "
                      << ddEdu0du0.mat[idim * 3 + kdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw0du0 " << kdim << " " << idim << " --> "
                      << (dEdw0_[idim] - dEdw0[idim]) / epsilon << " "
                      << ddEdu0dw0.mat[idim * 3 + kdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << "*ddEdu1du0 " << kdim << " " << idim << " --> "
                      << (dEdu1_[idim] - dEdu1[idim]) / epsilon << " "
                      << ddEdu0du1.mat[idim * 3 + kdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << "*ddEdw1du0 " << kdim << " " << idim << " --> "
                      << (dEdw1_[idim] - dEdw1[idim]) / epsilon << " "
                      << ddEdu0dw1.mat[idim * 3 + kdim] << std::endl;
        }
    }
    for (unsigned int jdim = 0; jdim < 3; jdim++) {
        CVector3D dEdu0_, dEdw0_, dEdu1_, dEdw1_;
        CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_,
            ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
        double energy_ = 0.0;
        CVector3D w(0, 0, 0);
        w[jdim] = epsilon;
        CMatrix3 dR;
        dR.SetRotMatrix_Cartesian(w);
        CMatrix3 R0_ = dR * R0;
        EdEddE_Joint(energy_, dEdu0_, dEdw0_, dEdu1_, dEdw1_, ddEdu0du0_,
                     ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_,
                     ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                     trans_stiff, rot_stiff, pj, cg0, u0, R0_, cg1, u1, R1);
        std::cout << "dEdw0: " << jdim << " -->  "
                  << (energy_ - energy) / epsilon << " " << dEdw0[jdim]
                  << std::endl;
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu0dw0 " << jdim << " " << idim << " --> "
                      << (dEdu0_[idim] - dEdu0[idim]) / epsilon << " "
                      << ddEdu0dw0.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw0dw0 " << jdim << " " << idim << " --> "
                      << (dEdw0_[idim] - dEdw0[idim]) / epsilon << " "
                      << ddEdw0dw0.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu1dw0 " << jdim << " " << idim << " --> "
                      << (dEdu1_[idim] - dEdu1[idim]) / epsilon << " "
                      << ddEdw0du1.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw1dw0 " << jdim << " " << idim << " --> "
                      << (dEdw1_[idim] - dEdw1[idim]) / epsilon << " "
                      << ddEdw0dw1.mat[idim * 3 + jdim] << std::endl;
        }
    }

    for (unsigned int jdim = 0; jdim < 3; jdim++) {
        CVector3D dEdu0_, dEdw0_, dEdu1_, dEdw1_;
        CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_,
            ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
        double energy_ = 0.0;
        CVector3D u1_ = u1;
        u1_[jdim] += epsilon;
        EdEddE_Joint(energy_, dEdu0_, dEdw0_, dEdu1_, dEdw1_, ddEdu0du0_,
                     ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_,
                     ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                     trans_stiff, rot_stiff, pj, cg0, u0, R0, cg1, u1_, R1);
        std::cout << "dEdu1: " << jdim << " -->  "
                  << (energy_ - energy) / epsilon << " " << dEdu1[jdim]
                  << std::endl;
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu0du1 " << jdim << " " << idim << " --> "
                      << (dEdu0_[idim] - dEdu0[idim]) / epsilon << " "
                      << ddEdu0du1.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw0du1 " << jdim << " " << idim << " --> "
                      << (dEdw0_[idim] - dEdw0[idim]) / epsilon << " "
                      << ddEdw0du1.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu1du1 " << jdim << " " << idim << " --> "
                      << (dEdu1_[idim] - dEdu1[idim]) / epsilon << " "
                      << ddEdu1du1.mat[idim * 3 + jdim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << "*ddEdw1du1 " << jdim << " " << idim << " --> "
                      << (dEdw1_[idim] - dEdw1[idim]) / epsilon << " "
                      << ddEdu1dw1.mat[idim * 3 + jdim] << std::endl;
        }
    }

    for (unsigned int jdim = 0; jdim < 3; jdim++) {
        CVector3D dEdu0_, dEdw0_, dEdu1_, dEdw1_;
        CMatrix3 ddEdu0du0_, ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_,
            ddEdw0du1_, ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_;
        double energy_ = 0.0;
        CVector3D w(0, 0, 0);
        w[jdim] = epsilon;
        CMatrix3 dR;
        dR.SetRotMatrix_Cartesian(w);
        CMatrix3 R1_ = dR * R1;
        EdEddE_Joint(energy_, dEdu0_, dEdw0_, dEdu1_, dEdw1_, ddEdu0du0_,
                     ddEdu0dw0_, ddEdu0du1_, ddEdu0dw1_, ddEdw0dw0_, ddEdw0du1_,
                     ddEdw0dw1_, ddEdu1du1_, ddEdu1dw1_, ddEdw1dw1_,
                     trans_stiff, rot_stiff, pj, cg0, u0, R0, cg1, u1, R1_);
        std::cout << "dEdw1: " << jdim << " -->  "
                  << (energy_ - energy) / epsilon << " " << dEdw1[jdim]
                  << std::endl;
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu0dw1 " << jdim << " " << idim << " --> "
                      << (dEdu0_[idim] - dEdu0[idim]) / epsilon << " "
                      << ddEdu0dw1.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw0dw1 " << jdim << " " << idim << " --> "
                      << (dEdw0_[idim] - dEdw0[idim]) / epsilon << " "
                      << ddEdw0dw1.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdu1dw1 " << jdim << " " << idim << " --> "
                      << (dEdu1_[idim] - dEdu1[idim]) / epsilon << " "
                      << ddEdu1dw1.mat[jdim * 3 + idim] << std::endl;
        }
        for (unsigned int idim = 0; idim < 3; idim++) {
            std::cout << " ddEdw1dw1 " << jdim << " " << idim << " --> "
                      << (dEdw1_[idim] - dEdw1[idim]) / epsilon << " "
                      << ddEdw1dw1.mat[idim * 3 + jdim] << std::endl;
        }
    }
}

void CPhysics::AddMatrix(Eigen::MatrixXd& M, unsigned int i0, unsigned int j0,
                         const CMatrix3& m, bool isnt_inverse)
{
    if (isnt_inverse) {
        for (unsigned int idim = 0; idim < 3; idim++) {
            for (unsigned int jdim = 0; jdim < 3; jdim++) {
                M(i0 + idim, j0 + jdim) += m.mat[idim * 3 + jdim];
            }
        }
    } else {
        for (unsigned int idim = 0; idim < 3; idim++) {
            for (unsigned int jdim = 0; jdim < 3; jdim++) {
                M(i0 + idim, j0 + jdim) += m.mat[jdim * 3 + idim];
            }
        }
    }
}

void CPhysics::EdEddE_Total(double& E, Eigen::VectorXd& dE,
                            Eigen::MatrixXd& ddE,
                            ////
                            const std::vector<CRigidBody>& aRigidBody,
                            const std::vector<CJoint>& aJoint,
                            ////
                            const CVector3D n, const CVector3D gravity,
                            const double cont_stiff, const double trans_stiff,
                            const double rot_stiff, bool is_friction)
{
    E = 0;
    ddE.setConstant(0);
    dE.setConstant(0);  // u0,w0, u1,w1, ....

    ////
    for (unsigned int irb = 0; irb < aRigidBody.size(); irb++) {
        const CRigidBody& rb = aRigidBody[irb];
        CVector3D cg = rb.cg;
        {
            double e = 0;
            CVector3D du, dw;
            EdEd_Potential(e, du, dw, cg, rb.u, rb.m, gravity);
            E += e;
            for (unsigned int idim = 0; idim < 3; idim++) {
                dE[irb * 6 + 0 + idim] += du[idim];
                dE[irb * 6 + 3 + idim] += dw[idim];
            }
        }
        for (unsigned int icp = 0; icp < rb.aCP.size(); icp++) {
            CVector3D cp = rb.aCP[icp];
            double e = 0;
            CVector3D du, dw;
            CMatrix3 ddu, ddw, dudw;
            if (is_friction) {
                EdEddE_ContactFriction(e, du, dw, ddu, ddw, dudw, cg, cp, rb.u,
                                       rb.R, cont_stiff);
            } else {
                EdEddE_Contact(e, du, dw, ddu, ddw, dudw, cg, cp, rb.u, rb.R,
                               cont_stiff, n);
            }
            E += e;
            for (unsigned int idim = 0; idim < 3; idim++) {
                dE[irb * 6 + 0 + idim] += du[idim];
                dE[irb * 6 + 3 + idim] += dw[idim];
            }
            for (unsigned int idim = 0; idim < 3; idim++) {
                for (unsigned int jdim = 0; jdim < 3; jdim++) {
                    ddE(irb * 6 + 0 + idim, irb * 6 + 0 + jdim) +=
                        ddu.mat[idim * 3 + jdim];
                }
            }
            AddMatrix(ddE, irb * 6 + 0, irb * 6 + 0, ddu, true);
            AddMatrix(ddE, irb * 6 + 0, irb * 6 + 3, dudw, false);
            ////
            AddMatrix(ddE, irb * 6 + 3, irb * 6 + 0, dudw, true);
            AddMatrix(ddE, irb * 6 + 3, irb * 6 + 3, ddw, true);
        }
        for (unsigned int iexf = 0; iexf < rb.aExForce.size(); iexf++) {
            CVector3D pex = rb.aExForce[iexf].first;
            CVector3D fex = rb.aExForce[iexf].second;
            double e = 0;
            CVector3D du, dw;
            CMatrix3 ddu, ddw, dudw;
            EdEddE_Exforce(e, du, dw, ddu, ddw, dudw, cg, pex, fex, rb.u, rb.R);
            E += e;
            for (unsigned int idim = 0; idim < 3; idim++) {
                dE[irb * 6 + 0 + idim] += du[idim];
                dE[irb * 6 + 3 + idim] += dw[idim];
            }
            for (unsigned int idim = 0; idim < 3; idim++) {
                for (unsigned int jdim = 0; jdim < 3; jdim++) {
                    ddE(irb * 6 + 0 + idim, irb * 6 + 0 + jdim) +=
                        ddu.mat[idim * 3 + jdim];
                }
            }
            AddMatrix(ddE, irb * 6 + 0, irb * 6 + 0, ddu, true);
            AddMatrix(ddE, irb * 6 + 0, irb * 6 + 3, dudw, false);
            ////
            AddMatrix(ddE, irb * 6 + 3, irb * 6 + 0, dudw, true);
            AddMatrix(ddE, irb * 6 + 3, irb * 6 + 3, ddw, true);
        }
    }
    for (unsigned int ij = 0; ij < aJoint.size(); ij++) {
        const CJoint& joint = aJoint[ij];
        CVector3D pj = joint.p;
        int irb0 = joint.irb0;
        int irb1 = joint.irb1;
        const CRigidBody& rb0 = aRigidBody[irb0];
        const CRigidBody& rb1 = aRigidBody[irb1];
        {
            double e = 0;
            CVector3D du0, dw0, du1, dw1;
            CMatrix3 du0du0, du0dw0, du0du1, du0dw1, dw0dw0, dw0du1, dw0dw1,
                du1du1, du1dw1, dw1dw1;
            EdEddE_Joint(e, du0, dw0, du1, dw1, du0du0, du0dw0, du0du1, du0dw1,
                         dw0dw0, dw0du1, dw0dw1, du1du1, du1dw1, dw1dw1,
                         trans_stiff, rot_stiff, pj, rb0.cg, rb0.u, rb0.R,
                         rb1.cg, rb1.u, rb1.R);
            E += e;
            for (unsigned int idim = 0; idim < 3; idim++) {
                dE[irb0 * 6 + 0 + idim] += du0[idim];
                dE[irb0 * 6 + 3 + idim] += dw0[idim];
            }
            for (unsigned int idim = 0; idim < 3; idim++) {
                dE[irb1 * 6 + 0 + idim] += du1[idim];
                dE[irb1 * 6 + 3 + idim] += dw1[idim];
            }
            AddMatrix(ddE, irb0 * 6 + 0, irb0 * 6 + 0, du0du0, true);
            AddMatrix(ddE, irb0 * 6 + 0, irb0 * 6 + 3, du0dw0, false);
            AddMatrix(ddE, irb0 * 6 + 0, irb1 * 6 + 0, du0du1, false);
            AddMatrix(ddE, irb0 * 6 + 0, irb1 * 6 + 3, du0dw1, false);
            ////
            AddMatrix(ddE, irb0 * 6 + 3, irb0 * 6 + 0, du0dw0, true);
            AddMatrix(ddE, irb0 * 6 + 3, irb0 * 6 + 3, dw0dw0, true);
            AddMatrix(ddE, irb0 * 6 + 3, irb1 * 6 + 0, dw0du1, false);
            AddMatrix(ddE, irb0 * 6 + 3, irb1 * 6 + 3, dw0dw1, false);
            ////
            AddMatrix(ddE, irb1 * 6 + 0, irb0 * 6 + 0, du0du1, true);
            AddMatrix(ddE, irb1 * 6 + 0, irb0 * 6 + 3, dw0du1, true);
            AddMatrix(ddE, irb1 * 6 + 0, irb1 * 6 + 0, du1du1, true);
            AddMatrix(ddE, irb1 * 6 + 0, irb1 * 6 + 3, du1dw1, false);
            ////
            AddMatrix(ddE, irb1 * 6 + 3, irb0 * 6 + 0, du0dw1, true);
            AddMatrix(ddE, irb1 * 6 + 3, irb0 * 6 + 3, dw0dw1, true);
            AddMatrix(ddE, irb1 * 6 + 3, irb1 * 6 + 0, du1dw1, true);
            AddMatrix(ddE, irb1 * 6 + 3, irb1 * 6 + 3, dw1dw1, true);
        }
    }
}

// solve one iteration
void CPhysics::SolveOneIteration()
/*
(
 std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 double damping_ratio,
 ////
 const CVector3D n,
 const CVector3D gravity,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff)*/
{
    int nDof = (int)aRigidBody.size() * 6;
    double E = 0;
    Eigen::VectorXd dE(nDof);
    Eigen::MatrixXd ddE(nDof, nDof);
    EdEddE_Total(E, dE, ddE, aRigidBody, aJoint, n, gravity, cont_stiff,
                 trans_stiff, rot_stiff, true);
    // std::cout << "energy : " << E << std::endl;
    /////
    Eigen::VectorXd b(nDof);
    Eigen::MatrixXd A(nDof, nDof);
    b = -dE;
    A = ddE;
    for (int i = 0; i < nDof; i++) {
        A(i, i) += damping_ratio;
    }
    Eigen::VectorXd x = A.partialPivLu().solve(b);
    for (unsigned int irb = 0; irb < aRigidBody.size(); irb++) {
        CRigidBody& rb = aRigidBody[irb];
        for (unsigned int idim = 0; idim < 3; idim++) {
            rb.u[idim] += x(irb * 6 + 0 + idim);
        }
        {
            CVector3D w;
            w.x = x(irb * 6 + 3 + 0);
            w.y = x(irb * 6 + 3 + 1);
            w.z = x(irb * 6 + 3 + 2);
            CMatrix3 dR;
            dR.SetRotMatrix_Cartesian(w);
            rb.R = dR * rb.R;
        }
    }
}

void CPhysics::Solve_InterPlane()
/*
(std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 double damping_ratio,
 int nitr,
 ////
 const CVector3D n,
 const CVector3D gravity,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff)
 */
{
    for (int itr = 0; itr < nitr; itr++) {
        // SolveOneIteration(aRigidBody, aJoint, damping_ratio,
        // n,gravity,cont_stiff,trans_stiff,rot_stiff);
        SolveOneIteration();
    }
    ComputeForces();
}

void CPhysics::ComputeForces()
/*
(std::vector<CRigidBody>& aRigidBody,
 std::vector<CJoint>& aJoint,
 ////
 const CVector3D n,
 const double cont_stiff,
 const double trans_stiff,
 const double rot_stiff)
 */
{
    for (CRigidBody& rb : aRigidBody) {
        const int ncp = rb.aCP.size();
        rb.aCForce.resize(ncp);
        for (int icp = 0; icp < ncp; icp++) {
            const CVector3D& cp = rb.aCP[icp];
            CVector3D Rv = rb.R * (cp - rb.cg);
            CVector3D cq = Rv + rb.cg + rb.u;
            rb.aCForce[icp] = ((cq * n) * cont_stiff) * n;
        }
    }
    for (CJoint& joint : aJoint) {
        CVector3D pj = joint.p;
        int irb0 = joint.irb0;
        int irb1 = joint.irb1;
        const CRigidBody& rb0 = aRigidBody[irb0];
        const CRigidBody& rb1 = aRigidBody[irb1];

        CVector3D cg0 = rb0.cg;
        CVector3D u0 = rb0.u;
        CMatrix3 R0 = rb0.R;

        CVector3D cg1 = rb1.cg;
        CVector3D u1 = rb1.u;
        CMatrix3 R1 = rb1.R;

        CVector3D trans_f;  // translation_force
        {
            CVector3D Rv0 = R0 * (pj - cg0);
            CVector3D qj0 =
                Rv0 + cg0 +
                u0;  // after deformation joint pos relative to rigid body 0
            CVector3D Rv1 = R1 * (pj - cg1);
            CVector3D qj1 =
                Rv1 + cg1 +
                u1;  // after deformation joint pos relative to rigid body 1
            trans_f = trans_stiff * (qj0 - qj1);
        }

        CVector3D torque_f;  // rotation_force
        {
            CVector3D av(0, 0, 0);
            for (int i = 0; i < 3; i++) {
                CVector3D r0(R0.mat[0 * 3 + i], R0.mat[1 * 3 + i],
                             R0.mat[2 * 3 + i]);
                CVector3D r1(R1.mat[0 * 3 + i], R1.mat[1 * 3 + i],
                             R1.mat[2 * 3 + i]);
                av += (r0 ^ r1);
            }
            av *= 0.5;
            torque_f = rot_stiff * av;
        }
        joint.linear = trans_f;
        joint.torque = torque_f;
    }
}

void CPhysics::PrintJointForce()
/*
(const std::vector<CRigidBody>& aRigidBody,
 const std::vector<CJoint>& aJoint,
 ////
 const double trans_stiff,
 const double rot_stiff)
 */
{
    std::cout << "force on joint" << std::endl;

    for (unsigned int ij = 0; ij < aJoint.size(); ij++) {
        const CJoint& joint = aJoint[ij];
        CVector3D pj = joint.p;
        int irb0 = joint.irb0;
        int irb1 = joint.irb1;
        const CRigidBody& rb0 = aRigidBody[irb0];
        const CRigidBody& rb1 = aRigidBody[irb1];

        CVector3D cg0 = rb0.cg;
        CVector3D u0 = rb0.u;
        CMatrix3 R0 = rb0.R;

        CVector3D cg1 = rb1.cg;
        CVector3D u1 = rb1.u;
        CMatrix3 R1 = rb1.R;

        CVector3D trans_f;  // translation_force
        {
            CVector3D Rv0 = R0 * (pj - cg0);
            CVector3D qj0 =
                Rv0 + cg0 +
                u0;  // after deformation joint pos relative to rigid body 0
            CVector3D Rv1 = R1 * (pj - cg1);
            CVector3D qj1 =
                Rv1 + cg1 +
                u1;  // after deformation joint pos relative to rigid body 1
            trans_f = trans_stiff * (qj0 - qj1);
        }

        CVector3D torque_f;  // rotation_force
        {
            CVector3D av(0, 0, 0);
            for (unsigned int i = 0; i < 3; i++) {
                CVector3D r0(R0.mat[0 * 3 + i], R0.mat[1 * 3 + i],
                             R0.mat[2 * 3 + i]);
                CVector3D r1(R1.mat[0 * 3 + i], R1.mat[1 * 3 + i],
                             R1.mat[2 * 3 + i]);
                av += (r0 ^ r1);
            }
            av *= 0.5;
            torque_f = rot_stiff * av;
        }

        std::cout << "force of joint: " << ij << std::endl;
        std::cout << "  trans_force:  " << trans_f.x << " " << trans_f.y << " "
                  << trans_f.z << std::endl;
        std::cout << "  torque_force: " << torque_f.x << " " << torque_f.y
                  << " " << torque_f.z << std::endl;
    }
}

void CPhysics::GetHeights(double& eh, double& hmin, unsigned int nH,
                          const std::vector<double>& aXY, double nx, double ny)
{
    const int np = (int)aXY.size() / 2;
    double hmax = 0.0;
    for (int ip = 0; ip < np; ip++) {
        double x = aXY[ip * 2 + 0];
        double y = aXY[ip * 2 + 1];
        double h = nx * x + ny * y;
        if (ip == 0) {
            hmin = h;
            hmax = h;
        } else {
            hmin = (h < hmin) ? h : hmin;
            hmax = (h > hmax) ? h : hmax;
        }
    }
    eh = (hmax - hmin) / (nH - 1);
    hmin += eh * 0.1;
    hmax -= eh * 0.1;
    eh = (hmax - hmin) / (nH - 1);
}

void CPhysics::GetPlateSection(std::vector<CPlateSection>& aPS, double eh,
                               double hmin, int nH,
                               const std::vector<double>& aXY, double nx,
                               double ny)
{
    aPS.clear();
    double tx = -ny;
    double ty = +nx;
    for (int ih = 0; ih < nH; ih++) {
        double h = hmin + ih * eh;
        const int ne = (int)aXY.size() / 2;
        std::set<CPlateSection::CEndPoint> setEP;
        for (int ie = 0; ie < ne; ie++) {
            int ip0 = ie;
            int ip1 = ip0 + 1;
            if (ip1 == ne) {
                ip1 = 0;
            }
            const double h0 = aXY[ip0 * 2 + 0] * nx + aXY[ip0 * 2 + 1] * ny;
            const double h1 = aXY[ip1 * 2 + 0] * nx + aXY[ip1 * 2 + 1] * ny;
            if ((h0 - h) * (h1 - h) > 0) continue;
            double r = (h - h0) / (h1 - h0);
            CPlateSection::CEndPoint ep;
            ep.ie = ie;
            ep.r = r;
            ep.x = (1 - r) * aXY[ip0 * 2 + 0] + r * aXY[ip1 * 2 + 0];
            ep.y = (1 - r) * aXY[ip0 * 2 + 1] + r * aXY[ip1 * 2 + 1];
            ep.t = ep.x * tx + ep.y * ty;
            setEP.insert(ep);
        }
        if (setEP.size() % 2 == 1) continue;
        int nps_h = (int)setEP.size() / 2;
        std::set<CPlateSection::CEndPoint>::iterator itr = setEP.begin();
        for (int ips_h = 0; ips_h < nps_h; ips_h++) {
            CPlateSection ps;
            ps.ep0 = *itr;
            itr++;
            ps.ep1 = *itr;
            itr++;
            ps.ih = ih;
            ps.length = sqrt((ps.ep0.x - ps.ep1.x) * (ps.ep0.x - ps.ep1.x) +
                             (ps.ep0.y - ps.ep1.y) * (ps.ep0.y - ps.ep1.y));
            ps.r = rand() / (RAND_MAX + 1.0);
            ps.g = rand() / (RAND_MAX + 1.0);
            ps.b = rand() / (RAND_MAX + 1.0);
            aPS.push_back(ps);
        }
    }
}

double CPhysics::GetMinDist_LineSegPoint2D(double& t,
                                           const double p[2],  // point
                                           const double s[2],  // source
                                           const double e[2])  // direction
{
    double d[2] = {e[0] - s[0], e[1] - s[1]};
    double ps[2] = {s[0] - p[0], s[1] - p[1]};
    double a = d[0] * d[0] + d[1] * d[1];
    double b = d[0] * ps[0] + d[1] * ps[1];
    t = -b / a;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double p1[2] = {ps[0] + t * d[0], ps[1] + t * d[1]};
    return sqrt(p1[0] * p1[0] + p1[1] * p1[1]);
}

void CPhysics::GetPlateLocalForce(
    std::vector<CPlate::CLocalForce>& alForce,
    ////
    int irb, const std::vector<CRigidBody>& aRigidBody,
    const std::vector<CJoint>& aJoint, const CVector3D& vec_p,
    const CVector3D& vec_t, const CVector3D& vec_b, const CVector3D& vec_n,
    const std::vector<double>& aXY)
{
    alForce.clear();
    const CRigidBody& rb = aRigidBody[irb];
    for (int icp = 0; icp < static_cast<int>(rb.aCP.size()); icp++) {
        CPlate::CLocalForce lf;
        CVector3D pos_cp = rb.aCP[icp];
        lf.lx = (pos_cp - vec_p) * vec_t;
        lf.ly = (pos_cp - vec_p) * vec_b;
        //    double n = (pos_cp-vec_p)*vec_n;
        //    std::cout << "n: " << (pos_cp-vec_p)*vec_n << std::endl;
        double lft = rb.aCForce[icp] * vec_t;
        double lfb = rb.aCForce[icp] * vec_b;
        double lfn = rb.aCForce[icp] * vec_n;
        lf.lfrc = CVector3D(lft, lfb, lfn);
        lf.is_trq = false;
        alForce.push_back(lf);
    }
    for (const auto& exForce : rb.aExForce) {
        const CVector3D p = exForce.first;
        const CVector3D f = exForce.second;
        CPlate::CLocalForce lf;
        lf.lx = (p - vec_p) * vec_t;
        lf.ly = (p - vec_p) * vec_b;
        //    double n = (pos_cp-vec_p)*vec_n;
        //    std::cout << "n: " << (pos_cp-vec_p)*vec_n << std::endl;
        double lft = f * vec_t;
        double lfb = f * vec_b;
        double lfn = f * vec_n;
        lf.lfrc = CVector3D(lft, lfb, lfn);
        lf.is_trq = false;
        alForce.push_back(lf);
    }
    for (const CJoint& j : aJoint) {
        if (j.irb0 != irb && j.irb1 != irb) {
            continue;
        }
        CPlate::CLocalForce lf;
        lf.lx = (j.p - vec_p) * vec_t;
        lf.ly = (j.p - vec_p) * vec_b;
        double lft = j.linear * vec_t;
        double lfb = j.linear * vec_b;
        double lfn = j.linear * vec_n;
        lf.lfrc = CVector3D(lft, lfb, lfn);
        ////
        double ltt = j.torque * vec_t;
        double ltb = j.torque * vec_b;
        double ltn = j.torque * vec_n;
        lf.ltrq = CVector3D(ltt, ltb, ltn);
        lf.is_trq = true;
        if (j.irb1 == irb) {
            lf.lfrc *= -1.0;
            lf.ltrq *= -1.0;
        }
        alForce.push_back(lf);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // if the force is on boundary, register its edge index and ratio on the
    // edge
    // std::cout << "############### " << std::endl;
    double area = fabs(AreaLoop2D(aXY));
    double rep_len = sqrt(area);  // use to judge nearness
    for (CPlate::CLocalForce& clf : alForce) {
        clf.is_boundary = false;
        const double p[2] = {clf.lx, clf.ly};
        int ie_min = -1;
        double r_min = 0;
        double dist_min = 100000;
        for (int ie = 0; ie < static_cast<int>(aXY.size() / 2); ie++) {
            int ip0 = ie;
            int ip1 = ie + 1;
            if (ip1 == static_cast<int>(aXY.size() / 2)) {
                ip1 = 0;
            };
            double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
            double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
            double r, dist = GetMinDist_LineSegPoint2D(r, p, p0, p1);
            if (ie_min == -1 || dist < dist_min) {
                ie_min = ie;
                r_min = r;
                dist_min = dist;
            }
        }
        if (dist_min < rep_len * 0.01) {
            //      std::cout << ilf << "   boundary force  " << dist_min << " "
            //      << rep_len << " " << r_min << "  " << clf.is_trq <<
            //      std::endl;
            clf.is_boundary = true;
            clf.iedge0 = ie_min;
            clf.ratio_edge0 = r_min;
        } else {
            //      std::cout << ilf << "   non boundary force  " << dist_min <<
            //      " " << rep_len << " " << r_min << "  " << clf.is_trq <<
            //      std::endl;
        }
    }
}

//! Area of the Triangle
double CPhysics::TriArea2D(const double v1[2], const double v2[2],
                           const double v3[2])
{
    return 0.5 * ((v2[0] - v1[0]) * (v3[1] - v1[1]) -
                  (v3[0] - v1[0]) * (v2[1] - v1[1]));
}

//! Normalize vector
void CPhysics::NormalizeVector2D(double v[2])
{
    double invlen = 1.0 / sqrt(v[0] * v[0] + v[1] * v[1]);
    v[0] *= invlen;
    v[1] *= invlen;
}

bool CPhysics::IsHitRay(double& r, const double org[2], const double end[2],
                        const double p0[2], const double p1[2])
{
    double a0 = TriArea2D(org, end, p0);
    double a1 = TriArea2D(end, org, p1);
    bool iflg = true;
    if (a0 * a1 < 0) {
        iflg = false;
    }
    r = a0 / (a0 + a1);
    double pm[2] = {(1 - r) * p0[0] + r * p1[0], (1 - r) * p0[1] + r * p1[1]};
    double dot = (pm[0] - org[0]) * (end[0] - org[0]) +
                 (pm[1] - org[1]) * (end[1] - org[1]);
    if (dot > 0 && iflg) {
        return true;
    }
    return false;
}

bool CPhysics::IsForceInsideHalfLoop2D(const std::vector<double>& aXY,
                                       const CPlate::CLocalForce& lf,
                                       const CPlateSection& ps)
{
    const int ie0 = ps.ep0.ie;
    const double re0 = ps.ep0.r;
    const int ie1 = ps.ep1.ie;
    const double re1 = ps.ep1.r;

    double org[2] = {lf.lx, lf.ly};
    double end[2] = {(ps.ep0.x + ps.ep1.x) * 0.5, (ps.ep0.y + ps.ep1.y) * 0.5};

    int ie = ie0;
    int icnt = 1;
    for (;;) {
        unsigned int ip0 = ie;
        unsigned int ip1 = ip0 + 1;
        if (ip1 == aXY.size() / 2) {
            ip1 = 0;
        }
        double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
        double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
        if (ie == ie0) {
            p0[0] = (1 - re0) * p0[0] + re0 * p1[0];
            p0[1] = (1 - re0) * p0[1] + re0 * p1[1];
            if (lf.is_boundary && lf.iedge0 == ie && lf.ratio_edge0 > re0) {
                return true;
            }
        } else if (ie == ie1) {
            p1[0] = (1 - re1) * p0[0] + re1 * p1[0];
            p1[1] = (1 - re1) * p0[1] + re1 * p1[1];
            if (lf.is_boundary && lf.iedge0 == ie && lf.ratio_edge0 < re0) {
                return true;
            }
        } else {
            if (lf.is_boundary && lf.iedge0 == ie) {
                return true;
            }
        }
        if (!lf.is_boundary) {
            double r;
            if (IsHitRay(r, org, end, p0, p1)) {
                icnt++;
            }
        }
        if (ie == ie1) break;
        ie = ie + 1;
        if (ie == static_cast<int>(aXY.size() / 2)) {
            ie = 0;
        }
    }

    if (!lf.is_boundary) {
        if (icnt % 2 == 1) {
            return true;
        }
        return false;
    }
    return false;
}

CVector3D CPhysics::MomentumTri(double& a, CVector3D lg, double p0[2],
                                double p1[2], double org[2])
{
    double cg[2] = {(p0[0] + p1[0] + org[0]) / 3.0,
                    (p0[1] + p1[1] + org[1]) / 3.0};
    double area = TriArea2D(p0, p1, org);
    CVector3D v(cg[0] - org[0], cg[1] - org[1], 0);
    a = area;
    return (v ^ lg) * area;
}

CVector3D CPhysics::MomentumHalfLoop2D(double lx, double ly,
                                       const std::vector<double>& aXY,
                                       const CVector3D& lg,
                                       const CPlateSection& ps)
{
    CVector3D vm(0, 0, 0);
    const int ie0 = ps.ep0.ie;
    const double re0 = ps.ep0.r;
    const int ie1 = ps.ep1.ie;
    const double re1 = ps.ep1.r;

    double org[2] = {lx, ly};

    int ie = ie0;
    // double area = 0;
    for (;;) {
        unsigned int ip0 = ie;
        unsigned int ip1 = ip0 + 1;
        if (ip1 == aXY.size() / 2) {
            ip1 = 0;
        }
        double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
        double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
        if (ie == ie0) {
            p0[0] = (1 - re0) * p0[0] + re0 * p1[0];
            p0[1] = (1 - re0) * p0[1] + re0 * p1[1];
        } else if (ie == ie1) {
            p1[0] = (1 - re1) * p0[0] + re1 * p1[0];
            p1[1] = (1 - re1) * p0[1] + re1 * p1[1];
        }
        ////
        double a;
        vm += MomentumTri(a, lg, p0, p1, org);
        // area += a;
        ////
        if (ie == ie1) break;
        ie = ie + 1;
        if (ie == static_cast<int>(aXY.size() / 2)) {
            ie = 0;
        }
    }
    {
        double p0[2] = {ps.ep1.x, ps.ep1.y};
        double p1[2] = {ps.ep0.x, ps.ep0.y};
        double a;
        vm += MomentumTri(a, lg, p0, p1, org);
        // area += a;
    }
    //  std::cout << area << std::endl;
    return vm;
}

CVector3D CPhysics::MomentumLoop2D(double lx, double ly,
                                   const std::vector<double>& aXY,
                                   const CVector3D& lg)
{
    // double area = 0;
    CVector3D vm(0, 0, 0);
    double org[2] = {lx, ly};
    for (int ie = 0; ie < static_cast<int>(aXY.size() / 2); ie++) {
        unsigned int ip0 = ie;
        unsigned int ip1 = ip0 + 1;
        if (ip1 == aXY.size() / 2) {
            ip1 = 0;
        }
        double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
        double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
        double a;
        vm += MomentumTri(a, lg, p0, p1, org);
        // area += a;
    }
    //  std::cout << area << std::endl;
    return vm;
}

double CPhysics::AreaLoop2D(const std::vector<double>& aXY)
{
    double org[2] = {0, 0};
    double area = 0;
    for (int ie = 0; ie < static_cast<int>(aXY.size() / 2); ie++) {
        unsigned int ip0 = ie;
        unsigned int ip1 = ip0 + 1;
        if (ip1 == aXY.size() / 2) {
            ip1 = 0;
        }
        double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
        double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
        area += TriArea2D(p0, p1, org);
    }
    //  std::cout << area << std::endl;
    return area;
}

bool CPhysics::FirstHitRay(int& ie0, double& r0, double org[2], double dir[2],
                           const std::vector<double>& aXY)
{
    const double end[2] = {org[0] + dir[0], org[1] + dir[1]};
    const double sqdir = sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    std::map<double, std::pair<int, double> > mapDistEdge;
    for (int ie = 0; ie < static_cast<int>(aXY.size() / 2); ie++) {
        unsigned int ip0 = ie;
        unsigned int ip1 = ip0 + 1;
        if (ip1 == aXY.size() / 2) {
            ip1 = 0;
        }
        double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
        double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
        if ((p0[0] - org[0]) * (p0[0] - org[0]) +
                (p0[1] - org[1]) * (p0[1] - org[1]) <
            sqdir * 1.0e-10)
            continue;
        if ((p1[0] - org[0]) * (p1[0] - org[0]) +
                (p1[1] - org[1]) * (p1[1] - org[1]) <
            sqdir * 1.0e-10)
            continue;
        double r0;
        bool res = IsHitRay(r0, org, end, p0, p1);
        if (!res) continue;
        double pm[2] = {(1 - r0) * p0[0] + r0 * p1[0],
                        (1 - r0) * p0[1] + r0 * p1[1]};
        double sqdist = (pm[0] - org[0]) * (pm[0] - org[0]) +
                        (pm[1] - org[1]) * (pm[1] - org[1]);
        mapDistEdge.insert(std::make_pair(sqdist, std::make_pair(ie, r0)));
    }
    if (mapDistEdge.empty()) return false;
    ie0 = mapDistEdge.begin()->second.first;
    r0 = mapDistEdge.begin()->second.second;
    return true;
}

void CPhysics::GetPlateSection_Slit(std::vector<CPlateSection>& aPS, int ie0,
                                    const std::vector<double>& aXY,
                                    bool is_loop_cc)
{
    int ip0 = ie0;
    int ip1 = ip0 + 1;
    if (ip1 == static_cast<int>(aXY.size() / 2)) {
        ip1 = 0;
    }
    double p0[2] = {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]};
    double p1[2] = {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]};
    double eu[2] = {p1[0] - p0[0], p1[1] - p0[1]};
    double len = sqrt(eu[0] * eu[0] + eu[1] * eu[1]);
    eu[0] /= len;
    eu[1] /= len;
    double ev[2];
    if (is_loop_cc) {
        ev[0] = -eu[1];
        ev[1] = +eu[0];
    } else {
        ev[0] = +eu[1];
        ev[1] = -eu[0];
    }
    const double PI = 3.1416;
    /*
    const int nSmpl = 3;
    const double aSmpl[nSmpl][2] = {
      {0.0,-PI*0.500},
      {0.5,+PI*0.000},
      {1.0,+PI*0.500},
    };
     */
    const int nSmpl = 11;
    const double aSmpl[nSmpl][2] = {
        {0.0, -PI * 0.625}, {0.0, -PI * 0.500}, {0.0, -PI * 0.375},
        {0.0, -PI * 0.250}, {0.0, -PI * 0.125}, {0.5, +PI * 0.000},
        {1.0, +PI * 0.125}, {1.0, +PI * 0.250}, {1.0, +PI * 0.375},
        {1.0, +PI * 0.500}, {1.0, +PI * 0.625}};
    ////
    for (int ismpl = 0; ismpl < nSmpl; ismpl++) {
        CPlateSection ps;
        {
            ps.ep0.ie = ie0;
            ps.ep0.r = aSmpl[ismpl][0];
            ps.ep0.x = (1 - ps.ep0.r) * p0[0] + ps.ep0.r * p1[0];
            ps.ep0.y = (1 - ps.ep0.r) * p0[1] + ps.ep0.r * p1[1];
            ps.ep0.t = 0;
        }
        {
            double tht = aSmpl[ismpl][1];
            double org[2] = {ps.ep0.x, ps.ep0.y};
            double dir[2] = {ev[0] * cos(tht) + eu[0] * sin(tht),
                             ev[1] * cos(tht) + eu[1] * sin(tht)};
            bool res = FirstHitRay(ps.ep1.ie, ps.ep1.r, org, dir, aXY);
            if (!res) continue;
            int jp0 = ps.ep1.ie;
            int jp1 = jp0 + 1;
            if (jp1 == static_cast<int>(aXY.size() / 2)) {
                jp1 = 0;
            }
            ps.ep1.x = (1.0 - ps.ep1.r) * aXY[jp0 * 2 + 0] +
                       ps.ep1.r * aXY[jp1 * 2 + 0];
            ps.ep1.y = (1.0 - ps.ep1.r) * aXY[jp0 * 2 + 1] +
                       ps.ep1.r * aXY[jp1 * 2 + 1];
            ps.ep1.t = 0;
        }
        ps.length = sqrt((ps.ep0.x - ps.ep1.x) * (ps.ep0.x - ps.ep1.x) +
                         (ps.ep0.y - ps.ep1.y) * (ps.ep0.y - ps.ep1.y));
        ps.r = rand() / (RAND_MAX + 1.0);
        ps.g = rand() / (RAND_MAX + 1.0);
        ps.b = rand() / (RAND_MAX + 1.0);
        aPS.push_back(ps);
    }
}

void CPhysics::AddWeakSection_Slit(CPlate& plt,  // (in,out)
                                                 ////
                                   const std::vector<double>& aXY,
                                   const CRigidBody& /* rb */,
                                   const CVector3D gravity,
                                   const std::vector<double>& /* aDir2D */,
                                   const int /* nH */, const double max_stress)
{
    CVector3D lg;
    {  // get gravity in local geometry
        double t = gravity * plt.t;
        double b = gravity * plt.b;
        double n = gravity * plt.n;
        lg = CVector3D(t, b, n);
        lg *= plt.rho * plt.thickness;
    }
    const double area = AreaLoop2D(aXY);
    for (int ilf = 0; ilf < static_cast<int>(plt.alForce.size()); ilf++) {
        const CPlate::CLocalForce& lfi = plt.alForce[ilf];
        if (!lfi.is_boundary || !lfi.is_trq) continue;
        const int ie0 = lfi.iedge0;
        // const double r0 = lfi.ratio_edge0;
        // assert( r0 > 0.1 && r0 < 0.9); // joint point should be middle.
        // qDebug() << "CPhysics::AddWeakWection_Slit() - Warning!! Joint point
        // should be middle";
        //    std::cout << ie0 << " " << r0 << std::endl;
        //    std::vector<CPlateSection> aPS;
        std::vector<CPlateSection> aPS;
        GetPlateSection_Slit(aPS, ie0, aXY, area > 0);
        for (CPlateSection& ps : aPS) {
            double lx = (ps.ep0.x + ps.ep1.x) * 0.5;
            double ly = (ps.ep0.y + ps.ep1.y) * 0.5;
            CVector3D bm0(0, 0, 0), bm1(0, 0, 0);
            {
                double dflg = 1;
                if (area > 0) {
                    dflg = -1;
                }
                CVector3D m01 = MomentumLoop2D(lx, ly, aXY, lg);
                CVector3D m0 = MomentumHalfLoop2D(lx, ly, aXY, lg, ps);
                CVector3D m1 = m01 - m0;
                bm0 += m0 * dflg;
                bm1 += m1 * dflg;
            }
            CVector3D ti;
            for (int jlf = 0; jlf < static_cast<int>(plt.alForce.size()); jlf++) {
                const CPlate::CLocalForce& lfj = plt.alForce[jlf];
                CVector3D v(lfj.lx - lx, lfj.ly - ly, 0);
                CVector3D t0 = -(v ^ lfj.lfrc);
                if (lfj.is_trq) {
                    t0 += lfj.ltrq;
                }
                //////
                if (ilf == jlf) {
                    ti = t0;
                } else {
                    bool is_inside0 = IsForceInsideHalfLoop2D(aXY, lfj, ps);
                    if (is_inside0) {
                        bm0 += t0;
                    } else {
                        bm1 += t0;
                    }
                }
            }
            ////
            for (int iside = 0; iside < 2; iside++) {
                CVector3D bm0a = bm0;
                CVector3D bm1a = bm1;
                if (iside == 0) {
                    bm0a += ti;
                } else {
                    bm1a += ti;
                }
                const CVector3D ltrq = (bm1a - bm0a) * 0.5;
                double len =
                    sqrt((ps.ep1.x - ps.ep0.x) * (ps.ep1.x - ps.ep0.x) +
                         (ps.ep1.y - ps.ep0.y) * (ps.ep1.y - ps.ep0.y));
                // double len_med = len;
                double dir_med[2];
                {
                    const int ne = static_cast<int>(aXY.size() / 2);
                    int ie0 = ps.ep0.ie;
                    int ie1 = ps.ep1.ie;
                    int ip0s = ie0;
                    int ip0e = ie0 + 1;
                    if (ip0e == ne) {
                        ip0e = 0;
                    }
                    int ip1s = ie1;
                    int ip1e = ie1 + 1;
                    if (ip1e == ne) {
                        ip1e = 0;
                    }
                    double t0[2] = {aXY[ip0e * 2 + 0] - aXY[ip0s * 2 + 0],
                                    aXY[ip0e * 2 + 1] - aXY[ip0s * 2 + 1]};
                    double t1[2] = {aXY[ip1e * 2 + 0] - aXY[ip1s * 2 + 0],
                                    aXY[ip1e * 2 + 1] - aXY[ip1s * 2 + 1]};
                    NormalizeVector2D(t0);
                    NormalizeVector2D(t1);
                    double vm[2] = {t0[0] - t1[0], t0[1] - t1[1]};
                    NormalizeVector2D(vm);  // medial axis direction
                    double n[2] = {ps.ep1.y - ps.ep0.y, ps.ep0.x - ps.ep1.x};
                    NormalizeVector2D(n);
                    double cos_theta = n[0] * vm[0] + n[1] * vm[1];
                    if (cos_theta < 0) {
                        cos_theta = -cos_theta;
                    }
                    if (cos_theta < cos(3.1415 * 0.25))
                        continue;  // more than sampling direction interval
                    // len_med = len * cos_theta;
                    dir_med[0] = -n[1];
                    dir_med[1] = +n[0];
                }
                double t = plt.thickness;
                ////
                double z0 = t * len * len / 6.0;  // section modulus
                double stress0 = fabs(ltrq.z / z0);
                ////
                double z1 = t * t * len / 6.0;  // section modulus
                double stress1 =
                    fabs(ltrq.x * dir_med[0] + ltrq.y * dir_med[1]) / z1;
                ////
                double stress = (stress0 > stress1) ? stress0 : stress1;
                //        std::cout << fabs(stress0) << " " << fabs(stress1) <<
                //        std::endl;
                ps.ltrq = ltrq;
                ps.weakness = stress / max_stress;
                if (stress > max_stress) {
                    plt.aPS.push_back(ps);
                }
            }
        }
    }
}

void CPhysics::AddWeakSection_Inside(CPlate& plt,  // (in,out)
                                                   ////
                                     const std::vector<double>& aXY,
                                     const CRigidBody& /* rb */,
                                     const CVector3D gravity,
                                     const std::vector<double>& aDir2D,
                                     const int nH, const double max_stress)
{
    /*
    std::cout << std::endl;
    {
      // check
      CVector3D lgf;
      {
        double t = gravity*plt.t;
        double b = gravity*plt.b;
        double n = gravity*plt.n;
        lgf = CVector3D(t,b,n)*rb.m;
      }
      CVector3D sumF(0,0,0);
      {
        for(int ifrc=0;ifrc<plt.alForce.size();ifrc++){
          const CPlate::CLocalForce& lf = plt.alForce[ifrc];
          sumF += lf.lfrc;
        }
        sumF += lgf;
        std::cout << "sumF: " << sumF.x << " " << sumF.y << " " << sumF.z <<
    std::endl;
      }
      CVector3D sumT(0,0,0);
      {
        double cgx = (rb.cg-plt.p)*plt.t;
        double cgy = (rb.cg-plt.p)*plt.b;
        for(int ifrc=0;ifrc<plt.alForce.size();ifrc++){
          const CPlate::CLocalForce& lf = plt.alForce[ifrc];
          CVector3D v(lf.lx-cgx,lf.ly-cgy,0);
          sumT += lf.lfrc^v;
          if( lf.is_trq ){ sumT += lf.ltrq; }
        }
        std::cout << "sumT: " << sumT.x << " " << sumT.y << " " << sumT.z <<
    std::endl;
      }
    }
     */
    //////
    CVector3D lg;
    {  // get gravity in local geometry
        double t = gravity * plt.t;
        double b = gravity * plt.b;
        double n = gravity * plt.n;
        lg = CVector3D(t, b, n);
        lg *= plt.rho * plt.thickness;
    }
    double area = AreaLoop2D(aXY);

    const int ndir = (int)aDir2D.size() / 2;
    for (int idir = 0; idir < ndir; idir++) {
        double nx = aDir2D[idir * 2 + 0];
        double ny = aDir2D[idir * 2 + 1];
        double eh;
        double hmin;
        GetHeights(eh, hmin, nH, aXY, nx, ny);
        std::vector<CPlateSection> aPS;
        GetPlateSection(aPS, eh, hmin, nH, aXY, nx, ny);
        for (CPlateSection& ps : aPS) {
            double lx = (ps.ep0.x + ps.ep1.x) * 0.5;
            double ly = (ps.ep0.y + ps.ep1.y) * 0.5;
            ///
            CVector3D bm0(0, 0, 0);
            CVector3D bm1(0, 0, 0);
            for (const CPlate::CLocalForce& lf : plt.alForce) {
                CVector3D v(lf.lx - lx, lf.ly - ly, 0);
                CVector3D t0 = -(v ^ lf.lfrc);
                if (lf.is_trq) {
                    t0 += lf.ltrq;
                }
                //////
                bool is_inside0 = IsForceInsideHalfLoop2D(aXY, lf, ps);
                if (is_inside0) {
                    bm0 += t0;
                } else {
                    bm1 += t0;
                }
            }
            {  // gravity
                double dflg = 1;
                if (area > 0) {
                    dflg = -1;
                }
                CVector3D m01 = MomentumLoop2D(lx, ly, aXY, lg);
                CVector3D m0 = MomentumHalfLoop2D(lx, ly, aXY, lg, ps);
                CVector3D m1 = m01 - m0;
                bm0 += m0 * dflg;
                bm1 += m1 * dflg;
            }
            //      std::cout << bm0.x << " " << bm0.y << " " << bm0.z << " ";
            //      std::cout << bm1.x << " " << bm1.y << " " << bm1.z <<
            //      std::endl;
            const CVector3D ltrq = (bm1 - bm0) * 0.5;
            double len = sqrt((ps.ep1.x - ps.ep0.x) * (ps.ep1.x - ps.ep0.x) +
                              (ps.ep1.y - ps.ep0.y) * (ps.ep1.y - ps.ep0.y));
            double len_med = len;
            double dir_med[2];
            {
                const int ne = (int)aXY.size() / 2;
                int ie0 = ps.ep0.ie;
                int ie1 = ps.ep1.ie;
                int ip0s = ie0;
                int ip0e = ie0 + 1;
                if (ip0e == ne) {
                    ip0e = 0;
                }
                int ip1s = ie1;
                int ip1e = ie1 + 1;
                if (ip1e == ne) {
                    ip1e = 0;
                }
                double t0[2] = {aXY[ip0e * 2 + 0] - aXY[ip0s * 2 + 0],
                                aXY[ip0e * 2 + 1] - aXY[ip0s * 2 + 1]};
                double t1[2] = {aXY[ip1e * 2 + 0] - aXY[ip1s * 2 + 0],
                                aXY[ip1e * 2 + 1] - aXY[ip1s * 2 + 1]};
                NormalizeVector2D(t0);
                NormalizeVector2D(t1);
                double vm[2] = {t0[0] - t1[0], t0[1] - t1[1]};
                NormalizeVector2D(vm);  // medial axis direction
                double n[2] = {ps.ep1.y - ps.ep0.y, ps.ep0.x - ps.ep1.x};
                NormalizeVector2D(n);
                double cos_theta = n[0] * vm[0] + n[1] * vm[1];
                if (cos_theta < 0) {
                    cos_theta = -cos_theta;
                }
                if (cos_theta < cos(3.1415 * 0.25))
                    continue;  // more than sampling direction interval
                len_med = len * cos_theta;
                dir_med[0] = -n[1];
                dir_med[1] = +n[0];
            }
            double t = plt.thickness;
            ////
            double z0 = t * len_med * len_med / 6.0;  // section modulus
            double stress0 = fabs(ltrq.z / z0);
            ////
            double z1 = t * t * len_med / 6.0;  // section modulus
            double stress1 =
                fabs(ltrq.x * dir_med[0] + ltrq.y * dir_med[1]) / z1;
            ////
            double stress = (stress0 > stress1) ? stress0 : stress1;
            //       std::cout << fabs(stress0) << " " << fabs(stress1) <<
            //       std::endl;
            ps.ltrq = ltrq;
            ps.weakness = stress / max_stress;
            if (stress < max_stress) {
                continue;
            }
            plt.aPS.push_back(ps);
        }
    }
}

void CPhysics::CutSlit()
{
    if (aRigidBody.size() != aPlate.size()) return;

    for (int irb = 0; irb < static_cast<int>(aRigidBody.size()); irb++) {
        const std::vector<double>& aXY = aPlate[irb].aXY;
        std::vector<CSlitData> aSlitData;
        double arealoop = this->AreaLoop2D(aXY);
        for (const CJoint& j : aJoint) {
            if (j.irb0 != irb && j.irb1 != irb) continue;
            const double* ep = (j.irb0 == irb) ? j.jp0 : j.jp1;
            double thickness = (j.irb0 == irb) ? aPlate[j.irb1].thickness
                                               : aPlate[j.irb0].thickness;
            {
                CSlitData sld;
                sld.eu[0] = ep[2] - ep[0];
                sld.eu[1] = ep[3] - ep[1];
                {
                    double len =
                        sqrt(sld.eu[0] * sld.eu[0] + sld.eu[1] * sld.eu[1]);
                    sld.eu[0] /= len;
                    sld.eu[1] /= len;
                }
                if (arealoop > 0) {
                    sld.ev[0] = -sld.eu[1];
                    sld.ev[1] = +sld.eu[0];
                } else {
                    sld.ev[0] = +sld.eu[1];
                    sld.ev[1] = -sld.eu[0];
                }
                sld.ci0[0] = ep[0] - sld.ev[0] * thickness * 0.5;
                sld.ci0[1] = ep[1] - sld.ev[1] * thickness * 0.5;
                sld.ci1[0] = ep[0] + sld.ev[0] * thickness * 0.5;
                sld.ci1[1] = ep[1] + sld.ev[1] * thickness * 0.5;

                FirstHitRay(sld.ie0, sld.r0, sld.ci0, sld.eu, aXY);
                FirstHitRay(sld.ie1, sld.r1, sld.ci1, sld.eu, aXY);

                aSlitData.push_back(sld);
            }
        }
        const int nv = static_cast<int>(aXY.size() / 2);
        std::vector<int> aFlgV0(
            nv,
            0);  // mark if this vertex is between the slit and will be removed
        for (unsigned int isl = 0; isl < aSlitData.size(); isl++) {
            const CSlitData& sld = aSlitData[isl];
            if (sld.ie0 == sld.ie1) continue;
            if (sld.ie0 < sld.ie1) {
                for (int iv = sld.ie0 + 1; iv <= sld.ie1; iv++) {
                    aFlgV0[iv] = 1;
                }
            } else {
                for (int iv = sld.ie0 + 1; iv < nv; iv++) {
                    aFlgV0[iv] = 1;
                }
                for (int iv = 0; iv <= sld.ie1; iv++) {
                    aFlgV0[iv] = 1;
                }
            }
        }
        std::set<CPointOnLoop> setPointOnLoop;
        for (unsigned int isd = 0; isd < aSlitData.size(); isd++) {
            {
                CPointOnLoop p0;
                p0.ie = aSlitData[isd].ie0;
                p0.r = aSlitData[isd].r0;
                p0.is_flont = true;
                p0.islit = isd;
                setPointOnLoop.insert(p0);
            }
            {
                CPointOnLoop p1;
                p1.ie = aSlitData[isd].ie1;
                p1.r = aSlitData[isd].r1;
                p1.is_flont = false;
                p1.islit = isd;
                setPointOnLoop.insert(p1);
            }
        }
        for (int iv = 0; iv < nv; iv++) {
            if (aFlgV0[iv] != 0) {
                continue;
            }
            CPointOnLoop p;
            p.ie = iv;
            p.r = 0;
            p.islit = -1;
            setPointOnLoop.insert(p);
        }
        std::vector<double> aXY1;
        std::set<CPointOnLoop>::iterator itr = setPointOnLoop.begin();
        for (; itr != setPointOnLoop.end(); itr++) {
            const CPointOnLoop& pol = *itr;
            int iv = pol.ie;
            if (pol.islit == -1) {
                aXY1.push_back(aXY[iv * 2 + 0]);
                aXY1.push_back(aXY[iv * 2 + 1]);
            } else {
                unsigned int jv = iv + 1;
                if (iv == nv - 1) {
                    jv = 0;
                }
                double pi[2] = {aXY[iv * 2 + 0], aXY[iv * 2 + 1]};
                double pj[2] = {aXY[jv * 2 + 0], aXY[jv * 2 + 1]};
                double r0 = pol.r;
                double pm[2] = {(1.0 - r0) * pi[0] + r0 * pj[0],
                                (1.0 - r0) * pi[1] + r0 * pj[1]};
                if (pol.is_flont) {
                    const double* ci0 = aSlitData[pol.islit].ci0;
                    aXY1.push_back(pm[0]);
                    aXY1.push_back(pm[1]);
                    aXY1.push_back(ci0[0]);
                    aXY1.push_back(ci0[1]);
                } else {
                    const double* ci1 = aSlitData[pol.islit].ci1;
                    aXY1.push_back(ci1[0]);
                    aXY1.push_back(ci1[1]);
                    aXY1.push_back(pm[0]);
                    aXY1.push_back(pm[1]);
                }
            }
        }
        aPlate[irb].aXY_slit = aXY1;
    }
}

// Setting problem here
void CPhysics::SetExample()
{
    aPlate.clear();
    aRigidBody.clear();
    aJoint.clear();

    {  // making plane 0
        CRigidBody rb;
        rb.cg = CVector3D(0, 1, 0);
        rb.m = 1.0;
        aRigidBody.push_back(rb);
    }
    {  // making plane 1
        CRigidBody rb;
        rb.cg = CVector3D(-1, 0.5, -1);
        rb.aCP.push_back(CVector3D(-1.1, 0, -1.1));
        rb.aCForce.resize(1);
        rb.m = 0.1;
        aRigidBody.push_back(rb);
    }
    {  // making plane 2
        CRigidBody rb;
        rb.cg = CVector3D(-1, 0.5, +1);
        rb.aCP.push_back(CVector3D(-1.1, 0, +1.1));
        rb.aCForce.resize(1);
        rb.m = 0.1;
        aRigidBody.push_back(rb);
    }
    {  // making plane 3
        CRigidBody rb;
        rb.cg = CVector3D(+1, 0.5, -1);
        rb.aCP.push_back(CVector3D(+1.1, 0, -1.1));
        rb.aCForce.resize(1);
        rb.m = 0.1;
        aRigidBody.push_back(rb);
    }
    {  // making plane 4
        CRigidBody rb;
        rb.cg = CVector3D(+1, 0.5, +1);
        rb.aCP.push_back(CVector3D(+1.1, 0, +1.1));
        rb.aCForce.resize(1);
        rb.m = 0.1;
        aRigidBody.push_back(rb);
    }

    {  // joint 0
        CJoint jt;
        jt.p = CVector3D(-1, +1, -1);
        jt.irb0 = 0;
        jt.irb1 = 1;
        aJoint.push_back(jt);
    }
    {  // joint 1
        CJoint jt;
        jt.p = CVector3D(-1, +1, +1);
        jt.irb0 = 0;
        jt.irb1 = 2;
        aJoint.push_back(jt);
    }
    {  // joint 2
        CJoint jt;
        jt.p = CVector3D(+1, +1, -1);
        jt.irb0 = 0;
        jt.irb1 = 3;
        aJoint.push_back(jt);
    }
    {  // joint 3
        CJoint jt;
        jt.p = CVector3D(+1, +1, +1);
        jt.irb0 = 0;
        jt.irb1 = 4;
        aJoint.push_back(jt);
    }
}

static void myGlVertex3d(const CVector3D& v) { ::glVertex3d(v.x, v.y, v.z); }

static void heatmap(float* color, double val)
{
    if (0 <= val && val <= 0.25) {
        color[0] = 0.0;
        color[1] = val * 4.0;
        color[2] = 1.0;
    } else if (0.25 < val && val <= 0.5) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 2.0 - val * 4.0;
    } else if (0.5 < val && val <= 0.75) {
        color[0] = val * 4 - 2.0;
        color[1] = 1.0;
        color[2] = 0.0;
    } else if (0.75 < val && val <= 1.0) {
        color[0] = 1.0;
        color[1] = 4.0 - val * 4.0;
        color[2] = 0.0;
    } else if (1.0 < val) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
    } else {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
    }
}

void CPhysics::DrawGL()
{
    const double small_rad = 0.005;
    const double big_rad = 0.1;
    GLUquadricObj* quadricSphere = gluNewQuadric();

    ::glDisable(GL_LIGHTING);
    for (int irb = 0; irb < static_cast<int>(aRigidBody.size()); irb++) {
        const CRigidBody& rb = aRigidBody[irb];
        CVector3D cg = rb.cg;
        if (is_draw_deformed) {
            cg += rb.u;
        }
        if (aPlate.size() == aRigidBody.size()) {
            const std::vector<double>& aXY = aPlate[irb].aXY_slit;
            const CVector3D& t = aPlate[irb].t;
            const CVector3D& b = aPlate[irb].b;
            const CVector3D& p = aPlate[irb].p;
            const int np = (int)aXY.size() / 2;
            if (irb_draw == -1 || irb_draw == irb) {  // draw_plate_outline
                ::glColor3d(0, 0, 0);
                ::glLineWidth(1);
                ::glBegin(GL_LINE_LOOP);
                for (int ip = 0; ip < np; ip++) {
                    double x1 = aXY[ip * 2 + 0];
                    double y1 = aXY[ip * 2 + 1];
                    CVector3D q = x1 * t + y1 * b + p;
                    myGlVertex3d(q);
                }
                ::glEnd();
            }
            for (const auto& exf : aRigidBody[irb].aExForce) {
                const CVector3D& p = exf.first;
                const CVector3D& f = exf.second;
                ::glColor3d(1, 0, 0);
                ::glLineWidth(3);
                ::glBegin(GL_LINES);
                ::myGlVertex3d(p - f);
                ::myGlVertex3d(p);
                ::glEnd();
            }
            if (is_draw_section) {
                const CPlate& plt = aPlate[irb];
                for (const CPlateSection& ps : plt.aPS) {
                    double x0 = ps.ep0.x;
                    double y0 = ps.ep0.y;
                    double x1 = ps.ep1.x;
                    double y1 = ps.ep1.y;
                    CVector3D q0 = x0 * t + y0 * b + p;
                    CVector3D q1 = x1 * t + y1 * b + p;
                    //          ::glColor3d(ps.r,ps.g,ps.b);
                    int lw = (int)((ps.weakness - 1.0) / 1.0 + 1.0);
                    ::glLineWidth(lw);
                    float color[4];
                    heatmap(color, (ps.weakness - 1.0) / 4.0);
                    ::glColor3fv(color);
                    ::glBegin(GL_LINES);
                    myGlVertex3d(q0);
                    myGlVertex3d(q1);
                    ::glEnd();
                }
            }

            if (is_draw_section_moment) {
                ::glLineWidth(1);
                ::glBegin(GL_LINES);
                const CPlate& plt = aPlate[irb];
                for (const CPlateSection& ps : plt.aPS) {
                    double x0 = ps.ep0.x;
                    double y0 = ps.ep0.y;
                    double x1 = ps.ep1.x;
                    double y1 = ps.ep1.y;
                    CVector3D q0 = x0 * t + y0 * b + p;
                    CVector3D q1 = x1 * t + y1 * b + p;
                    CVector3D qm = (q0 + q1) * 0.5;
                    CVector3D ltrq = ps.ltrq * scale_torque;
                    CVector3D trq =
                        ltrq.x * plt.t + ltrq.y * plt.b + ltrq.z * plt.n;
                    float color[4];
                    heatmap(color, (ps.weakness - 1.0) / 4.0);
                    ::glColor3fv(color);
                    ::myGlVertex3d(qm);
                    ::myGlVertex3d(qm + trq);
                }
                ::glEnd();
            }
        }
        if (is_draw_skeleton) {
            ::glColor3d(0, 1, 0);
            ::glPushMatrix();
            ::glTranslated(cg.x, cg.y, cg.z);
            //::glutWireSphere(0.1, 16, 16);
            gluSphere(quadricSphere, big_rad, 16, 16);
            ::glPopMatrix();
        }
        for (unsigned int icp = 0; icp < rb.aCP.size(); icp++) {
            CVector3D p = rb.aCP[icp];
            if (is_draw_deformed) {
                p = rb.R * (p - rb.cg) + rb.cg + rb.u;
            }
            ::glColor3d(1, 0, 0);
            ::glPushMatrix();
            ::glTranslated(p.x, p.y, p.z);
            //::glutWireSphere(0.02, 16, 16);
            gluSphere(quadricSphere, small_rad, 16, 16);
            ::glPopMatrix();

            if (is_draw_skeleton) {
                ::glLineWidth(5);
                ::glColor3d(0, 0, 0);
                ::glBegin(GL_LINES);
                myGlVertex3d(cg);
                myGlVertex3d(p);
                ::glEnd();
            }

            if (is_draw_force) {
                if (!rb.aCForce.empty()) {
                    CVector3D q = p + scale_force * rb.aCForce[icp];
                    ::glLineWidth(5);
                    ::glColor3d(1, 0, 0);
                    ::glBegin(GL_LINES);
                    myGlVertex3d(p);
                    myGlVertex3d(q);
                    ::glEnd();
                }
            }
        }
    }
    /////////////////////////////////////////////
    for (const CJoint& joint : aJoint) {
        int irb0 = joint.irb0;
        int irb1 = joint.irb1;
        const CRigidBody& rb0 = aRigidBody[irb0];
        const CRigidBody& rb1 = aRigidBody[irb1];

        CVector3D p0 = joint.p;
        CVector3D p1 = joint.p;
        if (is_draw_deformed) {
            p0 = rb0.R * (p0 - rb0.cg) + rb0.cg + rb0.u;
            p1 = rb1.R * (p1 - rb1.cg) + rb1.cg + rb1.u;
        }

        /*
         for(int i=0;i<2;i++){
         int irb = (i==0) ? irb0 : irb1;
         const double* ep = (i==0) ? joint.jp0 : joint.jp1;
         if( irb < 0 || irb >= aPlate.size() ) continue;
         const CVector3D& t = aPlate[irb].t;
         const CVector3D& b = aPlate[irb].b;
         const CVector3D& p = aPlate[irb].p;
         ////
         CVector3D q = ep[2]*t+ ep[3]*b + p;
         //      std::cout << ij << " " << i << "  " << ep[0] << " " << ep[1] <<
         " " << ep[2] << " " << ep[3] << std::endl;
         ::glPushMatrix();
         ::glTranslated(q.x,q.y,q.z);
         //::glutWireSphere(0.02, 16, 16);
         ::glColor3d(0,0,0);
         gluSphere(quadricSphere, small_rad, 16, 16);
         ::glPopMatrix();
         }
         */

        // joint point seen from rb0
        ::glColor3d(0, 0, 1);
        ::glPushMatrix();
        ::glTranslated(p0.x, p0.y, p0.z);
        //::glutWireSphere(0.02, 16, 16);
        gluSphere(quadricSphere, small_rad, 16, 16);
        ::glPopMatrix();

        // joint point seen from rb1
        ::glPushMatrix();
        ::glTranslated(p1.x, p1.y, p1.z);
        //::glutWireSphere(0.02, 16, 16);
        gluSphere(quadricSphere, small_rad, 16, 16);
        ::glPopMatrix();

        CVector3D cg0 = rb0.cg;
        CVector3D cg1 = rb1.cg;
        if (is_draw_deformed) {
            cg0 += rb0.u;
            cg1 += rb1.u;
        }
        if (is_draw_skeleton) {
            ::glLineWidth(5);
            ::glColor3d(0, 0, 0);
            ::glBegin(GL_LINES);
            myGlVertex3d(cg0);
            myGlVertex3d(p0);
            myGlVertex3d(cg1);
            myGlVertex3d(p1);
            ::glEnd();
        }

        if (is_draw_force) {
            ::glLineWidth(5);
            CVector3D q0 = p0 + scale_force * joint.linear;
            CVector3D q1 = p1 - scale_force * joint.linear;
            ::glColor3d(1, 0, 1);
            ::glBegin(GL_LINES);
            myGlVertex3d(p0);
            myGlVertex3d(q0);
            myGlVertex3d(p1);
            myGlVertex3d(q1);
            ::glEnd();

            CVector3D r0 = p0 + scale_torque * joint.torque;
            CVector3D r1 = p1 - scale_torque * joint.torque;
            ::glColor3d(0, 1, 1);
            ::glBegin(GL_LINES);
            myGlVertex3d(p0);
            myGlVertex3d(r0);
            myGlVertex3d(p1);
            myGlVertex3d(r1);
            ::glEnd();
        }
    }

    glLineWidth(1.0f);
    gluDeleteQuadric(quadricSphere);
}

void CPhysics::DrawEnhancedGL(const double scale_fac)
{
    const double small_rad = 0.05 / scale_fac;
    const double big_rad = 0.1 / scale_fac;

    float max_contact_linear = -FLT_MAX;
    float max_joint_linear = -FLT_MAX;
    float max_joint_torque = -FLT_MAX;

    // compute maximal forces
    for (const CRigidBody& rb : aRigidBody) {
        for (unsigned int icp = 0; icp < rb.aCP.size(); icp++) {
            CVector3D p = rb.aCP[icp];
            if (!rb.aCForce.empty()) {
                CVector3D q = p + scale_force * rb.aCForce[icp];
                QVector3D each_f = QVector3D(
                    rb.aCForce[icp].x, rb.aCForce[icp].y, rb.aCForce[icp].z);
                max_contact_linear =
                    qMax(float(each_f.length()), max_contact_linear);
            }
        }
    }

    for (const CJoint& joint : aJoint) {
        QVector3D each_f(joint.linear.x, joint.linear.y, joint.linear.z);
        QVector3D each_t(joint.torque.x, joint.torque.y, joint.torque.z);

        max_joint_linear = qMax(float(each_f.length()), max_joint_linear);
        max_joint_torque = qMax(float(each_t.length()), max_joint_torque);
    }

    // now draw
    for (int irb = 0; irb < static_cast<int>(aRigidBody.size()); irb++) {
        const CRigidBody rb = aRigidBody[irb];
        CVector3D cg = rb.cg;
        if (is_draw_deformed) {
            cg += rb.u;
        }
        if (aPlate.size() == aRigidBody.size()) {
            // const std::vector<double>& aXY = aPlate[irb].aXY_slit;
            const CVector3D& t = aPlate[irb].t;
            const CVector3D& b = aPlate[irb].b;
            const CVector3D& p = aPlate[irb].p;
            // const int np = (int)aXY.size() / 2;

            if (is_draw_section) {
                const CPlate& plt = aPlate[irb];
                for (const CPlateSection& ps : plt.aPS) {
                    double x0 = ps.ep0.x;
                    double y0 = ps.ep0.y;
                    double x1 = ps.ep1.x;
                    double y1 = ps.ep1.y;
                    CVector3D q0 = x0 * t + y0 * b + p;
                    CVector3D q1 = x1 * t + y1 * b + p;

                    float color[4];
                    heatmap(color, (ps.weakness - 1.0) / 4.0);
                    ::glColor3fv(color);

                    QVector3D offset =
                        QVector3D(plt.n.x, plt.n.y, plt.n.z) * plt.thickness;
                    GLutils::DrawCylinder(QVector3D(q0.x, q0.y, q0.z) + offset,
                                          QVector3D(q1.x, q1.y, q1.z) + offset,
                                          small_rad);
                    GLutils::DrawCylinder(QVector3D(q0.x, q0.y, q0.z) + offset,
                                          QVector3D(q0.x, q0.y, q0.z) - offset,
                                          small_rad);
                    GLutils::DrawCylinder(QVector3D(q1.x, q1.y, q1.z) + offset,
                                          QVector3D(q1.x, q1.y, q1.z) - offset,
                                          small_rad);
                    GLutils::DrawCylinder(QVector3D(q0.x, q0.y, q0.z) - offset,
                                          QVector3D(q1.x, q1.y, q1.z) - offset,
                                          small_rad);
                }
            }

            if (is_draw_section_moment) {
                const CPlate& plt = aPlate[irb];
                for (const CPlateSection& ps : plt.aPS) {
                    double x0 = ps.ep0.x;
                    double y0 = ps.ep0.y;
                    double x1 = ps.ep1.x;
                    double y1 = ps.ep1.y;
                    CVector3D q0 = x0 * t + y0 * b + p;
                    CVector3D q1 = x1 * t + y1 * b + p;
                    CVector3D qm = (q0 + q1) * 0.5;
                    CVector3D ltrq = ps.ltrq * scale_torque;
                    CVector3D trq =
                        ltrq.x * plt.t + ltrq.y * plt.b + ltrq.z * plt.n;
                    CVector3D qmtrq = qm + trq;
                    float color[4];
                    heatmap(color, (ps.weakness - 1.0) / 4.0);
                    ::glColor3fv(color);
                    // GLutils::DrawArrowFixedLength(QVector3D(qm.x, qm.y,
                    // qm.z), QVector3D(qmtrq.x, qmtrq.y, qmtrq.z), 0.25f +
                    // 0.75f * ps.weakness); qDebug() << ps.weakness;
                    GLutils::DrawArrowFixedLength(
                        QVector3D(qm.x, qm.y, qm.z),
                        QVector3D(qmtrq.x, qmtrq.y, qmtrq.z),
                        (0.25f +
                         0.75f * qMin(0.25 + (ps.weakness - 1.0) / 4.0, 1.0)) /
                            scale_fac);
                }
            }
        }
        if (is_draw_skeleton) {
            glColor3d(0, 1, 0);
            GLutils::DrawSphere(QVector3D(cg.x, cg.y, cg.z), big_rad);
        }
        for (unsigned int icp = 0; icp < rb.aCP.size(); icp++) {
            CVector3D p = rb.aCP[icp];
            if (is_draw_deformed) {
                p = rb.R * (p - rb.cg) + rb.cg + rb.u;
            }

            if (is_draw_skeleton) {
                // these are ground contacts
                glColor3d(1, 0, 0);
                GLutils::DrawSphere(QVector3D(p.x, p.y, p.z), big_rad);

                // these skeleton lines connect centroids to contact pts
                glColor3f(0.5f, 0.5f, 0.5f);
                GLutils::DrawCylinder(QVector3D(cg.x, cg.y, cg.z),
                                      QVector3D(p.x, p.y, p.z), small_rad);
            }

            if (is_draw_force) {
                if (!rb.aCForce.empty()) {
                    CVector3D q = p + scale_force * rb.aCForce[icp];
                    QVector3D each_f =
                        QVector3D(rb.aCForce[icp].x, rb.aCForce[icp].y,
                                  rb.aCForce[icp].z);
                    const float t = each_f.length() / max_contact_linear;

                    float color[4];
                    heatmap(color, t);
                    ::glColor3fv(color);

                    QVector3D a1(p.x, p.y, p.z);
                    QVector3D a2(q.x, q.y, q.z);
                    GLutils::DrawArrowFixedLength(a1, a2, 1.0f / scale_fac);
                }
            }
        }
    }

    for (const CJoint& joint : aJoint) {
        int irb0 = joint.irb0;
        int irb1 = joint.irb1;
        const CRigidBody& rb0 = aRigidBody[irb0];
        const CRigidBody& rb1 = aRigidBody[irb1];

        CVector3D p0 = joint.p;
        CVector3D p1 = joint.p;
        if (is_draw_deformed) {
            p0 = rb0.R * (p0 - rb0.cg) + rb0.cg + rb0.u;
            p1 = rb1.R * (p1 - rb1.cg) + rb1.cg + rb1.u;
        }

        if (is_draw_skeleton) {
            glColor3d(0, 0, 1);
            // joint point seen from rb0
            GLutils::DrawSphere(QVector3D(p0.x, p0.y, p0.z), big_rad);
            // joint point seen from rb1
            GLutils::DrawSphere(QVector3D(p1.x, p1.y, p1.z), big_rad);
        }

        CVector3D cg0 = rb0.cg;
        CVector3D cg1 = rb1.cg;
        if (is_draw_deformed) {
            cg0 += rb0.u;
            cg1 += rb1.u;
        }
        if (is_draw_skeleton) {
            glColor3f(0.5f, 0.5f, 0.5f);
            GLutils::DrawCylinder(QVector3D(cg0.x, cg0.y, cg0.z),
                                  QVector3D(p0.x, p0.y, p0.z), small_rad);
            glColor3f(0.5f, 0.5f, 0.5f);
            GLutils::DrawCylinder(QVector3D(cg1.x, cg1.y, cg1.z),
                                  QVector3D(p1.x, p1.y, p1.z), small_rad);
        }

        if (is_draw_force) {
            CVector3D q0 = p0 + scale_force * joint.linear;
            CVector3D q1 = p1 - scale_force * joint.linear;
            CVector3D r0 = p0 + scale_torque * joint.torque;
            CVector3D r1 = p1 - scale_torque * joint.torque;

            QVector3D a1(p0.x, p0.y, p0.z);
            QVector3D a2(q0.x, q0.y, q0.z);
            QVector3D a3(p1.x, p1.y, p1.z);
            QVector3D a4(q1.x, q1.y, q1.z);
            QVector3D a5(r0.x, r0.y, r0.z);
            QVector3D a6(r1.x, r1.y, r1.z);

            QVector3D each_f(joint.linear.x, joint.linear.y, joint.linear.z);
            QVector3D each_t(joint.torque.x, joint.torque.y, joint.torque.z);

            const float t1 = each_f.length() / max_joint_linear;
            const float t2 = each_t.length() / max_joint_torque;

            // linear
            // glColor3f(1.0f, 0.0f, 1.0f);
            // glColor3f(1.0f, 1.0f-t1, 1.0f);
            // glColor3f(1.0f, 1.0f-t1, 1.0f-t1);
            float color[4];
            heatmap(color, t1);
            ::glColor3fv(color);
            // GLutils::DrawArrowFixedLength(a1, a2,  t1*0.5f);
            // GLutils::DrawArrowFixedLength(a3, a4,  t1*0.5f);
            GLutils::DrawArrowFixedLength(a1, a2, 0.5f / scale_fac);
            GLutils::DrawArrowFixedLength(a3, a4, 0.5f / scale_fac);

            // torque
            // glColor3f(1.0f-t2, 1.0f, 1.0f);
            // glColor3f(1.0f, 1.0f-t2, 1.0f-t2);
            heatmap(color, t2);
            ::glColor3fv(color);
            // GLutils::DrawCylinderFixedLength(a1, a5, small_rad*0.5f,
            // t2*0.5f); GLutils::DrawCylinderFixedLength(a3, a6,
            // small_rad*0.5f, t2*0.5f);
            GLutils::DrawCylinderFixedLength(a1, a5, small_rad * 0.5f,
                                             0.5f / scale_fac);
            GLutils::DrawCylinderFixedLength(a3, a6, small_rad * 0.5f,
                                             0.5f / scale_fac);
        }
    }

    glLineWidth(1.0f);
}

void CPhysics::DrawFloorGL()
{
    if (is_draw_grid) {
        // draw floor
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glColor3d(0, 0, 0);
        double grid_x_min = -10;
        double grid_x_max = +10;
        double grid_z_min = -10;
        double grid_z_max = +10;
        int ndiv_grid = 30;
        for (int ix = 0; ix < ndiv_grid + 1; ix++) {
            double x0 = (grid_x_max - grid_x_min) / ndiv_grid * ix + grid_x_min;
            ::glVertex3d(x0, 0, grid_z_min);
            ::glVertex3d(x0, 0, grid_z_max);
        }
        for (int iz = 0; iz < ndiv_grid + 1; iz++) {
            double z0 = (grid_z_max - grid_z_min) / ndiv_grid * iz + grid_z_min;
            ::glVertex3d(grid_x_min, 0, z0);
            ::glVertex3d(grid_x_max, 0, z0);
        }
        ::glEnd();
    }
}

bool CPhysics::GetDrawDeformed() { return is_draw_deformed; }

bool CPhysics::GetDrawSkeleton() { return is_draw_skeleton; }

bool CPhysics::GetDrawForce() { return is_draw_force; }

bool CPhysics::GetDrawSection() { return is_draw_section; }

bool CPhysics::GetDrawSectionMoment() { return is_draw_section_moment; }

double CPhysics::GetMaximumStress() { return max_stress; }

void CPhysics::SetDrawDeformed(const bool b) { is_draw_deformed = b; }

void CPhysics::SetDrawSkeleton(const bool b) { is_draw_skeleton = b; }

void CPhysics::SetDrawForce(const bool b) { is_draw_force = b; }

void CPhysics::SetDrawSection(const bool b) { is_draw_section = b; }

void CPhysics::SetDrawSectionMoment(const bool b)
{
    is_draw_section_moment = b;
}

void CPhysics::SetMaximumStress(const double d) { max_stress = d; }
