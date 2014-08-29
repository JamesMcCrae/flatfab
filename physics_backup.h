//
//  rigidbody.h
//  inter_plane
//
//  Created by Nobuyuki Umetani on 12/15/13.
//  Copyright (c) 2013 Nobuyuki Umetani. All rights reserved.
//

#ifndef __inter_plane__rigidbody__
#define __inter_plane__rigidbody__

#include <iostream>
#include <set>
#include <cstdio>
#include <map>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
//#include <GL/glut.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#include "Eigen/Dense"
#else
#include <Eigen/Dense>
#endif

//#if defined(__LINUX__)
//#include <QtCore>
//#endif

#include "physics_vector3d.h"
#include "physics_matrix3d.h"

#include "glutils.h"

// class of rigid body
class CRigidBody
{
public:
  CRigidBody(){
    u = CVector3D(0,0,0);
    R.SetIdentity();
  }
public:
  CVector3D cg; // the center of gravity.
  double m;     // mass of this plate
  std::vector<CVector3D> aCP; // contact point
  std::vector< std::pair<CVector3D,CVector3D> > aExForce; // (position,magnitude) of external forces
  ////
  CVector3D u; // deformation
  CMatrix3  R; // rotation
  std::vector<CVector3D> aCForce;
};

class CJoint {
public:
  CVector3D p; // position
  int irb0;    // id of rigid body
  int irb1;    // id of rigid body
  ////
  double jp0[4]; // joint position
  double jp1[4]; // joint position
  ////
  CVector3D linear;
  CVector3D torque;
};


class CPlateSection
{
public:
  class CEndPoint{
  public:
    double x,y;
    int ie;
    double r;
    double t; // tangent
    bool operator < ( const CEndPoint& ep ) const{ return this->t < ep.t; }
  };
  public:
  CEndPoint ep0;
  CEndPoint ep1;
  ////
  int ih;
  double length;
  ////
  double r,g,b;
  double weakness; // stress scaled with breaking stress
  CVector3D ltrq;
};

class CPlate{
public:
  CVector3D t;
  CVector3D n;
  CVector3D b;
  CVector3D p;
  double thickness; // thickness of plate
  double rho;       // density of plate
  /////
  std::vector<double> aXY;
  std::vector<double> aXY_slit;
  /////
  class CLocalForce{ // force in 2D space
  public:
    double lx;
    double ly;
    CVector3D lfrc;
    ////
    bool is_trq;
    CVector3D ltrq;
    ////
    bool is_boundary;
    int iedge0;
    double ratio_edge0;
  };
  std::vector<CLocalForce> alForce; // local force
  std::vector<CPlateSection> aPS;   // plate section
};

class CSlitData
{
public:
  double ci0[2], ci1[2];
  double eu[2], ev[2];
  int ie0; double r0;
  int ie1; double r1;
  bool is_used;
public:
};

class CPointOnLoop
{
public:
  unsigned int ie;
  double r;
  int islit; // -1:loop vertex index, else:slit point index
  bool is_flont;
public:
  bool operator<(const CPointOnLoop& rhs) const {
    if( this->ie < rhs.ie ) return true;
    if( this->ie > rhs.ie ) return false;
    return this->r < rhs.r;
  }
};
      
class CPhysics {
public:
    
  CPhysics();
  
  //main class interface methods to set up problem
  void AddPlate(const double p[3], const double t[3], const double n[3], const double b[3], const double thickness, const double rho, const std::vector<double> & pts);
  void AddRigidBody(const double centre_of_mass[3],
                    const double mass,
                    const std::vector<double>& contact_points);
  void AddWeightToRigidBody(const double mass_pos[3], const double mass_magnitude[3]);
  void AddJoint(const double position[3], const double jp0[4], const double jp1[4], const int body_index1, const int body_index2);
  void ClearProblem();
  
  
  // Setting A example Problem
  void SetExample();    
  void ReadFile(const char* fname);
  void CutSlit();
  
  //solving methods
  void Solve(){
    Solve_InterPlane();
    Solve_InsidePlane();
  }
  void ComputeForces();

  void PrintJointForce();
  
  //display
  void DrawGL();
  void DrawEnhancedGL();
  void DrawFloorGL();

  //get/set's
    bool GetDrawDeformed();
    bool GetDrawSkeleton();
    bool GetDrawForce();
    bool GetDrawSection();
    bool GetDrawSectionMoment();

    void SetDrawDeformed(const bool b);
    void SetDrawSkeleton(const bool b);
    void SetDrawForce(const bool b);
    void SetDrawSection(const bool b);
    void SetDrawSectionMoment(const bool b);
  
private: // non-static private functions
  void Solve_InterPlane();
  void Solve_InsidePlane();
  void SolveOneIteration();
  //
private:  // static functions
  static void EdEd_Potential(double& energy, CVector3D& dEdu, CVector3D& dEdw, const CVector3D& cg, const CVector3D& u, const double mass, const CVector3D& g);
  static void EdEddE_Exforce(double& energy, CVector3D& dEdu, CVector3D& dEdw, CMatrix3& ddEddu,  CMatrix3& ddEddw, CMatrix3& ddEdudw, const CVector3D& cg,
                      const CVector3D& pex, const CVector3D& fex, const CVector3D& u, const CMatrix3& R);
  static void EdEddE_Contact(double& energy, CVector3D& dEdu, CVector3D& dEdw, CMatrix3& ddEddu, CMatrix3& ddEddw, CMatrix3&  ddEdudw, const CVector3D& cg,
                      const CVector3D& cp, const CVector3D& u, const CMatrix3& R, const double cont_stiff, const CVector3D& n);
  static void EdEddE_ContactFriction(double& energy, CVector3D& dEdu, CVector3D& dEdw, CMatrix3& ddEddu, CMatrix3& ddEddw, CMatrix3& ddEdudw, const CVector3D& cg,
                              const CVector3D& cp, const CVector3D& u, const CMatrix3& R, const double cont_stiff);
  static void EdEddE_Joint(double& energy, CVector3D& dEdu0, CVector3D& dEdw0, CVector3D& dEdu1, CVector3D& dEdw1,
                    CMatrix3& ddEdu0du0,  CMatrix3& ddEdu0dw0, CMatrix3& ddEdu0du1, CMatrix3& ddEdu0dw1, CMatrix3& ddEdw0dw0,
                    CMatrix3& ddEdw0du1,  CMatrix3& ddEdw0dw1, CMatrix3& ddEdu1du1, CMatrix3& ddEdu1dw1, CMatrix3& ddEdw1dw1,
                    const double trans_stiff, const double rot_stiff, const CVector3D& pj,
                    const CVector3D& cg0,  const CVector3D& u0,  const CMatrix3& R0,
                    const CVector3D& cg1,  const CVector3D& u1,  const CMatrix3& R1);
  static void EdEddE_Total(double& E,Eigen::VectorXd& dE, Eigen::MatrixXd& ddE,
                    const std::vector<CRigidBody>& aRigidBody, const std::vector<CJoint>& aJoint,
                    const CVector3D n, const CVector3D gravity,
                    const double cont_stiff, const double trans_stiff, const double rot_stiff, bool is_friction);
  static CVector3D rand_vec(double s);
  static CMatrix3 rand_rot();
  static void AddMatrix(Eigen::MatrixXd& M, unsigned int i0, unsigned int j0, const CMatrix3& m, bool isnt_inverse);
  static double GetMinDist_LineSegPoint2D(double& t, const double p[2],  const double s[2],  const double e[2]);
  static double TriArea2D(const double v1[2], const double v2[2], const double v3[2]);
  static bool IsHitRay(double& r, const double org[2], const double end[2], const double p0[2], const double p1[2]);
  static bool IsForceInsideHalfLoop2D(const std::vector<double>& aXY, const CPlate::CLocalForce& lf, const CPlateSection& ps);
  static CVector3D MomentumTri(double& a,CVector3D lg, double p0[2], double p1[2], double org[2]);
  static CVector3D MomentumHalfLoop2D(double lx, double ly, const std::vector<double>& aXY, const CVector3D& lg, const CPlateSection& ps);
  static CVector3D MomentumLoop2D(double lx, double ly, const std::vector<double>& aXY, const CVector3D& lg);
  static double AreaLoop2D(const std::vector<double>& aXY);
  static bool FirstHitRay(int& ie0, double& r, double org[2], double end[2], const std::vector<double>& aXY);

  static void GetHeights
  (double &eh,
   double& hmin,
   unsigned int nH,
   const std::vector<double>& aXY,
   double nx, double ny);

  static void GetPlateSection
  (std::vector<CPlateSection>& aPS,
   double eh,double hmin, int nH,
   const std::vector<double>& aXY,
   double nx, double ny);
  
  static void GetPlateSection_Slit
  (std::vector<CPlateSection>& aPS,
   int ie0,
   const std::vector<double>& aXY,
   bool is_loop_cc);  
  
  // solve one iteration
  static void SolveOneIteration
  (std::vector<CRigidBody>& aRigidBody,
   std::vector<CJoint>& aJoint,
   ////
   double damping_ratio,
   ////
   const CVector3D n,
   const CVector3D gravity,
   const double cont_stiff,
   const double trans_stiff,
   const double rot_stiff);
  
  
  static void GetPlateLocalForce
  (std::vector<CPlate::CLocalForce>& alForce,
   ////
   int irb,
   const std::vector<CRigidBody>& aRigidBody,
   const std::vector<CJoint>& aJoint,
   const CVector3D& vec_p,
   const CVector3D& vec_t,
   const CVector3D& vec_b,
   const CVector3D& vec_n,
   const std::vector<double>& aXY);
  
  
  static void AddWeakSection_Inside
  (CPlate& plate, // (in,out)
   ///
   const std::vector<double>& aXY,
   const CRigidBody& rb,
   const CVector3D gravity,
   ////
   const std::vector<double>& aDir2D,
   const int nH,
   const double max_stress);
    
  static void AddWeakWection_Slit
  (CPlate& plt, // (in,out)
   ////
   const std::vector<double>& aXY,
   const CRigidBody& rb,
   const CVector3D gravity,
   const std::vector<double>& aDir2D,
   const int nH,
   const double max_stress);

  static void CheckDiff_Contact();
  static void CheckDiff_Joint();
  static void CheckDiff_Exforce();
  static void CheckDiff_ContactFriction();  
public:
  //members
  std::vector<CRigidBody> aRigidBody; // array of rigid body
  std::vector<CJoint> aJoint; // array of joint
  std::vector<CPlate> aPlate;
  
  CVector3D n; // normal direction of floor (should be an unit vector)
  CVector3D gravity; // gravity
  double cont_stiff; // contact_stiffness (insensitive)
  double trans_stiff; // joint_translation_stiffness (insensitive)
  double rot_stiff; // joint_rotation_stiffness (insensitive)
  
  int nitr;
  double damping_ratio;
  ////
  bool is_draw_force;
  bool is_draw_skeleton;
  bool is_draw_deformed;
  bool is_draw_section;
  bool is_draw_grid;
  bool is_draw_section_moment;
  double scale_force;
  double scale_torque;
  int irb_draw; // if( irb_draw==-1 ) draw all the rigid bodies
  
  // inside plate analysis
  bool is_slit;
  std::vector<double> aDir2D; // number of sampling direction
  int nHeight; // number of sampling heights
  double max_stress; // maximum_tensil_stress
};
      



#endif /* defined(__inter_plane__rigidbody__) */
