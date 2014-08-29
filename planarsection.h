#ifndef PLANARSECTION_H
#define PLANARSECTION_H

#include <QtOpenGL>
#include <QVector2D>

#ifdef __APPLE__
    #include <OpenGL/glu.h>
#else
    #include <GL/glu.h>
#endif

#include "beziercurve.h"
#include "bezierfit.h"
#include "triangulate.h"
#include "glutils.h"


struct ExternalWeight {

    QVector2D p;
    float mass_kg;

    float rad;
    float height;

};

class PlanarSection
{

public:

    PlanarSection();

    //PlanarSection & operator=(const PlanarSection & rhs);

    BezierCurve & GetCurve(const int i);
    int GetNumCurves() const;
    void SetCurve(const int i, const BezierCurve & b);

    QList <BezierCurve> & GetCurves();
    void SetCurves(const QList <BezierCurve> & curves);

    void AddNewCurve();
    void RemoveCurve(const int i);
    void RemoveHoles();

    void SetNewBasis(const QVector3D & new_t, const QVector3D & new_n, const QVector3D & new_b, const QVector3D & new_p);

    //get/set for TNB frame and P
    void SetP(const QVector3D & v);    
    void SetT(const QVector3D & v);
    void SetN(const QVector3D & v);
    void SetB(const QVector3D & v);

    QVector3D P() const;
    QVector3D T() const;
    QVector3D N() const;
    QVector3D B() const;

    void CreateSquare(const float size);
    void CreateRectangle(QVector2D min_v, QVector2D max_v);
    void CreateCircle(QVector2D centre, const float radius);
    void CreateRadial(QVector2D centre, const float base_rad, const int num_sectors, const float radii[9]);
    void CreateRadialHoles(QVector2D centre, const float base_rad, const int num_sectors, const float radii[9]);

    void AddWeight(const QVector2D & pos_2d, const float mass);
    void AddWeightAtMousePos(const QVector2D & mouse_pos, const float mass);
    void RemoveWeights();
    int GetNumWeights();
    QVector3D GetWeightPosition(const int i);
    QVector3D GetWeightForce(const int i);
    float GetWeightMass(const int i);

    void MoveP(const QVector3D & v); //note: this updates the bezier curve so it remains fixed in space
    void Scale(const float x, const float y);
    void Rotate(const float angle_rad);

    void FlipN(); //note: this updates the bezier curve so it remains fixed in space
    void FlipTB(); //note: this updates the bezier curve so it remains fixed in space

    //drawing methods
    void Draw();
    void DrawForPicking(const int pickval);
    void DrawCurve();
    void DrawCurveControlPoints(const float cam_width);
    void DrawCurveControlPolygon();
    void DrawTris();
    void DrawSlab();
    void DrawSlabCurves();
    void DrawSketch();
    void DrawTNBFrame();
    void DrawShadow();
    void DrawXZPlaneLine();
    void DrawDeadzone(const QVector3D & slot_start, const QVector3D & slot_end, const float deadzone_radius);
    static void DrawConvexHull(const QList <QVector3D> & hull, const bool draw_green);
    void DrawWeights();

    //update method
    void Update(const int c, const double error_tolerance);
    void UpdateCurveTrisSlab();

    //file load/save methods
    static void SaveHeaderToFile(const QList <PlanarSection> & sections, QTextStream & ofs);
    void SaveToFile(QTextStream & ofs);
    static int LoadHeaderFromFile(QTextStream & ifs);
    void LoadFromFile(QTextStream & ifs);

    static void ComputePacking(const QList <PlanarSection> & sections, const double width_height_ratio, QVector <QVector2D> & pack_pts, QVector2D & bbox);
    static void SaveToSVG(const QList <PlanarSection> & sections, QTextStream & ofs, const double metres_per_unit, const int dpi, const double width_height_ratio, const bool use_numeric_labels, const QList <QString> & numeric_labels);
    static void SaveToDXF(const QList <PlanarSection> & sections, QTextStream & ofs, const double metres_per_unit, const int dpi, const double width_height_ratio, const bool use_numeric_labels, const QList <QString> & numeric_labels);
    static void DXFWriteLine(QTextStream & ofs, const QVector3D & p1, const QVector3D & p2, const int layer_index, const int colour_index);
    static void DXFWriteText(QTextStream & ofs, const QVector3D & p1, const QString & text, const int layer_index, const int colour_index, const float height);
    static void SaveConnectionsToOBJ(const QList <PlanarSection> & sections, QTextStream & ofs);

    //intersection and slot related
    void GetSlots(const PlanarSection & other, QList <QVector2D> & my_slots, QList <QList <QVector2D> > & my_slots_rect) const;

    bool GetIntersectionSlotLine(const PlanarSection & other, QVector2D & slot_start, QVector2D & slot_end) const; //1.  do the planes intersect, what axis?
    void GetContourIntersections(const PlanarSection & other, QList <QVector3D> & isecs, QList <bool> & isecs_which) const;  //2.  give me points where EITHER plane intersects contour

    bool LinesIntersectContour(const QList <QVector2D> & lines) const;
    void GetBoundingBox2D(QVector2D & min_v, QVector2D & max_v) const;
    void GetBoundingInterval(const QVector2D & direction, float & min_v, float & max_v) const;
    void ComputeContourLineIntersects(const QVector2D & lp, const QVector2D & ld, QList <QVector2D> & intersects) const;

    //physical feasibility
    static void ComputeIntersectors(const QList <PlanarSection> & sections, const int which_section, QVector <bool> & isecs);
    static void ComputeIntersectionGraph(const QList <PlanarSection> & sections, QVector <QVector <bool> > & graph);
    //static void ComputeCycles(const QVector <QVector <bool> > & graph, QList <QList <int> > & cycles, QList <QList <int> > & all_paths);
    static void TestConnectedness(QList <PlanarSection> & sections, QList <QList <int> > & all_paths);

    //skinning
    void SetPartOfCycle(const bool b);
    bool PartOfCycle() const;
    static bool CycleAssemblable(const QList <PlanarSection> & section, const QList <int> & cycle);

    //connected get/set's
    void SetConnected(const bool b);
    bool Connected() const;

    void SetQuality(const int i);

    void SketchClear(const int c);
    void SketchAdd(const int c, const QVector2D & v);
    void SketchSetEditing(const bool b);
    int SketchNumPoints() const;
    bool SketchEditing();
    void SnapSketchToTemplateCut(const QList <QVector3D> & cut_segments, const float snap_dist_3d);

    bool CurveClockwise(const int c);
    void SketchSymmetryTest();

    bool AddMouseRayIntersect(const int c, const QVector2D & v);
    bool MouseRayIntersect(const QVector2D & v, QVector3D & intersect);
    bool RayIntersect(const QVector3D & p0, const QVector3D & dir, QVector3D & intersect);
    bool MouseOutsideDeadzone(const QVector2D & v, const QVector3D & slot_start, const QVector3D & slot_end, const float deadzone_radius);

    void SelectMouseRayIntersect(const QVector2D & v, const float cam_width);
    void MoveWeightMouseRayIntersect(const QVector2D & v);
    void MoveCtrlPointMouseRayIntersect(const QVector2D & v, const bool keep_g1);

    bool IsCtrlPointSelected();
    int SelectedCtrlPoint();
    void UnselectCtrlPoint();
    void SetCtrlPointsAboveXZPlane();
    void InsertCtrlPoint();
    void DeleteCtrlPoint();

    bool IsWeightSelected();
    void InsertWeight(const QVector2D & mouse_pos, const float mass);
    void DeleteWeight();

    void MirrorPlanarSectionX(PlanarSection & new_section) const;
    void MirrorPlanarSectionZ(PlanarSection & new_section) const;
    void RotatePlanarSectionY(PlanarSection & new_section) const;

    QVector2D GetPoint2D(const QVector3D & v) const;
    QVector3D GetPoint3D(const QVector2D & v) const;

    QVector3D GetTangent3D(const QVector2D & v) const;

    void ContourVertices3D(QList <QVector3D> & contour_verts) const;

    const QList <QVector2D> & SliceTriangles() const;
    const QList <QVector3D> & SlabVertices() const;
    const QList <QVector3D> & SlabNormals() const;

    float SlabThickness() const;
    void SetSlabThickness(const float f);

    //will it stand related methods...
    static QVector3D Centroid(const QList <PlanarSection> & sections);
    QVector2D Centroid() const;
    float SignedArea() const;
    QVector3D GetCentroid3D() const;

    void GetXZPlanePoints(QList <QVector3D> & plane_pts) const;
    static void GetXZPlanePoints(const QList <PlanarSection> & sections, QList <QVector3D> & plane_pts);
    static void GetXZConvexHull(const QList <PlanarSection> & sections, QList <QVector3D> & hull);
    static bool GetXZPointInsideConvexHull(const QList <QVector3D> & hull, const QVector3D & p);

    void SplitAlongLine(const QVector2D & split_p, const QVector2D & split_d, QList <PlanarSection> & split_sections);
    bool IsIntersectingSection(const PlanarSection & other) const;

    void GetCurveBetweenSections(const int my_index, const PlanarSection & section1, const int section1_index, const PlanarSection & section2, const int section2_index, BezierCurve & curve);
    void GetCurveAroundSection(const PlanarSection & section, BezierCurve & curve);
    void GetSectionsAlongCurve(const PlanarSection & start_section, const BezierCurve & curve, const int num_frames, QList <PlanarSection> & sections);

    float GetPointDistance(const QVector3D & v);

private:  

    void WidenSlot(const float factor, QList <QVector2D> & slot);

    void ComputeCentroidAndArea();
    void ComputeSlab();

    static void ComputeSlitNumericLabels(const QList <PlanarSection> & sections, const QList <QString> & numeric_labels, QVector <QVector <QString> > & graph_labels);

    //plane spatial parameters
    QVector3D p; //origin
    QVector3D t; //tangent
    QVector3D b; //binormal
    QVector3D n; //normal

    bool editing_sketch;
    bool part_of_cycle;
    bool connected;

    //the input stroke    
    //there will be some "path" related members
    //probably some planar bezier thing, but maybe have control point "linkage"
    BezierFit bez_fit;
    QList <QList <QVector2D> > sketch_pts;
    QList <BezierCurve> bez_curve;

    QVector2D centroid;
    float area;

    //there will be some triangulated geometry members
    QList <QVector2D> tris; //planar surface
    QList <QVector3D> slab_vert; //slab
    QList <QVector3D> slab_norm;

    float slab_thickness;   
    int quality_samples;

    //there will be some topological "linking" or "grouping" parameters here
    static QList <QVector3D> convex_hull;

    //external weights attached to this plane
    QList <ExternalWeight> weights;
    int selected_weight;    

    //display lists for some drawing
    //bool update_slab_disp_list;
    //GLuint slab_disp_list;
    //bool update_slabcurves_disp_list;
    //GLuint slabcurves_disp_list;

};

#endif // PLANARSECTION_H
