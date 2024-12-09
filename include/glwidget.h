#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtGui>
#include <QtOpenGL>

#include "contourgraph.h"
#include "glutils.h"
#include "physics.h"
#include "pivotcamera.h"
#include "planarsection.h"
#include "transformwidget.h"
#include "tree.h"

enum GestureState
{                           // states (all with LMB held down, from 1 onward)
    STATE_NONE,             // 0.  Awaiting mouse drag
    STATE_SLOT,             // 1.  CHOOSING LINE ON EXISTING PLANE
    STATE_CAM_TRANSLATE,    // 2.  MOVING LINE ENDPOINT AROUND ON-SCREEN
    STATE_DEADZONE,         // 3  DEADZONE
    STATE_CURVE,            // 4.  DEFINING A CURVE
    STATE_ORBIT,            // 5.  orbiting camera
    STATE_MANIP_CTRLPOINT,  // 6. manipulating control point
    STATE_RECURSIVE_SETUP_SLOT,
    STATE_RECURSIVE_SLOT,
    STATE_MANIP_WEIGHT,
    STATE_RESKETCH_CURVE,
    STATE_TRANSFORM_WIDGET,  // this state replaces TRANSLATE, ROTATE, SCALE,
                             // TRANSLATE_NORMAL
    STATE_ADD_HOLE,
    STATE_DIMENSIONING_FIRST,
    STATE_DIMENSIONING_SECOND,
    STATE_PEN_POINT,
    STATE_PEN_DRAG
};

enum LastOperation
{
    OP_NONE,
    OP_SELECTION,
    OP_UNDO,
    OP_ADD_PLANE,
    OP_DELETE_PLANE,
    OP_ADD_WEIGHT,
    OP_DELETE_WEIGHT,
    OP_DELETE_WEIGHTS,
    OP_ADD_CTRLPOINT,
    OP_DELETE_CTRLPOINT,
    OP_ADD_HOLE,
    OP_DELETE_HOLES,
    OP_RESKETCH_CURVE,
    OP_MANIP_TRANSFORM,
    OP_MANIP_WEIGHT,
    OP_MANIP_CTRLPOINT,
    OP_MANIP_DIMENSIONING_TOOL,
    OP_MANIP_SNAP_MAJOR_AXIS,
    OP_COPY_MIRRORX,
    OP_COPY_ROTATEY,
    OP_COPY_MIRRORZ,
    OP_COPY_DUPLICATE,
    OP_GENERATE_BLEND,
    OP_GENERATE_BRANCHING,
    OP_GENERATE_SURFACE_FACETS,
    OP_GENERATE_GRID,
    OP_GENERATE_LINEAR,
    OP_GENERATE_MAKE_CIRCLE,
    OP_GENERATE_MAKE_RADIAL,
    OP_GENERATE_MAKE_RADIAL_HOLE,
    OP_GENERATE_MAKE_RECTANGLE,
    OP_GENERATE_REVOLVE,
    OP_GENERATE_SLICES
};

// Tool states
enum ToolState
{
    TOOLSTATE_DEFAULT,
    TOOLSTATE_TRANSFORMING,
    TOOLSTATE_GENERATE
};

// Generate states
enum GenerateState
{
    GENSTATE_BLEND,
    GENSTATE_LINEAR,
    GENSTATE_REVOLVE,
    GENSTATE_SLICES,
    GENSTATE_GRID,
    GENSTATE_BRANCH,
    GENSTATE_SURFACE_FACETS,
    GENSTATE_NUM  // note: enum declaration guaranteed to be in order: 0,
                  // 1, 2... therefore this enum is the # of previously defined
                  // enums
};

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget();
    ~GLWidget();

    QWidget *GetEditWidget();
    QWidget *GetGenerateWidget();
    QWidget *GetGuidesWidget();
    QWidget *GetPhysicsWidget();
    QWidget *GetViewWidget();

    CPhysics &GetPhysics();

    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

    bool GetDoMagneticCuts();
    bool GetDoLocalSymmetry();
    bool GetPenModeOn();
    bool GetDoCyclesTest();
    bool GetDoStabilityTest();
    bool GetDoPhysicsTest();
    bool GetDoConnectedTest();
    bool GetShowTNBFrames();
    bool GetShowShadow();
    bool GetShowTemplates();

    float GetGenerateGridSizeX();
    float GetGenerateGridSizeY();
    float GetGenerateGridStapleSize();
    float GetGenerateLinearSpacing();
    bool GetGenerateLinearScaleX();
    bool GetGenerateLinearScaleY();
    float GetGenerateSlicesSpacing();
    bool GetGenerateSlicesX();
    bool GetGenerateSlicesY();
    bool GetGenerateSlicesZ();
    int GetGenerateBlendSections();
    int GetGenerateRevolveSections();
    float GetGenerateBranchingScaleChild();
    int GetGenerateRadialSectors();
    bool GetGenerateSurfaceFacetsTeeth();

    double GetPhysicsNewWeightMass();
    double GetPhysicsMaterialDensity();
    double GetPhysicsMaximumStress();

    bool IsSectionSelected() const;

    QString GetOpenFilename() const;

protected:
    void initializeGL();
    void paintEvent(QPaintEvent *);
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

    void ResetCamera();
    void UpdateCamera();

    void DrawGroundPlane();
    void DrawSymmetryPlanes();
    void DrawSection(const int i);
    void DrawInstructions(QPainter &painter);
    void DrawMarkers();
    void DrawSlot(QVector3D start, QVector3D end);
    void DrawDimensionTool();
    void DrawTemplateSurface();
    void DrawTemplateImage();

    void DoCyclesConnectedTest();
    void DoStabilityTest();
    void DoPhysicsTest();

    int PickSection(const QVector2D &mouse_pos, const bool only_test_selected);
    int PickTransformWidgetElement(const QVector2D &mouse_pos);

signals:

public slots:

    void SetUnitsToInches();
    void SetUnitsToCentimetres();

    void SetViewIso1();
    void SetViewIso2();
    void SetViewIso3();
    void SetViewIso4();
    void SetViewX();
    void SetViewY();
    void SetViewZ();
    void SetViewPart();

    void SetXSymmetry(bool b);
    void SetYSymmetry(bool b);
    void SetZSymmetry(bool b);

    void SetTemplateImageX(const int i);
    void SetTemplateImageY(const int i);
    void SetTemplateImageRotate(const int i);
    void SetTemplateImageScale(const int i);
    void SetTemplateImageFlipX(const bool b);

    void SetDoMagneticCuts(const bool b);
    void SetDoLocalSymmetry(const bool b);
    void SetPenModeOn(const bool b);
    void SetDoCyclesTest(const bool b);
    void SetDoStabilityTest(const bool b);
    void ToggleDoPhysicsTest();
    void SetDoPhysicsTest(const bool b);
    void SetDoConnectedTest(const bool b);
    void SetShowTNBFrames(const bool b);
    void SetShowShadow(const bool b);
    void SetShowTemplates(const bool b);

    void ToggleShowTNBFrames();
    void ToggleDoCyclesTest();
    void ToggleDoStabilityTest();
    void ToggleShowShadow();
    void ToggleDoConnectedTest();

    void SetSlabThickness(const double d);
    void SetCalibrationFactor(const int i);
    void SetQualitySamples(const int i);
    void SetRotationDuration(const int i);
    void SetRotationAngle(const int i);
    void SetGenerateGridSizeX(const int i);
    void SetGenerateGridSizeY(const int i);
    void SetGenerateGridStapleSize(const int i);
    void SetGenerateLinearSpacing(const int i);
    void SetGenerateLinearScaleX(const bool b);
    void SetGenerateLinearScaleY(const bool b);
    void SetGenerateSlicesSpacing(const int i);
    void UpdateGenerateSlicesSpacing();
    void SetGenerateSlicesX(const bool b);
    void SetGenerateSlicesY(const bool b);
    void SetGenerateSlicesZ(const bool b);
    void SetGenerateBlendSections(const int i);
    void SetGenerateRevolveSections(const int i);
    void SetGenerateBranchingScaleChild(const int i);
    void SetGenerateRadialSectors(const int i);
    void SetGenerateSurfaceFacetsTeeth(const bool b);

    void ToggleDrawDeformed();
    void ToggleDrawSkeleton();
    void ToggleDrawForce();
    void ToggleDrawSection();
    void ToggleDrawMoments();

    void ToggleDrawTemplates();

    void Undo(const LastOperation op);
    void Redo();
    void DeleteSelected();
    void Transform();
    void CopyMirrorX();
    void CopyRotateY();
    void CopyMirrorZ();
    void CopyDuplicate();
    void SnapToMajorAxis();
    void AddHoleBoundary();
    void RemoveHolesBoundary();
    void ResketchCurve();
    void CreateSurfacePatches();
    void DeleteSurfacePatches();
    void StartDimensioningTool();

    void NewPlaneSketch();
    void LoadTemplateCurve();
    void LoadTemplateImage();
    void LoadTemplateOBJ();
    bool LoadPlaneSketch();
    bool SavePlaneSketch();
    void SaveSliceOBJ();
    void SaveSlabOBJ();
    void SaveFlattenedSlabOBJ();
    void SaveSurfaceOBJ();
    void SaveSVG();
    void SaveDXF();
    void SaveCalibration();
    void SavePhysicsOutput();

    // procedural modelling extras
    void DoGenerateBranchingSetRoot();
    void DoGenerateBranching();
    void DoGenerateMakeCircle();
    void DoGenerateMakeRectangle();
    void DoGenerateMakeRadialHole();
    void DoGenerateSurfaceFacets();

    bool ShowGenerate();
    void AcceptGenerate();
    void CancelGenerate();

    void AcceptPenCurve();
    void CancelPenCurve();

    void StartGenerateLinear();
    void StartGenerateBlend();
    void StartGenerateRevolve();
    void StartGenerateSlices();
    void StartGenerateGrid();
    //    void StartGenerateBranch();

    void SetSelectedAsRadial();
    void RemoveRadial();

    // physics stuff
    void DoPhysicsAddWeight();
    void DoPhysicsRemoveWeights();
    void SetPhysicsNewWeightMass(const double d);
    void SetPhysicsMaterialDensity(const double d);
    void SetPhysicsMaximumStress(const double d);

    void SetDrawDeformed(bool b);
    void SetDrawSkeleton(bool b);
    void SetDrawForce(bool b);
    void SetDrawSection(bool b);
    void SetDrawSectionMoment(bool b);

private:
    void SetupPlanarSection(PlanarSection &p);

    void SetMetresPerUnit(const double d);

    void ScalePlanarSections(const float s);

    void UpdateDraw();
    void UpdateAllTests();

    void AddToUndoList(const LastOperation op);
    void ClearAll();

    void SetSelected(int i);
    void ClearSelected();
    void UpdateMarkers();

    void ComputeTemplateCut(const QVector3D &n, const QVector3D &p,
                            QList<QVector3D> &cut_segments);
    void UpdateTemplateCut();
    void UpdateMagneticCutSurface(const QList<QVector3D> &plane_ns,
                                  const QList<QVector3D> &plane_ps);
    void DrawTemplateCut(const QList<QVector3D> &cut_segments);

    void DeleteTemplateImage();

    void QVector3DToArray(const QVector3D &p, double array[3]);

    // QTabWidget * sideWidget;

    QWidget *editWidget;
    QWidget *genWidget;
    QWidget *guidesWidget;
    QWidget *physicsWidget;
    QWidget *viewWidget;

    QWidget *widget;

    QTimer animate_timer;
    QDateTime animate_until;
    int animate_frames;
    int animate_mouseupdates;

    QVector2D mouse_pos;

    int grid_size;
    QVector3D default_lookat;

    // camera
    float cam_animate_duration;
    float cam_animate_angle;
    float cam_lookat_distance;
    PivotCamera cam;

    // transform widget
    TransformWidget transform_widget;
    bool transforming;

    // important quantity: metres per unit
    double metres_per_unit;

    // symmetry
    bool x_symmetry;
    bool y_symmetry;
    bool z_symmetry;

    float slab_thickness;
    float calibration_factor;
    int quality_samples;
    float deadzone_radius;

    int generate_blend_sections;
    int generate_revolve_sections;
    float generate_linear_spacing;
    bool generate_linear_scalex;
    bool generate_linear_scaley;
    float generate_slices_spacing;
    bool generate_slices_x;
    bool generate_slices_y;
    bool generate_slices_z;
    float generate_grid_sizex;
    float generate_grid_sizey;
    float generate_grid_staplesize;
    float generate_branching_scalechild;
    float generate_radial_sectors;
    float generate_radial_params[9];
    bool generate_surfacefacets_teeth;

    // toggles for various feasibility tests
    bool do_local_symmetry;
    bool do_cycles_test;
    bool do_stability_test;
    bool do_physics_test;
    bool do_connected_test;

    bool pen_mode;

    bool do_show_tnb_frames;
    bool do_show_shadow;
    bool do_show_templates;
    bool do_magnetic_cuts;
    float magnet_strength;

    double physics_new_weight_mass;
    double physics_material_density;
    double physics_max_stress;

    // gesture state
    GestureState state;
    LastOperation last_op;

    // sections - this is the main container
    PlanarSection active_section;
    PlanarSection active_section_symmetry;

    QList<PlanarSection> sections;
    float section_error_tolerance;
    float section_error_tolerance_template;

    // undo related stuff
    int max_undo_sections;  // max size of undo list
    int undo_index;         // current index within the undo list
    QList<QList<PlanarSection> > undo_sections;  // the undo list

    // geometric stability-related stuff
    QList<QVector3D> sections_convex_hull;
    QVector3D centre_of_mass;
    bool centre_of_mass_in_hull;

    // template surface stuff
    QVector<QVector3D> template_verts;
    QVector<QVector3D> template_magnetic_verts;
    QVector<QVector3D> template_norms;
    QVector<int> template_faces;
    QVector<QList<int> >
        template_poly_faces;  // these are the OBJ faces which haven't been
                              // subdivided into triangles
    QVector<int> template_facenorms;
    QList<QVector3D> template_cut;
    float template_cut_snap_distance_3d;

    // template image stuff
    QImage template_image;
    GLuint template_image_tex;
    QVector2D template_pos;
    float template_rotation;
    float template_scale;
    bool template_flipx;

    // markers for visual debugging
    QList<QVector3D> markers;
    QList<QVector3D> markers_col;
    QList<QVector3D> markers2;
    QList<QVector3D> markers_col2;
    QVector3D slot_start;
    QVector3D slot_end;

    // for recursive/procedurally defining stuff
    QVector3D recursive_slot_start;
    QVector3D recursive_slot_end;

    // for dimensioning tool
    QVector3D dimensiontool_start;
    QVector3D dimensiontool_end;

    // for skinning stuff
    ContourGraph contour_graph;

    // for physical simulation stuff
    CPhysics physics;

    // some side display widgets
    QLabel *thick_label;
    QLabel *calibration_label;
    QLabel *quality_label;
    QLabel *rotate_label;
    QLabel *rotate_angle_label;
    QLabel *linear_spacing_label;
    QLabel *slices_spacing_label;
    QLabel *num_blend_sections_label;
    QLabel *num_revolve_sections_label;
    QLabel *grid_sizex_label;
    QLabel *grid_sizey_label;
    QLabel *grid_staple_label;
    QLabel *branching_scalechild_label;
    QSlider *radial_slider[9];
    QLabel *radial_sectors_label;
    QLabel *new_weight_label;

    // selection related
    int selected;
    QList<int> last_selected;

    // translate/rotate/scale related
    QVector2D anchor_point;
    PlanarSection original_section;

    bool update_sections_disp_list;
    GLuint sections_disp_list;

    //---- Needed to connect to menus in mainWindow
    QCheckBox *deformed_checkbox;
    QCheckBox *skeleton_checkbox;
    QCheckBox *forces_checkbox;
    QCheckBox *section_checkbox;
    QCheckBox *moment_checkbox;

    QPushButton *testPhysicsButton;
    QPushButton *toggleTemplatesButton;

    QCheckBox *show_tnb_frames_checkbox;
    QCheckBox *show_cycles_test_checkbox;
    QCheckBox *show_stability_checkbox;
    QCheckBox *show_shadow_checkbox;
    QCheckBox *show_connectivity_checkbox;

    QString open_filename;

    ToolState current_tool_state;
    GenerateState gen_state;

    // Generate Selections
    QList<int> generate_selections;  // holds indices of selected sections for
                                     // the generate tool
    int num_generate_selected;       // the number selected
    int selections_per_gen_type[GENSTATE_NUM];  // the number of needed
                                                // selections for each type of
                                                // generate

    QList<PlanarSection> generated_sections;  // the temporary generate sections

    // Radial Tool
    PlanarSection pre_radial_section;
    PlanarSection radial_section;
    int radial_section_index;

    QErrorMessage errorMessage;
};

#endif  // GLWIDGET_H
