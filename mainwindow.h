#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>
#include <QWebView>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QUrl>

#include "glwidget.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    MainWindow();
    ~MainWindow();
    
protected:

    void keyPressEvent(QKeyEvent * event);
    void keyReleaseEvent(QKeyEvent * event);

private slots:   

    void ShowAppWidgets();
    void UpdateWindowTitle();

    void NewPlaneSketch();
    void LoadTemplateOBJ();
    void LoadTemplateCurve();
    void LoadTemplateImage();
    void LoadPlaneSketch();
    void SavePlaneSketch();
    void SaveSliceOBJ();
    void SaveSlabOBJ();
    void SaveSurfaceOBJ();
    void SaveSVG();
    void SaveDXF();
    void SavePhysicsOutput();
    void SaveCalibration();
    void Exit();

    void Undo();
    void Redo();
    void Delete();
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
    void ToggleUseMagneticCuts();
    void TogglePenMode();
    void ToggleLocalSymmetry();
    void StartDimensioningTool();

    void ViewIso1();
    void ViewIso2();
    void ViewIso3();
    void ViewIso4();
    void ViewX();
    void ViewY();
    void ViewZ();
    void ViewPart();
    void ToggleCyclesTest();
    void ToggleConnectedTest();
    void ToggleStabilityTest();
    void ToggleShowTNBFrames();
    void ToggleShowShadow();
    void ToggleShowTemplates();

    void SetMultisampling0();
    void SetMultisampling4();
    void SetMultisampling16();

    void GenerateBranchingSetRoot();
    void GenerateBranching();
    void GenerateLinear();
    void GenerateBlend();
    void GenerateGrid();
    void GenerateRevolve();
    void GenerateMakeCircle();
    void GenerateMakeRectangle();
    void GenerateMakeRadial();
    void GenerateMakeRadialHole();
    void GenerateSurfaceFacets();
    void GenerateSlices();

    void TogglePhysicsTest();
    void TogglePhysicsDeformed();
    void TogglePhysicsSkeleton();
    void TogglePhysicsForce();
    void TogglePhysicsSection();
    void TogglePhysicsSectionMoment();
    void PhysicsAddExternalMass();
    void PhysicsRemoveExternalMasses();

    void closeDialog();

    // New UI features

    void openEditWidget();
    void openGenerateWidget();
    void openGuidesWidget();
    void openPhysicsWidget();
    void openViewWidget();

    void setEditMenuChecks();
    void setViewMenuChecks();
    void setPhysicsMenuChecks();

private:   

    void createActions();
    void createMenus();

    void ShowWelcomePage();
    void SendTrackRequest();

    void SetMultisampling(const int i);

    virtual void resizeEvent(QResizeEvent * event);
    virtual void moveEvent(QMoveEvent * event);

    QMenu *fileMenu;
    QAction *newPlaneSketchAct;    
    QAction *loadCurveAct;
    QAction *loadImageAct;
    QAction *loadOBJAct;
    QAction *loadPlaneSketchAct;
    QAction *savePlaneSketchAct;
    QAction *saveSliceOBJAct;
    QAction *saveSlabOBJAct;
    QAction *saveSurfaceOBJAct;
    QAction *saveSVGAct;
    QAction *saveDXFAct;
    QAction *saveCalibrationAct;
    QAction *savePhysicsOutput;
    QAction *exitAct;

    QMenu *editMenu;
    QAction *useLocalSymmetryAct;
    QAction *copyMirrorXAct;
    QAction *copyRotateYAct;
    QAction *copyMirrorZAct;
    QAction *copyDuplicateAct;
    QAction *undoAct;
    QAction *redoAct;
    QAction *deleteAct;
    QAction *transformAct;
    QAction *resketchCurveAct;
    QAction *addHoleBoundaryAct;
    QAction *removeHolesBoundaryAct;
    QAction *dimensioningToolAct;
    QAction *snapMajorAxisAct;
    QAction *useMagneticCutsAct;
    QAction *createSurfacePatchesAct;
    QAction *deleteSurfacePatchesAct;

    QMenu *viewMenu;
    QAction *viewIso1Act;
    QAction *viewIso2Act;
    QAction *viewIso3Act;
    QAction *viewIso4Act;
    QAction *viewXAct;
    QAction *viewYAct;
    QAction *viewZAct;
    QAction *viewPartAct;
    QAction *viewTNBFramesAct;
    QAction *viewShadowAct;
    QAction *viewTemplatesAct;
    QAction *testCyclesAct;
    QAction *testConnectedAct;
    QAction *testStabilityAct;
    QAction *multisample0Act;
    QAction *multisample4Act;
    QAction *multisample16Act;

    QMenu *generateMenu;
    QAction *generateBranchingSetRootAct;
    QAction *generateBranchingAct;
    QAction *generateLinearAct;
    QAction *generateBlendAct;
    QAction *generateGridAct;
    QAction *generateRevolveAct;
    QAction *generateMakeCircleAct;
    QAction *generateMakeRectangleAct;
    QAction *generateMakeRadialAct;
    QAction *generateMakeRadialHoleAct;
    QAction *generateSlicesAct;
    QAction *generateSurfaceFacetsAct;

    QMenu *physicsMenu;
    QAction *testPhysicsAct;
    QAction *showPhysicsDeformedAct;
    QAction *showPhysicsSkeletonAct;
    QAction *showPhysicsForceAct;
    QAction *showPhysicsSectionAct;
    QAction *showPhysicsSectionMomentAct;
    QAction *addExternalWeightAct;
    QAction *removeExternalWeightsAct;

    GLWidget glWidget;
    //QWidget * sideWidget;
    QDockWidget * dockWidget;

    QWidget * bottomWidget;
    QDockWidget * bottomDockWidget;

    QWebView * webView;


    // New UI features

    void createSideBar();
    QToolBar *mainToolBar;
    QAction *openDock[5];
    QDockWidget *docks[5];

    void createQuickToolBar();
    QToolBar *quickToolBar;
    QDockWidget * toolWidget;

    QString window_title;
    QTimer window_title_timer;

};

#endif // MAINWINDOW_H
