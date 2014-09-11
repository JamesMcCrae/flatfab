#include "mainwindow.h"

MainWindow::MainWindow()
{   

    //title/window stuff
    setWindowTitle("FlatFab (Release 0.5 Beta - September 11, 2014)");

    //release 0.4 -
    //added "calibration shape" generation to make finding the right calibration setting a cinch
    //multisampling now a command line parameter, the default value is 4.  e.g: flatfab.exe -ms 0 would disable multisampling
    //smoother animations/frame updates
    //updated linux binary with RPATH set in the executable (rather than relying on a script)

    //release 0.3
    //updated name
    //updated program icon
    //added antialiasing/multisampling option

    //release 0.2
    // fixed rippling effect of bezier curve fit
    // made 80 degrees the default rotation angle        

    //side dock widget stuff
    sideWidget = glWidget.GetSideWidget();
    dockWidget = new QDockWidget();
    dockWidget->setWindowTitle("Settings");
    dockWidget->setWidget(sideWidget);

    //web view stuff
    webView = new QWebView(this);
    webView->setGeometry(0,0,800,400);
    webView->load(QUrl("http://flatfab.com/splash.html"));

    webView->show();

    createActions();
    createMenus();

    //set up bottom widget
    ShowWelcomePage();

    //track usage
    SendTrackRequest();
}

MainWindow::~MainWindow()
{

}

void MainWindow::ShowWelcomePage()
{

    QPushButton *button1 = new QPushButton("New FlatFab");
    QPushButton *button2 = new QPushButton("Open FlatFab...");
    QPushButton *button3 = new QPushButton("Quit");

    connect(button1, SIGNAL(clicked()), this, SLOT(NewPlaneSketch()));
    connect(button2, SIGNAL(clicked()), this, SLOT(LoadPlaneSketch()));
    connect(button3, SIGNAL(clicked()), this, SLOT(Exit()));

    connect(button1, SIGNAL(clicked()), this, SLOT(ShowAppWidgets()));
    connect(button2, SIGNAL(clicked()), this, SLOT(ShowAppWidgets()));

    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(button1);
    layout->addWidget(button2);
    layout->addWidget(button3);

    bottomWidget = new QWidget(this);
    bottomWidget->setLayout(layout);

    bottomDockWidget = new QDockWidget();
    bottomDockWidget->setFeatures(0);
    bottomDockWidget->setWindowTitle("Start using FlatFab:");
    bottomDockWidget->setWidget(bottomWidget);

    setCentralWidget(webView);
    addDockWidget(Qt::TopDockWidgetArea, bottomDockWidget);
}

void MainWindow::ShowAppWidgets()
{

    //main widget/window layout
    this->removeDockWidget(bottomDockWidget);

    addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    setCentralWidget(&glWidget);

}

void MainWindow::SendTrackRequest()
{

    //this code sends a track request to Google Anayltics for "splash.html" for the FFF account
    QNetworkAccessManager * nam = new QNetworkAccessManager(this);
    //QObject::connect(nam, SIGNAL(finished(QNetworkReply*)), this, SLOT(finishedSlot(QNetworkReply*)));

    QUrl collect_url("http://www.google-analytics.com/collect");

    QNetworkRequest request(collect_url);
    request.setHeader(QNetworkRequest::ContentTypeHeader, "application/x-www-form-urlencoded");

    QByteArray query_items("v=1&tid=UA-51900137-1&cid=555&t=pageview&dp=%2Fsplash.html");

    //QNetworkReply * reply = nam->post(request, query_items);
    nam->post(request, query_items);

}

void MainWindow::keyPressEvent(QKeyEvent *event)
{

    glWidget.keyPressEvent(event);

}

void MainWindow::keyReleaseEvent(QKeyEvent *event)
{

    glWidget.keyReleaseEvent(event);

}

void MainWindow::createActions()
{

    //file menu actions
    newPlaneSketchAct = new QAction(tr("&New FlatFab"), this);
    newPlaneSketchAct->setStatusTip(tr("Create a new FlatFab."));
    newPlaneSketchAct->setShortcuts(QKeySequence::New);
    connect(newPlaneSketchAct, SIGNAL(triggered()), this, SLOT(NewPlaneSketch()));

    loadPlaneSketchAct = new QAction(tr("&Open FlatFab..."), this);
    loadPlaneSketchAct->setStatusTip(tr("Open a FlatFab file (TXT)."));
    loadPlaneSketchAct->setShortcuts(QKeySequence::Open);
    connect(loadPlaneSketchAct, SIGNAL(triggered()), this, SLOT(LoadPlaneSketch()));

    loadCurveAct = new QAction(tr("Open &Curve Template..."), this);
    loadCurveAct->setStatusTip(tr("Open a template curve."));
    connect(loadCurveAct, SIGNAL(triggered()), this, SLOT(LoadTemplateCurve()));

    loadImageAct = new QAction(tr("Open &Image Template..."), this);
    loadImageAct->setStatusTip(tr("Open a template image."));
    connect(loadImageAct, SIGNAL(triggered()), this, SLOT(LoadTemplateImage()));

    loadOBJAct = new QAction(tr("Open Surface &Template..."), this);
    loadOBJAct->setStatusTip(tr("Open a template surface (OBJ)."));
    connect(loadOBJAct, SIGNAL(triggered()), this, SLOT(LoadTemplateOBJ()));

    savePlaneSketchAct = new QAction(tr("&Save FlatFab..."), this);
    savePlaneSketchAct->setStatusTip(tr("Save a FlatFab file (TXT)."));
    savePlaneSketchAct->setShortcuts(QKeySequence::Save);
    connect(savePlaneSketchAct, SIGNAL(triggered()), this, SLOT(SavePlaneSketch()));

    saveSliceOBJAct = new QAction(tr("Save Slice Geometry..."), this);
    saveSliceOBJAct->setStatusTip(tr("Save geometry of flat planar sections (OBJ)."));
    connect(saveSliceOBJAct, SIGNAL(triggered()), this, SLOT(SaveSliceOBJ()));

    saveSlabOBJAct = new QAction(tr("Save Slab Geometry..."), this);
    saveSlabOBJAct->setStatusTip(tr("Save geometry of thick planar sections (OBJ)."));
    connect(saveSlabOBJAct, SIGNAL(triggered()), this, SLOT(SaveSlabOBJ()));

    saveSurfaceOBJAct = new QAction(tr("Save Surface Patches..."), this);
    saveSurfaceOBJAct->setStatusTip(tr("Save geometry of created surface patches (OBJ)."));
    connect(saveSurfaceOBJAct, SIGNAL(triggered()), this, SLOT(SaveSurfaceOBJ()));

    saveSVGAct = new QAction(tr("Save SVG for Fabrication..."), this);
    saveSVGAct->setStatusTip(tr("Save fabrication template (SVG) file."));
    connect(saveSVGAct, SIGNAL(triggered()), this, SLOT(SaveSVG()));

    saveDXFAct = new QAction(tr("Save DXF for Fabrication..."), this);
    saveDXFAct->setStatusTip(tr("Save fabrication template (DXF) file."));
    connect(saveDXFAct, SIGNAL(triggered()), this, SLOT(SaveDXF()));

    saveCalibrationAct = new QAction(tr("Save Calibration Shape..."), this);
    saveCalibrationAct->setStatusTip(tr("Save calibration shape (SVG) file, used to find optimal slit width."));
    connect(saveCalibrationAct, SIGNAL(triggered()), this, SLOT(SaveCalibration()));

    savePhysicsOutput = new QAction(tr("Save PhysicsFile..."), this);
    savePhysicsOutput->setStatusTip(tr("Save model in format for external applications (TXT)."));
    connect(savePhysicsOutput, SIGNAL(triggered()), this, SLOT(SavePhysicsOutput()));

    exitAct = new QAction(tr("&Quit"), this);
    exitAct->setStatusTip(tr("Quit the application."));
    exitAct->setShortcut(QKeySequence::Quit);
    connect(exitAct, SIGNAL(triggered()), this, SLOT(Exit()));

    //edit menu actions
    useLocalSymmetryAct = new QAction(tr("&Local Symmetry"), this);
    useLocalSymmetryAct->setStatusTip(tr("Enable/disable the use of local symmetry for sketched planar sections."));
    useLocalSymmetryAct->setCheckable(true);
    useLocalSymmetryAct->setChecked(glWidget.GetDoLocalSymmetry());
    connect(useLocalSymmetryAct, SIGNAL(triggered()), this, SLOT(ToggleLocalSymmetry()));

    copyMirrorXAct = new QAction(tr("Copy (Mirror &X)"), this);
    copyMirrorXAct->setStatusTip(tr("Copy planar section, mirroring along X axis."));
    connect(copyMirrorXAct, SIGNAL(triggered()), this, SLOT(CopyMirrorX()));

    copyRotateYAct = new QAction(tr("Copy (Rotate &Y)"), this);
    copyRotateYAct->setStatusTip(tr("Copy planar section, rotation 90 degrees about Y axis."));
    connect(copyRotateYAct, SIGNAL(triggered()), this, SLOT(CopyRotateY()));

    copyMirrorZAct = new QAction(tr("Copy (Mirror &Z)"), this);
    copyMirrorZAct->setStatusTip(tr("Copy planar section, mirroring along Z axis."));
    connect(copyMirrorZAct, SIGNAL(triggered()), this, SLOT(CopyMirrorZ()));

    copyDuplicateAct = new QAction(tr("Copy (Duplicate)"), this);
    copyDuplicateAct->setStatusTip(tr("Copy planar section, duplicating it in place."));
    connect(copyDuplicateAct, SIGNAL(triggered()), this, SLOT(CopyDuplicate()));

    undoAct = new QAction(tr("&Undo"), this);
    undoAct->setStatusTip(tr("Undo the last change to the planar sections."));
    undoAct->setShortcut(QKeySequence::Undo);
    connect(undoAct, SIGNAL(triggered()), this, SLOT(Undo()));

    redoAct = new QAction(tr("&Redo"), this);
    redoAct->setStatusTip(tr("Redo the last undo to the planar sections."));
    redoAct->setShortcut(QKeySequence::Redo);
    connect(redoAct, SIGNAL(triggered()), this, SLOT(Redo()));

    deleteAct = new QAction(tr("&Delete"), this);
    deleteAct->setStatusTip(tr("Delete selected planar section."));
    deleteAct->setShortcut(QKeySequence::Delete);
    connect(deleteAct, SIGNAL(triggered()), this, SLOT(Delete()));

    transformAct = new QAction(tr("&Transform"), this);
    transformAct->setStatusTip(tr("Transform a planar section using translate/rotate/scale."));
    connect(transformAct, SIGNAL(triggered()), this, SLOT(Transform()));

    snapMajorAxisAct = new QAction(tr("Snap To Major Axis"), this);
    snapMajorAxisAct->setStatusTip(tr("Snaps the planar section normal to one of the 3 major axes."));
    connect(snapMajorAxisAct, SIGNAL(triggered()), this, SLOT(SnapToMajorAxis()));

    resketchCurveAct = new QAction(tr("Re-sketch &Boundary"), this);
    resketchCurveAct->setStatusTip(tr("Allows immediate re-sketching of an already-created planar section."));
    connect(resketchCurveAct, SIGNAL(triggered()), this, SLOT(ResketchCurve()));   

    addHoleBoundaryAct = new QAction(tr("Add &Hole"), this);
    addHoleBoundaryAct->setStatusTip(tr("Sketch a boundary for a hole into the selected planar section."));
    connect(addHoleBoundaryAct, SIGNAL(triggered()), this, SLOT(AddHoleBoundary()));

    removeHolesBoundaryAct = new QAction(tr("Remove Holes"), this);
    removeHolesBoundaryAct->setStatusTip(tr("Removes all holes from the selected planar section."));
    connect(removeHolesBoundaryAct, SIGNAL(triggered()), this, SLOT(RemoveHolesBoundary()));

    useMagneticCutsAct = new QAction(tr("Use Magnetic Cuts"), this);
    useMagneticCutsAct->setStatusTip(tr("Deform the template surface using existing planar sections."));
    useMagneticCutsAct->setCheckable(true);
    useMagneticCutsAct->setChecked(glWidget.GetDoMagneticCuts());
    connect(useMagneticCutsAct, SIGNAL(triggered()), this, SLOT(ToggleUseMagneticCuts()));

    createSurfacePatchesAct = new QAction(tr("Create Surface Patches"), this);
    createSurfacePatchesAct->setStatusTip(tr("Surface the model by creating a set of Coons patches."));
    connect(createSurfacePatchesAct, SIGNAL(triggered()), this, SLOT(CreateSurfacePatches()));

    deleteSurfacePatchesAct = new QAction(tr("Delete Surface Patches"), this);
    deleteSurfacePatchesAct->setStatusTip(tr("Delete a created set of Coons patches."));
    connect(deleteSurfacePatchesAct, SIGNAL(triggered()), this, SLOT(DeleteSurfacePatches()));

    dimensioningToolAct = new QAction(tr("Dimensioning Tool"), this);
    dimensioningToolAct->setStatusTip(tr("Scale the model by specifying dimensions along a defined line."));
    connect(dimensioningToolAct, SIGNAL(triggered()), this, SLOT(StartDimensioningTool()));

    //view menu actions
    viewIso1Act = new QAction(tr("Isometric &1"), this);
    viewIso1Act->setStatusTip(tr("Set camera to isometric view 1."));
    viewIso1Act->setShortcut(QKeySequence(tr("1")));
    connect(viewIso1Act, SIGNAL(triggered()), this, SLOT(ViewIso1()));

    viewIso2Act = new QAction(tr("Isometric &2"), this);
    viewIso2Act->setStatusTip(tr("Set camera to isometric view 2."));
    viewIso2Act->setShortcut(QKeySequence(tr("2")));
    connect(viewIso2Act, SIGNAL(triggered()), this, SLOT(ViewIso2()));

    viewIso3Act = new QAction(tr("Isometric &3"), this);
    viewIso3Act->setStatusTip(tr("Set camera to isometric view 3."));
    viewIso3Act->setShortcut(QKeySequence(tr("3")));
    connect(viewIso3Act, SIGNAL(triggered()), this, SLOT(ViewIso3()));

    viewIso4Act = new QAction(tr("Isometric &4"), this);
    viewIso4Act->setStatusTip(tr("Set camera to isometric view 4."));
    viewIso4Act->setShortcut(QKeySequence(tr("4")));
    connect(viewIso4Act, SIGNAL(triggered()), this, SLOT(ViewIso4()));

    viewXAct = new QAction(tr("Along &X"), this);
    viewXAct->setStatusTip(tr("Set camera to view along X axis."));
    viewXAct->setShortcut(QKeySequence(tr("Q")));
    connect(viewXAct, SIGNAL(triggered()), this, SLOT(ViewX()));

    viewYAct = new QAction(tr("Along &Y"), this);
    viewYAct->setStatusTip(tr("Set camera to view along Y axis."));
    viewYAct->setShortcut(QKeySequence(tr("W")));
    connect(viewYAct, SIGNAL(triggered()), this, SLOT(ViewY()));

    viewZAct = new QAction(tr("Along &Z"), this);
    viewZAct->setStatusTip(tr("Set camera to view along Z axis."));
    viewZAct->setShortcut(QKeySequence(tr("E")));
    connect(viewZAct, SIGNAL(triggered()), this, SLOT(ViewZ()));

    viewPartAct = new QAction(tr("&Selected"), this);
    viewPartAct->setStatusTip(tr("Set camera to view selected planar section."));
    viewPartAct->setShortcut(QKeySequence(tr("A")));
    connect(viewPartAct, SIGNAL(triggered()), this, SLOT(ViewPart()));

    viewTNBFramesAct = new QAction(tr("&TNB Frames"), this);
    viewTNBFramesAct->setStatusTip(tr("Show TNB frames for each planar section."));
    viewTNBFramesAct->setCheckable(true);
    viewTNBFramesAct->setChecked(glWidget.GetShowTNBFrames());
    connect(viewTNBFramesAct, SIGNAL(triggered()), this, SLOT(ToggleShowTNBFrames()));

    viewShadowAct = new QAction(tr("S&hadow"), this);
    viewShadowAct->setStatusTip(tr("Show shadow on the ground plane."));
    viewShadowAct->setCheckable(true);
    viewShadowAct->setChecked(glWidget.GetShowShadow());
    connect(viewShadowAct, SIGNAL(triggered()), this, SLOT(ToggleShowShadow()));

    viewTemplatesAct = new QAction(tr("T&emplates"), this);
    viewTemplatesAct->setStatusTip(tr("Show loaded surface or image templates."));
    viewTemplatesAct->setCheckable(true);
    viewTemplatesAct->setChecked(glWidget.GetShowTemplates());
    connect(viewTemplatesAct, SIGNAL(triggered()), this, SLOT(ToggleShowTemplates()));

    testCyclesAct = new QAction(tr("&Cycles"), this);
    testCyclesAct->setStatusTip(tr("Show cycles of the model."));
    testCyclesAct->setCheckable(true);
    testCyclesAct->setChecked(glWidget.GetDoCyclesTest());
    connect(testCyclesAct, SIGNAL(triggered()), this, SLOT(ToggleCyclesTest()));

    testConnectedAct = new QAction(tr("C&onnectivity"), this);
    testConnectedAct->setStatusTip(tr("Show connectivity of the model."));
    testConnectedAct->setCheckable(true);
    testConnectedAct->setChecked(glWidget.GetDoConnectedTest());
    connect(testConnectedAct, SIGNAL(triggered()), this, SLOT(ToggleConnectedTest()));

    testStabilityAct = new QAction(tr("&Stability"), this);
    testStabilityAct->setStatusTip(tr("Show stability of the model."));
    testStabilityAct->setCheckable(true);
    testStabilityAct->setChecked(glWidget.GetDoStabilityTest());
    connect(testStabilityAct, SIGNAL(triggered()), this, SLOT(ToggleStabilityTest()));

    multisample0Act = new QAction(tr("No Multisamping"), this);
    multisample0Act->setStatusTip(tr("Show stability of the model."));
    multisample0Act->setCheckable(true);
    multisample0Act->setChecked(glWidget.GetDoStabilityTest());
    connect(multisample0Act, SIGNAL(triggered()), this, SLOT(SetMultisampling0()));

    multisample4Act = new QAction(tr("4x Multisamping"), this);
    multisample4Act->setStatusTip(tr("Show stability of the model."));
    multisample4Act->setCheckable(true);
    multisample4Act->setChecked(glWidget.GetDoStabilityTest());
    connect(multisample4Act, SIGNAL(triggered()), this, SLOT(SetMultisampling4()));

    multisample16Act = new QAction(tr("16x Multisamping"), this);
    multisample16Act->setStatusTip(tr("Show stability of the model."));
    multisample16Act->setCheckable(true);
    multisample16Act->setChecked(glWidget.GetDoStabilityTest());
    connect(multisample16Act, SIGNAL(triggered()), this, SLOT(SetMultisampling16()));

    //generate menu actions
    generateBranchingSetRootAct = new QAction(tr("Branching (Set &Root Slit)"), this);
    generateBranchingSetRootAct->setStatusTip(tr("Set the slit for the root node to generate branching planar sections."));
    generateBranchingSetRootAct->setShortcut(QKeySequence("Ctrl+R"));
    connect(generateBranchingSetRootAct, SIGNAL(triggered()), this, SLOT(GenerateBranchingSetRoot()));

    generateBranchingAct = new QAction(tr("&Branching"), this);
    generateBranchingAct->setStatusTip(tr("Generate branching planar sections."));
    generateBranchingAct->setShortcut(QKeySequence("Ctrl+B"));
    connect(generateBranchingAct, SIGNAL(triggered()), this, SLOT(GenerateBranching()));

    generateLinearAct = new QAction(tr("&Linear"), this);
    generateLinearAct->setStatusTip(tr("Generate linearly arranged planar sections."));
    generateLinearAct->setShortcut(QKeySequence("Ctrl+L"));
    connect(generateLinearAct, SIGNAL(triggered()), this, SLOT(GenerateLinear()));

    generateBlendAct = new QAction(tr("Bl&end"), this);
    generateBlendAct->setStatusTip(tr("Generate an arrangement of blended planar sections."));
    generateBlendAct->setShortcut(QKeySequence("Ctrl+E"));
    connect(generateBlendAct, SIGNAL(triggered()), this, SLOT(GenerateBlend()));

    generateGridAct = new QAction(tr("&Grid"), this);
    generateGridAct->setStatusTip(tr("Create a grid of planar sections whose elements fit within a specified rectangular boundary."));
    generateGridAct->setShortcut(QKeySequence("Ctrl+G"));
    connect(generateGridAct, SIGNAL(triggered()), this, SLOT(GenerateGrid()));

    generateRevolveAct = new QAction(tr("Re&volve"), this);
    generateRevolveAct->setStatusTip(tr("Generate a sequence of planar sections which revolve around a base section."));
    connect(generateRevolveAct, SIGNAL(triggered()), this, SLOT(GenerateRevolve()));

    generateMakeCircleAct = new QAction(tr("Make &Circle"), this);
    generateMakeCircleAct->setStatusTip(tr("Mkae the boundary of an existing planar section a circle."));
    connect(generateMakeCircleAct, SIGNAL(triggered()), this, SLOT(GenerateMakeCircle()));

    generateMakeRectangleAct = new QAction(tr("Make Rec&tangle"), this);
    generateMakeRectangleAct->setStatusTip(tr("Mkae the boundary of an existing planar section a rectangle."));
    connect(generateMakeRectangleAct, SIGNAL(triggered()), this, SLOT(GenerateMakeRectangle()));

    generateMakeRadialAct = new QAction(tr("Make Radial"), this);
    generateMakeRadialAct->setStatusTip(tr("Make the boundary of an existing planar section a radial pattern."));
    connect(generateMakeRadialAct, SIGNAL(triggered()), this, SLOT(GenerateMakeRadial()));

    generateMakeRadialHoleAct = new QAction(tr("Make Radial Holes"), this);
    generateMakeRadialHoleAct->setStatusTip(tr("Make a radial set of holes in a selected planar section."));
    connect(generateMakeRadialHoleAct, SIGNAL(triggered()), this, SLOT(GenerateMakeRadialHole()));

    generateSurfaceFacetsAct = new QAction(tr("Surface &Facets"), this);
    generateSurfaceFacetsAct->setStatusTip(tr("Create sections from the planar facets of a template surface."));
    connect(generateSurfaceFacetsAct, SIGNAL(triggered()), this, SLOT(GenerateSurfaceFacets()));

    //physics menu actions
    testPhysicsAct = new QAction(tr("Test &Physics"), this);
    testPhysicsAct->setStatusTip(tr("Perform physically-based tests of the model."));
    testPhysicsAct->setCheckable(true);
    testPhysicsAct->setChecked(glWidget.GetDoPhysicsTest());
    connect(testPhysicsAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsTest()));

    CPhysics & physics = glWidget.GetPhysics();

    showPhysicsDeformedAct = new QAction(tr("Show &Deformed"), this);
    showPhysicsDeformedAct->setStatusTip(tr("Show deformed physics."));
    showPhysicsDeformedAct->setCheckable(true);
    showPhysicsDeformedAct->setChecked(physics.GetDrawDeformed());
    connect(showPhysicsDeformedAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsDeformed()));

    showPhysicsSkeletonAct = new QAction(tr("Show &Skeleton"), this);
    showPhysicsSkeletonAct->setStatusTip(tr("Show skeleton used in physics."));
    showPhysicsSkeletonAct->setCheckable(true);
    showPhysicsSkeletonAct->setChecked(physics.GetDrawSkeleton());
    connect(showPhysicsSkeletonAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsSkeleton()));

    showPhysicsForceAct = new QAction(tr("Show &Force"), this);
    showPhysicsForceAct->setStatusTip(tr("Show calcualted forces."));
    showPhysicsForceAct->setCheckable(true);
    showPhysicsForceAct->setChecked(physics.GetDrawForce());
    connect(showPhysicsForceAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsForce()));

    showPhysicsSectionAct = new QAction(tr("Show S&ections"), this);
    showPhysicsSectionAct->setStatusTip(tr("Show planar sections used by physics."));
    showPhysicsSectionAct->setCheckable(true);
    showPhysicsSectionAct->setChecked(physics.GetDrawSection());
    connect(showPhysicsSectionAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsSection()));

    showPhysicsSectionMomentAct = new QAction(tr("Show Section &Moments"), this);
    showPhysicsSectionMomentAct->setStatusTip(tr("Show within-plane forces calculated by physics."));
    showPhysicsSectionMomentAct->setCheckable(true);
    showPhysicsSectionMomentAct->setChecked(physics.GetDrawSectionMoment());
    connect(showPhysicsSectionMomentAct, SIGNAL(triggered()), this, SLOT(TogglePhysicsSectionMoment()));

    addExternalWeightAct = new QAction(tr("Add &Weight"), this);
    addExternalWeightAct->setStatusTip(tr("Add an external mass to the selected plane."));
    connect(addExternalWeightAct, SIGNAL(triggered()), this, SLOT(PhysicsAddExternalMass()));

    removeExternalWeightsAct = new QAction(tr("Remove Weights"), this);
    removeExternalWeightsAct->setStatusTip(tr("Removes external masses from the selected plane."));
    connect(removeExternalWeightsAct, SIGNAL(triggered()), this, SLOT(PhysicsRemoveExternalMasses()));

}

void MainWindow::createMenus()
{

    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(newPlaneSketchAct);
    fileMenu->addSeparator();
    fileMenu->addAction(loadPlaneSketchAct);
    fileMenu->addAction(loadCurveAct);
    fileMenu->addAction(loadImageAct);
    fileMenu->addAction(loadOBJAct);    
    fileMenu->addSeparator();
    fileMenu->addAction(savePlaneSketchAct);
    fileMenu->addAction(saveSliceOBJAct);
    fileMenu->addAction(saveSlabOBJAct);
    //fileMenu->addAction(saveSurfaceOBJAct);
    fileMenu->addAction(saveSVGAct);
    fileMenu->addAction(saveDXFAct);
    fileMenu->addAction(saveCalibrationAct);
    //fileMenu->addAction(savePhysicsOutput);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    editMenu = menuBar()->addMenu(tr("&Edit"));
    editMenu->addAction(undoAct);
    editMenu->addAction(redoAct);
    editMenu->addAction(deleteAct);
    editMenu->addSeparator();
    editMenu->addAction(transformAct);
    editMenu->addSeparator();
    editMenu->addAction(copyMirrorXAct);
    editMenu->addAction(copyRotateYAct);
    editMenu->addAction(copyMirrorZAct);
    editMenu->addAction(copyDuplicateAct);
    //editMenu->addSeparator();
    editMenu->addAction(createSurfacePatchesAct);
    editMenu->addAction(deleteSurfacePatchesAct);
    editMenu->addSeparator();
    editMenu->addAction(addHoleBoundaryAct);
    editMenu->addAction(removeHolesBoundaryAct);
    editMenu->addSeparator();
    editMenu->addAction(dimensioningToolAct);
    editMenu->addSeparator();
    editMenu->addAction(snapMajorAxisAct);
    editMenu->addAction(resketchCurveAct);
    editMenu->addAction(useMagneticCutsAct);
    editMenu->addAction(useLocalSymmetryAct);    

    viewMenu = menuBar()->addMenu(tr("&View"));
    viewMenu->addAction(viewIso1Act);
    viewMenu->addAction(viewIso2Act);
    viewMenu->addAction(viewIso3Act);
    viewMenu->addAction(viewIso4Act);
    viewMenu->addAction(viewXAct);
    viewMenu->addAction(viewYAct);
    viewMenu->addAction(viewZAct);
    viewMenu->addAction(viewPartAct);
    viewMenu->addSeparator();
    viewMenu->addAction(viewTNBFramesAct);
    viewMenu->addAction(viewShadowAct);
    viewMenu->addSeparator();
    viewMenu->addAction(viewTemplatesAct);
    viewMenu->addSeparator();
    viewMenu->addAction(testCyclesAct);
    viewMenu->addAction(testConnectedAct);
    viewMenu->addAction(testStabilityAct);
    //viewMenu->addAction(multisample0Act);
    //viewMenu->addAction(multisample4Act);
    //viewMenu->addAction(multisample16Act);

    generateMenu = menuBar()->addMenu(tr("&Generate"));        
    generateMenu->addAction(generateBlendAct);
    generateMenu->addAction(generateBranchingSetRootAct);
    generateMenu->addAction(generateBranchingAct);
    generateMenu->addAction(generateGridAct);    
    generateMenu->addAction(generateLinearAct);
    generateMenu->addAction(generateMakeCircleAct);
    generateMenu->addAction(generateMakeRadialAct);
    generateMenu->addAction(generateMakeRadialHoleAct);
    generateMenu->addAction(generateMakeRectangleAct);    
    generateMenu->addAction(generateRevolveAct);
    generateMenu->addAction(generateSurfaceFacetsAct);

    physicsMenu = menuBar()->addMenu(tr("&Physics"));
    physicsMenu->addAction(testPhysicsAct);
    physicsMenu->addSeparator();
    physicsMenu->addAction(showPhysicsDeformedAct);
    physicsMenu->addAction(showPhysicsSkeletonAct);
    physicsMenu->addAction(showPhysicsForceAct);
    physicsMenu->addAction(showPhysicsSectionAct);
    physicsMenu->addAction(showPhysicsSectionMomentAct);
    physicsMenu->addSeparator();
    physicsMenu->addAction(addExternalWeightAct);
    physicsMenu->addAction(removeExternalWeightsAct);

}

void MainWindow::NewPlaneSketch()
{
    glWidget.NewPlaneSketch();
}

void MainWindow::LoadTemplateOBJ()
{
    glWidget.LoadTemplateOBJ();
}

void MainWindow::LoadTemplateCurve()
{
    glWidget.LoadTemplateCurve();
}

void MainWindow::LoadTemplateImage()
{
    glWidget.LoadTemplateImage();
}

void MainWindow::LoadPlaneSketch()
{
    glWidget.LoadPlaneSketch();
}

void MainWindow::SavePlaneSketch()
{
    glWidget.SavePlaneSketch();
}

void MainWindow::SaveSliceOBJ()
{
    glWidget.SaveSliceOBJ();
}

void MainWindow::SaveSlabOBJ()
{
    glWidget.SaveSlabOBJ();
}

void MainWindow::SaveSurfaceOBJ()
{
    glWidget.SaveSurfaceOBJ();
}

void MainWindow::SaveSVG()
{
    glWidget.SaveSVG();
}

void MainWindow::SaveDXF()
{
    glWidget.SaveDXF();
}

void MainWindow::SavePhysicsOutput()
{
    glWidget.SavePhysicsOutput();
}

void MainWindow::SaveCalibration()
{
    glWidget.SaveCalibration();
}

void MainWindow::Exit()
{
    close();
    qApp->quit();
}

void MainWindow::ToggleUseMagneticCuts()
{
    glWidget.SetDoMagneticCuts(!glWidget.GetDoMagneticCuts());
}

void MainWindow::ToggleLocalSymmetry()
{
    glWidget.SetDoLocalSymmetry(!glWidget.GetDoLocalSymmetry());
}

void MainWindow::StartDimensioningTool()
{
    glWidget.StartDimensioningTool();
}

void MainWindow::Undo()
{
    glWidget.Undo(OP_UNDO);
    //glWidget.updateGL();
}

void MainWindow::Redo()
{
    glWidget.Redo();
}

void MainWindow::Delete()
{
    glWidget.DeleteSelected();
}

void MainWindow::Transform()
{
    glWidget.Transform();
}

void MainWindow::CopyMirrorX()
{
    glWidget.CopyMirrorX();
}

void MainWindow::CopyRotateY()
{
    glWidget.CopyRotateY();
}

void MainWindow::CopyMirrorZ()
{
    glWidget.CopyMirrorZ();
}

void MainWindow::CopyDuplicate()
{
    glWidget.CopyDuplicate();
}

void MainWindow::SnapToMajorAxis()
{
    glWidget.SnapToMajorAxis();
}

void MainWindow::AddHoleBoundary()
{
    glWidget.AddHoleBoundary();
}

void MainWindow::RemoveHolesBoundary()
{
    glWidget.RemoveHolesBoundary();
}

void MainWindow::ResketchCurve()
{
    glWidget.ResketchCurve();
}

void MainWindow::CreateSurfacePatches()
{
    glWidget.CreateSurfacePatches();
}

void MainWindow::DeleteSurfacePatches()
{
    glWidget.DeleteSurfacePatches();
}

void MainWindow::ViewIso1()
{
    glWidget.SetViewIso1();
}

void MainWindow::ViewIso2()
{
    glWidget.SetViewIso2();
}

void MainWindow::ViewIso3()
{
    glWidget.SetViewIso3();
}

void MainWindow::ViewIso4()
{
    glWidget.SetViewIso4();
}

void MainWindow::ViewX()
{
    glWidget.SetViewX();
}

void MainWindow::ViewY()
{
    glWidget.SetViewY();
}

void MainWindow::ViewZ()
{
    glWidget.SetViewZ();
}

void MainWindow::ViewPart()
{
    glWidget.SetViewPart();
}

void MainWindow::ToggleCyclesTest()
{
    glWidget.SetDoCyclesTest(!glWidget.GetDoCyclesTest());
}

void MainWindow::ToggleConnectedTest()
{
    glWidget.SetDoConnectedTest(!glWidget.GetDoConnectedTest());
}

void MainWindow::ToggleStabilityTest()
{
    glWidget.SetDoStabilityTest(!glWidget.GetDoStabilityTest());
}

void MainWindow::ToggleShowTNBFrames()
{
    glWidget.SetShowTNBFrames(!glWidget.GetShowTNBFrames());
}

void MainWindow::ToggleShowShadow()
{
    glWidget.SetShowShadow(!glWidget.GetShowShadow());
}

void MainWindow::ToggleShowTemplates()
{
    glWidget.SetShowTemplates(!glWidget.GetShowTemplates());
}

void MainWindow::GenerateBranchingSetRoot()
{
    glWidget.DoGenerateBranchingSetRoot();
}

void MainWindow::GenerateBranching()
{
    glWidget.DoGenerateBranching();
}

void MainWindow::GenerateLinear()
{
    glWidget.DoGenerateLinear();
}

void MainWindow::GenerateBlend()
{
    glWidget.DoGenerateBlend();
}

void MainWindow::GenerateGrid()
{
    glWidget.DoGenerateGrid();
}

void MainWindow::GenerateRevolve()
{
    glWidget.DoGenerateRevolve();
}

void MainWindow::GenerateMakeCircle()
{
    glWidget.DoGenerateMakeCircle();
}

void MainWindow::GenerateMakeRectangle()
{
    glWidget.DoGenerateMakeRectangle();
}

void MainWindow::GenerateMakeRadial()
{
    glWidget.DoGenerateMakeRadial();
}

void MainWindow::GenerateMakeRadialHole()
{
    glWidget.DoGenerateMakeRadialHole();
}

void MainWindow::GenerateSurfaceFacets()
{
    glWidget.DoGenerateSurfaceFacets();
}

void MainWindow::TogglePhysicsTest()
{
    glWidget.SetDoPhysicsTest(!glWidget.GetDoPhysicsTest());
}

void MainWindow::TogglePhysicsDeformed()
{
    CPhysics & physics = glWidget.GetPhysics();
    physics.SetDrawDeformed(!physics.GetDrawDeformed());
    //glWidget.updateGL();
}

void MainWindow::TogglePhysicsSkeleton()
{
    CPhysics & physics = glWidget.GetPhysics();
    physics.SetDrawSkeleton(!physics.GetDrawSkeleton());
    //glWidget.updateGL();
}

void MainWindow::TogglePhysicsForce()
{
    CPhysics & physics = glWidget.GetPhysics();
    physics.SetDrawForce(!physics.GetDrawForce());
    //glWidget.updateGL();
}

void MainWindow::TogglePhysicsSection()
{
    CPhysics & physics = glWidget.GetPhysics();
    physics.SetDrawSection(!physics.GetDrawSection());
    //glWidget.updateGL();
}

void MainWindow::TogglePhysicsSectionMoment()
{
    CPhysics & physics = glWidget.GetPhysics();
    physics.SetDrawSectionMoment(!physics.GetDrawSectionMoment());
    //glWidget.updateGL();
}

void MainWindow::PhysicsAddExternalMass()
{
    glWidget.DoPhysicsAddWeight();
    //glWidget.updateGL();
}

void MainWindow::PhysicsRemoveExternalMasses()
{
    glWidget.DoPhysicsRemoveWeights();
    //glWidget.updateGL();
}

void MainWindow::SetMultisampling0()
{
    SetMultisampling(0);
}

void MainWindow::SetMultisampling4()
{
    SetMultisampling(4);
}

void MainWindow::SetMultisampling16()
{
    SetMultisampling(16);
}

void MainWindow::SetMultisampling(const int i)
{
    //qDebug() << "MainWindow::SetMultisampling(const int i)" << i;
    // Setting up Multi-sampling
    QGLFormat glf = QGLFormat::defaultFormat();
    glf.setSampleBuffers(i > 0);
    glf.setSamples(i);
    QGLFormat::setDefaultFormat(glf);
}
