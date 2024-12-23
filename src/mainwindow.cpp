#include "mainwindow.h"
#include "version.h"

MainWindow::MainWindow()
{
    window_title = QString("FlatFab ") + PROJECT_VERSION;
    UpdateWindowTitle();
    createActions();
    createMenus();
    createSideBar();
    createQuickTools();
    setCentralWidget(&glWidget);
}

MainWindow::~MainWindow() {}

void MainWindow::ShowWelcomeMessageBox()
{
    QMessageBox mb;
    mb.setPalette(QPalette(QColor(230, 230, 230), QColor(255, 255, 255)));

    mb.setMinimumSize(500, 300);
    mb.setTextFormat(Qt::RichText);
    mb.setWindowTitle("Welcome to FlatFab");
    mb.setText(
        "<p><font size='6'><b>Getting Started</b></font></p>"

        "<p>Start by drawing your first shape with the left mouse button "
        "held</p>"

        "<p><font size='4'><b>Create a New Section</b></font>"
        "<table> <tr><td>- Left click on a section and drag to the section "
        "edge to form a slit<br>"
        "- Now draw the new section while holding the left mouse "
        "button</td></tr></table></p>"

        "<p><font size='4'><b>Selecting a Section</b></font>"
        "<table> <tr><td>- Sections must be selected to perform editing "
        "operation<br>"
        "- Click on a section with the right mouse button to select "
        "it</td></tr></table></p>"

        "<p><font size='4'><b>Delete a Section</b></font>"
        "<table> <tr><td>- Press <b>Delete</b>/<b>Backspace</b> to delete the "
        "selected section</td></tr></p>"

        "<p><font size='4'><b>Editing Section Control Points</b></font>"
        "<table> <tr> <td>- Right click and drag the control points to modify "
        "the section shape (hold ctrl to create sharp corners)<br>"
        "- While hovering over a control point press <b>-</b> to delete a "
        "handle or <b>+</b> to add one</td></tr></table></p>"

        "<p><font size='4'><b>Camera Controls</b></font>"
        "<table> <tr><td><b>orbit</b><br><b>zoom</b><br><b>dolly</b></td>"
        "<td>left-click (anywhere off the model)<br>"
        "middle-click + ctrl (or mouse wheel)<br>"
        "middle-click</td> </table>");

    mb.setStandardButtons(QMessageBox::Ok);
    mb.exec();
}

void MainWindow::createQuickTools()
{    
    // Local symmetry button
    QPushButton *symButton = new QPushButton(
        QIcon(":/icons/appbar.transform.flip.horizontal.png"), "", this);
    symButton->setCheckable(true);
    symButton->setChecked(false);
    symButton->setToolTip("Local Symmetry");
    symButton->setStyleSheet(
        "QPushButton {icon-size:35px; max-width:35px; max-height:35px;}"
        "QPushButton:closed {background-color: #22c024;}"
        "QPushButton:open {background-color: #fff;}");
    connect(symButton, SIGNAL(clicked()), this, SLOT(ToggleLocalSymmetry()));

    // Pen method
    QPushButton *penButton =
        new QPushButton(QIcon(":/icons/appbar.vector.pen.png"), "", this);
    penButton->setCheckable(true);
    penButton->setChecked(true);
    penButton->setToolTip("Pen Mode");
    penButton->setStyleSheet(
        "QPushButton {icon-size:35px; max-width:35px; max-height:35px;}"
        "QPushButton:closed {background-color: #22c024;}"
        "QPushButton:open {background-color: #fff;}");
    connect(penButton, SIGNAL(clicked()), this, SLOT(TogglePenMode()));

    // Setting up layout and widget
    QHBoxLayout *layout = new QHBoxLayout();
    layout->addWidget(symButton);
    layout->addWidget(penButton);
    QWidget *widget = new QWidget();
    widget->setLayout(layout);

    mainToolBar->addWidget(widget);
}

void MainWindow::createSideBar()
{
    QIcon icon[5];
    icon[0].addFile(":/icons/edit-normal.png", QSize(), QIcon::Normal,
                    QIcon::Off);
    icon[0].addFile(":/icons/edit-hover.png", QSize(), QIcon::Active);
    icon[0].addFile(":/icons/edit-checked.png", QSize(), QIcon::Normal,
                    QIcon::On);

    openDock[0] = new QAction(icon[0], tr("&Edit"), this);
    openDock[0]->setStatusTip(tr("Edit"));
    openDock[0]->setCheckable(true);
    connect(openDock[0], SIGNAL(triggered()), this, SLOT(openEditWidget()));

    icon[1].addFile(":/icons/generate-normal.png", QSize(), QIcon::Normal,
                    QIcon::Off);
    icon[1].addFile(":/icons/generate-hover.png", QSize(), QIcon::Active);
    icon[1].addFile(":/icons/generate-checked.png", QSize(), QIcon::Normal,
                    QIcon::On);

    openDock[1] = new QAction(icon[1], tr("&Generate"), this);
    openDock[1]->setStatusTip(tr("Generate"));
    openDock[1]->setCheckable(true);
    connect(openDock[1], SIGNAL(triggered()), this, SLOT(openGenerateWidget()));

    icon[2].addFile(":/icons/guides-normal.png", QSize(), QIcon::Normal,
                    QIcon::Off);
    icon[2].addFile(":/icons/guides-hover.png", QSize(), QIcon::Active);
    icon[2].addFile(":/icons/guides-checked.png", QSize(), QIcon::Normal,
                    QIcon::On);

    openDock[2] = new QAction(icon[2], tr("&Guides"), this);
    openDock[2]->setStatusTip(tr("Guides"));
    openDock[2]->setCheckable(true);
    connect(openDock[2], SIGNAL(triggered()), this, SLOT(openGuidesWidget()));

    icon[3].addFile(":/icons/physics-normal.png", QSize(), QIcon::Normal,
                    QIcon::Off);
    icon[3].addFile(":/icons/physics-hover.png", QSize(), QIcon::Active);
    icon[3].addFile(":/icons/physics-checked.png", QSize(), QIcon::Normal,
                    QIcon::On);

    openDock[3] = new QAction(icon[3], tr("&Physics"), this);
    openDock[3]->setStatusTip(tr("Physics"));
    openDock[3]->setCheckable(true);
    connect(openDock[3], SIGNAL(triggered()), this, SLOT(openPhysicsWidget()));

    icon[4].addFile(":/icons/views-normal.png", QSize(), QIcon::Normal,
                    QIcon::Off);
    icon[4].addFile(":/icons/views-hover.png", QSize(), QIcon::Active);
    icon[4].addFile(":/icons/views-checked.png", QSize(), QIcon::Normal,
                    QIcon::On);

    openDock[4] = new QAction(icon[4], tr("&View"), this);
    openDock[4]->setStatusTip(tr("View"));
    openDock[4]->setCheckable(true);
    connect(openDock[4], SIGNAL(triggered()), this, SLOT(openViewWidget()));

    mainToolBar = new QToolBar(tr("Tools"));
    for (int i = 0; i < 4; i++) {
        mainToolBar->addAction(openDock[i]);
        mainToolBar->widgetForAction(openDock[i])
            ->setMinimumSize(QSize(75, 75));
        mainToolBar->addSeparator();
    }
    mainToolBar->addAction(openDock[4]);
    mainToolBar->widgetForAction(openDock[4])->setMinimumSize(QSize(75, 75));

    mainToolBar->setMovable(false);
    mainToolBar->setContextMenuPolicy(Qt::PreventContextMenu);
    mainToolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
    mainToolBar->setIconSize(QSize(76, 76));

    addToolBar(Qt::LeftToolBarArea, mainToolBar);

    QWidget *widgets[5];
    widgets[0] = glWidget.GetEditWidget();
    widgets[1] = glWidget.GetGenerateWidget();
    widgets[2] = glWidget.GetGuidesWidget();
    widgets[3] = glWidget.GetPhysicsWidget();
    widgets[4] = glWidget.GetViewWidget();

    for (int i = 0; i < 5; i++) {
        docks[i] = new QDockWidget(this);

        if (i == 1) {
            QScrollArea *scrollArea = new QScrollArea;
            scrollArea->setMinimumWidth(245);
            scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
            widgets[i]->layout()->setSizeConstraint(QLayout::SetFixedSize);
            scrollArea->setWidget(widgets[i]);
            docks[i]->setWidget(scrollArea);
        } else {
            docks[i]->setWidget(widgets[i]);
        }

        docks[i]->setFeatures(QDockWidget::NoDockWidgetFeatures);

        addDockWidget(Qt::LeftDockWidgetArea, docks[i]);
        docks[i]->setVisible(false);
    }
    docks[0]->setWindowTitle("Edit Selected Section");
    docks[1]->setWindowTitle("Generate Sections");
    docks[2]->setWindowTitle("Guides and Dimensions");
    docks[3]->setWindowTitle("Physical Simulation");
    docks[4]->setWindowTitle("View");
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
    // file menu actions
    newPlaneSketchAct = new QAction(tr("&New FlatFab"), this);
    newPlaneSketchAct->setStatusTip(tr("Create a new FlatFab."));
    newPlaneSketchAct->setShortcuts(QKeySequence::New);
    connect(newPlaneSketchAct, SIGNAL(triggered()), this,
            SLOT(NewPlaneSketch()));

    loadPlaneSketchAct = new QAction(tr("&Open FlatFab..."), this);
    loadPlaneSketchAct->setStatusTip(tr("Open a FlatFab file (TXT)."));
    loadPlaneSketchAct->setShortcuts(QKeySequence::Open);
    connect(loadPlaneSketchAct, SIGNAL(triggered()), this,
            SLOT(LoadPlaneSketch()));

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
    connect(savePlaneSketchAct, SIGNAL(triggered()), this,
            SLOT(SavePlaneSketch()));

    saveSliceOBJAct = new QAction(tr("Save Slice Geometry..."), this);
    saveSliceOBJAct->setStatusTip(
        tr("Save geometry of flat planar sections (OBJ)."));
    connect(saveSliceOBJAct, SIGNAL(triggered()), this, SLOT(SaveSliceOBJ()));

    saveSlabOBJAct = new QAction(tr("Save Slab Geometry..."), this);
    saveSlabOBJAct->setStatusTip(
        tr("Save geometry of thick planar sections (OBJ)."));
    connect(saveSlabOBJAct, SIGNAL(triggered()), this, SLOT(SaveSlabOBJ()));

    saveFlattenedSlabOBJAct =
        new QAction(tr("Save Flattened Slab Geometry..."), this);
    saveFlattenedSlabOBJAct->setStatusTip(
        tr("Save geometry of thick planar sections onto a 2D plane (OBJ)."));
    connect(saveFlattenedSlabOBJAct, SIGNAL(triggered()), this,
            SLOT(SaveFlattenedSlabOBJ()));

    saveSurfaceOBJAct = new QAction(tr("Save Surface Patches..."), this);
    saveSurfaceOBJAct->setStatusTip(
        tr("Save geometry of created surface patches (OBJ)."));
    connect(saveSurfaceOBJAct, SIGNAL(triggered()), this,
            SLOT(SaveSurfaceOBJ()));

    saveSVGAct = new QAction(tr("Save SVG for Fabrication..."), this);
    saveSVGAct->setStatusTip(tr("Save fabrication template (SVG) file."));
    connect(saveSVGAct, SIGNAL(triggered()), this, SLOT(SaveSVG()));

    saveDXFAct = new QAction(tr("Save DXF for Fabrication..."), this);
    saveDXFAct->setStatusTip(tr("Save fabrication template (DXF) file."));
    connect(saveDXFAct, SIGNAL(triggered()), this, SLOT(SaveDXF()));

    saveCalibrationAct = new QAction(tr("Save Calibration Shape..."), this);
    saveCalibrationAct->setStatusTip(tr(
        "Save calibration shape (SVG) file, used to find optimal slit width."));
    connect(saveCalibrationAct, SIGNAL(triggered()), this,
            SLOT(SaveCalibration()));

    savePhysicsOutput = new QAction(tr("Save PhysicsFile..."), this);
    savePhysicsOutput->setStatusTip(
        tr("Save model in format for external applications (TXT)."));
    connect(savePhysicsOutput, SIGNAL(triggered()), this,
            SLOT(SavePhysicsOutput()));

    exitAct = new QAction(tr("&Quit"), this);
    exitAct->setStatusTip(tr("Quit the application."));
    exitAct->setShortcut(QKeySequence::Quit);
    connect(exitAct, SIGNAL(triggered()), this, SLOT(Exit()));

    // edit menu actions
    useLocalSymmetryAct = new QAction(tr("&Local Symmetry"), this);
    useLocalSymmetryAct->setStatusTip(
        tr("Enable/disable the use of local symmetry for sketched planar "
           "sections."));
    useLocalSymmetryAct->setCheckable(true);
    useLocalSymmetryAct->setChecked(glWidget.GetDoLocalSymmetry());
    connect(useLocalSymmetryAct, SIGNAL(triggered()), this,
            SLOT(ToggleLocalSymmetry()));

    clipToGroundPlaneAct = new QAction(tr("Clip to &Ground Plane"), this);
    clipToGroundPlaneAct->setStatusTip(
        tr("Enable/disable clipping new planar sections to the ground "
           "plane."));
    clipToGroundPlaneAct->setCheckable(true);
    clipToGroundPlaneAct->setChecked(glWidget.GetClipToGroundPlane());
    connect(clipToGroundPlaneAct, SIGNAL(triggered()), this,
            SLOT(ToggleClipToGroundPlane()));

    copyMirrorXAct = new QAction(tr("Copy (Mirror &X)"), this);
    copyMirrorXAct->setStatusTip(
        tr("Copy planar section, mirroring along X axis."));
    connect(copyMirrorXAct, SIGNAL(triggered()), this, SLOT(CopyMirrorX()));

    copyRotateYAct = new QAction(tr("Copy (Rotate &Y)"), this);
    copyRotateYAct->setStatusTip(
        tr("Copy planar section, rotation 90 degrees about Y axis."));
    connect(copyRotateYAct, SIGNAL(triggered()), this, SLOT(CopyRotateY()));

    copyMirrorZAct = new QAction(tr("Copy (Mirror &Z)"), this);
    copyMirrorZAct->setStatusTip(
        tr("Copy planar section, mirroring along Z axis."));
    connect(copyMirrorZAct, SIGNAL(triggered()), this, SLOT(CopyMirrorZ()));

    copyDuplicateAct = new QAction(tr("Copy (Duplicate)"), this);
    copyDuplicateAct->setStatusTip(
        tr("Copy planar section, duplicating it in place."));
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
    QList<QKeySequence> del_shortcuts;
    del_shortcuts.push_back(Qt::Key_Backspace);
    del_shortcuts.push_back(Qt::Key_Delete);
    deleteAct->setShortcuts(del_shortcuts);
    connect(deleteAct, SIGNAL(triggered()), this, SLOT(Delete()));
    addAction(deleteAct);

    transformAct = new QAction(tr("&Transform"), this);
    transformAct->setStatusTip(
        tr("Transform a planar section using translate/rotate/scale."));
    connect(transformAct, SIGNAL(triggered()), this, SLOT(Transform()));

    snapMajorAxisAct = new QAction(tr("Snap To Major Axis"), this);
    snapMajorAxisAct->setStatusTip(
        tr("Snaps the planar section normal to one of the 3 major axes."));
    connect(snapMajorAxisAct, SIGNAL(triggered()), this,
            SLOT(SnapToMajorAxis()));

    resketchCurveAct = new QAction(tr("Re-sketch &Boundary"), this);
    resketchCurveAct->setStatusTip(tr(
        "Allows immediate re-sketching of an already-created planar section."));
    connect(resketchCurveAct, SIGNAL(triggered()), this, SLOT(ResketchCurve()));

    addHoleBoundaryAct = new QAction(tr("Add &Hole"), this);
    addHoleBoundaryAct->setStatusTip(
        tr("Sketch a boundary for a hole into the selected planar section."));
    connect(addHoleBoundaryAct, SIGNAL(triggered()), this,
            SLOT(AddHoleBoundary()));

    removeHolesBoundaryAct = new QAction(tr("Remove Holes"), this);
    removeHolesBoundaryAct->setStatusTip(
        tr("Removes all holes from the selected planar section."));
    connect(removeHolesBoundaryAct, SIGNAL(triggered()), this,
            SLOT(RemoveHolesBoundary()));

    useMagneticCutsAct = new QAction(tr("Use Magnetic Cuts"), this);
    useMagneticCutsAct->setStatusTip(
        tr("Deform the template surface using existing planar sections."));
    useMagneticCutsAct->setCheckable(true);
    useMagneticCutsAct->setChecked(glWidget.GetDoMagneticCuts());
    connect(useMagneticCutsAct, SIGNAL(triggered()), this,
            SLOT(ToggleUseMagneticCuts()));

    createSurfacePatchesAct = new QAction(tr("Create Surface Patches"), this);
    createSurfacePatchesAct->setStatusTip(
        tr("Surface the model by creating a set of Coons patches."));
    connect(createSurfacePatchesAct, SIGNAL(triggered()), this,
            SLOT(CreateSurfacePatches()));

    deleteSurfacePatchesAct = new QAction(tr("Delete Surface Patches"), this);
    deleteSurfacePatchesAct->setStatusTip(
        tr("Delete a created set of Coons patches."));
    connect(deleteSurfacePatchesAct, SIGNAL(triggered()), this,
            SLOT(DeleteSurfacePatches()));

    dimensioningToolAct = new QAction(tr("Dimensioning Tool"), this);
    dimensioningToolAct->setStatusTip(
        tr("Scale the model by specifying dimensions along a defined line."));
    connect(dimensioningToolAct, SIGNAL(triggered()), this,
            SLOT(StartDimensioningTool()));

    // view menu actions
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
    viewPartAct->setStatusTip(
        tr("Set camera to view selected planar section."));
    viewPartAct->setShortcut(QKeySequence(tr("A")));
    connect(viewPartAct, SIGNAL(triggered()), this, SLOT(ViewPart()));

    viewTNBFramesAct = new QAction(tr("&TNB Frames"), this);
    viewTNBFramesAct->setStatusTip(
        tr("Show TNB frames for each planar section."));
    viewTNBFramesAct->setCheckable(true);
    viewTNBFramesAct->setChecked(glWidget.GetShowTNBFrames());
    connect(viewTNBFramesAct, SIGNAL(triggered()), this,
            SLOT(ToggleShowTNBFrames()));

    viewShadowAct = new QAction(tr("S&hadow"), this);
    viewShadowAct->setStatusTip(tr("Show shadow on the ground plane."));
    viewShadowAct->setCheckable(true);
    viewShadowAct->setChecked(glWidget.GetShowShadow());
    connect(viewShadowAct, SIGNAL(triggered()), this, SLOT(ToggleShowShadow()));

    viewTemplatesAct = new QAction(tr("T&emplates"), this);
    viewTemplatesAct->setStatusTip(
        tr("Show loaded surface or image templates."));
    viewTemplatesAct->setCheckable(true);
    viewTemplatesAct->setChecked(glWidget.GetShowTemplates());
    connect(viewTemplatesAct, SIGNAL(triggered()), this,
            SLOT(ToggleShowTemplates()));

    testCyclesAct = new QAction(tr("&Cycles"), this);
    testCyclesAct->setStatusTip(tr("Show cycles of the model."));
    testCyclesAct->setCheckable(true);
    testCyclesAct->setChecked(glWidget.GetDoCyclesTest());
    connect(testCyclesAct, SIGNAL(triggered()), this, SLOT(ToggleCyclesTest()));

    testConnectedAct = new QAction(tr("C&onnectivity"), this);
    testConnectedAct->setStatusTip(tr("Show connectivity of the model."));
    testConnectedAct->setCheckable(true);
    testConnectedAct->setChecked(glWidget.GetDoConnectedTest());
    connect(testConnectedAct, SIGNAL(triggered()), this,
            SLOT(ToggleConnectedTest()));

    testStabilityAct = new QAction(tr("&Stability"), this);
    testStabilityAct->setStatusTip(tr("Show stability of the model."));
    testStabilityAct->setCheckable(true);
    testStabilityAct->setChecked(glWidget.GetDoStabilityTest());
    connect(testStabilityAct, SIGNAL(triggered()), this,
            SLOT(ToggleStabilityTest()));   

    // generate menu actions
    generateBranchingSetRootAct =
        new QAction(tr("Branching (Set Root Slit)"), this);
    generateBranchingSetRootAct->setStatusTip(
        tr("Set the slit for the root node to generate branching planar "
           "sections."));
    connect(generateBranchingSetRootAct, SIGNAL(triggered()), this,
            SLOT(GenerateBranchingSetRoot()));

    generateBranchingAct = new QAction(tr("&Branching"), this);
    generateBranchingAct->setStatusTip(
        tr("Generate branching planar sections."));
    generateBranchingAct->setShortcut(QKeySequence("Ctrl+B"));
    connect(generateBranchingAct, SIGNAL(triggered()), this,
            SLOT(GenerateBranching()));

    generateLinearAct = new QAction(tr("&Linear"), this);
    generateLinearAct->setStatusTip(
        tr("Generate linearly arranged planar sections."));
    generateLinearAct->setShortcut(QKeySequence("Ctrl+L"));
    connect(generateLinearAct, SIGNAL(triggered()), this,
            SLOT(GenerateLinear()));

    generateBlendAct = new QAction(tr("Bl&end"), this);
    generateBlendAct->setStatusTip(
        tr("Generate an arrangement of blended planar sections."));
    generateBlendAct->setShortcut(QKeySequence("Ctrl+E"));
    connect(generateBlendAct, SIGNAL(triggered()), this, SLOT(GenerateBlend()));

    generateGridAct = new QAction(tr("&Grid"), this);
    generateGridAct->setStatusTip(
        tr("Create a grid of planar sections whose elements fit within a "
           "specified rectangular boundary."));
    generateGridAct->setShortcut(QKeySequence("Ctrl+G"));
    connect(generateGridAct, SIGNAL(triggered()), this, SLOT(GenerateGrid()));

    generateRevolveAct = new QAction(tr("&Revolve"), this);
    generateRevolveAct->setStatusTip(
        tr("Generate a sequence of planar sections which revolve around a base "
           "section."));
    generateRevolveAct->setShortcut(QKeySequence("Ctrl+R"));
    connect(generateRevolveAct, SIGNAL(triggered()), this,
            SLOT(GenerateRevolve()));

    generateMakeCircleAct = new QAction(tr("Set &Circle"), this);
    generateMakeCircleAct->setStatusTip(
        tr("Make the boundary of an existing planar section a circle."));
    connect(generateMakeCircleAct, SIGNAL(triggered()), this,
            SLOT(GenerateMakeCircle()));

    generateMakeRectangleAct = new QAction(tr("Set Rec&tangle"), this);
    generateMakeRectangleAct->setStatusTip(
        tr("Make the boundary of an existing planar section a rectangle."));
    connect(generateMakeRectangleAct, SIGNAL(triggered()), this,
            SLOT(GenerateMakeRectangle()));

    generateMakeRadialAct = new QAction(tr("Set R&adial"), this);
    generateMakeRadialAct->setStatusTip(tr(
        "Make the boundary of an existing planar section a radial pattern."));
    connect(generateMakeRadialAct, SIGNAL(triggered()), this,
            SLOT(GenerateMakeRadial()));

    generateMakeRadialHoleAct = new QAction(tr("Make Radial Holes"), this);
    generateMakeRadialHoleAct->setStatusTip(
        tr("Make a radial set of holes in a selected planar section."));
    connect(generateMakeRadialHoleAct, SIGNAL(triggered()), this,
            SLOT(GenerateMakeRadialHole()));

    generateSlicesAct = new QAction(tr("&Slices"), this);
    generateSlicesAct->setStatusTip(
        tr("Create regularly sampled sections from a template surface."));
    connect(generateSlicesAct, SIGNAL(triggered()), this,
            SLOT(GenerateSlices()));

    generateSurfaceFacetsAct = new QAction(tr("Surface &Facets"), this);
    generateSurfaceFacetsAct->setStatusTip(
        tr("Create sections from the planar facets of a template surface."));
    connect(generateSurfaceFacetsAct, SIGNAL(triggered()), this,
            SLOT(GenerateSurfaceFacets()));

    // physics menu actions
    testPhysicsAct = new QAction(tr("Test &Physics"), this);
    testPhysicsAct->setStatusTip(
        tr("Perform physically-based tests of the model."));
    testPhysicsAct->setCheckable(true);
    testPhysicsAct->setChecked(glWidget.GetDoPhysicsTest());
    connect(testPhysicsAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsTest()));

    CPhysics &physics = glWidget.GetPhysics();

    showPhysicsDeformedAct = new QAction(tr("Show &Deformed"), this);
    showPhysicsDeformedAct->setStatusTip(tr("Show deformed physics."));
    showPhysicsDeformedAct->setCheckable(true);
    showPhysicsDeformedAct->setChecked(physics.GetDrawDeformed());
    connect(showPhysicsDeformedAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsDeformed()));

    showPhysicsSkeletonAct = new QAction(tr("Show &Skeleton"), this);
    showPhysicsSkeletonAct->setStatusTip(tr("Show skeleton used in physics."));
    showPhysicsSkeletonAct->setCheckable(true);
    showPhysicsSkeletonAct->setChecked(physics.GetDrawSkeleton());
    connect(showPhysicsSkeletonAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsSkeleton()));

    showPhysicsForceAct = new QAction(tr("Show &Force"), this);
    showPhysicsForceAct->setStatusTip(tr("Show calculated forces."));
    showPhysicsForceAct->setCheckable(true);
    showPhysicsForceAct->setChecked(physics.GetDrawForce());
    connect(showPhysicsForceAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsForce()));

    showPhysicsSectionAct = new QAction(tr("Show S&ections"), this);
    showPhysicsSectionAct->setStatusTip(
        tr("Show planar sections used by physics."));
    showPhysicsSectionAct->setCheckable(true);
    showPhysicsSectionAct->setChecked(physics.GetDrawSection());
    connect(showPhysicsSectionAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsSection()));

    showPhysicsSectionMomentAct =
        new QAction(tr("Show Section &Moments"), this);
    showPhysicsSectionMomentAct->setStatusTip(
        tr("Show within-plane forces calculated by physics."));
    showPhysicsSectionMomentAct->setCheckable(true);
    showPhysicsSectionMomentAct->setChecked(physics.GetDrawSectionMoment());
    connect(showPhysicsSectionMomentAct, SIGNAL(triggered()), this,
            SLOT(TogglePhysicsSectionMoment()));

    addExternalWeightAct = new QAction(tr("Add &Weight"), this);
    addExternalWeightAct->setStatusTip(
        tr("Add an external mass to the selected plane."));
    connect(addExternalWeightAct, SIGNAL(triggered()), this,
            SLOT(PhysicsAddExternalMass()));

    removeExternalWeightsAct = new QAction(tr("Remove Weights"), this);
    removeExternalWeightsAct->setStatusTip(
        tr("Removes external masses from the selected plane."));
    connect(removeExternalWeightsAct, SIGNAL(triggered()), this,
            SLOT(PhysicsRemoveExternalMasses()));
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
    fileMenu->addAction(saveFlattenedSlabOBJAct);
    fileMenu->addAction(saveSVGAct);
    fileMenu->addAction(saveDXFAct);
    fileMenu->addAction(saveCalibrationAct);
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
    editMenu->addSeparator();
    editMenu->addAction(addHoleBoundaryAct);
    editMenu->addAction(removeHolesBoundaryAct);

    editMenu->addSeparator();
    editMenu->addAction(generateMakeCircleAct);
    editMenu->addAction(generateMakeRadialAct);
    editMenu->addAction(generateMakeRectangleAct);

    editMenu->addSeparator();
    editMenu->addAction(dimensioningToolAct);
    editMenu->addSeparator();
    editMenu->addAction(snapMajorAxisAct);
    editMenu->addAction(resketchCurveAct);
    editMenu->addAction(useMagneticCutsAct);
    editMenu->addAction(useLocalSymmetryAct);
    editMenu->addAction(clipToGroundPlaneAct);

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

    generateMenu = menuBar()->addMenu(tr("&Generate"));
    generateMenu->addAction(generateBlendAct);
    generateMenu->addAction(generateBranchingSetRootAct);
    generateMenu->addAction(generateBranchingAct);
    generateMenu->addAction(generateGridAct);
    generateMenu->addAction(generateLinearAct);
    generateMenu->addAction(generateRevolveAct);
    generateMenu->addAction(generateSlicesAct);
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

    connect(editMenu, SIGNAL(aboutToShow()), this, SLOT(setEditMenuChecks()));
    connect(viewMenu, SIGNAL(aboutToShow()), this, SLOT(setViewMenuChecks()));
    connect(physicsMenu, SIGNAL(aboutToShow()), this,
            SLOT(setPhysicsMenuChecks()));
}

void MainWindow::UpdateWindowTitle()
{
    QString open_filename = glWidget.GetOpenFilename();
    if (open_filename.length() > 0) {
        setWindowTitle(glWidget.GetOpenFilename() + " - " + window_title);
    } else {
        setWindowTitle(window_title);
    }
}

void MainWindow::NewPlaneSketch()
{
    glWidget.NewPlaneSketch();
    UpdateWindowTitle();
}

void MainWindow::LoadTemplateOBJ() { glWidget.LoadTemplateOBJ(); }

void MainWindow::LoadTemplateCurve() { glWidget.LoadTemplateCurve(); }

void MainWindow::LoadTemplateImage() { glWidget.LoadTemplateImage(); }

void MainWindow::LoadPlaneSketch()
{
    glWidget.LoadPlaneSketch();
    UpdateWindowTitle();
}

void MainWindow::SavePlaneSketch() {
    glWidget.SavePlaneSketch();
    UpdateWindowTitle();
}

void MainWindow::SaveSliceOBJ() { glWidget.SaveSliceOBJ(); }

void MainWindow::SaveSlabOBJ() { glWidget.SaveSlabOBJ(); }

void MainWindow::SaveFlattenedSlabOBJ() { glWidget.SaveFlattenedSlabOBJ(); }

void MainWindow::SaveSurfaceOBJ() { glWidget.SaveSurfaceOBJ(); }

void MainWindow::SaveSVG() { glWidget.SaveSVG(); }

void MainWindow::SaveDXF() { glWidget.SaveDXF(); }

void MainWindow::SavePhysicsOutput() { glWidget.SavePhysicsOutput(); }

void MainWindow::SaveCalibration() { glWidget.SaveCalibration(); }

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

void MainWindow::ToggleClipToGroundPlane()
{
    glWidget.SetClipToGroundPlane(!glWidget.GetClipToGroundPlane());
}

void MainWindow::TogglePenMode()
{
    glWidget.SetPenModeOn(!glWidget.GetPenModeOn());
}

void MainWindow::StartDimensioningTool() { glWidget.StartDimensioningTool(); }

void MainWindow::Undo()
{
    glWidget.Undo(OP_UNDO);
}

void MainWindow::Redo() { glWidget.Redo(); }

void MainWindow::Delete() { glWidget.DeleteSelected(); }

void MainWindow::Transform() { glWidget.Transform(); }

void MainWindow::CopyMirrorX() { glWidget.CopyMirrorX(); }

void MainWindow::CopyRotateY() { glWidget.CopyRotateY(); }

void MainWindow::CopyMirrorZ() { glWidget.CopyMirrorZ(); }

void MainWindow::CopyDuplicate() { glWidget.CopyDuplicate(); }

void MainWindow::SnapToMajorAxis() { glWidget.SnapToMajorAxis(); }

void MainWindow::AddHoleBoundary() { glWidget.AddHoleBoundary(); }

void MainWindow::RemoveHolesBoundary() { glWidget.RemoveHolesBoundary(); }

void MainWindow::ResketchCurve() { glWidget.ResketchCurve(); }

void MainWindow::CreateSurfacePatches() { glWidget.CreateSurfacePatches(); }

void MainWindow::DeleteSurfacePatches() { glWidget.DeleteSurfacePatches(); }

void MainWindow::ViewIso1() { glWidget.SetViewIso1(); }

void MainWindow::ViewIso2() { glWidget.SetViewIso2(); }

void MainWindow::ViewIso3() { glWidget.SetViewIso3(); }

void MainWindow::ViewIso4() { glWidget.SetViewIso4(); }

void MainWindow::ViewX() { glWidget.SetViewX(); }

void MainWindow::ViewY() { glWidget.SetViewY(); }

void MainWindow::ViewZ() { glWidget.SetViewZ(); }

void MainWindow::ViewPart() { glWidget.SetViewPart(); }

void MainWindow::ToggleCyclesTest() { glWidget.ToggleDoCyclesTest(); }

void MainWindow::ToggleConnectedTest() { glWidget.ToggleDoConnectedTest(); }

void MainWindow::ToggleStabilityTest() { glWidget.ToggleDoStabilityTest(); }

void MainWindow::ToggleShowTNBFrames() { glWidget.ToggleShowTNBFrames(); }

void MainWindow::ToggleShowShadow() { glWidget.ToggleShowShadow(); }

void MainWindow::ToggleShowTemplates() { glWidget.ToggleDrawTemplates(); }

void MainWindow::GenerateBranchingSetRoot()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.DoGenerateBranchingSetRoot();
}

void MainWindow::GenerateBranching()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.DoGenerateBranching();
}

void MainWindow::GenerateLinear()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.StartGenerateLinear();
}

void MainWindow::GenerateBlend()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.StartGenerateBlend();
}

void MainWindow::GenerateGrid()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.StartGenerateGrid();
}

void MainWindow::GenerateRevolve()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.StartGenerateRevolve();
}

void MainWindow::GenerateMakeCircle()
{
    if (!docks[0]->isVisible()) {
        openEditWidget();
    }
    glWidget.DoGenerateMakeCircle();
}

void MainWindow::GenerateMakeRectangle()
{
    if (!docks[0]->isVisible()) {
        openEditWidget();
    }
    glWidget.DoGenerateMakeRectangle();
}

void MainWindow::GenerateMakeRadial()
{
    if (!docks[0]->isVisible()) {
        openEditWidget();
    }
    glWidget.SetSelectedAsRadial();
}

void MainWindow::GenerateMakeRadialHole()
{
    if (!docks[0]->isVisible()) {
        openEditWidget();
    }
    glWidget.DoGenerateMakeRadialHole();
}

void MainWindow::GenerateSurfaceFacets()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.DoGenerateSurfaceFacets();
}

void MainWindow::GenerateSlices()
{
    if (!docks[1]->isVisible()) {
        openGenerateWidget();
    }
    glWidget.StartGenerateSlices();
}

void MainWindow::TogglePhysicsTest()
{
    glWidget.ToggleDoPhysicsTest();
}

void MainWindow::TogglePhysicsDeformed()
{
    glWidget.ToggleDrawDeformed();
}

void MainWindow::TogglePhysicsSkeleton()
{
    glWidget.ToggleDrawSkeleton();
}

void MainWindow::TogglePhysicsForce()
{
    glWidget.ToggleDrawForce();
}

void MainWindow::TogglePhysicsSection()
{
    glWidget.ToggleDrawSection();
}

void MainWindow::TogglePhysicsSectionMoment()
{
    glWidget.ToggleDrawMoments();
}

void MainWindow::PhysicsAddExternalMass()
{
    glWidget.DoPhysicsAddWeight();
}

void MainWindow::PhysicsRemoveExternalMasses()
{
    glWidget.DoPhysicsRemoveWeights();
}

// UI slots
void MainWindow::openEditWidget()
{
    for (int i = 0; i < 5; i++) {
        if (i != 0) {
            docks[i]->setVisible(false);
            openDock[i]->setChecked(false);
        }
    }
    docks[0]->toggleViewAction()->trigger();
}

void MainWindow::openGenerateWidget()
{
    for (int i = 0; i < 5; i++) {
        if (i != 1) {
            docks[i]->setVisible(false);
            openDock[i]->setChecked(false);
        }
    }
    docks[1]->toggleViewAction()->trigger();
}

void MainWindow::openGuidesWidget()
{
    for (int i = 0; i < 5; i++) {
        if (i != 2) {
            docks[i]->setVisible(false);
            openDock[i]->setChecked(false);
        }
    }
    docks[2]->toggleViewAction()->trigger();
}

void MainWindow::openPhysicsWidget()
{
    for (int i = 0; i < 5; i++) {
        if (i != 3) {
            docks[i]->setVisible(false);
            openDock[i]->setChecked(false);
        }
    }
    docks[3]->toggleViewAction()->trigger();
}

void MainWindow::openViewWidget()
{
    for (int i = 0; i < 5; i++) {
        if (i != 4) {
            docks[i]->setVisible(false);
            openDock[i]->setChecked(false);
        }
    }
    docks[4]->toggleViewAction()->trigger();
}

void MainWindow::setEditMenuChecks()
{
    useMagneticCutsAct->setCheckable(glWidget.GetDoMagneticCuts());
}

void MainWindow::setViewMenuChecks()
{
    viewTNBFramesAct->setChecked(glWidget.GetShowTNBFrames());
    viewShadowAct->setChecked(glWidget.GetShowShadow());
    viewTemplatesAct->setChecked(glWidget.GetShowTemplates());
    testCyclesAct->setChecked(glWidget.GetDoCyclesTest());
    testConnectedAct->setChecked(glWidget.GetDoConnectedTest());
    testStabilityAct->setChecked(glWidget.GetDoStabilityTest());
}

void MainWindow::setPhysicsMenuChecks()
{
    CPhysics &physics = glWidget.GetPhysics();
    showPhysicsDeformedAct->setChecked(physics.GetDrawDeformed());
    showPhysicsSkeletonAct->setChecked(physics.GetDrawSkeleton());
    showPhysicsForceAct->setChecked(physics.GetDrawForce());
    showPhysicsSectionAct->setChecked(physics.GetDrawSection());
    showPhysicsSectionMomentAct->setChecked(physics.GetDrawSectionMoment());
    testPhysicsAct->setChecked(glWidget.GetDoPhysicsTest());
}

