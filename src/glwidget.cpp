﻿#include "glwidget.h"

GLWidget::GLWidget() :
    editWidget(NULL),
    genWidget(NULL),
    guidesWidget(NULL),
    physicsWidget(NULL),
    viewWidget(NULL)
{

    setMouseTracking(true);
    setMinimumSize(800, 700);
    setAutoFillBackground(false);
    setFocusPolicy(Qt::ClickFocus);

    grid_size = 10;

    max_undo_sections = 1000;
    undo_index = -1;

    x_symmetry = false;
    y_symmetry = false;
    z_symmetry = false;

    do_local_symmetry = true;
    do_cycles_test = false;
    do_connected_test = false;
    do_stability_test = false;
    do_physics_test = false;

    pen_mode = false;

    do_show_tnb_frames = false;
    do_show_shadow = true;
    do_show_templates = true;

    do_magnetic_cuts = false;
    magnet_strength = 1.0f;

    deadzone_radius = 0.4f;

    template_image_tex = 0;
    template_pos = QVector2D(0, 0);
    template_rotation = 0.0f;
    template_scale = 5.0f;
    template_flipx = false;

    generate_grid_sizex = 2.0f;
    generate_grid_sizey = 2.0f;
    generate_grid_staplesize = 0.25f;
    generate_linear_spacing = 0.5f;
    generate_linear_scalex = false;
    generate_linear_scaley = false;
    generate_slices_spacing = 0.5f;
    generate_slices_x = true;
    generate_slices_y = false;
    generate_slices_z = true;
    generate_blend_sections = 10;    
    generate_revolve_sections = 10;
    generate_branching_scalechild = 0.8f;
    generate_radial_sectors = 10;
    for (int i=0; i<9; ++i) {
        generate_radial_params[i] = 1.0f;
    }
    generate_surfacefacets_teeth = false;

    metres_per_unit = 0.0254; //0.01 means that each unit is 1cm, 0.0254 that means each unit is 1 inch

    quality_samples = 15;

    //slab_thickness = 0.06f; //this is in distance of units
    slab_thickness = 0.125f; //this is for 1/8" acrylic
    calibration_factor = 1.0f;
    physics_new_weight_mass = 1.0f; //in kg
    physics_material_density = 1190.0f; //in kg/m^3
    //physics_max_stress = 50000.0f; //in Pa (acrylic is 60 MPa, or 60 million Pa)
    //physics_max_stress = 60000000.0f; //in Pa (acrylic is 60 MPa, or 60 million Pa)
    physics_max_stress = 6000000.0f; //this is 1/10 of the max... note that red will appear at 30MPa
    physics.SetMaximumStress(physics_max_stress);

    template_cut_snap_distance_3d = 0.5f;

    section_error_tolerance = 0.015f;
    //section_error_tolerance = 0.01f;
    section_error_tolerance_template = 0.005f;

    animate_until = QDateTime::currentDateTime();

    cam_animate_duration = 800.0f;
    cam_animate_angle = 80.0f;
    //cam_animate_angle = 90.0f;
    cam_lookat_distance = 100.0f;

    //brush_type = BOUNDARY;

    state = STATE_NONE;
    last_op = OP_NONE;

    update_sections_disp_list = false;
    sections_disp_list = 0;

    /*
    PlanarDisc first_disc;
    TNBFrame frame;
    frame.SetN(vec3(1, 0, 0));
    frame.SetCentre(vec3(0, 0, 0));
    first_disc.SetFrame(frame);
    first_disc.SetActiveEdit(false);

    discs.push_back(first_disc);
    disc_selected = 0;
    */

    ClearAll();

    //StartAnimation(cam_animate_duration);

    default_lookat = QVector3D(0, grid_size/3, 0);

    ResetCamera();

    animate_timer.setSingleShot(false);
    animate_timer.start(0);

    connect(&animate_timer, SIGNAL(timeout()), this, SLOT(update()));
    //LoadTemplateOBJ();

    //initialize the physics thing
    //physics.SetProblem();
    //physics.Solve(100,20);
    //physics.PrintJointForce();

    current_tool_state = TOOLSTATE_DEFAULT;

    selections_per_gen_type[GENSTATE_BLEND] = 3;
    selections_per_gen_type[GENSTATE_LINEAR] = 2;
    selections_per_gen_type[GENSTATE_REVOLVE] = 2;    
    selections_per_gen_type[GENSTATE_GRID] = 1;
    selections_per_gen_type[GENSTATE_BRANCH] = 1;
    selections_per_gen_type[GENSTATE_SLICES] = 0;
    selections_per_gen_type[GENSTATE_SURFACE_FACETS] = 0;

}

GLWidget::~GLWidget()
{

    if (editWidget != NULL) {
        delete editWidget;
        editWidget = NULL;
    }
    if (genWidget != NULL) {
        delete genWidget;
        genWidget = NULL;
    }
    if (guidesWidget != NULL) {
        delete guidesWidget;
        guidesWidget = NULL;
    }
    if (physicsWidget != NULL) {
        delete physicsWidget;
        physicsWidget = NULL;
    }
    if (viewWidget != NULL) {
        delete viewWidget;
        viewWidget = NULL;
    }   

    if (template_image_tex > 0) {
        glDeleteTextures(1, &template_image_tex);
        template_image_tex = 0;
    }

}

CPhysics & GLWidget::GetPhysics()
{
    return physics;
}

void GLWidget::initializeGL()
{

    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);

    glShadeModel(GL_SMOOTH);

    glClearColor(1, 1, 1, 1);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

}

QWidget * GLWidget::GetEditWidget()
{
    if (editWidget != NULL) {
        return editWidget;
    }

    grid_sizex_label = new QLabel(QString::number(GetGenerateGridSizeX()));
    grid_sizey_label = new QLabel(QString::number(GetGenerateGridSizeX()));
    grid_staple_label = new QLabel(QString::number(GetGenerateGridStapleSize()));
    linear_spacing_label = new QLabel(QString::number(GetGenerateLinearSpacing()));
    slices_spacing_label = new QLabel(QString::number(GetGenerateSlicesSpacing()));
    num_blend_sections_label = new QLabel(QString::number(GetGenerateBlendSections()));
    num_revolve_sections_label = new QLabel(QString::number(GetGenerateRevolveSections()));
    branching_scalechild_label = new QLabel(QString::number(GetGenerateBranchingScaleChild()));
    radial_sectors_label = new QLabel(QString::number(GetGenerateRadialSectors()));

    grid_sizex_label->setFixedWidth(30);
    grid_sizey_label->setFixedWidth(30);
    grid_staple_label->setFixedWidth(30);
    linear_spacing_label->setFixedWidth(30);
    slices_spacing_label->setFixedWidth(30);
    num_blend_sections_label->setFixedWidth(30);
    num_revolve_sections_label->setFixedWidth(30);
    branching_scalechild_label->setFixedWidth(30);
    radial_sectors_label->setFixedWidth(30);

    QFormLayout * editWidgetLayout = new QFormLayout();


    QPushButton * transformButton = new QPushButton("Transform");
    connect(transformButton, SIGNAL(clicked()), this, SLOT(Transform()));
    transformButton->setMinimumHeight(40);

    editWidgetLayout->addWidget(transformButton);

    QPushButton * snapButton = new QPushButton("Snap to Major Axis");
    connect(snapButton, SIGNAL(clicked()), this, SLOT(SnapToMajorAxis()));

    editWidgetLayout->addWidget(snapButton);

    QPushButton * resketchButton = new QPushButton("Re-sketch Boundary");
    connect(resketchButton, SIGNAL(clicked()), this, SLOT(ResketchCurve()));

    editWidgetLayout->addWidget(resketchButton);

    QPushButton * deleteButton = new QPushButton("Delete Section");
    deleteButton->setStyleSheet("QPushButton{ background-color: #40FA42; color:black}"); // for any delete button
    connect(deleteButton, SIGNAL(clicked()), this, SLOT(DeleteSelected()));

    editWidgetLayout->addWidget(deleteButton);



    // Copy group

    QPushButton * copyButton1 = new QPushButton("Duplicate");
    connect(copyButton1, SIGNAL(clicked()), this, SLOT(CopyDuplicate()));

    QPushButton * copyButton2 = new QPushButton("Mirror X");
    connect(copyButton2, SIGNAL(clicked()), this, SLOT(CopyMirrorX()));

    QPushButton * copyButton3 = new QPushButton("Rotate Y");
    connect(copyButton3, SIGNAL(clicked()), this, SLOT(CopyRotateY()));

    QPushButton * copyButton4 = new QPushButton("Mirror Z");
    connect(copyButton4, SIGNAL(clicked()), this, SLOT(CopyMirrorZ()));

    QGroupBox * copy_groupbox = new QGroupBox(tr("Copy Selected Section"));
    QGridLayout * copy_layout = new QGridLayout;
    copy_layout->addWidget(copyButton1,0,0);
    copy_layout->addWidget(copyButton2,0,1);
    copy_layout->addWidget(copyButton3,1,0);
    copy_layout->addWidget(copyButton4,1,1);
    copy_groupbox->setLayout(copy_layout);
    editWidgetLayout->addRow(copy_groupbox);


    // Hole group

    QPushButton * holeButton1 = new QPushButton("Add Hole");
    connect(holeButton1, SIGNAL(clicked()), this, SLOT(AddHoleBoundary()));

    QPushButton * holeButton2 = new QPushButton("Remove All Holes");
    connect(holeButton2, SIGNAL(clicked()), this, SLOT(RemoveHolesBoundary()));
    holeButton2->setStyleSheet("QPushButton{ background-color: #40FA42; color:black}"); // for any delete button

    QGroupBox * hole_groupbox = new QGroupBox(tr("Section Holes"));
    QGridLayout * hole_layout = new QGridLayout;
    hole_layout->addWidget(holeButton1,0,0);
    hole_layout->addWidget(holeButton2,1,0);
    hole_groupbox->setLayout(hole_layout);
    editWidgetLayout->addRow(hole_groupbox);




    // Radial group

    QPushButton * radialButton = new QPushButton("Set Radial");
    connect(radialButton, SIGNAL(clicked()), this, SLOT(SetSelectedAsRadial()));

    QPushButton * radialButton2 = new QPushButton("Remove Radial");
    connect(radialButton2, SIGNAL(clicked()), this, SLOT(RemoveRadial()));

    QSlider * radial_sectors_slider = new QSlider();
    radial_sectors_slider->setRange(1, 40);
    radial_sectors_slider->setOrientation(Qt::Horizontal);
    radial_sectors_slider->setValue(this->GetGenerateRadialSectors());
    connect(radial_sectors_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateRadialSectors(int)));


    QGroupBox * radial_groupbox = new QGroupBox(tr("Radial Operations"));
    QGridLayout * radial_layout = new QGridLayout;
    radial_layout->addWidget(radialButton,0,0,1,3);

    radial_layout->addWidget(radialButton2,1,0,1,3);
    radial_layout->addWidget(new QLabel("Sectors"), 2, 0);
    radial_layout->addWidget(radial_sectors_slider, 2, 1);
    radial_layout->addWidget(radial_sectors_label, 2, 2);
    radial_groupbox->setLayout(radial_layout);
    editWidgetLayout->addRow(radial_groupbox);

    editWidget = new QTabWidget();
    editWidget->setMinimumWidth(225);
    editWidget->setLayout(editWidgetLayout);

    return editWidget;

}

QWidget * GLWidget::GetGenerateWidget()
{

    if (genWidget != NULL) {
        return genWidget;
    }

    QFormLayout * generateWidgetLayout = new QFormLayout();

    // Blend group

    QPushButton * blendButton = new QPushButton("Blend");
    connect(blendButton, SIGNAL(clicked()), this, SLOT(StartGenerateBlend()));

    QSlider * num_blend_sections_slider = new QSlider();
    num_blend_sections_slider->setRange(3, 40);
    num_blend_sections_slider->setOrientation(Qt::Horizontal);
    num_blend_sections_slider->setValue(GetGenerateBlendSections());
    connect(num_blend_sections_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateBlendSections(int)));

    QGroupBox * blend_groupbox = new QGroupBox(tr("Blend Operation"));
    QGridLayout * blend_layout = new QGridLayout;
    blend_layout->addWidget(blendButton,0,0,1,3);
    blend_layout->addWidget(new QLabel("Sections"), 1, 0);
    blend_layout->addWidget(num_blend_sections_slider, 1, 1);
    blend_layout->addWidget(num_blend_sections_label, 1, 2);
    blend_groupbox->setLayout(blend_layout);
    generateWidgetLayout->addRow(blend_groupbox);



    // Branch group

    QPushButton * setRootButton = new QPushButton("Set Root");
    connect(setRootButton, SIGNAL(clicked()), this, SLOT(DoGenerateBranchingSetRoot()) );

    QPushButton * branchButton = new QPushButton("Add Branch");
    connect(branchButton, SIGNAL(clicked()), this, SLOT(DoGenerateBranching()) );

    QSlider * branching_scalechild_slider = new QSlider();
    branching_scalechild_slider->setRange(1, 200);
    branching_scalechild_slider->setOrientation(Qt::Horizontal);
    branching_scalechild_slider->setValue(this->GetGenerateBranchingScaleChild()*100);
    connect(branching_scalechild_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateBranchingScaleChild(int)));

    QGroupBox * branch_groupbox = new QGroupBox(tr("Branching Operation"));
    QGridLayout * branch_layout = new QGridLayout;
    branch_layout->addWidget(setRootButton, 0, 0,1,3);
    branch_layout->addWidget(branchButton, 0, 3,1,3);
    branch_layout->addWidget(new QLabel("Child Scale"), 1, 0);
    branch_layout->addWidget(branching_scalechild_slider, 1, 1,1,4);
    branch_layout->addWidget(branching_scalechild_label, 1, 5);
    branch_groupbox->setLayout(branch_layout);
    generateWidgetLayout->addRow(branch_groupbox);



    // Grid group

    QPushButton * gridButton = new QPushButton("Grid");
    connect(gridButton, SIGNAL(clicked()), this, SLOT(StartGenerateGrid()));

    QSlider * grid_sizex_slider = new QSlider();
    grid_sizex_slider->setRange(1, 50);
    grid_sizex_slider->setValue(GetGenerateGridSizeX() * 10);
    grid_sizex_slider->setOrientation(Qt::Horizontal);
    connect(grid_sizex_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateGridSizeX(int)));

    QSlider * grid_sizey_slider = new QSlider();
    grid_sizey_slider->setRange(1, 50);
    grid_sizey_slider->setOrientation(Qt::Horizontal);
    grid_sizey_slider->setValue(GetGenerateGridSizeY() * 10);
    connect(grid_sizey_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateGridSizeY(int)));

    QSlider * grid_staplesize_slider = new QSlider();
    grid_staplesize_slider->setRange(1, 40);
    grid_staplesize_slider->setOrientation(Qt::Horizontal);
    grid_staplesize_slider->setValue(GetGenerateGridStapleSize() * 20);
    connect(grid_staplesize_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateGridStapleSize(int)));

    QGroupBox * grid_groupbox = new QGroupBox(tr("Grid Operation"));
    QGridLayout * grid_layout = new QGridLayout;
    grid_layout->addWidget(gridButton, 0, 0, 1, 3);
    grid_layout->addWidget(new QLabel("Cell Width"), 1, 0);
    grid_layout->addWidget(grid_sizex_slider, 1, 1);
    grid_layout->addWidget(grid_sizex_label, 1, 2);
    grid_layout->addWidget(new QLabel("Cell Height"), 2, 0);
    grid_layout->addWidget(grid_sizey_slider, 2, 1);
    grid_layout->addWidget(grid_sizey_label, 2, 2);
    grid_layout->addWidget(new QLabel("Staple Size"), 3, 0);
    grid_layout->addWidget(grid_staplesize_slider, 3, 1);
    grid_layout->addWidget(grid_staple_label, 3, 2);
    grid_groupbox->setLayout(grid_layout);
    generateWidgetLayout->addRow(grid_groupbox);



    // Linear group

    QPushButton * linearButton = new QPushButton("Linear");
    connect(linearButton, SIGNAL(clicked()), this, SLOT(StartGenerateLinear()));

    QSlider * linear_spacing_slider = new QSlider();
    linear_spacing_slider->setRange(1, 50);
    linear_spacing_slider->setOrientation(Qt::Horizontal);
    linear_spacing_slider->setValue(GetGenerateLinearSpacing()*10);
    connect(linear_spacing_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateLinearSpacing(int)));

    QCheckBox * linear_scalex_checkbox = new QCheckBox("Scale X");
    linear_scalex_checkbox ->setChecked(GetGenerateLinearScaleX());
    connect(linear_scalex_checkbox , SIGNAL(toggled(bool)), this, SLOT(SetGenerateLinearScaleX(bool)));

    QCheckBox * linear_scaley_checkbox = new QCheckBox("Scale Y");
    linear_scaley_checkbox ->setChecked(GetGenerateLinearScaleY());
    connect(linear_scaley_checkbox , SIGNAL(toggled(bool)), this, SLOT(SetGenerateLinearScaleY(bool)));

    QGroupBox * linear_groupbox = new QGroupBox(tr("Linear Operation"));
    QGridLayout * linear_layout = new QGridLayout;
    linear_layout->addWidget(linearButton,0,0,1,4);
    linear_layout->addWidget(new QLabel("Spacing"), 1, 0);
    linear_layout->addWidget(linear_spacing_slider, 1, 1, 1, 2);
    linear_layout->addWidget(linear_spacing_label, 1, 3);
    linear_layout->addWidget(linear_scalex_checkbox, 2, 0, 1, 2);
    linear_layout->addWidget(linear_scaley_checkbox, 2, 2, 1, 2);
    linear_groupbox->setLayout(linear_layout);
    generateWidgetLayout->addRow(linear_groupbox);                


    // Revolve group

    QPushButton * revolveButton = new QPushButton("Revolve");
    connect(revolveButton, SIGNAL(clicked()), this, SLOT(StartGenerateRevolve()));

    QSlider * num_revolve_sections_slider = new QSlider();
    num_revolve_sections_slider->setRange(2, 40);
    num_revolve_sections_slider->setOrientation(Qt::Horizontal);
    num_revolve_sections_slider->setValue(GetGenerateRevolveSections());
    connect(num_revolve_sections_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateRevolveSections(int)));

    QGroupBox * revolve_groupbox = new QGroupBox(tr("Revolve Operation"));
    QGridLayout * revolve_layout = new QGridLayout;
    revolve_layout->addWidget(revolveButton,0,0,1,3);
    revolve_layout->addWidget(new QLabel("Sections"), 1, 0);
    revolve_layout->addWidget(num_revolve_sections_slider, 1, 1);
    revolve_layout->addWidget(num_revolve_sections_label, 1, 2);
    revolve_groupbox->setLayout(revolve_layout);
    generateWidgetLayout->addRow(revolve_groupbox);


    // Slices group

    QPushButton * slicesButton = new QPushButton("Slices");
    connect(slicesButton, SIGNAL(clicked()), this, SLOT(StartGenerateSlices()));

    QSlider * slices_spacing_slider = new QSlider();
    slices_spacing_slider->setRange(1, 50);
    slices_spacing_slider->setOrientation(Qt::Horizontal);
    slices_spacing_slider->setValue(GetGenerateSlicesSpacing()*10);
    connect(slices_spacing_slider, SIGNAL(sliderReleased()), this, SLOT(UpdateGenerateSlicesSpacing()));
    connect(slices_spacing_slider, SIGNAL(valueChanged(int)), this, SLOT(SetGenerateSlicesSpacing(int)));

    QCheckBox * slices_x_checkbox = new QCheckBox("X");
    slices_x_checkbox ->setChecked(GetGenerateSlicesX());
    connect(slices_x_checkbox , SIGNAL(toggled(bool)), this, SLOT(SetGenerateSlicesX(bool)));

    QCheckBox * slices_y_checkbox = new QCheckBox("Y");
    slices_y_checkbox ->setChecked(GetGenerateSlicesY());
    connect(slices_y_checkbox , SIGNAL(toggled(bool)), this, SLOT(SetGenerateSlicesY(bool)));

    QCheckBox * slices_z_checkbox = new QCheckBox("Z");
    slices_z_checkbox ->setChecked(GetGenerateSlicesZ());
    connect(slices_z_checkbox , SIGNAL(toggled(bool)), this, SLOT(SetGenerateSlicesZ(bool)));

    QGroupBox * slices_groupbox = new QGroupBox(tr("Slice Operation"));
    QGridLayout * slices_layout = new QGridLayout;
    slices_layout->addWidget(slicesButton,0,0,1,4);
    slices_layout->addWidget(new QLabel("Spacing"), 1, 0);
    slices_layout->addWidget(slices_spacing_slider, 1, 1, 1, 2);
    slices_layout->addWidget(slices_spacing_label, 1, 3);
    slices_layout->addWidget(slices_x_checkbox, 2, 0, 1, 1);
    slices_layout->addWidget(slices_y_checkbox, 2, 1, 1, 1);
    slices_layout->addWidget(slices_z_checkbox, 2, 3, 1, 1);
    slices_groupbox->setLayout(slices_layout);
    generateWidgetLayout->addRow(slices_groupbox);

    // Surface facets group
    QPushButton * surfaceFacetsButton = new QPushButton("Surface Facets");
    connect(surfaceFacetsButton, SIGNAL(clicked()), this, SLOT(DoGenerateSurfaceFacets()));

    QCheckBox * surfacefacets_teeth_checkbox = new QCheckBox("Teeth");
    surfacefacets_teeth_checkbox ->setChecked(GetGenerateSurfaceFacetsTeeth());
    connect(surfacefacets_teeth_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetGenerateSurfaceFacetsTeeth(bool)));

    QGroupBox * surfacefacet_groupbox = new QGroupBox(tr("Surface Facets Operation"));
    QGridLayout * surfacefacet_layout = new QGridLayout;
    surfacefacet_layout->addWidget(surfaceFacetsButton,0,0,1,4);
    surfacefacet_layout->addWidget(surfacefacets_teeth_checkbox, 2, 0, 1, 1);
    surfacefacet_groupbox->setLayout(surfacefacet_layout);
    generateWidgetLayout->addRow(surfacefacet_groupbox);

    //Finally, create the tabwidget

    genWidget = new QTabWidget();
    genWidget->setMinimumWidth(225);
    genWidget->setLayout(generateWidgetLayout);

    return genWidget;
}

QWidget * GLWidget::GetGuidesWidget()
{
    if (guidesWidget != NULL) {
        return guidesWidget;
    }   

    QFormLayout * guidesWidgetLayout = new QFormLayout();

    QPushButton * dimensioningToolButton = new QPushButton("Dimensioning Tool");
    connect(dimensioningToolButton, SIGNAL(clicked()), this, SLOT(StartDimensioningTool()));
    dimensioningToolButton->setMinimumHeight(40);
    guidesWidgetLayout->addRow(dimensioningToolButton);


    QDoubleSpinBox * thick_spinbox = new QDoubleSpinBox();
    thick_spinbox->setRange(0.001, 0.5);
    thick_spinbox->setDecimals(4);
    thick_spinbox->setSingleStep(0.001);
    thick_spinbox->setValue(slab_thickness);
    //thick_spinbox->setFocusPolicy(Qt::ClickFocus);
    connect(thick_spinbox, SIGNAL(valueChanged(double)), this, SLOT(SetSlabThickness(double)));

    QSlider * calibration_slider = new QSlider();
    calibration_slider->setRange(50, 100);
    calibration_slider->setOrientation(Qt::Horizontal);
    calibration_slider->setValue(calibration_factor * 100.0f);
    connect(calibration_slider, SIGNAL(valueChanged(int)), this, SLOT(SetCalibrationFactor(int)));

    QSlider * quality_slider = new QSlider();
    quality_slider->setRange(1, 30);
    quality_slider->setOrientation(Qt::Horizontal);
    quality_slider->setValue(quality_samples);
    connect(quality_slider, SIGNAL(valueChanged(int)), this, SLOT(SetQualitySamples(int)));

    thick_label = new QLabel(QString::number(thick_spinbox->value()) + QString(" units"));
    calibration_label = new QLabel(QString::number(calibration_slider->value()) + QString("%"));
    quality_label = new QLabel(QString::number(quality_slider->value()));

    QGroupBox * dim_groupbox = new QGroupBox(tr("Dimensions"));
    QGridLayout *dim_layout = new QGridLayout;
    QRadioButton *radio_inch = new QRadioButton(tr("&inch"));
    radio_inch->setChecked(true);
    connect(radio_inch, SIGNAL(clicked()), this, SLOT(SetUnitsToInches()));
    QRadioButton *radio_cm = new QRadioButton(tr("&cm"));
    connect(radio_cm, SIGNAL(clicked()), this, SLOT(SetUnitsToCentimetres()));
    dim_layout->addWidget(new QLabel("Units"), 0, 0);
    dim_layout->addWidget(radio_inch, 0, 1);
    dim_layout->addWidget(radio_cm, 0, 2);
    dim_layout->addWidget(new QLabel("Thickness"), 1, 0);
    dim_layout->addWidget(thick_spinbox, 1, 1, 1, 2);
    dim_layout->addWidget(new QLabel("Calibration"), 2, 0);
    dim_layout->addWidget(calibration_slider, 2, 1);
    dim_layout->addWidget(calibration_label, 2, 2);
    dim_layout->addWidget(new QLabel("Quality"), 3, 0);
    dim_layout->addWidget(quality_slider, 3, 1);
    dim_layout->addWidget(quality_label, 3, 2);
    dim_groupbox->setLayout(dim_layout);
    guidesWidgetLayout->addRow(dim_groupbox);




    QPushButton * loadObjButton = new QPushButton("Load OBJ Template");
    connect(loadObjButton, SIGNAL(clicked()), this, SLOT(LoadTemplateOBJ()));

    QCheckBox * magnetic_cuts_checkbox = new QCheckBox("Magnetic Cuts");
    magnetic_cuts_checkbox->setChecked(do_magnetic_cuts);
    connect(magnetic_cuts_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDoMagneticCuts(bool)));


    QGroupBox * template_groupbox = new QGroupBox(tr("3D Template"));
    QVBoxLayout *template_layout = new QVBoxLayout;
    template_layout->addWidget(loadObjButton);
    template_layout->addWidget(magnetic_cuts_checkbox);
    template_groupbox->setLayout(template_layout);
    guidesWidgetLayout->addRow(template_groupbox);



    QPushButton * loadImageButton = new QPushButton("Load Image Template");
    connect(loadImageButton, SIGNAL(clicked()), this, SLOT(LoadTemplateImage()));

    QSlider * imagex_spinbox = new QSlider();
    imagex_spinbox->setRange(-250, 250);
    imagex_spinbox->setValue(template_pos.x() * 20.0f);
    imagex_spinbox->setOrientation(Qt::Horizontal);
    connect(imagex_spinbox, SIGNAL(valueChanged(int)), this, SLOT(SetTemplateImageX(int)));

    QSlider * imagey_spinbox = new QSlider();
    imagey_spinbox->setRange(-250, 250);
    imagey_spinbox->setValue(template_pos.y() * 20.0f);
    imagey_spinbox->setOrientation(Qt::Horizontal);
    connect(imagey_spinbox, SIGNAL(valueChanged(int)), this, SLOT(SetTemplateImageY(int)));

    QSlider * imagerotate_spinbox = new QSlider();
    imagerotate_spinbox->setRange(0, 360);
    imagerotate_spinbox->setValue(template_rotation);
    imagerotate_spinbox->setOrientation(Qt::Horizontal);
    connect(imagerotate_spinbox, SIGNAL(valueChanged(int)), this, SLOT(SetTemplateImageRotate(int)));

    QSlider * imagescale_spinbox = new QSlider();
    imagescale_spinbox->setRange(1, 200);
    imagescale_spinbox->setValue(template_scale * 20.0f);
    imagescale_spinbox->setOrientation(Qt::Horizontal);
    connect(imagescale_spinbox, SIGNAL(valueChanged(int)), this, SLOT(SetTemplateImageScale(int)));

    QCheckBox * imageflipx_checkbox = new QCheckBox("Horizontal Flip");
    imageflipx_checkbox->setChecked(template_flipx);
    connect(imageflipx_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetTemplateImageFlipX(bool)));

    QGroupBox * image_groupbox = new QGroupBox(tr("Image Template"));
    QGridLayout *image_layout = new QGridLayout;
    image_layout->addWidget(loadImageButton, 0, 0, 1, 2);
    image_layout->addWidget(imageflipx_checkbox, 1, 0, 1, 2);
    image_layout->addWidget(new QLabel("X"), 2, 0);
    image_layout->addWidget(imagex_spinbox, 2, 1);
    image_layout->addWidget(new QLabel("Y"), 3, 0);
    image_layout->addWidget(imagey_spinbox, 3, 1);
    image_layout->addWidget(new QLabel("Rotate"), 4, 0);
    image_layout->addWidget(imagerotate_spinbox, 4, 1);
    image_layout->addWidget(new QLabel("Scale"), 5, 0);
    image_layout->addWidget(imagescale_spinbox, 5, 1);
    image_groupbox->setLayout(image_layout);
    guidesWidgetLayout->addRow(image_groupbox);


    toggleTemplatesButton = new QPushButton("Show Templates");
    toggleTemplatesButton->setCheckable(true);
    toggleTemplatesButton->setChecked(do_show_templates);
    toggleTemplatesButton->setMinimumHeight(40);
    connect(toggleTemplatesButton, SIGNAL(clicked()), this, SLOT(ToggleDrawTemplates()));

    guidesWidgetLayout->addRow(toggleTemplatesButton);


    guidesWidget = new QTabWidget();
    guidesWidget->setMinimumWidth(225);
    guidesWidget->setLayout(guidesWidgetLayout);

    return guidesWidget;
}

QWidget * GLWidget::GetPhysicsWidget()
{
    if (physicsWidget != NULL) {
        return physicsWidget;
    }   

    QFormLayout * physWidgetLayout = new QFormLayout();

    testPhysicsButton = new QPushButton("Test Physics");
    testPhysicsButton->setCheckable(true);
    testPhysicsButton->setMinimumHeight(40);
    connect(testPhysicsButton, SIGNAL(clicked()), this, SLOT(ToggleDoPhysicsTest()) );
    physWidgetLayout->addRow(testPhysicsButton);


    QPushButton * weightButton1 = new QPushButton("Add Weight");
    connect(weightButton1, SIGNAL(clicked()), this, SLOT(DoPhysicsAddWeight()));

    QPushButton * weightButton2 = new QPushButton("Remove All");
    connect(weightButton2, SIGNAL(clicked()), this, SLOT(DoPhysicsRemoveWeights()));
    weightButton2->setStyleSheet("QPushButton{ background-color: #40FA42; color:black}"); // for any delete button

    QDoubleSpinBox * new_weight_spinbox = new QDoubleSpinBox();
    new_weight_spinbox->setRange(0.001, 50);
    new_weight_spinbox->setDecimals(3);
    new_weight_spinbox->setSingleStep(0.1);
    new_weight_spinbox->setValue(physics_new_weight_mass);
    connect(new_weight_spinbox, SIGNAL(valueChanged(double)), this, SLOT(SetPhysicsNewWeightMass(double)));

    new_weight_label = new QLabel(QString::number(new_weight_spinbox->value()) + QString(" kilograms"));

    QGroupBox * weight_groupbox = new QGroupBox(tr("Weights"));
    QGridLayout *weight_layout = new QGridLayout;
    weight_layout->addWidget(weightButton1, 0, 0, 1, 2);
    weight_layout->addWidget(weightButton2, 1, 0, 1, 2);
    weight_layout->addWidget(new QLabel("Weight (kg)"),2,0);
    weight_layout->addWidget(new_weight_spinbox,2,1);
    weight_groupbox->setLayout(weight_layout);
    physWidgetLayout->addRow(weight_groupbox);


    //display checkboxes

    deformed_checkbox = new QCheckBox("Deformed");
    deformed_checkbox->setChecked(physics.GetDrawDeformed());
    connect(deformed_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDrawDeformed(bool)));

    skeleton_checkbox = new QCheckBox("Skeleton");
    skeleton_checkbox->setChecked(physics.GetDrawSkeleton());
    connect(skeleton_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDrawSkeleton(bool)));

    forces_checkbox = new QCheckBox("Forces");
    forces_checkbox->setChecked(physics.GetDrawForce());
    connect(forces_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDrawForce(bool)));

    section_checkbox = new QCheckBox("Section");
    section_checkbox->setChecked(physics.GetDrawSection());
    connect(section_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDrawSection(bool)));

    moment_checkbox = new QCheckBox("Section Moments");
    moment_checkbox->setChecked(physics.GetDrawSectionMoment());
    connect(moment_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDrawSectionMoment(bool)));



    QGroupBox * show_groupbox = new QGroupBox(tr("Display"));
    QGridLayout *show_layout = new QGridLayout;

    show_layout->addWidget(forces_checkbox, 0, 0);
    show_layout->addWidget(deformed_checkbox, 0, 1);
    show_layout->addWidget(skeleton_checkbox, 1, 0);
    show_layout->addWidget(section_checkbox, 1, 1);
    show_layout->addWidget(moment_checkbox, 2, 0, 1, 2);
    show_groupbox->setLayout(show_layout);
    physWidgetLayout->addRow(show_groupbox);


    QDoubleSpinBox * density_spinbox = new QDoubleSpinBox();
    density_spinbox->setRange(10, 2000);
    density_spinbox->setDecimals(1);
    density_spinbox->setSingleStep(10);
    density_spinbox->setValue(physics_material_density);
    connect(density_spinbox, SIGNAL(valueChanged(double)), this, SLOT(SetPhysicsMaterialDensity(double)));

    QDoubleSpinBox * stress_spinbox = new QDoubleSpinBox();
    stress_spinbox->setRange(1, 500);
    stress_spinbox->setDecimals(2);
    stress_spinbox->setSingleStep(1);
    stress_spinbox->setValue(physics_max_stress / 1000000.0);
    connect(stress_spinbox, SIGNAL(valueChanged(double)), this, SLOT(SetPhysicsMaximumStress(double)));

    QGroupBox * phys_groupbox = new QGroupBox(tr("Physical Values"));
    QGridLayout *phys_layout = new QGridLayout;
    phys_layout->addWidget(new QLabel("Density (kg/m^3)"),0,0,1,2);
    phys_layout->addWidget(density_spinbox,0,2);
    phys_layout->addWidget(new QLabel("Stress (MPa)"),1,0,1,2);
    phys_layout->addWidget(stress_spinbox,1,2);
    phys_groupbox->setLayout(phys_layout);
    physWidgetLayout->addRow(phys_groupbox);

    physicsWidget = new QTabWidget();
    physicsWidget->setMinimumWidth(225);
    physicsWidget->setLayout(physWidgetLayout);



    return physicsWidget;
}

QWidget * GLWidget::GetViewWidget()
{
    if (viewWidget != NULL) {
        return viewWidget;
    }

    QFormLayout * viewWidgetLayout = new QFormLayout();

    //view
    const int viewbtnwidth = 50;
    QPushButton * viewisobtn1 = new QPushButton("Iso1");
    viewisobtn1->setMaximumWidth(viewbtnwidth);
    connect(viewisobtn1, SIGNAL(clicked()), this, SLOT(SetViewIso1()));
    QPushButton * viewisobtn2 = new QPushButton("Iso2");
    viewisobtn2->setMaximumWidth(viewbtnwidth);
    connect(viewisobtn2, SIGNAL(clicked()), this, SLOT(SetViewIso2()));
    QPushButton * viewisobtn3 = new QPushButton("Iso3");
    viewisobtn3->setMaximumWidth(viewbtnwidth);
    connect(viewisobtn3, SIGNAL(clicked()), this, SLOT(SetViewIso3()));
    QPushButton * viewisobtn4 = new QPushButton("Iso4");
    viewisobtn4->setMaximumWidth(viewbtnwidth);
    connect(viewisobtn4, SIGNAL(clicked()), this, SLOT(SetViewIso4()));

    QPushButton * viewxbtn = new QPushButton("X");
    viewxbtn->setMaximumWidth(viewbtnwidth);
    connect(viewxbtn, SIGNAL(clicked()), this, SLOT(SetViewX()));
    QPushButton * viewybtn = new QPushButton("Y");
    viewybtn->setMaximumWidth(viewbtnwidth);
    connect(viewybtn, SIGNAL(clicked()), this, SLOT(SetViewY()));
    QPushButton * viewzbtn = new QPushButton("Z");
    viewzbtn->setMaximumWidth(viewbtnwidth);
    connect(viewzbtn, SIGNAL(clicked()), this, SLOT(SetViewZ()));
    QPushButton * viewpartbtn = new QPushButton("Part");
    viewpartbtn->setMaximumWidth(viewbtnwidth);
    connect(viewpartbtn, SIGNAL(clicked()), this, SLOT(SetViewPart()));


    QGroupBox * view_groupbox = new QGroupBox(tr("Views"));
    QGridLayout *view_layout = new QGridLayout;
    view_layout->setSpacing(0);
    view_layout->setHorizontalSpacing(0);
    view_layout->addWidget(viewisobtn1, 0, 0);
    view_layout->addWidget(viewisobtn2, 0, 1);
    view_layout->addWidget(viewisobtn3, 0, 2);
    view_layout->addWidget(viewisobtn4, 0, 3);
    view_layout->addWidget(viewxbtn, 1, 0);
    view_layout->addWidget(viewybtn, 1, 1);
    view_layout->addWidget(viewzbtn, 1, 2);
    view_layout->addWidget(viewpartbtn, 1, 3);
    view_groupbox->setLayout(view_layout);
    viewWidgetLayout->addRow(view_groupbox);



    QSlider * rotation_slider = new QSlider();
    rotation_slider->setRange(1, 30);
    rotation_slider->setOrientation(Qt::Horizontal);
    rotation_slider->setValue(cam_animate_duration / 100.0f);
    connect(rotation_slider, SIGNAL(valueChanged(int)), this, SLOT(SetRotationDuration(int)));

    QSlider * rotation_angle_slider = new QSlider();
    rotation_angle_slider->setRange(0, 90);
    rotation_angle_slider->setOrientation(Qt::Horizontal);
    rotation_angle_slider->setValue(cam_animate_angle);
    connect(rotation_angle_slider, SIGNAL(valueChanged(int)), this, SLOT(SetRotationAngle(int)));

    rotate_label = new QLabel(QString::number(rotation_slider->value()*100) + QString(" ms"));
    rotate_angle_label = new QLabel(QString::number(rotation_angle_slider->value()) + QString(" degrees"));

    QGroupBox * viewrot_groupbox = new QGroupBox(tr("View Rotation"));
    QGridLayout *viewrot_layout = new QGridLayout;
    viewrot_layout->addWidget(rotation_slider, 0, 1);
    viewrot_layout->addWidget(rotate_label, 0, 0);
    viewrot_layout->addWidget(rotation_angle_slider, 1, 1);
    viewrot_layout->addWidget(rotate_angle_label, 1, 0);
    viewrot_groupbox->setLayout(viewrot_layout);
    viewWidgetLayout->addRow(viewrot_groupbox);




    //display checkboxes

    show_tnb_frames_checkbox = new QCheckBox("TNB Frame");
    show_tnb_frames_checkbox->setChecked(do_show_tnb_frames);
    connect(show_tnb_frames_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetShowTNBFrames(bool)));

    show_cycles_test_checkbox = new QCheckBox("Cycles");
    show_cycles_test_checkbox->setChecked(do_cycles_test);
    connect(show_cycles_test_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDoCyclesTest(bool)));

    show_stability_checkbox = new QCheckBox("Stability");
    show_stability_checkbox->setChecked(do_stability_test);
    connect(show_stability_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDoStabilityTest(bool)));

    show_shadow_checkbox = new QCheckBox("Shadow");
    show_shadow_checkbox->setChecked(do_show_shadow);
    connect(show_shadow_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetShowShadow(bool)));

    show_connectivity_checkbox = new QCheckBox("Connectivity");
    show_connectivity_checkbox->setChecked(do_connected_test);
    connect(show_connectivity_checkbox, SIGNAL(toggled(bool)), this, SLOT(SetDoConnectedTest(bool)));

    QGroupBox * show_groupbox1 = new QGroupBox(tr("General Display"));
    QGridLayout *show_layout1 = new QGridLayout;

    show_layout1->addWidget(show_tnb_frames_checkbox, 0, 0);
    show_layout1->addWidget(show_shadow_checkbox, 0, 1);

    show_groupbox1->setLayout(show_layout1);
    viewWidgetLayout->addRow(show_groupbox1);

    QGroupBox * show_groupbox2 = new QGroupBox(tr("Display Tests"));
    QGridLayout *show_layout2 = new QGridLayout;

    show_layout2->addWidget(show_stability_checkbox, 0, 0);
    show_layout2->addWidget(show_cycles_test_checkbox, 0, 1);
    show_layout2->addWidget(show_connectivity_checkbox, 1, 0, 1, 2);

    show_groupbox2->setLayout(show_layout2);
    viewWidgetLayout->addRow(show_groupbox2);

    viewWidget = new QTabWidget();
    viewWidget->setMinimumWidth(225);
    viewWidget->setLayout(viewWidgetLayout);


    return viewWidget;
}


// This is needed to allow QPainter to work with standard GL functions
void GLWidget::paintEvent(QPaintEvent *)
{

    //makeCurrent(); //not needed - paintGL does this implicitly
    paintGL();

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    DrawInstructions(painter);
    painter.end();

}

void GLWidget::paintGL()
{         

    //keep animating if still in duration   
    //qDebug() << "GLWidget::paintGL()";

    deadzone_radius = 0.025f * cam.CamWidth();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    UpdateCamera();   

    GLfloat light_position[4];
    light_position[0]=cam.Eye().x();
    light_position[1]=cam.Eye().y();
    light_position[2]=cam.Eye().z();
    light_position[3]=1.0f;
    glLightfv (GL_LIGHT0, GL_POSITION, light_position);

    //this stuff packs nicely into a display list
    if (update_sections_disp_list) {
        if (sections_disp_list > 0) {
            glDeleteLists(sections_disp_list, 1);
            sections_disp_list = 0;
        }
        update_sections_disp_list = false;
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glColor3f(0.75f, 0.75f, 0.75f);

    DrawGroundPlane();

    if (IsSectionSelected() && do_show_shadow) {
        sections[selected].DrawShadow();
    }

    if (sections_disp_list > 0) {
        glCallList(sections_disp_list);
    }
    else {

        sections_disp_list = glGenLists(1);
        glNewList(sections_disp_list, GL_COMPILE_AND_EXECUTE);

        glEnable(GL_DEPTH_TEST);
        //DrawGroundPlane();

        if (do_show_shadow) {
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            glColor3f(0.75f, 0.75f, 0.75f);
            for (int i=0; i<sections.size(); ++i) {

/*                if(current_tool_state == TOOLSTATE_GENERATE)
                {
                    if ( !( i == selected || (generate_selections.size() > 0 && i == generate_selections[0]) ||
                            (generate_selections.size() > 1 && i == generate_selections[1]) ||
                           (generate_selections.size() == 3 && i == generate_selections[2]) ) ) {
                        sections[i].DrawShadow();
                    }
                }
                else */if (i != selected) {
                    sections[i].DrawShadow();
                }

            }
            glDisable(GL_CULL_FACE);
            glEnable(GL_DEPTH_TEST);
        }

        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);

        for (int i=0; i<sections.size(); ++i) {

            if(current_tool_state == TOOLSTATE_GENERATE)
            {
                if ( !( i == selected || (generate_selections.size() > 0 && i == generate_selections[0]) ||
                        (generate_selections.size() > 1 && i == generate_selections[1]) ||
                       (generate_selections.size() == 3 && i == generate_selections[2]) ) ) {
                    DrawSection(i);
                }
            }
            else if (i != selected) {
                DrawSection(i);
            }

        }

        glEndList();

    }

    //debug: show parameters of camera model
    //glColor3f(1, 0, 0);
    //GLutils::DrawArrow(QVector3D(0, 0, 0), cam.GetRightVector());
    //glColor3f(0, 1, 0);
    //GLutils::DrawArrow(QVector3D(0, 0, 0), cam.Up());
    //glColor3f(0, 0, 1);
    //GLutils::DrawArrow(QVector3D(0, 0, 0), cam.ViewDir());



    //this stuff changes, not for display list
    if (state == STATE_CURVE || state == STATE_PEN_POINT || state == STATE_PEN_DRAG) {
        glColor3f(0.25f, 0.60f, 0.25f);

        if (do_local_symmetry) {
            active_section_symmetry.DrawTris();
        }
        else {
            active_section.DrawTris();
        }

        //active_section.DrawInputPolyline();
    }


    if(current_tool_state == TOOLSTATE_GENERATE)
    {

        if(generate_selections.size() > 0 && gen_state != GENSTATE_GRID)
        {
            glColor3f(1.0f, 0.0f, 0.0f);
            sections[generate_selections[0]].DrawSlab();
            if(generate_selections.size() > 1)
            {
                glColor3f(0.0f, 1.0f, 0.0f);
                sections[generate_selections[1]].DrawSlab();
                if(generate_selections.size() > 2)
                {
                    glColor3f(0.0f, 0.0f, 1.0f);
                    sections[generate_selections[2]].DrawSlab();
                }
            }
        }


        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glColor4f(0.35f, 0.45f, 0.6f, 0.5f);
        for(int i = 0; i < generated_sections.size(); i++)
        {
            generated_sections[i].DrawSlab();
        }
        glDisable(GL_BLEND);

    }
    else if (IsSectionSelected()) {
        DrawSection(selected);
    }


    if (do_show_templates) {
        DrawTemplateSurface();
    }

    contour_graph.Draw();

    glDisable(GL_LIGHTING);
    if (do_stability_test) {
        PlanarSection::DrawConvexHull(sections_convex_hull, centre_of_mass_in_hull);
    }

    DrawSlot(slot_start, slot_end);

    if (recursive_slot_start != QVector3D(0, 0 ,0)) {
        DrawSlot(recursive_slot_start, recursive_slot_end);
    }

    //draw template cut only when:
    //1.  no sections exist,
    //2.  we are adding a new planar section,
    //3.  we are editing a planar section
    if (do_show_templates) {
        if (state == STATE_SLOT || //2.
                state == STATE_CAM_TRANSLATE ||
                state == STATE_DEADZONE ||
                state == STATE_CURVE ||
                state == STATE_PEN_POINT ||
                state == STATE_PEN_DRAG ||
                state == STATE_MANIP_CTRLPOINT || //3.
                sections.empty()) { //1.

            UpdateTemplateCut();
            DrawTemplateCut(template_cut);

        }
    }

    glEnable(GL_LIGHTING);

    if (do_physics_test) {

        const double scale_fac = 1.0 / metres_per_unit;

        //draw physics thing
        glPushMatrix();
        glScaled(scale_fac, scale_fac, scale_fac);
        physics.DrawEnhancedGL(scale_fac);
        glPopMatrix();
        /*
        glPushMatrix();
        glScaled(scale_fac, scale_fac, scale_fac);
        physics.DrawGL();
        glPopMatrix();
        */
        //physics_solved.DrawGL();
    }

    glClear(GL_DEPTH_BUFFER_BIT);
    UpdateMarkers();
    DrawMarkers();


    if (state == STATE_DIMENSIONING_SECOND) {
        DrawDimensionTool();
    }

    /*
    glDisable(GL_LIGHTING);
    if (do_physics_test) {
        //draw physics thing
        physics.DrawGL();
        //physics_solved.DrawGL();
    }
    */


    //draw the transform widget
    if (IsSectionSelected() && current_tool_state == TOOLSTATE_TRANSFORMING) {
        transform_widget.SetP(sections[selected].P());
        transform_widget.SetX(sections[selected].T());
        transform_widget.SetY(sections[selected].N());
        transform_widget.SetZ(sections[selected].B());
        transform_widget.DrawGL(cam.CamWidth(),height());
    }

    //draw TNB frame
    if (do_show_tnb_frames) {
        //if (IsSectionSelected()) {
        //    sections[selected].DrawTNBFrame();
        //}
        //draw all TNB frames
        for (int i=0; i<sections.size(); ++i) {
            sections[i].DrawTNBFrame();
        }
    }

    glDisable(GL_LIGHTING);

    //draw the sketch line
    if (state == STATE_DEADZONE) {
        active_section.DrawDeadzone(slot_start, slot_end, deadzone_radius);
    }
    else if (state == STATE_ADD_HOLE) {

        if (IsSectionSelected()) {

            glLineWidth(2.0f);
            glColor3f(0, 0, 0.3f);
            sections[selected].DrawCurve();
            glLineWidth(1.0f);

        }

    }
    else if (state == STATE_CURVE || state == STATE_PEN_POINT || state == STATE_PEN_DRAG) {

        //glColor3f(0.3, 0.3, 0.9);
        //sections.last().DrawSketch();                
        glLineWidth(2.0f);
        //glColor3f(0, 0, 0);
        glColor3f(0, 0, 0.3f);

        if (do_local_symmetry) {
            active_section_symmetry.DrawCurve();
        }
        else {
            active_section.DrawCurve();
        }
        //active_section.DrawSketch();
        //active_section.DrawCurveControlPoints();

        if (do_show_templates) {
            DrawTemplateCut(template_cut);
        }

        //draw control points/polygon
        //glColor3f(1.0, 0.4, 1.0);
        //sections.last().DrawCurveControlPoints();

        //glEnable(GL_LINE_STIPPLE);
        //glLineStipple(1, 0x00ff);
        //glColor3f(0.8, 0.8, 0.8);
        //sections.last().DrawCurveControlPolygon();
        //glDisable(GL_LINE_STIPPLE);

        glColor3f(0, 0, 0);
        glLineWidth(3.0f);
        active_section.DrawXZPlaneLine();
        glLineWidth(1.0f);

    }
    else if (state == STATE_CAM_TRANSLATE) {

    }
    else if (IsSectionSelected()) {

        //if (state == STATE_NONE || state == STATE_MANIP_CTRLPOINT) {
        if (current_tool_state == TOOLSTATE_DEFAULT && state != STATE_CURVE && state != STATE_RESKETCH_CURVE && state != STATE_CAM_TRANSLATE && state != STATE_DEADZONE && state != STATE_PEN_POINT && state != STATE_PEN_DRAG) {
//            glColor3f(1.0, 0.4, 1.0);
            sections[selected].DrawCurveControlPointsHandleStyle(cam.CamWidth(), cam.Eye());

//            glEnable(GL_LINE_STIPPLE);
//            glLineStipple(1, 0x00ff);
//            glColor3f(0.8, 0.8, 0.8);
//            sections[selected].DrawCurveControlPolygon();
//            glDisable(GL_LINE_STIPPLE);
        }

        //glColor3f(0.3, 0.3, 0.9);
        //sections.last().DrawSketch();

    }

    if (do_show_templates) {
        DrawTemplateImage();
    }

    if(state == STATE_PEN_POINT || state == STATE_PEN_DRAG)
        active_section.DrawCurveControlPointsHandleStyle(cam.CamWidth(), cam.Eye());


    // grab control point 2D positions

    QVector3D slot_start_projected, slot_end_projected;
    if(do_local_symmetry && (state == STATE_CURVE || state == STATE_RESKETCH_CURVE || state == STATE_DEADZONE || state == STATE_PEN_POINT || state == STATE_PEN_DRAG))
    {
        slot_start_projected = GLutils::ProjectPoint(slot_start);
        slot_end_projected = GLutils::ProjectPoint(slot_end);
    }


    cam.DrawGL_Ortho();

    // Draw symmetry line when using local symmetry with no planar sections
    if(do_local_symmetry)
    {
        if(sections.empty())
        {
            glEnable(GL_BLEND);
            glLineWidth(3);
            glColor4f(0.3f, 0.5f, 0.8f, 0.5f);
            glBegin(GL_LINES);
                glVertex2f(width()/2.0, 0);
                glVertex2f(width()/2.0, height());
            glEnd();
            glLineWidth(1);
        }
        else if(state == STATE_CURVE || state == STATE_RESKETCH_CURVE || state == STATE_DEADZONE || state == STATE_PEN_POINT || state == STATE_PEN_DRAG)
        {
            QVector3D dir = (slot_start_projected - slot_end_projected).normalized();
            float factor = ( 20.0 * 600.0/height() ) / cam.CamWidth();
            glEnable(GL_BLEND);
            glLineWidth(3);
            glColor4f(0.3f, 0.5f, 0.8f, 0.5f);
            glBegin(GL_LINES);
                glVertex2f(slot_start_projected.x() + factor*dir.x(), slot_start_projected.y() + factor*dir.y());
                glVertex2f(slot_start_projected.x() + 10000*dir.x(), slot_start_projected.y() + 10000*dir.y());
                glVertex2f(slot_end_projected.x() - factor*dir.x(), slot_end_projected.y() - factor*dir.y());
                glVertex2f(slot_end_projected.x() - 10000*dir.x(), slot_end_projected.y() - 10000*dir.y());
            glEnd();
            glLineWidth(1);

        }
    }

//    DrawInfo();

}



void GLWidget::DrawInstructions(QPainter & painter)
{
    if (current_tool_state != TOOLSTATE_GENERATE && state != STATE_PEN_POINT && state != STATE_PEN_DRAG) {
        return;
    }

    QString title;
    QString text;

    if (current_tool_state == TOOLSTATE_GENERATE) {

        switch (gen_state) {

        case GENSTATE_LINEAR:
            title = "Generate: Linear";
            text = "Instructions\n\n 1) Select the base section (right-click away from sections to reset selection)\n 2) Select a section that intersects the base\n 3) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_BLEND:
            title = "Generate: Blend";
            text = "Instructions\n\n 1) Select the base section (right-click away from sections to reset selection)\n 2) Select two sections that intersect the base\n 3) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_REVOLVE:
            title = "Generate: Revolve";
            text = "Instructions\n\n 1) Select the base section (right-click away from sections to reset selection)\n 2) Select a section that intersects the base\n 3) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_GRID:
            title = "Generate: Grid";
            text = "Instructions\n\n 1) Select a section to form the grid\n 2) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_SLICES:
            title = "Generate: Slices";
            text = "Instructions\n\n 1) A template obj must be loaded - This can be done in the Guides side panel\n 2) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_SURFACE_FACETS:
            title = "Generate: Surface Facets";
            text = "Instructions\n\n 1) A template obj must be loaded - This can be done in the Guides side panel\n 2) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_BRANCH: // branching isn't active yet
            title = "Generate: Branch";
            text = "Instructions\n\n 1) Set the root slot by clicking and dragging on a section\n 2) Add branches 3) Press Enter to accept or Escape to cancel";
            break;

        case GENSTATE_NUM:
            break;

        }

    }
    else if (state == STATE_PEN_POINT || state == STATE_PEN_DRAG) {
        title = "Pen Mode";
        text = "Left-click and drag to create a new point and handles\nPress Enter to accept or Escape to cancel";
    }

    QFontMetrics metrics = QFontMetrics(font());

    int border = qMax(10, metrics.leading());

    QRect rect = metrics.boundingRect(0, 0, width() - 2*border, int(height()*0.125),
                                      Qt::AlignLeft | Qt::TextWordWrap, text);
    rect.setHeight(rect.height() + 20); // This accounts for the title

    painter.setRenderHint(QPainter::TextAntialiasing);
    painter.fillRect(QRect(0, 0, width(), rect.height() + 2*border),
                     QColor(245, 245, 245, 127));


    QFont title_font = font();
    title_font.setPixelSize(14);
    painter.setFont(title_font);
    painter.setPen(QColor(100,100,100));
    painter.drawText(border, border,
                      rect.width(), rect.height(),
                      Qt::AlignLeft | Qt::TextWordWrap, title);

    painter.setFont(font());
    painter.setPen(QColor(136,136,136));
    painter.drawText(border, border + 20,
                      rect.width(), rect.height(),
                      Qt::AlignLeft | Qt::TextWordWrap, text);

}

void GLWidget::resizeGL(int width, int height)
{

    glViewport(0, 0, width, height);

}

void GLWidget::mousePressEvent(QMouseEvent *event)
{

/*

"algorithm"

on click, choose a plane to draw onto...
a) when there are no planes, make plane params for frontoparallel plane, PROCEED TO STATE CURVE
b) when there are planes, attempt to select one
  i) if selection fails, do nothing
  ii) otherwise, set that initial point as anchor for slot line, and PROCEED TO STATE SLOT


*/

    const bool ctrl_held = ((event->modifiers() & Qt::ControlModifier) > 0);
    const bool shift_held = ((event->modifiers() & Qt::ShiftModifier) > 0);

    UpdateCamera();

    mouse_pos = QVector2D(event->x(), height() - event->y()); 
    mouse_pos *= devicePixelRatio();

    if (event->button() == Qt::LeftButton) {

        if (!(pen_mode && state == STATE_PEN_POINT) && (ctrl_held || shift_held)) {
            state = STATE_ORBIT;
        }
        else {

            switch (state) {

            case STATE_TRANSFORM_WIDGET:

                if (IsSectionSelected()){

                    original_section = sections[selected];

                    QVector3D anchor_3d;
                    sections[selected].MouseRayIntersect(mouse_pos, anchor_3d);
    //                anchor_point = sections[selected].GetPoint2D(anchor_3d);
                    anchor_point = mouse_pos;

                    //figure out which widget element was clicked
                    const int index = PickTransformWidgetElement(mouse_pos);

                    if (index > 0) {
                        transform_widget.SetState(TransformWidgetState(index));
                        state = STATE_TRANSFORM_WIDGET;
                    }
                    else if(PickSection(mouse_pos, false) >= 0 ) {
                        transform_widget.SetState(NONE);
                        state = STATE_NONE;
                        mousePressEvent(event);
                    }
                    else{
                        state = STATE_ORBIT;
                    }


                }
                else {
                    state = STATE_NONE;
                }

                break;

            case STATE_DEADZONE:
                break;

            case STATE_PEN_POINT:
                active_section.AddCtrlPointPenPress(0, mouse_pos);
                active_section.SelectCtrlPoint(active_section.GetCurve(0).Points().size() - 3); // the last two points in the curve attach to the first point
                break;

            case STATE_RESKETCH_CURVE:
                break;

            case STATE_ADD_HOLE:


                if (IsSectionSelected()){
                    sections[selected].AddNewCurve();
                }

                break;

            default:

                if (sections.empty()) {

                    active_section = PlanarSection();
                    SetupPlanarSection(active_section);

                    active_section.SetP(cam.LookAt());
                    active_section.SetN(-cam.ViewDir());
                    active_section.SetT(cam.GetRightVector());
                    active_section.SetB(cam.Up());

                    active_section.SketchSetEditing(true);

                    // Make a pen point change the state
                    if (pen_mode) {
                        state = STATE_PEN_POINT;
                        active_section.AddCtrlPointPenPress(0, mouse_pos);
                        active_section.SelectCtrlPoint(active_section.GetCurve(0).Points().size() - 3); // the last two points in the curve attach to the first point
                    }
                    else {
                        state = STATE_CURVE;
                    }

                }
                else {

                    const int section_clicked = PickSection(mouse_pos, false);

                    if (section_clicked >= 0 && section_clicked < sections.size()) {

                        QVector3D p;
                        sections[section_clicked].MouseRayIntersect(mouse_pos, p);

                        // If selectiing sections for the generate tools
                        if (current_tool_state == TOOLSTATE_GENERATE && gen_state != GENSTATE_BRANCH) {

                            if (generate_selections.size() < selections_per_gen_type[gen_state] || gen_state == GENSTATE_GRID) {

                                if (!generate_selections.contains(section_clicked)) {

                                    if (gen_state == GENSTATE_GRID && generate_selections.size() == 1) {
                                        generate_selections[0] = section_clicked; // this allows the user to switch sections during grid generation
                                    }
                                    else {
                                        generate_selections.append(section_clicked);
                                    }

                                }

                            }

                            if (generate_selections.size() == selections_per_gen_type[gen_state]) {
                                ShowGenerate();
                            }

                            UpdateDraw();
                        }
                        else if (state == STATE_RECURSIVE_SETUP_SLOT) {
                            //modify properties of base section
                            recursive_slot_start = p;
                            recursive_slot_end = p;
                            state = STATE_RECURSIVE_SLOT;
                        }
                        else if (state == STATE_DIMENSIONING_FIRST) {
                            dimensiontool_start = p;
                            dimensiontool_end = p;
                            state = STATE_DIMENSIONING_SECOND;
                        }
                        else if(current_tool_state == TOOLSTATE_DEFAULT) {

                            active_section = PlanarSection();
                            SetupPlanarSection(active_section);
                            active_section.SetP(p);
                            active_section.SetN(sections[section_clicked].T());
                            active_section.SetT(sections[section_clicked].B());
                            active_section.SetB(sections[section_clicked].N());

                            slot_start = p;
                            slot_end = p;
                            state = STATE_SLOT;
                        }

                        if (current_tool_state != TOOLSTATE_GENERATE && section_clicked != selected) {
                            SetSelected(section_clicked);
                            UpdateDraw();

                        }

                    }
                    else {
                        state = STATE_ORBIT;
                    }

                }

                break;

            }

        }     

    }
    else if (event->button() == Qt::RightButton) {           

        //if something already selected, maybe we are trying to move it
        if (state == STATE_TRANSFORM_WIDGET) {
            transform_widget.SetState(NONE);
            state = STATE_NONE;
        }

        // If selectiing sections for the generate tools
        if (current_tool_state == TOOLSTATE_GENERATE) {

            const int new_select = PickSection(mouse_pos, false);

            //bugfix: only add the selection index to the list if there is an actual selection (new_select >= 0)
            if (new_select < 0) {
                generate_selections.clear();
                ShowGenerate();
            }
            else if (generate_selections.size() < selections_per_gen_type[gen_state]  || gen_state == GENSTATE_GRID) {

                if (!generate_selections.contains(new_select)) {
                    if (gen_state == GENSTATE_GRID && generate_selections.size() == 1) {
                        generate_selections[0] = new_select; // this allows the user to switch sections during grid generation
                    }
                    else {
                        generate_selections.append(new_select);
                    }
                }

            }

            if (generate_selections.size() == selections_per_gen_type[gen_state]) {
                ShowGenerate();
            }

            UpdateDraw();
        }
        else if (state == STATE_PEN_POINT) {
            if (active_section.IsCtrlPointSelected()) {
                state = STATE_PEN_DRAG;
            }
        }
        else if (IsSectionSelected()) {

            if (sections[selected].IsCtrlPointSelected()) {
                state = STATE_MANIP_CTRLPOINT;
                AddToUndoList(OP_MANIP_CTRLPOINT);
            }
            else if (sections[selected].IsWeightSelected()) {
                state = STATE_MANIP_WEIGHT;
                AddToUndoList(OP_MANIP_WEIGHT);
            }
            else {
                const int new_select = PickSection(mouse_pos, false);   
                if (new_select != selected) {
                    last_op = OP_SELECTION;
                    SetSelected(new_select);             
                    UpdateDraw();
                }

            }

        }
        else {         
            const int new_select = PickSection(mouse_pos, false);
            if (new_select != selected) {
                last_op = OP_SELECTION;
                SetSelected(new_select);
                UpdateDraw();
            }

        }

    }

    if (state != STATE_ORBIT && state != STATE_TRANSFORM_WIDGET && current_tool_state == TOOLSTATE_TRANSFORMING) {
        current_tool_state = TOOLSTATE_DEFAULT;
    }

}

void GLWidget::mouseMoveEvent(QMouseEvent * event)
{   

    const bool ctrl_held = ((event->modifiers() & Qt::ControlModifier) > 0);
    const bool shift_held = ((event->modifiers() & Qt::ShiftModifier) > 0);

    QVector2D p(event->x(), height() - event->y());
    p *= devicePixelRatio();

    //qDebug() << "GLWidget::mouseMoveEvent() - state " << state;
    QVector2D mouse_diff = mouse_pos -p;
    mouse_pos = p;

    UpdateCamera();

    if (event->buttons() & Qt::LeftButton) {       

        switch (state) {

        case STATE_TRANSFORM_WIDGET:

            if (IsSectionSelected()) {

                QVector3D p = sections[selected].P();
                QVector3D t = sections[selected].T();
                QVector3D n = sections[selected].N();
                QVector3D b = sections[selected].B();
                QVector3D cursor_p;
                QVector3D initial_cursor_p;

                switch (transform_widget.GetState()) {

                case TRANS_X:
                {
                    sections[selected].MouseRayIntersect(mouse_pos, cursor_p);
                    sections[selected].MouseRayIntersect(anchor_point, initial_cursor_p);
                    const float dot_prod = QVector3D::dotProduct(t, (cursor_p - initial_cursor_p));
                    sections[selected].SetP(p + t * dot_prod);
                }
                    break;

                case TRANS_Y:
                {
                    PlanarSection ps;
                    ps.SetP(p);
                    ps.SetN(t);
                    ps.MouseRayIntersect(mouse_pos, cursor_p);
                    ps.MouseRayIntersect(anchor_point, initial_cursor_p);
                    float dot_prod = QVector3D::dotProduct(n, (cursor_p - initial_cursor_p));
                    sections[selected].SetP(p + n * dot_prod);
                }
                    break;

                case TRANS_Z:
                {
                    sections[selected].MouseRayIntersect(mouse_pos, cursor_p);
                    sections[selected].MouseRayIntersect(anchor_point, initial_cursor_p);
                    const float dot_prod = QVector3D::dotProduct(b, (cursor_p - initial_cursor_p));
                    sections[selected].SetP(p + b * dot_prod);
                }
                    break;

                case ROT_X:
                {
                    PlanarSection ps;
                    ps.SetP(p);
                    ps.SetN(t);
                    ps.MouseRayIntersect(mouse_pos, cursor_p);
                    ps.MouseRayIntersect(anchor_point, initial_cursor_p);

                    QVector3D v1 = (initial_cursor_p - p).normalized();
                    QVector3D v2 = (cursor_p - p).normalized();
                    QVector3D new_n = GLutils::RotateVector(n,t,GLutils::SignedAngleBetweenRad(v1,v2,t));

                    sections[selected].SetN(new_n);
                    sections[selected].SetB(QVector3D::crossProduct(new_n, t).normalized());
                }
                    break;
                case ROT_Y:
                {
                    sections[selected].MouseRayIntersect(mouse_pos, cursor_p);
                    sections[selected].MouseRayIntersect(anchor_point, initial_cursor_p);

                    QVector3D v1 = (initial_cursor_p - p).normalized();
                    QVector3D v2 = (cursor_p - p).normalized();
                    QVector3D new_b = GLutils::RotateVector(b,n,GLutils::SignedAngleBetweenRad(v1,v2,n));

                    sections[selected].SetB(new_b);
                    sections[selected].SetT(QVector3D::crossProduct(new_b, n).normalized());
                }
                    break;
                case ROT_Z:
                {
                    PlanarSection ps;
                    ps.SetP(p);
                    ps.SetN(b);
                    ps.MouseRayIntersect(mouse_pos, cursor_p);
                    ps.MouseRayIntersect(anchor_point, initial_cursor_p);

                    QVector3D v1 = (initial_cursor_p - p).normalized();
                    QVector3D v2 = (cursor_p - p).normalized();
                    QVector3D new_t = GLutils::RotateVector(t,b,GLutils::SignedAngleBetweenRad(v1,v2,b));

                    sections[selected].SetT(new_t);
                    sections[selected].SetN(QVector3D::crossProduct(new_t, b).normalized());
                }
                    break;

                case SCALE_X:
                {
                    sections[selected].MouseRayIntersect(mouse_pos, cursor_p);
                    sections[selected].MouseRayIntersect(anchor_point, initial_cursor_p);
                    const float dot_prod = QVector3D::dotProduct(t, (cursor_p - initial_cursor_p));
                    sections[selected].Scale(1 + .5*dot_prod, 1.0f);
                }
                    break;

                case SCALE_Z:
                {
                    sections[selected].MouseRayIntersect(mouse_pos, cursor_p);
                    sections[selected].MouseRayIntersect(anchor_point, initial_cursor_p);
                    const float dot_prod = QVector3D::dotProduct(b, (cursor_p - initial_cursor_p));
                    sections[selected].Scale(1.0f, 1 + .5*dot_prod);
                }
                    break;

                default:
                    break;
                }

                sections[selected].UpdateCurveTrisSlab();
                anchor_point = mouse_pos;

            }

            break;

        case STATE_RECURSIVE_SLOT:

            if (PickSection(mouse_pos, true) == selected) {

                QVector3D p;
                sections[selected].MouseRayIntersect(mouse_pos, p);
                recursive_slot_end = p;

                sections[selected].MoveP(recursive_slot_end);
                //sections[selected].SetP(recursive_slot_end);
                sections[selected].UpdateCurveTrisSlab();

                //qDebug() << sections.size();

            }
            else {

                state = STATE_NONE;

            }

            break;

        case STATE_DIMENSIONING_SECOND:

            if (PickSection(mouse_pos, true) == selected) {

                QVector3D p;
                sections[selected].MouseRayIntersect(mouse_pos, p);
                dimensiontool_end = p;

            }
            break;

        case STATE_SLOT:

            if (PickSection(mouse_pos, true) == selected) {

                QVector3D p;
                sections[selected].MouseRayIntersect(mouse_pos, p);
                slot_end = p;

                const QVector3D new_t = (-sections[selected].N()).normalized();
                const QVector3D new_b = (p - active_section.P()).normalized();
                const QVector3D new_n = QVector3D::crossProduct(new_t, new_b).normalized();

                active_section.SetT(new_t); //the dir
                active_section.SetN(new_n); //cross prod
                active_section.SetB(new_b);

                UpdateTemplateCut();

            }
            else {

                state = STATE_CAM_TRANSLATE;

                //"slot_end" is the mouse camera pivot point
                //we need to initiate the rotation,  but keep the vertex position slot_end constant in screen space...

                sections[selected].MouseRayIntersect(mouse_pos, slot_end);

                QVector3D v;
                sections[selected].MouseRayIntersect(QVector2D(float(this->width())/2.0f, float(this->height()/2.0f)), v);   

                QVector3D vs[4];
                vs[0] = GLutils::RotateVector(active_section.T(), active_section.B(), cam_animate_angle * 3.14159f / 180.0f);
                vs[1] = GLutils::MirrorVector(vs[0], active_section.N());
                vs[2] = GLutils::MirrorVector(vs[0], active_section.T());
                vs[3] = GLutils::MirrorVector(vs[2], active_section.N());

                float rots[4];
                rots[0] = GLutils::AngleBetweenDeg(cam.ViewDir(), vs[0]);
                rots[1] = GLutils::AngleBetweenDeg(cam.ViewDir(), vs[1]);
                rots[2] = GLutils::AngleBetweenDeg(cam.ViewDir(), vs[2]);
                rots[3] = GLutils::AngleBetweenDeg(cam.ViewDir(), vs[3]);

                int best_index = 0;
                float best_val = FLT_MAX; //best val is a MINIMUM dot product

                //among 4 choices, pick the one where of the downward viewing ones, the dot product with current view is lowest
                for (int i=0; i<4; ++i) {

                    if (vs[i].y() > 0.0f) {
                        continue;
                    }

                    if (rots[i] < best_val) {
                        best_index = i;
                        best_val = rots[i];
                    }

                }

                //visual debugging
                /*
                markers2.push_back(vs[0]);
                markers2.push_back(vs[1]);
                markers2.push_back(vs[2]);
                markers2.push_back(vs[3]);
                markers2.push_back(cam.ViewDir());
                markers_col2.push_back(QVector3D(1,0,0));
                markers_col2.push_back(QVector3D(1,1,0));
                markers_col2.push_back(QVector3D(0,1,0));
                markers_col2.push_back(QVector3D(0,1,1));
                markers_col2.push_back(QVector3D(0,0,1));
                */

                cam.StartPivot(slot_end, v, vs[best_index], cam_animate_duration);                

            }

            break;

        case STATE_CAM_TRANSLATE:           

            //if (cam.InterpActive()) {
            //    cam.UpdatePivotOffset(mouse_pos);
            //}
            //cam.Update(mouse_pos);

            break;

        case STATE_DEADZONE:

            if (!pen_mode && active_section.MouseOutsideDeadzone(mouse_pos, slot_start, slot_end, deadzone_radius)) {
                state = STATE_CURVE;
                active_section.SketchSetEditing(true);
            }

            break;

        case STATE_RESKETCH_CURVE:

            if (IsSectionSelected()) {

                sections[selected].SetLocalSymmetryMode(do_local_symmetry);
                sections[selected].AddMouseRayIntersect(0, mouse_pos);

                //if we have a template, do snapping
                if (do_show_templates && !template_cut.empty()) {
                    sections[selected].SnapSketchToTemplateCut(template_cut, template_cut_snap_distance_3d);
                    sections[selected].Update(0, section_error_tolerance_template);
                }
                else {
                    sections[selected].Update(0, section_error_tolerance);
                }
            }

            break;

        case STATE_ADD_HOLE:

            if (IsSectionSelected()) {

                sections[selected].AddMouseRayIntersect(sections[selected].GetNumCurves()-1, mouse_pos);
                sections[selected].Update(sections[selected].GetNumCurves()-1, section_error_tolerance);

            }

            break;

        case STATE_PEN_POINT:

            active_section.SetLocalSymmetryMode(do_local_symmetry);
            active_section.MoveCtrlPointMouseRayIntersect(mouse_pos, !ctrl_held, true);            
            active_section.UpdateCurveTrisSlab();

            if (do_local_symmetry) {
                active_section_symmetry = active_section;
                active_section_symmetry.MirrorControlPoints();
                active_section_symmetry.UpdateCurveTrisSlab();
            }

            break;

        case STATE_CURVE:

            active_section.SetLocalSymmetryMode(do_local_symmetry);
            active_section.AddMouseRayIntersect(0, mouse_pos);            

            //if we have a template, do snapping
            if (do_show_templates && !template_cut.empty()) {
                active_section.SnapSketchToTemplateCut(template_cut, template_cut_snap_distance_3d);
                active_section.Update(0, section_error_tolerance_template);
            }
            else {
                active_section.Update(0, section_error_tolerance);
            }

            if (do_local_symmetry) {
                active_section_symmetry = active_section;

                active_section_symmetry.SketchSymmetryTest();
                active_section_symmetry.Update(0, section_error_tolerance);
                active_section_symmetry.CreateLocalSymmetry();

                active_section_symmetry.UpdateCurveTrisSlab();

            }

            break;

        case STATE_ORBIT:

            if (ctrl_held) {
                cam.Zoom(mouse_diff);                                
            }
            else if (shift_held) {
                cam.Dolly(mouse_diff);
            }
            else {
                cam.Orbit(mouse_diff);
            }

            break;

        default:
            break;

        }

    }
    else if (event->buttons() & Qt::RightButton) {

        if (state == STATE_PEN_DRAG) {
            active_section.MoveCtrlPointMouseRayIntersect(mouse_pos, !ctrl_held);
            active_section.UpdateCurveTrisSlab();

            if (do_local_symmetry) {
                active_section_symmetry = active_section;
                active_section_symmetry.MirrorControlPoints();
                active_section_symmetry.UpdateCurveTrisSlab();
            }
        }
        else if (IsSectionSelected()) {

            if (state == STATE_MANIP_WEIGHT) {
                sections[selected].MoveWeightMouseRayIntersect(mouse_pos);
            }
            else if (state == STATE_MANIP_CTRLPOINT) {
                sections[selected].MoveCtrlPointMouseRayIntersect(mouse_pos, !ctrl_held);
                sections[selected].UpdateCurveTrisSlab();
            }

            if (do_physics_test) {
                DoPhysicsTest();
            }

        }


    }
    else {  

        if(state == STATE_PEN_POINT)
        {
            active_section.SelectMouseRayIntersect(mouse_pos, cam.CamWidth());
        }
        else if (IsSectionSelected()) {
            sections[selected].SelectMouseRayIntersect(mouse_pos, cam.CamWidth());
        }



    }

    //updateGL();

}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{

    mouse_pos = QVector2D(event->x(), height() - event->y());
    mouse_pos *= devicePixelRatio();

    UpdateCamera();

    if (event->button() == Qt::LeftButton) {

        switch (state) {

        case STATE_SLOT:
        case STATE_CAM_TRANSLATE:

            state = STATE_NONE;
            cam.SetInterpActive(false);            

            break;

        case STATE_DEADZONE:

            //if the gesture wasn't completed, remove the section being added
            //sections.removeLast();
            //active_section = null_section;
            if (pen_mode) {
                state = STATE_PEN_POINT;
            }
            else if (selected != -1) {
                SetSelected(-1);
                UpdateDraw();
            }
            break;

        case STATE_PEN_POINT:
            active_section.UnselectCtrlPoint();
            break;

        case STATE_RESKETCH_CURVE:

            if (sections.size() >= 2 && do_local_symmetry) {
                //do not do symmetry for the very first plane
                sections[selected].SketchSymmetryTest();
            }

            //at the curve step, the gesture is considered "completed"
            if (!template_cut.empty()) {
                sections[selected].Update(0, section_error_tolerance_template);
            }
            else {
                sections[selected].Update(0, section_error_tolerance);
            }

            if (sections.size() >= 2 && do_local_symmetry) {
                //do not do symmetry for the very first plane
                sections[selected].CreateLocalSymmetry();
            }

            sections[selected].SetCtrlPointsAboveXZPlane();
            sections[selected].UpdateCurveTrisSlab();

            break;

        case STATE_CURVE:

            if (active_section.SketchNumPoints() >= 5) {

                //add existing planar section set to undo list
                AddToUndoList(OP_ADD_PLANE);

                active_section.SketchSetEditing(false);
                if (do_local_symmetry) {
                    active_section_symmetry.SketchSetEditing(false);
                }

                //enforce above-ground control points
                active_section.SetCtrlPointsAboveXZPlane();
                active_section.UpdateCurveTrisSlab();

                if (do_local_symmetry) {
                    active_section_symmetry.SetCtrlPointsAboveXZPlane();
                    active_section_symmetry.UpdateCurveTrisSlab();
                }

                if (do_local_symmetry && !active_section_symmetry.SliceTriangles().empty()) {
                    sections.push_back(active_section_symmetry);
                    UpdateAllTests();
                }
                else if (!do_local_symmetry && !active_section.SliceTriangles().empty()) {
                    //for the added curve, we now update our physical tests (TODO: checkbox/boolean to disable this)
                    sections.push_back(active_section);                   
                    UpdateAllTests();
                }

            }

            SetSelected(-1);

            break;

        case STATE_DIMENSIONING_SECOND:

        {
            bool ok;
            double new_length = QInputDialog::getDouble(this, QString("Dimensioning Tool"), QString("Specify a length for the defined line segment:"), 5.0, 0.1, 100.0, 1, &ok, Qt::Dialog);
            if (ok) {

                //scale the model by the ratio of the defined length and the entered value
                AddToUndoList(OP_MANIP_DIMENSIONING_TOOL);

                const double scale_val = new_length / (dimensiontool_end - dimensiontool_start).length();
                this->ScalePlanarSections(scale_val);

            }
        }

            break;

        default:
            break;

        }

        if (state == STATE_ORBIT && current_tool_state == TOOLSTATE_TRANSFORMING) {
            state = STATE_TRANSFORM_WIDGET;
        }

        if (state != STATE_TRANSFORM_WIDGET && state != STATE_PEN_POINT) {
            state = STATE_NONE;
        }

        if (state != STATE_PEN_POINT) {
            slot_end = slot_start;
        }

    }
    else if (event->button() == Qt::RightButton) {


        if (state == STATE_PEN_DRAG)
        {
            state = STATE_PEN_POINT;
        }
        else if (IsSectionSelected()) {

            sections[selected].SetCtrlPointsAboveXZPlane();
            sections[selected].UpdateCurveTrisSlab();

            sections[selected].SelectMouseRayIntersect(mouse_pos, cam.CamWidth());

            state = STATE_NONE;

        }


        UpdateAllTests();

    }

    //UpdateDraw();
    //updateGL();

}

void GLWidget::keyPressEvent(QKeyEvent *event)
{

    switch (event->key()) {

    case Qt::Key_Enter:
    case Qt::Key_Return:

        if (state == STATE_PEN_POINT) {
            AcceptPenCurve();
        }
        else if (current_tool_state == TOOLSTATE_GENERATE) {
            AcceptGenerate();
        }

        break;

    case Qt::Key_Escape:

        if (current_tool_state == TOOLSTATE_GENERATE) {
            CancelGenerate();
        }
        else if (state == STATE_PEN_POINT) {
            CancelPenCurve();
        }

        break;
    }

}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{

    switch (event->key()) {

    case Qt::Key_Minus:

        if(state == STATE_PEN_POINT) {
            active_section.DeleteCtrlPoint(0, active_section.GetCurve(0).GetNumControlPoints() - 3);
            active_section.UpdateCurveTrisSlab();

            if (do_local_symmetry) {
                active_section_symmetry = active_section;
                active_section_symmetry.MirrorControlPoints();
                active_section_symmetry.UpdateCurveTrisSlab();
            }
        }
        else if (IsSectionSelected()) {


            if (sections[selected].IsCtrlPointSelected()) {
                AddToUndoList(OP_DELETE_CTRLPOINT);
                sections[selected].DeleteSelectedCtrlPoint();
            }
            else if (sections[selected].IsWeightSelected()) {
                AddToUndoList(OP_DELETE_WEIGHT);
                sections[selected].DeleteWeight();
            }
            sections[selected].UpdateCurveTrisSlab();

            UpdateAllTests();

        }
        break;

    case Qt::Key_Equal:
    case Qt::Key_Plus:

        if (IsSectionSelected()) {

            if (sections[selected].IsCtrlPointSelected()) {
                AddToUndoList(OP_ADD_CTRLPOINT);
                sections[selected].InsertCtrlPoint();
            }
            else if (do_physics_test) {
                AddToUndoList(OP_ADD_WEIGHT);
                UpdateCamera();
                sections[selected].AddWeightAtMousePos(mouse_pos, physics_new_weight_mass);
            }

            sections[selected].UpdateCurveTrisSlab();

            UpdateAllTests();

        }
        break;  

    }

    UpdateDraw();

}

void GLWidget::UpdateCamera() {

    makeCurrent();

    cam.SetWinWidth(width());
    cam.SetWinHeight(height());
    cam.Update(mouse_pos);
    cam.DrawGL_3DOrtho();

    if (state == STATE_CAM_TRANSLATE && !cam.InterpActive()) {
            state = STATE_DEADZONE;
    }

}

void GLWidget::DrawGroundPlane()
{

    glColor3f(0.95f, 0.95f, 0.95f);
    glBegin(GL_LINES);
        for (int i=-grid_size; i<=grid_size; ++i) {

            if(i != 0)
            {
                glVertex3i(i, 0, -grid_size);
                glVertex3i(i, 0, grid_size);

                glVertex3i(-grid_size, 0, i);
                glVertex3i(grid_size, 0, i);
            }

        }

        glColor3f(0.8f, 0.8f, 0.8f);
        glVertex3i(0, 0, -grid_size);
        glVertex3i(0, 0, grid_size);

        glVertex3i(-grid_size, 0, 0);
        glVertex3i(grid_size, 0, 0);

    glEnd();
    glColor3i(255, 255, 255);

}

void GLWidget::DrawSymmetryPlanes()
{

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glDepthMask(GL_FALSE);

    glColor4f(0.0f, 0.0f, 0.0f, 0.03f);

    if (x_symmetry) {

        for (int i=0; i<2; ++i) {

            i == 0 ? glBegin(GL_LINE_LOOP) : glBegin(GL_QUADS);
            glVertex3i(0, 0, -grid_size);
            glVertex3i(0, 0, grid_size);
            glVertex3i(0, grid_size, grid_size);
            glVertex3i(0, grid_size, -grid_size);
            glEnd();

        }

    }

    if (y_symmetry) {

        for (int i=0; i<2; ++i) {

            i == 0 ? glBegin(GL_LINE_LOOP) : glBegin(GL_QUADS);
            glVertex3i(-grid_size, 0, -grid_size);
            glVertex3i(-grid_size, 0, grid_size);
            glVertex3i(grid_size, 0, grid_size);
            glVertex3i(grid_size, 0, -grid_size);
            glEnd();

        }

    }

    if (z_symmetry) {

        for (int i=0; i<2; ++i) {

            i == 0 ? glBegin(GL_LINE_LOOP) : glBegin(GL_QUADS);
            glVertex3i(-grid_size, 0, 0);
            glVertex3i(-grid_size, grid_size, 0);
            glVertex3i(grid_size, grid_size, 0);
            glVertex3i(grid_size, 0, 0);
            glEnd();

        }

    }

    glColor3i(255, 255, 255);

    glDepthMask(GL_TRUE);

    glDisable(GL_BLEND);

}

void GLWidget::DrawSection(const int i)
{   

    if (i == selected) {
        if( state == STATE_RESKETCH_CURVE)
        {
            glEnable(GL_BLEND);
            glColor4f(0.35f, 0.45f, 0.6f, 0.5f);
        }
        else
            glColor3f(0.35f, 0.45f, 0.6f);
    }
    else {
        if (sections[i].PartOfCycle() || !sections[i].Connected()) {
            glColor3f(0.5f, 0.25f, 0.25f);
        }
        else {
            glColor3f(0.25f, 0.60f, 0.25f);
        }
    }

    if (sections[i].SketchEditing()) {
        sections[i].DrawTris();
    }
    else {
        sections[i].DrawSlab();
    }   

    glColor3f(0, 0, 0);
    if (sections[i].SketchEditing()) {
        sections[i].DrawCurve();
    }
    else {
        sections[i].DrawSlabCurves();
    }

    if(do_physics_test)
        sections[i].DrawWeights();

    if( state == STATE_RESKETCH_CURVE)
    {
        glDisable(GL_BLEND);
    }

}

void GLWidget::DoCyclesConnectedTest()
{   

    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetPartOfCycle(false);
        sections[i].SetConnected(true);
    }

    //these tests are put together since both rely on computation of
    //the same topology info - graphs and paths
    if (!do_cycles_test && !do_connected_test) {
        return;
    }

    QVector <QVector <bool> > graph;    
    QList <QList <int> > cycles;

    PlanarSection::ComputeIntersectionGraph(sections, graph);

    Tree tree;
    tree.CreateFromGraph(graph, 0);

    if (do_cycles_test && tree.HasCycles()) {       

        QList <QList <int> > cycles;
        tree.GetCycles(cycles);

        //qDebug() << cycles;
        //for each cycle, we test if there are any intersection lines which are (about) parallel,
        //which means the cycle is in fact assemblable
        for (int i=0; i<cycles.size(); ++i) {

            if (!PlanarSection::CycleAssemblable(sections, cycles[i])) {
                for (int j=0; j<cycles[i].size(); ++j) {
                    sections[cycles[i][j]].SetPartOfCycle(true);
                }
            }

        }
    }

    if (do_connected_test && !tree.IsConnected()) {

        QList <int> disc_nodes;
        tree.GetDisconnectedNodes(disc_nodes);

        for (int i=0; i<disc_nodes.size(); ++i) {
            sections[disc_nodes[i]].SetConnected(false);
        }
    }

    /*
    //old inefficient code
    //QList <QList <int> > all_paths;
    PlanarSection::ComputeCycles(graph, cycles, all_paths);
    //qDebug() << graph;
    //qDebug() << "Detected" << cycles.size() << "cycles";

    if (do_cycles_test) {
        for (int i=0; i<sections.size(); ++i) {
            sections[i].SetPartOfCycle(false);
        }
        for (int i=0; i<cycles.size(); ++i) {
            for (int j=0; j<cycles[i].size(); ++j) {
                sections[cycles[i][j]].SetPartOfCycle(true);
            }
        }
    }

    if (do_connected_test) {
        //test the connectedness
        PlanarSection::TestConnectedness(sections, all_paths);
    }
    */

}

void GLWidget::DoStabilityTest()
{

    if (!do_stability_test) {
        return;
    }

    centre_of_mass = PlanarSection::Centroid(sections);

    PlanarSection::GetXZConvexHull(sections, sections_convex_hull);
    centre_of_mass_in_hull = PlanarSection::GetXZPointInsideConvexHull(sections_convex_hull, centre_of_mass);

    //markers.push_back(centre_of_mass);
    //markers_col.push_back(QVector3D(0.75, 0.75, 0.75));
    markers.push_back(QVector3D(centre_of_mass.x(), 0, centre_of_mass.z()));
    markers_col.push_back(QVector3D(0.75, 0.75, 0.75));

    for (int i=0; i<sections_convex_hull.size(); ++i) {
        markers.push_back(sections_convex_hull[i]);
        if (centre_of_mass_in_hull) {
            markers_col.push_back(QVector3D(0.25, 1, 0.25));
        }
        else {
            markers_col.push_back(QVector3D(1, 0.25, 0.25));
        }
    }

}

int GLWidget::PickTransformWidgetElement(const QVector2D & mouse_pos)
{

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    transform_widget.DrawSelectionGL(cam.CamWidth(), height());

    unsigned char r, g, b;
    int index;

    GLutils::ReadPixelColor_BackBuffer(mouse_pos.x(), mouse_pos.y(), r, g, b);
    GLutils::PixelColorToIndex(r, g, b, index);

    if (index < NUM_STATES) {
        return index;
    }
    else {
        return -1;
    }

}

int GLWidget::PickSection(const QVector2D & mouse_pos, const bool only_test_selected)
{   

    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (only_test_selected && IsSectionSelected()) {
        //only draw the one plane we are testing, is the user mouse still on it?
        sections[selected].DrawForPicking(selected);
    }
    else {
        //draw all planes... which of them do we want?        
        for (int i=0; i<sections.size(); ++i) {
            sections[i].DrawForPicking(i);
        }
    }

    unsigned char r, g, b;
    int index;

    GLutils::ReadPixelColor_BackBuffer(mouse_pos.x(), mouse_pos.y(), r, g, b);   
    GLutils::PixelColorToIndex(r, g, b, index);   

    if (index >= 0 && index < sections.size()) {
        return index;        
    }
    else {
        return -1;
    }

}

void GLWidget::SetUnitsToInches()
{

    SetMetresPerUnit(0.0254);

}

void GLWidget::SetUnitsToCentimetres()
{

    SetMetresPerUnit(0.01);

}

void GLWidget::SetupPlanarSection(PlanarSection & p)
{
    p.SetSlabThickness(slab_thickness);
    p.SetQuality(quality_samples);
}

void GLWidget::SetMetresPerUnit(const double d)
{

    metres_per_unit = d;
    //qDebug() << "GLWidget::SetMetresPerUnit() " << metres_per_unit;

    UpdateAllTests();
    //updateGL();
}

void GLWidget::ScalePlanarSections(const float s){

    if (s <= 0.0f) {
        return;
    }   

    for (int i=0; i<sections.size(); ++i) {
        sections[i].Scale(s, s);
        sections[i].SetP(sections[i].P() * s);
        sections[i].UpdateCurveTrisSlab();
    }

    UpdateAllTests();
    //updateGL();

}

void GLWidget::SetViewIso1()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(1, 1, 1).normalized() * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewIso2()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(-1, 1, 1).normalized() * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewIso3()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(-1, 1, -1).normalized() * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewIso4()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(1, 1, -1).normalized() * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewX()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(1, 0, 0) * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();
}

void GLWidget::SetViewY()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(0, 1, 0) * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewZ()
{

    ResetCamera();

    cam.SetLookAt(default_lookat);
    cam.SetEye(default_lookat + QVector3D(0, 0, 1) * cam_lookat_distance);

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetViewPart()
{  

    if (!sections.empty() && IsSectionSelected()) {

        QVector2D min_v, max_v;
        sections[selected].GetBoundingBox2D(min_v, max_v);

        const float desired_cam_width = qMax(max_v.x() - min_v.x(), max_v.y() - min_v.y()) * 1.2f;
        const QVector2D centre = (min_v + max_v) * 0.5f;
        const QVector3D centre_3d = sections[selected].GetPoint3D(centre);

        cam.SetLookAt(centre_3d);
        if (sections[selected].N().y() > 0.0f) {
            cam.SetEye(centre_3d + sections[selected].N() * cam_lookat_distance);
        }
        else {
            cam.SetEye(centre_3d - sections[selected].N() * cam_lookat_distance);
        }
        cam.SetCamWidth(desired_cam_width);

    }

    UpdateTemplateCut();

    //updateGL();

}

void GLWidget::SetXSymmetry(bool b)
{

    x_symmetry = b;
    //updateGL();

}

void GLWidget::SetYSymmetry(bool b)
{
    y_symmetry = b;
    //updateGL();
}

void GLWidget::SetZSymmetry(bool b)
{
    z_symmetry = b;
    //updateGL();
}

void GLWidget::SetTemplateImageX(const int i)
{
    template_pos.setX(float(i)/20.0f);
    //updateGL();
}

void GLWidget::SetTemplateImageY(const int i)
{
    template_pos.setY(float(i)/20.0f);
    //updateGL();
}

void GLWidget::SetTemplateImageRotate(const int i)
{
    template_rotation = i;
    //updateGL();
}

void GLWidget::SetTemplateImageScale(const int i)
{
    template_scale = float(i)/20.0f;
    //updateGL();
}

void GLWidget::SetTemplateImageFlipX(const bool b)
{
    template_flipx = b;
    //updateGL();
}

bool GLWidget::GetDoMagneticCuts()
{
    return do_magnetic_cuts;
}

bool GLWidget::GetPenModeOn()
{
    return pen_mode;
}

bool GLWidget::GetDoLocalSymmetry()
{
    return do_local_symmetry;
}

bool GLWidget::GetDoCyclesTest()
{
    return do_cycles_test;
}

bool GLWidget::GetDoStabilityTest()
{
    return do_stability_test;
}

bool GLWidget::GetDoPhysicsTest()
{
    return do_physics_test;
}

bool GLWidget::GetDoConnectedTest()
{
    return do_connected_test;
}

bool GLWidget::GetShowTNBFrames()
{
    return do_show_tnb_frames;
}

bool GLWidget::GetShowShadow()
{
    return do_show_shadow;
}

bool GLWidget::GetShowTemplates()
{
    return do_show_templates;
}

float GLWidget::GetGenerateGridSizeX()
{
    return generate_grid_sizex;
}

float GLWidget::GetGenerateGridSizeY()
{
    return generate_grid_sizey;
}

float GLWidget::GetGenerateGridStapleSize(){
    return generate_grid_staplesize;
}

float GLWidget::GetGenerateLinearSpacing()
{
    return generate_linear_spacing;
}

bool GLWidget::GetGenerateLinearScaleX()
{
    return generate_linear_scalex;
}

bool GLWidget::GetGenerateLinearScaleY()
{
    return generate_linear_scaley;
}

float GLWidget::GetGenerateSlicesSpacing()
{
    return generate_slices_spacing;
}

bool GLWidget::GetGenerateSlicesX()
{
    return generate_slices_x;
}

bool GLWidget::GetGenerateSlicesY()
{
    return generate_slices_y;
}

bool GLWidget::GetGenerateSlicesZ()
{
    return generate_slices_z;
}

int GLWidget::GetGenerateBlendSections()
{
    return generate_blend_sections;
}

int GLWidget::GetGenerateRevolveSections()
{
    return generate_revolve_sections;
}

float GLWidget::GetGenerateBranchingScaleChild()
{
    return generate_branching_scalechild;
}

int GLWidget::GetGenerateRadialSectors()
{
    return generate_radial_sectors;
}

bool GLWidget::GetGenerateSurfaceFacetsTeeth()
{
    return generate_surfacefacets_teeth;
}

double GLWidget::GetPhysicsNewWeightMass()
{
    return physics_new_weight_mass;
}

double GLWidget::GetPhysicsMaterialDensity()
{
    return physics_material_density;
}

double GLWidget::GetPhysicsMaximumStress()
{
    return physics_max_stress;
}

bool GLWidget::IsSectionSelected() const
{
    return (selected >= 0 && selected < sections.size());
}

QString GLWidget::GetOpenFilename() const
{
    return open_filename;
}

void GLWidget::SetDoMagneticCuts(const bool b)
{
    do_magnetic_cuts = b;

    UpdateTemplateCut();
}

void GLWidget::SetPenModeOn(const bool b)
{
    pen_mode = b;
    if (state == STATE_PEN_POINT) {
        CancelPenCurve();
    }
}

void GLWidget::SetDoLocalSymmetry(const bool b)
{
    do_local_symmetry = b;
}

void GLWidget::SetDoCyclesTest(const bool b)
{
    do_cycles_test = b;

    DoCyclesConnectedTest();
    update_sections_disp_list = true;
}

void GLWidget::SetDoStabilityTest(const bool b)
{
    do_stability_test = b;
    DoStabilityTest();
}

void GLWidget::SetDoPhysicsTest(const bool b)
{
    do_physics_test = b;

    DoPhysicsTest();

    update_sections_disp_list = true;

    testPhysicsButton->setChecked(do_physics_test);
}

void GLWidget::ToggleDoPhysicsTest()
{
    do_physics_test = !do_physics_test;

    if(do_physics_test)
        DoPhysicsTest();

    update_sections_disp_list = true;

    testPhysicsButton->setChecked(do_physics_test);

}

 void GLWidget::SetDoConnectedTest(const bool b)
 {
     do_connected_test = b;

     DoCyclesConnectedTest();

     update_sections_disp_list = true;
 }

 void GLWidget::SetShowTNBFrames(const bool b)
 {
     do_show_tnb_frames = b;
 }

 void GLWidget::SetShowShadow(const bool b)
 {
     do_show_shadow = b;
     update_sections_disp_list = true;
 }

 void GLWidget::SetShowTemplates(const bool b)
 {
     do_show_templates = b;
 }

void GLWidget::SetSlabThickness(const double d)
{
    slab_thickness = d;

    if (last_op == OP_GENERATE_SURFACE_FACETS) {
        Undo(OP_GENERATE_SURFACE_FACETS);
        DoGenerateSurfaceFacets();
    }

    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness);
        sections[i].UpdateCurveTrisSlab();
    }

    thick_label->setText(QString::number(slab_thickness) + QString(" units"));

    UpdateAllTests();
    UpdateDraw();
}

void GLWidget::SetCalibrationFactor(const int i)
{
    calibration_factor = float(i) / 100.0f;
    calibration_label->setText(QString::number(calibration_factor * 100.0f) + QString("%"));
}

void GLWidget::SetQualitySamples(const int i)
{
    quality_samples = i;

    active_section.SetQuality(quality_samples);
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetQuality(quality_samples);
        sections[i].UpdateCurveTrisSlab();
    }

    quality_label->setText(QString::number(quality_samples));

    UpdateDraw();
}

void GLWidget::SetRotationDuration(const int i)
{
    cam_animate_duration = float(i * 100);
    rotate_label->setText(QString::number(i * 100) + QString(" ms"));
}

void GLWidget::SetRotationAngle(const int i)
{
    cam_animate_angle = float(i);
    rotate_angle_label->setText(QString::number(i) + QString(" degrees"));
}

void GLWidget::SetGenerateGridSizeX(const int i)
{
    generate_grid_sizex = float(i) / 10.0f;
    grid_sizex_label->setText(QString::number(generate_grid_sizex));

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_GRID) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateGridSizeY(const int i)
{
    generate_grid_sizey = float(i) / 10.0f;
    grid_sizey_label->setText(QString::number(generate_grid_sizey));

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_GRID) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateGridStapleSize(const int i)
{
    generate_grid_staplesize = float(i) / 20.0f;
    grid_staple_label->setText(QString::number(generate_grid_staplesize));

    if(current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_GRID) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateLinearSpacing(const int i)
{
    generate_linear_spacing = float(i) / 10.0f;
    linear_spacing_label->setText(QString::number(generate_linear_spacing));

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_LINEAR) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateLinearScaleX(const bool b)
{
    generate_linear_scalex = b;

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_LINEAR) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateLinearScaleY(const bool b)
{
    generate_linear_scaley = b;

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_LINEAR) {
        ShowGenerate();
    }
}

void GLWidget::UpdateGenerateSlicesSpacing()
{
    if(current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_SLICES)
        ShowGenerate();
}

void GLWidget::SetGenerateSlicesSpacing(const int i)
{
    generate_slices_spacing = float(i) / 10.0f;
    slices_spacing_label->setText(QString::number(generate_slices_spacing));

}

void GLWidget::SetGenerateSlicesX(const bool b)
{
    generate_slices_x = b;

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_SLICES) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateSlicesY(const bool b)
{
    generate_slices_y = b;

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_SLICES) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateSlicesZ(const bool b)
{
    generate_slices_z = b;

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_SLICES) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateBranchingScaleChild(const int i)
{
    generate_branching_scalechild = float(i) / 100.0f;
    branching_scalechild_label->setText(QString::number(generate_branching_scalechild));

    if (last_op == OP_GENERATE_BRANCHING) {
        Undo(OP_GENERATE_BRANCHING);
        DoGenerateBranching();
    }
}

void GLWidget::SetGenerateRadialSectors(const int i)
{
    generate_radial_sectors = i;
    radial_sectors_label->setText(QString::number(generate_radial_sectors));

    if(IsSectionSelected() && sections[selected].IsRadial())
    {
        sections[selected].SetNumRadialSectors(generate_radial_sectors);
        sections[selected].UpdateRadial();
    }

    if (last_op == OP_GENERATE_MAKE_RADIAL_HOLE) {
        Undo(OP_GENERATE_MAKE_RADIAL_HOLE);
        DoGenerateMakeRadialHole();
    }
}

void GLWidget::SetGenerateSurfaceFacetsTeeth(const bool b)
{
    generate_surfacefacets_teeth = b;
}

void GLWidget::SetGenerateBlendSections(const int i)
{

    generate_blend_sections = i;
    num_blend_sections_label->setText(QString::number(generate_blend_sections));

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_BLEND) {
        ShowGenerate();
    }
}

void GLWidget::SetGenerateRevolveSections(const int i)
{
    generate_revolve_sections = i; 
    num_revolve_sections_label->setText(QString::number(generate_revolve_sections));

    if (current_tool_state == TOOLSTATE_GENERATE && gen_state == GENSTATE_REVOLVE) {
        ShowGenerate();
    }
}

void GLWidget::AddToUndoList(const LastOperation op)
{

    last_op = op;

    //if we're back in the undo list, remove the undo steps ahead of us
    while (undo_sections.size()-1 > undo_index) {
        undo_sections.pop_back();
    }

    //now add this config to the end of the list
    undo_sections.push_back(sections);
    ++undo_index;

    //and make sure that the undo list doesn't grow too large
    if (undo_sections.size() > max_undo_sections) {
        undo_sections.pop_front();
        --undo_index;
    }

    update_sections_disp_list = true;

    if (QString::compare(open_filename.right(1), "*") != 0) {
        open_filename += "*";
    }

}

void GLWidget::Undo(const LastOperation op)
{

    if (current_tool_state == TOOLSTATE_GENERATE) {
        return;
    }

    //add the existing sections to the undo list, but ONLY if we are at the end
    if (undo_index == undo_sections.size()-1) {
        AddToUndoList(op);
        --undo_index;
    }
    if (undo_index >= 0) {
        sections = undo_sections[undo_index];
        --undo_index;
    }

    if (!IsSectionSelected()) {
        selected = -1;
    }   

}

void GLWidget::Redo()
{
    if (current_tool_state == TOOLSTATE_GENERATE) {
        return;
    }

    if (undo_index < undo_sections.size()-2) {
        sections = undo_sections[undo_index+2];
        ++undo_index;

    }

    if (!IsSectionSelected()) {
        selected = -1;
    }

}

void GLWidget::DeleteSelected()
{

    if (state == STATE_PEN_POINT || !IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_DELETE_PLANE);

    sections.removeAt(selected);
    last_selected.clear();
    SetSelected(-1);

    UpdateAllTests();
    //updateGL();

}

void GLWidget::Transform()
{
    if(state == STATE_TRANSFORM_WIDGET)
    {
        state = STATE_NONE;
        current_tool_state = TOOLSTATE_DEFAULT;
        return;
    }

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_MANIP_TRANSFORM);

    state = STATE_TRANSFORM_WIDGET;
    current_tool_state = TOOLSTATE_TRANSFORMING;
}

void GLWidget::CopyMirrorX()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_COPY_MIRRORX);

    PlanarSection new_section;
    sections[selected].MirrorPlanarSectionX(new_section);
    new_section.UpdateCurveTrisSlab();
    sections.push_back(new_section);

    UpdateAllTests();
    //updateGL();

}

void GLWidget::CopyRotateY()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_COPY_ROTATEY);

    PlanarSection new_section;
    sections[selected].RotatePlanarSectionY(new_section);
    new_section.UpdateCurveTrisSlab();
    sections.push_back(new_section);

    UpdateAllTests();
    //updateGL();

}

void GLWidget::CopyMirrorZ()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_COPY_MIRRORZ);

    PlanarSection new_section;
    sections[selected].MirrorPlanarSectionZ(new_section);
    new_section.UpdateCurveTrisSlab();
    sections.push_back(new_section);

    UpdateAllTests();
    //updateGL();

}

void GLWidget::CopyDuplicate()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_COPY_DUPLICATE);

    PlanarSection new_section = sections[selected];
    sections.push_back(new_section);

    UpdateAllTests();
    //updateGL();

}

void GLWidget::SnapToMajorAxis()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_MANIP_SNAP_MAJOR_AXIS);

    //determine closest major axis
    QList <QVector3D> dirs;

    dirs.push_back(QVector3D(1, 0, 0));
    dirs.push_back(QVector3D(-1, 0, 0));
    dirs.push_back(QVector3D(0, 1, 0));
    dirs.push_back(QVector3D(0, -1, 0));
    dirs.push_back(QVector3D(0, 0, 1));
    dirs.push_back(QVector3D(0, 0, -1));

    int closest_index = -1;
    float closest_angle = FLT_MAX;

    for (int i=0; i<dirs.size(); ++i) {

        const float each_angle = GLutils::AngleBetweenRad(dirs[i], sections[selected].N());

        if (each_angle < closest_angle) {
            closest_angle = each_angle;
            closest_index = i;
        }

    }

    QVector3D rot_axis = QVector3D::crossProduct(dirs[closest_index], sections[selected].N()).normalized();

    QVector3D new_t = GLutils::RotateVector(sections[selected].T(), rot_axis, -closest_angle);
    QVector3D new_n = GLutils::RotateVector(sections[selected].N(), rot_axis, -closest_angle);
    QVector3D new_b = GLutils::RotateVector(sections[selected].B(), rot_axis, -closest_angle);

    sections[selected].SetT(new_t);
    sections[selected].SetN(new_n);
    sections[selected].SetB(new_b);

    sections[selected].UpdateCurveTrisSlab();

}

void GLWidget::AddHoleBoundary()
{
    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_ADD_HOLE);
    state = STATE_ADD_HOLE;
}

void GLWidget::RemoveHolesBoundary()
{
    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_DELETE_HOLES);

    sections[selected].RemoveHoles();
    sections[selected].UpdateCurveTrisSlab();
}

void GLWidget::ResketchCurve()
{
    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_RESKETCH_CURVE);

    sections[selected].SketchClear(0);
    sections[selected].Update(0, section_error_tolerance);

    state = STATE_RESKETCH_CURVE;
}

void GLWidget::NewPlaneSketch()
{
    ClearAll();
    ResetCamera();
    UpdateDraw();

    open_filename = "untitled";
}

void GLWidget::LoadTemplateCurve()
{

    QString filename = QFileDialog::getOpenFileName(this, tr("Open Curve Template"), QString(""), tr("SVG (*.svg)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return;
    }

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::LoadTemplateOBJ(): File " << filename << " can't be loaded";
        return;
    }

    QByteArray ba = file.readAll();
    file.close();

    QStringList str_list = QString(ba).split("\"");

    for (int i=0; i<str_list.size()-1; ++i) {

        if (str_list[i].contains("path") && str_list[i].contains("d=")) {

            BezierCurve curve;
            curve.LoadFromSVGData(str_list[i+1]);

            active_section = PlanarSection();
            SetupPlanarSection(active_section);

            active_section.SetP(cam.LookAt());
            active_section.SetN(-cam.ViewDir());
            active_section.SetT(cam.GetRightVector());
            active_section.SetB(cam.Up());

            active_section.SetCurve(0, curve);
            active_section.UpdateCurveTrisSlab();

            sections.push_back(active_section);

        }

    }

    UpdateDraw();

}

void GLWidget::LoadTemplateImage()
{

    DeleteTemplateImage();

    QString filename = QFileDialog::getOpenFileName(this, tr("Open Image Template"), QString(""), tr("All Files (*)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return;
    }

    template_image = QImage(filename);
    if (template_image.isNull()) {
        QMessageBox::information(this, tr("Error"), tr("Cannot load %1.").arg(filename));
    }

    makeCurrent();

    //get the OpenGL-friendly image
    QImage GL_formatted_image = QGLWidget::convertToGLFormat(template_image);

    //make sure its not null
    if(GL_formatted_image.isNull()) {
        qDebug() << "GLWidget::LoadTemplateImage() - GL formatted image is NULL";
        return;
    }

    //generate the texture name
    glGenTextures(1, &template_image_tex);

    //bind the texture ID
    glBindTexture(GL_TEXTURE_2D, template_image_tex);

    //texture parameters
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    //generate the texture
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, GL_formatted_image.width(),
                  GL_formatted_image.height(), 0, GL_RGBA,
                  GL_UNSIGNED_BYTE, GL_formatted_image.bits() );

    //bind the texture ID
    glBindTexture(GL_TEXTURE_2D, 0);

    //qDebug() << "texture" << template_image_tex;

}

void GLWidget::LoadTemplateOBJ()
{

    QString filename = QFileDialog::getOpenFileName(this, tr("Open Surface Template"), QString(""), tr("OBJ (*.obj);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".obj", Qt::CaseInsensitive) != 0) {
        filename += QString(".obj");
    }

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::LoadTemplateOBJ(): File " << filename << " can't be loaded";
        return;
    }

    template_verts.clear();
    template_magnetic_verts.clear();
    template_norms.clear();
    template_faces.clear();
    template_facenorms.clear();
    template_poly_faces.clear();
    template_cut.clear();

    QTextStream ifs(&file);

    while (!ifs.atEnd()) {

        const QStringList line = ifs.readLine().split(" ", QString::SkipEmptyParts);

        //skip any blank lines or lines without at least 4 entries (a command + 3 parameters)
        if (line.size() < 4) {
            continue;
        }

        if (QString::compare(line.first(), "v") == 0) {
            template_verts.push_back(QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat()));
        }
        else if (QString::compare(line.first(), "vn") == 0) {
            template_norms.push_back(QVector3D(line[1].toFloat(), line[2].toFloat(), line[3].toFloat()));
        }
        else if (QString::compare(line.first(), "f") == 0) {

            //populate poly faces
            QList <int> face;
            for (int j=0; j<line.size(); ++j) {
                const QStringList each_v = line[j].split("/");
                const int v = each_v.first().toInt();
                if (v > 0) { //indexes in OBJ should always be greater than 0, it's 0 if int conversion failed
                    face.push_back(v-1);
                }
            }
            if (face.size() >= 3) {
                template_poly_faces.push_back(face);
            }

            //for a face line, triangle face entries come from 1-2-3, 1-3-4, 1-4-5, ...
            const QStringList each_v1 = line[1].split("/");          

            //populate "faces" - which is just triangles
            for (int j=3; j<line.size(); ++j) {

                const QStringList each_v2 = line[j-1].split("/");
                const QStringList each_v3 = line[j].split("/");

                const int ind_v1 = each_v1.first().toInt()-1;
                const int ind_v2 = each_v2.first().toInt()-1;
                const int ind_v3 = each_v3.first().toInt()-1;

                template_faces.push_back(ind_v1); //note: 1-indexed in file, but 0-indexed in program
                template_faces.push_back(ind_v2);
                template_faces.push_back(ind_v3);

                if (each_v1.size() == 3) { //there is included normal information of format v1/vt1/vn1 or v1//vn1

                    const int ind_vn1 = each_v1.last().toInt()-1;
                    const int ind_vn2 = each_v1.last().toInt()-1;
                    const int ind_vn3 = each_v1.last().toInt()-1;

                    template_facenorms.push_back(ind_vn1);
                    template_facenorms.push_back(ind_vn2);
                    template_facenorms.push_back(ind_vn3);

                }
                else { //we compute a per-face normal, but shading won't be smooth across faces, TODO: normal smoothing

                    const QVector3D dir1 = (template_verts[ind_v2]-template_verts[ind_v1]).normalized();
                    const QVector3D dir2 = (template_verts[ind_v3]-template_verts[ind_v1]).normalized();
                    const QVector3D norm = QVector3D::crossProduct(dir1, dir2);

                    template_facenorms.push_back(template_norms.size());
                    template_facenorms.push_back(template_norms.size());
                    template_facenorms.push_back(template_norms.size());

                    //we make up a norm for this face
                    template_norms.push_back(norm);

                }

            }
        }

    }

    file.close();

    //compute smooth per-face normals
    //qDebug() << "Computing smooth face normals...";
    QVector <QVector3D> new_norms = QVector <QVector3D> (template_verts.size(), QVector3D(0, 0, 0));
    QVector <int> vert_degree = QVector <int> (template_verts.size(), 0);
    QVector <int> new_facenorms;

    for (int i=0; i<template_faces.size(); i+=3) {
        for (int j=0; j<3; ++j) {
            const int ind_v = template_faces[i+j];
            const int ind_vn = template_facenorms[i+j];
            new_norms[ind_v] += template_norms[ind_vn];
            ++vert_degree[ind_v];
            new_facenorms.push_back(ind_v);
        }
    }

    for (int i=0; i<new_norms.size(); ++i) {
        new_norms[i] *= 1.0f / float(vert_degree[i]);
    }

    template_norms = new_norms;
    template_facenorms = new_facenorms;

    //translate and scale the object
    QVector3D bbox_min, bbox_max;
    GLutils::GetBoundingBox(template_verts, bbox_min, bbox_max);
    const float bbox_diam = (bbox_max-bbox_min).length();
    //qDebug() << bbox_min << bbox_max << bbox_diam;
    for (int i=0; i<template_verts.size(); ++i) {
        //translate
        template_verts[i].setX(template_verts[i].x() - (bbox_max.x() + bbox_min.x()) * 0.5f); //move to X centroid
        template_verts[i].setY(template_verts[i].y() - bbox_min.y()); //put above ground plane on Y
        template_verts[i].setZ(template_verts[i].z() - (bbox_max.z() + bbox_min.z()) * 0.5f); //move to Z centroid

        //scale based on bbox diameter
        template_verts[i] *= 12.0f / bbox_diam;
    }

    template_magnetic_verts = template_verts;

    qDebug() << "GLWidget::LoadTemplateOBJ(): File " << filename << "loaded.  " << template_verts.size() << "vertices," << template_faces.size()/3 << "faces.";

    return;

}

bool GLWidget::LoadPlaneSketch()
{

    QString filename = QFileDialog::getOpenFileName(this, tr("Open FlatFab"), QString(), tr("FlatFab (*.txt);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    const QString extension(".txt");
    if (QString::compare(filename.right(extension.length()), extension) != 0) {
        filename += QString(".txt");
    }

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::LoadPlaneSketch(): File " << filename << " can't be loaded";
        return false;
    }

    open_filename = QFileInfo(filename).fileName();

    QTextStream ifs(&file);
    int nPlanes = PlanarSection::LoadHeaderFromFile(ifs);

    ClearAll();
    for (int i=0; i<nPlanes; ++i) {

        PlanarSection new_section;
        new_section.LoadFromFile(ifs);
        new_section.SetSlabThickness(slab_thickness);
        new_section.UpdateCurveTrisSlab();
        sections.push_back(new_section);

    }

    file.close();
    qDebug() << "GLWidget::LoadPlaneSketch(): File " << filename << "loaded.";   

    UpdateAllTests();
    UpdateDraw();

    SetViewIso1();

    return true;

}

bool GLWidget::SavePlaneSketch()
{

    //get filename
    QString filename = QFileDialog::getSaveFileName(this, tr("Save FlatFab"), QString(""), tr("FlatFab (*.txt);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return false;
    }

    //add extension if missing
    const QString extension(".txt");
    if (QString::compare(filename.right(extension.length()), extension) != 0) {
        filename += extension;
    }

    //open file for writing
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SavePlaneSketch(): File " << filename << " can't be saved";
        return false;
    }

    open_filename = QFileInfo(filename).fileName();

    //save out the data
    QTextStream ofs(&file);

    PlanarSection::SaveHeaderToFile(sections, ofs); //header
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SaveToFile(ofs); //each plane's data
    }

    //close file, report saving ok
    file.close();
    qDebug() << "GLWidget::SavePlaneSketch(): File " << filename << "saved.";

    return true;

}

void GLWidget::SaveSliceOBJ()
{

    QString filename = QFileDialog::getSaveFileName(this, tr("Save Slice OBJ"), QString(""), tr("OBJ (*.obj);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".obj") != 0) {
        filename += QString(".obj");
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveSliceOBJ(): File " << filename << " can't be exported";
        return;
    }

    QTextStream ofs(&file);

    ofs << "# Created with PlanarDisc\n";    
    PlanarSection::SaveConnectionsToOBJ(sections, ofs);

    int face_count = 0;
    for (int i=0; i<sections.size(); ++i) {

        const QList <QVector2D> & tris = sections[i].SliceTriangles();

        //group and shared normal
        ofs << "g plane" << i << "\n";
        ofs << "vn " << sections[i].N().x() << " " << sections[i].N().y() << " " << sections[i].N().z() << "\n";

        //triangle vertices
        for (int j=0; j<tris.size(); ++j) {
            const QVector3D v = sections[i].GetPoint3D(tris[j]);
            ofs << "v " << v.x() << " " << v.y() << " " << v.z() << "\n";
        }
        //triangle faces
        for (int j=0; j<tris.size(); j+=3) {
            ofs << "f " << face_count+1 << "//" << i+1 << " " << face_count+2 << "//" << i+1 << " " << face_count+3 << "//" << i+1 << "\n";
            face_count+=3;
        }

    }

    file.close();

    qDebug() << "GLWidget::SaveSliceOBJ(): File " << filename << "saved.";

    return;

}

void GLWidget::SaveSlabOBJ()
{

    QString filename = QFileDialog::getSaveFileName(this, tr("Save Slab OBJ"), QString(""), tr("OBJ (*.obj);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".obj") != 0) {
        filename += QString(".obj");
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveSlabOBJ(): File " << filename << " can't be exported";
        return;
    }

    QTextStream ofs(&file);

    ofs << "# Created with FlatFab\n";
    PlanarSection::SaveConnectionsToOBJ(sections, ofs);

    int vert_count = 0;
    int norm_count = 0;

    for (int i=0; i<sections.size(); ++i) {

        const QList <QVector3D> & vs = sections[i].SlabVertices();
        const QList <QVector3D> & ns = sections[i].SlabNormals();

        //group and shared normal
        ofs << "g plane" << i << "\n";

        //triangle vertices
        for (int j=0; j<vs.size(); ++j) {
            ofs << "v " << vs[j].x() << " " << vs[j].y() << " " << vs[j].z() << "\n";
        }
        //triangle normals
        for (int j=0; j<ns.size(); ++j) {
            ofs << "vn " << ns[j].x() << " " << ns[j].y() << " " << ns[j].z() << "\n";
        }
        //triangle faces
        for (int j=0; j<vs.size(); j+=3) {
            ofs << "f " << vert_count+1 << "//" << norm_count+1 << " " << vert_count+2 << "//" << norm_count+1 << " " << vert_count+3 << "//" << norm_count+1 << "\n";
            vert_count += 3;
            ++norm_count;
        }

    }

    file.close();

    qDebug() << "GLWidget::SaveSlabOBJ(): File " << filename << "saved.";

    return;

}

void GLWidget::SaveFlattenedSlabOBJ()
{

    QString filename = QFileDialog::getSaveFileName(this, tr("Save Flattened Slab OBJ"), QString(""), tr("OBJ (*.obj);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".obj") != 0) {
        filename += QString(".obj");
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveFlattenedSlabOBJ(): File " << filename << " can't be exported";
        return;
    }

    //1. compute the 2D packing so each piece has a location there
    QVector <QVector2D> pack_pts;
    QVector2D bbox;
    PlanarSection::ComputePacking(sections, 2.0f, pack_pts, bbox);

    //2. write out each slab, however map each slab to a coordinate frame that rests on the 2D ground plane

    QTextStream ofs(&file);

    ofs << "# Created with FlatFab\n";
    PlanarSection::SaveConnectionsToOBJ(sections, ofs);

    int vert_count = 0;
    int norm_count = 0;

    for (int i=0; i<sections.size(); ++i) {

        const QVector3D t = sections[i].T();
        const QVector3D n = sections[i].N();
        const QVector3D b = sections[i].B();
        const QVector3D p = sections[i].P();

        QVector2D min_v, max_v;
        sections[i].GetBoundingBox2D(min_v, max_v);

        sections[i].SetT(QVector3D(1,0,0));
        sections[i].SetN(QVector3D(0,0,1));
        sections[i].SetB(QVector3D(0,1,0));
        sections[i].SetP(pack_pts[i] - min_v);

        sections[i].UpdateCurveTrisSlab();

        const QList <QVector3D> & vs = sections[i].SlabVertices();
        const QList <QVector3D> & ns = sections[i].SlabNormals();

        //group and shared normal
        ofs << "g plane" << i << "\n";

        //triangle vertices
        for (int j=0; j<vs.size(); ++j) {
            ofs << "v " << vs[j].x() << " " << vs[j].y() << " " << vs[j].z() << "\n";
        }
        //triangle normals
        for (int j=0; j<ns.size(); ++j) {
            ofs << "vn " << ns[j].x() << " " << ns[j].y() << " " << ns[j].z() << "\n";
        }
        //triangle faces
        for (int j=0; j<vs.size(); j+=3) {
            ofs << "f " << vert_count+1 << "//" << norm_count+1 << " " << vert_count+2 << "//" << norm_count+1 << " " << vert_count+3 << "//" << norm_count+1 << "\n";
            vert_count += 3;
            ++norm_count;
        }

        sections[i].SetT(t);
        sections[i].SetN(n);
        sections[i].SetB(b);
        sections[i].SetP(p);

        sections[i].UpdateCurveTrisSlab();

    }

    file.close();

    qDebug() << "GLWidget::SaveFlattenedSlabOBJ(): File " << filename << "saved.";

    return;

}

void GLWidget::SaveSurfaceOBJ()
{

    QString filename = QFileDialog::getSaveFileName(this, tr("Save Surface OBJ"), QString(""), tr("OBJ (*.obj);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".obj") != 0) {
        filename += QString(".obj");
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveSurfaceOBJ(): File " << filename << " can't be exported";
        return;
    }

    if (contour_graph.GetNumPatches() == 0) {
        contour_graph.Create(sections);
    }

    QTextStream ofs(&file);

    ofs << "# Created with FlatFab\n";
    ofs << "# Mesh based on " << sections.size() << " planes and " << contour_graph.GetNumPatches() << "patches\n";

    contour_graph.SaveToOBJFile(ofs);

    file.close();

    qDebug() << "GLWidget::SaveSurfaceOBJ(): File " << filename << "saved.";

    return;
}

void GLWidget::SavePhysicsOutput()
{

    QString filename = QFileDialog::getSaveFileName(this, tr("Save Physics Output"), QString(""), tr("TXT (*.txt);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);

    if (filename.isNull()) {
        return;
    }

    if (QString::compare(filename.right(4), ".txt") != 0) {
        filename += QString(".txt");
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SavePhysicsOutput(): File " << filename << " can't be exported";
        return;
    }

    QTextStream ofs(&file);

    ofs << "PlaneSections " << sections.size() << "\n";

    // first add each rigid body
    for (int i=0; i<sections.size(); ++i) {

        ofs << "PlaneIndex " << i << "\n";

        //P and TNB frame
        const QVector3D p = sections[i].P();
        ofs << "PlaneP " << p.x() << " " << p.y() << " " << p.z() << "\n";

        const QVector3D t = sections[i].T();
        ofs << "PlaneT " << t.x() << " " << t.y() << " " << t.z() << "\n";

        const QVector3D n = sections[i].N();
        ofs << "PlaneN " << n.x() << " " << n.y() << " " << n.z() << "\n";

        const QVector3D b = sections[i].B();
        ofs << "PlaneB " << b.x() << " " << b.y() << " " << b.z() << "\n";

        const QVector3D mass_cent = sections[i].GetCentroid3D(); //centre of mass
        ofs << "PlaneCentreOfMass3D " << mass_cent.x() << " " << mass_cent.y() << " " << mass_cent.z() << "\n";

        const float area = fabsf(sections[i].SignedArea());
        ofs << "PlaneArea " << area << "\n";

        const float thick = sections[i].SlabThickness();
        ofs << "PlaneThickness " << thick << "\n";

        QList <QVector3D> contact_pts; //contact points
        sections[i].GetXZPlanePoints(contact_pts);
        ofs << "PlaneContactPoints3D " << contact_pts.size() << "\n";
        for (int j=0; j<contact_pts.size(); ++j) {
            ofs << contact_pts[j].x() << " " << contact_pts[j].y() << " " << contact_pts[j].z() << "\n";
        }

        ofs << "PlaneNumBoundaries " << sections[i].GetNumCurves() << "\n";
        for (int c=0; c<sections[i].GetNumCurves(); ++c) {
            const QList <QVector2D> & samples = sections[i].GetCurve(c).Samples();
            ofs << "PlaneBoundary2D " << samples.size() << "\n";
            for (int j=0; j<samples.size(); ++j) {
                ofs << samples[j].x() << " " << samples[j].y() << "\n";
            }
        }
    }

    //joints
    QList <int> joint_p1i;
    QList <int> joint_p2i;
    QList <QVector3D> joint_centre;
    QList <QVector3D> joint_axis;
    QList <QVector2D> joint_p1_pt1;
    QList <QVector2D> joint_p1_pt2;
    QList <QVector2D> joint_p2_pt1;
    QList <QVector2D> joint_p2_pt2;
    for (int i=0; i<sections.size(); ++i) {

        // then add the joint info
        for (int j=i+1; j<sections.size(); ++j) {

            if (i == j) {
                continue;
            }

            //int index, other_index;
            QList <QVector2D> my_slots;
            QList <QList <QVector2D> > my_slots_rect;

            //sections[i].GetSlots(sections[j], index, other_index, my_slots, my_slots_rect);
            sections[i].GetSlots(sections[j], my_slots, my_slots_rect);

            if (my_slots.size() >= 2) {

                joint_p1i.push_back(i);
                joint_p2i.push_back(j);
                joint_centre.push_back(sections[i].GetPoint3D(my_slots[1]));
                joint_axis.push_back((sections[i].GetPoint3D(my_slots[1]) - sections[i].GetPoint3D(my_slots[0])).normalized());
                joint_p1_pt1.push_back(my_slots[1]);
                joint_p1_pt2.push_back(my_slots[0]);

                joint_p2_pt1.push_back(sections[j].GetPoint2D(sections[i].GetPoint3D(my_slots[1])));
                joint_p2_pt2.push_back(sections[j].GetPoint2D(sections[i].GetPoint3D(my_slots[1]-my_slots[0]+my_slots[1])));

                //qDebug() << "adding joint at" << position << "between" << i << j;

            }
        }
    }
    ofs << "Joints " << joint_p1i.size() << "\n";
    for (int i=0; i<joint_p1i.size(); ++i) {
        ofs << "JointPlaneIndexes " << joint_p1i[i] << " " << joint_p2i[i] << "\n";
        ofs << "JointCentre3D " << joint_centre[i].x() << " " << joint_centre[i].y() << " " << joint_centre[i].z() << "\n";
        ofs << "JointAxis3D " << joint_axis[i].x() << " " << joint_axis[i].y() << " " << joint_axis[i].z() << "\n";
        ofs << "JointPlane1SlitEndpoints2D " << joint_p1_pt1[i].x() << " " << joint_p1_pt1[i].y() << " " << joint_p1_pt2[i].x() << " " << joint_p1_pt2[i].y() << "\n";;
        ofs << "JointPlane2SlitEndpoints2D " << joint_p2_pt1[i].x() << " " << joint_p2_pt1[i].y() << " " << joint_p2_pt2[i].x() << " " << joint_p2_pt2[i].y() << "\n";
    }

    //save out any external forces as well
    int total_forces = 0;
    for (int i=0; i<sections.size(); ++i) {
        total_forces += sections[i].GetNumWeights();
    }
    ofs << "ExternalForces " << total_forces << "\n";
    for (int i=0; i<sections.size(); ++i) {
        for (int j=0; j<sections[i].GetNumWeights(); ++j) {
            QVector3D pos = sections[i].GetWeightPosition(j);
            QVector3D force = sections[i].GetWeightForce(j);
            ofs << i << " " << pos.x() << " " << pos.y() << " " << pos.z() << " " << force.x() << " " << force.y() << " " << force.z() << "\n";
        }
    }

    file.close();

    qDebug() << "GLWidget::SavePhysicsOutput(): File " << filename << "saved.";

    return;

}

void GLWidget::SaveSVG()
{

    bool ok;

    int svg_dpi = QInputDialog::getInt(this, QString("Save SVG File - DPI"), QString("Enter DPI (dots per inch) for SVG.\n(Illustrator and many other programs use 72dpi (default), Inkscape uses 90dpi.)"), 72, 1, 1000, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    double width_height_ratio = QInputDialog::getDouble(this, QString("Save SVG File - Aspect Ratio"), QString("Enter cutting bed width to height ratio.\n(For example, a 24\"width X 12\"height ratio would be 2.0.)"), 1.5, 0.1, 10.0, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    bool use_numeric_labels;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, QString("Save SVG File - Numeric Labels"), QString("Do you want numeric labels to appear in the file?\n(Numeric labels help simplify the assembly process.)"), QMessageBox::No | QMessageBox::Yes);
    use_numeric_labels = (reply == QMessageBox::Yes);

    qDebug() << svg_dpi << width_height_ratio << use_numeric_labels;

    QString filename = QFileDialog::getSaveFileName(this, tr("Save SVG File"), QString(""), tr("SVG (*.svg);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return;
    }

    //add extension if missing
    const QString extension(".svg");
    if (QString::compare(filename.right(extension.length()), extension) != 0) {
        filename += extension;
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveSVG(): File " << filename << " can't be exported";
        return;
    }

    //quickly scale the thickness by the calibration factor
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness * calibration_factor);
    }

    QList <QString> labels;
    QTextStream ofs(&file);\
    PlanarSection::SaveToSVG(sections, ofs, metres_per_unit, svg_dpi, width_height_ratio, use_numeric_labels, labels);
    file.close();

    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness);
    }

    qDebug() << "GLWidget::SaveSVG(): File " << filename << "saved.";

    return;

}

void GLWidget::SaveCalibration()
{
    bool ok;

    int svg_dpi = QInputDialog::getInt(this, QString("Save Calibration SVG File - DPI"), QString("Enter DPI (dots per inch) for SVG.\n(Illustrator and many other programs use 72dpi (default), Inkscape uses 90dpi.)"), 72, 1, 1000, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    int calib_spacing = QInputDialog::getInt(this, QString("Save Calibration SVG File - Spacing value"), QString("Enter lower bound for calibration value."), 70, 1, 99, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    int calib_increment = QInputDialog::getInt(this, QString("Save Calibration SVG File - Increment value"), QString("Enter increment for calibration value."), 2, 1, 20, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    QString filename = QFileDialog::getSaveFileName(this, tr("Save SVG File"), QString(""), tr("SVG (*.svg);;All Files (*)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return;
    }

    //add extension if missing
    const QString extension(".svg");
    if (QString::compare(filename.right(extension.length()), extension, Qt::CaseInsensitive) != 0) {
        filename += extension;
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveCalibration(): File " << filename << " can't be exported";
        return;
    }

    //Since this is a calibration shape, scale all but the 1st
    /*
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness * calibration_factor);
    }
    */
    QList <PlanarSection> calib_sections;

    int nSecs = (100 - calib_spacing)/calib_increment + 2;

    PlanarSection p;
    p.SetT(QVector3D(1,0,0));
    p.SetB(QVector3D(0,1,0));
    p.SetN(QVector3D(0,0,1));
    p.SetP(QVector3D(0,0,0));

    BezierCurve b;
    b.AddPoint(QVector2D(0, 0)); //p1
    b.AddPoint(QVector2D(0, 1));
    b.AddPoint(QVector2D(0,nSecs-1));
    b.AddPoint(QVector2D(0,nSecs)); //p2
    b.AddPoint(QVector2D(-1,nSecs));
    b.AddPoint(QVector2D(-1,nSecs));
    b.AddPoint(QVector2D(-2,nSecs)); //p3
    b.AddPoint(QVector2D(-2,nSecs-1));
    b.AddPoint(QVector2D(-2,1));
    b.AddPoint(QVector2D(-2,0)); //p4
    b.AddPoint(QVector2D(-1,0));
    b.AddPoint(QVector2D(-1,0));
    b.AddPoint(QVector2D(0,0)); //return to p1

    p.SetCurve(0, b);
    p.SetSlabThickness(slab_thickness);

    calib_sections.push_back(p);
    calib_sections.last().UpdateCurveTrisSlab();

    int cury = 1;

    QList <QString> labels;
    labels.push_back("calibration_shape");

    for (int i=calib_spacing; i<=100; i+=calib_increment) {

        PlanarSection p;
        p.SetT(QVector3D(1,0,0));
        p.SetB(QVector3D(0,0,1));
        p.SetN(QVector3D(0,1,0));
        p.SetP(QVector3D(0,cury,0));

        BezierCurve b;
        b.AddPoint(QVector2D(-1,-1)); //p1
        b.AddPoint(QVector2D(0,-1));
        b.AddPoint(QVector2D(0,-1));
        b.AddPoint(QVector2D(1,-1)); //p2
        b.AddPoint(QVector2D(1,0));
        b.AddPoint(QVector2D(1,0));
        b.AddPoint(QVector2D(1,1)); //p3
        b.AddPoint(QVector2D(0,1));
        b.AddPoint(QVector2D(0,1));
        b.AddPoint(QVector2D(-1,1)); //p4
        b.AddPoint(QVector2D(-1,0));
        b.AddPoint(QVector2D(-1,0));
        b.AddPoint(QVector2D(-1,-1)); //return to p1

        p.SetCurve(0, b);
        p.SetSlabThickness(slab_thickness * float(i)/100.0f);

        calib_sections.push_back(p);
        calib_sections.last().UpdateCurveTrisSlab();

        labels.push_back(QString::number(i) + QString("%"));

        ++cury;
    }

    //labels.clear();

    QTextStream ofs(&file);\
    PlanarSection::SaveToSVG(calib_sections, ofs, metres_per_unit, svg_dpi, 1.0f, true, labels);
    file.close();

    /*
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness);
    }
    */

    qDebug() << "GLWidget::SaveCalibration(): File " << filename << "saved.";

    return;
}

void GLWidget::SaveDXF()
{

    bool ok;
    double width_height_ratio = QInputDialog::getDouble(this, QString("Save DXF File - Aspect Ratio"), QString("Enter cutting bed width to height ratio.\n(For example, a 24\"width X 12\"height ratio would be 2.0.)"), 1.5, 0.1, 10.0, 1, &ok, Qt::Dialog);
    if (!ok) {
        return;
    }

    const int dpi = 90;

    bool use_numeric_labels;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, QString("Save DXF File - Numeric Labels"), QString("Do you want numeric labels to appear in the file?\n(Numeric labels help simplify the assembly process.)"), QMessageBox::No | QMessageBox::Yes);
    use_numeric_labels = (reply == QMessageBox::Yes);

    qDebug() << dpi << width_height_ratio << use_numeric_labels;

    QString filename = QFileDialog::getSaveFileName(this, tr("Save DXF File"), QString(""), tr("DXF (*.dxf);; All Files (*)"), 0, QFileDialog::DontUseNativeDialog);
    if (filename.isNull()) {
        return;
    }

    //add extension if missing
    const QString extension(".dxf");
    if (QString::compare(filename.right(extension.length()), extension) != 0) {
        filename += extension;
    }

    QFile file(filename);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "GLWidget::SaveDXF(): File " << filename << " can't be exported";
        return;
    }   

    //quickly scale the thickness by the calibration factor
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness * calibration_factor);
    }

    QList <QString> labels;

    QTextStream ofs(&file);    
    PlanarSection::SaveToDXF(sections, ofs, metres_per_unit, dpi, width_height_ratio, use_numeric_labels, labels);
    file.close();

    //restore the thickness
    for (int i=0; i<sections.size(); ++i) {
        sections[i].SetSlabThickness(slab_thickness);
    }

    qDebug() << "GLWidget::SaveDXF(): File " << filename << "saved.";

    return;

}

void GLWidget::UpdateDraw()
{
    update_sections_disp_list = true;
    //updateGL();
}

void GLWidget::UpdateAllTests()
{

    DoCyclesConnectedTest();
    DoStabilityTest();
    DoPhysicsTest();

}

void GLWidget::ClearAll()
{

    ClearSelected();

    sections.clear();
    undo_sections.clear();

    contour_graph.Clear();

    template_verts.clear();
    template_magnetic_verts.clear();
    template_norms.clear();
    template_faces.clear();
    template_facenorms.clear();

    recursive_slot_start = QVector3D(0, 0, 0);
    recursive_slot_end = QVector3D(0, 0, 0);

    dimensiontool_start = QVector3D(0, 0, 0);
    dimensiontool_end = QVector3D(0, 0, 0);

    markers.clear();
    markers_col.clear();
    markers2.clear();
    markers_col2.clear();

    UpdateAllTests();

    undo_index = -1;

    state = STATE_NONE;
    current_tool_state = TOOLSTATE_DEFAULT;
    generated_sections.clear();
    generate_selections.clear();

}

void GLWidget::SetSelected(int i)
{

    selected = i;

    if (i != -1) {
        if (!last_selected.empty() && i != last_selected.last()) {
            last_selected.push_back(i);
        }
        else {
            last_selected.push_back(i);
        }
    }

    if (last_selected.size() > 10) {
        last_selected.pop_front();
    }   

}

void GLWidget::ClearSelected()
{
    selected = -1;
    last_selected.clear();
}

void GLWidget::UpdateMarkers()
{

    markers.clear();
    markers_col.clear();

    if (sections.size() >= 1 &&
            (state == STATE_SLOT || state == STATE_CAM_TRANSLATE || state == STATE_DEADZONE ||
             state == STATE_CURVE || state == STATE_PEN_POINT || state == STATE_PEN_DRAG) ) {
        markers.push_back(slot_start);
        markers_col.push_back(QVector3D(1, 0, 0));
        markers.push_back(slot_end);
        markers_col.push_back(QVector3D(1, 1, 0));
    }

    if (!sections.empty() &&
            (state == STATE_RECURSIVE_SETUP_SLOT || state == STATE_RECURSIVE_SLOT || recursive_slot_start != QVector3D(0, 0 ,0))) {
        markers.push_back(recursive_slot_start);
        markers_col.push_back(QVector3D(0, 0 ,1));
        markers.push_back(recursive_slot_end);
        markers_col.push_back(QVector3D(1, 0, 1));
    }   

    //points around convex hull
    if (do_stability_test) {
        DoStabilityTest();
    }

    //compute/display all plane-plane intersections on contours
    /*
    for (int i=0; i<sections.size(); ++i) {

        for (int j=0; j<sections.size(); ++j) {

            if (i == j) {
                continue;
            }

            QList <QVector3D> isecs;
            QList <bool> isecs_which;
            sections[i].GetContourIntersections(sections[j], isecs, isecs_which);

            for(int k=0; k<isecs.size(); ++k) {
                markers.push_back(isecs[k]);
                markers_col.push_back(QVector3D(0.5, 0.5, 1));
            }

        }

    }
    */


}

void GLWidget::DrawMarkers()
{

    for (int i=0; i<markers.size(); ++i) {
        glPushMatrix();
        glTranslatef(markers[i].x(), markers[i].y(), markers[i].z());
        glColor3f(markers_col[i].x(), markers_col[i].y(), markers_col[i].z());
        gluSphere(gluNewQuadric(), 0.1, 20, 20);
        glPopMatrix();
    }

    for (int i=0; i<markers2.size(); ++i) {
        glPushMatrix();
        glTranslatef(markers2[i].x(), markers2[i].y(), markers2[i].z());
        glColor3f(markers_col2[i].x(), markers_col2[i].y(), markers_col2[i].z());
        gluSphere(gluNewQuadric(), 0.15, 10, 10);
        glPopMatrix();
    }

}

void GLWidget::DrawSlot(QVector3D start, QVector3D end)
{

    //draw the slot as a dashed line
    glDisable(GL_DEPTH_TEST);
    glLineWidth(3.0f);
    glColor3f(1, 1, 1);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x00ff);
    glBegin(GL_LINES);

    glVertex3f(start.x(), start.y(), start.z());
    glVertex3f(end.x(), end.y(), end.z());

    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glLineWidth(1.0f);
    glEnable(GL_DEPTH_TEST);

}



void GLWidget::DrawDimensionTool()
{

    //draw the slot as a dashed line
    //glDisable(GL_DEPTH_TEST);
    /*
    glLineWidth(3.0f);
    glColor3f(1, 1, 1);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x00ff);
    glBegin(GL_LINES);

    glVertex3f(dimensiontool_start.x(), dimensiontool_start.y(), dimensiontool_start.z());
    glVertex3f(dimensiontool_end.x(), dimensiontool_end.y(), dimensiontool_end.z());

    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glLineWidth(1.0f);
    */
    //glEnable(GL_DEPTH_TEST);
    glColor3f(1, 1, 1);
    GLutils::DrawSphere(dimensiontool_start, cam.CamWidth() / 100.0f);
    GLutils::DrawSphere(dimensiontool_end, cam.CamWidth() / 100.0f);
    glColor3f(0.6f, 1, 0.6f);
    GLutils::DrawCylinder(dimensiontool_start, dimensiontool_end, cam.CamWidth() / 150.0f);

}

void GLWidget::DrawTemplateImage()
{

    if (template_image.isNull() || template_image_tex <= 0) {
        return;
    }

    const float w = 1.0f;
    const float h = float(template_image.height()) / float(template_image.width());

    QVector3D t, n, b, p;

    if (sections.empty()) {
        t = cam.GetRightVector();
        n = -cam.ViewDir();
        b = cam.Up();
        p = cam.LookAt();
    }
    else {

        if (last_selected.empty()) {
            return;
        }

        const int ls = last_selected.last();

        if (ls < 0 || ls >= sections.size()) {
            return;
        }

        t = sections[ls].T();
        n = sections[ls].N();
        b = sections[ls].B();
        p = sections[ls].P();

    }

    p = p + t * template_pos.x() + b * template_pos.y();

    t = GLutils::RotateVector(t, n, template_rotation * 3.14159f / 180.0f) * template_scale;
    b = GLutils::RotateVector(b, n, template_rotation * 3.14159f / 180.0f) * template_scale;

    const QVector3D p1 = t * (-w ) + b * (-h) + p;
    const QVector3D p2 = t * (w) + b * (-h) + p;
    const QVector3D p3 = t * (w ) + b * (h) + p;
    const QVector3D p4 = t * (-w) + b * (h) + p;

    //enable blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //enable 2d texturing
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, template_image_tex);

    glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    glBegin( GL_QUADS );
    if (!template_flipx) {
        glTexCoord2d(0.0,0.0);
        glVertex3f(p1.x(), p1.y(), p1.z());
        glTexCoord2d(1.0,0.0);
        glVertex3f(p2.x(), p2.y(), p2.z());
        glTexCoord2d(1.0,1.0);
        glVertex3f(p3.x(), p3.y(), p3.z());
        glTexCoord2d(0.0,1.0);
        glVertex3f(p4.x(), p4.y(), p4.z());
    }
    else {
        glTexCoord2d(1.0,0.0);
        glVertex3f(p1.x(), p1.y(), p1.z());
        glTexCoord2d(0.0,0.0);
        glVertex3f(p2.x(), p2.y(), p2.z());
        glTexCoord2d(0.0,1.0);
        glVertex3f(p3.x(), p3.y(), p3.z());
        glTexCoord2d(1.0,1.0);
        glVertex3f(p4.x(), p4.y(), p4.z());
    }
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);

}

void GLWidget::DrawTemplateSurface()
{

    //glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glColor4f(0.6f, 0.9f, 0.6f, 0.25f);
    glBegin(GL_TRIANGLES);
    for (int i=0; i<template_faces.size(); i+=3) {

        for (int j=0; j<3; ++j) {
            const int ind_v = template_faces[i+j];
            const int ind_vn = template_facenorms[i+j];
            glNormal3f(template_norms[ind_vn].x(), template_norms[ind_vn].y(), template_norms[ind_vn].z());
            glVertex3f(template_magnetic_verts[ind_v].x(), template_magnetic_verts[ind_v].y(), template_magnetic_verts[ind_v].z());
        }

    }
    glEnd();

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    //glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

}

void GLWidget::UpdateMagneticCutSurface(const QList <QVector3D> & plane_ns, const QList <QVector3D> & plane_ps)
{

    for (int i=0; i<template_magnetic_verts.size(); ++i) {

        const QVector3D & p = template_verts[i];
        const QVector3D & n = template_norms[i];

        //1.  get the deformed point due to EVERY existing section
        QVector <QVector3D> P_is = QVector <QVector3D> (plane_ps.size(), p);

        bool too_far = true;

        for (int j=0; j<plane_ps.size(); ++j) {

            const QVector3D & pl_n = plane_ns[j];
            const QVector3D & pl_p = plane_ps[j];

            const float dist = QVector3D::dotProduct(pl_n, p - pl_p);
            const float interp_val = (1.0f - (fabsf(dist) / magnet_strength)) * 4.0f; //this makes it smoother

            if (fabsf(dist) < magnet_strength) {

                too_far = false;

                //Rule: within magnet field, vertex should move onto side of plane that
                //the vert's normal points toward

                //Vector3 side = planarCut.params.norm * planarCut.params.norm.Dot(verts[i] - planarCut.params.cent);
                const QVector3D side = pl_n * dist;
                const QVector3D sidenorm = side.normalized();

                //qDebug() << "interpval" << interp_val << "dist" << dist << "magnetval" << magnet_strength;
                const QVector3D new_pt = (p - side) - sidenorm * QVector3D::dotProduct(sidenorm, n) * magnet_strength;

                if (interp_val > 0.0f && interp_val < 1.0f) { //we linterp close to the boundary of the field
                    P_is[j] = new_pt * interp_val + p * (1.0f - interp_val);
                }
                else {
                    P_is[j] = new_pt;
                }

            }

        }


        //2. use Karan's "wires" technique to blend each deformed points together into a single coherent deformation
        /*
        Overall deformation for any point P to be
        P' =  P + sum ( w_i * (P_i - P))
        i is running over all the planes.
        Where w_i = ||P_i-P ||^n  / sum_i(||P_i-P ||^n  )
        Usually n=1 or 1.5 should work well.
        -k
        */

        QVector <float> w_is = QVector <float> (plane_ps.size(), 0.0f);

        // w_i = ||P_i-P ||^n  / sum_i(||P_i-P ||^n  )
        float w_total = 0.0f;
        for (int j=0; j<P_is.size(); ++j) {
            w_is[j] = (P_is[j] - p).length();
            w_total += w_is[j];
        }
        for (int j=0; j<P_is.size(); ++j) {
            w_is[j] /= w_total;
        }

        //P' =  P + sum ( w_i * (P_i - P))
        QVector3D p_prime = p;
        if (w_total > 0.0f) {

            for (int j=0; j<P_is.size(); ++j) {
                p_prime += (P_is[j] - p) * w_is[j];
            }

            if (too_far) {
                qDebug() << p << p_prime << P_is << w_is;
            }

        }

        template_magnetic_verts[i] = p_prime;

    }

}

void GLWidget::ComputeTemplateCut(const QVector3D & n, const QVector3D & p, QList <QVector3D> & cut_segments)
{

    cut_segments.clear();

    //first see if we should do magnetic cuts
    template_magnetic_verts = template_verts;

    if (do_magnetic_cuts) {

        QList <QVector3D> plane_ns, plane_ps;
        for (int i=0; i<sections.size(); ++i) {
            plane_ns.push_back(sections[i].N());
            plane_ps.push_back(sections[i].P());
        }
        plane_ns.push_back(n);
        plane_ps.push_back(p);

        UpdateMagneticCutSurface(plane_ns, plane_ps);

        /*
        for (int i=0; i<template_magnetic_verts.size(); ++i) {

            //1.  compute each of the distances
            QVector <float> distances = QVector <float> (sections.size()+1); //note distances are signed!  (side of plane is stored)
            float total_dists = 0.0f;

            for (int j=0; j<sections.size(); ++j) {
                distances[j] = sections[j].GetPointDistance(template_verts[i]);
                total_dists += fabsf(distances[j]);
            }
            distances[sections.size()] = QVector3D::dotProduct(n, p - template_verts[i]);
            total_dists += fabsf(distances[sections.size()]);

            //2.  compute an offset vector, normalized by distances
            QVector3D offset(0, 0, 0);
            for (int j=0; j<distances.size(); ++j) {
                const float w_j = fabsf(distances[j]) / total_dists;
                if (j < sections.size()) {
                    const QVector3D Pi_P = sections[j].N() * distances[j];
                    offset += Pi_P * w_j;
                }
                else {
                    const QVector3D Pi_P = n * distances[j];
                    offset += Pi_P * w_j;
                }
            }

            template_magnetic_verts[i] += offset;

        }
        */

    }

    for (int i=0; i<template_faces.size(); i+=3) {

        const QVector3D & v1 = template_magnetic_verts[template_faces[i]];
        const QVector3D & v2 = template_magnetic_verts[template_faces[i+1]];
        const QVector3D & v3 = template_magnetic_verts[template_faces[i+2]];

        QVector3D int1, int2, int3;
        bool b1, b2, b3;

        b1 = GLutils::LineSegmentPlaneIntersection(p, n, v1, v2, int1);
        b2 = GLutils::LineSegmentPlaneIntersection(p, n, v2, v3, int2);
        b3 = GLutils::LineSegmentPlaneIntersection(p, n, v3, v1, int3);

        //intersection between plane and triangle, but remove zero-length segments
        if (b1 && b2 && int1 != int2) {
            cut_segments.push_back(int1);
            cut_segments.push_back(int2);
        }
        else if (b2 && b3 && int2 != int3) {
            cut_segments.push_back(int2);
            cut_segments.push_back(int3);
        }
        else if (b3 && b1 && int1 != int3) {
            cut_segments.push_back(int3);
            cut_segments.push_back(int1);
        }

        //JAMES EDIT
        /*
        //if ( GetPlaneTriIntersect(planarCut.params, eachTri, eachP1, eachP2, eachE1, eachE2) && (eachP1-eachP2).Length() > 0.00001f ) {
        if (GetPlaneTriIntersect(planarCut.params, eachTri, eachP1, eachP2, eachE1, eachE2)) {

            //if performing "aesthetic cuts", remove all segments whose
            //absolute value of dot product of surface and plane cut
            //is less than the specified threshold
            if (m_baesthetic) {

                Vector3 & norm = planarCut.params.norm;
                const Vector3 & surfnorm = facenorms[i];

                if (fabs(norm.Dot(surfnorm)) > m_baestheticthresh) {
                    continue;
                }
            }

            MeshSegment newSegment;
            newSegment.pt1 = eachP1;
            newSegment.pt2 = eachP2;
            newSegment.n = facenorms[i];
            newSegment.face = i;

            //save edge (edge defined as pair of vert indexes)
            newSegment.edge1 = MeshEdge(faces[i][eachE1], faces[i][(eachE1+1)%3]);
            newSegment.edge2 = MeshEdge(faces[i][eachE2], faces[i][(eachE2+1)%3]);

            planarCut.segments.push_back(newSegment);

        }
        */

    }

}

void GLWidget::DrawTemplateCut(const QList <QVector3D> & cut_segments)
{

    //draw the slot as a dashed line
    glLineWidth(2.0f);
    glColor3f(0.1, 0.5, 0.1);
    //glEnable(GL_LINE_STIPPLE);
    //glLineStipple(1, 0x00ff);
    glBegin(GL_LINES);

    for (int i=0; i<cut_segments.size(); ++i) {
        glVertex3f(cut_segments[i].x(), cut_segments[i].y(), cut_segments[i].z());
    }

    glEnd();
    //glDisable(GL_LINE_STIPPLE);
    glLineWidth(1.0f);

}

void GLWidget::DeleteTemplateImage()
{

    //deallocate previous texture if one already exists
    template_image = QImage();
    if (template_image_tex > 0) {
        glDeleteTextures(1, &template_image_tex);
        template_image_tex = 0;
    }

}

void GLWidget::QVector3DToArray(const QVector3D & p, double array[3])
{
    array[0] = p.x();
    array[1] = p.y();
    array[2] = p.z();
}

void GLWidget::UpdateTemplateCut()
{

    //do the template cut stuff
    /*
    if (!sections.empty()) {
        ComputeTemplateCut(sections.last().N(), sections.last().P(), template_cut);
    }
    else {
        ComputeTemplateCut(cam.ViewDir(), QVector3D(0, 0, 0), template_cut);
    }
    */
    //NOTE TODO you changed this
    if (sections.empty() && state == STATE_NONE) {
        //active_section = PlanarSection();
        active_section.SetP(cam.LookAt());
        active_section.SetN(-cam.ViewDir());
        active_section.SetT(cam.GetRightVector());
        active_section.SetB(cam.Up());
    }

    ComputeTemplateCut(active_section.N(), active_section.P(), template_cut);

}


void GLWidget::AcceptPenCurve()
{
    //add existing planar section set to undo list
    AddToUndoList(OP_ADD_PLANE);

    active_section.SketchSetEditing(false);

    if (do_local_symmetry) {
        active_section.MirrorControlPoints();
        active_section.UpdateCurveTrisSlab();
    }

    //enforce above-ground control points
    active_section.SetCtrlPointsAboveXZPlane();
    active_section.UpdateCurveTrisSlab();   

    if (!active_section.SliceTriangles().empty()) {
        sections.push_back(active_section);

        //for the added curve, we now update our physical tests (TODO: checkbox/boolean to disable this)
        UpdateAllTests();

    }

    SetSelected(-1);

    slot_end = slot_start;
    state = STATE_NONE;
}


void GLWidget::CancelPenCurve()
{
    active_section.SketchSetEditing(false);
    slot_end = slot_start;
    state = STATE_NONE;
}

void GLWidget::AcceptGenerate()
{

    current_tool_state = TOOLSTATE_DEFAULT;
    if(ShowGenerate())
    {

        switch (gen_state){

        case GENSTATE_LINEAR:
            AddToUndoList(OP_GENERATE_LINEAR);
            break;

        case GENSTATE_BLEND:
            AddToUndoList(OP_GENERATE_BLEND);
            break;

        case GENSTATE_REVOLVE:
            AddToUndoList(OP_GENERATE_REVOLVE);
            break;

        case GENSTATE_SLICES:
            AddToUndoList(OP_GENERATE_SLICES);
            break;

        case GENSTATE_GRID:
            AddToUndoList(OP_GENERATE_GRID);
            sections.removeAt(generate_selections[0]); // Note: this must come after AddToUndoList(OP_GENERATE_GRID); or else the undo will not work
            break;

        default:
            break;

        }
        sections.append(generated_sections);
    }



    generated_sections.clear();
    generate_selections.clear();

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::CancelGenerate()
{
//    ShowGenerate();
    current_tool_state = TOOLSTATE_DEFAULT;
    generated_sections.clear();
    generate_selections.clear();
    UpdateAllTests();
    UpdateDraw();
}


bool GLWidget::ShowGenerate()
{

    generated_sections.clear();

    switch (gen_state) {

    case GENSTATE_SLICES:
    {

        //qDebug() << "GLWidget::ShowGenerate() for GENSTATE_SLICES";

        if (template_verts.empty() || template_faces.empty()) {
            qDebug() << "GLWidget::ShowGenerate() - Template is empty (either no vertices or faces), aborting";

            UpdateAllTests();
            UpdateDraw();
            if (current_tool_state == TOOLSTATE_DEFAULT) { // this is should only be shown during AcceptGenerate()
                QMessageBox::warning ( this, "Generate Slices Error", "Can't Accept Generate Slices - Template is empty (either no vertices or faces).");
            }

            return false;
        }

        //1.  compute bounding box
        QVector3D bbox_min, bbox_max;
        GLutils::GetBoundingBox(template_verts, bbox_min, bbox_max);

        //2.  we sweep through each of X, Y, and Z adding planar sections along spacing defined by the slices_spacing value
        //float generate_slices_spacing;
        //bool generate_slices_x;
        //bool generate_slices_y;
        //bool generate_slices_z;
        for (int i=0; i<3; ++i) {

            //continue if this direction is not enabled
            if ((i == 0 && !generate_slices_x) ||
                    (i == 1 && !generate_slices_y) ||
                    (i == 2 && !generate_slices_z)) {
                continue;
            }

            //2a.  Assign normal sampling direction
            QVector3D t, n, b;

            switch (i) {
            case 0:
                t = QVector3D(0,0,1);
                n = QVector3D(1,0,0);
                b = QVector3D(0,-1,0);
                break;
            case 1:
                t = QVector3D(1,0,0);
                n = QVector3D(0,1,0);
                b = QVector3D(0,0,-1);
                break;
            case 2:
            default:
                t = QVector3D(0,1,0);
                n = QVector3D(0,0,1);
                b = QVector3D(-1,0,0);
                break;
            }

            //2b.  Compute interval range to iterate over (the bounding box dimensions along sampling direction n)
            const float dp_min = QVector3D::dotProduct(n, bbox_min);
            const float dp_max = QVector3D::dotProduct(n, bbox_max);

            const int j_min = floorf(dp_min/generate_slices_spacing);
            const int j_max = ceilf(dp_max/generate_slices_spacing);

            //2c.  iterate by the spacing amount, creating each planar section
            for (int j=j_min; j<=j_max; ++j) {

                if (i == 1 && j <= 0) { //skip slices that would be below the ground plane
                    continue;
                }

                QVector3D p = n * (generate_slices_spacing * j);

                PlanarSection new_sec;
                SetupPlanarSection(new_sec);
                new_sec.SetT(t);
                new_sec.SetN(n);
                new_sec.SetB(b);
                new_sec.SetP(p);


                //2c.1.  Compute the segment list with the plane intersected with triangles of the surface
                QList <QVector3D> segments;
                ComputeTemplateCut(n, p, segments);

                if (segments.empty()) {
                    continue;
                }

                //2c.2  Project each line segment endpoint onto the plane
                QList <QVector2D> segments_2d;
                for (int k=0; k<segments.size(); ++k) {
                    segments_2d.push_back(new_sec.GetPoint2D(segments[k]));
                }

                //qDebug() << "segments in 2d" << segments_2d;

                //2c.2.  Combine the segments together into one or more curves
                QList <QList <QVector2D> > curves;
                GLutils::CurvesFromLineSegments(segments_2d, curves);

                if (curves.empty()) {
                    continue;
                }

                BezierCurve bez_curve;
                bez_curve.SetPointsFromPolyline(curves.first());
                new_sec.SetCurve(0, bez_curve);

                //2c.3.  Assign the planar section the computed curves as the input polylines
                /*
                for (int k=0; k<curves.size(); ++k) {
                    new_sec.SketchSetCurve(k, curves[k]);
                    new_sec.Update(k, section_error_tolerance);

                    if (k < curves.size()-1) {
                        new_sec.AddNewCurve();
                    }

                }
                */

                //new_sec.CreateCircle(QVector2D(0,0), 1.0f);
                new_sec.UpdateCurveTrisSlab();

                generated_sections.push_back(new_sec);

            }

        }
    }

        break;


    case GENSTATE_LINEAR:
        {

            if (generate_selections.size() < 2) {
                qDebug() << "GLWidget::DoGenerateLinear() - Warning, 2 sections need to be selected";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Linear Error", "Can't Accept Generate Linear - 2 sections need to be selected.");

                return false;
            }

            //1.  test to make sure they differ
            const int ind1 = generate_selections.front();
            const int ind2 = generate_selections.back();

            if (ind1 == ind2) {
                qDebug() << "GLWidget::DoGenerateLinear() - Last two selections the same, aborting";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Linear Error", "Can't Accept Generate Linear - Last two selections the same.");

                return false;
            }

            //2.  grab the set of curves that should be copied, this involves getting all children
            //    of the revolving one's children (relative to tree with root at base planar section)
            QVector <QVector <bool> > graph;
            QList <QList <int> > cycles;
            PlanarSection::ComputeIntersectionGraph(sections, graph);

            QList <int> branch;
            Tree tree;
            tree.CreateFromGraph(graph, ind1);
            tree.GetBranch(ind2, branch);

            if (branch.contains(ind1)) { //node shouldn't have parent in its branch, this means a cycle, just use the 1 node
                branch.clear();
                branch.push_back(ind2);
            }

            //3.  compute lines of intersection
            PlanarSection & spine_section = sections[ind1];
            PlanarSection & copy_section = sections[ind2];

            //compute line of intersection
            QList <QVector3D> intersects;
            QList <bool> intersects_which;
            spine_section.GetContourIntersections(copy_section, intersects, intersects_which);

            if (intersects.size() < 2) {
                qDebug() << "GLWidget::DoGenerateLinear() - Aborting, need at least 2 contour intersections.";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Linear Error", "Can't Accept Generate Linear - Section 1 and Section 2 must intersect.");

                return false;
            }

            //compute 2d intersection ray points
            QVector2D isec_p1_2d = spine_section.GetPoint2D(intersects.first());
            QVector2D isec_p2_2d = spine_section.GetPoint2D(intersects.last());
            //compute 2d intersection ray direction
            QVector2D isec_dir_2d = spine_section.GetPoint2D(intersects.last()) - spine_section.GetPoint2D(intersects.first());

            //compute bounding region
            QVector2D copy_normal_2d = spine_section.GetPoint2D(copy_section.N()) - spine_section.GetPoint2D(QVector3D(0, 0, 0));
            float min_val, max_val;
            spine_section.GetBoundingInterval(copy_normal_2d, min_val, max_val);

            const float copy_val = QVector2D::dotProduct(copy_normal_2d, spine_section.GetPoint2D(copy_section.P()));

            //we need to perform a test to determine which side of the contour to attach to consistently when we cast rays
            //create a ray on the plane of spine_section
            QList <QVector2D> contour_isecs;
            spine_section.ComputeContourLineIntersects(isec_p1_2d, isec_p2_2d-isec_p1_2d, contour_isecs);
            //which contour_isec is closer to intersects[1]? // TODO: make this more robust, this choice is arbitrary, and depends on spline before point being defined
            const float contour_first_len = (spine_section.GetPoint3D(contour_isecs.first()) - intersects[1]).length();
            const float contour_last_len = (spine_section.GetPoint3D(contour_isecs.last()) - intersects[1]).length();
            const float contour_scale = (spine_section.GetPoint3D(contour_isecs.last()) - spine_section.GetPoint3D(contour_isecs.first())).length();
            const bool use_first_contour_isec = (contour_first_len < contour_last_len);

            //compute series of 2D intersections on spine's plane, along direction ld (projected to 2D)
            const int nBefore = -(int((copy_val - min_val) / generate_linear_spacing) + 1);
            const int nAfter = int((max_val - copy_val) / generate_linear_spacing) + 1;

            const QVector3D oldbasis_t = copy_section.T();
            const QVector3D oldbasis_n = copy_section.N();
            const QVector3D oldbasis_b = copy_section.B();
            const QVector3D oldbasis_p = copy_section.P();

            //4.  add the duplicates in a linear arrangement
            for (int i=nBefore; i<=nAfter; ++i) {

                //we skip the section at the copy_section position
                if (i == 0) {
                    continue;
                }

                const float d = copy_val + generate_linear_spacing * i;

                const QVector3D new_plane_p1 = spine_section.GetPoint3D(copy_normal_2d * d);
                const QVector3D new_plane_p2 = spine_section.GetPoint3D(copy_normal_2d * d + isec_dir_2d);

                //create a ray on the plane of spine_section
                QVector2D ray_p1 = spine_section.GetPoint2D(new_plane_p1);
                QVector2D ray_p2 = spine_section.GetPoint2D(new_plane_p2);

                //cast the ray and get the intersections with spine_section's contour
                QList <QVector2D> contour_isecs;
                spine_section.ComputeContourLineIntersects(ray_p1, ray_p2-ray_p1, contour_isecs);

                //if hits, we duplicate the copy_section there
                if (!contour_isecs.empty()) {

                    //convert intersection points to 3d
                    QList <QVector3D> contour_isecs_3d;
                    for (int j=0; j<contour_isecs.size(); ++j) {
                        contour_isecs_3d.push_back(spine_section.GetPoint3D(contour_isecs[j]));
                    }

                    //sort the intersection points
                    GLutils::SortPointsAlongDirection3D(intersects.last() - intersects.first(), contour_isecs_3d);

                    //visual debugging
                    /*
                    for (int j=0; j<contour_isecs_3d.size(); ++j) {
                        markers.push_back(contour_isecs_3d[j]);
                        if (j == 0) {
                            markers_col.push_back(QVector3D(0, 0, 1));
                        }
                        else {
                            markers_col.push_back(QVector3D(1, 0, 1));
                        }
                    }
                    */

                    //create a new planar section copy here, using the intersection value... try both, and use the one
                    //whose line of intersection distance is closest the piece we are cloning
                    const float each_contour_scale = (contour_isecs_3d.last() - contour_isecs_3d.first()).length() / contour_scale;

                    const QVector3D newbasis_t = copy_section.T();
                    const QVector3D newbasis_n = copy_section.N();
                    const QVector3D newbasis_b = copy_section.B();
                    QVector3D newbasis_p;
                    if (use_first_contour_isec) { //when scaling Y a translation is required so the section will still intersect at the slit
                        newbasis_p = contour_isecs_3d.first() + (copy_section.P() - intersects[1]) + (intersects[1] - intersects[2])*(1.0f - each_contour_scale);
                    }
                    else {
                        newbasis_p = contour_isecs_3d.last() + (copy_section.P() - intersects[2]) + (intersects[2] - intersects[1])*(1.0f - each_contour_scale);
                    }

                    for (int j=0; j<branch.size(); ++j) {

                        PlanarSection & copy_sec = sections[branch[j]];

                        const float scale_x = (generate_linear_scalex) ? each_contour_scale : 1.0f;
                        const float scale_y = (generate_linear_scaley) ? each_contour_scale : 1.0f;

                        const QVector3D new_t = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                           newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                           copy_sec.T(), 1.0f, 1.0f, 1.0f);
                        const QVector3D new_n = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                           newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                           copy_sec.N(), 1.0f, 1.0f, 1.0f);
                        const QVector3D new_b = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                           newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                           copy_sec.B(), 1.0f, 1.0f, 1.0f);
                        const QVector3D new_p = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, oldbasis_p,
                                                                          newbasis_t, newbasis_n, newbasis_b, newbasis_p,
                                                                          copy_sec.P(), scale_x, 1.0f, scale_y);

                        PlanarSection new_sec = copy_sec;
                        SetupPlanarSection(new_sec);
                        new_sec.SetT(new_t);
                        new_sec.SetN(new_n);
                        new_sec.SetB(new_b);
                        new_sec.SetP(new_p);
                        new_sec.Scale(scale_x, scale_y);
                        new_sec.UpdateCurveTrisSlab();

                        generated_sections.push_back(new_sec);

                    }

                }

            }
        }
        break;


    case GENSTATE_BLEND:
        {

            if (generate_selections.size() < 3) {
                qDebug() << "GLWidget::DoGenerateBlend() - Warning, 3 sections need to be selected";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Blend Error", "Can't Accept Generate Blend - 3 sections need to be selected.");

                return false;
            }

            //2.  need to ensure that sections 2 and 3 are connected to section 1
            //  2a.  test that connectivity
            PlanarSection & section1 = sections[generate_selections[0]];
            PlanarSection & section2 = sections[generate_selections[1]];
            PlanarSection & section3 = sections[generate_selections[2]];

            if (!section1.IsIntersectingSection(section2) || !section1.IsIntersectingSection(section3)) {
                qDebug() << "GLWidget::DoGenerateBlend() - Warning, section1 does not intersect both section2 and section3";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Blend Error", "Can't Accept Generate Blend - Section 1 does not intersect both Section 2 and Section 3.");

                return false;
            }

            AddToUndoList(OP_GENERATE_BLEND);

            //3.  we need to find the interval of section 1's boundary to regularly add sections to
            BezierCurve curve;
            section1.GetCurveBetweenSections(section2, section3, curve);
            /*
            const QList <QVector2D> & samples = curve.Samples();
            markers2.clear();
            markers_col2.clear();
            for (int i=0; i<samples.size(); ++i) {
                markers2.push_back(section1.GetPoint3D(samples[i]));
                markers_col2.push_back(QVector3D(float(i)/float(samples.size()), 1, 1));
            }
            */

            QList <PlanarSection> new_sections;
            section1.GetSectionsAlongCurve(section2, curve, generate_blend_sections-1, new_sections);

            //4.  set up first and last sections so their normals point tangential to the curve (roughly) and they are centered on the curve
            //  4a.  translate P of the first and last so they are on the curve (and the interpolated ones are translated correctly)
            section2.MoveP(new_sections.first().P());
            section3.MoveP(new_sections.last().P());

            //  4b.  set the normal of the first and last so that they are facing along the curve
            if (QVector3D::dotProduct(section2.N(), new_sections.first().N()) < 0.0f) {
                section2.FlipN();
                section2.UpdateCurveTrisSlab();
            }
            if (QVector3D::dotProduct(section3.N(), new_sections.last().N()) < 0.0f) {
                section3.FlipN();
                section3.UpdateCurveTrisSlab();
            }
            if (QVector3D::dotProduct(section2.T(), new_sections.first().T()) < 0.0f) {
                for (int i=0; i<new_sections.size(); ++i) {
                    new_sections[i].FlipTB();
                }
            }

            //  4c.  just update them
            section2.UpdateCurveTrisSlab();
            section3.UpdateCurveTrisSlab();

            //5.  create a correspondence between section boundaries of 2 and 3
            //  5a.  correspondence depends on same # of control points.  subdivide section with more, splitting at longest segments
            BezierCurve & curve2 = section2.GetCurve(0);
            BezierCurve & curve3 = section3.GetCurve(0);

            while (curve2.GetNumControlPoints() < curve3.GetNumControlPoints()) {
                curve2.SubdivideLongestSegment();
                //qDebug() << "subdivided curve2 - pts" << curve2.GetNumControlPoints();
            }
            while (curve3.GetNumControlPoints() < curve2.GetNumControlPoints()) {
                curve3.SubdivideLongestSegment();
                //qDebug() << "subdivided curve3 - pts" << curve3.GetNumControlPoints();
            }

            //  5b.  now find the control point correspondences (an index offset mapping ctrl points from 2 to 3)
            int corresp_offset;
            bool corresp_forward;
            curve2.GetCurvePointCorrespondence(curve3, corresp_offset, corresp_forward);
            //qDebug() << "using correspondence" << corresp_offset << corresp_forward;

            //6.  March along and do the blending interpolation
            float theta1 = GLutils::SignedAngleBetweenRad(section2.N(), new_sections.first().N(), section2.T());
            float theta2 = GLutils::SignedAngleBetweenRad(section3.N(), new_sections.last().N(), section2.T());

            for (int i=1; i<generate_blend_sections-1; ++i) {

                //t takes on values (0, 1) and not [0, 1]
                const float t = float(i) / float(generate_blend_sections);

                BezierCurve blend_curve;
                curve2.InterpolateCurvePointCorrespondence(curve3, corresp_offset, corresp_forward, t, blend_curve);

                const float theta_t = theta1 * (1.0f - t) + theta2 * t;
                //qDebug() << i << theta1 << theta2 << theta_t << t;

                QVector3D rotate_n = GLutils::RotateVector(new_sections[i].N(), section2.T(), -theta_t);
                QVector3D rotate_b = GLutils::RotateVector(new_sections[i].B(), section2.T(), -theta_t);

                new_sections[i].SetCurve(0, blend_curve);
                new_sections[i].SetN(rotate_n);
                new_sections[i].SetB(rotate_b);
                new_sections[i].SetSlabThickness(slab_thickness);
                new_sections[i].UpdateCurveTrisSlab();

                generated_sections.push_back(new_sections[i]);

            }
        }

        break;

    case GENSTATE_REVOLVE:
        {
            if (generate_selections.size() < 2) {
                qDebug() << "GLWidget::DoGenerateRevolve() - Warning, 2 sections need to be selected";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Revolve Error", "Can't Accept Generate Revolve - 2 sections need to be selected");

                return false;
            }

            //1.  need to ensure that section 1 is connected to section 2
            //  1a.  test that connectivity
            const int ind1 = generate_selections[0];
            const int ind2 = generate_selections[1];
            PlanarSection & section1 = sections[ind1];
            PlanarSection & section2 = sections[ind2];

            if (!section1.IsIntersectingSection(section2)) {
                qDebug() << "GLWidget::DoGenerateRevolve() - Warning, section1 does not intersect section2";

                UpdateAllTests();
                UpdateDraw();
                if(current_tool_state == TOOLSTATE_DEFAULT) // this is should only be shown during AcceptGenerate()
                    QMessageBox::warning ( this, "Generate Revolve Error", "Can't Accept Generate Revolve - Section 1 does not intersect Section 2.");

                return false;
            }

            //2.  grab the set of curves that should be copied, this involves getting all children
            //    of the revolving one's children (relative to tree with root at base planar section)
            QVector <QVector <bool> > graph;
            PlanarSection::ComputeIntersectionGraph(sections, graph);

            QList <int> branch;
            Tree tree;
            tree.CreateFromGraph(graph, ind1);
            tree.GetBranch(ind2, branch);

            if (branch.contains(ind1)) { //node shouldn't have parent in its branch, this means a cycle, just use the 1 node
                branch.clear();
                branch.push_back(ind2);
            }

            //3.  we need to find the interval of section 1's boundary to regularly add sections to
            BezierCurve curve;
            section1.GetCurveAroundSection(section2, curve);
            //section1.SetCurve(0, curve); //debugging

            QList <PlanarSection> new_sections;
            section1.GetSectionsAlongCurve(section2, curve, generate_revolve_sections, new_sections);

            //4.  set up first and last sections so their normals point tangential to the curve (roughly) and they are centered on the curve
            //  4a.  translate P of the first and last so they are on the curve (and the interpolated ones are translated correctly)
            section2.MoveP(new_sections.first().P());

            //  4b.  set the normal of the first and last so that they are facing along the curve
            if (QVector3D::dotProduct(section2.N(), new_sections.first().N()) < 0.0f) {
                section2.FlipN();
                section2.UpdateCurveTrisSlab();
                //qDebug() << "FLIPPED N2";
            }

            if (QVector3D::dotProduct(section2.T(), new_sections.first().T()) < 0.0f) {
                for (int i=0; i<new_sections.size(); ++i) {
                    new_sections[i].FlipTB();
                }
            }

            //  4c.  just update them
            section2.UpdateCurveTrisSlab();

            //5.  March along and add sections
            float theta = GLutils::SignedAngleBetweenRad(section2.N(), new_sections.first().N(), section2.T());

            const QVector3D oldbasis_t = section2.T();
            const QVector3D oldbasis_n = section2.N();
            const QVector3D oldbasis_b = section2.B();
            const QVector3D oldbasis_p = section2.P();

            for (int i=1; i<generate_revolve_sections; ++i) { //for each section sampled on the boundary...

                const QVector3D newbasis_t = new_sections[i].T();
                const QVector3D newbasis_n = GLutils::RotateVector(new_sections[i].N(), section2.T(), -theta);
                const QVector3D newbasis_b = GLutils::RotateVector(new_sections[i].B(), section2.T(), -theta);
                const QVector3D newbasis_p = new_sections[i].P();

                //qDebug() << "section2t" << section2.T() << "theta" << theta << "newsections[i].N" << new_sections[i].N() << "newsections[i].B" << new_sections[i].B();
                //qDebug() << "new section" << i << "t" << newbasis_t << "n" << newbasis_n << "b" << newbasis_b << "p" << newbasis_p;

                for (int j=0; j<branch.size(); ++j) { //and for every section which is part of the branch...

                    PlanarSection & copy_sec = sections[branch[j]];

                    const QVector3D new_t = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                       newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                       copy_sec.T(), 1.0f, 1.0f, 1.0f);
                    const QVector3D new_n = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                       newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                       copy_sec.N(), 1.0f, 1.0f, 1.0f);
                    const QVector3D new_b = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, QVector3D(0, 0, 0),
                                                                       newbasis_t, newbasis_n, newbasis_b, QVector3D(0, 0, 0),
                                                                       copy_sec.B(), 1.0f, 1.0f, 1.0f);
                    const QVector3D new_p = GLutils::GetVectorNewBasis(oldbasis_t, oldbasis_n, oldbasis_b, oldbasis_p,
                                                                       newbasis_t, newbasis_n, newbasis_b, newbasis_p,
                                                                       copy_sec.P(), 1.0f, 1.0f, 1.0f);

                    PlanarSection new_sec = copy_sec;
                    SetupPlanarSection(new_sec);
                    new_sec.SetT(new_t);
                    new_sec.SetN(new_n);
                    new_sec.SetB(new_b);
                    new_sec.SetP(new_p);
                    new_sec.UpdateCurveTrisSlab();
                    generated_sections.push_back(new_sec);

                }

            }
        }
        break;


    case GENSTATE_GRID:
    {

        if (generate_selections.size() < 1) {
            qDebug() << "GLWidget::DoGenerateGrid() - Warning, 1 section must be selected";

            UpdateAllTests();
            UpdateDraw();
            if (current_tool_state == TOOLSTATE_DEFAULT) { // this is should only be shown during AcceptGenerate()
                QMessageBox::warning ( this, "Generate Grid Error", "Can't Accept Generate Grid - 1 section must be selected");
            }

            return false;
        }

        //qDebug() << sections.size() << generate_selections.size() << generate_selections[0];
        PlanarSection & section1 = sections[generate_selections[0]];

        //TODO: this stuff should eventually be moved into planarsection class, and only parameters passed
        const float offset = -0.01f;

        //get bounding box
        QVector2D min_bb, max_bb;
        section1.GetBoundingBox2D(min_bb, max_bb);

        QList <PlanarSection> split_sections;
        split_sections.push_back(section1);

        QList <QVector2D> split_pts;
        QList <QVector2D> split_dirs;

        for (float x=min_bb.x()+(generate_grid_sizex*offset); x<=max_bb.x(); x += generate_grid_sizex) {
            split_pts.push_back(QVector2D(x, 0));
            split_dirs.push_back(QVector2D(0, 1));
        }

        for (float y=min_bb.y()+(generate_grid_sizey*offset); y<=max_bb.y(); y += generate_grid_sizey) {
            split_pts.push_back(QVector2D(0, y));
            split_dirs.push_back(QVector2D(1, 0));
        }

        //iterate through all cuts
        for (int i=0; i<split_pts.size(); ++i) {

            //this saves sections remaining following the cut
            QList <PlanarSection> split_sections_after_cut;

            for (int j=0; j<split_sections.size(); ++j) {

                QList <PlanarSection> each_split_sections;
                split_sections[j].SplitAlongLine(split_pts[i], split_dirs[i], each_split_sections);

                if (!each_split_sections.empty()) {
                    split_sections_after_cut += each_split_sections;
                }
                else {
                    split_sections_after_cut += split_sections[j];
                }

            }

            split_sections = split_sections_after_cut;

        }

        for (int i=0; i<split_sections.size(); ++i) {
            SetupPlanarSection(split_sections[i]);
        }

        //need to add more sections between neighbouring cuts
        //add a bunch of staples regularly, then, only keep the ones connected to at least 2 other planes


        QList <PlanarSection> staple_sections;
        //sections with normal which is vertical (wrt to the grid plane)
        for (float x=min_bb.x()+(generate_grid_sizex*(1.0f+offset)); x<=max_bb.x(); x += generate_grid_sizex) {
            for (float y=min_bb.y()+(generate_grid_sizey*(0.5+offset)); y<=max_bb.y(); y += generate_grid_sizey) {

                PlanarSection staple;
                SetupPlanarSection(staple);
                staple.SetP(section1.GetPoint3D(QVector2D(x, y)));
                staple.SetT(section1.N());
                staple.SetN(section1.B());
                staple.SetB(section1.T());
                staple.CreateSquare(generate_grid_staplesize);
                staple_sections.push_back(staple);

            }
        }

        //sections with normal which is horizontal (wrt to the grid plane)
        for (float y=min_bb.y()+(generate_grid_sizey*(1.0f+offset)); y<=max_bb.y(); y += generate_grid_sizey) {
            for (float x=min_bb.x()+(generate_grid_sizex*(0.5+offset)); x<=max_bb.x(); x += generate_grid_sizex) {

                PlanarSection staple;
                SetupPlanarSection(staple);
                staple.SetP(section1.GetPoint3D(QVector2D(x, y)));
                staple.SetT(section1.B());
                staple.SetN(section1.T());
                staple.SetB(section1.N());
                staple.CreateSquare(generate_grid_staplesize);
                staple_sections.push_back(staple);

            }
        }

        //for all staple sections now ensure they intersect at least 2 of the split sections
        for (int i=0; i<staple_sections.size(); ++i) {

            //count # of intersections
            int num_staple_split_ints = 0;
            for (int j=0; j<split_sections.size(); ++j) {

                if (staple_sections[i].IsIntersectingSection(split_sections[j])) {
                    ++num_staple_split_ints;
                    //qDebug() << i << "intersects" << j;
                }
            }

            //if less than 2, remove this staple
            //qDebug() << i << "intersects" << num_staple_split_ints << "sections";
            if (num_staple_split_ints < 2) {
                staple_sections.removeAt(i);
                --i;
            }
        }


        generated_sections += split_sections;
        generated_sections += staple_sections;
    }
        break;

    default:
        break;


    }

    //UpdateAllTests();
    UpdateDraw();
    return true;
}

void GLWidget::StartGenerateLinear()
{
    if (current_tool_state == TOOLSTATE_GENERATE) {
        CancelGenerate();
    }

    gen_state = GENSTATE_LINEAR;
    current_tool_state = TOOLSTATE_GENERATE;
    selected = -1;
    generate_selections.clear();
    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::StartGenerateBlend()
{
    if (current_tool_state == TOOLSTATE_GENERATE) {
        CancelGenerate();
    }

    gen_state = GENSTATE_BLEND;
    current_tool_state = TOOLSTATE_GENERATE;
    selected = -1;
    generate_selections.clear();
    UpdateAllTests();
    UpdateDraw();
}

void GLWidget::StartGenerateRevolve()
{
    if (current_tool_state == TOOLSTATE_GENERATE) {
        CancelGenerate();
    }

    gen_state = GENSTATE_REVOLVE;
    current_tool_state = TOOLSTATE_GENERATE;
    selected = -1;
    generate_selections.clear();
    UpdateAllTests();
    UpdateDraw();
}

void GLWidget::StartGenerateSlices()
{

    //sideWidget->setCurrentIndex(1);

    if (current_tool_state == TOOLSTATE_GENERATE) {
        CancelGenerate();
    }

    gen_state = GENSTATE_SLICES;
    current_tool_state = TOOLSTATE_GENERATE;
    selected = -1;
    generate_selections.clear();
    ShowGenerate(); //note: since no selections are required, perform operation immediately
    UpdateAllTests();
    UpdateDraw();
}

void GLWidget::StartGenerateGrid()
{
    if (current_tool_state == TOOLSTATE_GENERATE) {
        CancelGenerate();
    }

    gen_state = GENSTATE_GRID;
    current_tool_state = TOOLSTATE_GENERATE;
    selected = -1;
    generate_selections.clear();
    UpdateAllTests();
    UpdateDraw();
}



void GLWidget::SetSelectedAsRadial()
{
    if (IsSectionSelected() ) {
        sections[selected].MakeRadial(generate_radial_sectors,3);
    }
}

void GLWidget::RemoveRadial()
{
    if (IsSectionSelected() ) {
        sections[selected].SetRadial(false);
    }
}



void GLWidget::DoGenerateBranchingSetRoot()
{
    state = STATE_RECURSIVE_SETUP_SLOT;
}

void GLWidget::DoGenerateBranching()
{

    //sideWidget->setCurrentIndex(1);

    if (sections.size() < 2) {
        return;
    }

    AddToUndoList(OP_GENERATE_BRANCHING);

    QList <PlanarSection> last_sections = sections;


    //if there are n planar sections now, and root has k children, new tree consists of kn + 1 sections
    QList <PlanarSection> new_sections;
    QVector <bool> intersectors;
    PlanarSection::ComputeIntersectors(sections, 0, intersectors);

    int num_kids = 0;
    for (int i=0; i<intersectors.size(); ++i) {
        if (intersectors[i]) {
            ++num_kids;
        }
    }

    //qDebug() << "Recursing one level: " << sections.size() << "sections -> " << sections.size() * num_kids + 1 << "sections";

    for (int i=1; i<intersectors.size(); ++i) {

        //only add subtrees for immediate children of root node
        if (!intersectors[i]) {
            continue;
        }

        QList <PlanarSection> * sections_to_use = &sections;

        //randomly do stuff
        //if (qrand() % 100 < 50) {
        //    sections_to_use = &last_sections;
        //}

        //copy all of these sections to form subtree at child's TNB frame, incorporating transformation from parent's to child's TNB frame
        for (int j=0; j<sections_to_use->size(); ++j) {

            QVector3D new_t = GLutils::GetVectorNewBasis(sections[0].T(), sections[0].N(), sections[0].B(), QVector3D(0, 0, 0),
                    sections[i].T(), sections[i].N(), sections[i].B(), QVector3D(0, 0, 0),
                    (*sections_to_use)[j].T(), 1.0f, 1.0f, 1.0f);
            QVector3D new_n = GLutils::GetVectorNewBasis(sections[0].T(), sections[0].N(), sections[0].B(), QVector3D(0, 0, 0),
                    sections[i].T(), sections[i].N(), sections[i].B(), QVector3D(0, 0, 0),
                    (*sections_to_use)[j].N(), 1.0f, 1.0f, 1.0f);
            QVector3D new_b = GLutils::GetVectorNewBasis(sections[0].T(), sections[0].N(), sections[0].B(), QVector3D(0, 0, 0),
                    sections[i].T(), sections[i].N(), sections[i].B(), QVector3D(0, 0, 0),
                    (*sections_to_use)[j].B(), 1.0f, 1.0f, 1.0f);            
            QVector3D new_p = GLutils::GetVectorNewBasis(sections[0].T(), sections[0].N(), sections[0].B(), sections[0].P(),
                    sections[i].T(), sections[i].N(), sections[i].B(), sections[i].P(),
                    (*sections_to_use)[j].P(), generate_branching_scalechild, generate_branching_scalechild, generate_branching_scalechild);

            PlanarSection new_section = (*sections_to_use)[j];
            SetupPlanarSection(new_section);
            new_section.SetT(new_t);
            new_section.SetN(new_n);
            new_section.SetB(new_b);
            new_section.SetP(new_p);
            new_section.Scale(generate_branching_scalechild, generate_branching_scalechild);            
            new_section.UpdateCurveTrisSlab();
            new_sections.push_back(new_section);

        }

    }

    new_sections.push_front(sections.first());
    last_sections = sections;
    sections = new_sections;

    //qDebug() << "Done.  Sections:" << sections.size();

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::CreateSurfacePatches()
{
    contour_graph.Clear();
    contour_graph.Create(sections);
    UpdateDraw();
}

void GLWidget::DeleteSurfacePatches()
{
    contour_graph.Clear();
    UpdateDraw();
}

void GLWidget::StartDimensioningTool()
{

    state = STATE_DIMENSIONING_FIRST;

}

void GLWidget::DoPhysicsTest()
{

    if (!do_physics_test) {
        return;
    }

    physics.ClearProblem();

    //TODO/NOTE: physics expects everything in m (metres)

    // first add each rigid body
    double total_mass = 0.0;
    for (int i=0; i<sections.size(); ++i) {

        QVector3D mass_cent = sections[i].GetCentroid3D(); //centre of mass
        //double mass = fabsf(sections[i].SignedArea()) * sections[i].SlabThickness();
        double mass = physics_material_density *
                fabs(sections[i].SignedArea() * metres_per_unit * metres_per_unit) *
                sections[i].SlabThickness() * metres_per_unit;

        total_mass += mass;

        double mass_cent_array[3];
        mass_cent_array[0] = mass_cent.x() * metres_per_unit;
        mass_cent_array[1] = mass_cent.y() * metres_per_unit;
        mass_cent_array[2] = mass_cent.z() * metres_per_unit;

        QList <QVector3D> contact_pts; //contact points
        sections[i].GetXZPlanePoints(contact_pts);

        std::vector<double> contact_pts_stlvector;
        contact_pts_stlvector.resize(contact_pts.size() * 3);
        for (int j=0; j<contact_pts.size(); ++j) {
            contact_pts_stlvector[j*3+0] = contact_pts[j].x() * metres_per_unit;
            contact_pts_stlvector[j*3+1] = contact_pts[j].y() * metres_per_unit;
            contact_pts_stlvector[j*3+2] = contact_pts[j].z() * metres_per_unit;
        }

        physics.AddRigidBody(mass_cent_array, mass, contact_pts_stlvector);

        //add in the weights
        for (int j=0; j<sections[i].GetNumWeights(); ++j) {

            QVector3D pos = sections[i].GetWeightPosition(j);
            QVector3D force = sections[i].GetWeightForce(j);

            double mass_pos[3];
            double mass_force[3];

            mass_pos[0] = pos.x() * metres_per_unit;
            mass_pos[1] = pos.y() * metres_per_unit;
            mass_pos[2] = pos.z() * metres_per_unit;

            mass_force[0] = force.x();
            mass_force[1] = force.y();
            mass_force[2] = force.z();

            physics.AddWeightToRigidBody(mass_pos, mass_force);
        }

    }

    //qDebug() << "Total weight of model: " << total_mass << " kg";

    //then add each "plate" (each plate corresponds to a rigid body)
    for (int i=0; i<sections.size(); ++i) {

        double p[3];
        double t[3];
        double n[3];
        double b[3];

        QVector3DToArray(sections[i].P(), p);
        QVector3DToArray(sections[i].T(), t);
        QVector3DToArray(sections[i].N(), n);
        QVector3DToArray(sections[i].B(), b);

        p[0] *=  metres_per_unit;
        p[1] *=  metres_per_unit;
        p[2] *=  metres_per_unit;

        const double rho = 1.0;

        const QList <QVector2D> & samples = sections[i].GetCurve(0).Samples();

        std::vector<double> pts;
        pts.resize(samples.size() * 2);

        for (int j=0; j<samples.size(); ++j) {
            pts[j*2+0] = samples[j].x() * metres_per_unit;
            pts[j*2+1] = samples[j].y() * metres_per_unit;
        }

        physics.AddPlate(p, t, n, b, sections[i].SlabThickness() * metres_per_unit, rho, pts);

    }

    // then add the joint info
    for (int i=0; i<sections.size(); ++i) {
        for (int j=i+1; j<sections.size(); ++j) {

            //int index, other_index;
            QList <QVector2D> my_slots;
            QList <QList <QVector2D> > my_slots_rect;

            //sections[i].GetSlots(sections[j], index, other_index, my_slots, my_slots_rect);
            sections[i].GetSlots(sections[j], my_slots, my_slots_rect);

            if (my_slots.size() >= 2) {

                //QVector3D position = sections[i].GetPoint3D((my_slots[0] + my_slots[1]) * 0.5);

                double position_array[3];
                double jp0[4];
                double jp1[4];

                QVector3D position = sections[i].GetPoint3D(my_slots[1]);
                QVector2D p0_pt1 = my_slots[1];
                QVector2D p0_pt2 = my_slots[0];
                QVector2D p1_pt1 = sections[j].GetPoint2D(sections[i].GetPoint3D(my_slots[1]));
                QVector2D p1_pt2 = sections[j].GetPoint2D(sections[i].GetPoint3D(my_slots[1]-my_slots[0]+my_slots[1]));

                position_array[0] = position.x() * metres_per_unit;
                position_array[1] = position.y() * metres_per_unit;
                position_array[2] = position.z() * metres_per_unit;

                jp0[0] = p0_pt1.x() * metres_per_unit;
                jp0[1] = p0_pt1.y() * metres_per_unit;
                jp0[2] = p0_pt2.x() * metres_per_unit;
                jp0[3] = p0_pt2.y() * metres_per_unit;

                jp1[0] = p1_pt1.x() * metres_per_unit;
                jp1[1] = p1_pt1.y() * metres_per_unit;
                jp1[2] = p1_pt2.x() * metres_per_unit;
                jp1[3] = p1_pt2.y() * metres_per_unit;

                physics.AddJoint(position_array, jp0, jp1, i, j);
                //qDebug() << "adding joint at" << position << "between" << i << j;

            }
            else {
                //qDebug() << "tested" << i << j << "size was" << my_slots.size();
            }

        }
    }  

    physics.CutSlit();
    physics.Solve();

}

void GLWidget::ResetCamera()
{
    cam.SetEye(default_lookat + QVector3D(0, 0, cam_lookat_distance)); //looking along z axis
    cam.SetLookAt(default_lookat);
    cam.SetCamWidth(10.0f);
}

void GLWidget::DoPhysicsAddWeight()
{

    if (last_selected.empty() ) {
        return;
    }

    const int last = last_selected.last();

    if (last < 0 || last >= sections.size()) {
        return;
    }

    AddToUndoList(OP_ADD_WEIGHT);

    PlanarSection & section = sections[last_selected.last()];
    section.AddWeight(section.Centroid(), physics_new_weight_mass);

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::DoPhysicsRemoveWeights()
{

    if (last_selected.empty() ) {
        return;
    }

    const int last = last_selected.last();

    if (last < 0 || last >= sections.size()) {
        return;
    }

    if (sections[last].GetNumWeights() > 0) {

        AddToUndoList(OP_DELETE_WEIGHTS);
        sections[last].RemoveWeights();
        UpdateAllTests();
        UpdateDraw();

    }

}

void GLWidget::SetPhysicsNewWeightMass(const double d)
{
    physics_new_weight_mass = d;
    new_weight_label->setText(QString::number(physics_new_weight_mass) + QString(" kilograms"));
}

void GLWidget::SetPhysicsMaterialDensity(const double d)
{
    physics_material_density = d;

    UpdateAllTests();
    //updateGL();

}

void GLWidget::SetPhysicsMaximumStress(const double d)
{
    physics_max_stress = d * 1000000.0; //conversion to MPa (millions of Pascals)
    physics.SetMaximumStress(physics_max_stress);

    UpdateAllTests();
    //updateGL();
}

void GLWidget::DoGenerateMakeCircle()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_GENERATE_MAKE_CIRCLE);

    //1.  get existing points dimensions
    QVector2D min_v;
    QVector2D max_v;
    sections[selected].GetBoundingBox2D(min_v, max_v);

    //2.  parameterize circle based on dimensions
    const QVector2D centre = (min_v + max_v) * 0.5f;
    const float rad = (max_v.x() - min_v.x() + max_v.y() - min_v.y()) * 0.25f;
    sections[selected].CreateCircle(centre, rad);

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::DoGenerateMakeRectangle()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_GENERATE_MAKE_RECTANGLE);

    //1.  get existing points dimensions
    QVector2D min_v;
    QVector2D max_v;
    sections[selected].GetBoundingBox2D(min_v, max_v);

    //2.  parameterize rectangle based on dimensions
    sections[selected].CreateRectangle(min_v, max_v);

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::DoGenerateMakeRadialHole()
{

    if (!IsSectionSelected()) {
        return;
    }

    AddToUndoList(OP_GENERATE_MAKE_RADIAL_HOLE);

    //1.  get existing points dimensions
    QVector2D min_v;
    QVector2D max_v;
    sections[selected].GetBoundingBox2D(min_v, max_v);

    //2.  parameterize circle based on dimensions
    const QVector2D centre = (min_v + max_v) * 0.5f;
    const float rad = (max_v.x() - min_v.x() + max_v.y() - min_v.y()) * 0.25f;
    float small_radii[9];
    for (int i=0; i<9; ++i) {
        small_radii[i] = generate_radial_params[i] * 0.5f;
    }
    sections[selected].CreateRadialHoles(centre, rad, generate_radial_sectors, small_radii);

    UpdateAllTests();
    UpdateDraw();

}

void GLWidget::DoGenerateSurfaceFacets()
{

    //qDebug() << "GLWidget::DoGenerateSurfaceFacets()";
    if (template_verts.empty() || template_faces.empty()) {
        return;
    }

    AddToUndoList(OP_GENERATE_SURFACE_FACETS);

    sections.clear();

    //1.  iterate over each face, making a planar face
    for (int i=0; i<template_poly_faces.size(); ++i) {

        //qDebug() << "1.  Processing polyface" << i;
        PlanarSection new_section;
        SetupPlanarSection(new_section);

        //1a.  Set up the section's space
        const QList <int> & face = template_poly_faces[i];
        const QVector3D & v0 = template_verts[face[0]];
        const QVector3D & v1 = template_verts[face[1]];
        const QVector3D & v2 = template_verts[face[2]];

        QVector3D new_n = QVector3D::crossProduct((v1-v0).normalized(), (v2-v0).normalized()).normalized();
        QVector3D new_t = GLutils::GetOrthoVec(new_n);
        QVector3D new_b = QVector3D::crossProduct(new_n, new_t);

        new_section.SetT(new_t);
        new_section.SetN(new_n);
        new_section.SetB(new_b);

        //qDebug() << "Adding basis for polyface" << i << "with lengths" << new_n.length() << new_t.length() << new_b.length();

        //1b. Set up the section's 3D point (centroid of the poly verts)
        QVector3D new_p(0,0,0);
        for (int j=0; j<face.size(); ++j) {
            new_p += template_verts[face[j]];
        }
        new_p /= face.size();
        //also offset it by half the material thickness (so adjacent sections don't intersect one another)
        //new_p += new_n * slab_thickness * 0.5f;
        new_section.SetP(new_p);

        //1c.  Add control points to the section from the poly verts
        BezierCurve & curve = new_section.GetCurve(0);
        //curve.AddPoint(new_section.GetPoint2D(template_verts[face[0]])); //add initial point

        for (int j=0; j<face.size(); ++j) {

            const int i0 = j;
            const int i1 = (j+1) % face.size();

            QVector2D v0 = new_section.GetPoint2D(template_verts[face[i0]]);
            QVector2D v1 = new_section.GetPoint2D(template_verts[face[i1]]);

            //scale them in by the slab_thickness (allow a small gap so adjacent faces don't intersect, this is dependent on the slab_thickness)
            if (generate_surfacefacets_teeth) {

                //TODO: prevent intersections between adjacent sections
//                const float len_v0 = v0.length();
//                const float len_v1 = v1.length();
//                v0 *= (len_v0 - slab_thickness * 0.5f) / len_v0;
//                v1 *= (len_v1 - slab_thickness * 0.5f) / len_v1;

                float edge_len = (v1 - v0).length();
                const QVector2D x = (v1 - v0).normalized();
                const QVector2D y(-x.y(), x.x());

                int teethpairs = int(edge_len / slab_thickness / 8.0f);
                //const int teethpairs = 3; //TODO: the "ideal" toothpair should be a function of edge length and material thickness

                QVector2D v0_2 = v0;
                QVector2D v1_2 = v1;

                v0 += x * edge_len / float(teethpairs) / 2.0f;
                v1 -= x * edge_len / float(teethpairs) / 2.0f;

                if (j == 0) {
                    curve.AddPoint(v0_2);
                }
                curve.AddPoint(v0_2 * 0.75f + v0 * 0.25f);
                curve.AddPoint(v0_2 * 0.25f + v0 * 0.75f);
                curve.AddPoint(v0);

                edge_len = (v1 - v0).length();
                --teethpairs;

                for (int k=0; k<teethpairs; ++k) {

                    const float x0 = float(k)/float(teethpairs) * edge_len;
                    const float x1 = float(k+1)/float(teethpairs) * edge_len;
                    float y0 = slab_thickness;
                    float y1 = -slab_thickness;

                    const QVector2D p0(x0, y0);
                    const QVector2D p1(x0 + (x1-x0)*0.5f, y0);
                    const QVector2D p2(x0 + (x1-x0)*0.5f, y1);
                    const QVector2D p3(x1, y1);
                    const QVector2D p4(x1, 0);

                    const QVector2D _p0 = v0 + p0.x() * x + p0.y() * y;
                    const QVector2D _p1 = v0 + p1.x() * x + p1.y() * y;
                    const QVector2D _p2 = v0 + p2.x() * x + p2.y() * y;
                    const QVector2D _p3 = v0 + p3.x() * x + p3.y() * y;
                    const QVector2D _p4 = v0 + p4.x() * x + p4.y() * y;

                    const QVector2D _l = curve.Points().last();
                    curve.AddPoint(_l * 0.75f + _p0 * 0.25f);
                    curve.AddPoint(_l * 0.25f + _p0 * 0.75f);
                    curve.AddPoint(_p0);

                    curve.AddPoint(_p0 * 0.75f + _p1 * 0.25f);
                    curve.AddPoint(_p0 * 0.25f + _p1 * 0.75f);
                    curve.AddPoint(_p1);

                    curve.AddPoint(_p1 * 0.75f + _p2 * 0.25f);
                    curve.AddPoint(_p1 * 0.25f + _p2 * 0.75f);
                    curve.AddPoint(_p2);

                    curve.AddPoint(_p2 * 0.75f + _p3 * 0.25f);
                    curve.AddPoint(_p2 * 0.25f + _p3 * 0.75f);
                    curve.AddPoint(_p3);

                    if (k == teethpairs-1) {
                        curve.AddPoint(_p3 * 0.75f + _p4 * 0.25f);
                        curve.AddPoint(_p3 * 0.25f + _p4 * 0.75f);
                        curve.AddPoint(_p4);
                    }

                }

                curve.AddPoint(v1 * 0.75f + v1_2 * 0.25f);
                curve.AddPoint(v1 * 0.25f + v1_2 * 0.75f);
                curve.AddPoint(v1_2);

            }
            else {
                const float len_v0 = v0.length();
                const float len_v1 = v1.length();

                v0 = v0 * (len_v0 - slab_thickness * 1.5f) / len_v0;
                v1 = v1 * (len_v1 - slab_thickness * 1.5f) / len_v1;

                if (j == 0) {
                    curve.AddPoint(v0);
                }
                curve.AddPoint(v0 * 0.75f + v1 * 0.25f); //add 3 points for the Bezier segment
                curve.AddPoint(v0 * 0.25f + v1 * 0.75f);
                curve.AddPoint(v1);
            }

        }

        //1d.  Push the new planar section back
        //new_section.Scale(generate_branching_scalechild, generate_branching_scalechild);
        new_section.UpdateCurveTrisSlab();
        sections.push_back(new_section);

    }

    //this block of code adds the circular coins at edges to join adjacent faces
    if (!generate_surfacefacets_teeth) {

        QMap <QPair <int, int>, bool> edge_processed;

        //2.  iterate over each edge, adding a circular-shaped section at the midpoint
        for (int i=0; i<template_poly_faces.size(); ++i) {

            //qDebug() << "2.  Processing polyface" << i;
            const QList <int> & face = template_poly_faces[i];

            for (int j=0; j<face.size(); ++j) {

                const int i0 = face[j];
                const int i1 = face[(j+1) % face.size()];

                //2a.  ensure we only do this edge once
                QPair <int, int> edge_key(qMin(i0, i1), qMax(i0, i1));
                //qDebug() << "Considering face" << i << "edge" << j << "with key" << edge_key;
                if (edge_processed[edge_key]) {
                    continue;
                    //qDebug() << "Already processed.";
                }
                else {
                    edge_processed[edge_key] = true;
                    //qDebug() << "Processing for first time.";
                }

                //2b.  compute centroid point
                const QVector3D & v0 = template_verts[i0];
                const QVector3D & v1 = template_verts[i1];

                QVector3D new_n = (v1-v0).normalized();
                QVector3D new_t = GLutils::GetOrthoVec(new_n);
                QVector3D new_b = QVector3D::crossProduct(new_n, new_t);
                QVector3D new_p = (v0 + v1) * 0.5f;

                PlanarSection new_section;
                SetupPlanarSection(new_section);

                new_section.SetT(new_t);
                new_section.SetN(new_n);
                new_section.SetB(new_b);
                new_section.SetP(new_p);

                //2c.  create a disc-shaped piece for the boundary curve
                //new_section.CreateCircle(QVector2D(0,0), 1.0f);
                new_section.CreateCircle(QVector2D(0,0), slab_thickness * 4.0f);
                //new_section.UpdateCurveTrisSlab(); //CreateCircle() updates curvetrisslab
                sections.push_back(new_section);

            }

        }

    }

    UpdateAllTests();
    UpdateDraw();

}


void GLWidget::ToggleDrawDeformed()
{
    physics.SetDrawDeformed(!physics.GetDrawDeformed());
    deformed_checkbox->setChecked(physics.GetDrawDeformed());
}

void GLWidget::ToggleDrawSkeleton()
{
    physics.SetDrawSkeleton(!physics.GetDrawSkeleton());
    skeleton_checkbox->setChecked(physics.GetDrawSkeleton());
}

void GLWidget::ToggleDrawForce()
{
    physics.SetDrawForce(!physics.GetDrawForce());
    forces_checkbox->setChecked(physics.GetDrawForce());
}

void GLWidget::ToggleDrawSection()
{
    physics.SetDrawSection(!physics.GetDrawSection());
    section_checkbox->setChecked(physics.GetDrawSection());
}

void GLWidget::ToggleDrawMoments()
{
    physics.SetDrawSectionMoment(!physics.GetDrawSectionMoment());
    moment_checkbox->setChecked(physics.GetDrawSectionMoment());
}

void GLWidget::ToggleDrawTemplates()
{
    do_show_templates = !do_show_templates;
    toggleTemplatesButton->setChecked(do_show_templates); //keeps the menu option and button synced
}


void GLWidget::ToggleShowTNBFrames()
{
    do_show_tnb_frames = !do_show_tnb_frames;
    show_tnb_frames_checkbox->setChecked(do_show_tnb_frames);
}

void GLWidget::ToggleDoCyclesTest()
{
    do_cycles_test = !do_cycles_test;
    show_cycles_test_checkbox->setChecked(do_cycles_test);
    update_sections_disp_list = true;

    if (do_cycles_test) {
        DoCyclesConnectedTest();
    }
}

void GLWidget::ToggleDoStabilityTest()
{
    do_stability_test = !do_stability_test;
    show_stability_checkbox->setChecked(do_stability_test);
    update_sections_disp_list = true;

    if (do_stability_test) {
        DoStabilityTest();
    }
}

void GLWidget::ToggleShowShadow()
{
    do_show_shadow = !do_show_shadow;
    show_shadow_checkbox->setChecked(do_show_shadow);
    update_sections_disp_list = true;
}

void GLWidget::ToggleDoConnectedTest()
{
    do_connected_test = !do_connected_test;
    show_connectivity_checkbox->setChecked(do_connected_test);
    update_sections_disp_list = true;

    if (do_connected_test) {
        DoCyclesConnectedTest();
    }
}

void GLWidget::SetDrawDeformed(bool b)
{
    physics.SetDrawDeformed(b);
    deformed_checkbox->setChecked(b);
}

void GLWidget::SetDrawSkeleton(bool b)
{
    physics.SetDrawSkeleton(b);
    skeleton_checkbox->setChecked(b);
}

void GLWidget::SetDrawForce(bool b)
{
    physics.SetDrawForce(b);
    forces_checkbox->setChecked(b);
}

void GLWidget::SetDrawSection(bool b)
{
    physics.SetDrawSection(b);
    section_checkbox->setChecked(b);
}

void GLWidget::SetDrawSectionMoment(bool b)
{
    physics.SetDrawSectionMoment(b);
    moment_checkbox->setChecked(b);
}
