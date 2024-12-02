#-------------------------------------------------
#
# Project created by QtCreator 2012-11-10T11:13:21
#
#-------------------------------------------------

QT       += core gui opengl widgets webkit network webkitwidgets

TARGET = FlatFab
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    eigen.cpp \
    glutils.cpp \
    planarsection.cpp \
    beziercurve.cpp \
    pivotcamera.cpp \
    bezierfit.cpp \
    contourgraph.cpp \
    physics.cpp \
    physics_matrix3d.cpp \
    physics_vector3d.cpp \
    tree.cpp \
    transformwidget.cpp \
    triangulate2.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    eigen.h \
    glutils.h \
    planarsection.h \
    beziercurve.h \
    pivotcamera.h \
    bezierfit.h \
    contourgraph.h \
    physics.h \
    physics_matrix3d.h \
    physics_vector3d.h \
    tree.h \
    transformwidget.h \
    triangulate2.h

#linux specific build settings
unix:!macx:LIBS += -lGLU
unix:!macx:QMAKE_LFLAGS +=  '-Wl,-rpath,\'\$$ORIGIN/libs\'' #this sets the RPATH to "libs" local to the executable location
unix:!macx:INCLUDEPATH += /usr/include/eigen3

#windows specific build settings
win32:RC_FILE = flatfab.rc  #program icon


RESOURCES += \
    resources.qrc
