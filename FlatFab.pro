#-------------------------------------------------
#
# Project created by QtCreator 2012-11-10T11:13:21
#
#-------------------------------------------------

QT       += core gui opengl widgets webkit network webkitwidgets

TARGET = FlatFab
TEMPLATE = app


SOURCES += src/main.cpp\
    src/mainwindow.cpp \
    src/glwidget.cpp \
    src/eigen.cpp \
    src/glutils.cpp \
    src/planarsection.cpp \
    src/beziercurve.cpp \
    src/pivotcamera.cpp \
    src/bezierfit.cpp \
    src/contourgraph.cpp \
    src/physics.cpp \
    src/physics_matrix3d.cpp \
    src/physics_vector3d.cpp \
    src/tree.cpp \
    src/transformwidget.cpp \
    src/triangulate2.cpp

HEADERS  += include/mainwindow.h \
    include/glwidget.h \
    include/eigen.h \
    include/glutils.h \
    include/planarsection.h \
    include/beziercurve.h \
    include/pivotcamera.h \
    include/bezierfit.h \
    include/contourgraph.h \
    include/physics.h \
    include/physics_matrix3d.h \
    include/physics_vector3d.h \
    include/tree.h \
    include/transformwidget.h \
    include/triangulate2.h

#linux specific build settings
unix:!macx:LIBS += -lGLU
unix:!macx:QMAKE_LFLAGS +=  '-Wl,-rpath,\'\$$ORIGIN/libs\'' #this sets the RPATH to "libs" local to the executable location

#windows specific build settings
win32:RC_FILE = flatfab.rc  #program icon


RESOURCES += \
    resources.qrc
