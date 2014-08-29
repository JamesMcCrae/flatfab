#-------------------------------------------------
#
# Project created by QtCreator 2012-11-10T11:13:21
#
#-------------------------------------------------

QT       += core gui opengl widgets webkitwidgets

TARGET = FlatFab
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    eigen.cpp \
    glutils.cpp \
    triangle.c \
    triangulate.cpp \
    planarsection.cpp \
    beziercurve.cpp \
    pivotcamera.cpp \
    bezierfit.cpp \
    contourgraph.cpp \
    physics.cpp \
    physics_matrix3d.cpp \
    physics_vector3d.cpp \
    tree.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    eigen.h \
    glutils.h \
    triangle.h \
    triangulate.h \
    planarsection.h \
    beziercurve.h \
    pivotcamera.h \
    bezierfit.h \
    contourgraph.h \
    physics.h \
    physics_matrix3d.h \
    physics_vector3d.h \
    tree.h

unix:!macx:LIBS += -lGLU

unix:!macx:QMAKE_RPATHDIR += ./libs

# program icon (for Windows, anyway) Chris - can you figure out how to do it for Mac - perhaps the same/similar means?
win32:RC_FILE = flatfab.rc
