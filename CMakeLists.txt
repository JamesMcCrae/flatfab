cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0072 NEW)
project(FlatFab VERSION 0.8.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets Network WebEngineWidgets REQUIRED)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenGL REQUIRED)
find_library(GLU_LIBRARY GLU REQUIRED)

include_directories(include ${EIGEN3_INCLUDE_DIR})

set(SOURCES
    src/main.cpp
    src/mainwindow.cpp
    src/glwidget.cpp
    src/eigen.cpp
    src/glutils.cpp
    src/planarsection.cpp
    src/beziercurve.cpp
    src/pivotcamera.cpp
    src/bezierfit.cpp
    src/contourgraph.cpp
    src/physics.cpp
    src/physics_matrix3d.cpp
    src/physics_vector3d.cpp
    src/tree.cpp
    src/transformwidget.cpp
    src/triangulate2.cpp
)

# Windows-specific settings
if(WIN32)
    set_source_files_properties(flatfab.rc PROPERTIES HEADER_FILE_ONLY TRUE)
    list(APPEND SOURCES flatfab.rc)
endif()


# Add headers
set(HEADERS
    include/mainwindow.h
    include/glwidget.h
    include/eigen.h
    include/glutils.h
    include/planarsection.h
    include/beziercurve.h
    include/pivotcamera.h
    include/bezierfit.h
    include/contourgraph.h
    include/physics.h
    include/physics_matrix3d.h
    include/physics_vector3d.h
    include/tree.h
    include/transformwidget.h
    include/triangulate2.h
)

# Add resources
set(RESOURCES resources.qrc)

# Define the executable
add_executable(${PROJECT_NAME}
    ${SOURCES}
    ${HEADERS}
    ${RESOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} 
    Qt5::Core
    Qt5::Gui
    Qt5::OpenGL
    Qt5::Widgets
    Qt5::Network
    Qt5::WebEngineWidgets
    OpenGL::GL
    ${GLU_LIBRARY})

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}  # For binaries
)
