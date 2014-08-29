#ifndef CONTOURGRAPH_H
#define CONTOURGRAPH_H

#include <QtOpenGL>

#include "planarsection.h"
#include "glutils.h"

class ContourGraph
{

public:

    ContourGraph();

    void Create(const QList <PlanarSection> & sections);

    void AddVertex(const QVector3D & v);
    QVector3D Vertex(const int index) const;
    const QList <QVector3D> & Verts();
    int ClosestVertex(const QVector3D & p);
    void ClosestVertex(const QList <QVector3D> & ps, int & ps_index, int & vert_index);
    void AddEdge(const QPair <int, int> & edge, const QVector3D & norm, const QList <QVector3D> & edge_pts, const int edge_plane);
    int GetNumPatches();
    void Clear();    
    void Draw() const;
    void SaveToOBJFile(QTextStream & ofs);

private:

    //private methods
    void ComputeVerts(const QList <PlanarSection> & sections);
    void ComputeEdges(const QList <PlanarSection> & sections);
    void ComputeCycles(); //TODO: switch to tree implementation
    void ComputeCycleEdgePoints(const QList <PlanarSection> & sections);
    void ComputeCycleEdgeArcLengths();

    void PruneCycles(const QList <PlanarSection> & sections);

    QVector3D SampleCycleEdgeContour(const int c_ind, const int e_ind, const float f);

    void FitCoonsPatches();

    //private members
    QVector <QVector <bool> > plane_graph;

    //step 1 get the vertices
    QList <QVector3D> verts;

    //step 2 get the edges and edge info
    QList <QPair <int, int> > edges;
    QList <int> edge_resting_plane;
    QList <QVector3D> edge_normals;
    QList <QList <QVector3D> > edge_points;
    QVector <QVector <bool> > edge_matrix;

    //step 3 get the cycles (but only keep the ones that are patchable)
    int max_cycle_size;
    QList <QList <int> > cycles;
    QList <QList <int> > cycle_edges;

    //step 4 get edge parameterizations for making patches
    QList <QList <QVector3D> > cycle_edge_normals;
    QList <QList <QList <QVector3D> > > cycle_edge_points;
    QVector <QVector <float> > cycle_edge_lengths;
    QVector <QVector <QVector <float> > > cycle_edge_arclengths;

    //step 5 make the patches
    QVector <QVector <QVector <QVector3D> > > cycle_coons_patches;

};

#endif // CONTOURGRAPH_H
