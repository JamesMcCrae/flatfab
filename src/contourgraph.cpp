#include "contourgraph.h"

ContourGraph::ContourGraph() { max_cycle_size = 6; }

void ContourGraph::Create(const QList<PlanarSection> & sections)
{
    Clear();

    // there will be a list of plane-plane intersections, these are the contour
    // graph vertices and contours are the edges
    PlanarSection::ComputeIntersectionGraph(sections, plane_graph);

    // step 1
    ComputeVerts(sections);

    // step 2
    ComputeEdges(sections);

    // step 3
    ComputeCycles();
    PruneCycles(sections);

    // step 4
    ComputeCycleEdgePoints(sections);
    ComputeCycleEdgeArcLengths();

    // step 5
    FitCoonsPatches();

    qDebug() << "Created contour graph with" << verts.size() << "vertices,"
             << edges.size() << "edges, and" << cycles.size() << "cycles.";
}

void ContourGraph::AddVertex(const QVector3D & v) { verts.push_back(v); }

QVector3D ContourGraph::Vertex(const int index) const { return verts[index]; }

const QList<QVector3D> & ContourGraph::Verts() { return verts; }

int ContourGraph::ClosestVertex(const QVector3D & p)
{
    return GLutils::GetClosestPoint(verts, p);
}

void ContourGraph::ClosestVertex(const QList<QVector3D> & ps, int & ps_index,
                                 int & vert_index)
{
    GLutils::GetClosestPairOfPoints(ps, verts, ps_index, vert_index);
}

void ContourGraph::AddEdge(const QPair<int, int> & edge, const QVector3D & norm,
                           const QList<QVector3D> & edge_pts,
                           const int edge_plane)
{
    edges.push_back(edge);
    edge_normals.push_back(norm);
    edge_points.push_back(edge_pts);
    edge_matrix[edge.first][edge.second] = true;
    edge_matrix[edge.second][edge.first] = true;
    edge_resting_plane.push_back(edge_plane);
}

int ContourGraph::GetNumPatches() { return cycle_edge_points.size(); }

void ContourGraph::Clear()
{
    plane_graph.clear();
    verts.clear();
    edges.clear();
    edge_normals.clear();
    edge_points.clear();
    edge_matrix.clear();
    edge_resting_plane.clear();

    cycles.clear();
    cycle_edges.clear();
    cycle_edge_normals.clear();
    cycle_edge_points.clear();

    cycle_edge_lengths.clear();
    cycle_edge_arclengths.clear();

    cycle_coons_patches.clear();
}

void ContourGraph::Draw() const
{
    glEnable(GL_LIGHTING);

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < cycle_coons_patches.size(); ++i) {
        for (int s = 1; s < cycle_coons_patches[i].size(); ++s) {
            for (int t = 1; t < cycle_coons_patches[i][s].size(); ++t) {

                GLutils::ColorByIndex(i + (s + t) % 2);

                const QVector3D & p1 = cycle_coons_patches[i][s - 1][t - 1];
                const QVector3D & p2 = cycle_coons_patches[i][s][t - 1];
                const QVector3D & p3 = cycle_coons_patches[i][s][t];
                const QVector3D & p4 = cycle_coons_patches[i][s - 1][t];

                const QVector3D & n = QVector3D::crossProduct(
                    (p2 - p1).normalized(), (p4 - p1).normalized());

                glNormal3f(n.x(), n.y(), n.z());

                glVertex3f(p1.x(), p1.y(), p1.z());
                glVertex3f(p2.x(), p2.y(), p2.z());
                glVertex3f(p3.x(), p3.y(), p3.z());

                glVertex3f(p1.x(), p1.y(), p1.z());
                glVertex3f(p3.x(), p3.y(), p3.z());
                glVertex3f(p4.x(), p4.y(), p4.z());
            }
        }
    }
    glEnd();

    glDisable(GL_LIGHTING);    
}

void ContourGraph::SaveToOBJFile(QTextStream & ofs)
{
    // write out each patch's verts
    for (int i = 0; i < cycle_coons_patches.size(); ++i) {
        for (int s = 0; s < cycle_coons_patches[i].size(); ++s) {
            for (int t = 0; t < cycle_coons_patches[i][s].size(); ++t) {
                const QVector3D & p = cycle_coons_patches[i][s][t];
                ofs << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
        }
    }

    // write out each patches faces
    for (int i = 0; i < cycle_coons_patches.size(); ++i) {
        const int patch_size = cycle_coons_patches[i].size();
        // note the plus 1!  obj verts are 1-indexed!
        // assumed all pathces are square
        int base_patch_index = i * patch_size * patch_size + 1;

        ofs << "g surfacepatch" << i << "\n";
        for (int s = 1; s < patch_size; ++s) {
            for (int t = 1; t < patch_size; ++t) {
                const int index1 =
                    base_patch_index + (s - 1) * patch_size + (t - 1);
                const int index2 = base_patch_index + (s)*patch_size + (t - 1);
                const int index3 = base_patch_index + (s)*patch_size + (t);
                const int index4 =
                    base_patch_index + (s - 1) * patch_size + (t);

                ofs << "f " << index1 << " " << index2 << " " << index3 << " "
                    << index4 << "\n";
            }
        }
    }
}

void ContourGraph::ComputeVerts(const QList<PlanarSection> & sections)
{
    // every "true" in graph is a plane-plane intersection
    for (int i = 0; i < sections.size(); ++i) {
        for (int j = i + 1; j < sections.size(); ++j) {
            // skip if their planar sections don't intersect
            if (!plane_graph[i][j]) {
                continue;
            }

            // we now have to go through each contour,
            // and organize intervals of the contour into EDGES
            QList<QVector3D> isecs;
            QList<bool> isecs_which;
            sections[i].GetContourIntersections(sections[j], isecs,
                                                isecs_which);

            if (isecs.size() < 4) {
                continue;
            }

            for (int k = 0; k < isecs.size(); ++k) {
                if (isecs_which[k]) {
                    AddVertex(isecs[k]);
                }
            }
        }
    }
}

void ContourGraph::ComputeEdges(const QList<PlanarSection> & sections)
{
    if (verts.empty()) {
        return;
    }

    edge_matrix.resize(verts.size());
    for (int i = 0; i < verts.size(); ++i) {
        edge_matrix[i] = QVector<bool>(verts.size(), false);
    }

    // iterate through the section i's contour,
    // finding all the intersection points, and push them back
    for (int i = 0; i < sections.size(); ++i) {
        QList<int> split_index;
        QList<QVector3D> split_point;

        QList<QVector3D> contour_verts;
        sections[i].ContourVertices3D(contour_verts);

        for (int j = 0; j < contour_verts.size(); ++j) {
            const int ind1 = j;
            const int ind2 = (j + 1) % contour_verts.size();

            const QVector3D v1 = contour_verts[ind1];
            const QVector3D v2 = contour_verts[ind2];

            for (int k = 0; k < sections.size(); ++k) {
                // skip if their planar sections don't intersect
                if (!plane_graph[i][k]) {
                    continue;
                }

                // check for a contour-plane intersection with the current line
                // segment of the contour
                QVector3D intersect;
                if (GLutils::LineSegmentPlaneIntersection(
                        sections[k].P(), sections[k].N(), v1, v2, intersect)) {
                    // nice, we found a split point
                    split_index.push_back(j);
                    split_point.push_back(intersect);
                }
            }
        }

        // close it off (duplicate the last splitpoint)
        split_index.push_back(split_index.first());
        split_point.push_back(split_point.first());

        while (split_index.size() >= 2) {
            int cur_split_index = split_index.first();
            QVector3D cur_split_point = split_point.first();

            // take one off
            split_index.pop_front();
            split_point.pop_front();

            // generate our edge for this interval
            QList<QVector3D> each_edge;
            QPair<int, int> each_endpts;

            each_endpts.first = ClosestVertex(cur_split_point);
            each_endpts.second = ClosestVertex(split_point.first());

            const int start_index =
                (cur_split_index + 1) % contour_verts.size();
            const int end_index =
                (split_index.first() + 1) % contour_verts.size();
            int cur_index = start_index;

            // add points along the interval
            each_edge.push_back(cur_split_point);
            while (cur_index != end_index) {
                each_edge.push_back(contour_verts[cur_index]);
                cur_index = (cur_index + 1) % contour_verts.size();
            }
            each_edge.push_back(split_point.first());

            // add this subcontour to the graph (and the endpoint topology info)
            AddEdge(each_endpts, sections[i].N(), each_edge, i);
        }
    }
}

// TODO!  switch this to using the tree data structure and back edges
void ContourGraph::ComputeCycles()
{
    qDebug() << "Detecting all cycles in mesh...";

    QList<QList<int> > verts_walked;
    QList<QList<int> > edges_walked_list;
    QList<QSet<int> > edges_walked;

    // initialize a bunch of them to have 1 vert, no edges traversed
    for (int i = 0; i < edge_matrix.size(); ++i) {
        QList<int> each_cyc;
        each_cyc.push_back(i);
        verts_walked.push_back(each_cyc);

        QSet<int> each_walk;
        edges_walked.push_back(each_walk);

        QList<int> each_walked_list;
        edges_walked_list.push_back(each_walked_list);
    }

    // we do an iteration on our spreaders
    for (int i = 0; i < verts_walked.size(); ++i) {

        int first_vert = verts_walked[i].first();
        int end_vert = verts_walked[i].last();

        // is it a cycle?  if so we keep it
        if (!edges_walked[i].empty() && first_vert == end_vert) {
            cycles.push_back(verts_walked[i]);
            cycle_edges.push_back(edges_walked_list[i]);

            continue;
        }

        // for each spreader, we see the verts we can get to from here, and if
        // we took that edge already
        for (int j = 0; j < edges.size(); ++j) {
            // only follow edges of graph
            if (edges[j].first != end_vert && edges[j].second != end_vert) {
                continue;
            }

            int next_vert =
                (edges[j].first == end_vert) ? edges[j].second : edges[j].first;

            // new edge index
            const int new_edge_key = j;

            // only follow edge if we didn't yet
            if (!edges_walked[i].contains(new_edge_key)) {
                // we continue along, trying this edge
                QList<int> new_verts_walked = verts_walked[i];
                QSet<int> new_edges_walked = edges_walked[i];
                QList<int> new_edges_walked_list = edges_walked_list[i];

                new_verts_walked.push_back(next_vert);
                new_edges_walked.insert(new_edge_key);
                new_edges_walked_list.push_back(new_edge_key);

                // make sure we are not processing this exact edge list
                bool edgeset_not_yet_walked = true;
                for (int k = 0; k < edges_walked.size(); ++k) {
                    if (new_edges_walked == edges_walked[k]) {
                        edgeset_not_yet_walked = false;
                        break;
                    }
                }

                // add it to the list to process
                if (edgeset_not_yet_walked &&
                    new_edges_walked.size() <= max_cycle_size) {
                    verts_walked.push_back(new_verts_walked);
                    edges_walked.push_back(new_edges_walked);
                    edges_walked_list.push_back(new_edges_walked_list);
                }
            }
        }
    }
}

void ContourGraph::ComputeCycleEdgePoints(const QList<PlanarSection> & sections)
{
    for (int i = 0; i < cycle_edges.size(); ++i) {
        QList<QVector3D> each_cycle_edge_normals;
        QList<QList<QVector3D> > each_cycle_edge_points;

        for (int j = 0; j < cycle_edges[i].size(); ++j) {
            each_cycle_edge_normals.push_back(edge_normals[cycle_edges[i][j]]);

            int v1 = cycles[i][j];
            int v2 = cycles[i][(j + 1) % cycles[i].size()];

            QPair<int, int> & edge = edges[cycle_edges[i][j]];
            QList<QVector3D> & edge_pts = edge_points[cycle_edges[i][j]];

            if (edge.first == v1 && edge.second == v2) {
                each_cycle_edge_points.push_back(edge_pts);

            } else {
                // do it backwards
                QList<QVector3D> edge_backward;

                for (int l = 0; l < edge_pts.size(); ++l) {
                    edge_backward.push_back(edge_pts[edge_pts.size() - 1 - l]);
                }

                each_cycle_edge_points.push_back(edge_backward);
            }
        }

        // make edge normals consistent
        for (int j = 0; j < each_cycle_edge_normals.size(); ++j) {
            const int cur_edge_plane = edge_resting_plane[cycle_edges[i][j]];
            const int ind = (j + 1) % each_cycle_edge_normals.size();
            const QList<QVector3D> & edge = each_cycle_edge_points[ind];

            const QVector3D & edge_p = edge[edge.size() / 2];

            if (QVector3D::dotProduct(edge_p - sections[cur_edge_plane].P(),
                                      each_cycle_edge_normals[j]) < 0.0f) {
                each_cycle_edge_normals[j] = -each_cycle_edge_normals[j];
            }
        }

        cycle_edge_normals.push_back(each_cycle_edge_normals);
        cycle_edge_points.push_back(each_cycle_edge_points);
    }

    // for our cycles with just 2 edges, we need to split them both to create a
    // 4-cycle, for a "coons patch formulation"
    for (int i = 0; i < cycle_edge_points.size(); ++i) {
        // cycle is fine, leave it
        if (cycle_edge_points[i].size() == 4) {
            continue;
        }

        QList<QVector3D> new_cycle_edge_normals;
        QList<QList<QVector3D> > new_cycle_edge_points;
        for (int j = 0; j < 4; ++j) {
            QList<QVector3D> each_edge_pts;
            new_cycle_edge_points.push_back(each_edge_pts);
        }

        if (cycle_edge_points[i].size() == 2) {
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][0]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][0]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][1]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][1]);

            for (int j = 0; j < cycle_edge_points[i].size(); ++j) {
                for (int k = 0; k < cycle_edge_points[i][j].size(); ++k) {
                    const QVector3D & p = cycle_edge_points[i][j][k];

                    if (k < cycle_edge_points[i][j].size() / 2) {
                        new_cycle_edge_points[j * 2 + 0].push_back(p);
                    } else {
                        new_cycle_edge_points[j * 2 + 1].push_back(p);
                    }
                }
            }

        } else if (cycle_edge_points[i].size() == 3) {
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][0]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][1]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][2]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][2]);

            new_cycle_edge_points[0] = cycle_edge_points[i][0];
            new_cycle_edge_points[1] = cycle_edge_points[i][1];

            for (int k = 0; k < cycle_edge_points[i][2].size(); ++k) {
                const QVector3D & p = cycle_edge_points[i][2][k];

                if (k < cycle_edge_points[i][2].size() / 2) {
                    new_cycle_edge_points[2].push_back(p);
                } else {
                    new_cycle_edge_points[3].push_back(p);
                }
            }

        } else if (cycle_edge_points[i].size() >= 5) {
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][0]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][1]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][2]);
            new_cycle_edge_normals.push_back(cycle_edge_normals[i][3]);

            new_cycle_edge_points[0] = cycle_edge_points[i][0];
            new_cycle_edge_points[1] = cycle_edge_points[i][1];
            new_cycle_edge_points[2] = cycle_edge_points[i][2];

            for (int j = 3; j < cycle_edge_points[i].size(); ++j) {
                for (int k = 0; k < cycle_edge_points[i][j].size(); ++k) {
                    const QVector3D & p = cycle_edge_points[i][j][k];
                    new_cycle_edge_points[3].push_back(p);
                }
            }
        }

        cycle_edge_normals[i] = new_cycle_edge_normals;
        cycle_edge_points[i] = new_cycle_edge_points;
    }
}

void ContourGraph::ComputeCycleEdgeArcLengths()
{
    // compute edge lengths and arclengths
    cycle_edge_lengths.resize(cycle_edge_points.size());
    cycle_edge_arclengths.resize(cycle_edge_points.size());
    for (int i = 0; i < cycle_edge_points.size(); ++i) {
        // init edge lengths to zero for each cycle
        cycle_edge_lengths[i] =
            QVector<float>(cycle_edge_points[i].size(), 0.0f);
        cycle_edge_arclengths[i].resize(cycle_edge_points[i].size());

        // for each cycle's edge
        for (int j = 0; j < cycle_edge_points[i].size(); ++j) {
            cycle_edge_arclengths[i][j] =
                QVector<float>(cycle_edge_points[i][j].size(), 0.0f);

            // for each cycle's edge's points
            for (int k = 1; k < cycle_edge_points[i][j].size(); ++k) {
                // add segment length to edge length
                const float each_len = (cycle_edge_points[i][j][k] -
                                        cycle_edge_points[i][j][k - 1])
                                           .length();
                cycle_edge_lengths[i][j] += each_len;

                cycle_edge_arclengths[i][j][k] = cycle_edge_lengths[i][j];
            }

            // now we need to normalize arc lengths by total length of curve
            for (int k = 1; k < cycle_edge_points[i][j].size(); ++k) {
                cycle_edge_arclengths[i][j][k] /= cycle_edge_lengths[i][j];
            }
        }
    }

    qDebug() << cycle_edge_points.size() << cycle_edge_arclengths.size()
             << "these should be equal";
}

void ContourGraph::PruneCycles(const QList<PlanarSection> & sections)
{
    qDebug() << "Cycles before pruning:" << cycles.size() << "cycles.";

    for (int i = 0; i < cycles.size(); ++i) {
        bool prune_cycle = false;

        // pruning rule 1: alternate plane for every edge
        for (int j = 0; j < cycle_edges[i].size(); ++j) {
            int e1 = cycle_edges[i][j];
            int e2 = cycle_edges[i][(j + 1) % cycle_edges[i].size()];

            // resting plane between two consecutive edges must be different
            if (edge_resting_plane[e1] == edge_resting_plane[e2]) {
                prune_cycle = true;
            }
        }

        // pruning rule 2: remain within union of halfspaces defined by planes
        // taken
        for (int j = 0; j < cycle_edges[i].size(); ++j) {
            int edge_plane1 = edge_resting_plane[cycle_edges[i][j]];

            // for all edges not part of edge plane, compute their halfspace
            // side (they should be all negative or all positive)
            bool all_dot_prod_positive = true;
            bool all_dot_prod_negative = true;

            for (int k = 0; k < cycle_edges[i].size(); ++k) {
                int edge_plane2 = edge_resting_plane[cycle_edges[i][k]];

                if (edge_plane1 == edge_plane2) {
                    continue;
                }

                QList<QVector3D> & other_edge_pts =
                    edge_points[cycle_edges[i][k]];

                QVector3D edge_midpt =
                    other_edge_pts[other_edge_pts.size() / 2];

                const float each_dot_prod = QVector3D::dotProduct(
                    sections[edge_plane1].N(),
                    edge_midpt - sections[edge_plane1].P());

                if (each_dot_prod > 0.0f) {
                    all_dot_prod_negative = false;
                } else if (each_dot_prod < 0.0f) {
                    all_dot_prod_positive = false;
                }
            }

            if (!all_dot_prod_positive && !all_dot_prod_negative) {
                prune_cycle = true;
            }
        }

        // prune
        if (prune_cycle) {
            cycles.removeAt(i);
            cycle_edges.removeAt(i);
            cycle_edge_points.removeAt(i);
            --i;
        }
    }

    qDebug() << "Cycles after pruning:" << cycles.size() << "cycles.";
}

QVector3D ContourGraph::SampleCycleEdgeContour(const int c_ind, const int e_ind,
                                               const float f)
{
    for (int i = 1; i < cycle_edge_arclengths[c_ind][e_ind].size(); ++i) {
        const float arc_1 = cycle_edge_arclengths[c_ind][e_ind][i - 1];
        const float arc_2 = cycle_edge_arclengths[c_ind][e_ind][i];

        if (arc_1 <= f && arc_2 >= f) {
            const QVector3D & p1 = cycle_edge_points[c_ind][e_ind][i - 1];
            const QVector3D & p2 = cycle_edge_points[c_ind][e_ind][i];

            // f is within this segment, we linterp a precise value at f
            const float interp = (f - arc_1) / (arc_2 - arc_1);
            return p1 * (1.0f - interp) + p2 * interp;
        }
    }

    qDebug() << "ContourGraph::SampleCycleEdgeContour - Warning, didn't find "
                "suitable point along cycle"
             << c_ind << "edge" << e_ind << "given parameter" << f;
    return QVector3D(0, 0, 0);
}

void ContourGraph::FitCoonsPatches()
{
    const bool do_hermite = true;

    qDebug() << "Computing coons patches...";
    cycle_coons_patches.resize(cycles.size());

    for (int i = 0; i < cycle_edge_points.size(); ++i) {

        if (cycle_edge_points[i].size() != 4) {
            continue;
        }

        const int patch_samples = 100;

        cycle_coons_patches[i].resize(patch_samples + 1);

        for (int s = 0; s <= patch_samples; ++s) {
            cycle_coons_patches[i][s].resize(patch_samples + 1);

            const float s_f = float(s) / float(patch_samples);

            for (int t = 0; t <= patch_samples; ++t) {
                const float t_f = float(t) / float(patch_samples);

                const QVector3D C_00 = SampleCycleEdgeContour(i, 0, 0.0f);
                const QVector3D C_01 = SampleCycleEdgeContour(i, 0, 1.0f);
                const QVector3D C_10 = SampleCycleEdgeContour(i, 2, 1.0f);
                const QVector3D C_11 = SampleCycleEdgeContour(i, 2, 0.0f);
                const QVector3D C_s0 =
                    SampleCycleEdgeContour(i, 0, s_f);  // 1st
                const QVector3D D_1t =
                    SampleCycleEdgeContour(i, 1, t_f);  // 2nd
                const QVector3D C_s1 =
                    SampleCycleEdgeContour(i, 2, 1.0f - s_f);  // 3rd - flip
                const QVector3D D_0t =
                    SampleCycleEdgeContour(i, 3, 1.0f - t_f);  // 4th - flip

                if (do_hermite) {
                    const float s_f2 = s_f * s_f;
                    const float s_f3 = s_f2 * s_f;
                    const float t_f2 = t_f * t_f;
                    const float t_f3 = t_f2 * t_f;

                    float coeff_s[4];
                    float coeff_t[4];

                    coeff_s[0] = (2.0f * s_f3 - 3.0f * s_f2 + 1.0f);
                    coeff_s[1] = (s_f3 - 2.0f * s_f2 + s_f);
                    coeff_s[2] = (-2.0f * s_f3 + 3.0f * s_f2);
                    coeff_s[3] = (s_f3 - s_f2);

                    coeff_t[0] = (2.0f * t_f3 - 3.0f * t_f2 + 1.0f);
                    coeff_t[1] = (t_f3 - 2.0f * t_f2 + t_f);
                    coeff_t[2] = (-2.0f * t_f3 + 3.0f * t_f2);
                    coeff_t[3] = (t_f3 - t_f2);

                    // this one doesn't incorporate normals (but best one yet)
                    const QVector3D H_c = C_s0 * coeff_t[0] + C_s1 * coeff_t[2];
                    const QVector3D H_d = D_0t * coeff_s[0] + D_1t * coeff_s[2];                    

                    const QVector3D B = C_00 * coeff_s[0] * coeff_t[0] +
                                        C_01 * coeff_s[2] * coeff_t[0] +
                                        C_11 * coeff_s[2] * coeff_t[2] +
                                        C_10 * coeff_s[0] * coeff_t[2];                    

                    cycle_coons_patches[i][s][t] = H_c + H_d - B;

                    // s -> u, t -> v
                    // c0_0 -> P(0, 0)
                    // c0_1 -> P(0, 1)
                    // c1_1 -> P(1, 1)
                    // c1_0 -> P(1, 0)
                    // c0_s -> P(s, 0)
                    // c1_s -> P(s, 1)
                    // d1_t -> P(1, t)
                    // d0_t -> P(0, t)
                    // F1_u -> coeff_s[0]
                    // F2_u -> coeff_s[2]
                    // F3_u -> coeff_s[1]
                    // F4_u -> coeff_s[3]
                    // F1_w -> coeff_t[0]
                    // F2_w -> coeff_t[2]
                    // F3_w -> coeff_t[1]
                    // F4_w -> coeff_t[3]
                    // using ZERO for twist vectors

                } else {
                    const QVector3D L_c = C_s0 * (1.0f - t_f) + C_s1 * t_f;
                    const QVector3D L_d = D_0t * (1.0f - s_f) + D_1t * s_f;
                    const QVector3D B = C_00 * (1.0f - s_f) * (1.0f - t_f) +
                                        C_01 * s_f * (1.0f - t_f) +
                                        C_11 * s_f * t_f +
                                        C_10 * (1.0f - s_f) * t_f;

                    cycle_coons_patches[i][s][t] = L_c + L_d - B;
                }
            }
        }
    }

    qDebug() << "Done.";
}
