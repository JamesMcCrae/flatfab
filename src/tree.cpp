#include "tree.h"

Tree::Tree()
{
    is_connected = false;
    has_cycles = false;
}

void Tree::CreateFromGraph(const QVector<QVector<bool> > & graph,
                           const int root_node)
{
    is_connected = true;
    has_cycles = false;
    cycles.clear();
    nodes.clear();
    disconnected_nodes.clear();

    if (graph.size() < 1) {
        return;
    }

    nodes = QVector<TreeNode>(graph.size());

    for (int i = 0; i < nodes.size(); ++i) {
        nodes[i].index = i;
        nodes[i].visited = 0;  // white
    }

    // maintain a list of nodes that await being visited, use index 0 as root
    QList<int> visiting;
    nodes[0].visited = 1;  // grey
    visiting.push_back(root_node);

    int num_visited = 0;

    // the following performs BFS from the root
    // qDebug() << "doing BFS";
    while (!visiting.empty()) {
        const int cur_node = visiting.first();
        // qDebug() << "visiting" << cur_node;
        visiting.pop_front();

        nodes[cur_node].visited = 2;  // black

        for (int i = 0; i < graph.size(); ++i) {
            if (cur_node == i || !graph[cur_node][i]) {
                continue;
            }

            if (nodes[i].visited == 1) {  // node is grey (we have discovered
                                          // it)

                has_cycles = true;

                // add the cycle to this list
                // qDebug() << "found cycle via backedge between" << cur_node <<
                // i;
                QList<int> new_cycle;
                TraceBackCycle(cur_node, i, new_cycle);
                cycles.push_back(new_cycle);

            } else if (nodes[i].visited ==
                       0) {  // node is white (we have not discovered it)

                // set the topology now
                nodes[cur_node].children.push_back(i);

                nodes[i].visited = 1;  // set node to grey (discovered)
                nodes[i].parent = cur_node;

                visiting.push_back(i);
            }
        }

        ++num_visited;
    }

    if (num_visited != nodes.size()) {
        is_connected = false;

        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].visited != 2) {
                disconnected_nodes.push_back(i);
            }
        }
    }

    // qDebug() << "has cycles?" << has_cycles << "is connected?" <<
    // is_connected;
}

int Tree::GetNumNodes() { return nodes.size(); }

bool Tree::IsConnected() { return is_connected; }

bool Tree::HasCycles() { return has_cycles; }

int Tree::GetRoot() { return 0; }

bool Tree::HasParent(const int i) { return (nodes[i].parent != -1); }

int Tree::GetParent(const int i) { return nodes[i].parent; }

void Tree::GetChildren(const int i, QList<int> & children)
{
    children = nodes[i].children;
}

void Tree::GetBranch(const int i, QList<int> & branch)
{
    QSet<int> visited;
    QList<int> to_visit;
    to_visit.push_back(i);

    while (!to_visit.empty()) {
        const int j = to_visit.first();
        to_visit.pop_front();

        for (int k = 0; k < nodes[j].children.size(); ++k) {
            const int c = nodes[j].children[k];
            if (!visited.contains(c)) {
                to_visit.push_back(c);
            }
        }

        visited.insert(j);
    }

    branch = visited.values();
}

void Tree::GetCycles(QList<QList<int> > & cyc) { cyc = cycles; }

void Tree::GetDisconnectedNodes(QList<int> & disc_nodes)
{
    disc_nodes = disconnected_nodes;
}

void Tree::GetPathToRoot(const int i, QList<int> & path)
{
    int cur_node = i;
    while (cur_node != -1) {
        path.push_back(cur_node);
        cur_node = nodes[cur_node].parent;
    }
}

void Tree::TraceBackCycle(const int i1, const int i2, QList<int> & cycle)
{
    // 1.  get the parent lists of i1 and i2
    QList<int> path_i1;
    QList<int> path_i2;

    GetPathToRoot(i1, path_i1);
    GetPathToRoot(i2, path_i2);

    QSet<int> s1(path_i1.begin(), path_i1.end());
    QSet<int> s2(path_i2.begin(), path_i2.end());

    // 2.  for a tree, any i1 and i2 in the tree must share one common
    // (grand)parent
    QSet<int> common_parents = s1.intersect(s2);
    if (common_parents.empty()) {
        qDebug() << "Tree::TraceBackCycle - Problem!  Could not find common "
                    "parent for"
                 << i1 << i2;
        qDebug() << s1;
        qDebug() << s2;
        return;
    }

    // 3.  get the parent with smallest index in path_i1 (since multiple parents
    // can be shared)
    int common_parent = -1;
    for (int i = 0; i < path_i1.size(); ++i) {
        if (common_parents.contains(path_i1[i])) {
            common_parent = path_i1[i];
            break;
        }
    }

    if (common_parent == -1) {
        qDebug() << "Tree::TraceBackCycle - Problem!  Could not find common "
                    "parent for"
                 << i1 << i2 << "within sets";
        return;
    }

    // 4.  form a cycle which goes across the common parent
    //   4a.  loop up i1's path list
    const int ind1 = path_i1.indexOf(common_parent);
    const int ind2 = path_i2.indexOf(common_parent);

    for (int i = 0; i <= ind1; ++i) {
        cycle.push_back(path_i1[i]);
    }

    //  4b. loop down i2's path list
    for (int i = ind2 - 1; i >= 0; --i) {
        cycle.push_back(path_i2[i]);
    }
}
