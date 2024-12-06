#ifndef TREE_H
#define TREE_H

#include <QDebug>
#include <QList>
#include <QVector>

struct TreeNode
{
    TreeNode() : index(0), parent(-1), visited(0) {}

    int index;

    int parent;
    QList<int> children;

    int visited;
};

class Tree
{
public:
    Tree();

    void CreateFromGraph(const QVector<QVector<bool> > & graph,
                         const int root_node);

    int GetNumNodes();
    bool IsConnected();
    bool HasCycles();

    int GetRoot();
    bool HasParent(const int i);
    int GetParent(const int i);
    void GetChildren(const int i, QList<int> & children);
    void GetBranch(const int i, QList<int> & branch);

    void GetCycles(QList<QList<int> > & cycles);
    void GetDisconnectedNodes(QList<int> & disc_nodes);

private:
    void GetPathToRoot(const int i, QList<int> & path);
    void TraceBackCycle(const int i1, const int i2, QList<int> & cycle);

    bool is_connected;
    bool has_cycles;
    QVector<TreeNode> nodes;
    QList<QList<int> > cycles;
    QList<int> disconnected_nodes;
};

#endif  // TREE_H
