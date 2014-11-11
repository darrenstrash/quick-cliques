#ifndef VERTEX_SETS_H

#include <vector>
#include <list>
#include <vector>

class VertexSets {
public:
    VertexSets();
    ~VertexSets();

    void Recurse();
    void Backtrack();

    void MoveFromPToX(int const vertexInP);
    void MoveFromPToR(int const vertexInP);
    void MoveFromRToX(int const vertexInP);

    int  ChoosePivot(std::list<int> &pivotNonNeighbors) const;
    bool InP(int const vertex) const;

protected:
    class SetDelineator {
    public:
        SetDelineator(int const beginX, int const beginP, int const beginR) : m_BeginX(beginX), m_BeginP(beginP), m_BeginR(beginR) {}
        int m_BeginX;
        int m_BeginP;
        int m_BeginR;
    };

private: // members
    std::vector<int> m_VertexSets;
    std::vector<int> m_VertexLocation;
};

#endif //VERTEX_SETS_H
