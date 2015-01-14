#ifndef ARRAY_SETS_XPR_H
#define ARRAY_SETS_XPR_H

#include <iostream>
#include <vector>

class IteratorShim
{
public:
    IteratorShim(std::vector<int>::const_iterator _begin, std::vector<int>::const_iterator _end) : m_Begin(_begin), m_End(_end)
    {
    }
    ~IteratorShim() {}

    std::vector<int>::const_iterator begin() const { return m_Begin; }
    std::vector<int>::const_iterator end()   const { return m_End;   }

    std::vector<int>::const_iterator m_Begin;
    std::vector<int>::const_iterator m_End;
};


class ArraySetsXPR
{
protected:
    class SetDelineator {
    public:
        SetDelineator(int const beginX, int const beginP, int const beginR) : m_BeginX(beginX), m_BeginP(beginP), m_BeginR(beginR) {}
        int m_BeginX;
        int m_BeginP;
        int m_BeginR;
    };

public:
    ArraySetsXPR(size_t const size) : vertexSets(size,-1), vertexLookup(size,-1), m_lDelineators(), beginX(0), beginP(0), beginR(0)
    {
    }
    ~ArraySetsXPR() {}

    void PrintSummary() const
    {
        std::cout /*<< line*/ << ": X = [";
        std::set<int> setX;
        setX.insert(vertexSets.begin() + beginX, vertexSets.begin() + beginP);
        for (int const x : setX) {
            std::cout << x << " ";
        }

        std::set<int> setP;
        setP.insert(vertexSets.begin() + beginP, vertexSets.begin() + beginR);
        std::cout << "] P = [";
        for (int const p : setP) {
            std::cout << p << " ";
        }

////        std::set<int> setR;
////        setR.insert(R.begin(), R.end());
////        setX.insert(vertexSets.begin() + beginX, vertexSets.begin() + beginP);
////        std::cout << "] R = [";
////        for (int const r : setR) {
////            std::cout << r << " ";
////        }
        std::cout << "]" << std::endl;
    }

    bool XIsEmpty() const { return beginP <= beginX; }
    bool PIsEmpty() const { return beginR <= beginP; }

    size_t SizeOfP() const { return beginR - beginP; }
    size_t SizeOfX() const { return beginP - beginX; }

    void SaveState()
    {
        m_lDelineators.emplace_back(SetDelineator(beginX, beginP, beginR));
    }

    void RestoreState()
    {
        SetDelineator const &delineator = m_lDelineators.back();

        ////    std::cout << __LINE__ << ": " << CheckP() << std::endl;

        beginX = delineator.m_BeginX;
        beginP = delineator.m_BeginP;
        beginR = delineator.m_BeginR;
        m_lDelineators.pop_back();
    }

    // assumes no vertices are in R.
    void InsertIntoP(int const vertex)
    {
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;

        beginR++;
    }

    void MoveFromRToP(int const vertex)
    {
        int const vertexLocation = vertexLookup[vertex];
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;
        beginR++;
    }

    void MoveFromPToR(int const vertex)
    {
        if (vertex < 0 || vertex >= vertexLookup.size() || !InP(vertex)) {
            std::cout << "Attempting to move invalid vertex : " << vertex << std::endl << std::flush;
        }
        int const vertexLocation = vertexLookup[vertex];
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;
    }

    void MoveFromXToP(int const vertex)
    {
        int const vertexLocation = vertexLookup[vertex];
        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
    }

    void MoveFromPToX(int const vertex)
    {
        int const vertexLocation = vertexLookup[vertex];

        //swap vertex into X and increment beginP and beginR
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexLookup[vertexSets[beginP]] = vertexLocation;
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        beginP = beginP + 1;
    }

    bool InX(int const vertex) const {
        int const vertexLocation(vertexLookup[vertex]);
        return vertexLocation >= beginX && vertexLocation < beginP;
    }

    bool InP(int const vertex) const {
        int const vertexLocation(vertexLookup[vertex]);
        return vertexLocation >= beginP && vertexLocation < beginR;
    }

    void RemoveFromX(int const vertex)
    {
        int const vertexLocation = vertexLookup[vertex];
        vertexSets[vertexLocation] = vertexSets[beginX];
        vertexLookup[vertexSets[beginX]] = vertexLocation;
        vertexSets[beginX] = vertex;
        vertexLookup[vertex] = beginX;
        beginX++;
    }

    void RemoveFromP(int const vertex)
    {
        int const vertexLocation = vertexLookup[vertex];
        beginR--;
        vertexSets[vertexLocation] = vertexSets[beginR];
        vertexLookup[vertexSets[beginR]] = vertexLocation;
        vertexSets[beginR] = vertex;
        vertexLookup[vertex] = beginR;
    }

    IteratorShim GetP() const { return IteratorShim(vertexSets.begin() + beginP, vertexSets.begin() + beginR); }

protected:
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<SetDelineator> m_lDelineators;

    int beginX;
    int beginP;
    int beginR;
};

#endif //ARRAY_SETS_XPR_H
