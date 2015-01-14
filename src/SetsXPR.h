#ifndef SETS_XPR_H
#define SETS_XPR_H

#include "Set.h"

#include <iostream>
#include <set>

class SetsXPR
{

public:
    SetsXPR() : X(), P(), R() {}
    ~SetsXPR() {}

    void PrintSummary() const
    {
        std::cout /*<< line*/ << ": X = [";
        std::set<int> setX;
        setX.insert(X.begin(), X.end());
        for (int const x : setX) {
            std::cout << x << " ";
        }

        std::set<int> setP;
        setP.insert(P.begin(), P.end());
        std::cout << "] P = [";
        for (int const p : setP) {
            std::cout << p << " ";
        }

        std::set<int> setR;
        setR.insert(R.begin(), R.end());
        std::cout << "] R = [";
        for (int const r : setR) {
            std::cout << r << " ";
        }
        std::cout << "]" << std::endl;
    }

    bool XIsEmpty() const { return X.Empty(); }
    bool PIsEmpty() const { return P.Empty(); }

    size_t SizeOfP() const { return P.Size(); }
    size_t SizeOfX() const { return X.Size(); }

    void SaveState()
    {
        X.SaveState();
        P.SaveState();
        R.SaveState();
    }

    void RestoreState()
    {
        X.RestoreState();
        P.RestoreState();
        R.RestoreState();
    }

    void InsertIntoP(int const vertex)
    {
        P.Insert(vertex);
    }

    void MoveFromPToR(int const vertex)
    {
        P.MoveTo(vertex, R);
    }

    void MoveFromXToP(int const vertex)
    {
        X.MoveTo(vertex, P);
    }

    void MoveFromPToX(int const vertex)
    {
        P.MoveTo(vertex, X);
    }

    bool InX(int const vertex) const {
        return X.Contains(vertex);
    }

    bool InP(int const vertex) const {
        return P.Contains(vertex);
    }

    void RemoveFromX(int const vertex) {
        X.Remove(vertex);
    }

    void RemoveFromP(int const vertex) {
        P.Remove(vertex);
    }

    Set const &GetP() const { return P; }

protected:
    Set X;
    Set P;
    Set R;
};

#endif //SETS_XPR_H
