#ifndef ARRAY_SET_H
#define ARRAY_SET_H

#include <set>
#include <vector>
#include <iostream>

class ArraySet
{
public:
    ArraySet(size_t const size) : m_Lookup(size, -1), m_Elements(size, -1), m_iBegin(0), m_iEnd(-1)
    {
    }

    ~ArraySet() {}

    void PrintSummary() const
    {
        for (int const p : *this) {
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }

    bool Contains(int const x) const {
        return m_Lookup[x] != -1;
    }

    void Insert(int const x) {
        if (Contains(x)) return;
        m_iEnd++;
        m_Lookup[x] = m_iEnd;
        m_Elements[m_iEnd] = x;
    }

    void Remove(int const x) {
        if (!Contains(x)) return;
        int const locationX(m_Lookup[x]);
        m_Elements[locationX] = m_Elements[m_iEnd];
        m_Lookup[m_Elements[locationX]] = locationX;
        m_Lookup[x] = -1;
        m_iEnd--;
    }

    void MoveTo(int const x, ArraySet &other) {
        Remove(x);
        other.Insert(x);
    }

    size_t Size()  const { return m_iEnd - m_iBegin + 1; }
    bool   Empty() const { return m_iEnd < m_iBegin; }

    std::vector<int>::iterator begin() { return m_Elements.begin() + m_iBegin; }
    std::vector<int>::iterator end()   { return m_Elements.begin() + m_iEnd;   }

    std::vector<int>::const_iterator begin() const { return m_Elements.begin() + m_iBegin; }
    std::vector<int>::const_iterator end()   const { return m_Elements.begin() + m_iEnd;   }

    static bool Test() {
        std::cout << "ArraySet: ";
        ArraySet testSet(3);
        testSet.Insert(0);
        testSet.Insert(1);
        testSet.Insert(2);

        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: ArraySet failed ContainsTest" << std::endl;
            return false;
        }

        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Size() != 2) {
            std::cout << "FAILED: ArraySet failed RemoveFirst Test" << std::endl;
            return false;
        }

        testSet.Insert(0);
        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: ArraySet failed RemoveAdd Test" << std::endl;
            return false;
        }

        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Size() != 2) {
            std::cout << "FAILED: ArraySet failed RemoveLast Test" << std::endl;
            return false;
        }

        testSet.Insert(0);
        testSet.Remove(1);
        if (testSet.Contains(1) || testSet.Size() != 2) {
            std::cout << "FAILED: ArraySet failed RemoveMiddle Test" << std::endl;
            return false;
        }

        testSet.Insert(1);
        testSet.Remove(1);
        testSet.Remove(2);
        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Contains(1) || testSet.Contains(2) || testSet.Size() != 0 || !testSet.Empty()) {
            std::cout << "FAILED: ArraySet failed RemoveAll Test" << std::endl;
            return false;
        }

        std::cout << "PASSED!" << std::endl << std::flush;
        return true;
    }

private:
    std::vector<int> m_Lookup;
    std::vector<int> m_Elements;
    int m_iBegin;
    int m_iEnd;
};

#endif // ARRAY_SET_H
