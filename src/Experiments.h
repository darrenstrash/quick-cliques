#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <string>
#include <vector>

class Experiments
{
public:
    Experiments(std::string const sDataSetName, double const dTimeout, bool const bOutputLatex, bool const bPrintHeader, std::vector<std::vector<int>> const &adjacencyArray, std::vector<std::vector<char>> const &vAdjacencyMatrix);

    void RunKernelSize() const;
    void RunFastKernelSize() const;
    void KernelizeAndRunReductionSparseMISS() const;
    void KernelizeAndRunComponentWiseMISS() const;
    void RunComponentWiseMISS() const;
    void KernelizeAndRunComponentWiseReductionSparseMISS() const;
    void ComputeCriticalIndependentSet() const;
    void ComputeCriticalIndependentSetKernel() const;
    void ComputeMaximumCriticalIndependentSetKernel() const;
    void ComputeMaximumCriticalIndependentSet() const;
    void RunExactSearch() const;
    void RunStandardSearch() const;
    void RunComponentsMISS() const;
    void RunComponentsStandardSearch() const;
    void RunComponentsForwardSearch() const;
    void RunForwardSearch() const;

private:
    std::string m_sDataSetName;
    bool const m_bOutputLatex;
    bool const m_bPrintHeader;
    double const m_dTimeout;
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<std::vector<char>> const &m_AdjacencyMatrix;
};
#endif // EXPERIMENTS_H
