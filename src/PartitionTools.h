#ifndef PARTITION_TOOLS_H
#define PARTITION_TOOLS_H

#include <vector>
#include <string>

namespace PartitionTools
{
    void ReadPartitionFile(std::string const &filename, std::vector<int> &vVertexToPartition, std::vector<std::vector<int>> &vPartitions);
};
#endif //PARTITION_TOOLS_H
