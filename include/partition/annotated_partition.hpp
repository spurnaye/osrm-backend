#ifndef OSRM_PARTITION_ANNOTATE_HPP_
#define OSRM_PARTITION_ANNOTATE_HPP_

#include "partition/bisection_graph.hpp"
#include "util/typedefs.hpp"
#include <cstdint>

namespace osrm
{
namespace partition
{

// takes the result of a recursive bisection and turns it into an annotated partition for MLD. These
// annotated partitions provide a mapping from every node in the graph to a consecutively
// numbered cell in each level of the multi level partition. Instead of using the bisection directly
// (which can result in a unbalanced tree structure)
// 
//            _____o______
//           /            \
//          o          ____o____
//         / \        /         \
//        a   b       o         _o_
//                   / \       /   \
//                  c   d     o     o
//                           / \   / \
//                          e  f   g  h
// 
// we build a balanced structure that will result in a multi-cut on any level. We transform this
// layout into:
//
//            _____o__________
//           /        |       \
//          o         |        \
//         / \        |         \
//        a   b       o         _o_
//                   / \       /   \
//                  c   d     o     o
//                           / \   / \
//                          e  f   g  h
class AnnotatedPartition
{
  public:
    // Used to generate an implicit tree representation
    struct SizedID
    {
        BisectionID id;
        std::size_t count;

        bool operator<(const SizedID &other) const { return id < other.id; };
    };

    AnnotatedPartition(const BisectionGraph &graph,
                       const double balance,
                       const std::vector<BisectionID> &bisection_ids);

  private:
    // print distribution of level graph as it is
    void PrintLevels(const std::vector<SizedID> &implicit_tree);

    // find levels according to balance in the component tree
    void SearchLevels(double balance, const std::vector<SizedID> &implicit_tree);
};

} // namespace partition
} // namespace osrm

#endif // OSRM_PARTITION_ANNOTATE_HPP_
