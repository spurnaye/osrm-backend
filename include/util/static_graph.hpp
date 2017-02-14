#ifndef STATIC_GRAPH_HPP
#define STATIC_GRAPH_HPP

#include "util/integer_range.hpp"
#include "util/percent.hpp"
#include "util/shared_memory_vector_wrapper.hpp"
#include "util/static_graph_traits.hpp"
#include "util/typedefs.hpp"

#include <boost/assert.hpp>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace osrm
{
namespace util
{

namespace static_graph_details
{

using NodeIterator = NodeID;
using EdgeIterator = NodeID;

struct NodeArrayEntry
{
    // index of the first edge
    EdgeIterator first_edge;
};

template <typename EdgeDataT> struct EdgeArrayEntry;

template <typename EdgeDataT> struct EdgeArrayEntry
{
    NodeID target;
    EdgeDataT data;
};

template <> struct EdgeArrayEntry<void>
{
    NodeID target;
};

template <typename EdgeDataT> struct SortableEdgeWithData;

template <> struct SortableEdgeWithData<void>
{
    NodeIterator source;
    NodeIterator target;

    SortableEdgeWithData(NodeIterator source, NodeIterator target) : source(source), target(target)
    {
    }

    bool operator<(const SortableEdgeWithData &right) const
    {
        if (source != right.source)
        {
            return source < right.source;
        }
        return target < right.target;
    }
};

template <typename EdgeDataT> struct SortableEdgeWithData : SortableEdgeWithData<void>
{
    using Base = SortableEdgeWithData<void>;

    EdgeDataT data;

    template <typename... Ts>
    SortableEdgeWithData(NodeIterator source, NodeIterator target, Ts &&... data)
        : Base{source, target}, data{std::forward<Ts>(data)...}
    {
    }
};

} // namespace static_graph_details

template <typename EdgeDataT, bool UseSharedMemory = false> class StaticGraph
{
  public:
    using NodeIterator = static_graph_details::NodeIterator;
    using EdgeIterator = static_graph_details::EdgeIterator;
    using EdgeRange = range<EdgeIterator>;
    using NodeArrayEntry = static_graph_details::NodeArrayEntry;
    using EdgeArrayEntry = static_graph_details::EdgeArrayEntry<EdgeDataT>;

    EdgeRange GetAdjacentEdgeRange(const NodeID node) const
    {
        return irange(BeginEdges(node), EndEdges(node));
    }

    template <typename ContainerT> StaticGraph(const int nodes, const ContainerT &graph)
    {
        BOOST_ASSERT(std::is_sorted(const_cast<ContainerT &>(graph).begin(),
                                    const_cast<ContainerT &>(graph).end()));

        number_of_nodes = nodes;
        number_of_edges = static_cast<EdgeIterator>(graph.size());
        node_array.resize(number_of_nodes + 1);
        EdgeIterator edge = 0;
        EdgeIterator position = 0;
        for (const auto node : irange(0u, number_of_nodes + 1))
        {
            EdgeIterator last_edge = edge;
            while (edge < number_of_edges && graph[edge].source == node)
            {
                ++edge;
            }
            node_array[node].first_edge = position; //=edge
            position += edge - last_edge;           // remove
        }
        edge_array.resize(position); //(edge)
        edge = 0;
        for (const auto node : irange(0u, number_of_nodes))
        {
            EdgeIterator e = node_array[node + 1].first_edge;
            for (const auto i : irange(node_array[node].first_edge, e))
            {
                edge_array[i].target = graph[edge].target;
                CopyDataIfAvailable(
                    edge_array[i], graph[edge], traits::HasDataMember<EdgeArrayEntry>{});
                edge++;
            }
        }
    }

    StaticGraph(typename ShM<NodeArrayEntry, UseSharedMemory>::vector &nodes,
                typename ShM<EdgeArrayEntry, UseSharedMemory>::vector &edges)
    {
        number_of_nodes = static_cast<decltype(number_of_nodes)>(nodes.size() - 1);
        number_of_edges = static_cast<decltype(number_of_edges)>(edges.size());

        using std::swap;
        swap(node_array, nodes);
        swap(edge_array, edges);
    }

    unsigned GetNumberOfNodes() const { return number_of_nodes; }

    unsigned GetNumberOfEdges() const { return number_of_edges; }

    unsigned GetOutDegree(const NodeIterator n) const { return EndEdges(n) - BeginEdges(n); }

    inline NodeIterator GetTarget(const EdgeIterator e) const
    {
        return NodeIterator(edge_array[e].target);
    }

    auto &GetEdgeData(const EdgeIterator e) { return edge_array[e].data; }

    const auto &GetEdgeData(const EdgeIterator e) const { return edge_array[e].data; }

    EdgeIterator BeginEdges(const NodeIterator n) const
    {
        return EdgeIterator(node_array.at(n).first_edge);
    }

    EdgeIterator EndEdges(const NodeIterator n) const
    {
        return EdgeIterator(node_array.at(n + 1).first_edge);
    }

    // searches for a specific edge
    EdgeIterator FindEdge(const NodeIterator from, const NodeIterator to) const
    {
        for (const auto i : irange(BeginEdges(from), EndEdges(from)))
        {
            if (to == edge_array[i].target)
            {
                return i;
            }
        }
        return SPECIAL_EDGEID;
    }

    /**
     * Finds the edge with the smallest `.weight` going from `from` to `to`
     * @param from the source node ID
     * @param to the target node ID
     * @param filter a functor that returns a `bool` that determines whether an edge should be
     * tested or not.
     *   Takes `EdgeData` as a parameter.
     * @return the ID of the smallest edge if any were found that satisfied *filter*, or
     * `SPECIAL_EDGEID` if no
     *   matching edge is found.
     */
    template <typename FilterFunction>
    EdgeIterator
    FindSmallestEdge(const NodeIterator from, const NodeIterator to, FilterFunction &&filter) const
    {
        static_assert(traits::HasDataMember<EdgeArrayEntry>::value,
                      "Filtering on .data not possible without .data member attribute");

        EdgeIterator smallest_edge = SPECIAL_EDGEID;
        EdgeWeight smallest_weight = INVALID_EDGE_WEIGHT;
        for (auto edge : GetAdjacentEdgeRange(from))
        {
            const NodeID target = GetTarget(edge);
            const auto &data = GetEdgeData(edge);
            if (target == to && data.weight < smallest_weight &&
                std::forward<FilterFunction>(filter)(data))
            {
                smallest_edge = edge;
                smallest_weight = data.weight;
            }
        }
        return smallest_edge;
    }

    EdgeIterator FindEdgeInEitherDirection(const NodeIterator from, const NodeIterator to) const
    {
        EdgeIterator tmp = FindEdge(from, to);
        return (SPECIAL_NODEID != tmp ? tmp : FindEdge(to, from));
    }

    EdgeIterator
    FindEdgeIndicateIfReverse(const NodeIterator from, const NodeIterator to, bool &result) const
    {
        EdgeIterator current_iterator = FindEdge(from, to);
        if (SPECIAL_NODEID == current_iterator)
        {
            current_iterator = FindEdge(to, from);
            if (SPECIAL_NODEID != current_iterator)
            {
                result = true;
            }
        }
        return current_iterator;
    }

    const NodeArrayEntry &GetNode(const NodeID nid) const { return node_array[nid]; }
    const EdgeArrayEntry &GetEdge(const EdgeID eid) const { return edge_array[eid]; }

  private:
    template <typename OtherEdge>
    void CopyDataIfAvailable(EdgeArrayEntry &into, const OtherEdge &from, std::true_type)
    {
        into.data = from.data;
    }

    template <typename OtherEdge>
    void CopyDataIfAvailable(EdgeArrayEntry &into, const OtherEdge &from, std::false_type)
    {
        // Graph has no .data member, never copy even if `from` has a .data member.
        (void)into;
        (void)from;
    }

    NodeIterator number_of_nodes;
    EdgeIterator number_of_edges;

    typename ShM<NodeArrayEntry, UseSharedMemory>::vector node_array;
    typename ShM<EdgeArrayEntry, UseSharedMemory>::vector edge_array;
};

} // namespace util
} // namespace osrm

#endif // STATIC_GRAPH_HPP
