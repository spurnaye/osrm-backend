#include "partition/annotated_partition.hpp"

#include <algorithm>
#include <climits> // for CHAR_BIT
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <unordered_map>

namespace osrm
{
namespace partition
{

namespace
{

const constexpr std::uint32_t INVALID_CELLID = std::numeric_limits<std::uint32_t>::max();

auto masked(const BisectionID id, const std::uint32_t level)
{
    // 0.01.1 with 1 starting at the level+1_th most significant bit (level = 0 -> 01..1)
    const auto cut_below_level = (1 << (sizeof(BisectionID) * CHAR_BIT - 1 - level)) - 1;
    const auto mask = std::numeric_limits<BisectionID>::max() ^ cut_below_level;
    return id & mask;
}

// create a comparator for a given level
auto makeCompare(const std::uint32_t level)
{
    return [level](const AnnotatedPartition::SizedID lhs, const AnnotatedPartition::SizedID rhs) {
        return masked(lhs.id, level) < masked(rhs.id, level);
    };
}

// build a tree of cells from the IDs present:
auto leftChild(const BisectionID id_prefix, const std::uint32_t /*level*/) { return id_prefix; }

auto rightChild(const BisectionID id_prefix, const std::uint32_t level)
{
    return id_prefix | (1 << (sizeof(BisectionID) * CHAR_BIT - 1 - level));
}

// get the range of all children
auto getChildrenRange(const std::vector<AnnotatedPartition::SizedID> &implicit_tree,
                      const BisectionID id_prefix,
                      const std::uint32_t level)
{
    AnnotatedPartition::SizedID id = {id_prefix, 0};

    // find all elements of the same prefix as id_prefi
    auto range =
        std::equal_range(implicit_tree.begin(), implicit_tree.end(), id, makeCompare(level));

    // don't ever return our sentinel element as included
    if (range.second == implicit_tree.end())
        --range.second;

    return range;
};

auto getCellSize(const std::vector<AnnotatedPartition::SizedID> &implicit_tree,
                 const BisectionID id_prefix,
                 const std::uint32_t level)
{
    auto range = getChildrenRange(implicit_tree, id_prefix, level);
    return range.second->count - range.first->count;
}

bool hasChildren(const std::vector<AnnotatedPartition::SizedID> &implicit_tree,
                 const BisectionID id_prefix,
                 const std::uint32_t level)
{
    auto range = getChildrenRange(implicit_tree, id_prefix, level);
    return std::distance(range.first, range.second) > 1;
}

} // namespace

AnnotatedPartition::AnnotatedPartition(const BisectionGraph &graph,
                                       const double balance,
                                       const std::vector<BisectionID> &bisection_ids)
{
    // create a sorted vector of bisection ids that exist in the network
    std::vector<SizedID> implicit_tree = [&]() {
        std::map<BisectionID, SizedID> existing_ids;

        // insert an ID into the sized_id set or increase the count if the element should be already
        // present in the set of known ids
        const auto insert_or_augment = [&existing_ids](const BisectionID id) {
            SizedID sized_id = {id, 1};
            auto maybe_existing_id = existing_ids.find(id);
            if (maybe_existing_id == existing_ids.end())
                existing_ids[id] = sized_id;
            else
                maybe_existing_id->second.count++;
        };
        std::for_each(bisection_ids.begin(), bisection_ids.end(), insert_or_augment);

        std::vector<SizedID> result;
        result.resize(existing_ids.size() + 1);
        std::transform(existing_ids.begin(),
                       existing_ids.end(),
                       result.begin(),
                       [](const auto &pair) { return pair.second; });

        // sentinel
        result.back() = {std::numeric_limits<BisectionID>::max(), 0};

        return result;
    }();

    // calculate a prefix sum over all sorted IDs, this allows to get the size of any partition in
    // the array/level based on the prefix and lower bound on prefixes.
    // e.g 00,01,10,11 allow to search for (0) (1) to find (00) and (10) as lower bounds. The
    // difference in count is the size of all cells in the left part of the partition.
    std::transform(implicit_tree.begin(),
                   implicit_tree.end(),
                   implicit_tree.begin(),
                   [sum = std::size_t{0}](SizedID id) mutable {
                       const auto new_sum = sum + id.count;
                       id.count = sum;
                       sum = new_sum;
                       return id;
                   });

    PrintBisection(implicit_tree, graph, bisection_ids);
    SearchLevels(balance, implicit_tree, graph, bisection_ids);
}

void AnnotatedPartition::PrintBisection(const std::vector<SizedID> &implicit_tree,
                                        const BisectionGraph &graph,
                                        const std::vector<BisectionID> &bisection_ids) const
{
    // print some statistics on the bisection tree
    std::cout << "[Unmodified Bisection]:\n";

    std::queue<BisectionID> id_queue;
    id_queue.push(0);

    const auto add_child = [&id_queue, implicit_tree](const BisectionID prefix,
                                                      const std::uint32_t level) {
        const auto child_range = getChildrenRange(implicit_tree, prefix, level);
        if (std::distance(child_range.first, child_range.second) > 1)
            id_queue.push(prefix);
    };

    for (std::uint32_t level = 0; !id_queue.empty(); ++level)
    {
        auto level_size = id_queue.size();
        std::vector<std::pair<BisectionID, std::uint32_t>> current_level;
        while (level_size--)
        {
            const auto prefix = id_queue.front();
            id_queue.pop();
            if (level == 0 || hasChildren(implicit_tree, prefix, level-1))
            {
                current_level.push_back(
                    std::pair<BisectionID, std::uint32_t>(leftChild(prefix, level), level+1));
                current_level.push_back(
                    std::pair<BisectionID, std::uint32_t>(rightChild(prefix, level), level+1));
            }
            add_child(leftChild(prefix, level), level);
            add_child(rightChild(prefix, level), level);
        }
        if (!current_level.empty())
        {
            const auto cell_ids = ComputeCellIDs(current_level, graph, bisection_ids);
            const auto stats = AnalyseLevel(graph, cell_ids);
            stats.print(std::cout);
        }
    }
}

void AnnotatedPartition::SearchLevels(double balance,
                                      const std::vector<SizedID> &implicit_tree,
                                      const BisectionGraph &graph,
                                      const std::vector<BisectionID> &bisection_ids) const
{
    std::vector<std::pair<BisectionID, std::uint32_t>> current_level;
    std::cout << "[balanced via DFS]\n";

    // start searching with level 0 at prefix 0
    current_level.push_back({static_cast<BisectionID>(0), 0u});

    const auto print_level = [&]() {
        if (current_level.empty())
            return;
        const auto cell_ids = ComputeCellIDs(current_level, graph, bisection_ids);
        const auto stats = AnalyseLevel(graph, cell_ids);
        stats.print(std::cout);
    };

    while (!current_level.empty())
    {
        std::size_t total_size = 0;
        std::size_t count = 0;
        std::queue<std::pair<BisectionID, std::uint32_t>> id_queue;
        for (auto element : current_level)
        {
            // don't relax final cells
            if (element.second == 0 ||
                hasChildren(implicit_tree, element.first, element.second - 1))
            {
                total_size += getCellSize(
                    implicit_tree, leftChild(element.first, element.second), element.second);
                id_queue.push(std::pair<BisectionID, std::uint32_t>(
                    leftChild(element.first, element.second), element.second + 1));

                total_size += getCellSize(
                    implicit_tree, rightChild(element.first, element.second), element.second);
                id_queue.push(std::pair<BisectionID, std::uint32_t>(
                    rightChild(element.first, element.second), element.second + 1));
                count += 2;
            }
        }

        std::size_t max_size = balance * (total_size / static_cast<double>(count));

        current_level.clear();

        const auto relax = [&id_queue, implicit_tree, max_size, &current_level](
            const std::pair<BisectionID, std::uint32_t> &element) {
            const auto size = getCellSize(implicit_tree, element.first, element.second - 1);
            if (size <= max_size)
            {
                current_level.push_back(element);
            }
            else if (!hasChildren(implicit_tree, element.first, element.second))
            {
                current_level.push_back(element);
            }
            else
            {
                id_queue.push(std::pair<BisectionID, std::uint32_t>(
                    leftChild(element.first, element.second), element.second + 1));
                id_queue.push(std::pair<BisectionID, std::uint32_t>(
                    rightChild(element.first, element.second), element.second + 1));
            }
        };

        while (!id_queue.empty())
        {
            relax(id_queue.front());
            id_queue.pop();
        }
        print_level();
    }
}

AnnotatedPartition::LevelMetrics
AnnotatedPartition::AnalyseLevel(const BisectionGraph &graph,
                                 const std::vector<std::uint32_t> &cell_ids) const
{
    std::unordered_map<std::uint32_t, std::size_t> cell_sizes;
    std::unordered_map<std::uint32_t, std::size_t> border_nodes;
    std::unordered_map<std::uint32_t, std::size_t> border_arcs;

    // compute basic metrics of the level
    std::size_t border_nodes_total = 0;
    std::size_t border_arcs_total = 0;
    std::size_t contained_nodes = 0;

    for (const auto &node : graph.Nodes())
    {
        const auto cell_id = cell_ids[node.original_id];
        if (cell_id == INVALID_CELLID)
            continue;

        ++contained_nodes;
        const auto edge_range = graph.Edges(node);
        const auto border_arcs_at_node = std::count_if(
            edge_range.begin(), edge_range.end(), [&cell_id, &cell_ids, &graph](const auto &edge) {
                const auto target_cell_id = cell_ids[graph.Node(edge.target).original_id];
                return target_cell_id != cell_id;
            });

        cell_sizes[cell_id]++;
        border_arcs[cell_id] += border_arcs_at_node;
        border_arcs_total += border_arcs_at_node;
        if (border_arcs_at_node)
        {
            border_nodes[cell_id]++;
            ++border_nodes_total;
        }
    }

    const auto by_size = [](const std::pair<std::uint32_t, std::size_t> &lhs,
                            const std::pair<std::uint32_t, std::size_t> &rhs) {
        return lhs.second < rhs.second;
    };
    const auto max_nodes =
        border_nodes.empty()
            ? 0
            : std::max_element(border_nodes.begin(), border_nodes.end(), by_size)->second;
    const auto max_arcs =
        border_arcs.empty()
            ? 0
            : std::max_element(border_arcs.begin(), border_arcs.end(), by_size)->second;

    const auto squarded_size = [](const std::size_t accumulated,
                                  const std::pair<std::uint32_t, std::size_t> &element) {
        return accumulated + element.second * element.second;
    };

    const auto memory =
        4 * std::accumulate(border_arcs.begin(), border_nodes.end(), std::size_t(0), squarded_size);

    return {border_nodes_total,
            border_arcs_total,
            contained_nodes,
            border_nodes.size(),
            max_nodes,
            max_arcs,
            memory};
}

std::vector<std::uint32_t> AnnotatedPartition::ComputeCellIDs(
    const std::vector<std::pair<BisectionID, std::uint32_t>> &prefixes,
    const BisectionGraph &graph,
    const std::vector<BisectionID> &bisection_ids) const
{
    std::vector<std::uint32_t> cell_ids(graph.NumberOfNodes(), INVALID_CELLID);

    for (const auto &node : graph.Nodes())
    {
        // find the cell_id of node in the current levels
        const auto id = bisection_ids[node.original_id];

        const auto is_prefixed_by = [id](const auto &prefix) {
            return masked(id, prefix.second - 1) == prefix.first;
        };

        const auto prefix = std::find_if(prefixes.begin(), prefixes.end(), is_prefixed_by);

        if (prefix != prefixes.end())
            cell_ids[node.original_id] = std::distance(prefixes.begin(), prefix);
    }

    return cell_ids;
}

} // namespace partition
} // namespace osrm
