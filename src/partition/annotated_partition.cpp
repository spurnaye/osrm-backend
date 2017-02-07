#include "partition/annotated_partition.hpp"

#include <algorithm>
#include <climits> // for CHAR_BIT
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>

namespace osrm
{
namespace partition
{

namespace
{
// create a comparator for a given level
auto makeCompare(const std::uint32_t level)
{
    const auto masked = [](const BisectionID id, const std::uint32_t level) {
        // 0.01.1 with 1 starting at the level+1_th most significant bit (level = 0 -> 01..1)
        const auto cut_below_level = (1 << (sizeof(BisectionID) * CHAR_BIT - 1 - level)) - 1;
        const auto mask = std::numeric_limits<BisectionID>::max() ^ cut_below_level;
        return id & mask;
    };

    return [level, masked](const AnnotatedPartition::SizedID lhs,
                           const AnnotatedPartition::SizedID rhs) {
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

    PrintLevels(implicit_tree);
    SearchLevels(balance, implicit_tree);
}

void AnnotatedPartition::PrintLevels(const std::vector<SizedID> &implicit_tree)
{
    // print some statistics on the bisection tree
    std::cout << "Bisection Result:";
    std::queue<BisectionID> id_queue;
    id_queue.push(0);

    const auto add_child = [&id_queue, implicit_tree](const BisectionID prefix,
                                                      const std::uint32_t level) {
        const auto child_range = getChildrenRange(implicit_tree, prefix, level);
        std::cout << "\t" << std::bitset<32>(prefix)
                  << " - Size: " << (child_range.second->count - child_range.first->count)
                  << " Subcells: " << std::distance(child_range.first, child_range.second) << "\n";
        if (std::distance(child_range.first, child_range.second) > 1)
            id_queue.push(prefix);
    };

    for (std::uint32_t level = 0; !id_queue.empty(); ++level)
    {
        auto level_size = id_queue.size();
        std::cout << "Level: " << level << "\n";
        while (level_size--)
        {

            const auto prefix = id_queue.front();
            id_queue.pop();
            add_child(leftChild(prefix, level), level);
            add_child(rightChild(prefix, level), level);
        }
        std::cout << std::endl;
    }
}

void AnnotatedPartition::SearchLevels(double balance, const std::vector<SizedID> &implicit_tree)
{
    typedef std::pair<BisectionID, std::uint32_t> pbu;
    std::vector<pbu> current_level;

    // start searching with level 0 at prefix 0
    current_level.push_back({static_cast<BisectionID>(0), 0u});

    const auto print_level = [&current_level, &implicit_tree]() {
        std::cout << "Level\n";
        for (auto element : current_level)
        {
            std::cout << "\t Prefix: " << std::bitset<32>(element.first)
                      << " Depth: " << element.second
                      << " Size: " << getCellSize(implicit_tree, element.first, element.second - 1)
                      << std::endl;
        }
    };

    while (!current_level.empty())
    {

        std::size_t total_size = 0;
        std::size_t count = 0;
        std::queue<pbu> id_queue;
        for (auto element : current_level)
        {
            // don't relax final cells
            if (element.second == 0 || hasChildren(implicit_tree, element.first, element.second-1))
            {
                total_size += getCellSize(
                    implicit_tree, leftChild(element.first, element.second), element.second);
                id_queue.push(pbu(leftChild(element.first, element.second), element.second + 1));

                total_size += getCellSize(
                    implicit_tree, rightChild(element.first, element.second), element.second);
                id_queue.push(pbu(rightChild(element.first, element.second), element.second + 1));
                count += 2;
            }
        }

        std::size_t max_size = balance * (total_size / static_cast<double>(count));

        current_level.clear();

        const auto relax = [&id_queue, implicit_tree, max_size, &current_level](
            const pbu &element) {
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
                id_queue.push(pbu(leftChild(element.first, element.second), element.second + 1));
                id_queue.push(pbu(rightChild(element.first, element.second), element.second + 1));
            }
        };

        while (!id_queue.empty())
        {
            relax(id_queue.front());
            id_queue.pop();
        }

        print_level();
    };
}

} // namespace partition
} // namespace osrm
