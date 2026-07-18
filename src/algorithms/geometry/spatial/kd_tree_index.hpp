#pragma once

#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/geometry/spatial/kd_tree.hpp"

#include <vector>
#include <limits>
#include <stdexcept>
#include <concepts>

namespace algorithms::geometry::spatial
{

template<typename T>
concept KDTreeNode = requires(T node)
{
    { node.point } -> std::convertible_to<PointD>;
};

template<KDTreeNode Node>
class KDTreeIndex
{
private:

    KDTree tree;

    std::vector<KDPoint> buffer;

    std::size_t rebuildSize;


public:
    
    explicit KDTreeIndex(std::size_t rebuildSize_) : rebuildSize(rebuildSize_) {}
    KDTreeIndex() = delete;
    
    void insert(const std::vector<Node>& points,
                int value)
    {
        buffer.push_back(
            KDPoint{
                points[static_cast<std::size_t>(value)].point,
                value
            }
        );

        if(buffer.size() >= rebuildSize)
        {
            rebuild(points);
        }
    }


    int nearest(const PointD& point) const
    {
        bool found = false;

        int bestValue{};
        double bestDistance =
            std::numeric_limits<double>::max();


        if(!tree.empty())
        {
            auto result = tree.nearest(point);

            bestValue = result.value;
            bestDistance = result.distanceSquared;

            found = true;
        }


        for(const auto& item : buffer)
        {
            double dist =
                distanceSquared(
                    item.point,
                    point);

            if(dist < bestDistance)
            {
                bestDistance = dist;
                bestValue = item.value;
                found = true;
            }
        }


        if(!found)
            throw std::runtime_error(
                      "KDTreeIndex is empty");


        return bestValue;
    }


    void radiusSearch(const PointD& point,
                      double radius,
                      std::vector<int>& result) const
    {
        result.clear();


        // Основное дерево

        tree.radiusSearch(
            point,
            radius,
            result);



        // Буфер

        const double radiusSquared =
            radius * radius;


        for(const auto& item : buffer)
        {
            if(distanceSquared(item.point, point)
                <= radiusSquared)
            {
                result.push_back(item.value);
            }
        }
    }


private:


    void rebuild(const std::vector<Node>& points)
    {
        std::vector<KDPoint> temp;

        temp.reserve(
            points.size() + buffer.size()
        );


        for(std::size_t i = 0; i < points.size(); i++)
        {
            temp.push_back(
                {
                    points[i].point,
                    static_cast<int>(i)
                }
            );
        }


        temp.insert(
            temp.end(),
            buffer.begin(),
            buffer.end()
        );


        tree.build(std::move(temp));

        buffer.clear();
    }

};


} // namespace algorithms::geometry::spatial
