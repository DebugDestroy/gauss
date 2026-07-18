#pragma once
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/geometry/math.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstddef>

namespace algorithms::geometry::spatial
{

struct KDPoint
{
    PointD point;
    int value;
};

struct KDNearestResult
{
    int value;
    double distanceSquared;
};

class KDTree
{
public:

    ~KDTree()
    {
        clear(root);
    }
    
    KDTree() = default;
    
    KDTree(const KDTree&) = delete;
    KDTree& operator=(const KDTree&) = delete;

    KDTree(KDTree&& other) noexcept
    {
        root = other.root;
        other.root = nullptr;
    }


    KDTree& operator=(KDTree&& other) noexcept
    {
        if(this != &other)
        {
            clear(root);

            root = other.root;
            other.root = nullptr;
        }

        return *this;
    }
    
    bool empty() const
    {
        return root == nullptr;
    }
    
    void build(std::vector<KDPoint> points)
    {
        clear(root);

        if(points.empty())
        {
            root = nullptr;
            return;
        }

        root = buildRecursive(
            points,
            0,
            static_cast<int>(points.size()) - 1,
            0);
    }

    KDNearestResult nearest(const PointD& point) const
    {
        const Node* bestNode = nullptr;

        double bestDist2 =
            std::numeric_limits<double>::max();

        if (empty())
            throw std::runtime_error("KDTree is empty");
    
        nearestRecursive(
            root,
            point,
            0,
            bestNode,
            bestDist2);


        return {bestNode->point.value, bestDist2};
    }

    void radiusSearch(
        const PointD& point,
        double radius,
        std::vector<int>& result) const
    {
        result.clear();
           
        radiusSearchRecursive(
            root,
            point,
            radius * radius,
            0,
            result);
    }
    
private:

    struct Node
    {
        KDPoint point;

        Node* left = nullptr;
        Node* right = nullptr;
    };

    Node* root = nullptr;

    Node* buildRecursive(
        std::vector<KDPoint>& points,
        int left,
        int right,
        int depth)
    {
        if (left > right)
            return nullptr;

        const int axis = depth % 2;
        const int mid = left + (right - left) / 2;

        auto begin = points.begin() + static_cast<std::ptrdiff_t>(left);
        auto middle = points.begin() + static_cast<std::ptrdiff_t>(mid);
        auto end = points.begin() + static_cast<std::ptrdiff_t>(right + 1);

        std::nth_element(
            begin,
            middle,
            end,
            [axis](const KDPoint& a, const KDPoint& b)
            {
                if (axis == 0)
                    return a.point.x < b.point.x;

                return a.point.y < b.point.y;
            });

        Node* node = new Node;
        node->point = points[static_cast<std::size_t>(mid)];

        node->left =
            buildRecursive(
                points,
                left,
                mid - 1,
                depth + 1);

        node->right =
            buildRecursive(
                points,
                mid + 1,
                right,
                depth + 1);

        return node;
    }
    
    void nearestRecursive(
        const Node* node,
        const PointD& point,
        int depth,
        const Node*& bestNode,
        double& bestDist2) const
    {
        if (node == nullptr)
            return;


        // 1. Проверяем текущую вершину
        const double dist2 =
            distanceSquared(node->point.point, point);


        if (dist2 < bestDist2)
        {
            bestDist2 = dist2;
            bestNode = node;
        }


        // 2. Определяем ось разделения
        const int axis = depth % 2;


        // 3. Определяем ближнюю и дальнюю ветку

        const Node* nearChild;
        const Node* farChild;


        if ((axis == 0 && point.x < node->point.point.x) ||
            (axis == 1 && point.y < node->point.point.y))
        {
            nearChild = node->left;
            farChild = node->right;
        }
        else
        {
            nearChild = node->right;
            farChild = node->left;
        }


        // 4. Сначала идем туда, где может быть точка
        nearestRecursive(
            nearChild,
            point,
            depth + 1,
            bestNode,
            bestDist2);


        // 5. Нужно ли проверять вторую половину?

        double splitDistance;


        if (axis == 0)
            splitDistance =
                point.x - node->point.point.x;
        else
            splitDistance =
                point.y - node->point.point.y;


        if (splitDistance * splitDistance < bestDist2)
        {
            nearestRecursive(
                farChild,
                point,
                depth + 1,
                bestNode,
                bestDist2);
        }
    }

    void radiusSearchRecursive(
        const Node* node,
        const PointD& point,
        double radiusSquared,
        int depth,
        std::vector<int>& result) const
    {
        if (node == nullptr)
            return;


        // Проверяем текущую точку

        const double dist2 =
            distanceSquared(node->point.point, point);


        if (dist2 <= radiusSquared)
        {
            result.push_back(node->point.value);
        }


        const int axis = depth % 2;


        const Node* nearChild;
        const Node* farChild;


        if ((axis == 0 && point.x < node->point.point.x) ||
            (axis == 1 && point.y < node->point.point.y))
        {
            nearChild = node->left;
            farChild = node->right;
        }
        else
        {
            nearChild = node->right;
            farChild = node->left;
        }


        // Сначала проверяем ближнюю сторону

        radiusSearchRecursive(
            nearChild,
            point,
            radiusSquared,
            depth + 1,
            result);


        // Может ли круг пересекать вторую сторону?

        double splitDistance;


        if (axis == 0)
            splitDistance =
                point.x - node->point.point.x;
        else
            splitDistance =
                point.y - node->point.point.y;


        if (splitDistance * splitDistance <= radiusSquared)
        {
            radiusSearchRecursive(
                farChild,
                point,
                radiusSquared,
                depth + 1,
                result);
        }
    }
    
    void clear(Node*& node)
    {
        if(node == nullptr)
            return;

        clear(node->left);
        clear(node->right);

        delete node;
        node = nullptr;
    }
};
} // namespace algorithms::geometry::spatial
