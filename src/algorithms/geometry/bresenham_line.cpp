#include "algorithms/geometry/bresenham_line.hpp"
#include <cmath> // для std::abs

namespace algorithms::geometry {
 std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end) {
    
    std::vector<PointD> linePoints;
    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    int dx = std::abs(x1 - x0);
    int dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    size_t pointCount = 0;
    while (true) {
        linePoints.emplace_back(x0, y0);
        pointCount++;
        
        if (x0 == x1 && y0 == y1) {
            break;
        }
        
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
    
    return linePoints;
}
}
