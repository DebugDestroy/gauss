#include "algorithms/path/a_star/conditions.hpp"
#include "algorithms/kinematics/incline_angle.hpp"  // Для calculateWheelAngle
#include "algorithms/geometry/bresenham_line.hpp"

#include <cmath>
#include <sstream>
#include <limits>

namespace algorithms::path::a_star {

    Conditions::Conditions(const core::Config& cfg, core::Logger& lg) : config(cfg), logger(lg) {
        logger.trace("[Conditions] Инициализация проверок проходимости пути");
    }
    
bool Conditions::isVehicleRadiusValid(const algorithms::geometry::PointD& pixel, 
                         const std::vector<std::vector<double>>& binaryMap) const
{
    logger.trace("[PathFinder::isVehicleRadiusValid] Проверка радиуса в точке (" +
                 std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");

    const double x = pixel.x;
    const double y = pixel.y;

    const double effectiveRadius = config.vehicleRadius;
    const int radiusPixels = static_cast<int>(std::ceil(effectiveRadius));
    const double radiusSquared = effectiveRadius * effectiveRadius;

    size_t checkedPixels = 0;
    size_t collisionPixels = 0;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx * dx + dy * dy > radiusSquared) continue;

            const int nx = static_cast<int>(x) + dx;
            const int ny = static_cast<int>(y) + dy;
            ++checkedPixels;

            if (nx >= 0 && ny >= 0 &&
                ny < static_cast<int>(binaryMap.size()) &&
                nx < static_cast<int>(binaryMap[0].size()) &&
                std::fabs(binaryMap[ny][nx] - core::WHITE) < core::EPSILON) {

                ++collisionPixels;
                logger.trace("[PathFinder::isVehicleRadiusValid] Столкновение в (" +
                             std::to_string(nx) + "," + std::to_string(ny) + ")");
            }
        }
    }

    bool result = (collisionPixels == 0);
    logger.debug("[PathFinder::isVehicleRadiusValid] Результат: " + 
                 std::string(result ? "проходимо" : "столкновение") +
                 ", проверено пикселей: " + std::to_string(checkedPixels) + 
                 ", столкновений: " + std::to_string(collisionPixels));

    return result;
}

bool Conditions::isEdgeNavigable(const algorithms::geometry::Edge& edge, 
                    const std::unique_ptr<algorithms::gauss::Pole>& p,
                    const std::vector<std::vector<double>>& binaryMap) const {
    const auto line = algorithms::geometry::bresenhamLine(edge.a, edge.b);
    const double r = config.vehicleRadius;
    
    // Вектор пути и ортогональный ему
    const algorithms::geometry::PointD dir = {line.back().x - line.front().x, line.back().y - line.front().y};
    algorithms::geometry::PointD perp = {-dir.y, dir.x};
    double perpLen = std::hypot(perp.x, perp.y);
    
    // Защита от деления на ноль
        if (perpLen < core::EPSILON) {
            logger.error("Нулевая длина перпендикуляра");
            return false;
        }
        
    perp.x = perp.x / perpLen * r;
    perp.y = perp.y / perpLen * r;

    logger.debug("--- Начало проверки ребра ---");
    logger.debug("Параметры тележки: радиус=" + std::to_string(r) + 
                ", max_уклон=" + std::to_string(config.maxUpDownAngle) + 
                "°, max_крен=" + std::to_string(config.maxSideAngle) + "°");
    // Проверка на минимальную длину ребра
        if (line.size() < 2) {
            logger.error("Слишком короткое ребро для проверки");
            return false;
        }
        logger.trace("Длина ребра: " + std::to_string(line.size()) + " точек");
    
    for (size_t i = 0; i < line.size(); ++i) {
        const algorithms::geometry::PointD& center = line[i];
        double h_center; // Высота центра тележки
        // Проверка границ массива
            if (center.y < 0 || center.x < 0 || 
                center.y >= p->field.size() || 
                center.x >= p->field[0].size()) {
                logger.error("Точка за границами поля: (" + 
                           std::to_string(center.x) + "," + 
                           std::to_string(center.y) + ")");
                return false;
            }
            
        h_center = p->field[static_cast<int>(center.y)][static_cast<int>(center.x)];
        
        std::ostringstream pointHeader;
        pointHeader << "Точка [" << i << "/" << line.size()-1 << "] (" 
                   << center.x << "," << center.y << "): "
                   << "высота=" << h_center;
        logger.trace(pointHeader.str());

        // Проверка переднего колеса
        if (i + r < line.size()) {
            const algorithms::geometry::PointD& frontWheel = line[i + static_cast<size_t>(r)];
            double h_front = p->field[static_cast<int>(frontWheel.y)][static_cast<int>(frontWheel.x)];
            double frontAngle = algorithms::kinematics::calculateWheelAngle(center, frontWheel, p);
            
            logger.trace("  Переднее колесо (" + 
                        std::to_string(frontWheel.x) + "," + std::to_string(frontWheel.y) + 
                        "): угол=" + std::to_string(frontAngle) + 
                        "°, высота=" + std::to_string(h_front));
            
            if (std::fabs(frontAngle) > config.maxUpDownAngle) {
                logger.debug("  ! ПРЕВЫШЕНИЕ: передний угол " + std::to_string(frontAngle) + 
                           "° > допустимого " + std::to_string(config.maxUpDownAngle) + "°");
                return false;
            }
        }

        // Боковые колеса
        algorithms::geometry::PointD leftWheel = {std::round(center.x + perp.x), std::round(center.y + perp.y)};
        algorithms::geometry::PointD rightWheel = {std::round(center.x - perp.x), std::round(center.y - perp.y)};

        double leftAngle = 0, rightAngle = 0;
        if (leftWheel.x >= 0 && leftWheel.y >= 0 && 
            leftWheel.x < p->field[0].size() && leftWheel.y < p->field.size()) {
            double h_left = p->field[static_cast<int>(leftWheel.y)][static_cast<int>(leftWheel.x)];
            leftAngle = algorithms::kinematics::calculateWheelAngle(center, leftWheel, p);
            logger.trace("  Левое колесо (" + 
                        std::to_string(leftWheel.x) + "," + std::to_string(leftWheel.y) + 
                        "): угол=" + std::to_string(leftAngle) + 
                        "°, высота=" + std::to_string(h_left));
        }

        if (rightWheel.x >= 0 && rightWheel.y >= 0 && 
            rightWheel.x < p->field[0].size() && rightWheel.y < p->field.size()) {
            double h_right = p->field[static_cast<int>(rightWheel.y)][static_cast<int>(rightWheel.x)];
            rightAngle = algorithms::kinematics::calculateWheelAngle(center, rightWheel, p);
            logger.trace("  Правое колесо (" + 
                        std::to_string(rightWheel.x) + "," + std::to_string(rightWheel.y) + 
                        "): угол=" + std::to_string(rightAngle) + 
                        "°, высота=" + std::to_string(h_right));
        }

        if (std::fabs(leftAngle) > config.maxSideAngle || 
            std::fabs(rightAngle) > config.maxSideAngle) {
            logger.debug("  ! ПРЕВЫШЕНИЕ: боковые углы L=" + std::to_string(leftAngle) + 
                       "° R=" + std::to_string(rightAngle) + 
                       " > допустимого " + std::to_string(config.maxSideAngle) + "°");
            return false;
        }

        // Проверка коллизий
        bool collisionCheck = isVehicleRadiusValid(center, binaryMap);
        logger.trace("  Коллизия: " + std::string(collisionCheck ? "нет" : "есть"));
        
        if (!collisionCheck) {
            logger.debug("  ! СТОЛКНОВЕНИЕ в точке (" + 
                       std::to_string(center.x) + "," + std::to_string(center.y) + ")");
            return false;
        }
    }
    
    logger.debug("--- Ребро проходимо ---");
    return true;
}
}
