#pragma once
#include <vector>       // Для std::vector
#include <string>       // Для std::string
#include <cstdio>       // Для FILE, fprintf, popen, pclose
#include <utility>      // Для std::pair
#include <cmath>        // Для std::sqrt, std::abs
#include <memory>       // Для std::unique_ptr

class GnuplotInterface {
private:
    Logger& logger;

    double transformY(double y, int height) const {
        double transformed = height - y - 1;
        logger.trace(std::string("[GnuplotInterface::transformY] Преобразование Y: ") + 
                   std::to_string(y) + " -> " + std::to_string(transformed));
        return transformed;
    }

    void logPlotStart(const std::string& plotType, const std::string& filename) const {
        logger.info(std::string("[GnuplotInterface] Начало визуализации: ") + plotType);
        logger.debug(std::string("Файл вывода: ") + filename);
    }

    void logPlotEnd(const std::string& plotType) const {
        logger.info(std::string("[GnuplotInterface] Визуализация завершена: ") + plotType);
    }

public:
    GnuplotInterface(Logger& lg) : logger(lg) {
        logger.trace("[GnuplotInterface] Инициализация интерфейса Gnuplot");
    }

    void plotBinaryWithComponents(const std::vector<std::vector<double>>& CopyPole, 
                                const std::vector<Component>& components, 
                                const std::string& filename) {
        logPlotStart("BinaryWithComponents", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotBinaryWithComponents] Ошибка открытия gnuplot pipe");
            return;
        }

        const int height = CopyPole.size();
        const int width = CopyPole[0].size();
        logger.debug("Размер данных: " + std::to_string(width) + "x" + std::to_string(height) + 
                   ", компонентов: " + std::to_string(components.size()));

        // Настройки графика
        fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
        fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
        fprintf(gnuplotPipe, "set title 'Binary Image with Components Metadata'\n");
        fprintf(gnuplotPipe, "set size ratio -1\n");
        fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
        fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
        fprintf(gnuplotPipe, "unset key\n");

        // Многослойный график
        fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
        fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'red', \\\n");
        fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
        fprintf(gnuplotPipe, "'-' with vectors head filled lc 'green'\n");

        // 1. Данные бинарного изображения (с инверсией Y)
for (int y = height - 1; y >= 0; --y) {
    for (int x = 0; x < width; ++x) {
        fprintf(gnuplotPipe, "%f ", CopyPole[y][x]);
    }
            fprintf(gnuplotPipe, "\n");
        }
        fprintf(gnuplotPipe, "e\n");

         // 2. Границы компонент (красные прямоугольники)
        for (const auto& comp : components) {
            double y_min = transformY(comp.min_y, height);
            double y_max = transformY(comp.max_y, height);
            
            fprintf(gnuplotPipe, "%f %f\n", comp.min_x - 0.5, y_min - 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.max_x + 0.5, y_min - 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.max_x + 0.5, y_max + 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.min_x - 0.5, y_max + 0.5);
            fprintf(gnuplotPipe, "%f %f\n\n", comp.min_x - 0.5, y_min - 0.5);
        }
        fprintf(gnuplotPipe, "e\n");

        // 3. Центры компонент (синие точки)
        for (const auto& comp : components) {
            fprintf(gnuplotPipe, "%f %f\n", 
                    comp.center_x, 
                    transformY(comp.center_y, height));
        }
        fprintf(gnuplotPipe, "e\n");

         // 4. Собственные векторы (зеленые стрелки)
for (const auto& comp : components) {
    double cy = transformY(comp.center_y, height);
    
    // Масштабируем собственные вектора по собственным значениям
    double scale_factor = 1; // Общий масштаб для визуализации
    double vec1_scale = scale_factor * sqrt(comp.eigenvalue1);
    double vec2_scale = scale_factor * sqrt(comp.eigenvalue2);
    
    // Первый собственный вектор
    fprintf(gnuplotPipe, "%f %f %f %f\n", 
            comp.center_x,
            cy,
            comp.eigenvec1_x * vec1_scale,
            -comp.eigenvec1_y * vec1_scale);

    // Второй собственный вектор
    fprintf(gnuplotPipe, "%f %f %f %f\n", 
            comp.center_x,
            cy,
            comp.eigenvec2_x * vec2_scale,
            -comp.eigenvec2_y * vec2_scale);
}
fprintf(gnuplotPipe, "e\n");

        pclose(gnuplotPipe);
        logPlotEnd("BinaryWithComponents");
    }
      void gnuplot(std::unique_ptr<Pole>& p, const std::string& filename) {
        logPlotStart("HeightMap3D", filename);
        
        if (p == nullptr) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка: данные высот не инициализированы");
            return;
        }

        int rows = p->field.size();
        int cols = p->field[0].size();
        logger.debug(std::string("Размер сетки: ") + std::to_string(cols) + "x" + std::to_string(rows));
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // Настройки 3D графика
    fprintf(gnuplotPipe, "set terminal pngcairo enhanced size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set zlabel 'Height'\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1);
    fprintf(gnuplotPipe, "set zrange [*:*]\n"); // Автомасштабирование по Z
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n"); // Угол обзора
    
    // Формат данных: x y z
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Height Map'\n");
    
    // Записываем данные в правильном порядке
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            // Преобразуем координаты для правильной ориентации
            fprintf(gnuplotPipe, "%d %d %f\n", x, rows-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n"); // Пустая строка между слоями Y
    }
    
    fprintf(gnuplotPipe, "e\n"); // Конец данных
    pclose(gnuplotPipe);
    logPlotEnd("HeightMap3D");
}

void plotVoronoi(const std::unique_ptr<Pole>& p, 
                    const std::vector<VoronoiEdge>& edges, 
                    const std::vector<PointD>& sites, 
                    const std::string& filename) {
        logPlotStart("VoronoiDiagram", filename);
        
        if (!p || edges.empty()) {
            logger.warning("[GnuplotInterface::plotVoronoi] Нет данных для визуализации");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Диаграмма Вороного: ") + std::to_string(edges.size()) + " ребер, " + 
                   std::to_string(sites.size()) + " сайтов");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotVoronoi] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика с контрастными цветами для красного фона
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Voronoi Diagram on Red Field'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема:
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#00FF00' lw 2    # Ярко-зеленые ребра\n");
    fprintf(gnuplotPipe, "set style line 2 lc rgb '#FFFFFF' pt 7 ps 2 # Белые центры\n");

    // Многослойный график:
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines ls 1, \\\n");  // Ребра
    fprintf(gnuplotPipe, "'-' with points ls 2\n");     // Центры

    // 1. Данные поля (красный фон)
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Ребра Вороного (ярко-зеленые)
    for (const auto& edge : edges) {
        fprintf(gnuplotPipe, "%f %f\n%f %f\n\n", 
                edge.start.x, transformY(edge.start.y, height),
                edge.end.x, transformY(edge.end.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 3. Центры (белые точки)
    for (const auto& site : sites) {
        fprintf(gnuplotPipe, "%f %f\n", 
                site.x, 
                transformY(site.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("VoronoiDiagram");
}
    
    void plotDelaunay(const std::vector<Triangle>& triangles, 
                     std::unique_ptr<Pole>& p, 
                     const std::string& filename) {
        logPlotStart("DelaunayTriangulation", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotDelaunay] Ошибка открытия gnuplot pipe");
            return;
        }
        
        if (!p) {
            logger.error("[GnuplotInterface::plotDelaunay] Нет данных высот");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Триангуляция Делоне: ") + std::to_string(triangles.size()) + " треугольников");

    // Улучшенные настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Delaunay Triangulation'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#00FF00' lw 1.5\n");   // Ярко-зеленые линии
    
    // Многослойный график: фон + треугольники
    fprintf(gnuplotPipe, "plot '-' matrix with image, '-' with lines ls 1\n");

    // 1. Данные фона (с инверсией Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Данные треугольников (с инверсией Y)
    for (const auto& tri : triangles) {
        fprintf(gnuplotPipe, "%f %f\n", tri.a.x, transformY(tri.a.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.b.x, transformY(tri.b.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.c.x, transformY(tri.c.y, height));
        fprintf(gnuplotPipe, "%f %f\n\n", tri.a.x, transformY(tri.a.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("DelaunayTriangulation");
}

void plotPath(const std::vector<PointD>& path, 
                 const std::unique_ptr<Pole>& p, 
                 const std::string& filename, 
                 const DispatcherParams& params, 
                 const PathFinder& pathFinder, 
                 const double Radius) {
        logPlotStart("PathVisualization", filename);
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plotPath] Нет данных пути для визуализации");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Визуализация пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Path Visualization'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // 1. Собираем все пиксели пути
    std::vector<PointD> pathPixels;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto segment = pathFinder.bresenhamLine(path[i], path[i+1]);
        pathPixels.insert(pathPixels.end(), segment.begin(), segment.end());
    }

    // 2. Создаем команды для окружностей
    for (const auto& pixel : pathPixels) {
        const int x = static_cast<int>(pixel.x);
        const int y = static_cast<int>(pixel.y);
        const double currentHeight = p->field[y][x];
        const double heightDiff = std::abs(currentHeight - params.threshold);
        
        if (heightDiff < Radius) {
            const double radius = std::sqrt(
                Radius * Radius - 
                heightDiff * heightDiff
            );
            
            fprintf(gnuplotPipe, "set object circle at %d,%d size %f fc rgb '#004400' fs transparent solid 0.5 front\n",
                   static_cast<int>(pixel.x), static_cast<int>(transformY(pixel.y, height)), radius);
        }
    }
    
    // 3. Подписи START/END
    fprintf(gnuplotPipe, "set label 'START' at %d,%d front\n", 
            static_cast<int>(params.startPointX), static_cast<int>(transformY(params.startPointY, height)));
    fprintf(gnuplotPipe, "set label 'END' at %d,%d front\n",
            static_cast<int>(params.endPointX), static_cast<int>(transformY(params.endPointY, height)));
    
    // Многослойный график: фон + путь + точки
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines lw 2 lc rgb '#00FF0080', \\\n"); // Светло-зеленый 
    fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 9 ps 2 lc 'purple'\n");

    // 4. Данные фона (инверсия Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 5. Путь (инверсия Y)
    for (const auto& point : path) {
        fprintf(gnuplotPipe, "%f %f\n", 
                point.x, 
                transformY(point.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 6. Точка A (синяя)
    fprintf(gnuplotPipe, "%f %f\n", 
            params.startPointX, 
            transformY(params.startPointY, height));
    fprintf(gnuplotPipe, "e\n");

    // 7. Точка B (фиолетовая)
    fprintf(gnuplotPipe, "%f %f\n", 
            params.endPointX, 
            transformY(params.endPointY, height));
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("PathVisualization");
}

void plotInteractive3DPath(const std::vector<PointD>& path, 
                              const std::unique_ptr<Pole>& p, 
                              const PointD& start, 
                              const PointD& end, 
                              const PathFinder& pathFinder,  
                              const double sphereRadius) {
        logPlotStart("Interactive3DPath", "интерактивное окно");
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plotInteractive3DPath] Нет данных пути");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Интерактивный 3D путь: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotInteractive3DPath] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // 1. Настройки графика (должны идти первыми!)
    fprintf(gnuplotPipe, "set terminal qt size 1600,1200 enhanced\n");
    fprintf(gnuplotPipe, "set title '3D Path Visualization (Use mouse to rotate)'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set zlabel 'Height'\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    fprintf(gnuplotPipe, "set mouse\n");  // Включаем интерактивное управление
    fprintf(gnuplotPipe, "set key outside\n");
    
    // 2. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = pathFinder.bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX], sphereRadius);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX], sphereRadius);
    }
    
     fprintf(gnuplotPipe, "set label 'END' at %d,%d,%f front\n",
                endX, static_cast<int>(transformY(endY, height)), p->field[endY][endX] + 1);

    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "pause mouse close\n");
    pclose(gnuplotPipe);
    logPlotEnd("Interactive3DPath");
}

void plot3DPath(const std::vector<PointD>& path, 
                   const std::unique_ptr<Pole>& p, 
                   const std::string& filename, 
                   const PointD& start, 
                   const PointD& end, 
                   const PathFinder& pathFinder, 
                   const double sphereRadius) {
        logPlotStart("3DPathProjection", filename);
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plot3DPath] Нет данных пути");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("3D проекция пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plot3DPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // 1. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = pathFinder.bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }

    // 2. Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo enhanced size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title '3D Path Projection'\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX] , sphereRadius);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX], sphereRadius);
    }
    
     fprintf(gnuplotPipe, "set label 'END' at %d,%d,%f front\n",
                endX, static_cast<int>(transformY(endY, height)), p->field[endY][endX] + 1);

    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("3DPathProjection");
}

};
