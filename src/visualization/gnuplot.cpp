#include "gnuplot.hpp"

namespace visualization {

    GnuplotInterface::GnuplotInterface(core::Logger& lg) : logger(lg) {
        logger.trace("[GnuplotInterface] Инициализация интерфейса Gnuplot");
    }
    
    void GnuplotInterface::applyNiceStyle(FILE* pipe, const std::string& title, bool is3D, const std::string& terminal) {

    if (terminal == "qt") {
        fprintf(pipe,
            "set terminal qt enhanced font 'Arial,20'\n");
    }
    
    else if (terminal == "gif")
    {
        fprintf(pipe,
            "set terminal gif animate delay 5 optimize size 1600,1200\n");
    }
    
    else {
        fprintf(pipe,
            "set terminal pngcairo enhanced size 1600,1200 font 'Arial,20'\n");
    }

    fprintf(pipe,
        "set title '{/:Bold %s}'\n", title.c_str());

    fprintf(pipe,
        "set xlabel '{/:Bold X coordinate}' font 'Arial,18'\n");
    fprintf(pipe,
        "set ylabel '{/:Bold Y coordinate}' font 'Arial,18'\n");

    fprintf(pipe,
        "set xtics font 'Arial,18'\n");
    fprintf(pipe,
        "set ytics font 'Arial,18'\n");

    if (is3D) {
        fprintf(pipe,
            "set zlabel '{/:Bold Height}' font 'Arial,18'\n");
        fprintf(pipe,
            "set ztics font 'Arial,18'\n");
    }

   fprintf(pipe, "set grid\n"); // бонус
}

    void GnuplotInterface::logPlotStart(const std::string& plotType, const std::string& filename) const {
        logger.info(std::string("[GnuplotInterface] Начало визуализации: ") + plotType);
        logger.debug(std::string("Файл вывода: ") + filename);
    }

    void GnuplotInterface::logPlotEnd(const std::string& plotType) const {
        logger.info(std::string("[GnuplotInterface] Визуализация завершена: ") + plotType);
    }

    void GnuplotInterface::plotBinaryWithComponents(const std::vector<std::vector<double>>& binaryMap, 
                                const std::vector<algorithms::components::Component>& components, 
                                const std::string& filename) {
        logPlotStart("BinaryWithComponents", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotBinaryWithComponents] Ошибка открытия gnuplot pipe");
            return;
        }

        const int height = static_cast<int>(binaryMap.size());
        const int width = static_cast<int>(binaryMap[0].size());
        logger.debug("Число компонент: " + std::to_string(components.size()));

        // Настройки графика
        applyNiceStyle(gnuplotPipe, "Binary Image with Components");
        fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
        fprintf(gnuplotPipe, "set size ratio -1\n");
        fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
        fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
        fprintf(gnuplotPipe, "unset key\n");

        // Многослойный график
        fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
        fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'red', \\\n");
        fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
        fprintf(gnuplotPipe, "'-' with vectors head filled  lw 3 lc 'green'\n");

        // 1. Данные бинарного изображения (с инверсией Y)
for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
        fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
    }
            fprintf(gnuplotPipe, "\n");
        }
        fprintf(gnuplotPipe, "e\n");

         // 2. Границы компонент (красные прямоугольники)
        for (const auto& comp : components) {
            double y_min = comp.min_y;
            double y_max = comp.max_y;
            
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
                    comp.center_y);
        }
        fprintf(gnuplotPipe, "e\n");

         // 4. Собственные векторы (зеленые стрелки)
for (const auto& comp : components) {
    double cy = comp.center_y;
    
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
      void GnuplotInterface::gnuplot(const std::vector<std::vector<double>>& field, const std::string& filename) {
        logPlotStart("HeightMap3D", filename);
        
        if (field.empty()) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка: данные высот не инициализированы");
            return;
        }

        int height = static_cast<int>(static_cast<int>(field.size()));
        int width = static_cast<int>(static_cast<int>(field[0].size()));
        logger.debug(std::string("Размер сетки: ") + std::to_string(width) + "x" + std::to_string(height));
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // Настройки 3D графика
    applyNiceStyle(gnuplotPipe, "Height Map", true);
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "set zrange [*:*]\n"); // Автомасштабирование по Z
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n"); // Угол обзора
    
    // Формат данных: x y z
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Height Map'\n");
    
    // Записываем данные в правильном порядке
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, y, field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n"); // Пустая строка между слоями Y
    }
    
    fprintf(gnuplotPipe, "e\n"); // Конец данных
    pclose(gnuplotPipe);
    logPlotEnd("HeightMap3D");
}

void GnuplotInterface::plotKmeans(const std::vector<std::vector<double>>& binaryMap, const std::vector<algorithms::geometry::PointD>& centers, const std::string& filename) {
       logPlotStart("Kmeans", filename);
        
        if (binaryMap.empty() || centers.empty()) {
            logger.warning("[GnuplotInterface::plotKmeans] Нет данных для визуализации");
            return;
        }

        const int height = static_cast<int>(binaryMap.size());
        const int width = static_cast<int>(binaryMap[0].size());
        logger.debug(std::string("Kmeans: ") + std::to_string(centers.size()) + " центров");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotKmeans] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика с контрастными цветами
    applyNiceStyle(gnuplotPipe, "Kmeans");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема:
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#FFFFFF' pt 7 ps 2 # Белые центры\n");

    // Многослойный график:
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with points ls 1\n");     // Центры

    // 1. Данные поля
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Центры (белые точки)
    for (const auto& center : centers) {
        fprintf(gnuplotPipe, "%f %f\n", 
                center.x, 
                center.y);
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("Kmeans");
}

void GnuplotInterface::plotVoronoi(const std::vector<std::vector<double>>& binaryMap, 
                    const std::vector<algorithms::geometry::Edge>& edges, 
                    const std::vector<algorithms::geometry::PointD>& sites, 
                    const std::string& filename) {
        logPlotStart("VoronoiDiagram", filename);
        
        if (binaryMap.empty() || edges.empty()) {
            logger.warning("[GnuplotInterface::plotVoronoi] Нет данных для визуализации");
            return;
        }

        const int height = static_cast<int>(binaryMap.size());
        const int width = static_cast<int>(binaryMap[0].size());
        logger.debug(std::string("Диаграмма Вороного: ") + std::to_string(edges.size()) + " ребер, " + 
                   std::to_string(sites.size()) + " центров");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotVoronoi] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика с контрастными цветами для красного фона
    applyNiceStyle(gnuplotPipe, "Voronoi Diagram");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
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
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Ребра Вороного (ярко-зеленые)
    for (const auto& edge : edges) {
        fprintf(gnuplotPipe, "%f %f\n%f %f\n\n", 
                edge.a.x, edge.a.y,
                edge.b.x, edge.b.y);
    }
    fprintf(gnuplotPipe, "e\n");

    // 3. Центры (белые точки)
    for (const auto& site : sites) {
        fprintf(gnuplotPipe, "%f %f\n", 
                site.x, 
                site.y);
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("VoronoiDiagram");
}

void GnuplotInterface::plotGraph(
    const std::vector<std::vector<double>>& binaryMap,
    const std::unordered_map<
        algorithms::geometry::Pixel,
        std::vector<algorithms::geometry::Pixel>>& graph,
    const std::string& filename,
    std::optional<algorithms::geometry::Pixel> start,
    std::optional<algorithms::geometry::Pixel> goal)
{
    logPlotStart("NavigationGraph", filename);

    if (binaryMap.empty() || graph.empty()) {
        logger.warning(
            "[GnuplotInterface::plotGraph] Нет данных для визуализации");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width = static_cast<int>(binaryMap[0].size());

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        logger.error(
            "[GnuplotInterface::plotGraph] Ошибка открытия gnuplot pipe");
        return;
    }

    applyNiceStyle(gnuplotPipe, "Navigation Graph");

    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");
    
     // Подписи START/END
     if (start)
        fprintf(gnuplotPipe, "set label '{/:Bold START}' at %d,%d tc rgb 'blue' front font 'Arial,18'\n", 
                static_cast<int>(start->x), static_cast<int>(start->y));
    if (goal)
        fprintf(gnuplotPipe, "set label '{/:Bold END}' at %d,%d tc rgb 'purple' front font 'Arial,18'\n",
                static_cast<int>(goal->x), static_cast<int>(goal->y));
            
    fprintf(gnuplotPipe,
            "set style line 1 lc rgb '#00FFFF' lw 2\n"); // ребра
    fprintf(gnuplotPipe,
            "set style line 2 lc rgb '#FFFFFF' pt 7 ps 1.5\n"); // вершины
    fprintf(gnuplotPipe,
            "plot '-' matrix with image, \\\n"
            "'-' with lines ls 1, \\\n"
            "'-' with points ls 2");
if (start)
    fprintf(gnuplotPipe,
            ", \\\n'-' with points pt 7 ps 2 lc rgb 'blue'");

if (goal)
    fprintf(gnuplotPipe,
            ", \\\n'-' with points pt 9 ps 2 lc rgb 'purple'");

fprintf(gnuplotPipe, "\n");
    // ==================================================
    // Поле
    // ==================================================

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }

    fprintf(gnuplotPipe, "e\n");

    // ==================================================
    // Ребра графа 
    // ==================================================

    for (const auto& [node, neighbors] : graph)
    {
        for (const auto& neighbor : neighbors)
        {
            algorithms::geometry::PixelEdge e(node, neighbor);

            fprintf(gnuplotPipe,
                    "%d %d\n"
                    "%d %d\n\n",
                    e.a.x,
                    e.a.y,
                    e.b.x,
                    e.b.y);
        }
    }

    fprintf(gnuplotPipe, "e\n");

    // ==================================================
    // Вершины графа
    // ==================================================

    for (const auto& [node, _] : graph)
    {
        fprintf(gnuplotPipe,
                "%d %d\n",
                node.x,
                node.y);
    }

    fprintf(gnuplotPipe, "e\n");

if (start)
{
    fprintf(gnuplotPipe,
            "%d %d\n",
            start->x,
            start->y);
    fprintf(gnuplotPipe, "e\n");
}

if (goal)
{
    fprintf(gnuplotPipe,
            "%d %d\n",
            goal->x,
            goal->y);
    fprintf(gnuplotPipe, "e\n");
}
    
    pclose(gnuplotPipe);

    logPlotEnd("NavigationGraph");
}

void GnuplotInterface::plotGrid(
    const algorithms::path::common::Grid& grid,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename)
{
    logPlotStart("Grid", filename);

    if (binaryMap.empty() || grid.cells.empty()) {
        logger.warning("[plotGrid] пустые данные");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        logger.error("[plotGrid] Не удалось открыть gnuplot");
        return;
    }

    applyNiceStyle(gnuplotPipe, "Grid");

    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");

    fprintf(gnuplotPipe,
            "plot '-' matrix with image, \\\n"
            "'-' with lines lc rgb 'white' lw 1, \\\n"
            "'-' with lines lc rgb 'white' lw 1\n");

    // =====================================================
    // 1. Карта
    // =====================================================

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // =====================================================
    // 2. Вертикальные линии
    // =====================================================

    const int cs = grid.cellSize;

    for (int col = 0; col <= grid.cols; ++col)
    {
        const int x = col * cs;

        fprintf(gnuplotPipe,
                "%d %d\n"
                "%d %d\n\n",
                x, 0,
                x, static_cast<int>(height - 1));
    }

    fprintf(gnuplotPipe, "e\n");

    // =====================================================
    // 3. Горизонтальные линии
    // =====================================================

    for (int row = 0; row <= grid.rows; ++row)
    {
        int y = row * cs;

        fprintf(gnuplotPipe,
                "%d %d\n"
                "%d %d\n\n",
                0, y,
                static_cast<int>(width - 1), y);
    }

    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);

    logPlotEnd("Grid");
}

void GnuplotInterface::plotNavGrid(
    const algorithms::path::common::Grid& grid,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename,
    const std::optional<algorithms::path::common::GridCell> start,
    const std::optional<algorithms::path::common::GridCell> end)
{
    logPlotStart("NavGrid", filename);

    if (binaryMap.empty() || grid.cells.empty()) {
        logger.warning("[plotNavGrid] empty input");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());
    const int cs     = grid.cellSize;
    FILE* pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
        logger.error("[plotNavGrid] gnuplot pipe error");
        return;
    }
  
    applyNiceStyle(pipe, "NavGrid");
    fprintf(pipe, "set output '%s'\n", filename.c_str());
    fprintf(pipe, "unset key\n");
    fprintf(pipe, "set size ratio -1\n");
    fprintf(pipe, "set xrange [0:%d]\n", width - 1);
    fprintf(pipe, "set yrange [0:%d]\n", height - 1);

    int id = 1;

    for (const auto& cell : grid.cells)
    {
        double x1 = cell.col * cs;
        double x2 = x1 + cs;

        double y1 = cell.row * cs + cs;
        double y2 = cell.row * cs;

        std::string color;

        bool isStart = start && cell.row == start->row && cell.col == start->col;
        bool isEnd   = end   && cell.row == end->row   && cell.col == end->col;

        if (isStart || isEnd)
            color = "#A020F0"; // фиолетовый
        else if (!cell.traversable)
            color = "#FF0000"; // красный
        else
            color = "#00FF00"; // зелёный

        fprintf(pipe,
            "set object %d rect from %f,%f to %f,%f fc rgb '%s' fs solid border lc rgb 'black'\n",
            id++, x1, y1, x2, y2, color.c_str());
        
        if (isStart)
        {
            fprintf(pipe,
                "set label 'START' at %f,%f center tc rgb 'black' font ',10'\n",
                x1 + cs * 0.5, y1 + cs * 0.5);
        }

        if (isEnd)
        {
            fprintf(pipe,
                "set label 'END' at %f,%f center tc rgb 'black' font ',10'\n",
                x1 + cs * 0.5, y1 + cs * 0.5);
        }
    }

    fprintf(pipe, "plot NaN notitle\n");
    pclose(pipe);

    logPlotEnd("NavGrid");
}

void GnuplotInterface::plotGridPath(
    const algorithms::path::common::Grid& grid,
    std::vector<algorithms::path::common::GridCell>& gridPath,
    const algorithms::path::common::GridCell& start,
    const algorithms::path::common::GridCell& end,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename)
{
    logPlotStart("GridPath", filename);

    if (grid.cells.empty() || binaryMap.empty()) {
        logger.warning("[PlotGridPath] empty input");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());
    const int cs     = grid.cellSize;

    FILE* gnuplotPipe = popen("gnuplot", "w");
    if (!gnuplotPipe) {
        logger.error("[PlotGridPath] gnuplot pipe error");
        return;
    }

    applyNiceStyle(gnuplotPipe, "GridPath");
    
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");

    // =====================================================
    // 1. PATH lookup (быстро)
    // =====================================================
    auto isInPath = [&](const algorithms::path::common::GridCell& c)
    {
        return std::any_of(gridPath.begin(), gridPath.end(),
            [&](const algorithms::path::common::GridCell& p)
            {
                return p.row == c.row && p.col == c.col;
            });
    };

    // =====================================================
    // 2. КЛЕТКИ
    // =====================================================

    int objectId = 1;

    for (const auto& cell : grid.cells)
    {
        double x1 = cell.col * cs;
        double x2 = x1 + cs;

        double y1 = cell.row * cs + cs;
        double y2 = cell.row * cs;

        std::string color;

        bool isStart = (cell.row == start.row && cell.col == start.col);
        bool isEnd   = (cell.row == end.row && cell.col == end.col);
        bool inPath  = isInPath(cell);

        if (isStart || isEnd || inPath)
            color = "#A020F0"; // фиолетовый
        else if (!cell.traversable)
            color = "#FF0000"; // красный
        else
            color = "#00FF00"; // зелёный

        fprintf(gnuplotPipe,
            "set object %d rect from %f,%f to %f,%f "
            "fc rgb '%s' fs border lc rgb 'black'\n",
            objectId++, x1, y1, x2, y2, color.c_str());

        // =====================================================
        // 3. LABELS
        // =====================================================

        if (isStart)
        {
            fprintf(gnuplotPipe,
                "set label 'START' at %f,%f center tc rgb 'black' font ',10'\n",
                x1 + cs * 0.5, y1 + cs * 0.5);
        }

        if (isEnd)
        {
            fprintf(gnuplotPipe,
                "set label 'END' at %f,%f center tc rgb 'black' font ',10'\n",
                x1 + cs * 0.5, y1 + cs * 0.5);
        }
    }

    // =====================================================
    // 4. финальный render
    // =====================================================
    fprintf(gnuplotPipe, "plot NaN notitle\n");
    pclose(gnuplotPipe);

    logPlotEnd("GridPath");
}
               
    void GnuplotInterface::plotDelaunay(const std::vector<algorithms::geometry::Triangle>& triangles, 
                     const std::vector<std::vector<double>>& binaryMap, 
                     const std::string& filename) {
        logPlotStart("DelaunayTriangulation", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotDelaunay] Ошибка открытия gnuplot pipe");
            return;
        }
        
        if (binaryMap.empty()) {
            logger.error("[GnuplotInterface::plotDelaunay] Нет данных высот");
            return;
        }
        
         if (triangles.empty()) {
            logger.error("[GnuplotInterface::plotDelaunay] Нет треугольников");
            return;
        }

        const int height = static_cast<int>(binaryMap.size());
        const int width = static_cast<int>(binaryMap[0].size());
        logger.debug(std::string("Триангуляция Делоне: ") + std::to_string(triangles.size()) + " треугольников");

    // Улучшенные настройки графика
    applyNiceStyle(gnuplotPipe, "Delaunay Triangulation");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#00FF00' lw 1.5\n");   // Ярко-зеленые линии
    
    // Многослойный график: фон + треугольники
    fprintf(gnuplotPipe, "plot '-' matrix with image, '-' with lines ls 1\n");

    // 1. Данные фона (с инверсией Y)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Данные треугольников (с инверсией Y)
    for (const auto& tri : triangles) {
        fprintf(gnuplotPipe, "%f %f\n", tri.a.x, tri.a.y);
        fprintf(gnuplotPipe, "%f %f\n", tri.b.x, tri.b.y);
        fprintf(gnuplotPipe, "%f %f\n", tri.c.x, tri.c.y);
        fprintf(gnuplotPipe, "%f %f\n\n", tri.a.x, tri.a.y);
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("DelaunayTriangulation");
}

void GnuplotInterface::plotPath(const std::vector<algorithms::geometry::Pixel>& path, 
                 const std::vector<std::vector<double>>& field,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename, 
                 const command::DispatcherParams& params,
                 const int Radius) {
        logPlotStart("PathVisualization", filename);
        
        if (field.empty() || path.empty()) {
            logger.warning("[GnuplotInterface::plotPath] Нет данных пути для визуализации");
            return;
        }

        const int height = static_cast<int>(field.size());
        const int width = static_cast<int>(field[0].size());
        logger.debug(std::string("Визуализация пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика
    applyNiceStyle(gnuplotPipe, "Path Visualization");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // 1. Собираем все пиксели пути
    std::vector<algorithms::geometry::Pixel> pathPixels;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto segment = algorithms::geometry::bresenhamLine(path[i], path[i+1]);
        pathPixels.insert(pathPixels.end(), segment.begin(), segment.end());
    }

    // 2. Создаем команды для окружностей
    for (const auto& pixel : pathPixels) {
        const int x = static_cast<int>(pixel.x);
        const int y = static_cast<int>(pixel.y);
        const double currentHeight = field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)];
        const double heightDiff = std::abs(currentHeight - core::MID_GRAY - params.heightThresholdPixel);
        
        if (heightDiff < Radius) {
            const double radius = std::sqrt(
                Radius * Radius - 
                heightDiff * heightDiff
            );
            
            fprintf(gnuplotPipe, "set object circle at %d,%d size %f fc rgb '#004400' fs transparent solid 0.5 front\n",
                   static_cast<int>(pixel.x), static_cast<int>(pixel.y), radius);
        }
    }
    
    // 3. Подписи START/END
    fprintf(gnuplotPipe, "set label '{/:Bold START}' at %d,%d tc rgb 'blue' front font 'Arial,18'\n", 
            static_cast<int>(params.startPixelX), static_cast<int>(params.startPixelY));
    fprintf(gnuplotPipe, "set label '{/:Bold END}' at %d,%d tc rgb 'purple' front font 'Arial,18'\n",
            static_cast<int>(params.goalPixelX), static_cast<int>(params.goalPixelY));
    
    // Многослойный график: фон + путь + точки
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines lw 2 lc rgb '#00FF0080', \\\n"); // Светло-зеленый 
    fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 9 ps 2 lc 'purple'\n");

    // 4. Данные фона (инверсия Y)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 5. Путь (инверсия Y)
    for (const auto& point : path) {
        fprintf(gnuplotPipe, "%d %d\n", 
                point.x, 
                point.y);
    }
    fprintf(gnuplotPipe, "e\n");

    // 6. Точка A (синяя)
    fprintf(gnuplotPipe, "%d %d\n", 
            params.startPixelX, 
            params.startPixelY);
    fprintf(gnuplotPipe, "e\n");

    // 7. Точка B (фиолетовая)
    fprintf(gnuplotPipe, "%d %d\n", 
            params.goalPixelX, 
            params.goalPixelY);
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("PathVisualization");
}

void GnuplotInterface::plotRRT(
        const std::vector<algorithms::geometry::PointD>& path,
        const std::vector<algorithms::path::rrt::RRTNode>& treeRRT,
        const algorithms::geometry::PointD& start,
        const algorithms::geometry::PointD& end,
        const std::vector<std::vector<double>>& binaryMap,
        const std::string& filename)
{
    logPlotStart("RRT Animation", filename);

    if (binaryMap.empty() || treeRRT.empty())
    {
        logger.warning("[plotRRT] Empty data");
        return;
    }

    FILE* pipe = popen("gnuplot", "w");

    if (!pipe)
    {
        logger.error("[plotRRT] Cannot start gnuplot");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());

    applyNiceStyle(pipe,
                   "Rapidly Exploring Random Tree",
                   false,
                   "gif");

    fprintf(pipe,
        "set output '%s'\n",
        filename.c_str());

    fprintf(pipe, "unset key\n");
    fprintf(pipe, "set size ratio -1\n");

    fprintf(pipe,
        "set xrange [0:%d]\n"
        "set yrange [0:%d]\n",
        width - 1,
        height - 1);

    fprintf(pipe,
        "set label 1 'Start' at %f,%f tc rgb 'blue' front\n",
        start.x,
        start.y);

    fprintf(pipe,
        "set label 2 'Goal' at %f,%f tc rgb 'purple' front\n",
        end.x,
        end.y);

    // ============================
    // Постепенное построение дерева
    // ============================
    
    std::size_t step = std::max<std::size_t>(1, treeRRT.size() / 50);
        
    for (std::size_t frame = 1; frame <= treeRRT.size(); frame += step)
    {
        fprintf(pipe,
            "plot '-' matrix with image,\\\n"
            "'-' with lines lw 2 lc rgb '#00AA00',\\\n"
            "'-' with points pt 7 ps 2 lc rgb 'blue',\\\n"
            "'-' with points pt 9 ps 2 lc rgb 'purple'\n");

        //-------------------------
        // Карта
        //-------------------------

        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                fprintf(pipe,
                        "%f ",
                        binaryMap[static_cast<std::size_t>(y)]
                                 [static_cast<std::size_t>(x)]);
            }
            fprintf(pipe, "\n");
        }
        fprintf(pipe, "e\n");

        //-------------------------
        // Дерево
        //-------------------------

        for (std::size_t i = 1; i < frame; ++i)
        {
            const auto& node = treeRRT[i];

            if (node.parent < 0)
                continue;

            const auto& parent =
                treeRRT[static_cast<std::size_t>(node.parent)];

            fprintf(pipe,
                "%f %f\n",
                parent.point.x,
                parent.point.y);

            fprintf(pipe,
                "%f %f\n\n",
                node.point.x,
                node.point.y);
        }
        
        fprintf(pipe, "e\n");

        //-------------------------
        // Старт
        //-------------------------

        fprintf(pipe,
            "%f %f\n"
            "e\n",
            start.x,
            start.y);

        //-------------------------
        // Финиш
        //-------------------------

        fprintf(pipe,
            "%f %f\n"
            "e\n",
            end.x,
            end.y);
    }

    // =====================================
    // Последний кадр — рисуем путь поверх
    // =====================================
    for (int pause = 0; pause < 50; ++pause) 
    {
        fprintf(pipe,
            "plot '-' matrix with image,\\\n"
            "'-' with lines lw 2 lc rgb '#00AA00',\\\n"
            "'-' with lines lw 4 lc rgb 'red',\\\n"
            "'-' with points pt 7 ps 2 lc rgb 'blue',\\\n"
            "'-' with points pt 9 ps 2 lc rgb 'purple'\n");

        // карта
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
                fprintf(pipe,
                        "%f ",
                        binaryMap[static_cast<std::size_t>(y)]
                                 [static_cast<std::size_t>(x)]);
            fprintf(pipe, "\n");
        }
        fprintf(pipe, "e\n");

        // дерево
        for (std::size_t i = 1; i < treeRRT.size(); ++i)
        {
            if (treeRRT[i].parent < 0)
                continue;

            const auto& p =
                treeRRT[static_cast<std::size_t>(treeRRT[i].parent)];

            fprintf(pipe,
                "%f %f\n"
                "%f %f\n\n",
                p.point.x,
                p.point.y,
                treeRRT[i].point.x,
                treeRRT[i].point.y);
        }
        fprintf(pipe, "e\n");

        // путь
        for (const auto& p : path)
            fprintf(pipe,
                    "%f %f\n",
                    p.x,
                    p.y);
        fprintf(pipe, "e\n");

        // старт
        fprintf(pipe,
            "%f %f\n"
            "e\n",
            start.x,
            start.y);

        // финиш
        fprintf(pipe,
            "%f %f\n"
            "e\n",
            end.x,
            end.y);
    }
    pclose(pipe);

    logPlotEnd("RRT Animation");
}

void GnuplotInterface::plotInteractive3DPath(const std::vector<algorithms::geometry::Pixel>& path, 
                              const std::vector<std::vector<double>>& field, 
                              const algorithms::geometry::Pixel& start, 
                              const algorithms::geometry::Pixel& end,
                              const int sphereRadius) {
        logPlotStart("Interactive3DPath", "интерактивное окно");
        
        if (field.empty() || path.empty()) {
            logger.warning("[GnuplotInterface::plotInteractive3DPath] Нет данных пути");
            return;
        }

        const int height = static_cast<int>(field.size());
        const int width = static_cast<int>(field[0].size());
        logger.debug(std::string("Интерактивный 3D путь: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotInteractive3DPath] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // 1. Настройки графика (должны идти первыми!)
    applyNiceStyle(gnuplotPipe, "3D Path Visualization", true, "qt");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    fprintf(gnuplotPipe, "set mouse\n");  // Включаем интерактивное управление
    fprintf(gnuplotPipe, "set key outside\n");
    
    // 2. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = algorithms::geometry::bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(y), field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(startY), 
                field[static_cast<std::size_t>(startY)][static_cast<std::size_t>(startX)], sphereRadius);
    
    fprintf(gnuplotPipe, "set label '{/:Bold START}' at %d,%d,%f tc rgb '#FF0000' front font 'Arial,18'\n",
                startX, static_cast<int>(startY), field[static_cast<std::size_t>(startY)][static_cast<std::size_t>(startX)] + 1);
    }
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(endY), 
                field[static_cast<std::size_t>(endY)][static_cast<std::size_t>(endX)], sphereRadius);
    
     fprintf(gnuplotPipe, "set label '{/:Bold END}' at %d,%d,%f tc rgb '#0000FF' front font 'Arial,18'\n",
                endX, static_cast<int>(endY), field[static_cast<std::size_t>(endY)][static_cast<std::size_t>(endX)] + 1);
    }
    
    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "pause mouse close\n");
    pclose(gnuplotPipe);
    logPlotEnd("Interactive3DPath");
}

void GnuplotInterface::plot3DPath(const std::vector<algorithms::geometry::Pixel>& path, 
                   const std::vector<std::vector<double>>& field, 
                   const std::string& filename, 
                   const algorithms::geometry::Pixel& start, 
                   const algorithms::geometry::Pixel& end,
                   const int sphereRadius) {
        logPlotStart("3DPathProjection", filename);
        
        if (field.empty() || path.empty()) {
            logger.warning("[GnuplotInterface::plot3DPath] Нет данных пути");
            return;
        }

        const int height = static_cast<int>(field.size());
        const int width = static_cast<int>(field[0].size());
        logger.debug(std::string("3D проекция пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plot3DPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // 1. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = algorithms::geometry::bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }

    // 2. Настройки графика
    applyNiceStyle(gnuplotPipe, "3D Path Projection", true);
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(y), field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(startY), 
                field[static_cast<std::size_t>(startY)][static_cast<std::size_t>(startX)] , sphereRadius);
    
    fprintf(gnuplotPipe, "set label '{/:Bold START}' at %d,%d,%f tc rgb '#FF0000' front font 'Arial,18'\n",
                startX, static_cast<int>(startY), field[static_cast<std::size_t>(startY)][static_cast<std::size_t>(startX)] + 1);
    }
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %d fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(endY), 
                field[static_cast<std::size_t>(endY)][static_cast<std::size_t>(endX)], sphereRadius);
    
     fprintf(gnuplotPipe, "set label '{/:Bold END}' at %d,%d,%f tc rgb '#0000FF' front font 'Arial,18'\n",
                endX, static_cast<int>(endY), field[static_cast<std::size_t>(endY)][static_cast<std::size_t>(endX)] + 1);
    }
    
    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("3DPathProjection");
}

}
