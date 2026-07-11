#include "command/interface.hpp"
#include "command/validator.hpp"

#include <iostream>
#include <sstream>

namespace command {

  // Приватная функция для вывода справки
    void Interface::showHelp() {
    logger.info("Showing help information");
    
    // Содержимое справки
const std::string helpContent = R"(# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки и ее радиуса
  

## 🚀 Запуск программы

### Способ 1: Автоматический скрипт (рекомендуется)
```bash
# Даем права на выполнение (только при первом запуске)
chmod +x run.sh

# Запуск с интерфейсом командной строки
./run.sh

# Запуск с файлом команд (filename)
./run.sh filename
```

## 🛠 Команды управления (для командного файла command.txt)

| Команда              | Тип данных + *параметры*                                                                           | Описание                                                                 |
|----------------------|----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| help                 | -                                                                                                  | Создание файла с пояснением команд                                       |
| init                 | int *fieldWidth fieldHeight*                                                                 | Инициализация поля с заданой шириной и длиной. field[y][x] X-(width) Y-(height)|
| g                    | double *x y sx sy h*                                                                               | Создает гаусс с центром (x,y), размерами (sx,sy) и высотой h             |
| g_auto               | int *count_min count_max* double *xmin xmax ymin ymax ... h_min h_max*                             | Создает случайные гаусы с параметрами в промежутках                      |
| generate             | -                                                                                                  | Складывает все добавленные гауссы в итоговое поле                        |
| save_g               | string *filename.txt*                                                                              | Сохраняет параметры g в txt файл                                         |
| gnuplot              | string *filename.png*                                                                              | Сохраняет 3D-визуализацию поля в PNG файл                                |
| PlotMetedata         | string *filename.png*                                                                              | Визуализирует метаданные компонент с границами и центрами                |
| PlotKmeans           | string *filename.png*                                                                              | Визуализирует k_means                                                    |
| PlotGraph            | string *filename.png*                                                                              | Визуализирует граф                                                       |
| PlotVoronoi          | string *filename.png*                                                                              | Строит диаграмму Вороного по текущей триангуляции                        |
| PlotDelaunay         | string *filename.png*                                                                              | Визуализирует триангуляцию Делоне                                        |
| PlotGrid             | string *filename.png*                                                                              | Визуализирует сетку                                                      |
| PlotNavGrid          | string *filename.png*                                                                              | Визуализирует навигационную сетку                                        |
| PlotGridPath         | string *filename.png*                                                                              | Визуализирует путь на сетке                                              |
| PlotPath             | string *filename.png*                                                                              | Отображает найденный путь между точками A и B                            |
| PlotRRT              | string *filename.png*                                                                              | Делает gif изображение построения дерева и пути RRT                      |
| bmp_write            | string *filename.bmp [Full/Binary]*                                                                | Сохраняет поле в BMP: Full - полное, Binary - бинаризованное             |
| bmp_read             | string *filename.bmp*                                                                              | Загружает поле из BMP файла                                              |
| bin                  | int *slice*                                                                                        | Бинаризация с уровнем отклонения от равнины MID_GRAY                     |
| wave                 | int *noisy*                                                                                        | Удаляет компоненты размером ≤ noisy как шум                              |
| k_means              | int *k*                                                                                            | Кластеризует данные в k кластеров                                        |
| k_means_kern         | int *k p*                                                                                          | Kmeans с параметром k на ядрах размера p                                 |
| triangulate          | -                                                                                                  | Строит триангуляцию Делоне по центрам компонент                          |
| voronoi              | -                                                                                                  | Строит диаграмму Вороного                                                |
| build_nav_graph      | int *vehicleRadius* double *maxSideAngle maxUpDownAngle*                                           | Строит граф проходимости по диаграмме Вороного                           |
| grid                 | int *width*                                                                                        | Строит квадратную сетку на поле (ячейка размера width x width)           |
| build_nav_grid       | int *noisy*                                                                                        | Удаляет ячейки с числом опасных пикселей >= noisy                        |
| connect_to_grid      | int *Ax Ay Bx By*                                                                                  | Подключает стартовую и финишную ячейку к сетке                           |
| connect_to_graph     | int *Ax Ay Bx By* string *[Nearest/NearestK/All]* int *k (только для NearestK)*                    | Подключает старт и финиш к графу с разными режимами                      |
| astar_graph          | -                                                                                                  | A* ищет путь между точками A и B на графе                                |
| dekstra_graph        | -                                                                                                  | Dekstra ищет путь между точками A и B на графе                           |
| greedy_graph         | -                                                                                                  | Greedy ищет путь между точками A и B на графе                            |
| astar_grid           | -                                                                                                  | A* ищет путь между точками A и B на сетке                                |
| dekstra_grid         | -                                                                                                  | Dekstra ищет путь между точками A и B на сетке                           |
| greedy_grid          | -                                                                                                  | Greedy ищет путь между точками A и B на сетке                            |
| save_metrics         | string *filename.csv*                                                                              | Сохраняет метрики в файле csv                                            |
| Plot3DPath           | string *filename.png*                                                                              | Сохраняет 3D-визуализацию пути в PNG файл                                |
| plotInteractive3DPath| -                                                                                                  | Интерактвный 3D режим с путем                                            |
| end                  | -                                                                                                  | Завершает работу программы                                               |
|rrt|size_t *maxIterations* double *Ax Ay Bx By vehicleRadius heightThreshold maxSideAngle maxUpDownAngle interpEdge interpCollision interpAngle step goalRadius goalBias*|Строит путь соблюдая условия|


## ⚙️ Параметры конфигурационного файла (config.txt)

| Параметр                     | Тип данных                                     | Описание                                                                         |
|------------------------------|------------------------------------------------|----------------------------------------------------------------------------------|
| defaultfieldWidth            | int                                            | Ширина рабочего поля в пикселях по умолчанию                                     |  
| defaultfieldHeight           | int                                            | Высота рабочего поля в пикселях по умолчанию                                     |
| defaultCenterX               | double                                         | Стандартная X-координата центра гауссова распределения по умолчанию              |
| defaultCenterY               | double                                         | Стандартная Y-координата центра гауссова распределения по умолчанию              |
| defaultSigmaX                | double                                         | Стандартное отклонение по оси X по умолчанию                                     |
| defaultSigmaY                | double                                         | Стандартное отклонение по оси Y по умолчанию                                     |
| defaultHeight                | double                                         | Стандартная высота гауссова распределения по умолчанию                           |
| count_min                    | int                                            | Минимальное количество генерируемых гауссов                                      |
| count_max                    | int                                            | Максимальное количество генерируемых гауссов                                     |
| xmin                         | double                                         | Минимальная X-координата центра гаусса при g_auto                                |
| xmax                         | double                                         | Максимальная X-координата центра гаусса при g_auto                               |
| ymin                         | double                                         | Минимальная Y-координата центра гаусса при g_auto                                |
| ymax                         | double                                         | Максимальная Y-координата центра гаусса при g_auto                               |
| sx_min                       | double                                         | Минимальное отклонение по оси X при g_auto                                       |
| sx_max                       | double                                         | Максимальное отклонение по оси X при g_auto                                      |
| sy_min                       | double                                         | Минимальное отклонение по оси Y при g_auto                                       |
| sy_max                       | double                                         | Максимальное отклонение по оси Y при g_auto                                      |
| h_min                        | double                                         | Минимальная высота гаусса при g_auto                                             |
| h_max                        | double                                         | Максимальная высота гаусса при g_auto                                            |
| save_g                       | string                                         | Путь к файлу для сохранения параметров g                                         |
| defaultGnuplot               | string                                         | Путь к файлу для сохранения 3D-визуализации по умолчанию                         |
| defaultPlotMetedata          | string                                         | Путь к файлу для визуализации метаданных компонент по умолчанию                  |
| defaultPlotKmeans            | string                                         | Путь к файлу для k_means по умолчанию                                            |
| defaultPlotGraph             | string                                         | Путь к файлу для графа по умолчанию                                              |
| defaultPlotVoronoi           | string                                         | Путь к файлу для диаграммы Вороного по умолчанию                                 |
| defaultPlotDelaunay          | string                                         | Путь к файлу для триангуляции Делоне по умолчанию                                |
| defaultPlotGrid              | string                                         | Путь к файлу для сетки по умолчанию                                              |
| defaultPlotNavGrid           | string                                         | Путь к файлу для навигационной сетки по умолчанию                                |
| defaultPlotGridPath          | string                                         | Путь к файлу для пути на сетке по умолчанию                                      |
| defaultPlotPath              | string                                         | Путь к файлу для визуализации маршрута по умолчанию                              |
| PlotRRT                      | string                                         | Путь к файлу для визуализации rrt по умолчанию                                   |
| defaultWrite                 | string                                         | Путь к файлу для сохранения BMP-изображения по умолчанию                         |
| defaultWriteModeImage        | [Full/Binary]                                  | Режим сохранения BMP (Full/Binary) по умолчанию                                  |
| defaultRead                  | string                                         | Путь к файлу для загрузки BMP-изображения по умолчанию                           |
| heightThresholdPixel         | int                                            | Порог бинаризации по умолчанию                                                   |
| defaultWaveNoisy             | int                                            | Порог для удаления шумовых компонент по умолчанию                                |
| defaultKlaster               | int                                            | Количество кластеров для k-mean по умолчанию                                     |
| defaultKlasterKern           | int                                            | Размер ядра для кластеризации по умолчанию                                       |
| defaultgridWidth             | int                                            | Размер ячейки для сетки по умолчанию                                             |
| defaultgridNoisy             | int                                            | Шум для опасных пикселей в ячейки по умолчанию                                   |
| startPixelX                  | int                                            | X-координата точки A для поиска пути по умолчанию                                |
| startPixelY                  | int                                            | Y-координата точки A для поиска пути по умолчанию                                |
| goalPixelX                   | int                                            | X-координата точки B для поиска пути по умолчанию                                |
| goalPixelY                   | int                                            | Y-координата точки B для поиска пути по умолчанию                                |
| defaultConnectMode           | [Nearest/NearestK int/All]                     | Режимы подключения старта и финиша к графу (Nearest/NearestK k/All) по умолчанию |
| defaultPlot3DPath            | string                                         | Путь к файлу для 3D-визуализации маршрута по умолчанию                           |
| vehicleRadiusPixel           | int                                            | Радиус транспортного средства                                                    |
| maxSideAngle                 | double                                         | Максимальный угол поворота вбок (градусы)                                        |
| maxUpDownAngle               | double                                         | Максимальный угол наклона вверх/вниз (градусы)                                   |
| maxIterations                | size_t                                         | Максимальное число итераций по умолчанию                                         |
| startWorldX                  | double                                         | X-координата точки A для поиска пути по умолчанию                                |
| startWorldY                  | double                                         | Y-координата точки A для поиска пути по умолчанию                                |
| goalWorldX                   | double                                         | X-координата точки B для поиска пути по умолчанию                                |
| goalWorldY                   | double                                         | Y-координата точки B для поиска пути по умолчанию                                |
| vehicleRadiusWorld           | double                                         | Радиус транспортного средства                                                    |
| heightThresholdWorld         | double                                         | Порог для отклонения от равнины по умолчанию                                     |
| interpEdge                   | double                                         | Шаг интерполяции ребра                                                           |
| interpCollision              | double                                         | Шаг интерполяции радиуса                                                         |
| interpAngle                  | double                                         | Шаг интерполяции для окружности (расстояние между соседними точками окружности)  |
| step                         | double                                         | Шаг алгоритма                                                                    |
| goalRadius                   | double                                         | Радиус круга цели, внутри которого пробуем присоединить оказавшиеся там вершины  |
| goalBias                     | double                                         | Вероятность оказаться случайной точке у цели                                     |
| defaultsave_metrics          | string                                         | Путь к файлу по умолчанию для сохранения метрик                                  |
| logFileNameInterface         | string                                         | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | string                                         | Путь к лог-файлу управления                                                      |
| defaultHelp                  | string                                         | Путь где сохранить help файл                                                     |
| FiltrationLogLevelInterface  | [TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF]  | Уровень логирования интерфейса                                                   |
| FiltrationLogLevelControl    | [TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF]  | Уровень логирования управления                                                   |
| seedMode                     | [Random/Fixed size_t]                          | Режим случайной генерации (Random/Fixed seed) по умолчанию                       |


## ⚠️ Важно
1. Всегда начинайте с команды init
2. Маршрут будет найден не всегда!
3. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
4. Если пользуетесь программой, то важно использовать ту же файловую структуру!
5. Путь к файлу пишем относительный
6. Важен порядок команд, не забывайте делать картинки после команд
7. Для команды bmp_write не полагайтесь на значения по умолчанию
8. Для команды PlotKmeans не полагайтесь на значения по умолчанию (иначе вывод совпадет для kmeans и kmeans_with_kernels)
8. Точки A и B должны попадать в безопасную зону иначе путь не будет найден
9. Уровень "равнины" = 127, чтобы метод записи поля по гаусам согласовался с записью по картике bmp
11. Условия проходимости: 1) Если угол наклона в любом пикселе путя не превосходит допустимый (по направлению или вбок) 2) Расстояние до препятсвия на срезе меньше радиуса

## 📜 Командный файл (примеры)
1) Если нужно прочитать данные с файла Read.bmp
```
init
help
bmp_read results/visualizations/Read.bmp
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Polew.bmp Full
bin 147 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 5
PlotKmeans results/visualizations/kmeans.png
triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph
astar_graph
save_metrics
PlotPath results/visualizations/Path.png
end
```
2) Если данные вводятся с помощью гаусов
```
init
help
g 99 50 25 25 -25
g 50 50 20 20 25
g 200 50 20 20 -25
g 50 200 20 20 25
g 189 180 20 20 -25
generate
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full
bin 147 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 10
PlotKmeans results/visualizations/kmeans.png
triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph 60 130 150 135
astar_graph
save_metrics
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
```

## 📃️ Конфигурационный файл (пример)
```
defaultHelp var/help/help.txt


defaultfieldWidth 300
defaultfieldHeight 300


defaultCenterX 50.0
defaultCenterY 50.0
defaultSigmaX 5.0
defaultSigmaY 5.0
defaultHeight 10.0


count_min 50
count_max 100
xmin 0
xmax 300
ymin 0
ymax 300
sx_min 2
sx_max 5
sy_min 2
sy_max 5
h_min -120
h_max 120


defaultGnuplot results/visualizations/Gnuplot.png
defaultPlotMetedata results/visualizations/Metadata.png
defaultPlotKmeans results/visualizations/kmeans.png
defaultPlotGraph results/visualizations/Graph.png
defaultPlotVoronoi results/visualizations/Voronoi.png
defaultPlotDelaunay results/visualizations/Delaunay.png
defaultPlotGrid results/visualizations/Grid.png
defaultPlotNavGrid results/visualizations/NavGrid.png
defaultPlotGridPath results/visualizations/GridPath.png
defaultPlotPath results/visualizations/Path.png
PlotRRT results/visualizations/RRT.gif
defaultPlot3DPath results/visualizations/Plot3DPath.png

defaultWrite results/visualizations/Write.bmp 
defaultWriteModeImage Full
defaultRead results/visualizations/Read.bmp

save_g config/commands/gaussians.txt


heightThresholdPixel 3


defaultWaveNoisy 10


defaultKlaster 5
defaultKlasterKern 5


defaultgridWidth 1
defaultgridNoisy 0


startPixelX 150
startPixelY 150
goalPixelX 160
goalPixelY 160
defaultConnectMode All

vehicleRadiusPixel 1
maxSideAngle 90.0
maxUpDownAngle 90.0


maxIterations 10000
startWorldX 1.0
startWorldY 1.0
goalWorldX 100.0
goalWorldY 100.0
vehicleRadiusWorld 1.0
heightThresholdWorld 5.0
interpEdge 1.0
interpCollision 1.0
interpAngle 1.0
step 1.0
goalRadius 2.0
goalBias 0.2


defaultsave_metrics var/metrics/metrics.csv


logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt

FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO


seedMode Fixed 42
```

Полная документация: см. README.md
)";
    // 1. Вывод в консоль
    std::cout << helpContent << std::endl;

    // 2. Сохранение в файл
    std::ofstream helpFile(config.defaultHelp);
    
    if (helpFile) {
        helpFile << helpContent;
        logger.info("Help file created successfully at" + config.defaultHelp);
    } else {
        logger.error("Failed to create help file at " + config.defaultHelp);
        std::cerr << "Error: Could not create help file at " + config.defaultHelp << std::endl;
    }
}

    // Обработка команд из файла
    void Interface::processFileCommands(std::ifstream& file) {
        int commandCount = 0;
        
        while (file >> params.command) {
            commandCount++;
            logger.info(std::string("Processing command #") + std::to_string(commandCount) + ": " + params.command);
            
            if (params.command == "help") {
                showHelp();
                continue;
            }
            
            if (params.command == "end") {
                logger.info("Ending the program.");
                break;
            }
            
            if (!processCommand(file, false)) {
                break;
            }
        }
    }

    // Обработка команд с клавиатуры
    void Interface::processKeyboardCommands() {
        
        while (true) {
            const std::string commandshow = R"(Enter the command and its parameters immediately (help, init, g, g_auto, generate, save_g, gnuplot, PlotKmeans, PlotMetedata, PlotVoronoi, PlotDelaunay,
PlotGrid, PlotNavGrid, PlotGridPath, PlotPath, PlotRRT, bmp_write, bmp_read, bin, wave, k_means, k_means_kern, triangulate, voronoi, build_nav_graph, grid, build_nav_grid, connect_to_grid, 
connect_to_graph, astar_graph, dekstra_graph, greedy_graph, astar_grid, dekstra_grid, greedy_grid, save_metrics, Plot3DPath, plotInteractive3DPath, end, rrt):)";
            std::cout << commandshow;
            std::cin >> params.command;
            std::cout << "\n";
            logger.info(std::string("Received command: ") + params.command);
            
            if (params.command == "help") {
                showHelp();
                continue;
            }
            
            if (params.command == "end") {
                std::cout << "Ending the program" << std::endl;
                logger.info("Ending the program.");
                break;
            }
            
            if (!processCommand(std::cin, true)) {
                break;
            }
        }
    }

    // Общая обработка команды (для файла и клавиатуры)
    bool Interface::processCommand(std::istream& input, bool fromKeyboard) {
    std::string line;
    std::string showInfo;
    std::string modeWrite, modeConnect;
    
    if (params.command == "init") {
        if (initState) {
            throw std::runtime_error("Error: Multiple init commands.");
            return false;
        }
        
        params.fieldWidth = config.defaultfieldWidth;
        params.fieldHeight = config.defaultfieldHeight;
        
        std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.fieldWidth >> params.fieldHeight;
            
        Validator::validateFieldSize(params.fieldWidth, params.fieldHeight);
    
        initState = true;
            
        showInfo = std::string("Initializing field with size: ") + std::to_string(params.fieldWidth) + " x " + std::to_string(params.fieldHeight);
        
        if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Field initialized.");
    }
    else if (!initState) {
        throw std::runtime_error("The init command was not used.");
        return false;
    }
    else if (params.command == "g") {
        params.centerX = config.defaultCenterX;
        params.centerY = config.defaultCenterY;
        params.sigmaX = config.defaultSigmaX;
        params.sigmaY = config.defaultSigmaY;
        params.height = config.defaultHeight;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.centerX >> params.centerY >> params.sigmaX >> params.sigmaY >> params.height;
            
        Validator::validateGaussian(params);
         
        showInfo = std::string("Adding Gaussian: x=") + std::to_string(params.centerX) + 
                 ", y=" + std::to_string(params.centerY) + 
                 ", sx=" + std::to_string(params.sigmaX) + 
                 ", sy=" + std::to_string(params.sigmaY) + 
                 ", h=" + std::to_string(params.height);
                  
         if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info(showInfo);
        control.Dispetcher(params);
    }
    else if (params.command == "g_auto") {
        params.count_min = config.count_min;
        params.count_max = config.count_max;
        params.xmin = config.xmin;
        params.xmax = config.xmax;
        params.ymin = config.ymin;
        params.ymax = config.ymax;
        params.sx_min = config.sx_min;
        params.sx_max = config.sx_max;
        params.sy_min = config.sy_min;
        params.sy_max = config.sy_max;
        params.h_min = config.h_min;
        params.h_max = config.h_max;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.count_min >> params.count_max
            >> params.xmin >> params.xmax
            >> params.ymin >> params.ymax 
            >> params.sx_min >> params.sx_max
            >> params.sy_min >> params.sy_max
            >> params.h_min >> params.h_max;
                
      Validator::validateAutoGaussian(params);
        
         showInfo = std::string("Автогенерация гауссов: ") +
               "x=[" + std::to_string(params.xmin) + ", " + std::to_string(params.xmax) + "], " +
               "y=[" + std::to_string(params.ymin) + ", " + std::to_string(params.ymax) + "], " +
               "sx=[" + std::to_string(params.sx_min) + ", " + std::to_string(params.sx_max) + "], " +
               "sy=[" + std::to_string(params.sy_min) + ", " + std::to_string(params.sy_max) + "], " +
               "h=[" + std::to_string(params.h_min) + ", " + std::to_string(params.h_max) + "], " +
               "count=[" + std::to_string(params.count_min) + ", " + std::to_string(params.count_max) + "]";
                   
         if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info(showInfo);
        control.Dispetcher(params);
    }
    else if (params.command == "generate") {
        showInfo = "Generating field by summing all Gaussians";
        logger.info(showInfo);
        control.Dispetcher(params);
        if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info("Field generation completed");
    }
    else if (params.command == "save_g") {
        params.filename = config.save_g;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;

            showInfo = std::string("Save save_g to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("save_g save completed");
    }
    else if (params.command == "gnuplot") {
        params.filename = config.defaultGnuplot;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
       Validator::validateFileName(params.filename);
           
            showInfo = std::string("Calling gnuplot with filename: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Gnuplot visualization completed");
    }
    else if (params.command == "PlotMetedata") {
        params.filename = config.defaultPlotMetedata;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
            Validator::validateFileName(params.filename);
                
            showInfo = std::string("Plotting metadata to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Metadata plotting completed");
    }
    else if (params.command == "PlotKmeans") {
        params.filename = config.defaultPlotKmeans;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
            Validator::validateFileName(params.filename);
 
            showInfo = std::string("Plotting kmeans to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Kmeans plotting completed");
    }
    else if (params.command == "PlotGraph") {
        params.filename = config.defaultPlotGraph;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                        
            Validator::validateFileName(params.filename);
 
            showInfo = std::string("Plotting graph to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Graph plotting completed");
    }
    else if (params.command == "PlotVoronoi") {
        params.filename = config.defaultPlotVoronoi;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                        
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting Voronoi diagram to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }        

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Voronoi diagram plotting completed");
    }
    else if (params.command == "PlotDelaunay") {
        params.filename = config.defaultPlotDelaunay;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting Delaunay triangulation to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation plotting completed");
    }
    else if (params.command == "PlotGrid") {
        params.filename = config.defaultPlotGrid;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting grid to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Grid plotting completed");
    }
    else if (params.command == "PlotNavGrid") {
        params.filename = config.defaultPlotNavGrid;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting nav grid to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Nav grid plotting completed");
    }
    else if (params.command == "PlotGridPath") {
        params.filename = config.defaultPlotGridPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting grid path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Grid path plotting completed");
    }
    else if (params.command == "PlotPath") {
        params.filename = config.defaultPlotPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Path plotting completed");
    }
    else if (params.command == "PlotRRT") {
        params.filename = config.PlotRRT;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Plotting rrt to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("RRT plotting completed");
    }
    else if (params.command == "bmp_write") {
        params.filename = config.defaultWrite;
        modeWrite = config.defaultWriteModeImage;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename >> modeWrite;
            
            if (modeWrite == "Full") {
            params.bmpWriteMode = io::BmpWriteMode::Full;
        } 
       else if (modeWrite == "Binary") {
            params.bmpWriteMode = io::BmpWriteMode::Binary;
        } 
        
            Validator::validateBmpWriteMode(params.filename, modeWrite);
                    
        showInfo = std::string("Writing BMP file: ") + params.filename + " with mode: " + modeWrite;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("BMP writing completed");
    }
    else if (params.command == "bmp_read") {
        params.filename = config.defaultRead;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                        
            Validator::validateFileName(params.filename);
 
        showInfo = std::string("Reading BMP file: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("BMP reading completed");
    }
    else if (params.command == "bin") {
        params.heightThresholdPixel = config.heightThresholdPixel;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.heightThresholdPixel;
                        
            Validator::validateBinaryParameters(params.heightThresholdPixel);
                  
        showInfo = std::string("Applying binary filter with slice: ") + std::to_string(params.heightThresholdPixel);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Binary filtering completed");
    }
    else if (params.command == "wave") {
        params.waveNoisy = config.defaultWaveNoisy;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.waveNoisy;
                                
            Validator::validateWaveNoiseSize(params.waveNoisy);
           
        showInfo = std::string("Applying wave filter with noisy level: ") + std::to_string(params.waveNoisy);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info(std::string("Wave filtering completed"));
    }
    else if (params.command == "k_means") {
        params.clusterCount = config.defaultKlaster;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.clusterCount;
                                         
            Validator::validateKMeans(params.clusterCount);
          
        showInfo = std::string("Running k-means with k: ") + std::to_string(params.clusterCount);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("k-means clustering completed");
    }
    else if (params.command == "k_means_kern") {
        params.clusterCount = config.defaultKlaster;
        params.kernelSize = config.defaultKlasterKern;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.clusterCount >> params.kernelSize;
                                        
            Validator::validateKernelSize(params.clusterCount, params.kernelSize);
           
        showInfo = std::string("Running k-means with k:") + std::to_string(params.clusterCount) + 
                   std::string(", kernel size: ") + std::to_string(params.kernelSize);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("k-means with kernels completed");
    }
    else if (params.command == "triangulate") {
        showInfo = "Starting Delaunay triangulation";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation completed");
    }
    else if (params.command == "voronoi") {
        showInfo = "Building Voronoi diagram";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Voronoi diagram construction completed");
    }
    else if (params.command == "build_nav_graph") {
        params.vehicleRadiusPixel = config.vehicleRadiusPixel;
        params.maxSideAngle = config.maxSideAngle;
        params.maxUpDownAngle = config.maxUpDownAngle;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.vehicleRadiusPixel >> params.maxSideAngle >> params.maxUpDownAngle;
        
                                        
            Validator::validateNavigationParameters(params.vehicleRadiusPixel, params.maxSideAngle, params.maxUpDownAngle);
                  
        showInfo = std::string("Building navigation graph with vehicleRadiusPixel = ") + std::to_string(params.vehicleRadiusPixel) + 
                                std::string(", maxSideAngle = ") + std::to_string(params.maxSideAngle) +
                                std::string(", maxUpDownAngle = ") + std::to_string(params.maxUpDownAngle);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Navigation graph construction completed");
    }
    else if (params.command == "grid") {
        params.gridWidth = config.defaultgridWidth;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.gridWidth;
                                         
            Validator::validateGrid(params.fieldWidth, params.fieldHeight, params.gridWidth);
          
        showInfo = std::string("Building grid with cell widht = ") + std::to_string(params.gridWidth);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Grid completed");
    }
    else if (params.command == "build_nav_grid") {
        params.gridNoisy = config.defaultgridNoisy;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.gridNoisy;
                                
            Validator::validateGridNoiseSize(params.gridNoisy);
           
        showInfo = std::string("Applying wave filter with noisy level: ") + std::to_string(params.gridNoisy);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Navigation grid construction completed");
    }
    else if (params.command == "connect_to_grid") {
        params.startPixelX = config.startPixelX;
        params.startPixelY = config.startPixelY;
        params.goalPixelX = config.goalPixelX;
        params.goalPixelY = config.goalPixelY;        
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.startPixelX >> params.startPixelY >> params.goalPixelX >> params.goalPixelY;
                
            Validator::validateConnectParameters(params);
 
            showInfo =
                std::string("Connecting points (") +
                std::to_string(params.startPixelX) + ", " +
                std::to_string(params.startPixelY) + ") and (" +
                std::to_string(params.goalPixelX) + ", " +
                std::to_string(params.goalPixelY) + ") to navigation grid";

            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Connection completed");
    }
    else if (params.command == "connect_to_graph") {
        params.startPixelX = config.startPixelX;
        params.startPixelY = config.startPixelY;
        params.goalPixelX = config.goalPixelX;
        params.goalPixelY = config.goalPixelY;       
        modeConnect = config.defaultConnectMode;
        params.nearestVerticesCount = config.defaultNearestVerticesCount;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.startPixelX >> params.startPixelY >> params.goalPixelX >> params.goalPixelY >> modeConnect;
            
            if (modeConnect == "Nearest")
            params.connectMode =
                algorithms::path::common::ConnectMode::Nearest;

        else if (modeConnect == "All")
            params.connectMode =
                algorithms::path::common::ConnectMode::All;

        else if (modeConnect == "NearestK") {
            params.connectMode =
                algorithms::path::common::ConnectMode::NearestK;
            iss >> params.nearestVerticesCount;
        }
                
            Validator::validateConnectParameters(params, modeConnect);
 
            showInfo =
                std::string("Connecting points (") +
                std::to_string(params.startPixelX) + ", " +
                std::to_string(params.startPixelY) + ") and (" +
                std::to_string(params.goalPixelX) + ", " +
                std::to_string(params.goalPixelY) + ") to navigation graph" + ", mode [" + modeConnect + "]";
    
            if (modeConnect == "NearestK")
            showInfo +=", k = " + std::to_string(params.nearestVerticesCount);

            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Connection completed");
    }
    else if (params.command == "astar_graph") {
        showInfo = "Starting A* graph path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("A* graph path search completed");
    }
    else if (params.command == "dekstra_graph") {
        showInfo = "Starting dekstra graph path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Dekstra graph path search completed");
    }
    else if (params.command == "greedy_graph") {
        showInfo = "Starting greedy graph path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Greedy graph path search completed");
    }
    else if (params.command == "astar_grid") {
        showInfo = "Starting A* grid path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("A* grid path search completed");
    }
    else if (params.command == "dekstra_grid") {
        showInfo = "Starting dekstra grid path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Dekstra grid path search completed");
    }
    else if (params.command == "greedy_grid") {
        showInfo = "Starting greedy grid path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Greedy grid path search completed");
    }
    else if (params.command == "rrt") {
        params.maxIterations = config.maxIterations;
        params.startWorldX = config.startWorldX;
        params.startWorldY = config.startWorldY;
        params.goalWorldX = config.goalWorldX;
        params.goalWorldY = config.goalWorldY;
        params.vehicleRadiusWorld = config.vehicleRadiusWorld;
        params.heightThresholdWorld = config.heightThresholdWorld;
        params.maxSideAngle = config.maxSideAngle;
        params.maxUpDownAngle = config.maxUpDownAngle;
        params.interpEdge = config.interpEdge;
        params.interpCollision = config.interpCollision;
        params.interpAngle = config.interpAngle;
        params.step = config.step;
        params.goalRadius = config.goalRadius;
        params.goalBias = config.goalBias;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.maxIterations
            >> params.startWorldX >> params.startWorldY
            >> params.goalWorldX >> params.goalWorldY 
            >> params.vehicleRadiusWorld >> params.heightThresholdWorld
            >> params.maxSideAngle >> params.maxUpDownAngle
            >> params.interpEdge >> params.interpCollision >> params.interpAngle
            >> params.step
            >> params.goalRadius
            >> params.goalBias;
                
      Validator::validateRRT(params);
        
         showInfo = std::string("RRT: ") +
             "maxIterations=" + std::to_string(params.maxIterations) +
             ", start=(" + std::to_string(params.startWorldX) + ", " +
                   std::to_string(params.startWorldY) + ")" +
             ", goal=(" + std::to_string(params.goalWorldX) + ", " +
                  std::to_string(params.goalWorldY) + ")" +
             ", vehicleRadius=" + std::to_string(params.vehicleRadiusWorld) +
             ", heightThreshold=" + std::to_string(params.heightThresholdWorld) +
             ", maxSideAngle=" + std::to_string(params.maxSideAngle) +
             ", maxUpDownAngle=" + std::to_string(params.maxUpDownAngle) +
             ", interpEdge=" + std::to_string(params.interpEdge) +
             ", interpCollision=" + std::to_string(params.interpCollision) +
             ", interpAngle=" + std::to_string(params.interpAngle) +
             ", step=" + std::to_string(params.step) +
             ", goalRadius=" + std::to_string(params.goalRadius) +
             ", goalBias=" + std::to_string(params.goalBias);
                   
         if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info(showInfo);
        control.Dispetcher(params);
    }
    else if (params.command == "save_metrics") {
        params.filename = config.defaultsave_metrics;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                               
            Validator::validateFileName(params.filename);
     
        showInfo = std::string("Metrics saving to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("save_metrics completed");
    }
    else if (params.command == "Plot3DPath") {
        params.filename = config.defaultPlot3DPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                               
            Validator::validateFileName(params.filename);
     
        showInfo = std::string("Plotting 3D path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("3D path plotting completed");
    }
    else if (params.command == "plotInteractive3DPath") {
        showInfo = "Starting interactive 3D path visualization";
        
        if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Interactive 3D visualization completed");
    }
    else {
        throw std::runtime_error(
            "Unknown command: '" + params.command + "'"
        );
    }
    
    return true;
}

    Interface::Interface(core::Config& cfg, core::Logger& log, Control& c) 
        : config(cfg), logger(log), control(c) {
        logger.info("Interface initialized");
    }
    
    void Interface::print() {
        logger.info("Starting user interface");
        
        bool useFileInput;
        std::string commandfilename;
        
        std::cout << "Hello, dear user, this program builds Gaussians.\n"
                  << "Enter commands from a text file (PRESS 1) or from the keyboard (PRESS 0)?" 
                  << std::endl;
        std::cin >> useFileInput;
        
        logger.info(std::string("User chose input method: ") + 
                  std::string(useFileInput ? "file" : "keyboard"));
        
        if (useFileInput) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> commandfilename;
            
            logger.info(std::string("Attempting to open command file: ") + commandfilename);
            std::ifstream file(commandfilename);
            
            if (!file) {
                logger.error(std::string("Failed to open command file: ") + commandfilename);
                std::cout << "File not found" << std::endl;
                return;
            }
            
            processFileCommands(file);
            file.close();
            logger.info("Closed input file.");
        } else {
            logger.info("User selected keyboard input mode");
            std::cout << "You will enter commands from the keyboard" << std::endl;
            processKeyboardCommands();
        }
    }
}
