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

| Команда              | Параметры                                                                                          | Описание                                                                 |
|----------------------|----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| help                 | -                                                                                                  | Создание файла с пояснением команд                                       |
| init                 | fieldWidth fieldHeight                                                                             | Инициализация поля с заданой шириной и длиной                            |
| g                    | x y sx sy h                                                                                        | Создает гаусс с центром (x,y), размерами (sx,sy) и высотой h             |
| g_auto               | xmin, xmax, ymin, ymax, ... , h_min, h_max, [Random/Fixed seed]                                    | Создает случайные гаусы с параметрами в промежутках + режим генерации    |
| generate             | -                                                                                                  | Складывает все добавленные гауссы в итоговое поле                        |
| save_g               | filename.png                                                                                       | Сохраняет параметры g в txt файл                                         |
| gnuplot              | filename.png                                                                                       | Сохраняет 3D-визуализацию поля в PNG файл                                |
| PlotMetedata         | filename.png                                                                                       | Визуализирует метаданные компонент с границами и центрами                |
| PlotKmeans           | filename.png                                                                                       | Визуализирует k_means                                                    |
| PlotGraph            | filename.png                                                                                       | Визуализирует граф                                                       |
| PlotVoronoi          | filename.png                                                                                       | Строит диаграмму Вороного по текущей триангуляции                        |
| PlotDelaunay         | filename.png                                                                                       | Визуализирует триангуляцию Делоне                                        |
| PlotPath             | filename.png                                                                                       | Отображает найденный путь между точками A и B                            |
| bmp_write            | filename.bmp [Full/Binary]                                                                         | Сохраняет поле в BMP: Full - полное, Binary - бинаризованное             |
| bmp_read             | filename.bmp                                                                                       | Загружает поле из BMP файла                                              |
| bin                  | slice [Peaks/Valleys/All]                                                                          | Бинаризация: Peaks - только пики, Valleys - впадины, All - по модулю     |
| wave                 | noisy                                                                                              | Удаляет компоненты размером ≤ noisy как шум                              |
| k_means              | k                                                                                                  | Кластеризует данные в k кластеров                                        |
| k_means_kern         | kk                                                                                                 | Кластеризация с ядрами размера kk                                        |
| triangulate          | -                                                                                                  | Строит триангуляцию Делоне по центрам компонент                          |
| voronoi              | -                                                                                                  | Строит диаграмму Вороного                                                |
| build_nav_graph      | vehicleRadius maxSideAngle maxUpDownAngle                                                          | Строит граф проходимости по диаграмме Вороного, учитывая наклон и радиус |
| connect_to_graph     | Ax Ay Bx By [Nearest/NearestK k/All]                                                               | Подключает старт и финиш к графу с разными режимами                      |
| find_path_astar      | -                                                                                                  | A* ищет путь между точками A и B                                         |
| find_path_dekstra    | -                                                                                                  | Dekstra ищет путь между точками A и B                                    |
| find_path_greedy     | -                                                                                                  | Greedy ищет путь между точками A и B                                     |
| save_metrics         | filename.csv                                                                                       | Сохраняет метрики в файле csv                                            |
| Plot3DPath           | filename.png                                                                                       | Сохраняет 3D-визуализацию путя в PNG файл                                |
| plotInteractive3DPath| -                                                                                                  | Интерактвный 3D режим с путем                                            |
| end                  | -                                                                                                  | Завершает работу программы                                               |


## ⚙️ Параметры конфигурационного файла (config.txt)

| Параметр                     | Значение                                       | Описание                                                                         |
|------------------------------|------------------------------------------------|----------------------------------------------------------------------------------|
| defaultfieldWidth            | `fieldWidth`                                   | Ширина рабочего поля в пикселях по умолчанию                                     |  
| defaultfieldHeight           | `fieldHeight`                                  | Высота рабочего поля в пикселях по умолчанию                                     |
| defaultCenterX               | `defaultCenterX`                               | Стандартная X-координата центра гауссова распределения по умолчанию              |
| defaultCenterY               | `defaultCenterY`                               | Стандартная Y-координата центра гауссова распределения по умолчанию              |
| defaultSigmaX                | `defaultSigmaX`                                | Стандартное отклонение по оси X по умолчанию                                     |
| defaultSigmaY                | `defaultSigmaY`                                | Стандартное отклонение по оси Y по умолчанию                                     |
| defaultHeight                | `defaultHeight`                                | Стандартная высота гауссова распределения по умолчанию                           |
| xmin                         | `xmin`                                         | Минимальная X-координата центра гаусса при g_auto                                |
| xmax                         | `xmax`                                         | Максимальная X-координата центра гаусса при g_auto                               |
| ymin                         | `ymin`                                         | Минимальная Y-координата центра гаусса при g_auto                                |
| ymax                         | `ymax`                                         | Максимальная Y-координата центра гаусса при g_auto                               |
| sx_min                       | `sx_min`                                       | Минимальное отклонение по оси X при g_auto                                       |
| sx_max                       | `sx_max`                                       | Максимальное отклонение по оси X при g_auto                                      |
| sy_min                       | `sy_min`                                       | Минимальное отклонение по оси Y при g_auto                                       |
| sy_max                       | `sy_max`                                       | Максимальное отклонение по оси Y при g_auto                                      |
| h_min                        | `h_min`                                        | Минимальная высота гаусса при g_auto                                             |
| h_max                        | `h_max`                                        | Максимальная высота гаусса при g_auto                                            |
| count_min                    | `count_min`                                    | Минимальное количество генерируемых гауссов                                      |
| count_max                    | `count_max`                                    | Максимальное количество генерируемых гауссов                                     |
| defaultGAutoMode             | 'GAutoMode'                                    | Режим генерации (Random/Fixed seed) по умолчанию                                 |
| save_g                       | `filename_save_g`                              | Путь к файлу для сохранения параметров g                                         |
| defaultGnuplot               | `filename_gnuplot.png`                         | Путь к файлу для сохранения 3D-визуализации по умолчанию                         |
| defaultPlotMetedata          | `filename_metadata.png`                        | Путь к файлу для визуализации метаданных компонент по умолчанию                  |
| defaultPlotKmeans            | `filename_kmeans.png`                          | Путь к файлу для k_means по умолчанию                                            |
| defaultPlotGraph             | `filename_graph.png`                           | Путь к файлу для графа по умолчанию                                              |
| defaultPlotVoronoi           | `filename_voronoi.png`                         | Путь к файлу для диаграммы Вороного по умолчанию                                 |
| defaultPlotDelaunay          | `filename_delaunay.png`                        | Путь к файлу для триангуляции Делоне по умолчанию                                |
| defaultPlotPath              | `filename_path.png`                            | Путь к файлу для визуализации маршрута по умолчанию                              |
| defaultWrite                 | `filename_write.bmp`                           | Путь к файлу для сохранения BMP-изображения по умолчанию                         |
| defaultWriteModeImage        | `writeMode`                                    | Режим сохранения BMP (Full/Binary) по умолчанию                                  |
| defaultRead                  | `filename_read.bmp`                            | Путь к файлу для загрузки BMP-изображения по умолчанию                           |
| defaultThreshold             | `defaultThreshold`                             | Порог бинаризации по умолчанию                                                   |
| defaultBinMode               | `binMode`                                      | Режим бинаризации (Peaks/Valleys/All) по умолчанию                               |
| defaultNoisy                 | `defaultNoisy`                                 | Порог для удаления шумовых компонент по умолчанию                                |
| defaultKlaster               | `defaultKlaster`                               | Количество кластеров для k-mean по умолчанию                                     |
| defaultKlasterKern           | `defaultKlasterKern`                           | Размер ядра для кластеризации по умолчанию                                       |
| defaultpointA_x              | `pointA_x`                                     | X-координата точки A для поиска пути по умолчанию                                |
| defaultpointA_y              | `pointA_y`                                     | Y-координата точки A для поиска пути по умолчанию                                |
| defaultpointB_x              | `pointB_x`                                     | X-координата точки B для поиска пути по умолчанию                                |
| defaultpointB_y              | `pointB_y`                                     | Y-координата точки B для поиска пути по умолчанию                                |
| defaultConnectMode           | `connectMode`                                  | Режимы подключения старта и финиша к графу (Nearest/NearestK k/All) по умолчанию |
| defaultPlot3DPath            | `filename_plot3dpath.png`                      | Путь к файлу для 3D-визуализации маршрута по умолчанию                           |
| vehicleRadius                | `vehicleRadius`                                | Радиус транспортного средства                                                    |
| maxSideAngle                 | `maxSideAngle`                                 | Максимальный угол поворота вбок (градусы)                                        |
| maxUpDownAngle               | `maxUpDownAngle`                               | Максимальный угол наклона вверх/вниз (градусы)                                   |
| defaultsave_metrics          | `filename.csv`                                 | Путь к файлу по умолчанию для сохранения метрик                                  |
| logFileNameInterface         | `filename_log_interface.txt`                   | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | `filename_log_control.txt`                     | Путь к лог-файлу управления                                                      |
| defaultHelp                  | `/home/log/Gauss/results/docs/help.txt`        | Путь где сохранить help файл                                                     |
| FiltrationLogLevelInterface  | `logLevelInterface`                            | Уровень логирования интерфейса (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |
| FiltrationLogLevelControl    | `logLevelControl`                              | Уровень логирования управления (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |


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
find_path_astar
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
find_path_astar
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
count_min 50
count_max 100
defaultGAutoMode Fixed 42


defaultGnuplot results/visualizations/Gnuplot.png
defaultPlotMetedata results/visualizations/Metadata.png
defaultPlotKmeans results/visualizations/kmeans.png
defaultPlotGraph results/visualizations/Graph.png
defaultPlotVoronoi results/visualizations/Voronoi.png
defaultPlotDelaunay results/visualizations/Delaunay.png
defaultPlotPath results/visualizations/Path.png
defaultPlot3DPath results/visualizations/Plot3DPath.png

defaultWrite results/visualizations/Write.bmp 
defaultWriteModeImage Full
defaultRead results/visualizations/Read.bmp

save_g config/commands/gaussians.txt


defaultThreshold 130
defaultBinMode All


defaultNoisy 10


defaultKlaster 5
defaultKlasterKern 5


defaultstartPointX 150
defaultstartPointY 150
defaultendPointX 160
defaultendPointY 160
defaultConnectMode All

vehicleRadius 1
maxSideAngle 90.0
maxUpDownAngle 90.0


defaultsave_metrics var/metrics/metrics.csv


logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt

FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO
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
PlotPath, bmp_write, bmp_read, bin, wave, k_means, k_means_kern, triangulate, voronoi, build_nav_graph, connect_to_graph, find_path_astar, find_path_dekstra, find_path_greedy, save_metrics, 
Plot3DPath, plotInteractive3DPath, end):)";
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
    std::string modeGAuto, modeWrite, modeBin, modeConnect;
    
    if (params.command == "init") {
        if (initState) {
            std::cout << "The init command has already been called.\nError\n";
            logger.error("Error: Multiple init commands.");
            return false;
        }
        
        params.fieldWidth = config.defaultfieldWidth;
        params.fieldHeight = config.defaultfieldHeight;
        
        std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.fieldWidth >> params.fieldHeight;
            
        if (!Validator::validateFieldSize(params.fieldWidth, params.fieldHeight, logger))
            return false;
    
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
        std::cout << "The init command was not used.\nError\n";
        logger.error("Error: The init command was not used.");
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
            
        if (!Validator::validateGaussian(params, logger))
            return false;
         
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
        params.count_min = config.count_min;
        params.count_max = config.count_max;
        modeGAuto = config.defaultGAutoMode;
        params.seedGAuto = config.defaultSeedGAuto;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.xmin >> params.xmax >> params.ymin >> params.ymax 
            >> params.sx_min >> params.sx_max >> params.sy_min >> params.sy_max
            >> params.h_min >> params.h_max
            >> params.count_min >> params.count_max
            >> modeGAuto;
      
      if (modeGAuto == "Random") {
            params.gAutoMode = algorithms::gauss::GAutoMode::Random;
        } 
 else if (modeGAuto == "Fixed") {
            params.gAutoMode = algorithms::gauss::GAutoMode::Fixed;
            iss >> params.seedGAuto;
        }
                
      if (!Validator::validateAutoGaussian(params, modeGAuto, logger))
          return false;
        
         showInfo = std::string("Автогенерация гауссов: ") +
               "x=[" + std::to_string(params.xmin) + ", " + std::to_string(params.xmax) + "], " +
               "y=[" + std::to_string(params.ymin) + ", " + std::to_string(params.ymax) + "], " +
               "sx=[" + std::to_string(params.sx_min) + ", " + std::to_string(params.sx_max) + "], " +
               "sy=[" + std::to_string(params.sy_min) + ", " + std::to_string(params.sy_max) + "], " +
               "h=[" + std::to_string(params.h_min) + ", " + std::to_string(params.h_max) + "], " +
               "count=[" + std::to_string(params.count_min) + ", " + std::to_string(params.count_max) + "], " +
               "gAutoMode=" + modeGAuto;
               
         if (modeGAuto == "Fixed")
          showInfo += ", seedGAuto=" + std::to_string(params.seedGAuto);
                   
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
            
       if (!Validator::validateFileName(params.filename, logger))
           return false;
           
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
            
            if (!Validator::validateFileName(params.filename, logger))
                return false;
                
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
            
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
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
                        
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
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
                        
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
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
                         
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
        showInfo = std::string("Plotting Delaunay triangulation to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation plotting completed");
    }
    else if (params.command == "PlotPath") {
        params.filename = config.defaultPlotPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                         
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
        showInfo = std::string("Plotting path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Path plotting completed");
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
        
            if (!Validator::validateBmpWriteMode(params.filename, modeWrite, logger))
                return false;
                    
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
                        
            if (!Validator::validateFileName(params.filename, logger))
                return false;
 
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
        params.threshold = config.defaultThreshold;
        modeBin = config.defaultBinMode;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.threshold >> modeBin;
                        
            if (!Validator::validateBinaryParameters(params.threshold, modeBin, logger))
                return false;
 
        if (modeBin == "Peaks") {
            params.thresholdMode = algorithms::components::ThresholdMode::Peaks;
        } else if (modeBin == "Valleys") {
            params.thresholdMode = algorithms::components::ThresholdMode::Valleys;
        } else {
            params.thresholdMode = algorithms::components::ThresholdMode::All;
        }
                  
        showInfo = std::string("Applying binary filter with slice: ") + std::to_string(params.threshold) + " and mode: " + modeBin;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Binary filtering completed");
    }
    else if (params.command == "wave") {
        params.noiseLevel = config.defaultNoisy;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.noiseLevel;
                                
            if (!Validator::validateNoiseSize(params.noiseLevel, logger))
                return false;
           
        showInfo = std::string("Applying wave filter with noisy level: ") + std::to_string(params.noiseLevel);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info(std::string("Wave filtering completed. Components count: ") + 
                   std::to_string(state.components.size()));
    }
    else if (params.command == "k_means") {
        params.clusterCount = config.defaultKlaster;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.clusterCount;
                                         
            if (!Validator::validateKMeans(params.clusterCount, logger))
                return false;
          
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
                                        
            if (!Validator::validateKernelSize(params.clusterCount, params.kernelSize, logger))
                return false;
           
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
        params.vehicleRadius = config.vehicleRadius;
        params.maxSideAngle = config.maxSideAngle;
        params.maxUpDownAngle = config.maxUpDownAngle;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.vehicleRadius >> params.maxSideAngle >> params.maxUpDownAngle;
        
                                        
            if (!Validator::validateNavigationParameters(params.vehicleRadius, params.maxSideAngle, params.maxUpDownAngle, logger))
                return false;
                  
        showInfo = std::string("Building navigation graph with vehicleRadius = ") + std::to_string(params.vehicleRadius) + 
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
    else if (params.command == "connect_to_graph") {
        params.startPointX = config.defaultstartPointX;
        params.startPointY = config.defaultstartPointY;
        params.endPointX = config.defaultendPointX;
        params.endPointY = config.defaultendPointY;        
        modeConnect = config.defaultConnectMode;
        params.nearestVerticesCount = config.defaultNearestVerticesCount;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.startPointX >> params.startPointY >> params.endPointX >> params.endPointY >> modeConnect;
            
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
                
            if (!Validator::validateConnectParameters(params, modeConnect, logger))
                return false;
 
            showInfo =
                std::string("Connecting points (") +
                std::to_string(params.startPointX) + ", " +
                std::to_string(params.startPointY) + ") and (" +
                std::to_string(params.endPointX) + ", " +
                std::to_string(params.endPointY) + ") to navigation graph" + ", mode [" + modeConnect + "]";
    
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
    else if (params.command == "find_path_astar") {
        showInfo = "Starting A* path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("A* path search completed");
    }
    else if (params.command == "find_path_dekstra") {
        showInfo = "Starting dekstra path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Dekstra path search completed");
    }
    else if (params.command == "find_path_greedy") {
        showInfo = "Starting greedy path search";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Greedy path search completed");
    }
    else if (params.command == "save_metrics") {
        params.filename = config.defaultsave_metrics;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
                               
            if (!Validator::validateFileName(params.filename, logger))
                return false;
     
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
                               
            if (!Validator::validateFileName(params.filename, logger))
                return false;
     
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
        std::cout << "Unknown command: " << params.command << std::endl;
        logger.warning(std::string("Unknown command received: ") + params.command);
        return false;
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
