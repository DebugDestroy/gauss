#pragma once
#include <vector>       // Для std::vector
#include <string>       // Для std::string
#include <fstream>      // Для std::ifstream, std::ofstream
#include <sstream>      // Для std::istringstream
#include <iostream>     // Для std::cin, std::cout
#include <unordered_map> // Для std::unordered_map (если потребуется)

// Локальные заголовки
#include "core/Config.hpp"          // Для Config
#include "core/Logger.hpp"          // Для Logger
#include "command/DispatcherParams.hpp"  // Новый include
#include "io/BmpHandler.hpp"        // Для BmpWriteMode
#include "algorithms/Component.hpp" // Для Component
#include "core/Geometry.hpp"        // Для PointD
#include "algorithms/PathFinder.hpp" // Для PathFinder

// Предварительное объявление класса Control
class Control;


class Interface {
private:
    Config& config;
    Logger& logger;
    Control& control;
    DispatcherParams params;
    int n = 0; // Флаг для команды init

    // Приватная функция для вывода справки
    void showHelp() {
    logger.info("Showing help information");
    
    // Содержимое справки
    const std::string helpContent = R"(
# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки и ее радиуса
  

## 💻️ Системные требования

### Обязательные компоненты
1. **Компилятор C++17**  
   - `g++` (GCC) или `clang++`  
   - *Зачем*: Для сборки исходного кода программы

2. **Gnuplot 5.4+**  
   - *Что это*: Программа для построения графиков  
   - *Зачем*: Для визуализации карт высот, триангуляции и маршрутов  
   - Установка:  
     ```bash
     sudo apt install gnuplot  # Linux (Debian/Ubuntu)
     brew install gnuplot      # macOS (Homebrew)
     ```

### Опциональные компоненты
3. **CMake 3.12+**  
   - *Что это*: Система управления сборкой  
   - *Зачем*: Для упрощённой компиляции в разных ОС (если не используете прямой вызов g++)  
   - Установка:  
     ```bash
     sudo apt install cmake  # Linux
     ```

### Совместимость ОС
✅ **Полная поддержка**:  
- Linux (Ubuntu/Debian/Arch)  
- macOS (Intel/Apple Silicon)  

⚠️ **Ограниченная поддержка**:  
- Windows (требуется WSL2 или Cygwin)  
  - Рекомендуемый способ:  
    ```bash
    wsl --install -d Ubuntu
    ```

### Проверка установки
```bash
# Проверить версии компонентов
g++ --version
gnuplot --version
cmake --version
```

## 🚀 Запуск программы

### Способ 1: Автоматический скрипт (рекомендуется)
```bash
# Даем права на выполнение (только при первом запуске)
chmod +x run.sh

# Запуск с интерфейсом командной строки
./run.sh

# Запуск с файлом команд (commandsGauss.cmd)
./run.sh commands
```

### Способ 2: С CMake (опционально)
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel $(nproc)

# Запуск из корня проекта
./run.sh commands
```

## 📂 Файлы проекта
```
~/gauss/                                # Корневая папка проекта
├── bin/                                # Ресурсы программы
│   └── etc/
│       ├── commands/                   # Файлы команд
│       │   ├── commandsGauss.cmd       # Пример с гауссовыми холмами
│       │   └── commandsRead.cmd        # Пример с загрузкой BMP
│       └── config.conf                 # Основные параметры
├── build/                              # Сборка (генерируется CMake)
├── src/                                # Исходный код
│   ├── algorithms/                     # Алгоритмы обработки
│   ├── app/                            # Главный исполняемый модуль
│   ├── command/                        # Система команд
│   ├── core/                           # Ядро программы
│   ├── io/                             # Ввод-вывод
│   ├── utils/                          # Вспомогательные утилиты
│   └── visualization/                  # Визуализация
├── results/                            # Результаты работы
│   ├── help.txt                        # Справочник по командам
│   └── visualizations/                 # Графические результаты
│       ├── *.png                       # 2D визуализации
│       ├── *.bmp                       # Бинарные карты
│       └── Plot3DPath.png              # 3D маршрут
├── var/                                # Системные файлы
│   └── logs/                           # Логи программы
├── .gitignore                          # Игнорируемые файлы
├── CMakeLists.txt                      # Конфигурация сборки
├── LICENSE                             # Лицензия
└── run.sh                              # Автоматический скрипт запуска
```


## 🛠 Команды управления (для командного файла command.txt)

| Команда              | Параметры                      | Описание                                                                 |
|----------------------|--------------------------------|--------------------------------------------------------------------------|
| help                 | -                              | Создание файла с пояснением команд                                       |
| init                 | -                              | Инициализация поля                                                       |
| g                    | x y sx sy h                    | Создает гаусс с центром (x,y), размерами (sx,sy) и высотой h             |
| generate             | -                              | Складывает все добавленные гауссы в итоговое поле                        |
| gnuplot              | filename.png                   | Сохраняет 3D-визуализацию поля в PNG файл                                |
| PlotMetedata         | filename.png                   | Визуализирует метаданные компонент с границами и центрами                |
| PlotVoronoi          | filename.png                   | Строит диаграмму Вороного по текущей триангуляции                        |
| PlotDelaunay         | filename.png                   | Визуализирует триангуляцию Делоне                                        |
| PlotPath             | filename.png                   | Отображает найденный путь между точками A и B                            |
| bmp_write            | filename.bmp [Full/Binary]     | Сохраняет поле в BMP: Full - полное, Binary - бинаризованное             |
| bmp_read             | filename.bmp                   | Загружает поле из BMP файла                                              |
| bin                  | slice [Peaks/Valleys/All]      | Бинаризация: Peaks - только пики, Valleys - впадины, All - по модулю     |
| wave                 | noisy                          | Удаляет компоненты размером ≤ noisy как шум                              |
| k_means              | k                              | Кластеризует данные в k кластеров                                        |
| k_means_kern         | kk                             | Кластеризация с ядрами размера kk                                        |
| triangulate          | -                              | Строит триангуляцию Делоне по центрам компонент                          |
| find_path            | Ax Ay Bx By                    | Ищет путь между точками A и B через триангуляцию                         |
| Plot3DPath           | filename.png                   | Сохраняет 3D-визуализацию путя в PNG файл                                |
| plotInteractive3DPath| -                              | Интерактвный 3D режим с путем                                            |
| end                  | -                              | Завершает работу программы                                               |


## ⚙️ Параметры конфигурационного файла (config.txt)

| Параметр                     | Значение                                       | Описание                                                                         |
|------------------------------|------------------------------------------------|----------------------------------------------------------------------------------|
| fieldWidth                   | `fieldWidth`                                   | Ширина рабочего поля в пикселях                                                  |  
| fieldHeight                  | `fieldHeight`                                  | Высота рабочего поля в пикселях                                                  |
| defaultX                     | `defaultX`                                     | Стандартная X-координата центра гауссова распределения по умолчанию              |
| defaultY                     | `defaultY`                                     | Стандартная Y-координата центра гауссова распределения по умолчанию              |
| defaultSx                    | `defaultSx`                                    | Стандартное отклонение по оси X по умолчанию                                     |
| defaultSy                    | `defaultSy`                                    | Стандартное отклонение по оси Y по умолчанию                                     |
| defaultH                     | `defaultH`                                     | Стандартная высота гауссова распределения по умолчанию                           |
| defaultGnuplot               | `filename_gnuplot.png`                         | Путь к файлу для сохранения 3D-визуализации по умолчанию                         |
| defaultPlotMetedata          | `filename_metadata.png`                        | Путь к файлу для визуализации метаданных компонент по умолчанию                  |
| defaultPlotVoronoi           | `filename_voronoi.png`                         | Путь к файлу для диаграммы Вороного по умолчанию                                 |
| defaultPlotDelaunay          | `filename_delaunay.png`                        | Путь к файлу для триангуляции Делоне по умолчанию                                |
| defaultPlotPath              | `filename_path.png`                            | Путь к файлу для визуализации маршрута по умолчанию                              |
| defaultWrite                 | `filename_write.bmp`                           | Путь к файлу для сохранения BMP-изображения по умолчанию                         |
| defaultWriteModeImage        | `writeMode`                                    | Режим сохранения BMP (Full/Binary) по умолчанию                                  |
| defaultRead                  | `filename_read.bmp`                            | Путь к файлу для загрузки BMP-изображения по умолчанию                           |
| defaultSlice                 | `defaultSlice`                                 | Порог бинаризации по умолчанию                                                   |
| defaultBinMode               | `binMode`                                      | Режим бинаризации (Peaks/Valleys/All) по умолчанию                               |
| defaultNoisy                 | `defaultNoisy`                                 | Порог для удаления шумовых компонент по умолчанию                                |
| defaultKlaster               | `defaultKlaster`                               | Количество кластеров для k-mean по умолчанию                                     |
| defaultKlasterKern           | `defaultKlasterKern`                           | Размер ядра для кластеризации по умолчанию                                       |
| defaultpointA_x              | `pointA_x`                                     | X-координата точки A для поиска пути по умолчанию                                |
| defaultpointA_y              | `pointA_y`                                     | Y-координата точки A для поиска пути по умолчанию                                |
| defaultpointB_x              | `pointB_x`                                     | X-координата точки B для поиска пути по умолчанию                                |
| defaultpointB_y              | `pointB_y`                                     | Y-координата точки B для поиска пути по умолчанию                                |
| defaultPlot3DPath            | `filename_plot3dpath.png`                      | Путь к файлу для 3D-визуализации маршрута по умолчанию                           |
| vehicleRadius                | `vehicleRadius`                                | Радиус транспортного средства                                                    |
| maxSideAngle                 | `maxSideAngle`                                 | Максимальный угол поворота вбок (градусы)                                        |
| maxUpDownAngle               | `maxUpDownAngle`                               | Максимальный угол наклона вверх/вниз (градусы)                                   |
| logFileNameInterface         | `filename_log_interface.txt`                   | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | `filename_log_control.txt`                     | Путь к лог-файлу управления                                                      |
| defaultHelp                  | `/home/log/Gauss/results/docs/help.txt`        | Путь где сохранить help файл                                                     |
| FiltrationLogLevelInterface  | `logLevelInterface`                            | Уровень логирования интерфейса (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |
| FiltrationLogLevelControl    | `logLevelControl`                              | Уровень логирования управления (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |


## ⚠️ Важно
1. Всегда начинайте с команды init
2. Точки A/B задаются в config.txt
3. Маршрут будет найден не всегда!
4. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
5. Если пользуетесь программой, то важно использовать ту же файловую структуру!
6. Путь к файлу пишем относительный
7. Важен порядок команд, не забывайте делать картинки после команд
8. Для команды bmp_write не полагайтесь на значения по умолчанию
9. Точки A и B должны попадать в триангуляцию
10. Уровень "равнины" = 127, чтобы метод записи поля по гаусам согласовался с записью по картике
11. Условия проходимости: 1)Если угол наклона в любом пикселе путя превосходит допустимый (по направлению или вбок) 2) Расстояние до препятсвия на срезе меньше радиуса

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
bmp_write results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path
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
bmp_write results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path 60 130 150 135
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
```

## 📃️ Конфигурационный файл (пример)
```
fieldWidth 250
fieldHeight 250
defaultX 50.0
defaultY 50.0
defaultSx 20.0
defaultSy 20.0
defaultH 200.0
defaultGnuplot results/visualizations/Gnuplot.png
defaultPlotMetedata results/visualizations/Metadata.png
defaultPlotVoronoi results/visualizations/Voronoi.png
defaultPlotDelaunay results/visualizations/Delaunay.png
defaultPlotPath results/visualizations/Path.png
defaultWrite results/visualizations/Write.bmp 
defaultWriteModeImage Full
defaultRead results/visualizations/Read.bmp
defaultSlice 127
defaultBinMode All
defaultNoisy 10
defaultKlaster 5
defaultKlasterKern 5
defaultpointA_x 150.0
defaultpointA_y 150.0
defaultpointB_x 160.0
defaultpointB_y 160.0
defaultPlot3DPath results/visualizations/Plot3DPath.png
vehicleRadius 5
maxSideAngle 90.0
maxUpDownAngle 90.0
logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt
defaultHelp results/help.txt
FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO
```

##### ======================================
## 🔄 SAFE WORKFLOW v1.0
##### ======================================
#### Философия: "Мои изменения — священны 😤️, main обновляется без боли 😇️"

### 1️⃣ Начало работы (без опасного pull!)
```
git checkout main
git checkout -b feature/improved-structure  # Создаем ветку ОТ локального main
```

### 2️⃣ Работа с файлами
#### Редактируете файлы -> сохраняете -> проверяете:
```
git status
```

### 3️⃣ Фиксация изменений в ветке
```
git add .
git commit -m "Мои изменения ..."
git push origin feature/improved-structure
```

### 4️⃣ Лучший способ обновить main:

#### --- Способ 1: Идеальный мир (без конфликтов) ---
```
git checkout main
git merge --squash feature/improved-structure #Лень решать конфликты? Пишем "git reset --merge" после и идем ко 2 способу!
git commit -m "РЕЛИЗ: Новая структура проекта"
git push origin main
```

#### --- Способ 2: Жёсткий reset (когда конфликты лень решать) ---
```
git checkout main
git reset --hard feature/improved-structure
git push origin main --force-with-lease
```

#### --- Способ 3: Аккуратный мерж (если main меняли другие) ---
```
git fetch origin
git merge origin/main --no-commit
[ручное разрешение конфликтов]
git commit -m "Мерж улучшенной структуры"
git push origin main
```

### 5️⃣ Создание тега
```
git tag -a v2.0.0 -m "Обновленная структура проекта"
git push origin v2.0.0
```

### 6️⃣ Очистка
```
git branch -D feature/improved-structure
git push origin --delete feature/improved-structure
```


## 📄 Лицензия
Этот проект лицензирован под MIT License. Вы можете свободно использовать, изменять и распространять код, при условии, что вы укажете автора.

Разрешенные действия:

   1. Использование кода в коммерческих и некоммерческих проектах
   2. Модификация кода
   3. Распространение кода

Обязательное условие: при использовании кода, пожалуйста, укажите ссылку на автора.

Developed with ❤️ by **DebugDestroy**  
[GitHub Profile](https://github.com/DebugDestroy)
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
    void processFileCommands(std::ifstream& file) {
        int commandCount = 0;
        
        while (file >> params.s) {
            commandCount++;
            logger.info(std::string("Processing command #") + std::to_string(commandCount) + ": " + params.s);
            
            if (params.s == "help") {
                showHelp();
                continue;
            }
            
            if (params.s == "end") {
                logger.info("Ending the program.");
                break;
            }
            
            if (!processCommand(file, false)) {
                break;
            }
        }
    }

    // Обработка команд с клавиатуры
    void processKeyboardCommands() {
        
        while (true) {
            const std::string commandshow = R"(Enter the command and its parameters immediately (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, wave,
             PlotMetedata, PlotVoronoi, PlotDelaunay, PlotPath, k_means, k_means_kern, triangulate, find_path, Plot3DPath, plotInteractive3DPath, end):)";
            std::cout << commandshow;
            std::cin >> params.s;
            std::cout << "\n";
            logger.info(std::string("Received command: ") + params.s);
            
            if (params.s == "help") {
                showHelp();
                continue;
            }
            
            if (params.s == "end") {
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
    bool processCommand(std::istream& input, bool fromKeyboard) {
    std::string line;
    std::string showInfo;
    std::string modeWrite, modeBin;
    
    if (params.s == "init") {
        if (n != 0) {
            std::cout << "The init command has already been called.\nError\n";
            logger.error("Error: Multiple init commands.");
            return false;
        }
        n = 1;
        showInfo = std::string("Initializing field with size: ") + std::to_string(params.A) + " x " + std::to_string(params.B);
        
        if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Field initialized.");
    }
    else if (n != 1) {
        std::cout << "The init command was not used.\nError\n";
        logger.error("Error: The init command was not used.");
        return false;
    }
    else if (params.s == "g") {      
        params.x = config.defaultX;
        params.y = config.defaultY;
        params.sx = config.defaultSx;
        params.sy = config.defaultSy;
        params.h = config.defaultH;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.x >> params.y >> params.sx >> params.sy >> params.h;
      
         showInfo = std::string("Adding Gaussian: x=") + std::to_string(params.x) + 
                  ", y=" + std::to_string(params.y) + 
                  ", sx=" + std::to_string(params.sx) + 
                  ", sy=" + std::to_string(params.sy) + 
                  ", h=" + std::to_string(params.h);
                  
         if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info(showInfo);
        control.Dispetcher(params);
    }
    else if (params.s == "generate") {
        showInfo = "Generating field by summing all Gaussians";
        logger.info(showInfo);
        control.Dispetcher(params);
        if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        logger.info("Field generation completed");
    }
    else if (params.s == "gnuplot") {
        params.filename = config.defaultGnuplot;

            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;

            showInfo = std::string("Calling gnuplot with filename: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Gnuplot visualization completed");
    }
    else if (params.s == "PlotMetedata") {
        params.filename = config.defaultPlotMetedata;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
            showInfo = std::string("Plotting metadata to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Metadata plotting completed");
    }
    else if (params.s == "PlotVoronoi") {
        params.filename = config.defaultPlotVoronoi;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
        showInfo = std::string("Plotting Voronoi diagram to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }        

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Voronoi diagram plotting completed");
    }
    else if (params.s == "PlotDelaunay") {
        params.filename = config.defaultPlotDelaunay;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
             
        showInfo = std::string("Plotting Delaunay triangulation to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation plotting completed");
    }
    else if (params.s == "PlotPath") {
        params.filename = config.defaultPlotPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
             
        showInfo = std::string("Plotting path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Path plotting completed");
    }
    else if (params.s == "bmp_write") {
        params.filename = config.defaultWrite;
        modeWrite = config.defaultWriteModeImage;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename >> modeWrite;
         
        if (modeWrite == "Full") {
            params.bmp_mode = BmpWriteMode::Full;
        } else {
            params.bmp_mode = BmpWriteMode::Binary;
        }
                    
        showInfo = std::string("Writing BMP file: ") + params.filename + " with mode: " + modeWrite;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("BMP writing completed");
    }
    else if (params.s == "bmp_read") {
        params.filename = config.defaultRead;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
        showInfo = std::string("Reading BMP file: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("BMP reading completed");
    }
    else if (params.s == "bin") {
        params.slice = config.defaultSlice;
        modeBin = config.defaultBinMode;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.slice >> modeBin;

        if (modeBin == "Peaks") {
            params.bin_mode = ThresholdMode::Peaks;
        } else if (modeBin == "Valleys") {
            params.bin_mode = ThresholdMode::Valleys;
        } else {
            params.bin_mode = ThresholdMode::All;
        }
                  
        showInfo = std::string("Applying binary filter with slice: ") + std::to_string(params.slice) + " and mode: " + modeBin;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
        
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Binary filtering completed");
    }
    else if (params.s == "wave") {
        params.noisy = config.defaultNoisy;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.noisy;
                  
        showInfo = std::string("Applying wave filter with noisy level: ") + std::to_string(params.noisy);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info(std::string("Wave filtering completed. Components count: ") + 
                   std::to_string(control.componenti.size()));
    }
    else if (params.s == "k_means") {
        params.k = config.defaultKlaster;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.k;
                  
        showInfo = std::string("Running k-means with k: ") + std::to_string(params.k);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("k-means clustering completed");
    }
    else if (params.s == "k_means_kern") {
        params.kk = config.defaultKlasterKern;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.kk;
                  
        showInfo = std::string("Running k-means with kernel size: ") + std::to_string(params.kk);
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("k-means with kernels completed");
    }
    else if (params.s == "triangulate") {              
        showInfo = "Starting Delaunay triangulation";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }
         
        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation completed");
    }
    else if (params.s == "find_path") {
        params.pointA_x = config.defaultpointA_x;
        params.pointA_y = config.defaultpointA_y;
        params.pointB_x = config.defaultpointB_x;
        params.pointB_y = config.defaultpointB_y;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.pointA_x >> params.pointA_y >> params.pointB_x >> params.pointB_y;
            
        showInfo = std::string("Finding path from (") + std::to_string(params.pointA_x) + "," + 
                   std::to_string(params.pointA_y) + ") to (" + 
                   std::to_string(params.pointB_x) + "," + 
                   std::to_string(params.pointB_y) + ")";
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("Path finding completed");
    }
    else if (params.s == "Plot3DPath") {
        params.filename = config.defaultPlot3DPath;
        
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
            
        showInfo = std::string("Plotting 3D path to: ") + params.filename;
            
            if (fromKeyboard) {
            std::cout << "\n";
            std::cout << showInfo << std::endl;
         }

        logger.info(showInfo);
        control.Dispetcher(params);
        logger.info("3D path plotting completed");
    }
    else if (params.s == "plotInteractive3DPath") {
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
        std::cout << "Unknown command: " << params.s << std::endl;
        logger.warning(std::string("Unknown command received: ") + params.s);
        return false;
    }
    
    return true;
}

public:
    Interface(Config& cfg, Logger& log, Control& c) 
        : config(cfg), logger(log), control(c) {
        logger.info("Interface initialized");
        logger.debug(std::string("Field dimensions: ") + std::to_string(config.fieldWidth) + 
                   "x" + std::to_string(config.fieldHeight));
        
        params.A = config.fieldWidth;
        params.B = config.fieldHeight;
    }
    
    void print() {
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
};
