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
~/gauss/                      # Корневая папка проекта
.
├── CMakeLists.txt            # Главный файл сборки CMake
├── LICENSE                   # Лицензия проекта
├── README.md                 # Этот файл
├── run.sh                    # Скрипт для запуска проекта
├── User Guide.pdf            # Руководство пользователя
├── config/                   # Конфигурационные файлы
│   ├── commands/             # Примеры команд
│   │   ├── commandsGauss.cmd
│   │   └── commandsRead.cmd
│   └── config.conf           # Параметры конфигурации
├── results/                  # Результаты работы
│   ├── help.txt              # Подсказки по командам
│   └── visualizations/       # Графические выходные данные
│       ├── *.png             # Визуализации диаграмм, маршрутов и т.п.
│       └── *.bmp             # Бинарные изображения
├── src/                      # Исходный код
│   ├── algorithms/           # Алгоритмы и обработка данных
│   │   ├── components/       # Компонентный анализ, кластеризация
│   │   ├── gauss/            # Генерация гауссовых холмов
│   │   ├── geometry/         # Геометрические алгоритмы (триангуляция, Вороной и т.д.)
│   │   ├── kinematics/       # Расчёт углов наклона и крена
│   │   └── path/             # Алгоритмы поиска пути (A*, Dekstra, Greedy)
│   ├── app/                  # Точка входа (main.cpp)
│   ├── command/              # Интерпретация команд пользователя
│   ├── core/                 # Базовые сущности и логика (конфиг, логгирование и т.д.)
│   ├── io/                   # Ввод/вывод, работа с изображениями
│   ├── utils/                # Утилиты, хеши и т.п.
│   └── visualization/        # Модули визуализации (gnuplot, цвета)
├── var/                      # Временные/системные файлы
│   └── logs/                 # Логи работы программы

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
| find_path_astar      | Ax Ay Bx By                    | A* ищет путь между точками A и B через триангуляцию                      |
| find_path_dekstra    | Ax Ay Bx By                    | Dekstra ищет путь между точками A и B через триангуляцию                 |
| find_path_greedy     | Ax Ay Bx By                    | Greedy ищет путь между точками A и B через триангуляцию                  |
| Plot3DPath           | filename.png                   | Сохраняет 3D-визуализацию путя в PNG файл                                |
| plotInteractive3DPath| -                              | Интерактвный 3D режим с путем                                            |
| end                  | -                              | Завершает работу программы                                               |


## ⚙️ Параметры конфигурационного файла (config.txt)

| Параметр                     | Значение                                       | Описание                                                                         |
|------------------------------|------------------------------------------------|----------------------------------------------------------------------------------|
| fieldWidth                   | `fieldWidth`                                   | Ширина рабочего поля в пикселях                                                  |  
| fieldHeight                  | `fieldHeight`                                  | Высота рабочего поля в пикселях                                                  |
| defaultCenterX               | `defaultCenterX`                               | Стандартная X-координата центра гауссова распределения по умолчанию              |
| defaultCenterY               | `defaultCenterY`                               | Стандартная Y-координата центра гауссова распределения по умолчанию              |
| defaultSigmaX                | `defaultSigmaX`                                | Стандартное отклонение по оси X по умолчанию                                     |
| defaultSigmaY                | `defaultSigmaY`                                | Стандартное отклонение по оси Y по умолчанию                                     |
| defaultHeight                | `defaultHeight`                                | Стандартная высота гауссова распределения по умолчанию                           |
| defaultGnuplot               | `filename_gnuplot.png`                         | Путь к файлу для сохранения 3D-визуализации по умолчанию                         |
| defaultPlotMetedata          | `filename_metadata.png`                        | Путь к файлу для визуализации метаданных компонент по умолчанию                  |
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
find_path_astar
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
find_path_astar 60 130 150 135
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
```
3) Крутое гауссово поле для размера 300 x 300
```
init

g 150 180 5 6 10
g 60 50 9 12 18
g 50 130 7 6 16
g 10 120 10 20 -10
g 80 75 10 15 -8
g 200 250 7 6 16
g 100 100 5 5 16
g 120 130 17 7 12
g 50 200 10 10 20
g 250 20 10 10 -29
g 170 50 20 10 20
g 10 10 10 10 7
g 290 290 10 10 7
g 290 10 10 10 7
g 10 290 10 30 7
g 230 90 20 10 7
g 180 130 10 20 10
g 200 190 20 10 7
g 230 90 10 10 7
g 230 90 10 10 7
g 30 200 10 10 7
g 100 250 8 8 10
g 200 100 10 10 -10
g 130 220 10 10 10
g 170 170 8 8 12
g 250 250 8 8 -20
g 270 150 6 6 25
g 220 30 10 10 10
g 140 40 6 10 -15
g 90 180 8 8 10
g 160 220 10 10 -12
g 70 250 6 6 10
g 270 70 10 10 -15
g 130 70 8 8 10
g 200 60 6 6 10
g 60 160 8 8 -10
g 110 30 6 6 10
g 180 30 8 8 8
g 250 120 6 6 10
g 40 90 6 6 -10
g 220 220 10 10 15


generate
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full
bin 132 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 10
bmp_write results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path_astar 20 27 100 298
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
```

## 📃️ Конфигурационный файл (пример)
```
fieldWidth 250
fieldHeight 250
defaultCenterX 50.0
defaultCenterY 50.0
defaultSigmaX 20.0
defaultSigmaY 20.0
defaultHeight 200.0
defaultGnuplot results/visualizations/Gnuplot.png
defaultPlotMetedata results/visualizations/Metadata.png
defaultPlotVoronoi results/visualizations/Voronoi.png
defaultPlotDelaunay results/visualizations/Delaunay.png
defaultPlotPath results/visualizations/Path.png
defaultWrite results/visualizations/Write.bmp 
defaultWriteModeImage Full
defaultRead results/visualizations/Read.bmp
defaultThreshold 130
defaultBinMode All
defaultNoisy 10
defaultKlaster 5
defaultKlasterKern 5
defaultstartPointX 150.0
defaultstartPointY 150.0
defaultendPointX 160.0
defaultendPointY 160.0
defaultPlot3DPath results/visualizations/Plot3DPath.png
vehicleRadius 1
maxSideAngle 90.0
maxUpDownAngle 90.0
logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt
defaultHelp results/help.txt
FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO
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
