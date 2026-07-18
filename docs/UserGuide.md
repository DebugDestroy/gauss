# User Guide

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

## Как скачать программу
- С гитхаба загрузите последнюю версию
- Распакуйте и перетащите gauss в домашнюю папку

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

## 📂 Файлы проекта
```
~/gauss/                      # Корневая папка проекта
.
├── CMakeLists.txt            # Главный файл сборки CMake
├── LICENSE                   # Лицензия проекта
├── notebooks                 # Исследования
├── README.md                 # Витрина проекта
├── docs                      # Документы
│   ├── UserGuide.md          # Гайд пользователю (этот файл)
│   ├── WIP.md                # План и незакомиченные изменения
│   └──DeveloperGuide.md      # Гайд разработчику
├── run.sh                    # Скрипт для запуска проекта
├── config/                   # Конфигурационные файлы
│   ├── commands/             # Примеры команд
│   └── config.conf           # Параметры конфигурации
├── results/                  # Результаты работы
│   ├── showcase/             # Красивые картинки
│   └── visualizations/       # Результат работы
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
│   ├── statistics/           # Сбор статистики
│   ├── utils/                # Утилиты, хеши и т.п.
│   └── visualization/        # Модули визуализации (gnuplot, цвета)
├── var/                      # Временные/системные файлы
│   └── logs/                 # Логи работы программы

```


## 🛠 Команды управления (для командного файла command.txt)

| Команда              | Тип данных + *параметры*                                                                           | Описание                                                                 |
|----------------------|----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| help                 | -                                                                                                  | Создание файла с пояснением команд                                       |
| init                 | int *fieldWidth fieldHeight*                                                                 | Инициализация поля с заданой шириной и длиной. field[y][x] X-(width) Y-(height)|
| g                    | double *x y sx sy h*                                                                               | Создает гаусс с центром (x,y), размерами (sx,sy) и высотой h             |
| g_auto               | int *count_min count_max* double *xmin xmax ymin ymax ... h_min h_max*                             | Создает случайные гаусы с параметрами в промежутках                      |
| g_grid               | int *g_cell_size*                                                                                  |Строит пространственную сетку для ускорения вычислений высот в непр случае|
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
| PlotPathDiscrete     | string *filename.png*                                                                              | Отображает найденный путь между точками A и B для дискретного путя       |
| PlotRRT              | string *filename.png*                                                                              | Делает gif изображение построения дерева и пути RRT                      |
| PlotRRTStar          | string *filename.png*                                                                              | Делает gif изображение построения дерева и пути RRT*                     |
| PlotPathContinuous   | string *filename.png*                                                                              | Отображает найденный путь между точками A и B для непрерывного путя      |
| bmp_write            | string *filename.bmp [Full/Binary]*                                                                | Сохраняет поле в BMP: Full - полное, Binary - бинаризованное             |
| bmp_read             | string *filename.bmp*                                                                              | Загружает поле из BMP файла                                              |
| bin                  | int *slice*                                                                                        | Бинаризация с уровнем отклонения от равнины MID_GRAY                     |
| wave                 | int *noisy*                                                                                        | Удаляет компоненты размером ≤ noisy как шум                              |
| k_means              | int *k*                                                                                            | Кластеризует данные в k кластеров                                        |
| k_means_kern         | int *k p*                                                                                          | Kmeans с параметром k на ядрах размера p                                 |
| triangulate          | -                                                                                                  | Строит триангуляцию Делоне по центрам компонент                          |
| voronoi              | -                                                                                                  | Строит диаграмму Вороного                                                |
| build_nav_graph      | int *vehicleRadius* double *maxSideAngle maxUpDownAngle*                                           | Строит граф проходимости по диаграмме Вороного                           |
| grid                 | int *grid_cell_size*                                                                        | Строит квадратную сетку на поле (ячейка размера grid_cell_size x grid_cell_size)|
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
|rrt|size_t *rebuildSize maxIterations* double *Ax Ay Bx By vehicleRadius heightThreshold maxSideAngle maxUpDownAngle interpEdge interpCollision interpAngle step goalRadius goalBias*|Строит путь соблюдая условия|
|rrt_star|size_t *rebuildSize maxIterations* double *Ax Ay Bx By vehicleRadius heightThreshold maxSideAngle maxUpDownAngle interpEdge interpCollision interpAngle step maxFindRadius gammaConstant goalRadius goalBias*|Строит путь соблюдая условия|
| shortcut_discrete    | -                                                                                                  | Удаляет лишние вершины в дискретном пути, сокращая путь                  |
| shortcut_continuous  | -                                                                                                  | Удаляет лишние вершины в непрерывном пути, сокращая путь                 |
| spline_discrete      |size_t *samplesPerSegment* double *vehicleRadius heightThreshold maxSideAngle maxUpDownAngle interpEdge interpCollision interpAngle*| Делает путь плавным                      |
| spline_continuous    |size_t *samplesPerSegment*                                                                          | Делает путь плавным                                                      |


### Замечания к командному файлу
1. Обязательно полностью настроить конфиг файл
2. Важно строго соблюдать последовательность команд
```
Для графовых алгоритмов поиска пути:
init -> g несколько раз или один раз g_auto или bmp_read ->  generate -> bin -> wave -> kmeans/k_means_kern (не влияет на работу) -> triangulate -> voronoi -> build_nav_graph -> connect_to_graph -> 
astar_graph

Для алгоритмов поиска пути на сетке:
init -> g несколько раз или один раз g_auto или bmp_read ->  generate -> bin -> wave -> kmeans/k_means_kern (не влияет на работу) -> grid -> build_nav_grid -> connect_to_grid -> astar_grid

Для RRT:
init -> g несколько раз или один раз g_auto -> rrt (Но для визуализации нужно всю равно  сделать generate -> bin, чтобы была дискретная матрица для визуализации)
```
3. Не забывайте визуализировать после алгоритмов
4. Всегда начинайте с команды init
5. Для команды bmp_write не полагайтесь на значения по умолчанию
6. Для команды PlotKmeans не полагайтесь на значения по умолчанию (иначе вывод совпадет для kmeans и kmeans_with_kernels)
7. Команда connect_to_graph имеет режим NearestK k (команда+число). Она имеет свой численный параметр, его нужно указать в конфиге после соотв команды либо в командном файле после режима
8. Не забывайте в конце писать команду end для завершения программы


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
| g_cell_size                  | int                                            | Размер сетки по умолчанию                                                        |
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
| PlotPathDiscrete             | string                                         | Путь к файлу для визуализации дискретного маршрута по умолчанию                  |
| PlotRRT                      | string                                         | Путь к файлу для визуализации rrt по умолчанию                                   |
| PlotRRTStar                  | string                                         | Путь к файлу для визуализации rrt* по умолчанию                                  |
| PlotPathContinuous           | string                                         | Путь к файлу для визуализации непрерывного маршрута по умолчанию                 |
| defaultWrite                 | string                                         | Путь к файлу для сохранения BMP-изображения по умолчанию                         |
| defaultWriteModeImage        | [Full/Binary]                                  | Режим сохранения BMP (Full/Binary) по умолчанию                                  |
| defaultRead                  | string                                         | Путь к файлу для загрузки BMP-изображения по умолчанию                           |
| heightThresholdPixel         | int                                            | Порог бинаризации по умолчанию                                                   |
| defaultWaveNoisy             | int                                            | Порог для удаления шумовых компонент по умолчанию                                |
| defaultKlaster               | int                                            | Количество кластеров для k-mean по умолчанию                                     |
| defaultKlasterKern           | int                                            | Размер ядра для кластеризации по умолчанию                                       |
| grid_cell_size               | int                                            | Размер ячейки для сетки по умолчанию                                             |
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
| rebuildSize                  | size_t                                         | Частота обновления kd-tree                                                       |
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
| maxFindRadius                | double                                         | Максимальный радиус для присоединения соседей к новой вершине по умолчанию       |
| gammaConstant                | double                                         | Константа RRT* для пересчета радиуса подключения к соседям по умолчанию          |
| goalRadius                   | double                                         | Радиус круга цели, внутри которого пробуем присоединить оказавшиеся там вершины  |
| goalBias                     | double                                         | Вероятность оказаться случайной точке у цели                                     |
| samplesPerSegment            | size_t                                         | На сколько маленьких кусочков разбить один участок сплайна                       |
| defaultsave_metrics          | string                                         | Путь к файлу по умолчанию для сохранения метрик                                  |
| logFileNameInterface         | string                                         | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | string                                         | Путь к лог-файлу управления                                                      |
| defaultHelp                  | string                                         | Путь где сохранить help файл                                                     |
| FiltrationLogLevelInterface  | [TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF]  | Уровень логирования интерфейса                                                   |
| FiltrationLogLevelControl    | [TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF]  | Уровень логирования управления                                                   |
| seedMode                     | [Random/Fixed size_t]                          | Режим случайной генерации (Random/Fixed seed) по умолчанию                       |


### Замечания к конфигурационному файлу
1. Конфиг файл обязательно должен быть полностью задан иначе программа не будет работать
2. Команды seedMode в режиме Fixed и defaultConnectMode в режиме Nearest имеют доп численный параметр который нужно обязательно указать

## ⚠️ Важно
1. Маршрут будет найден не всегда!
2. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
3. Если пользуетесь программой, то важно использовать ту же файловую структуру!
4. Путь к файлу пишем относительный
5. Важен порядок команд, не забывайте делать картинки после команд
6. Точки A и B должны попадать в безопасную зону иначе путь не будет найден
7. Уровень "равнины" = 127, чтобы метод записи поля по гаусам согласовался с записью по картике bmp
8. Условия проходимости: 1) Если угол наклона в любом пикселе путя не превосходит допустимый (по направлению или вбок) 2) Расстояние до препятсвия на срезе меньше радиуса
9. Если что-то не работает смотрите логи. Возможно опечатались в командном или конфиг файле

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
PlotPathDiscrete results/visualizations/Path.png
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
PlotPathDiscrete results/visualizations/Path.png
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
PlotKmeans results/visualizations/kmeans.png
triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph 20 27 100 298
astar_graph
save_metrics
PlotPathDiscrete results/visualizations/Path.png
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


g_cell_size 10


defaultGnuplot results/visualizations/Gnuplot.png
defaultPlotMetedata results/visualizations/Metadata.png
defaultPlotKmeans results/visualizations/kmeans.png
defaultPlotGraph results/visualizations/Graph.png
defaultPlotVoronoi results/visualizations/Voronoi.png
defaultPlotDelaunay results/visualizations/Delaunay.png
defaultPlotGrid results/visualizations/Grid.png
defaultPlotNavGrid results/visualizations/NavGrid.png
defaultPlotGridPath results/visualizations/GridPath.png
PlotPathDiscrete results/visualizations/PathDiscrete.png
PlotPathContinuous results/visualizations/PathContinuous.png
PlotRRT results/visualizations/RRT.gif
PlotRRTStar results/visualizations/RRTStar.gif
defaultPlot3DPath results/visualizations/Plot3DPath.png

defaultWrite results/visualizations/Write.bmp 
defaultWriteModeImage Full
defaultRead results/visualizations/Read.bmp

save_g config/commands/gaussians.txt


heightThresholdPixel 3


defaultWaveNoisy 10


defaultKlaster 5
defaultKlasterKern 5


grid_cell_size 1
defaultgridNoisy 0


startPixelX 150
startPixelY 150
goalPixelX 160
goalPixelY 160
defaultConnectMode All

vehicleRadiusPixel 1
maxSideAngle 90.0
maxUpDownAngle 90.0


rebuildSize 20
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
maxFindRadius 20.0
gammaConstant 100.0
goalRadius 2.0
goalBias 0.2


samplesPerSegment 20


defaultsave_metrics var/metrics/metrics.csv


logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt

FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO


seedMode Fixed 42
```
