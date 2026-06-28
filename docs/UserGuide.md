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
├── README.md                 # Витрина проекта
├── docs
│   ├── UserGuide.md          # Гайд пользователю (этот файл)
│   ├── WIP.md                # План и незакомиченные изменения
│   └──DeveloperGuide.md      # Гайд разработчику
├── run.sh                    # Скрипт для запуска проекта
├── config/                   # Конфигурационные файлы
│   ├── commands/             # Примеры команд
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

| Команда              | Параметры                                                                                          | Описание                                                                 |
|----------------------|----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| help                 | -                                                                                                  | Создание файла с пояснением команд                                       |
| init                 | -                                                                                                  | Инициализация поля                                                       |
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
| k_means_kern         | k p                                                                                                | Kmeans с параметром k на ядрах размера p                                 |
| triangulate          | -                                                                                                  | Строит триангуляцию Делоне по центрам компонент                          |
| voronoi              | -                                                                                                  | Строит диаграмму Вороного                                                |
| build_nav_graph      | vehicleRadius maxSideAngle maxUpDownAngle                                                          | Строит граф проходимости по диаграмме Вороного, учитывая наклон и радиус |
| connect_to_graph     | Ax Ay Bx By [Nearest/NearestK k/All]                                                               | Подключает старт и финиш к графу с разными режимами                      |
| find_path_astar      | -                                                                                                  | A* ищет путь между точками A и B                                         |
| find_path_dekstra    | -                                                                                                  | Dekstra ищет путь между точками A и B                                    |
| find_path_greedy     | -                                                                                                  | Greedy ищет путь между точками A и B                                     |
| Plot3DPath           | filename.png                                                                                       | Сохраняет 3D-визуализацию пути в PNG файл                                |
| plotInteractive3DPath| -                                                                                                  | Интерактвный 3D режим с путем                                            |
| end                  | -                                                                                                  | Завершает работу программы                                               |


### Замечания к командному файлу
1. Если настроен конфиг то многие параметры можно задать там один раз и не придется в каждом командном файле их указывать
2. Важно строго соблюдать последовательность команд
```
init -> g несколько раз или один раз g_auto или bmp_read ->  generate -> bin -> wave -> kmeans/k_means_kern -> triangulate -> voronoi -> build_nav_graph -> connect_to_graph -> алгоритмы поиска пути
```
3. Не забывайте визуализировать после алгоритмов
4. Всегда начинайте с команды init
5. Для команды bmp_write не полагайтесь на значения по умолчанию
6. Для команды PlotKmeans не полагайтесь на значения по умолчанию (иначе вывод совпадет для kmeans и kmeans_with_kernels)
7. Команды g_auto и connect_to_graph имеют режимы Fixed seed и NearestK k. Они имеют свой численный параметр, его нужно указать либо в конфиге либо в командном файле после этих режимов
8. Не забывайте в конце писать команду end для завершения программы


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
| logFileNameInterface         | `filename_log_interface.txt`                   | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | `filename_log_control.txt`                     | Путь к лог-файлу управления                                                      |
| defaultHelp                  | `/home/log/Gauss/results/docs/help.txt`        | Путь где сохранить help файл                                                     |
| FiltrationLogLevelInterface  | `logLevelInterface`                            | Уровень логирования интерфейса (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |
| FiltrationLogLevelControl    | `logLevelControl`                              | Уровень логирования управления (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |


### Замечания к конфигурационному файлу
1. Конфиг файл обязательно должен быть полностью задан иначе программа не будет работать


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
PlotKmeans results/visualizations/kmeans.png
triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph 60 130 150 135
find_path_astar
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
PlotKmeans results/visualizations/kmeans.png
triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph 20 27 100 298
find_path_astar
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
```

## 📃️ Конфигурационный файл (пример)
```
defaultHelp results/help.txt


fieldWidth 300
fieldHeight 300


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


logFileNameInterface var/logs/log_interface.txt
logFileNameControl var/logs/logcontrol.txt

FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO
```
