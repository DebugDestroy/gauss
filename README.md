# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.


## 🧭 Pipeline обработки данных
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
![field](results/visualizations/gnuplot.png)

- Бинаризация и компоненты
![binary](results/visualizations/Metadata.png)

- Кластеризация объектов методом k-means
![kmeans](results/visualizations/kmeans.png)

- Триангуляция Делоне с учетом высот
![delaunay](results/visualizations/Triangulation_Delone.png)

- Построение диаграммы Вороного
![voronoi](results/visualizations/Diagramma_Voronova.png)

- Навигационный граф, построенный по диаграмме Вороного с ограничениями по углу наклона тележки и ее радиуса
![graph](results/visualizations/Graph.png)

- Поиск пути
![path](results/visualizations/AstarPlot3DPath.png)


## 📚 Документация
- 📘 [User Guide](docs/UserGuide.md) — описание всех команд и сценариев использования  
- 🛠 [Developer Guide](docs/DeveloperGuide.md) — архитектура и расширение системы
- 🚧 [WIP / Research Notes](docs/WIP.md) — текущие задачи, исследовательские направления и последние изменения проекта

## 📊 Исследования
- [Project 1 — Graph-based Path Planning](notebooks/project1.ipynb) — сравнение алгоритмов A*, Dijkstra и Greedy на навигационном графе, построеным на основе диаграммы Вороного

## 📄 Лицензия
Этот проект лицензирован под MIT License. Вы можете свободно использовать, изменять и распространять код, при условии, что вы укажете автора.

Разрешенные действия:

   1. Использование кода в коммерческих и некоммерческих проектах
   2. Модификация кода
   3. Распространение кода

Обязательное условие: при использовании кода, пожалуйста, укажите ссылку на автора.

Developed with ❤️ by **DebugDestroy**  
[GitHub Profile](https://github.com/DebugDestroy)
