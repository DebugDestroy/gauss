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
![delaunay](results/visualizations/Delaunay.png)

- Построение диаграммы Вороного
![voronoi](results/visualizations/Voronoi.png)

- Навигационный граф, построенный по диаграмме Вороного с ограничениями по углу наклона тележки и ее радиуса
![graph](results/visualizations/Graph.png)

- Поиск пути
![path](results/visualizations/AstarPlot3DPath.png)


## 🚀 Запуск
```bash
chmod +x run.sh
./run.sh
```

## 📚 Документация
- 📘 [User Guide](docs/User%20Guide.md) — описание всех команд и сценариев использования  
- 🛠 [Developer Guide](docs/Developer%20Guide.md) — архитектура и расширение системы
- 🚧 [WIP / Research Notes](docs/WIP.md) — текущие задачи, исследовательские направления и последние изменения проекта

## 📄 Лицензия
Этот проект лицензирован под MIT License. Вы можете свободно использовать, изменять и распространять код, при условии, что вы укажете автора.

Разрешенные действия:

   1. Использование кода в коммерческих и некоммерческих проектах
   2. Модификация кода
   3. Распространение кода

Обязательное условие: при использовании кода, пожалуйста, укажите ссылку на автора.

Developed with ❤️ by **DebugDestroy**  
[GitHub Profile](https://github.com/DebugDestroy)
