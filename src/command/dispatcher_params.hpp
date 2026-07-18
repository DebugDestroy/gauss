#pragma once
#include <string>
#include "io/bmp_handler.hpp"  // Для BmpWriteMode
#include "algorithms/path/common/graph.hpp" // Для СonnectMode
#include "algorithms/gauss/gauss_builder.hpp" // Для GAutoMode

namespace command {
struct DispatcherParams {
    std::string command;                                      // Команда для выполнения (help, init, g, generate и т.д.)
    
    int fieldWidth, fieldHeight;                              // Размеры поля (ширина × высота в пикселях)
    
    // ----- addgauss -----
    double height;                                            // Высота гауссова распределения (в условных единицах)
    double centerX, centerY;                                  // Координаты центра (x, y)
    double sigmaX, sigmaY;                                    // Разброс по осям (σx, σy)
    
    // ----- Параметры для g_auto -----
    std::size_t count_min, count_max;                         // Диапазон количества генерируемых гауссов
    double xmin, xmax;                                        // Диапазон X координаты центра
    double ymin, ymax;                                        // Диапазон Y координаты центра
    double sx_min, sx_max;                                    // Диапазон σx
    double sy_min, sy_max;                                    // Диапазон σy
    double h_min, h_max;                                      // Диапазон высоты
    
    // -----  g_grid ----- 
    int g_cell_size;                                          // Размер ячейки
    
    // ----- Ввод / вывод -----
    std::string filename;                                     // Имя файла для операций ввода/вывода
    io::BmpWriteMode bmpWriteMode;                            // Режим записи BMP (Full/Binary)
    
    // ----- binary -----
    int heightThresholdPixel;                                 // Бинаризиция на уровне отклонения от равнины  
    
    // ----- wave -----
    std::size_t waveNoisy;                                    // Максимальный размер шумовых компонент (в пикселях)
    
    // ----- kmeans -----
    std::size_t clusterCount;                                 // Количество кластеров для k-means
    std::size_t kernelSize;                                   // Размер ядра для кластеризации
    
    // ----- grid -----
    int grid_cell_size;                                       // Размер ячейки grid_cell_size x grid_cell_size
    std::size_t gridNoisy;                                    // Число допустимого шума в ячейке
    
    // ----- Путь дискретный-----
    int startPixelX, startPixelY;                             // Координаты начальной точки маршрута (Ax, Ay)
    int goalPixelX, goalPixelY;                               // Координаты конечной точки маршрута (Bx, By)
    algorithms::path::common::ConnectMode connectMode;        // Режим присоединению старта и финиша к графу
    std::size_t nearestVerticesCount;                         // используется только для режима NearestK
    
    int vehicleRadiusPixel;                                   // Радиус тележки
    double maxSideAngle, maxUpDownAngle;                      // Углы максимального наклона
    
    // ----- Путь непрерывный -----
    double startWorldX, startWorldY;                          // Координаты начальной точки маршрута (Ax, Ay)
    double goalWorldX, goalWorldY;                            // Координаты конечной точки маршрута (Bx, By)
    double vehicleRadiusWorld;                                // Радиус тележки
    
    // ----- rrt -----
    std::size_t rebuildSize;                                  // Частота перестройки kd-tree
    std::size_t maxIterations;                                // Предел числа итераций
    double heightThresholdWorld;                              // Допустимый уровень отклонения высоты от core::MID_GRAY = 127
    double interpEdge;                                        // С каким шагом проверять углы колес на ребре
    double interpCollision;                                   // С каким шагом проверять радиус в круге
    double interpAngle;                                       // С каким шагом проверять дугу окружности (не угол, а расстояние)
    double step;                                              // Шаг
    double goalRadius;                                        // Как близко нужно подойти к цели чтобы попробовать к ней присоединиться
    double goalBias;                                          // Вероятность генерации точки у цели
    
    // ----- rrt_star -----
    double maxFindRadius;                                     // Максимальный радиус для присоединения соседей к новой вершине
    double gammaConstant;                                     // Константа RRT* для пересчета радиуса подключения к соседям
    
    // ----- spline -----
    std::size_t samplesPerSegment;                            // На сколько маленьких кусочков разбить один участок сплайна
};
}
