/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память
3) Лит обзор добавить
4) Также при желании добавить правила проекта
5) Автогенерация гауссов
6) Повернуть колеса в  condition
7) Разные способы подключения старта/финиша к графу

Изменения:
1) Добавил новые алгоритмы пути (A*, Dekstra, Greedy)
2) Теперь стартовая/конечная точки присоединяются ко всем видимым вершинам диаграммы воронова (а не только к ближайшим)
3) Небольшой рефакторинг
4) В логах теперь есть основные метрики для алгоритмов путей (по этим параметрам можно сравнить алгоритмы)
                          
*/
#include "command/control.hpp"
#include "command/interface.hpp"

int main() {
    // Все пути считаются относительно корня (текущей директории)
    core::Config config("config/config.conf");
    core::Logger loggerinterface(config.logFileNameInterface, "Interface", config.FiltrationLogLevelInterface);
    core::Logger loggercontrol(config.logFileNameControl, "Control", config.FiltrationLogLevelControl);
    
    command::Control c(config, loggercontrol);
    command::Interface i(config, loggerinterface, c);
    
    i.print();
    return 0;
}
