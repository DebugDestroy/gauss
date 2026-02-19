/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Добавить новые алгоритмы пути и сравнить их
3) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память
4) Лит обзор добавить
5) Также при желании добавить правила проекта

Изменения:
1) Исправил A*, чтобы был реализован по канону (добавил closed set)
                          
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
