/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Возможно новые требования
3) Оптимизация
4) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память

Изменения:

                           Программа готова!!! Но всю равно следите за обновлениями на гитхаб [GitHub Profile](https://github.com/DebugDestroy)
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
