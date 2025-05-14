/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Возможно новые требования
3) Оптимизация
4) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память

Изменения:
1) Благодаря тестам понял, что предыдущая проеврка на проходимость проверяла окажется ли тележка внутри гаусса, что просто смешно! Теперь компоненты на бинарном срезе это опасная область и если тележка окажется выше или ниже, то данное ребро не проходимо, что логично! Также добавлены проверки для первого и последнего ребер!

                           Программа готова!!! Но всю равно следите за обновлениями на гитхаб [GitHub Profile](https://github.com/DebugDestroy)
*/
#include "command/Control.hpp"
#include "command/Interface.hpp"

int main() {
    // Все пути считаются относительно корня (текущей директории)
    Config config("bin/etc/config.conf");
    Logger loggerinterface(config.logFileNameInterface, "Interface", config.FiltrationLogLevelInterface);
    Logger loggercontrol(config.logFileNameControl, "Control", config.FiltrationLogLevelControl);
    
    Control c(config, loggercontrol);
    Interface i(config, loggerinterface, c);
    
    i.print();
    return 0;
}
