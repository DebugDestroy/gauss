/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Возможно новые требования
3) Оптимизация


Изменения: 
1) Благодаря тестам обнаружена ошибка при построении супер треугольника в триангуляции, плохая инкапсуляция в этом методе + ошибка в методе проверящем принадлежность точки к кругу. Все это исправлено! 
2) Добавлены доп проверки в PathFinder в метод проверки проходимости (защита от деления на ноль + защита если неправильно отработал метод поиска диаграммы Воронова)
3) Диаграмма Воронова строилась неверно при сложном поле! Заметил плохую реализацию в этом методе и исправил + убрал лишнее + убрал структуру Воронова + добавляются обрезанные точки на границе? +
 улучшено логирование
4) PathFinder использовал не ребра Воронова, а их приближенную версию. Пришлось изменить логику и явно передавать ребра Воронова!

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
