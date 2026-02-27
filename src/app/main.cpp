/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память
3) Лит обзор добавить
4) Также при желании добавить правила проекта
5) Повернуть колеса в  condition
6) Разные способы подключения старта/финиша к графу
7) Все пути на одной карте визуально
8) Читаемость кода улучшить + оптимизация

Изменения:
1) Добавил новую команду автогенерацию гауссов
2) Сделал тестовые карты
3) Упроситил запуск через скрипт
                          
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
