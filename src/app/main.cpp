/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2) Ошибка если уровень среза 127, алгоритм wave работает рекурсивно и при компонете размерером с поле рекурсия может переполнить память (также и для большого поля) нужно исправить рекурсию и хранить только пиксели коомпонент (структура пиксель x y из geometry) и сделать vector<pixel> components, тогда будет быстрее
3) Лит обзор добавить
4) Также при желании добавить правила проекта
5) Повернуть колеса в  condition на 45 градусов
6) Разные способы подключения старта/финиша к графу
7) Все пути на одной карте визуально
8) Читаемость кода улучшить + оптимизация
9) Сделать базу данных (state) отдельно а не в контрол (еще возможно стоит поле задать как vector<Pole> а не через умный указатель) и добавить вместо копии поля binary map
10) Добавить функцию grid которая создает сетку из поля, проходимость ячейки зависит от того есть ли в ней пиксель из компоненты а можно высоту хранить в ней
11) Использовать метод distance в алгоритмах поиска пути вместо std::hypot

порядок важности: 3 9 2 4 10 11 8...

Изменения:
1) Улучшил визуализацию
2) Исправил скрипт, теперь работает как задумывалось
                          
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
