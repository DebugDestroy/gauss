/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2)Возможно новые требования
3) Оптимизация

6)В алгоритме пути сначала проверки для ребер а уже потом алгоритм
7)Считать углы не между пискселями а между центром и пикселем пути удалееным на радиус и также для двух перпендикулярных направлений

8)Написать методичку (пикасо для обрезки фото) + псевдокод к алгоритму пути


Изменения: компилируется в гз +  параметр в конфиг для help + автодополнение команд при вводе с клавиатуры + Убрал список пути + cmake + супер файловая структура!!!

                           Программа готова!!! Но всю равно следите за обновлениями на гитхаб [GitHub Profile](https://github.com/DebugDestroy)
*/
#include "command/Control.hpp" 
#include "command/Interface.hpp"

int main(int argc, char* argv[]) {
    // Берём корень проекта из аргументов (argv[2])
    std::string root_path = (argc > 2) ? argv[2] : ".";
    
    // Используем в путях
    Config config(root_path + "/bin/etc/config.conf");
    Logger loggerinterface(root_path + "/" + config.logFileNameInterface, "Interface", config.FiltrationLogLevelInterface);
    Logger loggercontrol(root_path + "/" +  config.logFileNameControl, "Control", config.FiltrationLogLevelControl);
    // Создаем интерфейс
    Control c(config, loggercontrol);
    Interface i(config, loggerinterface, c);
    
    // Вызываем метод print() интерфейса
    i.print();

    return 0;
}
