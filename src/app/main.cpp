#include "command/control.hpp"
#include "command/interface.hpp"
#include <iostream>

int main() {
    try {
        core::Config config("config/config.conf");

        core::Logger loggerinterface(
            config.logFileNameInterface,
            "Interface",
            config.FiltrationLogLevelInterface);

        core::Logger loggercontrol(
            config.logFileNameControl,
            "Control",
            config.FiltrationLogLevelControl);

        command::Control c(loggercontrol, config.seedMode, config.seed);
        command::Interface i(config, loggerinterface, c);

        i.print();
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: "
                  << e.what()
                  << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
