#pragma once

#include <string>
#include <fstream>
#include "command/dispatcher_params.hpp"
#include "command/application_state.hpp"
#include "command/control.hpp"

namespace command {

class Interface {
public:
    Interface(core::Config& cfg, core::Logger& log, Control& c);

    void print();

private:
    core::Config& config;
    core::Logger& logger;
    DispatcherParams params;
    ApplicationState state;
    Control& control;
    bool initState = false;

    void showHelp();
    void processFileCommands(std::ifstream& file);
    void processKeyboardCommands();
    bool processCommand(std::istream& input, bool fromKeyboard);
};

}
