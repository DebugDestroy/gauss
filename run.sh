#!/bin/bash

# ======================================
# Terrain Navigator - Build & Run Script
# ======================================
# Автоматически:
# 1. Собирает проект через CMake
# 2. Запускает программу
# 3. Поддерживает оба режима ввода (файл/клавиатура)
#
# Использование:
#   ./run.sh          # Запуск с клавиатуры
#   ./run.sh commands # Запуск с файлом commandsGauss.cmd
# ======================================

set -euo pipefail

# --- Конфигурация ---
BUILD_DIR="build"
COMMAND_FILE="bin/etc/commands/commandsGauss.cmd"
EXEC_NAME="terrain_navigator"

# --- Цвета для вывода ---
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# --- Функции ---
error() {
    echo -e "${RED}ERROR: $1${NC}" >&2
    exit 1
}

success() {
    echo -e "${GREEN}$1${NC}"
}

info() {
    echo -e "${YELLOW}$1${NC}"
}

# --- Проверка зависимостей ---
check_dependencies() {
    if ! command -v cmake &>/dev/null; then
        error "CMake not found! Install with:\n  Ubuntu: sudo apt install cmake\n  macOS: brew install cmake"
    fi

    if (( $(cmake --version | grep -oE '[0-9]+\.[0-9]+' | head -1 | awk '{print ($1 < 3.12)}') )); then
        error "CMake 3.12+ required"
    fi
}

# --- Сборка проекта ---
build_project() {
    info "🔧 Configuring project..."
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release
    
    CPU_CORES=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
    info "🏗  Building with $CPU_CORES cores..."
    cmake --build . --parallel "$CPU_CORES"
    
    cd ..
}

# --- Запуск программы ---
run_program() {
    local mode=$1
    
    info "🚀 Launching program..."
    
    # Получаем абсолютный путь к корню проекта
    PROJECT_ROOT="$(dirname "$(realpath "$0")")"
    
    if [[ "$mode" == "file" ]]; then
        CMD_FILE="$PROJECT_ROOT/$COMMAND_FILE"
        if [[ ! -f "$CMD_FILE" ]]; then
            error "Command file $CMD_FILE not found!"
        fi
        # Передаём абсолютные пути в программу
        echo -e "1\n$CMD_FILE\n$PROJECT_ROOT" | "$PROJECT_ROOT/$BUILD_DIR/bin/$EXEC_NAME"
    else
        echo -e "0\n$PROJECT_ROOT" | "$PROJECT_ROOT/$BUILD_DIR/bin/$EXEC_NAME"
    fi
    
    if [[ $? -eq 0 ]]; then
        success "✅ Program finished successfully"
        info "📁 Results saved in $PROJECT_ROOT/results/"
    else
        error "Program failed"
    fi
}

# --- Главный скрипт ---
main() {
    check_dependencies
    
    # Сборка если нужно
    if [[ ! -f "$BUILD_DIR/bin/$EXEC_NAME" ]]; then
        build_project
    else
        info "✅ Build already exists - skipping"
    fi
    
    # Определение режима запуска
    if [[ $# -gt 0 && "$1" == "commands" ]]; then
        run_program "file"
    else
        run_program "interactive"
    fi
}

main "$@"
