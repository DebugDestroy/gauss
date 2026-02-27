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
#   ./run.sh filename # Запуск с файлом filename
# ======================================


set -euo pipefail

# --- Конфигурация ---
BUILD_DIR="build"
EXEC_PATH="$BUILD_DIR/bin/terrain_navigator"

# --- Цвета ---
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

error() { echo -e "${RED}ERROR: $1${NC}" >&2; exit 1; }
info() { echo -e "${YELLOW}$1${NC}"; }

# --- Определение корня ---
get_project_root() {
    local script_dir="$(dirname "$(realpath "$0")")"
    while [[ "$script_dir" != "/" ]]; do
        [[ -d "$script_dir/config" ]] && { echo "$script_dir"; return; }
        script_dir="$(dirname "$script_dir")"
    done
    error "Корень проекта не найден!"
}

PROJECT_ROOT="$(get_project_root)"
info "📁 Корень проекта: $PROJECT_ROOT"

# --- Сборка ---
check_and_build() {
    if ! command -v cmake &>/dev/null; then
        error "CMake не установлен. Установите:\n  Ubuntu: sudo apt install cmake\n  macOS: brew install cmake"
    fi

    if [[ ! -f "$PROJECT_ROOT/$EXEC_PATH" ]]; then
        info "🔧 Сборка проекта..."
        mkdir -p "$PROJECT_ROOT/$BUILD_DIR"
        cd "$PROJECT_ROOT/$BUILD_DIR"
        cmake .. -DCMAKE_BUILD_TYPE=Release
        cmake --build . --parallel $(nproc --all 2>/dev/null || echo 4)
        cd "$PROJECT_ROOT"
    else
        info "✅ Программа уже собрана"
    fi
}

# --- Запуск ---
run_program() {
    info "🚀 Запуск программы..."
    local binary="$PROJECT_ROOT/$EXEC_PATH"

    if [[ $# -gt 0 ]]; then
        local user_file="$1"

        # Если путь относительный — делаем его от корня проекта
        if [[ ! "$user_file" = /* ]]; then
            user_file="$PROJECT_ROOT/$user_file"
        fi

        [[ ! -f "$user_file" ]] && error "Файл команд не найден: $user_file"

        info "📄 Используется файл команд: $user_file"
        echo -e "1\n$user_file" | "$binary"
    else
        info "⌨️  Режим ввода с клавиатуры"
        "$binary"
    fi
}

main() {
    check_and_build
    run_program "$@"
}

main "$@"
