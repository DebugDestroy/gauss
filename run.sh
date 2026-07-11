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
#
# Параметры:
# BUILD_TYPE="Release" или "Debug"  # Режим сборки CMake
# ======================================

set -euo pipefail

# ---------- Параметры ----------
BUILD_TYPE="Debug"

# ---------- Конфигурация ----------
BUILD_DIR="build"
EXEC_PATH="$BUILD_DIR/bin/terrain_navigator"

# ---------- Цвета ----------
RED='\033[1;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

error() { echo -e "${RED}ERROR: $1${NC}" >&2; exit 1; }
info()  { echo -e "${YELLOW}$1${NC}"; }

# ---------- Поиск корня проекта ----------
get_project_root() {
    local dir
    dir="$(dirname "$(realpath "$0")")"

    while [[ "$dir" != "/" ]]; do
        [[ -d "$dir/config" ]] && {
            echo "$dir"
            return
        }
        dir="$(dirname "$dir")"
    done

    error "Корень проекта не найден."
}

PROJECT_ROOT="$(get_project_root)"

# ---------- Поиск файла команд ----------
resolve_command_file() {
    local file="$1"

    if [[ "$file" = /* ]]; then
        [[ -f "$file" ]] || error "Файл не найден: $file"
        echo "$file"
        return
    fi

    if [[ -f "$PROJECT_ROOT/$file" ]]; then
        echo "$PROJECT_ROOT/$file"
        return
    fi

    if [[ -f "$PROJECT_ROOT/config/commands/$file" ]]; then
        echo "$PROJECT_ROOT/config/commands/$file"
        return
    fi

    error "Файл не найден: $file"
}

# ---------- Сборка ----------
check_and_build() {

    command -v cmake >/dev/null ||
        error "CMake не установлен."

    info "📁 Корень проекта: $PROJECT_ROOT"
    info "🔧 Сборка ($BUILD_TYPE)..."

    cmake \
        -S "$PROJECT_ROOT" \
        -B "$PROJECT_ROOT/$BUILD_DIR" \
        -DCMAKE_BUILD_TYPE="$BUILD_TYPE"

    cmake \
        --build "$PROJECT_ROOT/$BUILD_DIR" \
        --parallel "$(nproc --all 2>/dev/null || echo 4)"
}

# ---------- Запуск ----------
run_program() {

    local binary="$PROJECT_ROOT/$EXEC_PATH"

    info "🚀 Запуск программы"

    if [[ $# -eq 0 ]]; then
        info "⌨️  Режим ввода с клавиатуры"
        "$binary"
        return
    fi

    local command_file
    command_file="$(resolve_command_file "$1")"

    info "📄 Файл команд: $command_file"

    echo -e "1\n$command_file" | "$binary"
}

main() {
    check_and_build
    run_program "$@"
}

main "$@"
