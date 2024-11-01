cmake_minimum_required(VERSION 3.20)  # Проверка версии CMake.

project(hello_grid)  # Название проекта

# Указываем папку, где находятся заголовочные файлы
include_directories(include)

set(SOURCE_EXE main.cpp)  # Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB grid.cpp)  # Установка переменной со списком исходников для библиотеки

# Создание статической библиотеки
add_library(grid STATIC src/${SOURCE_LIB}) 

# Создание исполняемого файла
add_executable(main ${SOURCE_EXE})  

# Линковка исполняемого файла с библиотекой
target_link_libraries(main grid)  