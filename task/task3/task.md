# Задание, неделя 1

## Цели

- Научится работать со сборкой проектов
- Понять различия типов линковки библиотек
- Попробовать подключить сторонню библиотеку для работы

## Сборка проектов

Сборка проектов на С++ осуществляется с помощь утилиты [CMake](https://github.com/Kitware/CMake/releases/tag/v3.26.4).
CMake не лучший инструмент для сборки проектов, но самый распространенный.

В дополнение к CMake нужно установить Make ([MSYS2](https://packages.msys2.org/package/make)
или [альтернатива](https://stackoverflow.com/questions/32127524/how-to-install-and-use-make-in-windows))
или [ninja](https://github.com/ninja-build/ninja/releases/tag/v1.12.1).

## Задача

Требуется создать библиотеку по работе с тригонометрическими функциями. Предлагается воспользоваться
разложением функции в [ряд Тейлора](https://ru.wikipedia.org/wiki/%D0%A2%D1%80%D0%B8%D0%B3%D0%BE%D0%BD%D0%BE%D0%BC%D0%B5%D1%82%D1%80%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B5_%D1%84%D1%83%D0%BD%D0%BA%D1%86%D0%B8%D0%B8)
для реализации функций.

Тригонометрическая функция принимает на вход переменную `double` и должна возвращать структуру
`AnswerHandler`. Структура содержит два поля: значение `x` и значение функции в точке `x`.

Структура проекта следующая:
week01/
 ├─project/
 ├───include/
 ├─────ans_handler.hpp    <- Contains definition of AnswerHandler
 ├─────my_math.hpp        <- Contains definition of functions
 ├───src/
 ├─────ans_handler.cpp    <- Contains realization of AnswerHandler (possible empty)
 ├─────my_math.cpp        <- Contains realization of functions
 ├───CMakeLists.txt          that are defined in my_math.hpp
 └───main.cpp             <- include "my_math.hpp" and use functions

Шаблоны для реализации указаны ниже:


#### CMakeLists.txt

# Please, replace ??? and %PROJECTNAME% with proper things

cmake_minimum_required(VERSION 3.10)
# ^ Minimal of version of CMake you need to begin the build

set(CMAKE_C_COMPILER gcc)
# ^ Variable of path to C compiler
set(CMAKE_CXX_COMPILER g++)
# ^ Variable of path to C++ compiler
project(%YOURPROJECTNAME%Project)
# ^ Name of project

message(INFO:\ ${%YOURPROJECTNAME%_SOURCE_DIR})
# ^ Log message
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${%YOURPROJECTNAME%_SOURCE_DIR}/bin)
# ^ Directory where store build of project

include_directories(???/)
# ^ Directory with headers to be included (like -I option of compiler)

add_executable(%YOURPROJECTNAME% main.cpp src/???) 
# ^ Create executable with sources

#### ans_handler.hpp
// Should be something included here?
// Should be added macro variable?
struct AnswerHandler {
  ...;
};

#### ans_handler.cpp
// Should be something included here?

#### my_math.hpp
// Should be something included here?
// Should be added macro variable?
??? my_sin(double x);

#### my_math.cpp
// Should be something included here?
??? my_sin(double x) {
  ...;
}
