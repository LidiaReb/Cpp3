cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_COMPILER "g++")

project(Ilovecpp LANGUAGES CXX)	

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(maincraft)
target_sources(maincraft PRIVATE main.cpp)

# Возможно, у тебя есть директории с заголовочными файлами - добавим их
# include_directories(include)

# Для сборки проекта
# Чтобы построить проект, просто запустите в терминале:
# mkdir build
# cd build
# cmake ..
# make
