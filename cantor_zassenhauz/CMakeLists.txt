cmake_minimum_required(VERSION 3.14)
project(cantor_zassenhauz)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib")

SET(TARGET_H
        cantor_zassenhaus.h
        test.h)

SET(TARGET_SRC
        cantor_zassenhaus.cpp
        test.cpp)

add_executable(cantor_zassenhauz main.cpp ${TARGET_H} ${TARGET_SRC})
target_link_libraries(cantor_zassenhauz givaro gmp gmpxx)
target_compile_options(cantor_zassenhauz PUBLIC -fabi-version=6 -I/usr/local/include)
