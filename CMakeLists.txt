cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(pixel LANGUAGES CXX)

# ROOT 6.40 用 C++20 编译，整个工程跟着用它
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Geant4 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/geant4/include ${ROOT_INCLUDE_DIRS})

add_executable(pixel
  geant4/src/main.cc
  geant4/src/DetectorConstruction.cc
  geant4/src/PrimaryGeneratorAction.cc
  geant4/src/RunAction.cc
  geant4/src/EventAction.cc
  geant4/src/SteppingAction.cc
  geant4/src/ActionInitialization.cc
)

target_link_libraries(pixel PRIVATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

target_compile_options(pixel PRIVATE -O2 -Wall -Wextra -Wno-unused-parameter)

# -----------------------------------------------------------------------------
#  digitize：hits(Geant4) -> GasDigital -> digi.root
#
#  只依赖 ROOT。若还想用 Garfield++ 提供气体参数（W 等），加
#      -DUSE_GARFIELD=ON -DGARFIELD_ROOT=<Garfield++ 安装前缀>
#  找不到就自动退回自带的 gauss 后端，不影响编译。
# -----------------------------------------------------------------------------
option(USE_GARFIELD "用 Garfield++ 的 MediumMagboltz 提供气体参数" OFF)
set(GARFIELD_ROOT "" CACHE PATH "Garfield++ 安装前缀")

add_executable(digitize
  gasdigital/app/digitize.cc
  gasdigital/src/GasDigital.cc
)
target_include_directories(digitize PRIVATE
  ${PROJECT_SOURCE_DIR}/gasdigital/include
  ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(digitize PRIVATE ${ROOT_LIBRARIES})
target_compile_options(digitize PRIVATE -O2 -Wall -Wextra -Wno-unused-parameter)

if(USE_GARFIELD)
  find_path(GARFIELD_INCLUDE_DIR Garfield/MediumMagboltz.hh
            HINTS ${GARFIELD_ROOT}/include NO_DEFAULT_PATH)
  find_library(GARFIELD_LIBRARY NAMES Garfield
               HINTS ${GARFIELD_ROOT}/lib NO_DEFAULT_PATH)

  if(GARFIELD_INCLUDE_DIR AND GARFIELD_LIBRARY)
    target_include_directories(digitize PRIVATE ${GARFIELD_INCLUDE_DIR})
    target_link_libraries(digitize PRIVATE ${GARFIELD_LIBRARY})
    target_compile_definitions(digitize PRIVATE
      HAVE_GARFIELD=1
      GARFIELD_HEED_DB="${GARFIELD_ROOT}/Heed/heed++/database"
    )
    # Heed 的数据文件和 libGarfield.so 运行时都得找得到
    set_target_properties(digitize PROPERTIES
      BUILD_RPATH "${GARFIELD_ROOT}/lib"
      INSTALL_RPATH "${GARFIELD_ROOT}/lib")

    # Garfield 里 Magboltz 是 Fortran，静态库要补 libgfortran
    find_library(GFORTRAN_LIBRARY NAMES gfortran
                 HINTS ${CMAKE_PREFIX_PATH}/lib /usr/lib/x86_64-linux-gnu)
    if(GFORTRAN_LIBRARY)
      target_link_libraries(digitize PRIVATE ${GFORTRAN_LIBRARY})
    endif()

    message(STATUS "digitize: Garfield++ 后端已启用 (${GARFIELD_LIBRARY})")
  else()
    message(WARNING
      "USE_GARFIELD=ON 但没找到 Garfield++（GARFIELD_ROOT=${GARFIELD_ROOT}），"
      "digitize 将只用自带的 gauss 后端")
  endif()
endif()
