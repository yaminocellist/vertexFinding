cmake_minimum_required(VERSION 3.10...3.22)
project(Para_Project)

find_package(nlohmann_json 3.11.3 QUIET REQUIRED)
find_package(ROOT REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(Para_dPhi para_dPhi.cpp)
target_link_libraries(Para_dPhi PUBLIC ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_include_directories(Para_dPhi PUBLIC ${ROOT_INCLUDE_DIRS})
target_compile_options(Para_dPhi PUBLIC
    -std=c++17
    -m64
    -pthread
    -stdlib=libc++
)

add_executable(Para_dEta para_dEta.cpp)
target_link_libraries(Para_dEta PUBLIC ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_include_directories(Para_dEta PUBLIC ${ROOT_INCLUDE_DIRS})
target_compile_options(Para_dEta PUBLIC
    -std=c++17
    -m64
    -pthread
    -stdlib=libc++
)

find_package(ROOT REQUIRED COMPONENTS Spectrum)
add_executable(para_INTT_zFinding para_INTT_zFinding.cpp)
target_link_libraries(para_INTT_zFinding PUBLIC ${ROOT_LIBRARIES} Spectrum OpenMP::OpenMP_CXX)
target_include_directories(para_INTT_zFinding PUBLIC ${ROOT_INCLUDE_DIRS})
target_compile_options(para_INTT_zFinding PUBLIC
    -std=c++17
    -m64
    -pthread
    -stdlib=libc++
)

add_executable(para_angulars_AI para_angulars_AI.cpp)
target_link_libraries(para_angulars_AI PUBLIC ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_include_directories(para_angulars_AI PUBLIC ${ROOT_INCLUDE_DIRS})
target_compile_options(para_angulars_AI PUBLIC
    -std=c++17
    -m64
    -pthread
    -stdlib=libc++
)
