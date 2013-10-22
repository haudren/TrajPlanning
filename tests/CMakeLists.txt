ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(TrajTest TrajTest.cpp)
target_link_libraries(TrajTest ${Boost_LIBRARIES} TrajPlanning)
PKG_CONFIG_USE_DEPENDENCY(TrajTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(TrajTest RBDyn)
add_test(TrajTestUnit TrajTest)

file(COPY "empty1.map" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "obs1.map" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "obs2.map" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "obs3.map" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
