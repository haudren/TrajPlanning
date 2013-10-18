set(SOURCES ObsPen.cpp SmoothnessMatrix.cpp Cost.cpp)
set(HEADERS ObsPen.h SmoothnessMatrix.h Cost.h)

include_directories(${Boost_INCLUDE_DIRS})

add_library(TrajPlanning SHARED ${SOURCES} ${HEADERS})
target_link_libraries(TrajPlanning ${Boost_LIBRARIES})

PKG_CONFIG_USE_DEPENDENCY(TrajPlanning SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(TrajPlanning RBDyn)

set(INSTALL_PATH include/TrajPlanning)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})

