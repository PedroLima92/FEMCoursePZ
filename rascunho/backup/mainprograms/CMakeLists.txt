
add_executable(FEMCoursePZ main.cpp ${SOURCES})

if(USING_PZ)
 target_link_libraries(FEMCoursePZ pz)

endif(USING_PZ)

if(USING_ARMADILLO)
 target_link_libraries(FEMCoursePZ ${ARMADILLO_LIBRARIES})

endif(USING_ARMADILLO)

if(USING_EIGEN3)
 target_link_libraries(FEMCoursePZ ${EIGEN_LIBRARIES})

endif(USING_EIGEN3)
