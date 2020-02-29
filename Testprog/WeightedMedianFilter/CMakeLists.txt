project(WeightedMedianFilter)

file(GLOB sources *.cc)
add_executable(${PROJECT_NAME} ${sources})

#add_subdirectory(twoD)
add_subdirectory(wmfTest)
add_subdirectory(wmfTime)
