project(stereoTime)

add_definitions("-DPROFILE -DUSE_TBB")

file(GLOB sources *.cc)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} TUTools tbb)
