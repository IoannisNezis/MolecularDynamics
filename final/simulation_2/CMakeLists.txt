# Defining the executable target
add_executable(simulation_2 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(simulation_2 PUBLIC my_md_lib)
