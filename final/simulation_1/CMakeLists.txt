# Defining the executable target
add_executable(simulation_1_Timestep find_timestep.cpp)
add_executable(simulation_1_LatticeConstant find_lattice_constant.cpp)
add_executable(simulation_1_EAMCutoff find_EAMcutoff.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(simulation_1_LatticeConstant PUBLIC my_md_lib)
target_link_libraries(simulation_1_Timestep PUBLIC my_md_lib)
target_link_libraries(simulation_1_EAMCutoff PUBLIC my_md_lib)

