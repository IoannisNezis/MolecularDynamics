# Defining the executable target
add_executable(085 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(085 PUBLIC my_md_lib)
