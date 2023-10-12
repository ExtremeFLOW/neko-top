# Main cmake definition to construct examples
# Add each example directory

set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Examples directory")

add_subdirectory(neko_examples/tgv)