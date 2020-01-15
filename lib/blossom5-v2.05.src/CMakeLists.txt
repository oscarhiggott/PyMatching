file(GLOB_RECURSE BLOSSOM_SOURCES "*.cpp")

add_library(Blossom5 ${BLOSSOM_SOURCES})
if(UNIX AND NOT APPLE)
    target_link_libraries(Blossom5 rt)
endif()
