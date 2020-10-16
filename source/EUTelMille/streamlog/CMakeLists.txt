set(LIB_SRC src/prefix.cc src/logstream.cc)
add_library(streamlog STATIC ${LIB_SRC})

target_include_directories(streamlog
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

# set(LIB_PUBLIC_HEADERS include/EUTelMille.h)
# set_target_properties(streamlog PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

install(TARGETS streamlog
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
