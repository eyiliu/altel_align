add_subdirectory(streamlog)


set(LIB_SRC src/EUTelMille.cc src/Mille.cc)
add_library(altel-eutelmille STATIC ${LIB_SRC})

target_include_directories(altel-eutelmille
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )


target_link_libraries(altel-eutelmille
  PUBLIC streamlog
  )


set(LIB_PUBLIC_HEADERS include/EUTelMille.h)
set_target_properties(altel-eutelmille PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

install(TARGETS altel-eutelmille
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
