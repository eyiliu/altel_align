find_package (Eigen3 REQUIRED NO_MODULE)

set(LIB_SRC src/TelMille.cc src/Mille.cc src/exampleUtil.cpp)

add_library(altel-telmille STATIC ${LIB_SRC})

target_include_directories(altel-telmille
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(altel-telmille
  PUBLIC Eigen3::Eigen
  PRIVATE mycommon
  )

set(LIB_PUBLIC_HEADERS include/EUTelMille.h)
set_target_properties(altel-telmille PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

install(TARGETS altel-telmille
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )


add_subdirectory(exe)
