cmake_minimum_required(VERSION 3.16)

project(SatelliteTransponderModel VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets PrintSupport)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets PrintSupport)

add_library(qcustomplot STATIC
    qcustomplot.cpp
    qcustomplot.h
)

target_link_libraries(qcustomplot
    PUBLIC
        Qt${QT_VERSION_MAJOR}::Widgets
        Qt${QT_VERSION_MAJOR}::PrintSupport
)

set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/qcustomplot/qcustomplot.cpp
    ${CMAKE_SOURCE_DIR}/qcustomplot/qcustomplot.h
    PROPERTIES SKIP_AUTOMOC ON
)

target_include_directories(qcustomplot
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(SatelliteTransponderModel
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
        graphplotting.cpp graphplotting.h
        signalprocessing.h signalprocessing.cpp
        ofdm.h ofdm.cpp
        fft.h fft.cpp
        HelpfullStructs.h
        pamodels.h pamodels.cpp
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET SatelliteTransponderModel APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(SatelliteTransponderModel SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(SatelliteTransponderModel
            ${PROJECT_SOURCES}
        )
    endif()
endif()

# ---- FFTW ----
set(FFTW_ROOT "C:/Users/denlo/Documents/Qt/SatelliteTransponderModel/fftw")

add_library(fftw3 SHARED IMPORTED)

set_target_properties(fftw3 PROPERTIES
    IMPORTED_LOCATION "${FFTW_ROOT}/lib/libfftw3-3.dll"
    IMPORTED_IMPLIB  "${FFTW_ROOT}/lib/libfftw3-3.a"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW_ROOT}/include"
)

target_link_libraries(SatelliteTransponderModel
    PRIVATE
        Qt${QT_VERSION_MAJOR}::Widgets
        Qt${QT_VERSION_MAJOR}::PrintSupport
        qcustomplot
        fftw3
)

add_custom_command(TARGET SatelliteTransponderModel POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${FFTW_ROOT}/lib/libfftw3-3.dll"
    $<TARGET_FILE_DIR:SatelliteTransponderModel>
)
# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.SatelliteTransponderModel)
endif()
set_target_properties(SatelliteTransponderModel PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS SatelliteTransponderModel
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(SatelliteTransponderModel)
endif()
