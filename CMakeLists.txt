cmake_minimum_required(VERSION 3.5)

project(QAnalySeriesWASM VERSION 0.1 LANGUAGES CXX)


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "WASI")
    set(WASI TRUE)
endif()


set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 REQUIRED COMPONENTS Widgets PrintSupport)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets PrintSupport Sql Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Network)

include_directories(${CMAKE_SOURCE_DIR}/eigen-3.4.0)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        qcustomplot.cpp
        qcustomplot.h
        mainwindow.h
        mainwindow.ui
        mainwindow.cpp

        globals.cpp
        globals.h

        mynr.cpp
        mynr.h

        correlation.cpp
        correlation.h
        correlation.ui

        sedrate.cpp
        sedrate.h
        sedrate.ui

        tuning.cpp
        tuning.h
        tuning.ui

        readcsv.cpp
        readcsv.h
        readcsv.ui

        laskar.cpp
        laskar.h
        laskar.ui

        dataproc.cpp
        dataproc.h
        dataproc.ui

        insolation.cpp
        insolation.h
        insolation.ui

        ssa.cpp
        ssa.h
        ssa.ui

        etp.cpp
        etp.h
        etp.ui

        #qtwebutils.cpp
        #qtwebutils.h

        myMath.h
        math1.cpp
        spectral.h
        spectral.cpp

        spectraldialog.h
        spectraldialog.cpp
        spectraldialog.ui

        icons/QAnalySeries.icns
        Resources/exampleProject.qaproj
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(QAnalySeriesWASM
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET QAnalySeriesWASM APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(QAnalySeriesWASM SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(QAnalySeriesWASM
            ${PROJECT_SOURCES}
        )
    endif()
endif()

#qt_add_resources(QAnalySeriesWASM "Resources"
#    PREFIX "/"
#    #BIG_RESOURCES
#    FILES
#        Resources/#.db3
#)

qt_add_resources(QAnalySeriesWASM "Resources"
    PREFIX "/"
    FILES
        icons/QAnalySeries.icns
        Resources/exampleProject.qaproj
)


target_link_libraries(QAnalySeriesWASM PRIVATE ${CODE_LIB_FILES} Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Sql Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::Network)

#if (MACOSX)
  #qt_add_executable(QAnalySeriesWASM MACOSX_BUNDLE main.cpp ${app_icon_macos})

  set_target_properties(QAnalySeriesWASM PROPERTIES
      MACOSX_BUNDLE_GUI_IDENTIFIER org.palaeowiggles
      MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
      MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
      MACOSX_BUNDLE TRUE
      WIN32_EXECUTABLE TRUE
  )

  set(MACOSX_BUNDLE_ICON_FILE QAnalySeries.icns)
  set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/icons/QAnalySeries.icns")
  set_source_files_properties(${app_icon_macos} PROPERTIES
         MACOSX_PACKAGE_LOCATION "Resources")

#elseif(WASI)
  set_target_properties(QAnalySeriesWASM PROPERTIES
      QT_WASM_INITIAL_MEMORY "128MB"
  )
#endif()

#target_link_options(QAnalySeriesWASM PUBLIC
#    -sASYNCIFY
#    -Os
#    -sLLD_REPORT_UNDEFINED
#    -sERROR_ON_UNDEFINED_SYMBOLS=0
#)

#target_link_options(QAnalySeriesWASM PUBLIC -sASYNCIFY -Os -sASSERTIONS)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(QAnalySeriesWASM)
endif()
