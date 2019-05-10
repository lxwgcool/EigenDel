TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -pthread

SOURCES += \
        main.cpp \
    clsparsebam.cpp \
    ../../../ShareLibrary/clsbasealgorithm.cpp \
    ../../../ShareLibrary/clsvcf1000genome.cpp \
    clsconfig.cpp \
    ../../../ShareLibrary/clsreadconfigini.cpp \
    clsdrawimage.cpp \
    clssvdeldetect.cpp \
    clsdebug.cpp \
    clslearning.cpp \
    clscomparison.cpp \
    ../../../ShareLibrary/clsfastareader.cpp

DISTFILES += \
    Workflow_Structure \
    useful_script.txt

HEADERS += \
    clsparsebam.h \
    ../../../ShareLibrary/clsbasealgorithm.h \
    ../../../ShareLibrary/clsvcf1000genome.h \
    clsconfig.h \
    ../../../ShareLibrary/clsreadconfigini.h \
    clsparsebam.h \
    clsdrawimage.h \
    clssvdeldetect.h \
    clsdebug.h \
    clslearning.h \
    clscomparison.h \
    ../../../ShareLibrary/clsfastareader.h

unix:!macx: LIBS += -L$$PWD/../../../bamtools/install/lib/ -lbamtools

INCLUDEPATH += $$PWD/../../../bamtools/install/include/bamtools
DEPENDPATH += $$PWD/../../../bamtools/install/include/bamtools


unix:!macx: LIBS += -L$$PWD/../../../../Software/opencv/install/lib/ -lopencv_core

unix:!macx: LIBS += -L$$PWD/../../../../Software/opencv/install/lib/ -lopencv_highgui

unix:!macx: LIBS += -L$$PWD/../../../../Software/opencv/install/lib/ -lopencv_imgproc

unix:!macx: LIBS += -L$$PWD/../../../../Software/opencv/install/lib/ -lopencv_imgcodecs

INCLUDEPATH += $$PWD/../../../../Software/opencv/install/include
DEPENDPATH += $$PWD/../../../../Software/opencv/install/include
