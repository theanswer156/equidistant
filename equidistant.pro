#-------------------------------------------------
#
# Project created by QtCreator 2024-07-22T11:06:22
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = equidistant
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11
INCLUDEPATH += D:\VSCode_Program/archheight
INCLUDEPATH += D:\VSCode_Program/segcut
SOURCES += \
        main.cpp \
        mainwindow.cpp \
        D:/VSCode_Program/archheight/archheight.cpp \
        D:\VSCode_Program/segcut/segment.cpp
HEADERS += \
        mainwindow.h \
        D:/VSCode_Program/archheight/archheight.h \
        D:\VSCode_Program/segcut/segment.h
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
