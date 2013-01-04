#-------------------------------------------------
#
# Project created by QtCreator 2012-08-05T14:51:52
#
#-------------------------------------------------

QT       += core gui

TARGET = QualitySNPng
TEMPLATE = app


SOURCES += main.cpp\
    core/trunk/Variation.cpp \
    core/trunk/SeqRead.cpp \
    core/trunk/SAMRead.cpp \
    core/trunk/SAMFile.cpp \
    core/trunk/SAMContig.cpp \
    core/trunk/Logger.cpp \
    core/trunk/HaploType.cpp \
    core/trunk/CSVWriter.cpp \
    core/trunk/ContigProvider.cpp \
    core/trunk/ContigPrinter.cpp \
    core/trunk/Contig.cpp \
    core/trunk/Configuration.cpp \
    core/trunk/ACEFile.cpp \
    core/trunk/CSVReader.cpp \
    haplotypemodel.cpp \
    rundialog.cpp \
    mainwindow.cpp \
    runqsnp.cpp \
    settingsdialog.cpp \
    markerlist.cpp \
    alignmentpicture.cpp \
    contiglistmodel.cpp \
    readgroupmodel.cpp \
    core/trunk/ReadGroup.cpp

HEADERS  += \
    core/trunk/Variation.h \
    core/trunk/SeqRead.h \
    core/trunk/SAMRead.h \
    core/trunk/SAMFile.h \
    core/trunk/SAMContig.h \
    core/trunk/QualitySNPpp.h \
    core/trunk/Logger.h \
    core/trunk/HaploType.h \
    core/trunk/CSVWriter.h \
    core/trunk/ContigProvider.h \
    core/trunk/ContigPrinter.h \
    core/trunk/ContigFile.h \
    core/trunk/Contig.h \
    core/trunk/Configuration.h \
    core/trunk/ACEFile.h \
    core/trunk/CSVReader.h \
    haplotypemodel.h \
    rundialog.h \
    mainwindow.h \
    runqsnp.h \
    settingsdialog.h \
    markerlist.h \
    alignmentpicture.h \
    contiglistmodel.h \
    readgroupmodel.h \
    core/trunk/ReadGroup.h

FORMS    += \
    rundialog.ui \
    markerlist.ui \
    settingsdialog.ui
