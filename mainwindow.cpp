#include <QDockWidget>
#include <QPushButton>
#include <QFileDialog>
#include <QTableWidget>
#include <QHeaderView>
#include <QDesktopWidget>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QApplication>
#include <QFormLayout>
#include <QLineEdit>
#include <QSpinBox>
#include <QLabel>
#include <QUrl>
#include <QDesktopServices>
#include <sstream>
#include "core/trunk/Variation.h"
#include "haplotypemodel.h"
#include "readgroupmodel.h"
#include "mainwindow.h"
#include "rundialog.h"
#include "markerlist.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    _contigListModel = NULL;
    _haploTypeView  = NULL;
    _readGroupView  = NULL;
    _contigList     = NULL;
    _alignmentDock  = NULL;
    _haplotypeDock  = NULL;
    _readGroupDock  = NULL;

            resize(900,600);
    createMenuBar();

    createContigListDock();
    fillContigListDock();

    createAlignmentDock();
    moveToCenterOfScreen();
    setWindowTitle("QualitySNPng");
}

void MainWindow::moveToCenterOfScreen()
{
    QDesktopWidget *desktop = QApplication::desktop();
    int screenWidth, width;
    int screenHeight, height;
    int x, y;
    QSize windowSize;

    screenWidth = desktop->width(); // get width of screen
    screenHeight = desktop->height(); // get height of screen

    windowSize = size(); // size of our application window
    width = windowSize.width();
    height = windowSize.height();

    // little computations
    x = (screenWidth - width) / 2;
    y = (screenHeight - height) / 2;
    y -= 50;

    // move window to desired coordinates
    move ( x, y );

}

void MainWindow::showContigDetails(QModelIndex clickedIndex)
{
    QModelIndex index = _contigList->model()->index(clickedIndex.row(),0);
    QString qcontigName = index.data().toString();
    string contigName = qcontigName.toStdString();

    Contig* pContig = _reader.getContig(contigName);

    if(pContig == NULL) {
        return;
    }

#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    _haplotypeDock->setWindowTitle(qcontigName + " - haplotypes");      
    _readGroupDock->setWindowTitle(qcontigName + " - groups");
    _alignmentDock->setWindowTitle(qcontigName + " - reads");


    _alignmentPicture = new AlignmentPicture(pContig, this);
    _alignmentDock->setWidget(_alignmentPicture);
    _alignmentPicture->zoomToFit();

    showHaploTypeView(pContig);
    bool bReadGroups = showReadGroupView(pContig);
    _readGroupDock->setVisible(bReadGroups);

#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
}

void MainWindow::showHaploTypeView(Contig *pContig)
{
    if(_haploTypeView != NULL) {
        delete _haploTypeView;
    }
    _haploTypeView = new QTableView(_haplotypeDock);
    _haploTypeView->setSelectionMode(QAbstractItemView::NoSelection);
    HaploTypeModel* model = new HaploTypeModel(pContig, this);
    _haploTypeView->setModel(model);
    _haploTypeView->resizeColumnsToContents();
    _haplotypeDock->setWidget(_haploTypeView);
    connect(_haploTypeView, SIGNAL(clicked(QModelIndex)), this, SLOT(scrollToAlignment(QModelIndex)));
}

bool MainWindow::showReadGroupView(Contig *pContig)
{
    ReadGroupModel* model = new ReadGroupModel(pContig, this);
    if(model->rowCount() == 0) {
        return false;
    }

    if(_readGroupView != NULL) {
        delete _readGroupView;
    }
    _readGroupView = new QTableView(_readGroupDock);
    _readGroupView->setSelectionMode(QAbstractItemView::NoSelection); 
    _readGroupView->setModel(model);
    _readGroupView->resizeColumnsToContents();
    _readGroupDock->setWidget(_readGroupView);
    connect(_readGroupView, SIGNAL(clicked(QModelIndex)), this, SLOT(scrollToAlignment(QModelIndex)));

    return true;
}

void MainWindow::runQSNP()
{
    RunDialog run(this);
    run.exec();
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
    fillContigListDock();
}

void MainWindow::openFileDialog()
{
    QString dirName = QFileDialog::getExistingDirectory(this, "Choose a directory");

    Configuration* pConfig = Configuration::getConfig(true);
    pConfig->readFile("config.cfg", dirName.toStdString());
    pConfig->setString("inputDirectory", dirName.toStdString());
    fillContigListDock();
}

void MainWindow::createContigListDock()
{
    QDockWidget*  contigDock = new QDockWidget(tr("Contigs"), this);
    contigDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    contigDock->setAllowedAreas(Qt::LeftDockWidgetArea);
    contigDock->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    QWidget* filterWidget = new QWidget();
    filterWidget->setWindowTitle(tr("Filter"));

    QGridLayout *gridLayout = new QGridLayout;
    QLabel* hapLabel = new QLabel(tr("Haplotypes"));
    hapLabel->setToolTip(tr("max haplotypes per SNP"));
    QLabel* readLabel = new QLabel(tr("Reads"));
    QLabel* SNPLabel = new QLabel(tr("SNPs"));
    _hapMinLineEdit = new QSpinBox();
    _hapMaxLineEdit = new QSpinBox();
    _readMinLineEdit = new QSpinBox();
    _readMaxLineEdit = new QSpinBox();
    _SNPMinLineEdit = new QSpinBox();
    _SNPMaxLineEdit = new QSpinBox();

    _contigCountLabel = new QLabel(tr("0"));

    QLabel* contigSearchLabel = new QLabel(tr("Contig name:"));
    _ContigSearchLineEdit = new QLineEdit();

    QLabel* filterLabel = new QLabel(tr("Filter"));
    QLabel* minLabel = new QLabel(tr("Min #"));
    QLabel* maxLabel = new QLabel(tr("Max #"));
    gridLayout->addWidget(filterLabel, 0, 0);
    gridLayout->addWidget(minLabel, 0, 1);
    gridLayout->addWidget(maxLabel, 0, 2);

    gridLayout->addWidget(hapLabel, 1, 0);
    gridLayout->addWidget(_hapMinLineEdit, 1, 1);
    gridLayout->addWidget(_hapMaxLineEdit, 1, 2);

    gridLayout->addWidget(readLabel, 2, 0);
    gridLayout->addWidget(_readMinLineEdit, 2, 1);
    gridLayout->addWidget(_readMaxLineEdit, 2, 2);

    gridLayout->addWidget(SNPLabel, 3, 0);
    gridLayout->addWidget(_SNPMinLineEdit, 3, 1);
    gridLayout->addWidget(_SNPMaxLineEdit, 3, 2);
    gridLayout->addWidget(contigSearchLabel, 4, 0);
    gridLayout->addWidget(_ContigSearchLineEdit, 4, 1, 1, 2);
    gridLayout->addWidget(_contigCountLabel, 4, 3, 1, 1, Qt::AlignRight);

    filterWidget->setLayout(gridLayout);

    _contigList = new QTableView(contigDock);
    _contigList->setEditTriggers(NULL);
    _contigList->verticalHeader()->hide();

    _contigList->setSelectionBehavior(QAbstractItemView::SelectRows);
    _contigList->setSortingEnabled(true);

    connect(_contigList, SIGNAL(clicked(QModelIndex)), this, SLOT(showContigDetails(QModelIndex)));

    QWidget* contigWidget = new QWidget();
    QVBoxLayout* boxLayout = new QVBoxLayout();

    boxLayout->addWidget(filterWidget);
    boxLayout->addWidget(_contigList);
    contigWidget->setLayout(boxLayout);
    contigDock->setWidget(contigWidget);

    addDockWidget(Qt::LeftDockWidgetArea, contigDock);
}

void MainWindow::createAlignmentDock()
{
    _haplotypeDock = new QDockWidget(tr("haplotypes"), this);
    _haplotypeDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    _haplotypeDock->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    addDockWidget(Qt::RightDockWidgetArea, _haplotypeDock);

    _readGroupDock = new QDockWidget(tr("readgroups"), this);
    _readGroupDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    _readGroupDock->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    addDockWidget(Qt::RightDockWidgetArea, _readGroupDock);

    _alignmentDock = new QDockWidget(tr("alignment"), this);
    _alignmentDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    _alignmentDock->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    addDockWidget(Qt::RightDockWidgetArea, _alignmentDock);
}

void MainWindow::createMenuBar()
{
    QAction* openDirAction = new QAction(tr("&Load saved run"), this);
    openDirAction->setStatusTip(tr("Load previous results"));
    connect(openDirAction, SIGNAL(triggered()), this, SLOT(openFileDialog()));
    QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openDirAction);

    QAction* markersAction = new QAction(tr("&Export marker list"), this);
    markersAction->setStatusTip(tr("Export marker list"));
    connect(markersAction, SIGNAL(triggered()), this, SLOT(exportMarkers()));
    fileMenu->addAction(markersAction);

    QAction* exitAction = new QAction(tr("&Quit"), this);
    exitAction->setStatusTip(tr("Quit program"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
    fileMenu->addAction(exitAction);

    QAction* runAction = new QAction(tr("&Run QualitySNPng"), this);
    runAction->setStatusTip(tr("Run QualitySNPng"));
    connect(runAction, SIGNAL(triggered()), this, SLOT(runQSNP()));
    QMenu* runMenu = menuBar()->addMenu(tr("Run"));
    runMenu->addAction(runAction);

    QAction* helpAction = new QAction(tr("&Manual"), this);
    helpAction->setStatusTip(tr("Opens the on-line manual"));
    connect(helpAction, SIGNAL(triggered()), this, SLOT(helpURL()));
    QMenu* helpMenu = menuBar()->addMenu(tr("Help"));
    helpMenu->addAction(helpAction);
}

void MainWindow::helpURL() {
    QDesktopServices::openUrl(QUrl("http://www.bioinformatics.nl/QualitySNPng/", QUrl::TolerantMode));
}

void MainWindow::fillContigListDock()
{
    _contigInfoList = _reader.getContigList();

    if(_contigListModel != NULL) {
        delete _contigListModel;
    }

    _contigListModel = new ContigListModel(_contigInfoList, this);
    _contigList->setModel(_contigListModel);

    int maxSNP = _contigListModel->getMaxSNPCount();
    int maxHap = _contigListModel->getMaxHapCount();
    int maxReads = _contigListModel->getMaxReadsCount();

    _hapMinLineEdit->setRange(0,maxHap);
    _hapMinLineEdit->setValue(0);
    _hapMaxLineEdit->setRange(0,maxHap);
    _hapMaxLineEdit->setValue(maxHap);

    _readMinLineEdit->setRange(0, maxReads);
    _readMinLineEdit->setValue(0);
    _readMaxLineEdit->setRange(0, maxReads);
    _readMaxLineEdit->setValue(maxReads);

    _SNPMinLineEdit->setRange(0, maxSNP);
    _SNPMinLineEdit->setValue(1);
    _SNPMaxLineEdit->setRange(0, maxSNP);
    _SNPMaxLineEdit->setValue(maxSNP);

    connect(_SNPMinLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_SNPMaxLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_readMinLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_readMaxLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_hapMinLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_hapMaxLineEdit, SIGNAL(valueChanged(int)), this, SLOT(filterContigList()));
    connect(_ContigSearchLineEdit, SIGNAL(textEdited(QString)), this, SLOT(filterContigList()));

    filterContigList();

    _contigList->resizeColumnsToContents();

}

void MainWindow::resetDocks() {

    if(_alignmentDock != NULL) {
        _alignmentDock->setWidget(NULL);
        _alignmentDock->setWindowTitle("alignment");
    }

    if(_haplotypeDock != NULL) {
        _haplotypeDock->setWidget(NULL);
        _haplotypeDock->setWindowTitle("haplotypes");
    }

    if(_readGroupDock != NULL) {
        _readGroupDock->setWidget(NULL);
        _readGroupDock->setWindowTitle("readgroups");
    }
}

void MainWindow::filterContigList()
{
    int minSNP = _SNPMinLineEdit->value();
    int maxSNP = _SNPMaxLineEdit->value();
    int minReads = _readMinLineEdit->value();
    int maxReads = _readMaxLineEdit->value();
    int minHap = _hapMinLineEdit->value();
    int maxHap = _hapMaxLineEdit->value();

    QString contigSearch = _ContigSearchLineEdit->text();
    int nRowCount = _contigListModel->filter(minSNP, maxSNP, minHap, maxHap, minReads, maxReads, contigSearch);

    _contigCountLabel->setText(QString("%1 contigs").arg(nRowCount));

    resetDocks();
}

void MainWindow::exportMarkers()
{
    list<Contig*> rgContigs;
    QList<QString> contigNames = _contigListModel->getVisibleContigNames();
    QList<QString>::iterator itNames;
    for(itNames = contigNames.begin(); itNames != contigNames.end(); itNames++) {
        string contigName = (*itNames).toStdString();
        Contig* pContig = _reader.getContig(contigName);
        rgContigs.push_back(pContig);
    }

    MarkerList markerList(rgContigs, this);
    markerList.exec();
}

void MainWindow::scrollToAlignment(QModelIndex index)
{
    int position = static_cast<HaploTypeModel*>(_haploTypeView->model())->getPositionForColumn(index.column());
    _alignmentPicture->scrollToPosition(position);
}



