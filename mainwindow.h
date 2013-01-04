#ifndef CONTIGLISTWINDOW_H
#define CONTIGLISTWINDOW_H

#include <QMainWindow>
#include <QTableView>
#include <QTableWidget>
#include <QSpinBox>
#include <QLabel>
#include <vector>
#include "core/trunk/Contig.h"
#include "core/trunk/CSVReader.h"
#include "alignmentpicture.h"
#include "contiglistmodel.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);

private slots:
    void showContigDetails(QModelIndex);
    void openFileDialog();
    void runQSNP();
    void helpURL();
    void filterContigList();
    void exportMarkers();
    void scrollToAlignment(QModelIndex);

private:
    void createContigListDock();
    void createAlignmentDock();
    void createMenuBar();
    void fillContigListDock();
    void resetDocks();

    void moveToCenterOfScreen();
    void showAlignmentView(Contig* pContig);
    void showHaploTypeView(Contig* pContig);
    bool showReadGroupView(Contig* pContig);

    QDockWidget*    _alignmentDock;
    QDockWidget*    _haplotypeDock;
    QDockWidget*    _readGroupDock;
    QTableView*     _contigList;
    QTableView*     _haploTypeView;
    QTableView*     _readGroupView;
    ContigListModel*    _contigListModel;

    vector<ContigInfo> _contigInfoList;

    QSpinBox* _hapMinLineEdit;
    QSpinBox* _hapMaxLineEdit;
    QSpinBox* _readMinLineEdit;
    QSpinBox* _readMaxLineEdit;
    QSpinBox* _SNPMinLineEdit;
    QSpinBox* _SNPMaxLineEdit;
    QLineEdit* _ContigSearchLineEdit;
    QLabel*    _contigCountLabel;

    AlignmentPicture* _alignmentPicture;

    CSVReader   _reader;
};

#endif // CONTIGLISTWINDOW_H
