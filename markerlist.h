#ifndef MARKERLIST_H
#define MARKERLIST_H

#include <QDialog>
#include <list>
#include "core/trunk/Contig.h"

namespace Ui {
class MarkerList;
}

class MarkerList : public QDialog
{
    Q_OBJECT
    
public:
    explicit MarkerList(QWidget *parent);
    MarkerList(list<Contig*>&, QWidget *parent);
    ~MarkerList();

private slots:
    void on_previewButton_clicked();
    void on_exportButton_clicked();

private:
    Ui::MarkerList *ui;
    list<Contig*>   _rgContigs;
};

#endif // MARKERLIST_H
