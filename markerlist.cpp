#include <sstream>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include "markerlist.h"
#include "ui_markerlist.h"
#include "core/trunk/Contig.h"
#include "core/trunk/Variation.h"

MarkerList::MarkerList(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MarkerList) {
}

MarkerList::MarkerList(list<Contig*>& rgContigs, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MarkerList),
    _rgContigs(rgContigs)
{
    ui->setupUi(this);

    QValidator* intValidator = new QIntValidator();
    ui->flankSizeLineEdit->setValidator(intValidator);
    ui->markerListTableWidget->setColumnCount(3);
    ui->markerListTableWidget->setHorizontalHeaderLabels(QString("Contig;Position;Sequence").split(";"));
    ui->exportButton->setEnabled(false);
}

MarkerList::~MarkerList()
{
    delete ui;
}

void MarkerList::on_exportButton_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                    "markers.csv",
                                                    tr("*.csv"));
    QFile f( fileName );

    if( f.open( QIODevice::WriteOnly ) ){
        QTextStream ts( &f );  // #include <QtCore/QTextStream>
        QStringList strList;

        for( int r = 0; r < ui->markerListTableWidget->rowCount(); ++r ){
            strList.clear();
            for( int c = 0; c < ui->markerListTableWidget->columnCount(); ++c ){
                strList << "\""+ui->markerListTableWidget->item( r, c )->text()+"\"";
            }
            ts << strList.join( "," )+"\n";
        }
        f.close();
    }
}

void MarkerList::on_previewButton_clicked()
{
    int flankLength = ui->flankSizeLineEdit->text().toInt();
    ui->markerListTableWidget->clearContents();
    ui->markerListTableWidget->setRowCount(0);

    for(list<Contig*>::const_iterator itContig = _rgContigs.begin(); itContig != _rgContigs.end(); itContig++) {
        const vector<Variation*>& variations = (*itContig)->getVariations();
        for(vector<Variation*>::const_iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
            if((*itVar)->isReliable() && (*itVar)->getFlankLength() >= flankLength) {
                int pos = (*itVar)->getPos();
                string contigName = (*itContig)->getName();
                string sequence = (*itContig)->getSequenceIUPAC();
                stringstream ss;
                int iSeq = pos - 1;
                int iFlank = flankLength - 1;
                string lFlank(flankLength,' ');
                while(iFlank >= 0) {
                    if(sequence[iSeq] != '*') {
                        lFlank[iFlank] = sequence.at(iSeq);
                        iFlank--;
                    }
                    iSeq--;
                }

                ss << lFlank;
                ss << '[' <<  (*itVar)->getMajorAllele() << '/' <<  (*itVar)->getMinorAllele() << ']';

                iSeq = pos + 1;
                iFlank = 0;
                string rFlank(flankLength,' ');
                while(iFlank < flankLength) {
                    if(sequence[iSeq] != '*') {
                        rFlank[iFlank] = sequence.at(iSeq);
                        iFlank++;
                    }
                    iSeq++;
                }
                ss << rFlank;

                int new_row = ui->markerListTableWidget->rowCount();
                ui->markerListTableWidget->insertRow(new_row);
                QTableWidgetItem* contigNameItem = new QTableWidgetItem(QString::fromStdString(contigName));
                QTableWidgetItem* positionItem = new QTableWidgetItem(QString(tr("%1").arg(pos)));
                QTableWidgetItem* sequenceItem = new QTableWidgetItem(QString::fromStdString(ss.str()));
                sequenceItem->setFont(QFont("Courier"));

                ui->markerListTableWidget->setItem(new_row, 0, contigNameItem);
                ui->markerListTableWidget->setItem(new_row, 1, positionItem);
                ui->markerListTableWidget->setItem(new_row, 2, sequenceItem);
            }
        }
    }
    ui->markerListTableWidget->resizeColumnsToContents();
    ui->exportButton->setEnabled(true);
}



