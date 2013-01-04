#ifndef CONTIGLISTMODEL_H
#define CONTIGLISTMODEL_H

#include <QAbstractTableModel>
#include <QStringList>
#include <QFont>
#include <vector>
#include "core/trunk/CSVReader.h"

struct ContigRow {
    QString name;
    int     length;
    int     cRelSNP;
    int     maxHapSNP;
    int     cReads;
};

class ContigListModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit ContigListModel(vector<ContigInfo>& contigInfoList, QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

    int getMaxSNPCount() { return _maxSNPCount; }
    int getMaxHapCount() { return _maxHapCount; }
    int getMaxReadsCount() { return _maxReadsCount; }
    int filter(int minSNP, int maxSNP, int minHap, int maxHap, int minRead, int maxRead, QString name);
    QList<QString>  getVisibleContigNames();
    void sort ( int column, Qt::SortOrder order = Qt::AscendingOrder );

signals:
    
public slots:

private:
    QList<ContigRow*>    _visibleRowList;
    QList<ContigRow>     _rowList;
    vector<ContigInfo>&  _contigInfoList;
    int                  _maxSNPCount;
    int                  _maxHapCount;
    int                  _maxReadsCount;
    int                  _cRow;
    int                  _cColumn;
    QList<QString>       _headers;
    QFont                _font;
    
};

#endif // CONTIGLISTMODEL_H
