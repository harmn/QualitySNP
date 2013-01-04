#ifndef READGROUPMODEL_H
#define READGROUPMODEL_H

#include <QStringList>
#include <QAbstractTableModel>
#include <QFont>
#include <QColor>
#include "core/trunk/Contig.h"

class ReadGroupModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit ReadGroupModel(Contig* pContig, QObject *parent = 0);
    
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

signals:

public slots:

private:
    QList<QStringList>      _rowList;
    QList<QStringList>      _rowDetailsList;
    int                     _cRow;
    int                     _cColumn;
    QMap<QString, QColor>   _colorMap;
    QFont                   _font;
    QList<int>              _positions;
    QList<QString>          _rowLabels;
    QList<QString>          _rowLabelTooltips;
    QList<int>              _readCounts;
};

#endif // READGROUPMODEL_H
