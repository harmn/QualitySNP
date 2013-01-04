#ifndef HAPLOTYPEMODEL_H
#define HAPLOTYPEMODEL_H

#include <QStringList>
#include <QAbstractTableModel>
#include <QFont>
#include <QColor>
#include "core/trunk/Contig.h"

class HaploTypeModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit HaploTypeModel(Contig* pContig, QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

    int getPositionForColumn(int);

signals:
    
public slots:

private:
    QList<QStringList>      _rowList;
    int                     _cRow;
    int                     _cColumn;
    QMap<QString, QColor>   _colorMap;
    QFont                   _font;
    QList<int>              _positions;
    QList<int>              _flanks;
    QList<QString>          _haploTypeIDs;
    QMap<QString, int>      _readCounts;
};

#endif // HAPLOTYPEMODEL_H
