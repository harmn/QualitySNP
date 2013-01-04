#include <list>
#include <vector>
#include "core/trunk/Variation.h"
#include "haplotypemodel.h"

bool haploTypeSortFunction(HaploType* pHap1, HaploType* pHap2) {
    return pHap1->getID() < pHap2->getID();
}

HaploTypeModel::HaploTypeModel(Contig* pContig, QObject *parent) :
    QAbstractTableModel(parent), _cRow(0), _cColumn(0), _font(QFont("Courier", 10, QFont::Bold))
{
    list<HaploType*> haploTypes = pContig->getHaploTypes();
    haploTypes.sort(haploTypeSortFunction);
    vector<Variation*> variations = pContig->getVariations();

    vector<Variation*>::iterator itVar = variations.begin();
    while(itVar != variations.end()) {
        if(!(*itVar)->isReliable()) {
            itVar = variations.erase(itVar);
        } else {
            itVar++;
        }
    }

    _cColumn = variations.size();
    _cRow = haploTypes.size();

    for(list<HaploType*>::iterator itHap = haploTypes.begin(); itHap != haploTypes.end(); itHap++) {
        QString hapID = QString::number((*itHap)->getID());
        _haploTypeIDs.push_back(hapID);
        QString nucleotideList;
        _readCounts[hapID] = (*itHap)->getReadCount();
        for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
            int pos = (*itVar)->getPos();
            nucleotideList.append((*itHap)->getNucleotideAt(pos));
        }
        _rowList.push_back(nucleotideList.split("", QString::SkipEmptyParts));
    }

    for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
        _positions.push_back((*itVar)->getPos());
        _flanks.push_back((*itVar)->getFlankLength());
    }

    _colorMap["A"] = QColor(255, 0, 0, 127);
    _colorMap["C"] = QColor(0, 255, 0, 127);
    _colorMap["G"] = QColor(255, 255, 0, 127);
    _colorMap["T"] = QColor(0, 0, 255, 127);
}

int HaploTypeModel::rowCount(const QModelIndex &parent) const
{
    return _cRow;
}

int HaploTypeModel::columnCount(const QModelIndex &parent) const
{
    return _cColumn;
}

QVariant HaploTypeModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role == Qt::DisplayRole) {
        return _rowList[index.row()][index.column()];
    } else if ( role == Qt::BackgroundRole) {
        return _colorMap.value(_rowList[index.row()][index.column()], Qt::white);
    } else if ( role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    } else if ( role == Qt::FontRole) {
        return _font;
    } else
        return QVariant();
}

QVariant HaploTypeModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole) {
         if (orientation == Qt::Horizontal){
             return QString("%1").arg(_positions[section] + 1);
         } else {
             QString hapID = _haploTypeIDs[section];
             return QString("HaploType %1 (%2 reads)").arg( hapID).arg(_readCounts[hapID]);
         }
    } else if ( role == Qt::FontRole) {
        return _font;
    } else if ( role == Qt::ToolTipRole && orientation == Qt::Horizontal) {
        return QString("%1").arg(_flanks[section]);
    }else {
        return QVariant();
    }
}

int HaploTypeModel::getPositionForColumn(int column)
{
    return _positions[column];
}
