#include "contiglistmodel.h"

ContigListModel::ContigListModel(vector<ContigInfo>& contigInfoList, QObject *parent) :
    QAbstractTableModel(parent), _contigInfoList(contigInfoList), _font(QFont("Courier", 10, QFont::Bold))
{
    _cColumn = 5;
    _cRow = 0;
    _maxSNPCount = 0;
    _maxHapCount = 0;
    _maxReadsCount = 0;

    _headers << "Name" << "Length" << "SNPs" << "Haplotypes" << "Reads";

    vector<ContigInfo>::iterator itContigs;
    for(itContigs = _contigInfoList.begin(); itContigs != _contigInfoList.end(); itContigs++) {
        QString qcontigName = QString(QString::fromUtf8((*itContigs).name.data(), (*itContigs).name.size()));

        int sequenceLength = (*itContigs).sequence.length();

        ContigRow row;
        row.name = qcontigName;
        row.length = sequenceLength;
        row.cRelSNP = (*itContigs).relSNP;
        row.maxHapSNP = (*itContigs).maxHapSNP;
        row.cReads = (*itContigs).reads;

        _maxSNPCount = (row.cRelSNP < _maxSNPCount) ? _maxSNPCount : row.cRelSNP;
        _maxHapCount = (row.maxHapSNP < _maxHapCount) ? _maxHapCount : row.maxHapSNP;
        _maxReadsCount = (row.cReads < _maxReadsCount) ? _maxReadsCount : row.cReads;

        _rowList << row;
        _visibleRowList << &_rowList.back();
        _cRow++;
    }
}

int ContigListModel::rowCount(const QModelIndex &parent) const
{
    return _cRow;
}

int ContigListModel::columnCount(const QModelIndex &parent) const
{
    return _cColumn;
}

QVariant ContigListModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role == Qt::DisplayRole) {
        ContigRow* pRow = _visibleRowList[index.row()];
        switch(index.column()) {
        case 0:
            return pRow->name;
        case 1:
            return pRow->length;
        case 2:
            return pRow->cRelSNP;
        case 3:
            return pRow->maxHapSNP;
        case 4:
            return pRow->cReads;
        default:
            return QVariant();
        }

    } else if ( role == Qt::FontRole) {
        return _font;
    } else
        return QVariant();
}

QVariant ContigListModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole) {
         if (orientation == Qt::Horizontal){
             return _headers[section];
         }
    } else if ( role == Qt::FontRole) {
        return _font;
    }

    return QVariant();
}


QList<QString> ContigListModel::getVisibleContigNames()
{
    QList<QString> contigNameList;
    QList<ContigRow*>::const_iterator itRow;
    for(itRow = _visibleRowList.constBegin(); itRow != _visibleRowList.constEnd(); itRow++) {
        contigNameList << (*itRow)->name;
    }
    return contigNameList;
}

static int s_sortColumn = -1;
static Qt::SortOrder s_order = Qt::AscendingOrder;

bool contigSortFunction(ContigRow* pRow1, ContigRow* pRow2) {
    if(s_order == Qt::DescendingOrder) {
        ContigRow* tmp = pRow1;
        pRow1 = pRow2;
        pRow2 = tmp;
    }

    switch(s_sortColumn) {
    case 0:
        return pRow1->name < pRow2->name;
    case 1:
        return pRow1->length < pRow2->length;
    case 2:
        return pRow1->cRelSNP < pRow2->cRelSNP;
    case 3:
        return pRow1->maxHapSNP < pRow2->maxHapSNP;
    case 4:
        return pRow1->cReads < pRow2->cReads;
    }

    return true;
}

int ContigListModel::filter(int minSNP, int maxSNP, int minHap, int maxHap, int minRead, int maxRead, QString name)
{
    _cRow = 0;
    QList<ContigRow>::iterator itRow;
    _visibleRowList.clear();
    for(itRow = _rowList.begin(); itRow != _rowList.end(); itRow++) {
        if( (*itRow).cRelSNP >= minSNP && (*itRow).cRelSNP <= maxSNP &&
            (*itRow).cReads >= minRead  && (*itRow).cReads <= maxRead &&
            (*itRow).maxHapSNP >= minHap && (*itRow).maxHapSNP <= maxHap &&
            (*itRow).name.contains(name, Qt::CaseInsensitive) ) {
            _visibleRowList << &(*itRow);
            _cRow++;
        }
    }

    if(s_sortColumn != -1) {
        qSort(_visibleRowList.begin(), _visibleRowList.end(), contigSortFunction);
    }

    reset();
    return _cRow;
}

void ContigListModel::sort(int column, Qt::SortOrder order)
{
    s_sortColumn = column;
    s_order = order;
    qSort(_visibleRowList.begin(), _visibleRowList.end(), contigSortFunction);
    reset();
}
