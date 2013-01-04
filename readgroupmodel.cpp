#include <list>
#include <vector>
#include "core/trunk/Variation.h"
#include "core/trunk/SeqRead.h"
#include "readgroupmodel.h"

int compareRows(QStringList l1, QStringList l2) {
    int minSize = (l1.size() < l2.size()) ? l1.size() : l2.size();

    for(int iList = 0; iList != minSize; iList++) {
        if(l1[iList] != " " && l2[iList] != " " && l1[iList] != l2[iList]) {
            return 1;
        }
    }

    return 0;
}

ReadGroupModel::ReadGroupModel(Contig* pContig, QObject *parent) :
    QAbstractTableModel(parent), _font(QFont("Courier", 10, QFont::Bold))
{
    list<HaploType*> haploTypes = pContig->getHaploTypes();
    vector<Variation*> variations = pContig->getVariations();
    list<SeqRead*> reads = pContig->getReads();

    // remove all non-reliable SNPs
    vector<Variation*>::iterator itVar = variations.begin();
    while(itVar != variations.end()) {
        if(!(*itVar)->isReliable()) {
            itVar = variations.erase(itVar);
        } else {
            itVar++;
        }
    }

    _cColumn = variations.size();

    map<string, list<SeqRead*> > groupMap;

    for(list<SeqRead*>::iterator itRead = reads.begin(); itRead != reads.end(); itRead++) {
        string group = (*itRead)->getGroup();
        if(!group.empty()) {
            groupMap[group].push_back(*itRead);
        }
    }

    for(map<string, list<SeqRead*> >::iterator itGroup = groupMap.begin(); itGroup != groupMap.end(); itGroup++) {
        QStringList nucleotideList;
        QStringList detailList;
        QString groupID = QString::fromStdString((*itGroup).first);

        int cReads = (*itGroup).second.size();
        for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
            int pos = (*itVar)->getPos();

            map<char, int> nucCount;
            for(list<SeqRead*>::iterator itRead = (*itGroup).second.begin(); itRead != (*itGroup).second.end(); itRead++) {
                char nuc = (*itRead)->getNucleotideAt(pos);
                nucCount[nuc]++;
            }
            int cMajorAllele = 0;
            int cMinorAllele = 0;
            char majorAllele = ' ';
            char minorAllele = ' ';
            for(map<char, int>::iterator itNuc = nucCount.begin(); itNuc != nucCount.end(); itNuc++) {
                if((*itNuc).first == ' ') {

                } else if((*itNuc).second >= cMajorAllele) {
                    minorAllele = majorAllele;
                    majorAllele = (*itNuc).first;
                    cMinorAllele = cMajorAllele;
                    cMajorAllele = (*itNuc).second;
                } else if ((*itNuc).second >= cMinorAllele) {
                    minorAllele = (*itNuc).first;
                    cMinorAllele = (*itNuc).second;
                }
            }

            nucleotideList.push_back(QString(tr("%1")).arg(majorAllele));

            if(cMinorAllele > 0) {
                detailList.push_back(QString(tr("%1/%2 (%3/%4)")).arg(majorAllele).arg(minorAllele).arg(cMajorAllele).arg(cMinorAllele));
            } else if(cMajorAllele > 0) {
                detailList.push_back(QString(tr("%1 (%2)")).arg(majorAllele).arg(cMajorAllele));
            } else {
                detailList.push_back(QString(tr("-")));
            }
        }

        bool bFound = false;
//        int iRow = 0;
//        while(iRow != _rowList.size() && !bFound) {
//            if(compareRows(_rowList[iRow], nucleotideList) == 0) {
//                _rowLabelTooltips[iRow].append(QString(tr(", %1")).arg(groupID));
//                _rowLabels[iRow] = "multiple";
//                _readCounts[iRow] += cReads;
//                 bFound = true;
//            }
//            iRow++;
//        }

        if(!bFound) {
            _rowLabelTooltips.push_back(groupID);
            _rowLabels.push_back(groupID);
            _rowList.push_back(nucleotideList);
            _rowDetailsList.push_back(detailList);
            _readCounts.push_back(cReads);
        }
    }


    _cRow = _rowList.size();

    for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
        _positions.push_back((*itVar)->getPos());
    }

    _colorMap["A"] = QColor(255, 0, 0, 127);
    _colorMap["C"] = QColor(0, 255, 0, 127);
    _colorMap["G"] = QColor(255, 255, 0, 127);
    _colorMap["T"] = QColor(0, 0, 255, 127);
}

int ReadGroupModel::rowCount(const QModelIndex &parent) const
{
    return _cRow;
}

int ReadGroupModel::columnCount(const QModelIndex &parent) const
{
    return _cColumn;
}

QVariant ReadGroupModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role == Qt::DisplayRole) {
        return _rowDetailsList[index.row()][index.column()];
    } else if ( role == Qt::BackgroundRole) {
        return _colorMap.value(_rowList[index.row()][index.column()], Qt::white);
    } else if ( role == Qt::TextAlignmentRole) {
        return Qt::AlignCenter;
    } else if ( role == Qt::FontRole) {
        return _font;
    } else
        return QVariant();
}

QVariant ReadGroupModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole) {
         if (orientation == Qt::Horizontal){
             return QString("%1").arg(_positions[section] + 1);
         } else {
             return QString("%1 (%2 reads)").arg(_rowLabels[section]).arg(_readCounts[section]);
         }
    } else if ( role == Qt::FontRole) {
        return _font;
    } else if ( role == Qt::ToolTipRole && orientation == Qt::Vertical) {
        return _rowLabelTooltips[section];
    }else {
        return QVariant();
    }
}
