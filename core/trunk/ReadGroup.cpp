#include "Configuration.h"
#include "SeqRead.h"
#include "Variation.h"
#include "Contig.h"
#include "ReadGroup.h"

ReadGroup::ReadGroup(string name, Contig* pContig): _name(name), _pContig(pContig)
{
}

void ReadGroup::addRead(SeqRead* pRead)
{
    _reads.push_back(pRead);
}

void ReadGroup::addReads(list<SeqRead *> & reads)
{
    _reads.insert(_reads.end(), reads.begin(), reads.end());
}

string ReadGroup::toCSV(int nPos, char majorAllele, char minorAllele)
{
    stringstream csv;

    int cMajorAllele = 0;
    int cMinorAllele = 0;
    list<SeqRead*>::iterator itReads;
    for (itReads=_reads.begin(); itReads != _reads.end(); itReads++) {
        char nuc = (*itReads)->getNucleotideAt(nPos);
        if(nuc == majorAllele) {
            cMajorAllele++;
        } else if (nuc == minorAllele) {
            cMinorAllele++;
        }
    }

    csv << majorAllele << ":" << cMajorAllele << "/" << minorAllele << ":" << cMinorAllele;

    return csv.str();
}
