#ifndef __CSVREADER_H__
#define __CSVREADER_H__

#include "Contig.h"

using namespace std;

struct ContigInfo{
 string name;
 int potSNP;
 int hqSNP;
 int relSNP;
 double Dvalue;
 int reads;
 int haplotypes;
 int maxHapSNP;
 string sequence;
};

class CSVReader
{
public:
    CSVReader();

    Contig* getContig(string contigName);
    const vector<ContigInfo> getContigList();
    void getVariations(Contig* pContig);

private:
    void getReads(Contig* pContig, map<int, vector<SeqRead*> >&);
    void getHaploTypes(Contig* pContig, map<int, vector<SeqRead*> >&);
    void readIndices();

private:
    map<string, streampos>  _contigsIndex;
    map<string, streampos>  _haploTypesIndex;
    map<string, streampos>  _readsIndex;
    map<string, streampos>  _variationsIndex;
};

#endif // __CSVREADER_H__
