#ifndef READGROUP_H
#define READGROUP_H

#include <list>
#include <string>

class Contig;
class SeqRead;

using namespace std;

class ReadGroup
{
public:
    ReadGroup(string name, Contig* pContig);

    void addRead(SeqRead*);
    void addReads(list<SeqRead*>&);
    string toCSV(int nPos, char majorAllele, char minorAllele);

private:
    string                  _name;
    Contig*					_pContig;
    list<SeqRead*>			_reads;
};

#endif // READGROUP_H
