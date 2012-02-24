#ifndef __CONTIGPRINTER_H__
#define __CONTIGPRINTER_H__

class ContigPrinter
{
public:
	ContigPrinter(Contig*);
	~ContigPrinter(void);
	string toString();

private:
	void printSummaryLine(stringstream& result);
	void printSNPs(stringstream& result, bool bHQOnly);
	void printHaploTypes(stringstream& result);
	void printMarkerSNPs(stringstream& result, unsigned int flankSize);
	void printAlignment(stringstream& result);
	Contig* _pContig;
};

#endif
