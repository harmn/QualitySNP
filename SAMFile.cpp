#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include "HaploType.h"
#include "SeqRead.h"
#include "Contig.h"
#include "Variation.h"
#include "ACEFile.h"
#include "SAMRead.h"
#include "SAMContig.h"
#include "SAMFile.h"

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

SAMFile::SAMFile(const string& contigFilename): _contigFilename(contigFilename), _pRead(NULL)
{
}

SAMFile::~SAMFile(void)
{
}


bool SAMFile::openFile() {
	if(_ifSAMFile.is_open()) {
		return true;
	}

	_ifSAMFile.open(_contigFilename.c_str(),ios::in);
	if (_ifSAMFile.fail()) {
		Logger::getLogger()->log(QSNP_ERROR, "could not open file: " + _contigFilename);
		return false;
	}

    return true;
}

bool SAMFile::isValid() {
	if(!openFile()) {
		return false;
	}
		
    _ifSAMFile.seekg (0, ios::beg);

	int cLines = 0;
	string line;
	bool foundAt = false;

	// scan 100 lines to find a line starting with @ indicating a header line
	while(cLines < 100 && getline(_ifSAMFile,line) && !foundAt) {
		if (!emptyLine(line) && line[0] == '@') {
            foundAt = true;
		}
		cLines++;
	}

    _ifSAMFile.clear();	
    _ifSAMFile.seekg (0, ios::beg);
    return foundAt;
}

// return the next contig in the SAM file, or NULL if there are no more contigs
Contig* SAMFile::nextContig() {
    Logger::getLogger()->log(QSNP_INFO, "SAMFile::nextContig()");
    openFile();
	string line;

	SAMContig contig;

	if(_pRead != NULL) {
		contig.setName(_pRead->getContigName());
		contig.addRead(_pRead);
	}

	int counter = 0;
	while(getline(_ifSAMFile, line)) {
		counter++;
		if (!line.empty() && line[0] != '@') {
			SAMRead* pRead = parseReadLine(line);
			if (pRead != NULL) {
				if(contig.getName().empty()) {
					contig.setName(pRead->getContigName());
				}

				if(pRead->getContigName() == contig.getName()) {
					contig.addRead(pRead);
				} else {
					_pRead = pRead;
					contig.stretchReads();
					contig.mergeReadPairs();
					contig.constructReferenceSequence();
					return contig.toContig();
				}
			}
		}
	}

	_pRead = NULL;
	contig.stretchReads();
	contig.mergeReadPairs();
	contig.constructReferenceSequence();
	return contig.toContig();
}

SAMRead* SAMFile::parseReadLine(const string& line) {
//SPWG_7257172	35	Contig_23	1	60	19S80M	=	3	120	CATTTATCGCCTCAGAGTTCATGGCTATGAAACCAAAACAGAGAGTGTAGAAGATGAAAGCATGAGATGATGATATTTGAATTTGCTGCTTCATATTGT	GHIIGIIIIDIIIIIIGIIIGIIIIIIIIHIIIIIIIIIFIIIHIDGFIGFBIGHIDIIIIHIIBIFHHHGIHGHGIHBHFIHHGHEGDBEEGE>EC>C	NH:i:1

	string readName, contigName, cigar, mrnm, seq, qual, tag, vtype, value, optional;
	int flag, nStart, nMapq, nMpos, nIsize;
	std::stringstream ssLine(line);
	ssLine >> readName >> flag >> contigName >> nStart >> nMapq >> cigar >> mrnm >> nMpos >> nIsize >> seq >> qual >> optional;

	SAMRead* pRead = new SAMRead(readName);
    Logger::getLogger()->log(QSNP_INFO, "parseReadLine for read: " + readName + " (" + contigName + ")");

	pRead->setSequence(seq);
	pRead->setContigName(contigName);
	pRead->setStartPosition(nStart);
	pRead->setQuality(qual);

	list<operation>& ops = pRead->getOperations();
	if(parseCigar(cigar, ops) && pRead->processOperations()) {
		return pRead;
	}

	// something when wrong reading the readline, so delete read and return NULL
	delete pRead;
	return NULL;
}

bool SAMFile::parseCigar(const string& cigar, list<operation>& ops) {
    if(cigar == "*") {
		Logger::getLogger()->log(QSNP_INFO, "skipping read with cigar *");
		return false;
    }

	size_t prevFound = 0;
	size_t found = cigar.find_first_of ("MIDNSHP", prevFound);

	while(found != string::npos) {
		char cOperation = cigar[found];
		int iSize;

		if(from_string<int>(iSize, cigar.substr(prevFound, found - prevFound), std::dec)) {
			operation op = {cOperation, iSize};
			ops.push_back(op);
		} else {
			Logger::getLogger()->log(QSNP_ERROR, "Parsing of cigar line failed: " + cigar);
			return false;
		}
		prevFound = found + 1;
		found = cigar.find_first_of ("MIDNSHP", prevFound);
	} 

	return true;
}


bool SAMFile::emptyLine(const string& line) {
	return line.empty() || (line.size() == 1 && line[0] == '\r');
}
