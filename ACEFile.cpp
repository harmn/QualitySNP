#include <iostream>
#include <sstream>
#include "SeqRead.h"
#include "Contig.h"
#include "ACEFile.h"

using namespace std;

ACEFile::ACEFile(const string& contigFilename): _contigFilename(contigFilename), _pContig(NULL)
{
}

ACEFile::~ACEFile(void)
{
}

bool ACEFile::openFile() {
	if(_ifACEFile.is_open()) {
		return true;
	}

	_ifACEFile.open(_contigFilename.c_str(),ios::in);
	if (_ifACEFile.fail()) {
		Logger::getLogger()->log(QSNP_ERROR, "could not open contig file: " + _contigFilename);
		return false;
	}
	return true;
}

bool ACEFile::isValid() {
	if(!openFile()) {
		return false;
	}
		
    _ifACEFile.seekg (0, ios::beg);

	int cLines = 1;
	string line;
	bool foundCO = false;

	// scan 100 lines to find a CO line indicating a new contig block start
	while(cLines < 100 && getline(_ifACEFile,line) && !foundCO) {
		if (!emptyLine(line) && line[0] == 'C' && line[1] == 'O') {
            foundCO = true;
		}
		cLines++;
	};

    _ifACEFile.clear();	
    _ifACEFile.seekg (0, ios::beg);
    return foundCO;
}

Contig* ACEFile::nextContig() {
	//char line[MAX_LINE_LENGTH];
	string line;

	string contig_sequence;
	string read_sequence;

	Contig* pContig = _pContig;
	SeqRead* pRead = NULL;
	_pContig = NULL;

	enum ReadState {
		readingContig,
		readingRead,
		readingQuality,
		lookingForBlock
	};

	ReadState readState = readingContig;

	// find the CO line indicating a new contig block start
	while(pContig == NULL) {
		if(!getline(_ifACEFile, line)) {
			return NULL;
		}

		if (line[0] == 'C' && line[1] == 'O' && line[2] <= ' ') {
			pContig = parseContigLine(line);
		}
	}

	string data;

	// start parsing the contig text block
	while(getline(_ifACEFile,line)) {
		if(emptyLine(line)) {
			switch(readState) {
				case readingContig:
					pContig->setSequence(data);
					break;
				case readingRead:
					pRead->setSequence(data);
					break;
				case readingQuality:
					pContig->setQuality(data);
					break;
                default:
                    break;
			}
			data.clear();
			readState = lookingForBlock;
		} else if (readState == lookingForBlock) {
			// looking for new sub block
			switch (getBlockLabel(line)){
				case CO:
					// new contig found, store the new contig line
					// and return the current contig object
					_pContig = parseContigLine(line);
					return pContig;
					break;
				case RD:
					pRead = parseReadLine(line, pContig);
					if (pRead == NULL) {
						Logger::getLogger()->log(QSNP_ERROR, "Invalid read line" + line);
						return NULL;
					}
					readState = readingRead;
					break;
				case AF:
					pRead = parseReadLocationLine(line, pContig);
					pContig->addRead(pRead);
					break;
				case QA:
					parseQualitySegmentLine(line,pRead);
					break;
				case BQ:
					readState = readingQuality;
					break;
                default:
                    break;
			}
		} else {
			data += line;
		}
	}

	return pContig;
}

bool ACEFile::emptyLine(const string& line) {
	return line.empty() || (line.size() == 1 && line[0] == '\r');
}

BLOCKLABEL ACEFile::getBlockLabel(const string& line) {
	if(line[0] == 'C' && line[1] == 'O' && line[2] <= ' ') {
		return CO;
	}
	if(line[0] == 'R' && line[1] == 'D' && line[2] <= ' ') {
		return RD;
	}
	if(line[0] == 'A' && line[1] == 'F' && line[2] <= ' ') {
		return AF;
	}
	if(line[0] == 'Q' && line[1] == 'A' && line[2] <= ' ') {
		return QA;
	}
	if(line[0] == 'B' && line[1] == 'Q' && line[2] <= ' ') {
		return BQ;
	}

	return UnknownTag;
}

Contig* ACEFile::parseContigLine(const string& line) {
	string name, tag;
	std::stringstream ssLine(line);
	ssLine >> tag >> name;
	return new Contig(name);
}

SeqRead* ACEFile::parseReadLine(const string& line, Contig* pContig) {
	if (pContig == NULL) {
		Logger::getLogger()->log(QSNP_ERROR, "Found a read line, but no active contig " + string(line));

		return NULL;
	}

	string name, tag;
	std::stringstream ssLine(line);
	ssLine >> tag >> name;
	stripTrailingFields(name);
	SeqRead* pRead = pContig->findRead(name);

	return pRead;
}

SeqRead* ACEFile::parseReadLocationLine(const string& line, Contig* pContig) {
	string tag, readName, CorU;
	int position;
	std::stringstream ssLine(line);
	ssLine >> tag >> readName >> CorU >> position;
	stripTrailingFields(readName);
	SeqRead* pRead = new SeqRead(readName, pContig);
	pRead->setStartPosition(position);

	return pRead;
}


void ACEFile::parseQualitySegmentLine(const string& line, SeqRead* pRead) {
	string tag, readName, CorU;
	int startPosition;
	int endPosition;
	std::stringstream ssLine(line);
	ssLine >> tag >> startPosition >> endPosition;

	pRead->setQualClip(startPosition, endPosition);
}

void ACEFile::stripTrailingFields(string& strIn) {
	size_t stDelimiterPos = strIn.find_first_of("|");
	if (stDelimiterPos > 0) {
		strIn = strIn.substr(0,stDelimiterPos);
	}
}
