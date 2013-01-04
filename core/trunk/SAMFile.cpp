/*
* This File is part of QualitySNP; a program to detect Single Nucleotide Variations
* https://trac.nbic.nl/qualitysnp/
*
*   Copyright (C) 2012 Harm Nijveen
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

SAMFile::SAMFile(const string& contigFilename): _contigFilename(contigFilename), _pRead(NULL), _cReads(-1)
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
	Configuration* pConfig = Configuration::getConfig();

	if(_pRead != NULL) {
		contig.setName(_pRead->getContigName());
		contig.addRead(_pRead);
	}

	int counter = 0;
	while(getline(_ifSAMFile, line)) {
		counter++;
		if (!line.empty() && line[0] != '@') {
            SAMRead* pRead = parseReadLine(line, false);
			if (pRead != NULL) {
				if (pRead->getMapQuality() < pConfig->getInt("minimalMappingQuality")) {
					Logger::getLogger()->log(QSNP_INFO, "Skipping read (" + pRead->getName() + ") with low mapping quality");
                    continue;
				}
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

SAMRead* SAMFile::parseReadLine(const string& line, bool bCollectStatics) {
//SPWG_7257172	35	Contig_23	1	60	19S80M	=	3	120	CATTTATCGCCTCAGAGTTCATGGCTATGAAACCAAAACAGAGAGTGTAGAAGATGAAAGCATGAGATGATGATATTTGAATTTGCTGCTTCATATTGT	GHIIGIIIIDIIIIIIGIIIGIIIIIIIIHIIIIIIIIIFIIIHIDGFIGFBIGHIDIIIIHIIBIFHHHGIHGHGIHBHFIHHGHEGDBEEGE>EC>C	NH:i:1

    string readName, contigName, cigar, mrnm, seq, qual, optional;
	int flag, nStart, nMapq, nMpos, nIsize;
	std::stringstream ssLine(line);
    ssLine >> readName;
    ssLine >> flag;
    ssLine >> contigName;
    ssLine >> nStart;
    ssLine >> nMapq;
    ssLine >> cigar;
    ssLine >> mrnm;
    ssLine >> nMpos;
    ssLine >> nIsize;
    ssLine >> seq;
    ssLine >> qual;

    if(contigName.empty() || contigName == "*" ) {
        return NULL;
    }

    string group;
    while(ssLine >> optional) {
        if(optional.substr(0,3).compare("RG:") == 0) {
            size_t valueStart = optional.find(":", 4);
            group = optional.substr(valueStart + 1);
        }
    }

	SAMRead* pRead = new SAMRead(readName);
    pRead->setGroup(group);
    pRead->setContigName(contigName);

    if(bCollectStatics) {
        return pRead;
    }

    Logger::getLogger()->log(QSNP_DEBUG, "parseReadLine for read: " + readName + " (" + contigName + ")");

    pRead->setSequence(seq);
    pRead->setStartPosition(nStart);
    pRead->setQuality(qual);
    pRead->setMapQuality(nMapq);

	list<operation>& ops = pRead->getOperations();
	if(parseCigar(cigar, ops) && pRead->processOperations()) {
		return pRead;
	}

    // something went wrong reading the readline, so delete read and return NULL
	delete pRead;
    return NULL;
}

bool SAMFile::collectStatistics()
{
    if(!openFile()) {
        return false;
    }

    _ifSAMFile.seekg (0, ios::beg);

    _cReads = 0;
    _cContigs = 0;
    _readGroups.clear();

    string line;
    string contigName;
    map<string,int> groupMap;

    while(getline(_ifSAMFile,line)) {
        if (!emptyLine(line) && line[0] != '@') {
            SAMRead* pRead = parseReadLine(line, true);
            if (pRead != NULL) {
                _cReads++;
                if(!pRead->getGroup().empty()) {
                    groupMap[pRead->getGroup()]++;
                }
                if(contigName != pRead->getContigName()) {
                    contigName = pRead->getContigName();
                    _cContigs++;
                }
                delete pRead;
            }
        }
    }


    map<string,int>::iterator itGroupMap;
    for(itGroupMap = groupMap.begin(); itGroupMap != groupMap.end(); itGroupMap++) {
        _readGroups.push_back((*itGroupMap).first);
    }

    _ifSAMFile.clear();
    _ifSAMFile.seekg (0, ios::beg);
    return true;
}

bool SAMFile::parseCigar(const string& cigar, list<operation>& ops) {
    if(cigar == "*") {
		Logger::getLogger()->log(QSNP_DEBUG, "skipping read with cigar *");
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

int SAMFile::readCount()
{
    if(_cReads == -1) {
        collectStatistics();
    }

    return _cReads;
}

int SAMFile::contigCount()
{
    if(_cContigs == -1) {
        collectStatistics();
    }

    return _cContigs;
}

vector<string> &SAMFile::readGroups()
{
    if(_cReads == -1) {
        collectStatistics();
    }

    return _readGroups;
}
