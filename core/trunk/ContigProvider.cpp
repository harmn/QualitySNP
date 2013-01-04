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

#include "HaploType.h"
#include "Variation.h"
#include "SeqRead.h"
#include "Contig.h"
#include "SAMFile.h"
#include "ContigProvider.h"

ContigProvider::ContigProvider(void): _pContigFile(NULL)
{
}

ContigProvider::~ContigProvider(void)
{
	delete _pContigFile;
}

bool ContigProvider::init() {
	ios_base::sync_with_stdio(false);
	string contigFilename = Configuration::getConfig()->getString("contigFileName");

	ifstream contigFile;
	contigFile.open(contigFilename.c_str(),ios::in);
	if (contigFile.fail()) {
		Logger::getLogger()->log(QSNP_ERROR, "could not open contig file: " + contigFilename);
		return false;
	}
	contigFile.close();


	ACEFile* pACEFile = new ACEFile(contigFilename);
	
	if(pACEFile->isValid()) {
		Logger::getLogger()->log(QSNP_INFO, "Starting to read ACE file: " + contigFilename);
		_pContigFile = pACEFile;
		return true;
	} 

    if(contigFilename.compare(contigFilename.size() -3, 3, "ace") == 0) {
		Logger::getLogger()->log(QSNP_INFO, "Starting to read ACE file: " + contigFilename);
		_pContigFile = pACEFile;
		return true;
    }
	
	delete pACEFile;

	SAMFile* pSAMFile = new SAMFile(contigFilename);

	if(pSAMFile->isValid()) {
		Logger::getLogger()->log(QSNP_INFO, "Starting to read SAM file: " + contigFilename);
		_pContigFile = pSAMFile;
		return true;
	}

    
    if(contigFilename.compare(contigFilename.size() -3, 3, "sam") == 0) {
		Logger::getLogger()->log(QSNP_INFO, "Starting to read SAM file: " + contigFilename);
		_pContigFile = pSAMFile;
        cerr << "SAM" << endl;
		return true;
    }

    delete pSAMFile;

	Logger::getLogger()->log(QSNP_ERROR, "No valid file format: " + contigFilename);

    return false;
}

int ContigProvider::getReadCount()
{
    return _pContigFile->readCount();
}

int ContigProvider::getContigCount()
{
    return _pContigFile->contigCount();
}

vector<string> &ContigProvider::getReadGroups()
{
    return _pContigFile->readGroups();
}

Contig* ContigProvider::nextContig() {
	return _pContigFile->nextContig();
}
