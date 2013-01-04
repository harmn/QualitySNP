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

#include "CSVWriter.h"

CSVWriter::CSVWriter(){
    _bReadGroups = false;
    _bFail = false;
}

CSVWriter::CSVWriter(vector<string>& readGroupNames){
    // constructor with readgroups
    _readGroupNames = readGroupNames;
    _bReadGroups = true;
    _bFail = false;
}

CSVWriter::~CSVWriter(void)
{
	vector<ofstream*>::iterator itOutputs;
	for(itOutputs = _rgOutputs.begin(); itOutputs != _rgOutputs.end(); itOutputs++) {
		delete *itOutputs;
	}

	vector<ofstream*>::iterator itIndexFiles;
	for(itIndexFiles = _rgIndexFiles.begin(); itIndexFiles != _rgIndexFiles.end(); itIndexFiles++) {
		delete *itIndexFiles;
	}
}

bool CSVWriter::init() {
    char sep = Configuration::getConfig()->getChar("fieldSeparator");

    _rgOutputTypes.push_back("contigsFile");
    stringstream labels;
    labels << "name" << sep;
    labels << "Potential SNP" << sep;
    labels << "HighQuality SNP" << sep;
    labels << "Reliable SNP" << sep;
    labels << "Dvalue" << sep;
    labels << "reads" << sep;
    labels << "haplotypes" << sep;
    labels << "max haplotypes SNP" << sep << "sequence";

    _rgColumnNames.push_back(labels.str());

    _rgOutputTypes.push_back("readsFile");
    labels.str("");
    labels << "name" << sep << "hapid" << sep;
    labels << "contig" << sep << "start" << sep;
    labels <<"group" <<  sep << "sequence";
    _rgColumnNames.push_back(labels.str());

    _rgOutputTypes.push_back("haploTypesFile");
    labels.str("");
    labels << "contig" << sep << "hapid";
    _rgColumnNames.push_back(labels.str());

    _rgOutputTypes.push_back("variationsFile");
    labels.str("");
    labels << "contig" << sep << "position" << sep;
    labels << "major allele" << sep << "minor allele" << sep;
    labels << "high quality" << sep << "reliable" << sep;
    labels << "defining" << sep << "high quality flank";
    _rgColumnNames.push_back(labels.str());

    if(_bReadGroups) {
        _rgOutputTypes.push_back("readGroupsFile");

        vector<string>::iterator itRGNames;
        labels.str("");
        labels << "contig" << sep << "position";

        for(itRGNames = _readGroupNames.begin(); itRGNames != _readGroupNames.end(); itRGNames++) {
            labels << sep << (*itRGNames);
        }

        _rgColumnNames.push_back(labels.str());
    }


	vector<string>::iterator itOutput;
	int iFile = 0;
    for ( itOutput=_rgOutputTypes.begin(); itOutput < _rgOutputTypes.end(); itOutput++) {
		if (!openOutputFile(*itOutput)) {
            _bFail = true;
			return false;
		}
		
		(*_rgOutputs[iFile]) << _rgColumnNames[iFile] << NEWLINE;
		iFile++;
	}

	return true;
}

bool CSVWriter::writeContig(Contig* pContig) {
	if(_bFail) {
		return false;
	}

    string contigName = pContig->getName();

    (*_rgIndexFiles[CONTIGSFILE]) << contigName << '\t' << _rgOutputs[CONTIGSFILE]->tellp() << endl;

    (*_rgOutputs[CONTIGSFILE]) << pContig->toCSV();
	_rgOutputs[CONTIGSFILE]->flush();

    (*_rgIndexFiles[READSFILE]) << contigName << '\t' << _rgOutputs[READSFILE]->tellp() << endl;

	(*_rgOutputs[READSFILE]) << pContig->reads2CSV();
	_rgOutputs[READSFILE]->flush();

    (*_rgIndexFiles[HAPLOTYPESFILE]) << contigName << '\t' << _rgOutputs[HAPLOTYPESFILE]->tellp() << endl;

    (*_rgOutputs[HAPLOTYPESFILE]) << pContig->haploTypes2CSV();
    _rgOutputs[HAPLOTYPESFILE]->flush();

    (*_rgIndexFiles[VARIATIONSFILE]) << contigName << '\t' << _rgOutputs[VARIATIONSFILE]->tellp() << endl;

	(*_rgOutputs[VARIATIONSFILE]) << pContig->variations2CSV();
	_rgOutputs[VARIATIONSFILE]->flush();

    if(_bReadGroups) {
        (*_rgIndexFiles[READGROUPSFILE]) << contigName << '\t' << _rgOutputs[READGROUPSFILE]->tellp() << endl;

        (*_rgOutputs[READGROUPSFILE]) << pContig->readGroups2CSV(_readGroupNames);
        _rgOutputs[READGROUPSFILE]->flush();
    }

	return true;
}

bool CSVWriter::openOutputFile(const string& output) {
	Configuration* pConfig = Configuration::getConfig();
	string filename = pConfig->getString("outputDirectory") + "/" + pConfig->getString(output);

	ofstream* pOs = new ofstream();
	pOs->open(filename.c_str(),ios::out);
	if (pOs->fail()) {
        Logger::getLogger()->log(QSNP_ERROR, "Could not create output file: " + filename);
        delete pOs;
		return false;
	}

	_rgOutputs.push_back(pOs);

    ofstream* pIndex = new ofstream();
    filename += ".inx";
    pIndex->open(filename.c_str(),ios::out);
    if (pIndex->fail()) {
        Logger::getLogger()->log(QSNP_ERROR, "Could not create index file: " + filename);
        return false;
    }
    _rgIndexFiles.push_back(pIndex);

	return true;
}
