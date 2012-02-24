#include "CSVWriter.h"

CSVWriter::CSVWriter(void)
{
	_rgOutPutTypes.push_back("contigsFile");
	_rgColumnNames.push_back("name\tPotential SNP\tHighQuality SNP\tReliable SNP\tDvalue\treads\thaplotypes\tsingle read haplotypes\tmax haplotypes SNP\tsequence");

	_rgOutPutTypes.push_back("readsFile");
	_rgColumnNames.push_back("name\thapid\tcontig\tstart\tsequence");

	_rgOutPutTypes.push_back("haploTypesFile");
	_rgColumnNames.push_back("contig\thapid");
	
	_rgOutPutTypes.push_back("variationsFile");
	_rgColumnNames.push_back("contig\tposition\tmajor allele\tminor allele\thigh quality\treliable\tdefining\thigh quality flank");
}

CSVWriter::~CSVWriter(void)
{
	vector<ofstream*>::iterator itOutputs;
	for(itOutputs = _rgOutputs.begin(); itOutputs != _rgOutputs.end(); itOutputs++) {
		delete *itOutputs;
	}
}

bool CSVWriter::init() {
	_bFail = true;


	vector<string>::iterator itOutput;
	int iFile = 0;
	for ( itOutput=_rgOutPutTypes.begin(); itOutput < _rgOutPutTypes.end(); itOutput++) {
		if (!openOutputFile(*itOutput)) {
			return false;
		}
		
		(*_rgOutputs[iFile]) << _rgColumnNames[iFile] << NEWLINE;
		iFile++;
	}

	_bFail = false;

	return true;
}

bool CSVWriter::writeContig(Contig* pContig) {
	if(_bFail) {
		return false;
	}
	(*_rgOutputs[CONTIGSFILE]) << pContig->toCSV();
	_rgOutputs[CONTIGSFILE]->flush();

	(*_rgOutputs[READSFILE]) << pContig->reads2CSV();
	_rgOutputs[READSFILE]->flush();

	(*_rgOutputs[HAPLOTYPESFILE]) << pContig->haploTypes2CSV();
	_rgOutputs[HAPLOTYPESFILE]->flush();

	(*_rgOutputs[VARIATIONSFILE]) << pContig->variations2CSV();
	_rgOutputs[VARIATIONSFILE]->flush();


	return true;
}

bool CSVWriter::openOutputFile(const string& output) {
	Configuration* pConfig = Configuration::getConfig();
	string filename = pConfig->getString("outputDirectory") + "/" + pConfig->getString(output);

	ofstream* pOs = new ofstream();
	pOs->open(filename.c_str(),ios::out);
	if (pOs->fail()) {
		Logger::getLogger()->log(QSNP_ERROR, "Could not open output file: " + filename);
		return false;
	}

	_rgOutputs.push_back(pOs);

	return true;
}
