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

Contig* ContigProvider::nextContig() {
	return _pContigFile->nextContig();
}
