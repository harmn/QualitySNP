#include "Contig.h"
#include "ACEFile.h"
#include "Configuration.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string.h>
using namespace std;

void usage(char* programName) {
	cout << programName ;
	cout << " [-acefile <acefile>]";
	cout << " [-fastafile <fastafile>]";
	cout << " [-logLevel n]";
	cout << endl;
	cout << programName << " -help" << endl;

	cout << endl;

	cout << "   -acefile            a valid ace file as produced by CAP3 and other tools" << endl;
	cout << "   -fastafile          the name of the FASTA file" << endl;
	cout << "   -logLevel           logging level, 1 (only errors), 2 (warnings) or 3 (info) (default 1)" << endl;

	return;
}

bool parseOpt(int argc, char* argv[]) {
	map<string, string> optionMap;
	optionMap["acefile"]	= "ACEFileName";
	optionMap["fastafile"]	= "FASTAFileName";
	optionMap["logLevel"]	= "logLevel";

	int i = 1;
	Configuration* pConfig = Configuration::getConfig();
	while (i < argc) {
		if(argv[i][0] == '-') {
			char* currentOption = argv[i]+1;
			if(strcmp(currentOption,"help") == 0) {
				usage(argv[0]);
				return true;
			} else if(optionMap.count(currentOption) == 0) {
				cerr << "Unknown option: " << currentOption << endl;
				usage(argv[0]);
				return false;
			}

			string value;
			if (i + 1 < argc && argv[i + 1][0] != '-') {
				i++;
				value = argv[i];
			}

			bool bSuccess = pConfig->setValue(optionMap[currentOption], value);
			if (!bSuccess) {
				cerr << "Could not set value \"" << value << "\"";
				cerr << " for option \"" << currentOption << "\"" << endl;
				cerr << endl;
				usage(argv[0]);
				return false;
			}
		} else {
			pConfig->setString("ACEFileName", argv[i]);
		}
		i++;
	}

	return true;
}

int main(int argc, char* argv[]) {

if (argc > 1) {
		bool bParseSuccess = parseOpt(argc, argv);
		if (!bParseSuccess) {
			return 1;
		}
		Configuration* pConfig = Configuration::getConfig();
		// Open ACE file
		string ACEFilename	= pConfig->getString("ACEFileName");
		ACEFile af(ACEFilename);
		if(!af.openFile()) {
			cerr <<  "Could not open ACE file: " + ACEFilename;
			cerr << endl;
			return 1;
		}
		// Open FASTA file
		string FASTAFileName	= pConfig->getString("FASTAFileName");
		ofstream FASTAFile;
		FASTAFile.open((FASTAFileName.c_str()), ios::out);
		if(!FASTAFile.is_open()){
			cerr <<  "Could not open FASTA file: " + FASTAFileName;
			cerr << endl;
			return 1;
		}
		// Loop through contigs and print names and sequences
		Contig* pContig = af.nextContig();
		while (pContig != NULL) {
			// Create FASTA header from the name
			string Name = ">";
			Name.append(pContig->getName());
			string Sequence = pContig->getSequence();
			FASTAFile << Name << endl;
			FASTAFile << Sequence << endl;
			delete pContig;
			pContig = af.nextContig();
		}
		delete pContig;
		// Cleanup
		FASTAFile.close();
	} else {
		usage(argv[0]);
	}

	return EXIT_SUCCESS;
}


