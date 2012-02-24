#include <iostream>
#include <iomanip>
#include <fstream>
#include "HaploType.h"
#include "Variation.h"
#include "SeqRead.h"
#include "Contig.h"
#include "ContigProvider.h"
#include "QualitySNPpp.h"
#include "Configuration.h"
#include "ContigPrinter.h"
#include "CSVWriter.h"
#include <stdlib.h>
#include <string.h>

using namespace std;


void usage(char* programName) {
	cout << programName << " [-outdir <outputdir>] [-minAlleles n]";
	cout << " [-lq5 n] [-lq3 n] [-lq3p n] [-lqWeight n] [-minQualScore n]";
	cout << " [-minConf n] [-simPol n] [-simAll n] [-minHQReads n] [-reliableMarkers F/T]";
	cout << " [-maxSNPsFlank n] [-indelLQ T/F] [-printSummaryLine T/F]";
	cout << " [-printAlignment F/T] [-printHaploTypes F/T] [-useIUPACCodes F/T]";
	cout << " [-printSNPs F/T] [-printMarkers F/T] [-logLevel n]";
	cout << " <acefile> " << endl;
	cout << programName << " -config <configurationfile> <contigfile>" << endl;
	cout << programName << " -help" << endl;

	cout << endl;

	cout << "   contigfile          a valid contig file in ACE or SAM format as produced by CAP3 and other tools" << endl;
	cout << "   -outdir             the directory for the log file and the csv data files (default: /tmp)" << endl;
	cout << "   -lq5                the number of nucleotides at the 5' end of each read that should be marked as low quality (default: 0)" << endl;
	cout << "   -lq3                the number of nucleotides at the 3' end of each read that should be marked as low quality (default: 0)" << endl;
	cout << "   -lq3p               the fraction (0-1) of nucleotides from the 3' end of each read that should be marked as low quality (default: 0)" << endl;
	cout << "   -minReads           minimal number of reads required for a valid allele (default: 2)" << endl;
	cout << "   -minReadsp          minimal number of reads required for a valid allele, as a fraction (0-1) (default: 0.0)" << endl;
	cout << "   -lqWeight           weight of the low quality nucleotides (default 0.5)" << endl;
	cout << "   -minQualScore       minimal contig quality score for high confidence (default: 20)" << endl;
	cout << "   -minConf            minimal score for high confidence (default: 2)" << endl;
	cout << "   -simPol             minimal similarity score per polymorphic site (default: 0.75)" << endl;
	cout << "   -simAll             minimal similarity score over all polymorphic sitea (default: 0.8)" << endl;
	cout << "   -minHQReads         minimal number of high quality reads for a SNP to be of high confidence" << endl;
	cout << "   -reliableMarkers	only use reliable SNPs for markers, not the high confidence ones (default: F)" << endl;
	cout << "   -maxSNPsFlank       maximal number of SNPs allowed in flanking regions for marker SNPs" << endl;
	cout << "   -useIUPACCodes      return the reference sequence with IUPAC codes for the variations" << endl;
	cout << "   -indelLQ            mark indels SNPs as low quality (default: T)" << endl;
	cout << "   -printAlignment     print an alignment of the contig and reads (default: F)" << endl;
	cout << "   -printHaploTypes    print the different haplotypes (default: F)" << endl;
	cout << "   -printSNPs          print the SNPs (default: F)" << endl;
	cout << "   -printMarkers       print the marker region (default: F)" << endl;
	cout << "   -printSummaryLine   print a summary line per contig (default: T)" << endl;
	cout << "   -logLevel           logging level, 1 (only errors), 2 (warnings) or 3 (info) (default 1)" << endl;
	cout << "   -config             load configuration file" << endl;
	
	return;
}

bool parseOpt(int argc, char* argv[]) {
	map<string, string> optionMap;
	optionMap["config"]				= "configurationFile";
	optionMap["minReads"]			= "minimalNumberOfReadsPerAllele";
	optionMap["minReadsp"]			= "minimalNumberOfReadsPerAllelep";
	optionMap["outdir"]				= "outputDirectory";
	optionMap["lq5"]				= "lowQualityRegion5prime";
	optionMap["lq3"]				= "lowQualityRegion3prime";
	optionMap["lq3p"]				= "lowQualityRegion3primePerc";
	optionMap["lqWeight"]			= "weightLowQualityRegion";
	optionMap["minQualScore"]		= "minSNPQualityScore";
	optionMap["minConf"]			= "minimalConfidenceScore";
	optionMap["simPol"]				= "similarityPerPolymorphicSite";
	optionMap["simAll"]				= "similarityAllPolymorphicSites";
	optionMap["minHQReads"]			= "minNumberOfHighQualityReads";
	optionMap["maxSNPsFlank"]		= "maxNumberOfSNPsInFlanks";
	optionMap["indelLQ"]			= "indelLowQuality";
	optionMap["printAlignment"]		= "printAlignment";
	optionMap["printHaploTypes"]	= "printHaploTypes";
	optionMap["printSNPs"]			= "printSNPs";
	optionMap["printMarkers"]		= "printMarkers";
	optionMap["printSummaryLine"]	= "printSummaryLine";
	optionMap["logLevel"]			= "logLevel";
	optionMap["reliableMarkers"]	= "onlyReliableMarkers";
	optionMap["useIUPACCodes"]		= "useIUPACCodes";
	
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
			pConfig->setString("contigFileName", argv[i]);
		}
		i++;
	}

	return true;
}

int main(int argc, char* argv[]) {
	if (argc <= 1) {
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	bool bParseSuccess = parseOpt(argc, argv);
	if (!bParseSuccess) {
		return EXIT_FAILURE;
	}

	Configuration* pConfig = Configuration::getConfig();
	if (pConfig->getString("configurationFile") != "")  {
		pConfig->readFile(pConfig->getString("configurationFile"));
	}
	Logger* pLogger = Logger::getLogger();
	pLogger->log(QSNP_ALWAYS, "QualitySNP started with settings: " + pConfig->toString());

	ContigProvider contigProvider;
	if (!contigProvider.init()) {
		return EXIT_FAILURE;	
	}

	CSVWriter csvWriter;
	if (!csvWriter.init()) {
		return EXIT_FAILURE;
	}

	Contig* pContig = contigProvider.nextContig();
	while (pContig != NULL) {
		pContig->sortReads();
		pContig->calculateProperties();
		if (pContig->getPotentialSNPCount() > 0) {
			ContigPrinter contigPrinter(pContig);
			cout << contigPrinter.toString();
			csvWriter.writeContig(pContig);
		}
		delete pContig;
		pContig = contigProvider.nextContig();
	}

	return EXIT_SUCCESS;
}


