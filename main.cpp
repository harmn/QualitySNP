#include <QtGui/QApplication>
#include <QThread>
#include <QDir>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "mainwindow.h"
#include "core/trunk/HaploType.h"
#include "core/trunk/Variation.h"
#include "core/trunk/SeqRead.h"
#include "core/trunk/Contig.h"
#include "core/trunk/ContigProvider.h"
#include "core/trunk/QualitySNPpp.h"
#include "core/trunk/Configuration.h"
#include "core/trunk/ContigPrinter.h"
#include "core/trunk/CSVWriter.h"
#include "runqsnp.h"


using namespace std;

void usage(char* programName) {
    cout << programName << " [-outdir <outputdir>] [-minAlleles n]";
    cout << " [-lq5 n] [-lq3 n] [-lq3p n] [-lqWeight n] [-minQualScore n] [-minMapQuality n ]";
    cout << " [-minConf n] [-simPol n] [-simAll n] [-minHQReads n] [-reliableMarkers F/T]";
    cout << " [-maxSNPsFlank n] [-indelLQ T/F] [-printSummaryLine T/F]";
    cout << " [-printAlignment F/T] [-printHaploTypes F/T] [-useIUPACCodes F/T]";
    cout << " [-printSNPs F/T] [-printMarkers F/T] [-logLevel n] [-servermode F/T]";
    cout << " <acefile> " << endl;
    cout << programName << " -config <configurationfile> <contigfile>" << endl;
    cout << programName << " -help" << endl;

    cout << endl;

    cout << "   contigfile              a valid contig file in ACE or SAM format as produced by CAP3 and other tools" << endl;
    cout << "   -outdir                 the directory for the log file and the csv data files (default: /tmp)" << endl;
    cout << "   -lq5                    the number of nucleotides at the 5' end of each read that should be marked as low quality (default: 0)" << endl;
    cout << "   -lq3                    the number of nucleotides at the 3' end of each read that should be marked as low quality (default: 0)" << endl;
    cout << "   -lq3p                   the fraction (0-1) of nucleotides from the 3' end of each read that should be marked as low quality (default: 0)" << endl;
    cout << "   -minReads               minimal number of reads required for a valid allele (default: 2)" << endl;
    cout << "   -minReadsp              minimal number of reads required for a valid allele, as a fraction (0-1) (default: 0.0)" << endl;
    cout << "   -lqWeight               weight of the low quality nucleotides (default 0.5)" << endl;
    cout << "   -minQualScore           minimal contig quality score for high confidence (default: 20)" << endl;
    cout << "   -minMapQuality          minimal mapping quality score (default: 0)" << endl;
    cout << "   -minConf                minimal score for high confidence (default: 2)" << endl;
    cout << "   -simPol                 minimal similarity score per polymorphic site (default: 0.75)" << endl;
    cout << "   -simAll                 minimal similarity score over all polymorphic sites (default: 0.8)" << endl;
    cout << "   -minHQReads             minimal number of high quality reads for a SNP to be of high confidence" << endl;
    cout << "   -reliableMarkers        only use reliable SNPs for markers, not the high confidence ones (T/F, default: F)" << endl;
    cout << "   -maxSNPsFlank           maximal number of SNPs allowed in flanking regions for marker SNPs" << endl;
    cout << "   -useIUPACCodes          return the reference sequence with IUPAC codes for the variations" << endl;
    cout << "   -indelLQ                mark indels SNPs as low quality (T/F, default: T)" << endl;
    cout << "   -printAlignment         print an alignment of the contig and reads (T/F, default: F)" << endl;
    cout << "   -printHaploTypes        print the different haplotypes (T/F, default: F)" << endl;
    cout << "   -printSNPs              print the SNPs (T/F, default: F)" << endl;
    cout << "   -printMarkers           print the marker region (T/F, default: F)" << endl;
    cout << "   -printSummaryLine       print a summary line per contig (T/F, default: T)" << endl;
    cout << "   -logLevel               logging level, 1 (only errors), 2 (warnings) or 3 (info) (default 1)" << endl;
    cout << "   -showContigsWithoutSNP  report on all contigs, also the ones without SNPs (T/F, default; F)" << endl;
    cout << "   -config                 load configuration file" << endl;
    cout << "   -servermode             run in servermode (without graphical interface) (T/F, default F)" << endl;

    return;
}

bool parseOpt(int argc, char* argv[]) {
    map<string, string> optionMap;
    optionMap["config"]					= "configurationFile";
    optionMap["minReads"]				= "minimalNumberOfReadsPerAllele";
    optionMap["minReadsp"]				= "minimalNumberOfReadsPerAllelep";
    optionMap["outdir"]					= "outputDirectory";
    optionMap["lq5"]					= "lowQualityRegion5prime";
    optionMap["lq3"]					= "lowQualityRegion3prime";
    optionMap["lq3p"]					= "lowQualityRegion3primePerc";
    optionMap["lqWeight"]				= "weightLowQualityRegion";
    optionMap["minQualScore"]			= "minSNPQualityScore";
    optionMap["minConf"]				= "minimalConfidenceScore";
    optionMap["simPol"]					= "similarityPerPolymorphicSite";
    optionMap["simAll"]					= "similarityAllPolymorphicSites";
    optionMap["minHQReads"]				= "minNumberOfHighQualityReads";
    optionMap["minMapQuality"]			= "minimalMappingQuality";
    optionMap["maxSNPsFlank"]			= "maxNumberOfSNPsInFlanks";
    optionMap["indelLQ"]				= "indelLowQuality";
    optionMap["printAlignment"]			= "printAlignment";
    optionMap["printHaploTypes"]		= "printHaploTypes";
    optionMap["printSNPs"]				= "printSNPs";
    optionMap["printMarkers"]			= "printMarkers";
    optionMap["printSummaryLine"]		= "printSummaryLine";
    optionMap["logLevel"]				= "logLevel";
    optionMap["reliableMarkers"]		= "onlyReliableMarkers";
    optionMap["useIUPACCodes"]			= "useIUPACCodes";
    optionMap["showContigsWithoutSNP"]	= "showContigsWithoutSNP";
    optionMap["servermode"]             = "servermode";

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
                return false;
            }
        } else {
            pConfig->setString("contigFileName", argv[i]);
        }
        i++;
    }

    return true;
}

class QThreadEx : public QThread
{
protected:
    void run() { exec(); }
};

int main(int argc, char* argv[]) {
    Configuration* pConfig = Configuration::getConfig();

    bool bParseSuccess = parseOpt(argc, argv);
    if (!bParseSuccess && pConfig->getBool("servermode")) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (pConfig->getString("configurationFile") != "")  {
        pConfig->readFile(pConfig->getString("configurationFile"), QDir::currentPath().toStdString());
        parseOpt(argc, argv);
    }

    QApplication a(argc, argv, !pConfig->getBool("servermode"));

    if(pConfig->getBool("servermode")) {
        if(pConfig->getString("contigFileName").empty()) {
            cerr << "No input file given" <<endl;
            usage(argv[0]);
            return EXIT_FAILURE;
        }

        QThreadEx runThread;
        RunQSNP runqsnp;
        runqsnp.moveToThread(&runThread);

        runqsnp.connect(&runThread, SIGNAL(started()), &runqsnp, SLOT(run()), Qt::DirectConnection);
        runThread.connect(&runqsnp, SIGNAL(done()), &runThread, SLOT(quit()),Qt::DirectConnection);

        runThread.start();

        runThread.wait();

        Logger* pLogger = Logger::getLogger();
        pLogger->log(QSNP_ALWAYS, "QualitySNP done");

        return EXIT_SUCCESS;
    } else {
        MainWindow w;
        w.show();

        return a.exec();
    }
}
