#include <iostream>
#include <fstream>
#include "SeqRead.h"
#include "Variation.h"
#include "HaploType.h"
#include "CSVReader.h"

CSVReader::CSVReader()
{
}

Contig* CSVReader::getContig(string contigName) {
    ifstream	ifContigFile;

    Configuration* pConfig = Configuration::getConfig();
    string contigPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("contigsFile");

    ifContigFile.open(contigPath.c_str(),ios::in);

    if(_contigsIndex.count(contigName) != 0) {
        ifContigFile.seekg(_contigsIndex[contigName], ios::beg);
    }

    Contig* pContig = NULL;
    string line;
    while(getline(ifContigFile, line) && pContig == NULL) {
        string name, sequence;
        int potSNP, hqSNP, relSNP, reads, haplotypes, maxHapSNP;
        double Dvalue;
        char tab;

        std::stringstream ssLine(line);
        // trick with noskipsws and tab to get the complete sequencing including any leading or internal whitespace in the sequence
        ssLine >> name >> potSNP >> hqSNP >> relSNP >> Dvalue >> reads >> haplotypes >> maxHapSNP >> noskipws >> tab;
        getline(ssLine,sequence);
        if(name == contigName) {
            pContig = new Contig(name);
            pContig->setSequence(sequence);

            map<int, vector<SeqRead*> > hapmap;
            getReads(pContig, hapmap);
            getVariations(pContig);
            getHaploTypes(pContig, hapmap);
        }
    }

    ifContigFile.close();

    return pContig;
}

const vector<ContigInfo> CSVReader::getContigList()
{
    Configuration* pConfig = Configuration::getConfig();

    readIndices();

    ifstream	ifContig;
    string contigPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("contigsFile");
    ifContig.open(contigPath.c_str(),ios::in);

    vector<ContigInfo> contigInfoList;
    string line;

    while(getline(ifContig, line)) {
        string name, sequence;
        int potSNP, hqSNP, relSNP, reads, haplotypes, maxHapSNP;
        double Dvalue;

        std::stringstream ssLine(line);
        ssLine >> name >> potSNP >> hqSNP >> relSNP >> Dvalue >> reads >> haplotypes >> maxHapSNP >> sequence;
        if(!sequence.empty()) {
            ContigInfo contigInfo;
            contigInfo.name = name;
            contigInfo.potSNP = potSNP;
            contigInfo.hqSNP = hqSNP;
            contigInfo.relSNP = relSNP;
            contigInfo.Dvalue = Dvalue;
            contigInfo.reads = reads;
            contigInfo.haplotypes = haplotypes;
            contigInfo.maxHapSNP = maxHapSNP;
            contigInfo.sequence = sequence;
            contigInfoList.push_back(contigInfo);
        }
    }

    ifContig.close();

    return contigInfoList;
}



void CSVReader::getReads(Contig* pContig, map<int, vector<SeqRead*> >& hapmap) {
    Configuration* pConfig = Configuration::getConfig();
    string readsPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("readsFile");

    ifstream	ifReads;

    ifReads.open(readsPath.c_str(),ios::in);

    if(_readsIndex.count(pContig->getName())  != 0) {
        ifReads.seekg(_readsIndex[pContig->getName()], ios::beg);
    }

    string line;
    bool bFound = false;

    while(getline(ifReads, line)) {
        string name, contigName, group, sequence;
        int hapid, start;

        std::stringstream ssLine(line);
        ssLine >> name >> hapid >> contigName >> start >> group >> ws;
        getline(ssLine, sequence);

        if(contigName == pContig->getName()) {
            bFound = true;
            SeqRead* pRead = new SeqRead(name, pContig);
            pRead->setSequence(sequence);
            pRead->setStartPosition(start + 1);

            if(group != "-") {
                pRead->setGroup(group);
            }

            hapmap[hapid].push_back(pRead);

            pContig->addRead(pRead);
        } else if(bFound) {
            ifReads.close();
            return;
        }
    }

    ifReads.close();
}


void CSVReader::getHaploTypes(Contig* pContig, map<int, vector<SeqRead*> >& hapmap) {
    Configuration* pConfig = Configuration::getConfig();
    string haploTypesPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("haploTypesFile");

    ifstream	ifHaploTypes;

    ifHaploTypes.open(haploTypesPath.c_str(),ios::in);

    if(_haploTypesIndex.count(pContig->getName()) != 0) {
        ifHaploTypes.seekg(_haploTypesIndex[pContig->getName()], ios::beg);
    }

    bool bFound = false;

    string line;
    while(getline(ifHaploTypes, line)) {
        string contigName;
        int hapid = -1;

        std::stringstream ssLine(line);
        ssLine >> contigName >> hapid;

        if(contigName == pContig->getName()) {
            bFound = true;
            HaploType* pHaploType = new HaploType(pContig, hapid);
            pContig->addHaploType(pHaploType);
            if(hapmap.count(hapid) != 0) {
                vector<SeqRead*>::iterator itRead;
                for(itRead = hapmap[hapid].begin(); itRead != hapmap[hapid].end(); itRead++) {
                    pHaploType->addRead(*itRead);
                }
            }
        } else if (bFound) {
            ifHaploTypes.close();
        }
    }

    ifHaploTypes.close();
}

void CSVReader::readIndices()
{
    Configuration* pConfig = Configuration::getConfig();

    string contigsPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("contigsFile") + ".inx";
    string variationsPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("variationsFile") + ".inx";
    string readsPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("readsFile") + ".inx";
    string haploTypesPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("haploTypesFile") + ".inx";

    ifstream ifContigIndex;
    ifstream ifVariationsIndex;
    ifstream ifReadsIndex;
    ifstream ifHaploTypesIndex;

    ifContigIndex.open(contigsPath.c_str(),ios::in);
    string line;
    while(getline(ifContigIndex, line)) {
        string contigName;
        long position;
        stringstream ssLine(line);
        ssLine >> contigName >> position;
        _contigsIndex[contigName] = position;
    }
    ifContigIndex.close();

    ifVariationsIndex.open(variationsPath.c_str(),ios::in);
    while(getline(ifVariationsIndex, line)) {
        string contigName;
        long position;
        stringstream ssLine(line);
        ssLine >> contigName >> position;
        _variationsIndex[contigName] = position;
    }
    ifVariationsIndex.close();

     ifReadsIndex.open(readsPath.c_str(),ios::in);
     while(getline(ifReadsIndex, line)) {
        string contigName;
        long position;
        stringstream ssLine(line);
        ssLine >> contigName >> position;
        _readsIndex[contigName] = position;
    }
    ifReadsIndex.close();

    ifHaploTypesIndex.open(haploTypesPath.c_str(),ios::in);
    while(getline(ifHaploTypesIndex, line)) {
        string contigName;
        long position;
        stringstream ssLine(line);
        ssLine >> contigName >> position;
        _haploTypesIndex[contigName] = position;
    }
    ifHaploTypesIndex.close();
}



void CSVReader::getVariations(Contig* pContig) {
    Configuration* pConfig = Configuration::getConfig();
    string variationsPath = pConfig->getString("inputDirectory") + "/" + pConfig->getString("variationsFile");

    ifstream	ifVariations;

    ifVariations.open(variationsPath.c_str(),ios::in);

    if(_variationsIndex.count(pContig->getName()) != 0) {
        ifVariations.seekg(_variationsIndex[pContig->getName()], ios::beg);
    }

    bool bFound = false;

    string line;
    while(getline(ifVariations, line)) {
        string contigName;
        int position;
        char majAllele;
        char minAllele;
        int highConfidence;
        int reliable;
        int defining;
        int hqFlank;
        std::stringstream ssLine(line);

        //major allele	minor allele	high quality	reliabldefining	high quality flank

        ssLine >> contigName >> position >> majAllele >> minAllele >> highConfidence >> reliable >> defining >> hqFlank;

        if(contigName == pContig->getName()) {
            bFound = true;
            Variation* pVariation = new Variation(pContig, position);
            pVariation->setMajorAllele(majAllele);
            pVariation->setMinorAllele(minAllele);
            pVariation->setIsDefining(defining);
            pVariation->setIsHighConfidence(highConfidence);
            pVariation->setIsReliable(reliable);
            pVariation->setFlankLength(hqFlank);
            pContig->addVariation(pVariation);
        } else if (bFound) {
            ifVariations.close();
            return;
        }
    }

    ifVariations.close();
}
