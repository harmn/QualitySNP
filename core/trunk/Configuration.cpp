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
#include <string>
#include <algorithm>
#include <QDir>
#include <dirent.h>
#include "Configuration.h"

Configuration* Configuration::s_pTheConfiguration = 0;

std::string trim(std::string& s,const std::string& drop = " \t\"'")
{
    std::string r=s.erase(s.find_last_not_of(drop)+1);
    return r.erase(0,r.find_first_not_of(drop));
}

Configuration::Configuration(void)
{
    _name = "default";

    intMap["minimalNumberOfReadsPerAllele"]		= 2;
    doubleMap["minimalNumberOfReadsPerAllelep"] = 0.1;
    intMap["minimalNumberOfReadsPerHaploType"]	= 2;
    intMap["lowQualityRegion5prime"]			= 0;
    intMap["lowQualityRegion3prime"]			= 0;
    intMap["maxNumberOfSNPsInFlanks"]			= 0;
    doubleMap["lowQualityRegion3primePerc"]		= 0.0;
    doubleMap["weightHighQualityRegion"]		= 1.0;
    doubleMap["weightLowQualityRegion"]			= 0.5;
    intMap["minSNPQualityScore"]				= 20;
    intMap["minimalConfidenceScore"]			= 5;
    intMap["minNumberOfHighQualityReads"]		= 2;
    intMap["logLevel"]							= QSNP_WARNING;
    intMap["minimalMappingQuality"]				= 0;
    intMap["maxNumberOfReads"]                  = 0;
    intMap["lowComplexityRegionSize"]           = 6;
    intMap["lowComplexityRepeatCount"]          = 5;
    doubleMap["similarityPerPolymorphicSite"]	= 0.75;
    doubleMap["similarityAllPolymorphicSites"]	= 0.8;
    doubleMap["alleleMajorityThreshold"]		= 0.75;
    stringMap["outputDirectory"]				= "/tmp";
    stringMap["inputDirectory"]                 = "/tmp";
    stringMap["contigFileName"]                 = "";
    stringMap["logFileName"]					= "QSNP.log";
    stringMap["contigsFile"]					= "contigs.csv";
    stringMap["readsFile"]						= "reads.csv";
    stringMap["haploTypesFile"]					= "haplotypes.csv";
    stringMap["variationsFile"]					= "variations.csv";
    stringMap["readGroupsFile"]                 = "readgroups.csv";
    stringMap["readNameGroupSeparator"]         = "";
    stringMap["configurationFile"]              = "";

    charMap["fieldSeparator"]					= '\t';
    boolMap["indelLowQuality"]					= true;
    boolMap["lowComplexityLowQuality"]          = true;
    boolMap["printAlignment"]					= false;
    boolMap["printHaploTypes"]					= false;
    boolMap["printSNPs"]						= false;
    boolMap["printMarkers"]						= false;
    boolMap["printSummaryLine"]					= true;
    boolMap["useIUPACCodes"]					= false;
    boolMap["onlyReliableMarkers"]				= true;
    boolMap["showContigsWithoutSNP"]			= true;
    boolMap["collectStatistics"]                = true;
    boolMap["outputReadGroups"]                 = true;
    boolMap["servermode"]                       = false;

    nucMap = new int[256];
    for (int i = 0; i < 256; i++) {
        nucMap[i] = -1;
    }

    nucMap['A'] = 0;
    nucMap['C'] = 1;
    nucMap['G'] = 2;
    nucMap['T'] = 3;
    nucMap['*'] = 4;
    rnucMap = "ACGT*";

    IUPACCode[03] = 'W';
    IUPACCode[12] = 'S';
    IUPACCode[02] = 'R';
    IUPACCode[13] = 'Y';
    IUPACCode[01] = 'M';
    IUPACCode[23] = 'K';
}

Configuration::~Configuration(void)
{
    delete[] nucMap;
}

map<string, string> Configuration::findConfigurations() {
    map<string, string> configurations;
    configurations["default"] = "";

    QDir confDir = QDir::home();
    if(!confDir.exists(CONF_DIR)) {
        Logger::getLogger()->log(QSNP_INFO, "No configuration directory, first time?" + confDir.path().toStdString());
        return configurations;
    }

    confDir.cd(CONF_DIR);
    QStringList filters;
    filters << "*.cfg";

    QStringList configurationFiles = confDir.entryList(filters);

    QStringListIterator itConfigurationFiles(configurationFiles);
    while(itConfigurationFiles.hasNext()) {
        Configuration conf;
        string filename = itConfigurationFiles.next().toStdString();
        conf.readFile(filename, "");
        string confName = conf.getName();
        configurations[confName] = filename;
    }

    return configurations;
}

bool Configuration::setInt(const string& key, int value) {
    if(intMap.count(key) == 0) {
        return false;
    }

    intMap[key] = value;
    return true;
}

bool Configuration::setChar(const string& key, char value) {
    if(charMap.count(key) == 0) {
        return false;
    }

    charMap[key] = value;
    return true;
}

bool Configuration::setDouble(const string& key, double value) {
    if(doubleMap.count(key) == 0) {
        return false;
    }

    doubleMap[key] = value;
    return true;
}

bool Configuration::setString(const string& key, const string& value) {
    if(stringMap.count(key) == 0) {
        return false;
    }

    stringMap[key] = value;
    return true;
}

bool Configuration::setBool(const string& key, bool value) {
    if(boolMap.count(key) == 0) {
        return false;
    }

    boolMap[key] = value;
    return true;
}


const string Configuration::getValue(const string &key)
{
    if(stringMap.count(key) > 0) {
        return getString(key);
    }

    stringstream converter;

    if(intMap.count(key) > 0) {
        converter << getInt(key);
    } else if(charMap.count(key) > 0) {
        converter << getChar(key);
    } else if (doubleMap.count(key) > 0) {
        converter << getDouble(key);
    } else if (boolMap.count(key) > 0) {
        converter << getBool(key);
    } else {
        return "";
    }

    return converter.str();
}

bool Configuration::setValue(const string& label, const string& value) {
    if(label == "name") {
        _name = value;
        return true;
    }

    istringstream converter(value);
    if(stringMap.count(label) > 0) {
        setString(label,value);
    } else if(intMap.count(label) > 0) {
        int nValue;
        if(!(converter >> nValue)) {
            return false;
        }
        setInt(label,nValue);
    } else if(charMap.count(label) > 0) {
        char cValue;
        if (value == "\\t") {
            cValue = '\t';
        } else {
            if((converter >> cValue)) {
                setChar(label,cValue);
            } else if(!value.empty()) {
                return false;
            }
        }
    } else if (doubleMap.count(label) > 0) {
        double dValue;
        if(!(converter >> dValue)) {
            return false;
        }
        setDouble(label,dValue);
    } else if (boolMap.count(label) > 0) {
        bool bValue;
        string tmpValue = value;
        std::transform(tmpValue.begin(), tmpValue.end(), tmpValue.begin(), ::tolower);
        bValue = (tmpValue == "true" || tmpValue == "t" || tmpValue == "1" || tmpValue.empty());
        setBool(label,bValue);
    } else {
        return false;
    }

    return true;
}

bool Configuration::writeFile(const string& strFileName, const string& strPath) {
    Logger* pLogger = Logger::getLogger();

    string strFilePathName;

    if(strPath.empty()) {
        QDir confPath = QDir::home();
        if(!confPath.exists(CONF_DIR)) {
            if(!confPath.mkdir(CONF_DIR)) {
                pLogger->log(QSNP_ERROR,"could not create config directory: " + confPath.path().toStdString() + CONF_DIR);
                return false;
            }
        }
        strFilePathName = confPath.path().toStdString() + "/" + CONF_DIR
                + "/" + strFileName;
    } else {
        strFilePathName = strPath + "/" + strFileName;
    }

    // configuration files should have the "*.cfg" extension
    if(strFilePathName.substr(strFilePathName.length() - 4, 4) != ".cfg") {
        strFilePathName += ".cfg";
    }

    ofstream	confFile;

    confFile.open(strFilePathName.c_str(),ios::out);
    if (confFile.fail()) {
        pLogger->log(QSNP_ERROR,"could not create config file: " + strFilePathName);
        return false;
    }

    confFile << toString();

    confFile.close();

    return true;
}

// advisable to do this on a new configuration object
bool Configuration::readFile(const string& strFileName, const string& strPath) {
    Logger* pLogger = Logger::getLogger();

    string strFilePathName;

    if(strPath.empty()) {
        QDir confDir = QDir::home();
        if(confDir.exists(CONF_DIR)) {
            strFilePathName = confDir.path().toStdString() + "/" + CONF_DIR
                    + "/" + strFileName;
        } else {
            pLogger->log(QSNP_ERROR, "Could not open config path: " + confDir.path().toStdString() + CONF_DIR);
            return false;
        }
    } else {
        strFilePathName = strPath + "/" + strFileName;
    }

    // configuration files should have the "*.cfg" extension
    if(strFilePathName.substr(strFilePathName.length() - 4, 4) != ".cfg") {
        strFilePathName += ".cfg";
    }

    ifstream	confFile;
    confFile.open(strFilePathName.c_str(),ios::in);
    if (confFile.fail()) {
        pLogger->log(QSNP_ERROR,"could not open config file: " + strFilePathName);
        return false;
    }

    char line[CONFIG_MAX_LINE_LENGTH];

    while(confFile.getline(line,CONFIG_MAX_LINE_LENGTH)) {
        stringstream stream(line);
        string label, value;
        getline(stream, label, '=');
        getline(stream, value);
        label = trim(label);
        value = trim(value);
        if(label.length() > 0 && !setValue(label,value)) {
            Logger::getLogger()->log(QSNP_ERROR, "could not set parameter " + label + " to value: " + value);
        }
    }

    confFile.close();

    return true;
}

const string Configuration::toString() {
    stringstream result;
    map<string,string>::iterator	strIt;
    map<string,double>::iterator	doubleIt;
    map<string,int>::iterator		intIt;
    map<string,char>::iterator		charIt;
    map<string,bool>::iterator		boolIt;

    result << "name" << "=" << _name << endl;

    for (strIt = stringMap.begin(); strIt != stringMap.end(); strIt++) {
        result << (*strIt).first << "=" << (*strIt).second << endl;
    }

    for (doubleIt = doubleMap.begin(); doubleIt != doubleMap.end(); doubleIt++) {
        result << (*doubleIt).first << "=" << (*doubleIt).second << endl;
    }

    for (intIt = intMap.begin(); intIt != intMap.end(); intIt++) {
        result << (*intIt).first << "=" << (*intIt).second << endl;
    }

    for (charIt = charMap.begin(); charIt != charMap.end(); charIt++) {
        result << (*charIt).first << "=" << (*charIt).second << endl;
    }

    for (boolIt =boolMap.begin(); boolIt != boolMap.end(); boolIt++) {
        result << (*boolIt).first << "=" << (*boolIt).second << endl;
    }

    return result.str();
}

list<string> Configuration::getSettings()
{
    list<string> settings;

    map<string,string>::iterator	strIt;

    for (strIt = stringMap.begin(); strIt != stringMap.end(); strIt++) {
        settings.push_back((*strIt).first);
    }

    map<string,double>::iterator	doubleIt;

    for (doubleIt = doubleMap.begin(); doubleIt != doubleMap.end(); doubleIt++) {
        settings.push_back((*doubleIt).first);
    }

    map<string,int>::iterator		intIt;

    for (intIt = intMap.begin(); intIt != intMap.end(); intIt++) {
        settings.push_back((*intIt).first);
    }

    map<string,char>::iterator		charIt;

    for (charIt = charMap.begin(); charIt != charMap.end(); charIt++) {
        settings.push_back((*charIt).first);
    }

    map<string,bool>::iterator		boolIt;

    for (boolIt = boolMap.begin(); boolIt != boolMap.end(); boolIt++) {
        settings.push_back((*boolIt).first);
    }

    return settings;
}

