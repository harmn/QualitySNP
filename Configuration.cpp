#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include "Configuration.h"

Configuration* Configuration::s_pTheConfiguration = 0;

typedef enum dataType {
	INT,
	DOUBLE,
	BOOL,
	CHAR,
	STRING
} DATATYPE;

std::string trim(std::string& s,const std::string& drop = " \t\"'")
{
	std::string r=s.erase(s.find_last_not_of(drop)+1);
	return r.erase(0,r.find_first_not_of(drop));
}

Configuration::Configuration(void)
{
	intMap["minimalNumberOfReadsPerAllele"] = 2;
	doubleMap["minimalNumberOfReadsPerAllelep"] = 0.0;
	intMap["lowQualityRegion5prime"]		= 0;
	intMap["lowQualityRegion3prime"]		= 0;
	intMap["maxNumberOfSNPsInFlanks"]		= 0;
	doubleMap["lowQualityRegion3primePerc"]	= 0.0;
	doubleMap["weightHighQualityRegion"]	= 1.0;
	doubleMap["weightLowQualityRegion"]		= 0.5;
	intMap["minSNPQualityScore"]			= 20;
	intMap["minimalConfidenceScore"]		= 2;
	intMap["minNumberOfHighQualityReads"]	= 2;
	intMap["logLevel"]				= QSNP_WARNING;
	doubleMap["similarityPerPolymorphicSite"]	= 0.75;
	doubleMap["similarityAllPolymorphicSites"]	= 0.8;
	stringMap["configurationFile"]	= "";
    stringMap["contigFileName"]		= "";
	stringMap["outputDirectory"]	= "/tmp";
	stringMap["logFileName"]		= "QSNP.log";
	stringMap["contigsFile"]		= "contigs.csv";
	stringMap["readsFile"]			= "reads.csv";
	stringMap["haploTypesFile"]		= "haplotypes.csv";
	stringMap["variationsFile"]		= "variations.csv";
	stringMap["FASTAFileName"]		= "";

	charMap["fieldSeparator"]		= '\t';
	boolMap["indelLowQuality"]		= true;
	boolMap["printAlignment"]		= false;
	boolMap["printHaploTypes"]		= false;
	boolMap["printSNPs"]			= false;
	boolMap["printMarkers"]			= false;
	boolMap["printSummaryLine"]		= true;
	boolMap["useIUPACCodes"]		= false;
	boolMap["onlyReliableMarkers"]	= true;

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

bool Configuration::setValue(const string& label, const string& value) {
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
			if(!(converter >> cValue)) {
				return false;
			}
		}
		setChar(label,cValue);
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
		bValue = (tmpValue == "true" || tmpValue == "t");
		setBool(label,bValue);
	} else {
		return false;
	}

	return true;
}

bool Configuration::readFile(const string& strFileName) {
	ifstream	confFile;
	confFile.open(strFileName.c_str(),ios::in);
	if (confFile.fail()) {
		Logger::getLogger()->log(QSNP_ERROR,"could not open config file: " + strFileName);
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
			Logger::getLogger()->log(QSNP_ERROR, "could not set parameter " + label + " to value + " + value);
			confFile.close();
			return false;
		}
	}
	
 	confFile.close();

    return true;
}

string Configuration::toString() {
	stringstream result;
	map<string,string>::iterator	strIt;
	map<string,double>::iterator	doubleIt;
	map<string,int>::iterator		intIt;
	map<string,char>::iterator		charIt;
	map<string,bool>::iterator		boolIt;

	for (strIt = stringMap.begin(); strIt != stringMap.end(); strIt++) {
		result << (*strIt).first << ": " << (*strIt).second << endl;
	}

	for (doubleIt = doubleMap.begin(); doubleIt != doubleMap.end(); doubleIt++) {
		result << (*doubleIt).first << ": " << (*doubleIt).second << endl;
	}

	for (intIt = intMap.begin(); intIt != intMap.end(); intIt++) {
		result << (*intIt).first << ": " << (*intIt).second << endl;
	}

	for (charIt = charMap.begin(); charIt != charMap.end(); charIt++) {
		result << (*charIt).first << ": " << (*charIt).second << endl;
	}

	for (boolIt =boolMap.begin(); boolIt != boolMap.end(); boolIt++) {
		result << (*boolIt).first << ": " << (*boolIt).second << endl;
	}

	return result.str();
}
