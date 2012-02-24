#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

#include <string>
#include <map>
#include "Logger.h"

#define NEWLINE '\n';
static const int CONFIG_MAX_LINE_LENGTH = 4000;

using namespace std;


class Configuration
{
public:
	~Configuration(void);
	int getInt(const string& key)		{ return intMap[key]; }
	int getChar(const string& key)		{ return charMap[key]; }
	double getDouble(const string& key) { return doubleMap[key]; }
	const string& getString(const string& key) { return stringMap[key]; }
	bool getBool(const string& key)		{ return boolMap[key]; 	}
	bool setValue(const string& key, const string& value);
	bool setInt(const string& key, int value); 
	bool setChar(const string& key, char value); 
	bool setDouble(const string& key, double value);
	bool setString(const string& key, const string& value); 
	bool setBool(const string& key, bool value);
	int getNumberOfNucs()				{ return rnucMap.size(); }
	int nuc2int(char nuc)				{ return nucMap[nuc]; }

	char int2nuc(int nuc) { 
		return (nuc >= 0 && nuc < (int) rnucMap.size()) ? rnucMap[nuc] : ' ';
	}

	char getIUPACCode(int nucleotides) {
		return (IUPACCode.count(nucleotides) > 0) ? IUPACCode[nucleotides] : ' ';
	}

	bool readFile(const string& strFileName);

	static Configuration* getConfig() {
		if (s_pTheConfiguration == 0) {
			s_pTheConfiguration = new Configuration();
		}
		return s_pTheConfiguration;
	}

	string toString();


private:
	Configuration(void);
	map<string,string>		stringMap;
	map<string,double>		doubleMap;
	map<string,int>			intMap;
	map<string,char>		charMap;
	map<string,bool>		boolMap;
	map<int,char>			IUPACCode;
	int*					nucMap;
	string					rnucMap;
	static Configuration*	s_pTheConfiguration;
};

#endif
