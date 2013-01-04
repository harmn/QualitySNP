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

#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

#include <string>
#include <map>
#include <list>
#include "Logger.h"

#define NEWLINE '\n';
static const int CONFIG_MAX_LINE_LENGTH = 4000;
static const char* CONF_DIR = ".QualitySNPng";

using namespace std;


class Configuration
{
public:
	~Configuration(void);
    const string getName()                    { return _name; }
    void setName(string name)           { _name = name; }
	int getInt(const string& key)		{ return intMap[key]; }
	int getChar(const string& key)		{ return charMap[key]; }
    map<string, string> findConfigurations();
	double getDouble(const string& key) { return doubleMap[key]; }
	const string& getString(const string& key) { return stringMap[key]; }
	bool getBool(const string& key)		{ return boolMap[key]; 	}
    const string getValue(const string& key);
	bool setValue(const string& key, const string& value);
	bool setInt(const string& key, int value); 
	bool setChar(const string& key, char value); 
	bool setDouble(const string& key, double value);
	bool setString(const string& key, const string& value); 
	bool setBool(const string& key, bool value);
	int getNumberOfNucs()				{ return rnucMap.size(); }
	int nuc2int(char nuc)				{ return nucMap[nuc]; }

	char int2nuc(int nuc) { 
        return (nuc >= 0 && nuc < (int) rnucMap.size()) ? rnucMap[nuc] : 'N';
	}

	char getIUPACCode(int nucleotides) {
        return (IUPACCode.count(nucleotides) > 0) ? IUPACCode[nucleotides] : 'N';
	}

    bool readFile(const string& strFileName, const string& strPath);
    bool writeFile(const string& strFileName, const string& strPath);

    static Configuration* getConfig(bool bNew=false) {
        if(bNew) {
             delete s_pTheConfiguration;
            s_pTheConfiguration = 0;
        }

        if (s_pTheConfiguration == 0) {
			s_pTheConfiguration = new Configuration();
        }
		return s_pTheConfiguration;
	}

	const string toString();

    list<string>        getSettings();

private:
	Configuration(void);
    string                  _name;
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
