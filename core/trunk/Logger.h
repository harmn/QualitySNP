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

#ifndef __LOG_H__
#define __LOG_H__

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <QMutex>
#include "Configuration.h"

#define QSNP_ALWAYS 0
#define QSNP_ERROR 1
#define QSNP_WARNING 2
#define QSNP_INFO 3
#define QSNP_DEBUG 4

using namespace std;

class Logger
{
public:
	~Logger(void);
	bool openFile();
	void log(int level, const string& message);


	static Logger* getLogger() {
		if (s_pTheLogger == 0) {
			s_pTheLogger = new Logger();
			s_pTheLogger->openFile();
		}
		return s_pTheLogger;
	}

private:
	Logger(void);
	string			_logFileName;
	ofstream		_logFile;
	vector<string>	_logLevels;
	int				_level;
	static Logger*	s_pTheLogger;
    QMutex          _mutex;
};

#endif
