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

#include <QApplication>
#include "Logger.h"

Logger* Logger::s_pTheLogger = 0;

Logger::Logger(void)
{
    _level = QSNP_ALWAYS;
	Configuration* pConfig = Configuration::getConfig();

    _logFileName = pConfig->getString("outputDirectory") + "/" + pConfig->getString("logFileName");
	_logLevels.push_back("Always");
	_logLevels.push_back("Error");
	_logLevels.push_back("Warning");
	_logLevels.push_back("Info");
	_logLevels.push_back("Debug");
}

Logger::~Logger(void)
{
	_logFile.flush();
}

bool Logger::openFile() {
	_logFile.open(_logFileName.c_str(),ios::out);
    bool bSuccess = true;

	if (_logFile.fail()) {
        log(QSNP_ERROR, "Could not create log file: " + _logFileName);
        bSuccess = false;
	}

    log(QSNP_ALWAYS, "Logger started");

	_level = Configuration::getConfig()->getInt("logLevel");
    if (_level > _logLevels.size() - 1){
        _level = _logLevels.size() - 1;
    }
    log(QSNP_ALWAYS,"Log level: " + _logLevels[_level]);

    return bSuccess;
}

void Logger::log(int level, const string& message) {
	if(_level < level) {
		return;
	}

    _mutex.lock();
	if(_logFile.is_open()) {
		_logFile << _logLevels[level] << ": " << message << endl;
	} else {
		cerr << _logLevels[level] << ": " << message << endl;
	}
    _mutex.unlock();
}
