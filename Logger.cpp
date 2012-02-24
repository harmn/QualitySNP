#include "Logger.h"

Logger* Logger::s_pTheLogger = 0;

Logger::Logger(void)
{
	Configuration* pConfig = Configuration::getConfig();
	_logFileName = pConfig->getString("outputDirectory") + "/" + pConfig->getString("logFileName");
	_logLevels.push_back("Always");
	_logLevels.push_back("Error");
	_logLevels.push_back("Warning");
	_logLevels.push_back("Info");
}

Logger::~Logger(void)
{
	_logFile.flush();
}

bool Logger::openFile() {
	_logFile.open(_logFileName.c_str(),ios::out);
	if (_logFile.fail()) {
		cerr << "Could not open log file" << _logFileName << NEWLINE;
		return false;
	}

	_logFile << "Logger started" << endl;

	_level = Configuration::getConfig()->getInt("logLevel");
    if (_level > _logLevels.size() - 1){
        _level = _logLevels.size() - 1;
    }
	_logFile << "Log level: " << _logLevels[_level] << endl;

    return true;
}

void Logger::log(int level, const string& message) {
	if(_level < level) {
		return;
	}

	if(_logFile.is_open()) {
		_logFile << _logLevels[level] << ": " << message << endl;
	} else {
		cerr << _logLevels[level] << ": " << message << endl;
	}
}
