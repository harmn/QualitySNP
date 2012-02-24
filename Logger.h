#ifndef __LOG_H__
#define __LOG_H__

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "Configuration.h"

#define QSNP_ALWAYS 0
#define QSNP_ERROR 1
#define QSNP_WARNING 2
#define QSNP_INFO 3

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
};

#endif
