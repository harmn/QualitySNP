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

#ifndef __CONTIGFILE_H__
#define __CONTIGFILE_H__

class ContigFile
{
public:
	virtual Contig* nextContig() = 0;
	virtual bool isValid() = 0;
	virtual ~ContigFile() {}
    virtual int contigCount() = 0;
    virtual int readCount() = 0;
    virtual vector<string>& readGroups() = 0;
};

#endif
