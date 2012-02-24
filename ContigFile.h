#ifndef __CONTIGFILE_H__
#define __CONTIGFILE_H__

class ContigFile
{
public:
	virtual Contig* nextContig() = 0;
	virtual bool isValid() = 0;
	virtual ~ContigFile() {}
};

#endif
