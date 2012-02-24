#include <list>
#include <sstream>
#include <cstdlib>
#include "Logger.h"
#include "SAMRead.h"

using namespace std;

SAMRead::~SAMRead(void)
{
}

void SAMRead::setSequence(const string& sequence) {
	_sequence = sequence;
}

void SAMRead::setQuality(const string& quality) {
	_quality = quality;

	if(_quality == "*") {
		// when there is no quality information the quality score for each nucleotide is set to the highest
		string qual(_sequence.size(), '~'); //~ is the highest quality score
		_quality = qual;
	}
}

bool SAMRead::processOperations() {
	list<operation>::iterator itOps;

	_opsequence.assign(_sequence.size(),' ');

	//remove hard clip operations from front and back
	if (_ops.front().op == 'H') {
		_ops.pop_front();
	}
	if(_ops.back().op == 'H') {
		_ops.pop_back();
	}

	//trim soft clipped ends
	if (_ops.front().op == 'S') {
		if (_sequence.size() < _ops.front().size) {
			Logger::getLogger()->log(QSNP_ERROR, "Sequence shorter than expected from cigar operations for read:  " + _name);
			return false;
		}
		_sequence = _sequence.substr(_ops.front().size);
		_ops.pop_front();
	}
	if(_ops.back().op == 'S') {
		int newSize = _sequence.size() - _ops.back().size;
		if (newSize < 0) {
			Logger::getLogger()->log(QSNP_ERROR, "Sequence shorter than expected from cigar operations for read:  " + _name);
			return false;
		}

		_sequence.resize(newSize);
		_ops.pop_back();
	}

	int nPosition = 0;
	for ( itOps=_ops.begin() ; itOps != _ops.end(); itOps++ ) {
		switch ((*itOps).op) {
			case 'S':
				Logger::getLogger()->log(QSNP_ERROR, "Internal S in SAM cigar line for read: " + _name);
				return false;
				break;
			case 'M':
			case 'I':
				for(int i = nPosition; i < nPosition + (*itOps).size; i++) {
					_opsequence[i] = (*itOps).op;
				}
				break;
			case 'P':
			case 'D':
			case 'N':
				_sequence.insert(nPosition, (*itOps).size,'*');
				_opsequence.insert(nPosition, (*itOps).size,(*itOps).op);
				_quality.insert(nPosition, (*itOps).size,'~'); //~ is the highest quality score
				break;
			case '=':
			case 'X':
				for(int i = nPosition; i < nPosition + (*itOps).size; i++) {
					_opsequence[i] = (*itOps).op;
				}
			default:
				break;
		} 
		nPosition += (*itOps).size;
	}

	return true;
}

char SAMRead::getNucleotideAt(unsigned int pos) {
	pos--;
	if(_startPosition >= 0 && pos < static_cast<unsigned int>(_startPosition)) {
		return ' ';
	}

	pos -= _startPosition;
	return (pos > _sequence.size() - 1) ? ' ' : _sequence[pos];
}

char SAMRead::getOperationAt(unsigned int pos) {
	pos--;
	if(_startPosition >= 0 && pos < static_cast<unsigned int>(_startPosition)) {
		return '<';
	}

	pos -= _startPosition;
	return (pos > _opsequence.size() - 1) ? '>' : _opsequence[pos];
}

void SAMRead::insertGapAt(unsigned int pos) {
	pos--;
	if(_startPosition >= 0 && pos < static_cast<unsigned int>(_startPosition)) {
		_startPosition++;
		return;
	}

	char chGap = '*';

	pos -= _startPosition;
	if (pos < _sequence.size() && _opsequence[pos] != 'I') {
		_sequence.insert(pos, 1, chGap);
		_opsequence.insert(pos, 1, 'I');
		_quality.insert(pos, 1, '~'); // ~ is the highest quality score
	}
}

void SAMRead::merge(SAMRead* pRead) {
	SAMRead* pRead1;
	SAMRead* pRead2;

	if(getStartPosition() < pRead->getStartPosition()) {
		pRead1 = this;
		pRead2 = pRead;
	} else {
		pRead1 = pRead;
		pRead2 = this;
	}

	stringstream ss;
	stringstream qual;

	int spacing = pRead2->getStartPosition() - pRead1->getStartPosition() - pRead1->getSequence().size();
	if(spacing < 0) {
		// overlapping sequences
		int limit = (abs(spacing) > pRead2->getQuality().size() ) ? pRead2->getQuality().size() + spacing : 0;

		for(int i = spacing; i < limit; i++) {
			if(pRead2->getQuality()[i - spacing] < pRead1->getQuality()[pRead1->getQuality().size() + i]) {
				pRead2->_sequence[i - spacing] = pRead1->getSequence()[pRead1->getSequence().size() + i];
				pRead2->_quality[i - spacing] = pRead1->getQuality()[pRead1->getQuality().size() + i];
			}
		}
		pRead1->setSequence(pRead1->getSequence().substr(0,pRead1->getSequence().size() + spacing));
		pRead1->setQuality(pRead1->getQuality().substr(0,pRead1->getQuality().size() + spacing));

		ss << pRead1->getSequence() << pRead2->getSequence();
		qual << pRead1->getQuality() << pRead2->getQuality();
	} else {
		ss << pRead1->getSequence() << string(spacing, ' ') << pRead2->getSequence();
		qual << pRead1->getQuality() << string(spacing, '~') << pRead2->getQuality();
	}

	_startPosition = pRead1->getStartPosition();
	_sequence = ss.str();
	_quality= qual.str();
}

SeqRead* SAMRead::toSeqRead() {
	SeqRead* pRead = new SeqRead(_name);
	pRead->setSequence(_sequence);
	pRead->setStartPosition(_startPosition + 1); // in SeqRead this will be corrected to zero-based
	pRead->setQualitySanger(_quality);
	return pRead;
}

