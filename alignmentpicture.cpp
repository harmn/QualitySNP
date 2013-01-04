#include <QMap>
#include <qmath.h>
#include <QAction>
#include <QMouseEvent>
#include <QApplication>
#include <QTime>
#include "core/trunk/Variation.h"
#include "core/trunk/SeqRead.h"
#include "alignmentpicture.h"

bool sortByStartPositionFunction(SeqRead* s1, SeqRead* s2) {
    return (s1->getStartPosition() < s2->getStartPosition());
}

bool sortByHaploTypeFunction(SeqRead* s1, SeqRead* s2) {
    HaploType* h1 = s1->getHaploType();
    HaploType* h2 = s2->getHaploType();

    if(h1 == NULL && h2 == NULL) {
        return (s1->getStartPosition() < s2->getStartPosition());
    } else if(h1 == NULL) {
        return false;
    } else if(h2 == NULL) {
        return true;
    } else if (s1->getHaploType()->getID() == s2->getHaploType()->getID()) {
        return (s1->getStartPosition() < s2->getStartPosition());
    } else {
        return s1->getHaploType()->getID() < s2->getHaploType()->getID();
    }
}

AlignmentPicture::AlignmentPicture(Contig* pContig, QWidget *parent ) :
    QWidget(parent), _pict(NULL), _overviewPict(NULL), _bPixelMode(true),
    _nW(1), _nH(1), _x(0), _y(0), _prevX(0), _prevY(0), _pContig(pContig),
    _coveragePlotHeight(50), _coveragePlot(NULL), _zoomX(1.0), _zoomY(1.0),
    _prevZoomX(1.0), _prevZoomY(1.0)
{
    initColorMap();
    calculateTextDimensions();

    _sortAction = new QAction(QString(), this);
    connect(_sortAction, SIGNAL(triggered()), this, SLOT(sortReads()));
    addAction(_sortAction);
    setContextMenuPolicy(Qt::ActionsContextMenu);

    QAction* fullSizeAction = new QAction(tr("show overview"), this);
    connect(fullSizeAction, SIGNAL(triggered()), this, SLOT(showOverview()));
    addAction(fullSizeAction);

    _reads = pContig->getReads();

    _bSortbyHaploType = false;
    _cRows = pContig->getReadCount();
    _cCols = pContig->getSequenceLength();
    sortReads();

    QPalette p;
    p.setColor(QPalette::Window, Qt::white);
    setPalette(p);
    setAutoFillBackground(true);
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);

    drawCoveragePlot();
}

const double AlignmentPicture::ZOOM_FACTOR = 0.8;

AlignmentPicture::~AlignmentPicture()
{
    delete _pContig;
}

// setting the colors for the different nucleotides
void AlignmentPicture::initColorMap() {
    _colorMap["*"] = QColor(255, 255, 255 ); // white
    _colorMap["N"] = QColor(255, 255, 255 ); // white
    _colorMap["A"] = QColor(255, 0, 0, 127);
    _colorMap["C"] = QColor(0, 255, 0, 127);
    _colorMap["G"] = QColor(255, 255, 0, 127);
    _colorMap["T"] = QColor(0, 0, 255, 127);
}

void AlignmentPicture::calculateTextDimensions() {
    QPainter painter;

    QPixmap mp(100,100);
    painter.begin(&mp);
    _charRect = painter.boundingRect(rect(), "G");
    painter.end();
    _charWidth =  _charRect.width() + 1.0;
    _charHeight =  _charRect.height() + 1.0;
    updateCharDimensions();
}

void AlignmentPicture::updateCharDimensions() {
    _nW = int(_charWidth*_zoomX) + 1;
    _nH = int(_charHeight*_zoomX) + 1;

    if(_nW < 3) {
        _bPixelMode = true;
        // in pixelmode all nucleotides are 1 pixel in size
        _nW = 1;
        _nH = 1;
    } else {
        _bPixelMode = false;
    }
}

void AlignmentPicture::paintEvent(QPaintEvent *)
{
    updateCharDimensions();

    _painter.begin(this);
    _painter.setBackground(Qt::white);

    if(_x > 0) {
        _x = 0;
    }

    if(_y > 0) {
        _y = 0;
    }

    if(_bPixelMode) {
        drawOverview();

        qreal scaleX = _zoomX * _charWidth;
        qreal scaleY = _zoomY * _charWidth;

        _painter.save();
        _painter.translate(_x,_y);
        _painter.scale(scaleX, scaleY);
        _painter.drawPixmap(0, _coveragePlotHeight/scaleY, *_overviewPict);
        _painter.restore();
        _painter.translate(_x,0);
        _painter.scale(scaleX, 1.0);
        _painter.drawPixmap(0, 0, *_coveragePlot);
    } else {
        drawRegion(_x, _y, width(),height());

        _painter.drawPixmap(0, _coveragePlotHeight, *_pict);

        qreal scaleX = int(_zoomX * _charWidth) + 1;
        _painter.translate(_x,0);
        _painter.scale(scaleX, 1.0);
        _painter.drawPixmap(0, 0, *_coveragePlot);
    }
    _painter.end();
}

void AlignmentPicture::mouseDoubleClickEvent(QMouseEvent *event)
{
    int x = event->x();
    int y = event->y();

    zoomIn(x, y);
}

void AlignmentPicture::mousePressEvent(QMouseEvent *event)
{
    _prevX = event->x();
    _prevY = event->y();
}

void AlignmentPicture::mouseReleaseEvent(QMouseEvent *event)
{
}

void AlignmentPicture::mouseMoveEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::LeftButton) {
        int movedX = event->x() - _prevX;
        int movedY = event->y() - _prevY;

        _prevX = event->x();
        _prevY = event->y();

        _x += movedX;
        _y += movedY;

        if(_x > 0)
            _x = 0;
        if(_y > 0)
            _y = 0;

        update();
    }

    int x = event->x() - _x;
    int y = event->y() - _y - _coveragePlotHeight;
    int col = 0;
    int row = 0;

    if (_bPixelMode) {
        col = x/(_zoomX*_charWidth) + 1;
        row = y/(_zoomY*_charWidth);
        if(_cRows > MAX_SCREEN_HEIGHT) {
            row *= 1.0 * _cRows/MAX_SCREEN_HEIGHT;
        }
    } else {
        col = int(x/_nW) + 1;
        row = int(y/_nH);
    }
    QString strToolTip;
    if(event->y() < 50) {
        if(col > 0 && col <= _coverage.size()) {
            strToolTip = QString("pos: %1, coverage: %2").arg(col).arg(_coverage[col - 1]);
        }
    } else if(row < _cRows && row >= 0) {
        ReadInfo read = _readInfo[row];
        if(col > read.start && col <= read.end) {
            if(read.haplotypeID == -1) {
                strToolTip = QString("pos: %1\nread: %2\nhaplotype: none").arg(col).arg(read.name);
            } else {
                strToolTip = QString("pos: %1\nread: %2\nhaplotype: %3").arg(col).arg(read.name).arg(read.haplotypeID);
            }
            if(!read.readgroup.isEmpty()) {
                strToolTip += QString("\ngroup: %1").arg(read.readgroup);
            }
        } else {
            strToolTip = QString("pos: %1, row: %2").arg(col).arg(row + 1);
        }
    } else {
        strToolTip = QString("pos: %1, row: %2").arg(col).arg(row + 1);
    }

    setToolTip(strToolTip);
}

void AlignmentPicture::keyPressEvent(QKeyEvent *event)
{
    int moveStep = 100;

    if(event->key() == Qt::Key_Minus) {
        zoomOut();

    } else if(event->key() == Qt::Key_Equal || event->key() == Qt::Key_Plus) {
        zoomIn(width()/2, height()/2);

    } else if(event->key() == Qt::Key_Z) {
        zoomInX();
    } else if(event->key() == Qt::Key_X) {
        zoomOutX();
    } else if(event->key() == Qt::Key_A) {
        zoomInY();
    } else if(event->key() == Qt::Key_S) {
        zoomOutY();
    } else if(event->key() == Qt::Key_Left) {
        _x += moveStep;
        if(_x > 0) // moving rightways is limited
            _x = 0;
        update();
    } else if(event->key() == Qt::Key_Right) {
        _x -= moveStep;
        update();
    } else if(event->key() == Qt::Key_Up) {
        _y += moveStep;
        if(_y > 0)
            _y = 0;
        update();
    } else if(event->key() == Qt::Key_Down) {
        _y -= moveStep;
        update();
    } else if(event->key() == Qt::Key_PageUp) {
        _y += height();
        if(_y > 0)
            _y = 0;
        update();
    } else if(event->key() == Qt::Key_PageDown) {
        _y -= height();
        update();
    } else if(event->key() == Qt::Key_F) {
        zoomToFit();
        update();
    } else if(event->key() == Qt::Key_H) {
        _x = 0;
        _y = 0;
        update();
    } else if(event->key() == Qt::Key_C) {
        _x = _cCols / 2;
        _y = _cRows / 2;
        update();
    }
}

void AlignmentPicture::zoomIn(int x, int y) {
    if ( _zoomX < 1.0) {
        _zoomX /= ZOOM_FACTOR;
    }

    if ( _zoomY < 1.0) {
        _zoomY /= ZOOM_FACTOR;
    }

    rescale(x,y);
}

void AlignmentPicture::zoomInX() {
    if ( _zoomX > 1.0) {
        return;
    }

    _zoomX /= ZOOM_FACTOR;
    rescale(width()/2, height()/2);
}

void AlignmentPicture::zoomInY() {
    if ( _zoomY > 1.0) {
        return;
    }

    _zoomY /= ZOOM_FACTOR;
    rescale(width()/2, height()/2);
}

void AlignmentPicture::zoomOut() {
    _zoomX *= ZOOM_FACTOR;
    _zoomY *= ZOOM_FACTOR;
    rescale(width()/2, height()/2);
}

// horizontal stretch
void AlignmentPicture::zoomOutX() {
    _zoomX *= ZOOM_FACTOR;
    rescale(width()/2, height()/2);
}

// vertical stretch
void AlignmentPicture::zoomOutY() {
    _zoomY *= ZOOM_FACTOR;
    rescale(width()/2, height()/2);
}


void AlignmentPicture::rescale(int x, int y) {
    qreal cX = 0.0;
    qreal cY = 0.0;

    if(_bPixelMode) {
        qreal scaleX = _prevZoomX * _charWidth;
        qreal scaleY = _prevZoomY * _charWidth;
        cX = (x - _x)/scaleX;
        cY = (y - _y - _coveragePlotHeight)/scaleY;
        if(_cRows > MAX_SCREEN_HEIGHT) {
            cY *= 1.0 * _cRows/MAX_SCREEN_HEIGHT;
        }
    } else {
        cX = (x -_x)/_nW + 1;
        cY = (y -_y - _coveragePlotHeight)/_nH;
    }

    // after a zoom operation the character height and width are recalculated
    updateCharDimensions();

    if(_bPixelMode) {
        qreal scaleX = _zoomX * _charWidth;
        qreal scaleY = _zoomY * _charWidth;

        if(_cRows > MAX_SCREEN_HEIGHT) {
            cY /= 1.0 * _cRows/MAX_SCREEN_HEIGHT;
        }

        _x = x - cX*scaleX;
        _y = y - cY*scaleY - _coveragePlotHeight;
    } else {
        _x = x -cX*_nW;
        _y = y- cY*_nH - _coveragePlotHeight;
    }


    update();

    _prevZoomX = _zoomX;
    _prevZoomY = _zoomY;
}

// load the sequence information into the datastructures for visualization
void AlignmentPicture::loadSequences() {
    _coverage.resize(_pContig->getSequenceLength());
    _readInfo.clear();
    _readInfo.reserve(_cRows);
    for(list<SeqRead*>::const_iterator itReads = _reads.begin(); itReads != _reads.end(); itReads++) {
        string sequence =  (*itReads)->getSequence();
        int startPos =  (*itReads)->getStartPosition();
        if(startPos < 0) {
            sequence = sequence.substr(-startPos);
            startPos = 0;
        }
        QString qsequence(startPos, QChar(' '));
        qsequence += QString(QString::fromUtf8(sequence.data(), sequence.size()));

        ReadInfo readInfo;
        readInfo.name = QString::fromStdString((*itReads)->getName());
        readInfo.start = startPos;
        readInfo.end = qsequence.length();
        readInfo.readgroup = QString::fromStdString((*itReads)->getGroup());

        // record the coverage per position
        for(int i = readInfo.start; i < readInfo.end; i ++) {
            _coverage[i]++;
        }
        readInfo.sequence = qsequence;
        HaploType* pHap = (*itReads)->getHaploType();
        readInfo.haplotypeID = (pHap != NULL) ? pHap->getID() : -1;
        _readInfo.push_back(readInfo);
    }

    // normalize the coverage
    int i = 0;
    int maxCoverage = 0;
    for(i = 0; i < _coverage.size(); i++) {
        maxCoverage = (maxCoverage < _coverage[i]) ? _coverage[i] : maxCoverage;
    }

    _coverageNorm.resize(_coverage.size());
    for(i = 0; i < _coverage.size(); i++) {
        _coverageNorm[i] = 1.0 * _coverage[i] / maxCoverage * _coveragePlotHeight + 1.0;
    }
}

void AlignmentPicture::drawCoveragePlot() {
    if(_coveragePlot != NULL) {
        return;
    }

    _coveragePlot = new QPixmap(_cCols, _coveragePlotHeight);
    _coveragePlot->fill(Qt::white);

    QPainter painter;
    painter.begin(_coveragePlot);

    QVector<QLine> coverageLines(_coverageNorm.size());
    for(int iCoverage = 0; iCoverage < _coverageNorm.size(); iCoverage++) {
        coverageLines[iCoverage] = QLine(iCoverage,_coveragePlotHeight,iCoverage,_coveragePlotHeight -_coverageNorm[iCoverage]);
    }

    QPen penCoverage(Qt::gray);
    penCoverage.setWidth(1);
    painter.setPen(penCoverage);
    painter.drawLines(coverageLines);

    painter.end();
}

void AlignmentPicture::drawOverview() {
    if(_overviewPict != NULL) {
        // overview picture was already drawn, no need to do that again
        return;
    }

#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif

    QMap<QString, QVector<QPoint> > nucMap;

    nucMap[tr("A")];
    nucMap[tr("C")];
    nucMap[tr("G")];
    nucMap[tr("T")];

    double increment = 1.0;
    int height = _cRows;

    // limit the screen height
    if(_cRows > MAX_SCREEN_HEIGHT) {
        increment = 1.0 * _cRows/MAX_SCREEN_HEIGHT;
        height = MAX_SCREEN_HEIGHT;
    }

    _overviewPict = new QPixmap(_cCols, height);
    _overviewPict->fill(Qt::white);

    QPainter painter;
    painter.begin(_overviewPict);

    int y = 0;
    for(double iRead = 0; iRead < _cRows; iRead+= increment) {
        QString qsequence = _readInfo[iRead].sequence;
        QStringList seqList = qsequence.split("", QString::SkipEmptyParts);

        for (int iSeq = 0; iSeq < seqList.size(); iSeq++) {
            if(_colorMap.contains(seqList[iSeq])) {
                nucMap[seqList[iSeq]].append(QPoint(iSeq,y));
            }
        }
        y++;
    }

    QMap<QString, QVector<QPoint> >::iterator itNucMap;
    for(itNucMap = nucMap.begin(); itNucMap != nucMap.end(); itNucMap++) {
        painter.setPen(_colorMap[itNucMap.key()]);
        painter.drawPoints(itNucMap.value().data(), itNucMap.value().size());
    }

    QPen pen(Qt::black);
    pen.setWidth(3);
    painter.setPen(pen);
    vector<Variation*> variations = _pContig->getVariations();
    for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
        if((*itVar)->isReliable()) {
            int pos = (*itVar)->getPos();
            painter.drawLine(pos,0, pos, height);
        }
    }

    painter.end();

#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
}

void AlignmentPicture::drawRegion(int x, int y, int width, int height)
{
    QPainter painter;
    QMap<QString, QPixmap> nucMap;

    QMap<QString,QColor>::const_iterator itColMap;
    for(itColMap = _colorMap.constBegin(); itColMap != _colorMap.constEnd(); itColMap++) {
        nucMap[itColMap.key()] = QPixmap(_nW,_nH);
        nucMap[itColMap.key()].fill(_colorMap[itColMap.key()]);
        painter.begin(&nucMap[itColMap.key()]);
        painter.scale(_zoomX,_zoomX);
        painter.setBackground(_colorMap[itColMap.key()]);
        painter.drawText(_charRect, itColMap.key());
        painter.end();
    }

    if(_pict != NULL) {
        delete _pict;
    }

    _pict = new QPixmap(width, height);
    _pict->fill(Qt::white);

    painter.begin(_pict);

    int startRead = -y/_nH;
    if(startRead < 0) {
        startRead = 0;
    }
    int endRead = height/_nH + startRead;
    if(endRead >= _readInfo.size()) {
        endRead = _readInfo.size();
    }

    int yPos = 0;

    for(int iRead = startRead; iRead < endRead; iRead++) {
        QString qsequence = _readInfo[iRead].sequence;
        QStringList seqList = qsequence.split("", QString::SkipEmptyParts);

        int startPos = -x/_nW;
        if(startPos < 0) {
            startPos = 0;
        }
        int endPos = width/_nW + startPos;
        if(endPos >= seqList.size()) {
            endPos = seqList.size();
        }
        int xPos = 0;
        for (int iSeq = startPos; iSeq < endPos; iSeq++) {
            if(seqList[iSeq] != " ") {
                if(nucMap.contains(seqList[iSeq])) {
                    painter.drawPixmap(xPos*_nW, yPos*_nH,nucMap[seqList[iSeq]]);
                } else {
                    painter.drawText(xPos*_nW, yPos*_nH + _nH * 0.8,seqList[iSeq]);
                }
            }
            xPos++;
        }
        yPos++;
    }

    int startPos = -x/_nW;
    int endPos = width/_nW + startPos;
    vector<Variation*> variations = _pContig->getVariations();
    for(vector<Variation*>::iterator itVar = variations.begin(); itVar != variations.end(); itVar++) {
        int pos = (*itVar)->getPos();
        if((*itVar)->isReliable() && pos >= startPos && pos <= endPos) {
            painter.setPen(Qt::black);
            QRect varRect((pos - startPos)*_nW, 0, _nW, _cRows * _nH);
            painter.setBrush(QColor(0, 0, 0, 75));
            painter.drawRect(varRect);
        }
    }

    painter.end();
}

void AlignmentPicture::sortReads()
{
    if(_bSortbyHaploType) {
        _reads.sort(sortByHaploTypeFunction);
        _sortAction->setText(tr("sort by start position"));
    } else {
        _reads.sort(sortByStartPositionFunction);
        _sortAction->setText(tr("sort by haplotype"));
    }
    loadSequences();
    if(_overviewPict != NULL) {
        delete _overviewPict;
        _overviewPict = NULL;
    }

    update();
    _bSortbyHaploType = !_bSortbyHaploType;
}

void AlignmentPicture::zoomToFit()
{
    int overviewHeight = (_cRows < MAX_SCREEN_HEIGHT) ? _cRows : MAX_SCREEN_HEIGHT;
    _zoomX = 1.0*width()/_cCols/_charWidth;
    _zoomY = 1.0*(height() - _coveragePlotHeight)/overviewHeight/_charWidth;

    _x = 0;
    _y = 0;
}

void AlignmentPicture::showOverview()
{
    zoomToFit();
    update();
}

void AlignmentPicture::scrollToPosition(int pos)
{
    int iRow = 0;
    while(iRow < _cRows && _readInfo[iRow].end < pos) {
        iRow++;
    }

    _y = 0;
    if(_bPixelMode) {
        _x = -(pos - width() / 2);
        _y = -iRow*_zoomY*_charWidth;
    } else {
        _x = -(pos*_nW - width() / 2);
        _y = -iRow*_nH;
    }

    update();
}
