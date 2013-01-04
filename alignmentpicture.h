#ifndef ALIGNMENTPICTURE_H
#define ALIGNMENTPICTURE_H

#include <QWidget>
#include <QPainter>
#include <QPicture>
#include <QPixmap>
#include <QMap>
#include <list>
#include "core/trunk/Contig.h"

struct ReadInfo{
 QString name;
 int start;
 int end;
 int haplotypeID;
 QString sequence;
 QString readgroup;
};

class AlignmentPicture : public QWidget
{
    Q_OBJECT
public:
    explicit AlignmentPicture(Contig* pContig, QWidget *parent = 0);
    ~AlignmentPicture();
    void        scrollToPosition(int pos);
    void        zoomToFit();

protected:
    void        zoomIn(int x, int y);
    void        zoomOut();
    void        zoomInX();
    void        zoomInY();
    void        zoomOutX();
    void        zoomOutY();
    void        rescale(int x, int y);
    void        paintEvent( QPaintEvent * );
    void        mouseDoubleClickEvent(QMouseEvent *event);
    void        mousePressEvent ( QMouseEvent * event );
    void        mouseReleaseEvent (QMouseEvent* event );
    void        mouseMoveEvent ( QMouseEvent * event );
    void        keyPressEvent(QKeyEvent* event);

private:
    void drawRegion(int x, int y, int width, int height);
    void drawOverview();
    void drawCoveragePlot();
    void loadSequences();
    void initColorMap();
    void calculateTextDimensions();
    void updateCharDimensions();
    void calculateCoverage();
    static const int MAX_SCREEN_HEIGHT = 5000;
    static const double ZOOM_FACTOR;

signals:
    
public slots:
    void sortReads();
    void showOverview();


private:
    list<SeqRead*>      _reads;
    Contig*             _pContig;
    QVector<ReadInfo>   _readInfo;
    QVector<int>        _coverage;
    QVector<int>        _coverageNorm;
    QPixmap*            _pict;
    QPixmap*            _overviewPict;
    QPixmap*            _coveragePlot;
    QPainter            _painter;
    qreal               _prevZoomX;
    qreal               _prevZoomY;
    qreal               _zoomX;
    qreal               _zoomY;
    QAction*            _sortAction;
    QMap<QString, QColor> _colorMap;
    QRectF              _charRect;
    int                 _cRows;
    int                 _cCols;
    int                 _x;
    int                 _y;
    int                 _prevX;
    int                 _prevY;
    int                 _nW;
    int                 _nH;
    int                 _charHeight;
    int                 _charWidth;
    int                 _coveragePlotHeight;
    int                 _maxCoverage;
    bool                _bSortbyHaploType;
    bool                _bPixelMode;
};

#endif // ALIGNMENTPICTURE_H
