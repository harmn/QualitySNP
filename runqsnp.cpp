#include <QMessageBox>
#include <QApplication>
#include <QThreadPool>
#include <QMutex>
#include "runqsnp.h"
#include "core/trunk/Configuration.h"
#include "core/trunk/Contig.h"
#include "core/trunk/ContigProvider.h"
#include "core/trunk/CSVWriter.h"
#include "core/trunk/ContigPrinter.h"

class I : public QThread
{
public:
    static void sleep(unsigned long secs) {
        QThread::sleep(secs);
    }
    static void msleep(unsigned long msecs) {
        QThread::msleep(msecs);
    }
    static void usleep(unsigned long usecs) {
        QThread::usleep(usecs);
    }
};

class QSNPRunTask : public QRunnable
{
    Contig* _pContig;
    CSVWriter* _pWriter;
    bool _bShowContigsWithoutSNP;
    QMutex* _pMutex;

public:
    void setContig(Contig* pContig) {
        _pContig = pContig;
    }

    void setWriter(CSVWriter* pWriter) {
        _pWriter = pWriter;
    }

    void setShowContigsWithoutSNP(bool bShowContigsWithoutSNP) {
        _bShowContigsWithoutSNP = bShowContigsWithoutSNP;
    }

    void setMutex(QMutex* pMutex) {
        _pMutex = pMutex;
    }

    void run() {
        _pContig->sortReads();
        _pContig->calculateProperties();
        if (_bShowContigsWithoutSNP || _pContig->getPotentialSNPCount() > 0) {
            _pMutex->lock();
            _pWriter->writeContig(_pContig);
            _pMutex->unlock();
        }
        delete _pContig;
    }
};

RunQSNP::RunQSNP(): _bCancelled(false)
{
}


void RunQSNP::run()
{
    Configuration* pConfig = Configuration::getConfig();

    Logger* pLogger = Logger::getLogger();
    pLogger->log(QSNP_ALWAYS, "QualitySNP started with settings: " + pConfig->toString());

    ContigProvider contigProvider;
    if (!contigProvider.init()) {
        emit reportError("Could read contig file");
        return;
    }

    vector<string> readGroups;
    if(pConfig->getBool("outputReadGroups")) {
        readGroups = contigProvider.getReadGroups();
    }

    CSVWriter csvWriter(readGroups);
    if (!csvWriter.init()) {
        emit reportError("Could not write results");
        emit done();
        return;
    }

    pConfig->writeFile("config.cfg", pConfig->getString("outputDirectory"));

    QString message;

    if(pConfig->getBool("collectStatistics")) {
        message = QString("Processing a total of %1 contigs with %2 reads").arg(contigProvider.getContigCount()).arg(contigProvider.getReadCount());
        emit printMessage(message);
    }

    QMutex mutex;
    Contig* pContig = contigProvider.nextContig();
    while (pContig != NULL && !_bCancelled) {
        message = "Processing: " + QString::fromStdString(pContig->getName());
        emit printMessage(message);

        QSNPRunTask *runTask = new QSNPRunTask();
        runTask->setContig(pContig);
        runTask->setWriter(&csvWriter);
        runTask->setShowContigsWithoutSNP(pConfig->getBool("showContigsWithoutSNP"));
        runTask->setMutex(&mutex);

        while(!QThreadPool::globalInstance()->tryStart(runTask)) {
            I::msleep(1); // sleep for 1 millisecond before trying again
        }

        message = "Done";
        emit printMessage(message);

        pContig = contigProvider.nextContig();
    }

    QThreadPool::globalInstance()->waitForDone();

    emit done();
}

void RunQSNP::cancel()
{
    _bCancelled = true;
}
