#ifndef RUNQSNP_H
#define RUNQSNP_H

#include <QObject>

class RunQSNP : public QObject
{
    Q_OBJECT
public:
    explicit RunQSNP();
    void setName(QString name) { _name = name; }

    ~RunQSNP() {
    }
    
signals:
    void done();
    void printMessage(QString);
    void reportError(QString);
    
public slots:
    void run();
    void cancel();

private:
    QString _name;
    bool _bCancelled;
    
};

#endif // RUNQSNP_H
