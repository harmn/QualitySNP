#ifndef RUNDIALOG_H
#define RUNDIALOG_H

#include <QDialog>
#include <string>
#include "runqsnp.h"

namespace Ui {
class RunDialog;
}

class RunDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit RunDialog(QWidget *parent = 0);
    ~RunDialog();
    

signals:
    void cancelRun();

public slots:
    void runDone();
    void printMessage(QString message);
    void showError(QString);

private slots:
    void on_browseButton_clicked();

    void on_stopButton_clicked();

    void on_runButton_clicked();

    void on_editButton_clicked();

    void on_outputDirButton_clicked();

private:
    void fillComboBox();
    void activateConfiguration();
    Ui::RunDialog *ui;
    RunQSNP* _pRun;
};

#endif // RUNDIALOG_H
