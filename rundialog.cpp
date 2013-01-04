#include <QFileDialog>
#include <QMessageBox>
#include <QThread>
#include <QDir>
#include "rundialog.h"
#include "ui_rundialog.h"
#include "settingsdialog.h"
#include "core/trunk/Configuration.h"
#include "core/trunk/Contig.h"
#include "core/trunk/ContigProvider.h"
#include "core/trunk/CSVWriter.h"
#include "core/trunk/ContigPrinter.h"

RunDialog::RunDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RunDialog), _pRun(NULL)
{
    ui->setupUi(this);

    Configuration* pConfig = Configuration::getConfig();
    string contigFileName = pConfig->getString("contigFileName");
    if(!contigFileName.empty()) {
        ui->contigFileName->setText(QString::fromStdString(contigFileName));
        ui->runButton->setEnabled(true);
    }

    ui->runName->setText(QString().fromStdString(pConfig->getString("outputDirectory")));
    fillComboBox();
}

RunDialog::~RunDialog()
{
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
    delete ui;
}

void RunDialog::fillComboBox() {
    ui->configurationComboBox->clear();
    Configuration* pConfig = Configuration::getConfig();
    string currentConfigName = pConfig->getName();

    map<string, string> configs = pConfig->findConfigurations();
    int currentItem = 0;
    for(map<string, string>::const_iterator itConfs = configs.begin(); itConfs != configs.end(); itConfs++) {
        string configName = (*itConfs).first;
        ui->configurationComboBox->addItem(QString::fromStdString(configName));

        if(configName == currentConfigName) {
            currentItem = ui->configurationComboBox->count() - 1;
        }
    }
    ui->configurationComboBox->setCurrentIndex(currentItem);
}

void RunDialog::on_browseButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Choose a contig file");
    ui->contigFileName->setText(fileName);
    ui->runButton->setEnabled(true);
}

void RunDialog::on_stopButton_clicked()
{
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
    emit cancelRun();
    this->close();
}

class QThreadEx : public QThread
{
protected:
    void run() { exec(); }
};

void RunDialog::activateConfiguration() {
    Configuration* pConfig = Configuration::getConfig();
    string selectedConfig = ui->configurationComboBox->currentText().toStdString();
    if(pConfig->getName() != selectedConfig) {
        pConfig = Configuration::getConfig(true);
        if(selectedConfig != "default") {
            map<string, string> configs = pConfig->findConfigurations();
            pConfig->readFile(configs[selectedConfig], "");
        }
    }
}

void RunDialog::on_runButton_clicked()
{
    activateConfiguration();

    Configuration* pConfig = Configuration::getConfig();

    string outputDirectory = ui->runName->text().toStdString();

    pConfig->setString("contigFileName", ui->contigFileName->text().toStdString());

    pConfig->setString("outputDirectory", outputDirectory);
    pConfig->setString("inputDirectory", outputDirectory);

    QThreadEx* pThread = new QThreadEx();
    _pRun = new RunQSNP();
    _pRun->moveToThread(pThread);
    connect(this, SIGNAL(cancelRun()), _pRun, SLOT(cancel()), Qt::DirectConnection);
    connect(pThread, SIGNAL(started()), _pRun, SLOT(run()));
    connect(_pRun, SIGNAL(done()), pThread, SLOT(quit()));
    connect(_pRun, SIGNAL(done()), this, SLOT(runDone()));
    connect(_pRun, SIGNAL(printMessage(QString)), this, SLOT(printMessage(QString)));

    pThread->start();
    printMessage(tr("Started"));

#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif

    ui->runButton->setEnabled(false);
    ui->stopButton->setEnabled(true);
    ui->browseButton->setEnabled(false);
    ui->editButton->setEnabled(false);
    ui->outputDirButton->setEnabled(false);
    ui->configurationComboBox->setEnabled(false);
}

void RunDialog::runDone()
{
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif
    ui->stopButton->setText("Show");
}

void RunDialog::printMessage(QString message)
{
    ui->logText->appendPlainText(message);
}

void RunDialog::showError(QString message)
{
    QMessageBox messageBox;
    messageBox.setIcon(QMessageBox::Critical);
    messageBox.setText(message);
    messageBox.exec();
    ui->runButton->setEnabled(true);
    ui->logText->clear();
}

void RunDialog::on_editButton_clicked()
{
    activateConfiguration();
    Configuration* pConfig = Configuration::getConfig();
    pConfig->setString("contigFileName", ui->contigFileName->text().toStdString());

    SettingsDialog settingsDialog(this);
    settingsDialog.exec();

    if(settingsDialog.result() == QDialog::Accepted) {
        fillComboBox();
        int index = ui->configurationComboBox->findText(QString::fromStdString(pConfig->getName()));
        ui->configurationComboBox->setCurrentIndex(index);
    }
}

void RunDialog::on_outputDirButton_clicked()
{
    Configuration* pConfig = Configuration::getConfig();
    QString dirName = QFileDialog::getExistingDirectory(this, tr("Choose the output directory"), QString().fromStdString(pConfig->getString("outputDirectory")));
    ui->runName->setText(dirName);
}
