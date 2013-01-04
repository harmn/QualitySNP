#include <list>
#include <QFormLayout>
#include <QDialogButtonBox>
#include <QDir>
#include "settingsdialog.h"
#include "ui_settingsdialog.h"
#include "core/trunk/Configuration.h"

SettingsDialog::SettingsDialog(QWidget *parent) :
    QDialog(parent),
      ui(new Ui::SettingsDialog)
{
     ui->setupUi(this);

    Configuration* pConfig = Configuration::getConfig();

    QString name = (pConfig->getName().empty()) ? tr("Default") : QString::fromStdString(pConfig->getName());

    ui->name->setText(name);

    ui->minimalNumberOfReadsPerAllele->setValue(pConfig->getInt("minimalNumberOfReadsPerAllele"));
    ui->minimalNumberOfReadsPerAllelep->setValue(pConfig->getDouble("minimalNumberOfReadsPerAllelep")*100);
    ui->minimalNumberOfReadsPerHaploType->setValue(pConfig->getInt("minimalNumberOfReadsPerHaploType"));
    ui->lowQualityRegion5prime->setValue(pConfig->getInt("lowQualityRegion5prime"));
    ui->lowQualityRegion3prime->setValue(pConfig->getInt("lowQualityRegion3prime"));
    ui->lowQualityRegion3primePerc->setValue(pConfig->getDouble("lowQualityRegion3primePerc")*100);
    ui->lowComplexityRegionSize->setValue(pConfig->getInt("lowComplexityRegionSize"));
    ui->lowComplexityRepeatCount->setValue(pConfig->getInt("lowComplexityRepeatCount"));
    ui->maxNumberOfSNPsInFlanks->setValue(pConfig->getInt("maxNumberOfSNPsInFlanks"));
    ui->weightHighQualityRegion->setValue(pConfig->getDouble("weightHighQualityRegion"));
    ui->weightLowQualityRegion->setValue(pConfig->getDouble("weightLowQualityRegion"));
    ui->minSNPQualityScore->setValue(pConfig->getInt("minSNPQualityScore"));
    ui->minimalConfidenceScore->setValue(pConfig->getInt("minimalConfidenceScore"));
    ui->minNumberOfHighQualityReads->setValue(pConfig->getInt("minNumberOfHighQualityReads"));
    ui->minimalMappingQuality->setValue(pConfig->getInt("minimalMappingQuality"));
    ui->similarityPerPolymorphicSite->setValue(pConfig->getDouble("similarityPerPolymorphicSite"));
    ui->similarityAllPolymorphicSites->setValue(pConfig->getDouble("similarityAllPolymorphicSites"));
    ui->alleleMajorityThreshold->setValue(pConfig->getDouble("alleleMajorityThreshold"));
    ui->indelLowQuality->setChecked(pConfig->getBool("indelLowQuality"));
    ui->lowComplexityLowQuality->setChecked(pConfig->getBool("lowComplexityLowQuality"));
    ui->useIUPACCodes->setChecked(pConfig->getBool("useIUPACCodes"));
    ui->onlyReliableMarkers->setChecked(pConfig->getBool("onlyReliableMarkers"));
    ui->readNameGroupSeparator->setText(QString::fromStdString(pConfig->getString("readNameGroupSeparator")));
 }

SettingsDialog::~SettingsDialog()
{
}

void SettingsDialog::on_buttonBox_accepted()
{
    Configuration* pConfig = Configuration::getConfig();

    pConfig->setInt("minimalNumberOfReadsPerAllele", ui->minimalNumberOfReadsPerAllele->value());
    pConfig->setDouble("minimalNumberOfReadsPerAllelep",1.0 * ui->minimalNumberOfReadsPerAllelep->value()/100);
    pConfig->setInt("minimalNumberOfReadsPerHaploType", ui->minimalNumberOfReadsPerHaploType->value());
    pConfig->setInt("lowQualityRegion5prime", ui->lowQualityRegion5prime->value());
    pConfig->setInt("lowQualityRegion3prime",ui->lowQualityRegion3prime->value());
    pConfig->setDouble("lowQualityRegion3primePerc",1.0 * ui->lowQualityRegion3primePerc->value()/100);
    pConfig->setInt("lowComplexityRegionSize",ui->lowComplexityRegionSize->value());
    pConfig->setInt("lowComplexityRepeatCount",ui->lowComplexityRepeatCount->value());
    pConfig->setInt("maxNumberOfSNPsInFlanks",ui->maxNumberOfSNPsInFlanks->value());
    pConfig->setDouble("weightHighQualityRegion",ui->weightHighQualityRegion->value());
    pConfig->setDouble("weightLowQualityRegion", ui->weightLowQualityRegion->value());
    pConfig->setInt("minSNPQualityScore",ui->minSNPQualityScore->value());
    pConfig->setInt("minimalConfidenceScore",ui->minimalConfidenceScore->value());
    pConfig->setInt("minNumberOfHighQualityReads",ui->minNumberOfHighQualityReads->value());
    pConfig->setInt("minimalMappingQuality",ui->minimalMappingQuality->value());
    pConfig->setDouble("similarityPerPolymorphicSite",ui->similarityPerPolymorphicSite->value());
    pConfig->setDouble("similarityAllPolymorphicSites",ui->similarityAllPolymorphicSites->value());
    pConfig->setInt("alleleMajorityThreshold", ui->alleleMajorityThreshold->value());
    pConfig->setBool("indelLowQuality",ui->indelLowQuality->isChecked());
    pConfig->setBool("lowComplexityLowQuality",ui->lowComplexityLowQuality->isChecked());
    pConfig->setBool("useIUPACCodes", ui->useIUPACCodes->isChecked());
    pConfig->setBool("onlyReliableMarkers", ui->onlyReliableMarkers->isChecked());
    pConfig->setString("readNameGroupSeparator",ui->readNameGroupSeparator->text().toStdString());

    string configName = ui->name->text().toStdString();
    if(!configName.empty() && configName.compare("Default") != 0) {
        pConfig->setName(configName);
        pConfig->writeFile(configName, "");
    }

    QDialog::accept();
}

