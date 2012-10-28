#include "opendatadiag.h"



OpenDataDiag::OpenDataDiag(QWidget *parent): QDialog(parent)
{
    this->setFixedHeight(DIAG_H);
    this->setFixedWidth(DIAG_W);
    layout = new QVBoxLayout(this);
    btnLayout = new QHBoxLayout();
    ok = new QPushButton("OK");
    connect(ok, SIGNAL(clicked()),this, SLOT(loadData()));
    cancel = new QPushButton("Cancel");
    connect(cancel, SIGNAL(clicked()),this, SLOT(close()));
    this->setWindowTitle("This is a diag");
    tabWidget = new QTabWidget();
    tabWidget->addTab(new VCFTab(), tr("VCF"));
    tabWidget->addTab(new VCFTab(), tr("PED"));
    tabWidget->setTabPosition(QTabWidget::West);
    layout->addWidget(tabWidget);
    btnLayout->addWidget(ok);
    btnLayout->addWidget(cancel);
    layout->addLayout(btnLayout);
}

OpenDataDiag::~OpenDataDiag()
{
    delete ok;
    delete cancel;
    delete tabWidget;
    delete layout;
}

void OpenDataDiag::loadData()
{

}

VCFTab::VCFTab()
{
    browseF = new QHBoxLayout(this);
    filePath = new QLineEdit();
    browseFile = new QPushButton("Browse...");
    connect(browseFile, SIGNAL(clicked()),this, SLOT(showOpen()));
    browseF->addWidget(filePath);
    browseF->addWidget(browseFile);
}

VCFTab::~VCFTab()
{
    delete filePath;
    delete browseFile;
    delete browseF;
}

void VCFTab::showOpen()
{
    const QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "", tr("All files (*.*);;VCF (*.vcf"));
    if (fileName != NULL)
        filePath->setText(fileName);
}

