#ifndef OPENDATADIAG_H
#define OPENDATADIAG_H

#include <QTabWidget>
#include <QDialog>
#include<QPushButton>
#include <QGridLayout>
#include <QLineEdit>
#include <QFileDialog>

#define DIAG_W 800
#define DIAG_H 400

class OpenDataDiag: public QDialog
{
    Q_OBJECT
public:
    OpenDataDiag(QWidget *parent);
    virtual ~OpenDataDiag();

protected:
    QTabWidget *tabWidget;
    QPushButton *ok;
    QPushButton *cancel;
    QVBoxLayout *layout;
    QHBoxLayout *btnLayout;
protected slots:
    void loadData();

};

class VCFTab : public QWidget
{
    Q_OBJECT

public:
    VCFTab();
    virtual ~VCFTab();
protected:
    QLineEdit *filePath;
    QPushButton *browseFile;
    QHBoxLayout *browseF;

protected slots:
    void showOpen();
};

#endif // OPENDATADIAG_H
