#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QAction>
#include <QTabWidget>
#include <QMenu>
#include <QWidget>
#include <QMenuBar>
#include <QApplication>
#include <QPixmap>
#include <QPainter>
#include <QPoint>
#include <QFont>
#include <QPen>
#include <qpixmapcache.h>
#include <Qt>
#include <QPoint>

#include <inttypes.h>
#include <stdint.h>

extern "C" {
#include "globals.h"
}
#include "opendatadiag.h"

#define MAINW_W 800
#define MAINW_H 600

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private:

    QMenu *fileBar;

    QAction* openAct;
    QAction* hideAct;
    QAction* exitAct;
    OpenDataDiag *openDataDiag;
    QTabWidget *mainTab;

    void createMainM();
public slots:
    void openData();
};


/********************
  Main tabs' classes
  ******************/
/**Margin for the drawing to the window's borders */
#define MAINLDTAB_MARGIN 50

#define INTER_SNP_SPACE 20

#define HPOS_SNP_NAMES 150

//drawing taken from http://doc.qt.digia.com/qq/qq12-qpixmapcache.html
class MainLDTab : public QWidget
{
    Q_OBJECT

public:
    MainLDTab(QWidget *parent);
    virtual ~MainLDTab();
protected:

    QPixmap *pixmapLD;

    void paintEvent(QPaintEvent *e);
    void generatePixmap(marker* markerArr, size_t markerArrLen);
protected slots:
};

class MainHaploTTab : public QWidget
{
    Q_OBJECT

public:
    MainHaploTTab();
    virtual ~MainHaploTTab();
protected:

protected slots:
};

class MainCheckMTab : public QWidget
{
    Q_OBJECT

public:
    MainCheckMTab();
    virtual ~MainCheckMTab();
protected:

protected slots:
};

class MainTaggerTab : public QWidget
{
    Q_OBJECT

public:
    MainTaggerTab();
    virtual ~MainTaggerTab();
protected:

protected slots:
};

#endif // MAINWINDOW_H
