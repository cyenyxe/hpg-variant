#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    // Set size
    this->setFixedHeight(MAINW_H);
    this->setFixedWidth(MAINW_W);

    this->createMainM();

    // Create main tabs
    mainTab = new QTabWidget();
    mainTab->addTab(new MainLDTab(mainTab), tr("LD Plot"));
    mainTab->addTab(new MainHaploTTab(), tr("Haplotypes"));
    mainTab->addTab(new MainCheckMTab(), tr("Check Markers"));
    mainTab->addTab(new MainTaggerTab(), tr("Tagger"));
    this->setCentralWidget(mainTab);
}

void MainWindow::createMainM()
{
//    showAct = new QAction(("&Show"), this);
//    connect(showAct, SIGNAL(triggered()),label1, SLOT(show()));
    openAct = new QAction(("&Open new data"), this);
    connect(openAct, SIGNAL(triggered()),this, SLOT(openData()));
    exitAct = new QAction(("&Exit"), this);
    connect(exitAct, SIGNAL(triggered()),qApp, SLOT(quit()));
    fileBar = menuBar()->addMenu("&File");
//    mnFile->addAction(showAct);
    fileBar->addAction(openAct);
    fileBar->addSeparator();
    fileBar->addAction(exitAct);
}

void MainWindow::openData()
{
//    QMessageBox msgBox;
//     msgBox.setText("Test open");
//     msgBox.setInformativeText("Test text");
//     msgBox.setStandardButtons(QMessageBox::Ok);
//     msgBox.setDefaultButton(QMessageBox::Ok);
//     int ret = msgBox.exec();
    openDataDiag = new OpenDataDiag(this);
    openDataDiag->open();
}

MainWindow::~MainWindow()
{
   //delete ui;
   delete fileBar;
//   delete showAct;
//   delete hideAct;
   delete exitAct;
}

MainLDTab::MainLDTab(QWidget *parent): QWidget(parent)
{
    // Test data
    marker *test = (marker *) malloc(sizeof(*test) * 3);


    generatePixmap(test, 3);
}

void MainLDTab::paintEvent(QPaintEvent *e)
{
  Q_UNUSED(e);
    QPainter qp(this);
    qp.drawPixmap(0, 0, pixmapLD->width(), pixmapLD->height(), *pixmapLD);
//    QString key = QString("lights:%1:%2")
//                                 .arg(m_color.name())
//                                 .arg(m_diameter);
//           QPixmap pixmap;

//           if (!QPixmapCache::find(key, pixmap)) {
//               pixmap = generatePixmap();
//               QPixmapCache::insert(key, pixmap);
//           }
//           bitBlt(this, 0, 0, &pixmap);
}

void MainLDTab::generatePixmap(marker* markerArr, size_t markerArrLen)
   {

    // Add wait sign and trigger the cleaning
    if (pixmapLD != NULL)
        delete pixmapLD;

    // Processing can be stopped by sending null
    if (markerArr != NULL)
    {
        // Determine the width and height
        unsigned int w = 600, h = 600;

        QPixmap *pixmap = new QPixmap(w, h);
        pixmap->fill(this, 0, 0);
        QPainter painter(pixmap);
        painter.setBrush(Qt::black);

        // Assume that the markers are sorted and the first has the lowest pos in the genome, and with the last having the highest
        uint64_t distance = markerArr[markerArrLen-1].position - markerArr[0].position;

        // Draw relative lines showing the position in the Genome for the markers
        unsigned int x = MAINLDTAB_MARGIN;
        const uint64_t initPosSnp = markerArr[0].position;
        for (size_t idx=0; idx<markerArrLen; idx++)
        {
            painter.drawLine(x, MAINLDTAB_MARGIN, x, MAINLDTAB_MARGIN+10);
            x +=  (markerArr[1].position - initPosSnp) / markerArrLen;
        }




        painter.setBrush(Qt::white);
        const unsigned int cellSnpWidth = painter.fontMetrics().height() + INTER_SNP_SPACE;
        painter.drawRect(MAINLDTAB_MARGIN, MAINLDTAB_MARGIN, ((markerArrLen-1)*cellSnpWidth) , 10);

        //
        unsigned int posTextW = MAINLDTAB_MARGIN +painter.fontMetrics().height()/2 ;
        painter.setPen(Qt::black);
        //painter.begin(pixmap);
        //painter.rotate(-60);



        for (size_t idx=0; idx<markerArrLen; idx++)
        {
            painter.save();
            painter.translate(posTextW, 150);// to correct the text's position
            painter.rotate(-90);
            painter.drawText(0, 0, "SNP 1");
            posTextW += cellSnpWidth;
            // call the basic class's function to draw the text, however rotated.
            painter.restore();
        }
        //painter.end();

//       painter.setBrush(
//               m_color == Qt::red ? Qt::red
//                                  : Qt::lightGray);
//       painter.drawEllipse(0, 0, w, w);

//       painter.setBrush(
//               m_color == Qt::yellow ? Qt::yellow
//                                     : Qt::lightGray);
//       painter.drawEllipse(0, w, w, w);

//       painter.setBrush(
//               m_color == Qt::green ? Qt::green
//                                    : Qt::lightGray);
//       painter.drawEllipse(0, 2 * w, w, w);
        pixmapLD = pixmap;
    } else
           pixmapLD = new QPixmap(1, 1);

   }

MainLDTab::~MainLDTab()
{
    delete pixmapLD;
}

MainHaploTTab::MainHaploTTab() {}
MainHaploTTab::~MainHaploTTab() {}

MainCheckMTab::MainCheckMTab() {}
MainCheckMTab::~MainCheckMTab() {}

MainTaggerTab::MainTaggerTab() {}
MainTaggerTab::~MainTaggerTab() {}
