#include "mainwindow.h"
//#include "segment.h"
#include "archheight.h"
#include <QApplication>
#include <QGraphicsEllipseItem>
#include <QRandomGenerator>
#include <QPainterPath>
#include <qboxlayout.h>
#include <QDebug>
#include <QtMath>
#include <limits.h>
#include <QThread>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QGraphicsView *graphicsView = new QGraphicsView(this);
    QGraphicsScene *scene = new QGraphicsScene(this);
    graphicsView->setScene(scene);

    scene->setBackgroundBrush(Qt::white);

    graphicsView->fitInView(scene->sceneRect(),Qt::KeepAspectRatio);

    graphicsView->setRenderHint(QPainter::Antialiasing);
    setCentralWidget(graphicsView);

    qDebug()<<std::thread::hardware_concurrency();

    initial();
    ArchHeight archHeightSeg(5,2);
    int index = 0;
    for(const vector<Point>& points:archHeightSeg.outSrcPoint()){
        for(const Point& point:points){
            srcdata[index%2].append(point);
        }
        ++index;
    }
    for(const vector<Point>& points:archHeightSeg.outDesPoint()){
        for(const Point& point:points){
            desdata[index%2].append(point);
        }
        ++index;
    }
    for(const vector<Point>& points:archHeightSeg.outPiecePoint()){
        for(const Point& point:points){
            archdata[index%2].append(point);
        }
        ++index;
    }

    for(const vector<Point>& points:archHeightSeg.outAdaptPoint()){
        for(const Point& point:points){
            adaptdata[index%2].append(point);
        }
        ++index;
    }
//    int i = 0;
//    segment seg1(20);
//    for(const vector<Point>& points:seg1.outSrcData()){
//        for(const Point& point:points){
//            srcdata[i%2].append(QPointF(point.x,point.y));
//        }
//        ++i;
//    }
//    for(const vector<Point>& points:seg1.outDesData()){
//        for(const Point& point:points){
//            desdata[i%2].append(QPointF(point.x,point.y));
//        }
//        ++i;
//    }
//    for(const vector<Point>& points:seg1.outSegData()){
//        for(const Point& point:points){
//            segdata[i%2].append(QPointF(point.x,point.y));
//        }
//        ++i;
//    }
//    for(const Point& crosspoint:seg1.outSelfCrossPoint()){
//        selfcrosspoint.append(QPointF(crosspoint.x,crosspoint.y));
//    }
//    for(const Point& crosspoint1:seg1.outCrossPoint()){
//        crosspoint.append(QPointF(crosspoint1.x,crosspoint1.y));
//    }
//    for(const Rectangle& rect:seg1.outMinRectangle()){
//        QPointF topLeft{rect.LeftUp.x,rect.LeftUp.y};
//        QPointF bottomRight{rect.RightDown.x,rect.RightDown.y};
//        minRectangle.append(QRectF{topLeft,bottomRight});
//    }
    drawContralPoint(scene);
    drawCurate(scene);
    drawArchHeight(scene);
    drawAdaptSampling(scene);
//    drawSegment(scene);

//    drawCrossPoint(scene,selfcrosspoint);
//    drawCrossPoint(scene,crosspoint);
//    drawMinRectangle(scene);

}

MainWindow::~MainWindow()
{

}

void MainWindow::initial()
{
    srcdata.resize(2);
    desdata.resize(2);
    segdata.resize(2);
    segtime.resize(2);
    timingseq.resize(2);
    crosstime.resize(2);
    pointseq.resize(2);
    coefficient.resize(2);
    coefBPrime.resize(2);
    coefficientG.resize(2);
    archdata.resize(2);
    adaptdata.resize(2);
}




//!     画贝塞尔曲线
void MainWindow::drawCurate(QGraphicsScene *scene)
{
    if(desdata.empty()) return;
    QPainterPath path;
    int i = 0;
    for(const QVector<QPointF> &points:desdata){
        path.moveTo(desdata[i].at(0));
        for(const QPointF& point:points){
            path.lineTo(point);
        }
        scene->addPath(path,QPen(Qt::red,2,Qt::SolidLine));
        ++i;
    }
}
//!     画贝塞尔曲线的控制点
void MainWindow::drawContralPoint(QGraphicsScene *scene)
{
    for(const QVector<QPointF> &points:srcdata){
        qDebug()<<points<<"\n";
        for(const QPointF& point:points){
            QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
            item->setBrush(Qt::black);
            item->setPos(point);
            scene->addItem(item);
        }
    }
}
//!     画坐标栅格图
void MainWindow::drawGrid(QGraphicsScene *scene)
{
    qreal width = scene->width();
    qreal height = scene->height();
    QPainterPath path;
    QPen pen;
    pen.setWidthF(1.0);

    pen.setStyle(Qt::DotLine);

    for(int t = 0;t<width;t+=scale){
        path.lineTo(QPointF(t,height));
        path.moveTo(QPointF(t+scale,0));
    }
    for(int t = 0;t<height;t+=scale){
        path.lineTo(QPointF(width,t));
        path.moveTo(QPointF(0,t+scale));
    }
    scene->addPath(path,pen);
    qDebug()<<"Grid drawed";
}
//!     画等距分段贝塞尔曲线
void MainWindow::drawSegment(QGraphicsScene *scene)
{
    if(segdata.empty()) return;
    for(const QVector<QPointF> &points:segdata){
//        qDebug()<<points;
        for(const QPointF& point:points){
            QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,3,3);
            item->setBrush(Qt::green);
            item->setPos(point);
            scene->addItem(item);
        }

    }
    qDebug("\n");
}

void MainWindow::drawArchHeight(QGraphicsScene *scene)
{
    if(archdata.empty()) return;
    QPainterPath path;
    int i = 0;
    for(const QVector<QPointF> &points:archdata){
        path.moveTo(archdata[i%2].at(0));
        for(const QPointF& point:points){
            path.lineTo(point);
        }
        scene->addPath(path,QPen(Qt::blue,2,Qt::DashDotLine));
        ++i;

    }
    qDebug("\n");
}


void MainWindow::drawAdaptSampling(QGraphicsScene *scene)
{
    if(adaptdata.empty()) return;
    for(const QVector<QPointF> &points:adaptdata){
//        qDebug()<<points;
        for(const QPointF& point:points){
            QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,3,3);
            item->setBrush(Qt::green);
            item->setPos(point);
            scene->addItem(item);
        }
    }
}


void MainWindow::drawCrossPoint(QGraphicsScene *scene, const QVector<QPointF>& selfcrosspoint)
{
    for(const QPointF& crosspoint:selfcrosspoint){
        qDebug()<<crosspoint;
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
        item->setBrush(Qt::yellow);
        item->setPos(crosspoint);
        scene->addItem(item);
    }
    qDebug("\n");

}

void MainWindow::drawMinRectangle(QGraphicsScene *scene)
{
    for(const QRectF& rect:minRectangle){
//        qDebug()<<rect;
        scene->addRect(rect,QPen(Qt::blue));

    }
    qDebug("MinRectangle drawed \n");

}
