#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QApplication>

#include <QList>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QEvent>
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow()override;
    void initial();
protected:

private:
    void drawCurate(QGraphicsScene *scene);
    void drawContralPoint(QGraphicsScene *scene);
    void drawGrid(QGraphicsScene *scene);
    void drawSegment(QGraphicsScene *scene);
    void drawSelfCrossPoint(QGraphicsScene *scene,const QVector<QPointF>& selfcrosspoint);

private:
//    QVector<int> pascaTri = {1,3,3,1};
    QVector<QVector<QPointF>> srcdata;
    QVector<QVector<QPointF>> desdata;

    QVector<QVector<QPointF>> segdata;         //      存储等距数据点
    QVector<QVector<QPointF>> coefficient;     //      存储贝塞尔曲线关于时间 t 的方差的参数
    QVector<QVector<QPointF>> coefBPrime;
    QVector<QVector<qreal>> coefficientG;

    QVector<QPointF> selfcrosspoint;
    QVector<QVector<qreal>> segtime;
    QVector<QVector<QPointF>> pointseq;
    QVector<QVector<qreal>> timingseq;
    QVector<qreal> crosstime;
    QVector<int> pascaTri{0,1,3,3,1};
    qreal arclength = 0;            //      弧长初始化为零
    qreal precis = 1e-2;            //      设置矩阵法计算的贝塞尔曲线精度
    qreal delta = 1e-5;             //      设置牛顿迭代法参数
    qreal epslion = 1e-5;
    qreal tolerancrerror = 1e-7;    //      设置计算弧长时的容忍误差
    int scale = 50;                 //      设置坐标的刻度
    int segcount = 20;              //      设置等距分段的数目

};

#endif // MAINWINDOW_H
