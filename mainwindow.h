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
    void setSrcdata();                      //      随机设置贝塞尔曲线的控制点
    void getDesdata();                      //      利用矩阵方法计算贝塞尔曲线
    void setCoeficient();                   //      设置贝塞尔曲线关于时间t的参数方程的参数
    void setCoeficientG();                  //      设置函数 g 的参数
    void setCoeficientBPrime();

    void computeArcLength();                //      计算贝塞尔曲线的弧长

    void getIsoTime();                      //      计算贝塞尔曲线等距点的时间值
    void NewtonIterator();                  //      牛顿迭代方法计算等距点的时间值
    void computeIsoPoint();

    qreal compute_f(const qreal &t);        //      计算 f 对应的函数值
    qreal compute_gradf(const qreal &t);    //      计算 f 的导数
    qreal computeSubArcLength(const qreal &t);
    QPointF computePoint(qreal &t);

protected:

private:
    void drawCurate(QGraphicsScene *scene);
    void drawContralPoint(QGraphicsScene *scene);
    void drawGrid(QGraphicsScene *scene);
    void drawSegment(QGraphicsScene *scene);
private:
//    QVector<int> pascaTri = {1,3,3,1};
    QList<QPointF> srcdata;
    QList<QPointF> desdata;

    QList<QPointF> segdata;         //      存储等距数据点
    QList<QPointF> coefficient;     //      存储贝塞尔曲线关于时间 t 的方差的参数
    QList<qreal> coefficientG;

    QList<QPointF> coefBPrime;
    QList<qreal> segtime;
    QVector<int> pascaTri{0,1,3,3,1};
    qreal arclength = 0;            //      弧长初始化为零
    qreal precis = 1e-2;            //      设置矩阵法计算的贝塞尔曲线精度
    qreal delta = 1e-5;             //      设置牛顿迭代法参数
    qreal epslion = 1e-5;
    qreal tolerancrerror = 1e-7;    //      设置计算弧长时的容忍误差
    int scale = 50;                 //      设置坐标的刻度
    int segcount = 50;              //      设置等距分段的数目

};

#endif // MAINWINDOW_H
