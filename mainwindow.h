#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QList>
#include <QGraphicsView>
#include <QGraphicsScene>
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow()override;
    void setSrcdata();
    void getDesdata();
    void setCoeficient();
<<<<<<< HEAD
    void computeArcLength();
    QPointF computePoint(qreal &t);
    QPointF isoPoint(QPointF &beginPoint);
=======
    QPointF computePoint(qreal &t);
    QPointF isoPoint(QPointF &beginPoint);
    qreal computeArcLength();
>>>>>>> 8c39eb3122a5d63f309499200dd7ace477280e65
protected:

private:
    void drawCurate(QGraphicsScene *scene);
    void drawPoint(QGraphicsScene *scene);
    void drawGrid(QGraphicsScene *scene);
    void drawSegment(QGraphicsScene *scene);
private:
//    QVector<int> pascaTri = {1,3,3,1};
    QList<QPointF> srcdata;
    QList<QPointF> desdata;
<<<<<<< HEAD
    QList<QPointF> equidata;        //      存储等距数据点
    QList<QPointF> coefficient;     //      设置贝塞尔曲线关于时间 t 的方差的参数
    qreal arclenght = 0;            //      弧长初始化为零
    qreal precis = 1e-2;            //      设置矩阵法计算的贝塞尔曲线精度
    qreal delta = 1e-3;             //      设置牛顿迭代法参数
    qreal epslion = 1e-5;
    qreal tolerancrerror = 1e-7;    //      设置计算弧长时的容忍误差
    int scale = 50;                 //      设置坐标的刻度
=======
    QList<QPointF> equidata;
    QList<QPointF> coefficient;
    qreal precis = 1e-2;
    qreal delta = 1e-3;
    qreal tolerancrerror = 1e-7;
    int scale = 50;
>>>>>>> 8c39eb3122a5d63f309499200dd7ace477280e65

};

#endif // MAINWINDOW_H
