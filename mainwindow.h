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
    QPointF computePoint(qreal &t);
    QPointF isoPoint(QPointF &beginPoint);
    qreal computeArcLength();
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
    QList<QPointF> equidata;
    QList<QPointF> coefficient;
    qreal precis = 1e-2;
    qreal delta = 1e-3;
    qreal tolerancrerror = 1e-7;
    int scale = 50;

};

#endif // MAINWINDOW_H
