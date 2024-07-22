#include "mainwindow.h"
#include <QApplication>
#include <QGraphicsEllipseItem>
#include <QRandomGenerator>
#include <qboxlayout.h>
#include <QDebug>
#include <QtMath>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QGraphicsView *graphicsView = new QGraphicsView(this);
    QGraphicsScene *scene = new QGraphicsScene(this);
    graphicsView->setScene(scene);

    // 设置视图的一些属性
    graphicsView->setRenderHint(QPainter::Antialiasing);
//    graphicsView->setFixedSize(400, 400); // 设置视图固定大小
//    graphicsView->setSizeAdjustPolicy(QSizePolicy::Expanding);
    setCentralWidget(graphicsView);


    if(!srcdata.isEmpty()) return;
    setSrcdata();
    if(!desdata.isEmpty()) return;
    getDesdata();
    drawPoint(scene);
    drawCurate(scene);
    drawGrid(scene);

}

MainWindow::~MainWindow()
{

}

void MainWindow::setSrcdata()
{
    QRandomGenerator random;
    quint32 seed = static_cast<quint32>(QRandomGenerator::global()->bounded(501));
    random.seed(seed);
    for(int i = 0;i<4;++i){
        qreal value_x = static_cast<qreal>(random.bounded(this->width()));
        qreal value_y = static_cast<qreal>(random.bounded(this->height()));
        srcdata.append(QPointF(value_x,value_y));
    }
    std::sort(srcdata.begin(),srcdata.end(),[](const QPointF &a,const QPointF &b){
        return a.x()==b.x()?a.y()<b.y():a.x()<b.x();
    });
    for(const QPointF &point:srcdata)
        qDebug()<<point;
}

void MainWindow::getDesdata()
{
    int size = srcdata.size();
    QVector<int> pascaTri(size+1,0);
    pascaTri[1] = 1;
    for(int i = 1;i<=size;++i){
        int temp1 = pascaTri[0];
        int temp2 = 0;
        for(int j = 1;j<=i;++j){
            temp2 = pascaTri[j];
            pascaTri[j]=temp1+temp2;
            temp1 = temp2;
        }
    }

    for(qreal t = 0;t<1.0000;t+=precis){
        QPointF resultPoint={0,0};
        QVector<qreal> bernstein1(size,1),bernstein2(size,1);

        for(int i = 1;i<=size-1;++i){
            bernstein1[i]*=bernstein1[i-1]*t;
            bernstein2[size-i-1]*=bernstein2[size-i]*(1-t);
        }
        for(int i = 0;i<size;++i){
            resultPoint+=srcdata[i]*pascaTri[i+1]*bernstein1[i]*bernstein2[i];
        }
        desdata.append(resultPoint);
    }
}

void MainWindow::setCoeficient()
{
    if(!coefficient.isEmpty()) return;
    //!     coefficient[0]-[3]
    //!     贝塞尔曲线的四个系数   t in [0,1]
    coefficient.append(-srcdata[0]+3*srcdata[1]-3*srcdata[2]+srcdata[3]);
    coefficient.append(3*srcdata[0]-6*srcdata[1]+3*srcdata[2]);
    coefficient.append(-3*srcdata[0]+3*srcdata[1]);
    coefficient.append(srcdata[0]-srcdata[3]);
}

QPointF MainWindow::computePoint(qreal &t)
{
    int size = srcdata.size();
    QPointF resultPoint;
    QVector<int> pascaTri(size+1,0);
    pascaTri[1] = 1;
    for(int i = 1;i<=size;++i){
        int temp1 = pascaTri[0];
        int temp2 = 0;
        for(int j = 1;j<=i;++j){
            temp2 = pascaTri[j];
            pascaTri[j]=temp1+temp2;
            temp1 = temp2;
        }
    }
    for(int i = 0;i<size;++i){
        resultPoint+=srcdata[i]*pascaTri[i+1]*qPow(t,i)*qPow(1-t,size-i);
    }

    return resultPoint;
}

QPointF MainWindow::isoPoint(QPointF &beginPoint)
{
    QPointF endPoint;
    //!     采用牛顿迭代法计算得到下一个点的  t  坐标
    //!     x_{k+1} = x_{k}-f(x_{k})/f'(x_{k})
    //! 当|x_{k+1}-x_{k}|<epslion or |f(t)|足够小 几乎没变化时结束
    //!     f严格单调递增  所以f'(t)严格大于0
    //!     可以采用理查森外推的四阶中心差分公式计算f'(t)
    //!     这里的f(t)是弧长的积分不是贝塞尔曲线的方程
    //!     合理的方式是  在计算desdata时顺便计算弧长公式  这里是以时间  t  均分的
    //!     不行  这样会误差累积         将误差设置比较小  使得累积的误差也比较小
    return endPoint;
}

qreal MainWindow::computeArcLength()
{
    qreal arclength;
    setCoeficient();
    QVector<qreal> g_coeficient(5,0);
    g_coeficient[0] = 9*(coefficient.at(0).x()*coefficient.at(0).x()+coefficient.at(0).y()*coefficient.at(0).y());
    g_coeficient[1] = 12*(coefficient.at(0).x()*coefficient.at(1).x()+coefficient.at(0).y()*coefficient.at(1).y());

    g_coeficient[2] = 4*(coefficient.at(1).x()*coefficient.at(1).x()+coefficient.at(0).y()*coefficient.at(0).y());
    g_coeficient[2]+= 6*(coefficient.at(0).x()*coefficient.at(2).x()+coefficient.at(0).y()*coefficient.at(2).y());

    g_coeficient[3] = 4*(coefficient.at(1).x()*coefficient.at(2).x()+coefficient.at(1).y()*coefficient.at(2).y());
    g_coeficient[4] = coefficient.at(2).x()*coefficient.at(2).x()+coefficient.at(2).y()*coefficient.at(2).y();
    //!     利用龙贝格积分计算弧长

    return arclength;
}



void MainWindow::drawCurate(QGraphicsScene *scene)
{
    QPainterPath path;
    path.moveTo(desdata.at(0));
    for(const QPointF &point:desdata){
        //!     不画点  只画线
        path.lineTo(point);
    }
    scene->addPath(path,QPen(Qt::red));
}

void MainWindow::drawPoint(QGraphicsScene *scene)
{
    for(const QPointF &point:srcdata){
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
        item->setBrush(Qt::blue);
        item->setPos(point);
        scene->addItem(item);
    }
}

void MainWindow::drawGrid(QGraphicsScene *scene)
{
    qreal width = this->width();
    qreal height = this->height();
    QPainterPath path;
    QPen pen;
    pen.setWidthF(0.5);
    pen.setStyle(Qt::DashDotLine);
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

void MainWindow::drawSegment(QGraphicsScene *scene)
{
    for(const QPointF &point:equidata){
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
        item->setBrush(Qt::green);
        item->setPos(point);
        scene->addItem(item);
    }
}
