#include "mainwindow.h"
#include <QApplication>
#include <QGraphicsEllipseItem>
#include <QRandomGenerator>
#include <qboxlayout.h>
#include <QDebug>
#include <QtMath>
#include <limits.h>
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
    drawContralPoint(scene);
    drawCurate(scene);
    drawGrid(scene);
    setCoeficient();
    setCoeficientG();
    setCoeficientBPrime();
    computeArcLength();
    getIsoTime();
    computeIsoPoint();
    drawSegment(scene);


}

MainWindow::~MainWindow()
{

}

void MainWindow::setSrcdata()
{
    QRandomGenerator random;
    quint32 seed = static_cast<quint32>(QRandomGenerator::global()->bounded(101));
    random.seed(seed);
    for(int i = 0;i<4;++i){
        qreal value_x = static_cast<qreal>(random.bounded(this->width())/2e1);
        qreal value_y = static_cast<qreal>(random.bounded(this->height()/2e1));
        srcdata.append(QPointF(value_x,value_y));
    }
    std::sort(srcdata.begin(),srcdata.end(),[](const QPointF &a,const QPointF &b){
        qreal tolerence = 1e-9;
        return abs(a.x()-b.x())<tolerence?a.y()<b.y():a.x()<b.x();
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
//!     贝塞尔曲线关于时间t的参数方程的四个系数
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

//!     根据计算得到时间t时  曲线上的点
QPointF MainWindow::computePoint(qreal &t)
{
    int size = srcdata.size();
    QPointF resultPoint;
    //!     帕斯卡三角算一次点就要用一次  直接当作成员把吧

    for(int i = 0;i<size;++i){
        resultPoint+=srcdata[i]*pascaTri[i+1]*qPow(t,i)*qPow(1-t,size-i);
    }

    return resultPoint;
}
//!     用牛顿迭代法，得到曲线上的等分点的时间
void MainWindow::getIsoTime()
{
    //!     选择一个初始猜测值x0,它应该尽可能接近方程的根，以提高迭代的收敛速度。
    //!     采用牛顿迭代法计算得到下一个点的  t  坐标
    //!     x_{k+1} = x_{k}-f(x_{k})/f'(x_{k})
    //! 当|x_{k+1}-x_{k}|< δ or |f(t)|< ε 几乎没变化时结束
    //!     f严格单调递增  所以f'(t)严格大于0
    //!     可以采用理查森外推的四阶中心差分公式计算f'(t)
    //!     这里的f(t)是弧长的积分不是贝塞尔曲线的方程
    //!     合理的方式是  在计算desdata时顺便计算弧长公式  这里是以时间  t  均分的
    //!     不行  这样会误差累积         将误差设置比较小  使得累积的误差也比较小
    //!     记F_{i}(τ) = \int_{0}^{τ}  f(t) dt - (i/segcount)*arclength;
    //!     则  F'(τ) = f(τ)
    int MaxIterations = 1e4;
    qreal damping = 5e-4;       //!     由于函数值与导数值的比很大  所以设置阻尼系数
    segtime.append(precis);
    for(int j = 1;j<segcount;++j){
        //!     这里  t=0 被插入了两次
        QList<qreal>::reverse_iterator it = segtime.rbegin();
        qreal x0  = *it;
        for(int i = 0;i<MaxIterations;++i){
            qreal dF = compute_f(x0);
            if(abs(dF)<delta){
                qDebug()<<"The derivative is close to zero, stopping the iteration";
                qDebug()<<"stop to compute the "<<j<<"th segment";
                qDebug()<<"The time of "<<j<<"th segment time is \t"<<x0;
//                segtime.append(x0);
                break;
            }
            //!     一直往负方向跑  什么东西呀
            //!     要么computeSubArcLength(x0) 要么arclength没有  为什么
            qreal subArcLength = computeSubArcLength(x0);
            qDebug("The sub ArcLength from 0 to %llf is %llf",x0,subArcLength);
            qDebug()<<"j:"<<j<<"\t"<<"segcount:"<<segcount<<"\t"<<"arclength:"<<arclength;
            //!     这就很奇怪了   rest为什么不计算
            qreal rest = static_cast<qreal>(j*0.02)*arclength;
            qDebug()<<"The value of rest is"<<rest;
            qreal F = subArcLength-(rest);
            qDebug()<<"F/dF:"<<F/dF;
            qreal x1 = x0-damping*(F/dF);
            //!     采用反射策略
            if(x1<0){
                x1=-x1;
            }else if(x1>1){
                x1=2-x1;
            }
            if(abs(x0-x1)<epslion){
                qDebug()<<"abs(x0-x1)<epslion,stop the iteration";
                qDebug()<<"The time of "<<j<<"th segment time is \t"<<x0;
                //                segtime.append(x1);
                break;
            }

            x0 = x1;
        }
        qDebug()<<"Reaching the MaxIterations,"<<"stop to compute the"<<j<<"th segment";
        qDebug()<<"The time of "<<j<<"th segment time is \t"<<x0;
        segtime.append(x0);
    }
    segtime.append(1);

}

void MainWindow::computeIsoPoint()
{
    for(auto &time:segtime){
        segdata.append(computePoint(time));
    }

}


//!     计算从0——t，贝塞尔曲线的弧长
//!     这里如果是begintime——endtime会不会产生误差累积
//!     难道类似龙贝格积分 再算一次？？  只是选取 M  N  小一些？？？
qreal MainWindow::computeSubArcLength(const qreal &t)
{
    //!     t=0时返回0
    if(abs(t)<1e-9) return 0;
    int n = 5;
    int m = 5;
    QVector<qreal> romberg(n+1,0);
    //!     h = t   表示计算从0——t   贝塞尔曲线的弧长
    qreal h = t;
    //!     为什么这里总会是NAN
    //!     compute_f(*) 我们经常用  不如直接先计算
    //!     但是这里的  h  会不同  所以不能先计算
    int size = static_cast<int>(qPow(2,n)+1e-4);
    QVector<qreal> f(size,0);
    for(int i = 0;i<size;++i){
        f[i] = compute_f(i*h/size);
    }
    romberg[0] = (compute_f(0)+compute_f(h))/2;//!      这里为什么不计算值
    for(int i = 1;i<=n;++i){
        h/=2;
        qreal rest = 0;
        for(qreal j = 1;j<=qPow(2,i-1)+1e-5;++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i] = romberg[i-1]+rest;
    }
    //!     计算m阶龙贝格积分R(n,m)
    for(int i = 1;i<=m;++i){
        qreal temp1 = romberg[i-1];
        qreal temp2 = 0;
        for(int j = i;j<=n;++j){
            temp2 = romberg[j];
            romberg[j] = temp2 + qPow((qPow(4,i)-1),-1)*(temp2-temp1);
            temp1 = temp2;
        }
    }
    return romberg[n];
}


//!     利用数值积分计算得到曲线弧长
void MainWindow::computeArcLength()
{
    int n = 10;
    int m = 10;
//!    qreal interval_num = qPow(2,n);     //!      设置等分区间的数目
    //!      计算零阶龙贝格积分R(n,0)
    QVector<qreal> romberg(n+1,0);
    qreal h = 1;
    //!     还是一样   总是compute_f(0)+compute_f(h)没有
    //!     可能是因为在调用compute_f时  算出来的值不在定义域内
    //!     所以返回NAN
    romberg[0] = (compute_f(0)+compute_f(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        qreal rest = 0;
        for(qreal j = 1;j<=qPow(2,i-1)+1e-5;++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i] = romberg[i-1]+rest;
    }
    //!     计算m阶龙贝格积分R(n,m)
    for(int i = 1;i<=m;++i){
        qreal temp1 = romberg[i-1];
        qreal temp2 = 0;
        for(int j = i;j<=n;++j){
            temp2 = romberg[j];
            romberg[j] = temp2 + qPow((qPow(4,i)-1),-1)*(temp2-temp1);
            temp1 = temp2;
        }
    }
    //!     因为需要计算  2^{M}+1 个函数值，所以通常只选取一个适度的 M 值
    //!     更为精致的算法应该包含一个自动终止程序，当达到指定的误差标准时停止计算
    arclength = romberg[m];

}

//!     设置函数G的系数
void MainWindow::setCoeficientG()
{
    if(!coefficientG.isEmpty()) return;
    coefficientG.append(9*(coefficient.at(0).x()*coefficient.at(0).x()+coefficient.at(0).y()*coefficient.at(0).y()));
    coefficientG.append(12*(coefficient.at(0).x()*coefficient.at(1).x()+coefficient.at(0).y()*coefficient.at(1).y()));
    coefficientG.append(4*(coefficient.at(1).x()*coefficient.at(1).x()+coefficient.at(0).y()*coefficient.at(0).y())+
                        6*(coefficient.at(0).x()*coefficient.at(2).x()+coefficient.at(0).y()*coefficient.at(2).y()));
    coefficientG.append(4*(coefficient.at(1).x()*coefficient.at(2).x()+coefficient.at(1).y()*coefficient.at(2).y()));
    coefficientG.append(coefficient.at(2).x()*coefficient.at(2).x()+coefficient.at(2).y()*coefficient.at(2).y());
}

void MainWindow::setCoeficientBPrime()
{
    coefBPrime.append(3*(-srcdata[0]+3*srcdata[1]-3*srcdata[2]+srcdata[3]));
    coefBPrime.append(6*(srcdata[0]-2*srcdata[1]+srcdata[2]));
    coefBPrime.append(3*(-srcdata[0]+srcdata[1]));
}
qreal MainWindow::compute_f(const qreal &t)
{
    //!     这里加上一个  const   就解决了
    //!     这样   1   0  都可以用了
    //!     还要注意sqrt接受的函数是否可能小于零
    //!     这里result可能小于0  因此我们要将这种
    //!     情况特殊处理
    //!     所以这里返回sqrt(result)是不行的
    //!     只能是返回result  然后根据正负做特殊处理
    if(t<0) return 1e-9;
    double result = 0;
    //!     我直接判断不就好了？？？
    qreal xvalue = coefBPrime[0].x()*t*t+coefBPrime[1].x()*t+coefBPrime[2].x();
    qreal yvalue = coefBPrime[0].y()*t*t+coefBPrime[1].y()*t+coefBPrime[2].y();
    result = qPow(xvalue,2)+qPow(yvalue,2);
//    for(int i = 0;i<5;++i){
//        result+=coefficientG.at(i)*qPow(t,4-i);
//    }
    result = sqrt(result);
    qDebug("f(%llf)=%llf",t,result);
    return result;
}

qreal MainWindow::compute_gradf(const qreal &t)
{
    qreal result = 0;
    for(int i = 0;i<4;++i){
        result+=(4-i)*coefficientG.at(i)*qPow(t,3-i);
    }
    return 0.5*qPow(compute_f(t),-1)*result;
}


//!     画贝塞尔曲线
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
//!     画贝塞尔曲线的控制点
void MainWindow::drawContralPoint(QGraphicsScene *scene)
{
    for(const QPointF &point:srcdata){
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
        item->setBrush(Qt::blue);
        item->setPos(point);
        scene->addItem(item);
    }
}
//!     画坐标栅格图
void MainWindow::drawGrid(QGraphicsScene *scene)
{
    qreal width = this->width();
    qreal height = this->height();
    QPainterPath path;
    QPen pen;
    pen.setWidthF(0.5);

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
    for(const QPointF &point:segdata){
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,5,5);
        item->setBrush(Qt::green);
        item->setPos(point);
        scene->addItem(item);
    }
}
