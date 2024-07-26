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

    scene->setBackgroundBrush(Qt::white);

    graphicsView->fitInView(scene->sceneRect(),Qt::KeepAspectRatio);

    // 设置视图的一些属性
    graphicsView->setRenderHint(QPainter::Antialiasing);
    setCentralWidget(graphicsView);


    if(!srcdata.isEmpty()) return;
    setSrcdata();
    if(!desdata.isEmpty()) return;
    getDesdata();
    drawContralPoint(scene);
    drawCurate(scene);
//    drawGrid(scene);
    setCoeficient();
    setCoeficientBPrime();
    setCoeficientG();
    computeArcLength();

    getIsoTime();
//    NewtonIterator();
    computeIsoPoint();
    drawSegment(scene);


}

MainWindow::~MainWindow()
{

}

void MainWindow::setSrcdata()
{
    QRandomGenerator random;
    quint32 seed = static_cast<quint32>(QRandomGenerator::global()->bounded(2001));
    random.seed(seed);
    for(int i = 0;i<4;++i){
        qreal value_x = static_cast<qreal>(random.bounded(this->width())*2);
        qreal value_y = static_cast<qreal>(random.bounded(this->height())*2);
        srcdata.append(QPointF(value_x,value_y));
    }
//    std::sort(srcdata.begin(),srcdata.end(),[](const QPointF &a,const QPointF &b){
//        qreal tolerence = 1e-9;
//        return abs(a.x()-b.x())<tolerence?a.y()<b.y():a.x()<b.x();
//    });

//    srcdata.push_back(QPointF(215.0,76.0));
//    srcdata.push_back(QPointF(302.0,469.0));
//    srcdata.push_back(QPointF(359.0,60.0));
//    srcdata.push_back(QPointF(505.0,207.0));
    for(const QPointF &point:srcdata)
        qDebug()<<point;
}

void MainWindow::getDesdata()
{
    int size = srcdata.size();
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



//!     用二分法，得到曲线上的等分点的时间
void MainWindow::getIsoTime()
{
    //!     选择一个初始猜测值x0,它应该尽可能接近方程的根，以提高迭代的收敛速度。
    //!     采用牛顿迭代法计算得到下一个点的  t  坐标
    //!     x_{k+1} = x_{k}-f(x_{k})/f'(x_{k})
    //!     当|x_{k+1}-x_{k}|< δ or |f(t)|< ε 几乎没变化时结束
    //!     f严格单调递增  所以f'(t)严格大于0
    //!     可以采用理查森外推的四阶中心差分公式计算f'(t)
    //!     这里的f(t)是弧长的积分不是贝塞尔曲线的方程
    //!     合理的方式是  在计算desdata时顺便计算弧长公式  这里是以时间  t  均分的
    //!     不行  这样会误差累积         将误差设置比较小  使得累积的误差也比较小
    //!     记F_{i}(τ) = \int_{0}^{τ}  f(t) dt - (i/segcount)*arclength;
    //!     则  F'(τ) = f(τ)*f'(τ)

    qreal tol = 1e-1;
    segtime.append(0);
    for(qreal i = 1;i<segcount-1e-5;++i){
        qreal left = segtime.at(segtime.size()-1);
        qreal right = 1;
        qreal subArcLength = static_cast<qreal>(i/segcount)*arclength;
        qreal biCutPrircis = 1e-5;
        while(right-left > biCutPrircis){
            qreal mid = (left+right)/2;
            //!     计算中间的函数值和
            //!     采用这种方法的逻辑不对  这样的话只能
            //!     动右边而且要小心提防Right变化过大
            //!     找不到想要的地方
            qreal F_mid = computeSubArcLength(mid);
            qDebug()<<"Compute the "<<i<<"th segment";

            qDebug("The arclength from 0 to %lf is %lf",mid,F_mid);
            if(abs(F_mid-subArcLength)<tol){
                segtime.append(mid);
                break;
            }
            if(F_mid>subArcLength){
                right = mid;
            }else{
                left = mid;
            }
        }

//!     segtime.append((left+right)/2);
//!     默认不会到这一步
    }
    segtime.append(1);

}
//!     牛顿迭代法寻找等分点
void MainWindow::NewtonIterator()
{
    int MaxIterator = 1e4;
    qreal tol = 1e-4;
    qreal damping = 1;
    segtime.append(0);
    for(qreal i = 1;i<segcount;++i){
        qreal subarclength = (i/segcount)*arclength;
        qreal x0 = segtime.at(segtime.size()-1)+1.00/segcount;
        for(int j = 0;j<MaxIterator;++j){
            qreal f0 = computeSubArcLength(x0)-subarclength;
            if(abs(f0)<epslion){
                qDebug("abs(f0): %lf <epslion,stop this loop",f0);
                qDebug("The %lf th segtime is %lf",i,x0);
                segtime.append(x0);
                break;
            }
            qreal df = compute_gradf(x0);
            if(abs(df)<delta){
                qDebug("abs(df): %lf <epslion,stop this loop",df);
                qDebug("The %lf th segtime is %lf",i,x0);
                segtime.append(x0);
                break;
            }
            qreal x1=x0 - damping*(f0/df);
            if(abs(x1-x0)<tol){
                qDebug("x0: %lf, x1: %lf .abs(x1-x0)<tol,stop this loop",x0,x1);
                qDebug("The %lf th segtime is %lf",i,x0);
                segtime.append(x0);
                break;
            }
            if(j == MaxIterator-1){
                qDebug("reaching the MaxIterator , stop this loop");
                qDebug("The %lf th segtime is %lf",i,x0);
                segtime.append(x0);
                break;
            }
            x0 = x1;
        }

    }
    segtime.append(1);
}

void MainWindow::computeIsoPoint()
{
    for(auto &time:segtime){
        segdata.append(computePoint(time));
    }

}
//!     根据计算得到时间t时  曲线上的点
QPointF MainWindow::computePoint(qreal &t)
{

    QPointF resultPoint;

    int size = srcdata.size();
    QVector<qreal> coefficient(size, 0);
    coefficient[0] = 1.000;
    qreal u1 = 1.0 - t;
    //!
    //! De Casteljau递推算法
    //!
    for (int j = 1; j <= size - 1; j++) {
        qreal saved = 0.0;
        for (int k = 0; k < j; k++){
            qreal temp = coefficient[k];
            coefficient[k] = saved + u1 * temp;
            saved = t * temp;
        }
        coefficient[j] = saved;
    }
    for (int i = 0; i < size; i++) {
        QPointF point = srcdata.at(i);
        resultPoint = resultPoint + point * coefficient[i];
    }

    return resultPoint;
}

//!     计算从0——t，贝塞尔曲线的弧长
qreal MainWindow::computeSubArcLength(const qreal &t)
{
    //!     t=0时返回0
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    qreal tol = 1e-1;
    //!     这里容忍误差的设定与坐标的范围相关
    qreal result = 0;
    QVector<QVector<qreal>> romberg(n+1,QVector<qreal>(m+1,0));
    qreal h = t;
    romberg[0][0] = h*(compute_f(0)+compute_f(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        qreal rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=qPow(2,i-1);++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            qreal error = qPow(qPow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    qDebug()<<"\n";
    qDebug("The arclength from 0 to %lf is %lf",t,result);
//    for(int i = 0;i<romberg.size();++i){
//        if(romberg[i][0]<1e-5) break;
//        qDebug()<<romberg[i];
//    }


    return result;
}
//!     计算从begin到end时间曲线的弧长
qreal MainWindow::computeSubArcLength(const qreal &begin, const qreal &end)
{
    //!     t=0时返回0
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    qreal tol = 1e-1;
    //!     这里容忍误差的设定与坐标的范围相关
    qreal result = 0;
    QVector<QVector<qreal>> romberg(n+1,QVector<qreal>(m+1,0));
    qreal t1 = begin;
    qreal t2 = end;
    qreal h = t2-t1;
    romberg[0][0] = h*(compute_f(t1)+compute_f(t2))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        qreal rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=qPow(2,i-1);++j){
            rest+=compute_f(t1+(2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            qreal error = qPow(qPow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    qDebug()<<"\n";
    qDebug("The arclength from %lf to %lf is %lf",t1,t2,result);
    for(int i = 0;i<romberg.size();++i){
        if(romberg[i][0]<1e-5) break;
        qDebug()<<romberg[i];
    }


    return result;
}


//!     利用数值积分计算得到曲线弧长
void MainWindow::computeArcLength()
{
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    qreal tol = 1e-1;
    //!     这里容忍误差的设定与坐标的范围相关
    QVector<QVector<qreal>> romberg(n+1,QVector<qreal>(m+1,0));
    qreal h = 1;
    romberg[0][0] = h*(compute_f(0)+compute_f(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        qreal rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=qPow(2,i-1);++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            arclength = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            qreal error = qPow(qPow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                arclength = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    for(int i = 0;i<romberg.size();++i){
        if(romberg[i][0]<1e-5) break;
        qDebug()<<romberg[i];
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
    coefficient.append(srcdata[0]);
}
void MainWindow::setCoeficientBPrime()
{
    coefBPrime.append(3*coefficient[0]);
    coefBPrime.append(2*coefficient[1]);
    coefBPrime.append(coefficient[2]);
}
//!     设置函数G的系数
//!
//!
//!     注意根号前还有个系数    3
void MainWindow::setCoeficientG()
{
    if(!coefficientG.isEmpty()) return;
    coefficientG.append(coefBPrime.at(0).x()*coefBPrime.at(0).x()+coefBPrime.at(0).y()*coefBPrime.at(0).y());
    coefficientG.append(4*(coefBPrime[0].x()*coefBPrime[1].x()+coefBPrime[0].y()*coefBPrime[1].y()));
    coefficientG.append(2*(coefBPrime[0].x()*coefBPrime[2].x()+coefBPrime[0].y()*coefBPrime[2].y())
                        +4*(coefBPrime[1].x()*coefBPrime[1].x()+coefBPrime[1].y()*coefBPrime[1].y()));
    coefficientG.append(4*(coefBPrime[1].x()*coefBPrime[2].x()+coefBPrime[1].y()*coefBPrime[2].y()));
    coefficientG.append(coefBPrime[2].x()*coefBPrime[2].x()+coefBPrime[2].y()*coefBPrime[2].y());
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
    qreal xvalue = coefBPrime[0].x()*t*t+coefBPrime[1].x()*t+coefBPrime[2].x();
    qreal yvalue = coefBPrime[0].y()*t*t+coefBPrime[1].y()*t+coefBPrime[2].y();
    result = qPow(xvalue,2)+qPow(yvalue,2);

//    qreal xvalue = 3*coefficient[0].x()*t*t+2*coefficient[1].x()*t+coefficient[2].x();
//    qreal yvalue = 3*coefficient[0].y()*t*t+2*coefficient[1].y()*t+coefficient[2].y();
//    result = qPow(xvalue,2)+qPow(yvalue,2);
//    result = 3*sqrt(abs(result));

//!      用下面这个也太离谱了

    result = sqrt(abs(result));
//    qDebug("f(%lf)=%lf",t,result);
    return result;
}

qreal MainWindow::compute_gradf(const qreal &t)
{
    qreal xvalue = coefBPrime.at(0).x()*t*t+2*coefBPrime.at(1).x()*t+coefBPrime.at(2).x();
    qreal xsign = xvalue<0?-1:1;
    xvalue = abs(xvalue)*xsign*2*(coefBPrime.at(0).x()*t+coefBPrime.at(1).x());
    qreal yvalue = coefBPrime.at(0).y()*t*t+2*coefBPrime.at(1).y()*t+coefBPrime.at(2).y();
    qreal ysign = yvalue<0?-1:1;
    yvalue = abs(yvalue)*ysign*2*(coefBPrime.at(0).y()*t+coefBPrime.at(1).y());
    return xvalue+yvalue;
}


//!     画贝塞尔曲线
void MainWindow::drawCurate(QGraphicsScene *scene)
{
    QPainterPath path;
    path.moveTo(desdata.at(0));
    for(const QPointF &point:desdata){
        path.lineTo(point);
    }
    scene->addPath(path,QPen(Qt::red,2,Qt::SolidLine));
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
    for(const QPointF &point:segdata){
        qDebug()<<point;
        QGraphicsEllipseItem *item = new QGraphicsEllipseItem(0,0,3,3);
        item->setBrush(Qt::green);
        item->setPos(point);
        scene->addItem(item);
    }
}
