#include "segment.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <algorithm>
using namespace std;
segment::segment(const int& count)
{
   segcount = count;
   if(!srcdata.empty()) return;
   setSrcData();
   if(!desdata.empty()) return;
   getDesData();
   setCoeficient();
   setCoeficientBPrime();
   setCoeficientG();
//    computeArcLength();

//    getIsoTime();
//    NewtonIterator();
//    computeIsoPoint();
   computeCrossTime();

}
segment::segment()
{
   if(!srcdata.empty()) return;
   setSrcData();
   if(!desdata.empty()) return;
   getDesData();
   setCoeficient();
   setCoeficientBPrime();
   setCoeficientG();
//    computeArcLength();

//    getIsoTime();
//    NewtonIterator();
//    computeIsoPoint();
   computeCrossTime();

}


segment::~segment()
{
}

void segment::setSrcData()
{
    // srand(time(0));
    // for (int i = 0; i < 4; ++i)
    // {
    //     srcdata.emplace_back(rand()/4e1,rand()/4e1);
    // }
   srcdata.emplace_back(Point(500.0,0.0));
   srcdata.emplace_back(Point(0.0,500.0));
   srcdata.emplace_back(Point(500,500.0));
   srcdata.emplace_back(Point(60.0,22.0));
}

void segment::getDesData()
{
    size_t size = srcdata.size();
    for (double t = 0; t < 1.0000; t += precis)
    {
        Point resultPoint = {0, 0};
        vector<double> bernstein1(size, 1), bernstein2(size, 1);
        for (int i = 1; i <= size - 1; ++i)
        {
            bernstein1[i] *= bernstein1[i - 1] * t;
            bernstein2[size - i - 1] *= bernstein2[size - i] * (1 - t);
        }
        for (int i = 0; i < size; ++i)
        {
            resultPoint += srcdata[i] * pascaTri[i + 1] * bernstein1[i] * bernstein2[i];
        }
        desdata.emplace_back(resultPoint);
    }
}

void segment::setCoeficient()
{
    if (!coefficient.empty()) return;
    //     coefficient[0]-[3]
    //     贝塞尔曲线的四个系数   t in [0,1]
    coefficient.emplace_back(-srcdata[0] + 3 * srcdata[1] - 3 * srcdata[2] + srcdata[3]);
    coefficient.emplace_back(3 * srcdata[0] - 6 * srcdata[1] + 3 * srcdata[2]);
    coefficient.emplace_back(-3 * srcdata[0] + 3 * srcdata[1]);
    coefficient.emplace_back(1*srcdata[0]);
}

void segment::setCoeficientG()
{
    if(!coefficientG.empty()) return;
    coefficientG.emplace_back(coefBPrime[0].x*coefBPrime[0].x+coefBPrime[0].y*coefBPrime[0].y);
    coefficientG.emplace_back(4*(coefBPrime[0].x*coefBPrime[1].x+coefBPrime[0].y*coefBPrime[1].y));
    coefficientG.emplace_back(2*(coefBPrime[0].x*coefBPrime[2].x+coefBPrime[0].y*coefBPrime[2].y)
                        +4*(coefBPrime[1].x*coefBPrime[1].x+coefBPrime[1].y*coefBPrime[1].y));
    coefficientG.emplace_back(4*(coefBPrime[1].x*coefBPrime[2].x+coefBPrime[1].y*coefBPrime[2].y));
    coefficientG.emplace_back(coefBPrime[2].x*coefBPrime[2].x+coefBPrime[2].y*coefBPrime[2].y);
}

void segment::setCoeficientBPrime()
{
    coefBPrime.emplace_back(3 * coefficient[0]);
    coefBPrime.emplace_back(2 * coefficient[1]);
    coefBPrime.emplace_back(1 * coefficient[2]);
}

void segment::computeArcLength()
{
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;
    //!     这里容忍误差的设定与坐标的范围相关
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = 1;
    romberg[0][0] = h*(compute_f(0)+compute_f(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            arclength = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
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
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),ostream_iterator<double>(std::cout," "));
    // }

}

void segment::getIsoTime()
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

    double tol = 1e-1;
    segtime.emplace_back(0);
    for(int i = 1;i<segcount-1e-5;++i){
        double left = segtime.at(segtime.size()-1);
        double right = 1;
        double subArcLength = static_cast<double>(i *1.0/segcount)*arclength;
        double biCutPrircis = 1e-5;
        while(right-left > biCutPrircis){
            double mid = (left+right)/2;
            //!     计算中间的函数值和
            //!     采用这种方法的逻辑不对  这样的话只能
            //!     动右边而且要小心提防Right变化过大
            //!     找不到想要的地方
            double F_mid = computeSubArcLength(mid);
            std::cout<<std::endl;
            std::cout<<"Compute the "<<i<<"th segment"<<std::endl;
            std::cout<<"The arclength from 0 to "<<std::setprecision(6)<<mid<<" is "
                    <<std::setprecision(10)<<F_mid<<std::endl;
            if(abs(F_mid-subArcLength)<tol){
                segtime.emplace_back(mid);
                break;
            }
            if(F_mid>subArcLength){
                right = mid;
            }else{
                left = mid;
            }
        }

//!     segtime.emplace_back((left+right)/2);
//!     默认不会到这一步
    }
    segtime.emplace_back(1);
}

void segment::NewtonIterator()
{
    int MaxIterator = 1e2;
    double tol = 1e-4;
    double damping = 1;
    segtime.emplace_back(0);
    for(double i = 1;i<segcount;++i){
        double subarclength = (i/segcount)*arclength;
        double x0 = segtime.at(segtime.size()-1)+1.00/segcount;
        for(int j = 0;j<MaxIterator;++j){
            double f0 = computeSubArcLength(x0)-subarclength;
            if(abs(f0)<epslion){
                printf("abs(f0): %lf <epslion,stop this loop",f0);
                printf("The %lf th segtime is %lf",i,x0);
                segtime.emplace_back(x0);
                break;
            }
            double df = compute_gradf(x0);
            if(abs(df)<delta){
                printf("abs(df): %lf <epslion,stop this loop",df);
                printf("The %lf th segtime is %lf",i,x0);
                segtime.emplace_back(x0);
                break;
            }
            double x1=x0 - damping*(f0/df);
            if(abs(x1-x0)<tol){
                printf("x0: %lf, x1: %lf .abs(x1-x0)<tol,stop this loop",x0,x1);
                printf("The %lf th segtime is %lf",i,x0);
                segtime.emplace_back(x0);
                break;
            }
            if(j == MaxIterator-1){
                printf("reaching the MaxIterator , stop this loop");
                printf("The %lf th segtime is %lf",i,x0);
                segtime.emplace_back(x0);
                break;
            }
            x0 = x1;
        }

    }
    segtime.emplace_back(1);
}

void segment::computeIsoPoint()
{
    for(auto &time:segtime){
    segdata.emplace_back(computePoint(time));
    }
}
//      计算贝塞尔曲线自相交的交点情况
//      如果A、B线性无关  交点很好算  能够根据公式直接求
//      如果A、B线性相关  这样要怎么说明呢
//      判断A、C的线性相关性  如果A、C线性无关  那么
//      C无法由A、B线性表示  从而曲线没有自相交的点
//      如果A、C也线性相关   那么A、B、C三个向量是等价的
//      即A、B、C三个向量在同一条直线上(A、B、C三个向量不能有零向量)
//      这时贝塞尔曲线的四个控制点在同一条直线上  曲线没有自相交的点
//      如果B或C为零向量  表明四个控制点有三个在同一直线上(重合的点也算)
//      这时曲线没有自相交的点
//      如果A为零向量     
//      说明三阶贝塞尔曲线变成了二阶贝塞尔曲线  只不过这时
//      是由四个控制点控制而已   相当于是升阶(或者说是降阶)了    这时也没有自相交的点
//      综上所述   只有当A、B线性无关的时候曲线才可能会有相交的点  只要在这个时候求交点就好了
void segment::computeCrossTime()
{
    if((coefficient.at(0) == Point{0,0})|| (coefficient.at(1) == Point{0,0})
                                        ||(coefficient.at(2) == Point{0,0})) return;
    double judgement1 = coefficient.at(0).x*coefficient.at(1).y-coefficient.at(0).y*coefficient.at(1).x;
    if(abs(judgement1)<1e-10) return;
    double k1 = (coefficient.at(2).x*coefficient.at(1).y-coefficient.at(2).y*coefficient.at(1).x)/judgement1;
    double k2 = (coefficient.at(0).x*coefficient.at(2).y-coefficient.at(0).y*coefficient.at(2).x)/judgement1;
    if(k1>0 || k2>0) return;
    double delta = -4*k1-3*k2*k2;
    //!     当delta判别式等于0的时候  time1==time2  我们把这种情况剔除
    if(delta<=0) return;
    double time1 = (-k2+sqrt(delta))/2;
    double time2 = (-k2-sqrt(delta))/2;
    //!     如果有一个不在范围内的话   另外就不能与之称为一对有效的解
    //!     从而只要有一个不对  那么另外一个也不需要了
    if(time1>=0 && time1 <=1 && time2>=0 && time2<=1){
        selfcrosstime = time1;
    }
    if(selfcrosstime>=0){
        printf("The bezier curve self-cross at time %llf",selfcrosstime);
    }else{
        printf("The bezier curve is not self-cross");
    }
}

void segment::computeCrossTime1()
{
    //! 计算两个贝塞尔曲线交点  还是只能采用包络框的方式 由computTimeSeq函数计算
    //! 得到的两条曲线单调性变化的位置最后由二分法计算   可以得到所有的交点
    //! 所以我们要做的就是能够快速找到有交的包络框
    //! 由于timeseq1、timeseq2的定义  我们可以知道在两个时间段内曲线要么是严格凹的
    //! 要么是严格凸的  

    //      先做退出去的判断   这里不能对crosstime做判断  因为我们是要做递归调用 
    //      采用类似于剪枝的方法去做的
    //      其实可以暴力的使用这种方法去计算交点  比如将0-1以0.001分为1000份
    //      进行一百万次计算就可以精确到0.001  很是不错
    //      如果通过先暴力  平均分为10份   总共分成100个区域  结合曲线的单调性来
    //      计算  将在区间内严格单调的
    //      不是  我直接根据计算而来的两个贝塞尔曲线单调性变化的数据timeseq将
    //      两个贝塞尔曲线直接分成这些曲线段  然后用牛顿迭代法去找不就是了  这样岂不是更快
    //      这样最多要执行25次牛顿迭代  但是牛顿迭代是二次收敛的   怕什么
    //      
    
    //      double eps = 1e-2;
    for(double i = 0;i<1.0;i+=0.01){
        double end1 = i+0.01;
        Point point1 = computePoint(i);
        Point point2 = computePoint(end1);
        for(double j = 0;j<1.0;j+=0.01){
            double end2 = j+0.01;
            Point point3 = computePoint1(j);
            Point point4 = computePoint1(end2);
            if(isRectangleIntersecting(point1,point2,point3,point4)){
                crosstime.emplace_back(i);
                crosstime.emplace_back(end1);
                crosstime.emplace_back(j);
                crosstime.emplace_back(end2);
            }
        }
    }
}

void segment::computeTimeSeq()
{   //  计算 x 参数方程单调性变化的两个时间点  我们只要单调性变化的时间点  
    //  恒为常数或者一直递增或递减可以不存入  这些情况都没有单调性的变化
    //  但是这样一直用if 确实不太好
    if(true){
        if(coefBPrime.at(0).x != 0){
            double delta = coefBPrime.at(1).x*coefBPrime.at(1).x-4*coefBPrime.at(0).x*coefBPrime.at(2).x;
            if(delta>=0){
                double time1 = (-coefBPrime.at(1).x+sqrt(delta))/(2*coefBPrime.at(0).x);
                double time2 = (-coefBPrime.at(1).x-sqrt(delta))/(2*coefBPrime.at(0).x);
                if(time1>=0 && time1 <=1){
                    timeseq.emplace_back(time1);
                }
                if(time2>=0 && time2 <=1){
                    timeseq.emplace_back(time2);
                }
            }
        }else if(coefBPrime.at(1).x != 0){
            double time1 = -coefBPrime.at(2).x/coefBPrime.at(1).x;
            if(time1>=0 && time1 <=1){
                timeseq.emplace_back(time1);
            }
        }
    }
    //  计算 y 参数方程单调性变化的两个时间点
    if(true){
        if(coefBPrime.at(0).y != 0){
            double delta = coefBPrime.at(1).y*coefBPrime.at(1).y-4*coefBPrime.at(0).y*coefBPrime.at(2).y;
            if(delta>=0){
                double time1 = (-coefBPrime.at(1).y+sqrt(delta))/(2*coefBPrime.at(0).y);
                double time2 = (-coefBPrime.at(1).y-sqrt(delta))/(2*coefBPrime.at(0).y);
                if(time1>=0 && time1 <=1){
                    timeseq.emplace_back(time1);
                }
                if(time2>=0 && time2 <=1){
                    timeseq.emplace_back(time2);
                }
            }
        }else if(coefBPrime.at(1).y != 0){
            double time1 = -coefBPrime.at(2).y/coefBPrime.at(1).y;
            if(time1>=0 && time1 <=1){
                timeseq.emplace_back(time1);
            }
        }
    }
    //      对单调性变化的时间进行逆排序   从1到0的逆序排序
    std::sort(timeseq.begin(),timeseq.end());
    reverse(timeseq.begin(),timeseq.end());
}

bool segment::isRectangleIntersecting(const Point &p1, const Point &p2, const Point &p3, const Point &p4)
{
    return !(p2.x < p3.x || p1.x > p4.x || p2.y > p3.y || p1.y < p4.y);
}

double segment::compute_f(const double &t)
{
        if(t<0) return 1e-9;
    double result = 0;
    double xvalue = coefBPrime[0].x*t*t+coefBPrime[1].x*t+coefBPrime[2].x;
    double yvalue = coefBPrime[0].y*t*t+coefBPrime[1].y*t+coefBPrime[2].y;
    result = pow(xvalue,2)+pow(yvalue,2);


    result = sqrt(abs(result));
//    printf("f(%lf)=%lf",t,result);
    return result;
}

double segment::compute_gradf(const double &t)
{
    double xvalue = coefBPrime.at(0).x*t*t+2*coefBPrime.at(1).x*t+coefBPrime.at(2).x;
    double xsign = xvalue<0?-1:1;
    xvalue = abs(xvalue)*xsign*2*(coefBPrime.at(0).x*t+coefBPrime.at(1).x);
    double yvalue = coefBPrime.at(0).y*t*t+2*coefBPrime.at(1).y*t+coefBPrime.at(2).y;
    double ysign = yvalue<0?-1:1;
    yvalue = abs(yvalue)*ysign*2*(coefBPrime.at(0).y*t+coefBPrime.at(1).y);
    return xvalue+yvalue;
}

double segment::computeSubArcLength(const double &t)
{
    //!     t=0时返回0
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;
    //!     这里容忍误差的设定与坐标的范围相关
    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = t;
    romberg[0][0] = h*(compute_f(0)+compute_f(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
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
    std::cout<<"\n";
    printf("The arclength from 0 to %lf is %lf",t,result);
//    for(int i = 0;i<romberg.size();++i){
//        if(romberg[i][0]<1e-5) break;
//        std::cout<<romberg[i];
//    }


    return result;
}

double segment::computeSubArcLength(const double &begin, const double &end)
{
        //!     t=0时返回0
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;
    //     这里容忍误差的设定与坐标的范围相关
    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double t1 = begin;
    double t2 = end;
    double h = t2-t1;
    romberg[0][0] = h*(compute_f(t1)+compute_f(t2))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f(t1+(2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
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
    std::cout<<"\n";
    printf("The arclength from %lf to %lf is %lf",t1,t2,result);
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),std::ostream_iterator<double>(std::cout," "));

    // }
    return result;
}

double segment::dist(const Point &p1, const Point &p2)
{
    return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
}

Point segment::computePoint(double &t)
{
    Point resultPoint;

    int size = srcdata.size();
    vector<double> coefficient(size, 0);
    coefficient[0] = 1.000;
    double u1 = 1.0 - t;
    //!
    //! De Casteljau递推算法
    //!
    for (int j = 1; j <= size - 1; j++) {
        double saved = 0.0;
        for (int k = 0; k < j; k++){
            double temp = coefficient[k];
            coefficient[k] = saved + u1 * temp;
            saved = t * temp;
        }
        coefficient[j] = saved;
    }
    for (int i = 0; i < size; i++) {
        Point point = srcdata.at(i);
        resultPoint = resultPoint + point * coefficient[i];
    }

    return resultPoint;
}

Point segment::computePoint1(double &t)
{
        Point resultPoint;

    int size = srcdata1.size();
    vector<double> coefficient(size, 0);
    coefficient[0] = 1.000;
    double u1 = 1.0 - t;
    //!
    //! De Casteljau递推算法
    //!
    for (int j = 1; j <= size - 1; j++) {
        double saved = 0.0;
        for (int k = 0; k < j; k++){
            double temp = coefficient[k];
            coefficient[k] = saved + u1 * temp;
            saved = t * temp;
        }
        coefficient[j] = saved;
    }
    for (int i = 0; i < size; i++) {
        Point point = srcdata1.at(i);
        resultPoint = resultPoint + point * coefficient[i];
    }

    return resultPoint;
}

vector<Point> segment::outSegData()
{
    return segdata;
}

vector<Point> segment::outSrcData()
{
    return srcdata;
}

vector<Point> segment::outDesData()
{
    return desdata;
}

double segment::outCrossTime()
{
    return selfcrosstime;
}
