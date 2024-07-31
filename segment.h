#include <vector>

using namespace std;

struct Point{
    double x;
    double y;
    //      在这里可以直接定义点的模长以及叉乘了
    Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
    Point operator-()const{
        return Point(-x,-y);
    }
    friend Point operator*(const Point& p,double scalar){
        return {p.x*scalar,p.y*scalar};
    }
    friend Point operator*(double scalar,const Point& p){
        return {scalar*p.x,scalar*p.y};
    }
    friend Point operator*(const Point& p1,const Point& p2){
        return {p1.x*p2.x,p1.y*p2.y};
    }
    friend Point operator+(double scalar,const Point& p){
        return {p.x+scalar,p.y+scalar};
    }
    friend Point operator+(const Point& p,double scalar){
        return scalar+p;
    }
    friend Point operator+(const Point& p1,const Point& p2){
        return {p1.x+p2.x,p1.y+p2.y};
    }
    friend Point operator-(const Point& p1,const Point& p2){
        return {p1.x-p2.x,p1.y-p2.y};
    }
    bool operator==(const Point& other) const {
        return (x == other.x) && (y == other.y);
    }
    Point& operator+=(const Point& p1){
        this->x += p1.x;
        this->y += p1.y;
        return *this;
    }
};
class segment{
public:
    segment();
    segment(const int& count);
    ~segment();
    void setSrcData();                      //      随机设置贝塞尔曲线的控制点
    void getDesData();                      //      利用矩阵方法计算贝塞尔曲线
    void setCoeficient();                   //      设置贝塞尔曲线关于时间t的参数方程的参数
    void setCoeficientG();                  //      设置函数 g 的参数
    void setCoeficientBPrime();

    void computeArcLength();                //      计算贝塞尔曲线的弧长

    void getIsoTime();                      //      计算贝塞尔曲线等距点的时间值
    void NewtonIterator();                  //      牛顿迭代方法计算等距点的时间值
    void computeIsoPoint();
    void computeCrossTime();                //      计算贝塞尔曲线自交的交点
    double compute_f(const double &t);      //      计算 f 对应的函数值
    double compute_gradf(const double &t);  //      计算 f 的导数
    double computeSubArcLength(const double &begin);//计算从0到t时间的曲线弧长
    double computeSubArcLength(const double &begin,const double &end);//计算从begin到end时间的曲线弧长
    Point computePoint(double &t);
    vector<Point> outSegData();//    提供公共接口访问私有数据
    vector<Point> outSrcData();//    提供公共接口访问私有数据
    vector<Point> outDesData();//    提供公共接口访问私有数据
    double outCrossTime();//提供公共接口访问 crosstime 变量
private:
    vector<Point> srcdata;
    vector<Point> desdata;
    vector<Point> segdata;

    vector<double> segtime;

    vector<Point> coefficient;
    vector<Point> coefBPrime;
    vector<double> coefficientG;
    vector<int> pascaTri{0,1,3,3,1};

    bool corsslabel = false;
    double crosstime = -1;

    double arclength = 0;
    double precis = 1e-2;
    double delta = 1e-5;
    double epslion = 1e-5;
    double tolerenceerror = 1e-7;

    int segcount = 50;
    


};
