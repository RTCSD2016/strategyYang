#include <iostream> 
#include <string>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include <algorithm>
#include <cmath>
using namespace std;

double xa, ya, va, a;  //输入冰球信息

double xx1, yy1;   //冰球轨迹与上下壁交点坐标
double xx2, yy2;   //冰球轨迹与左壁交点坐标

double xhit, yhit;     //击球点
double xhit1, yhit1;   //轨迹与击球区域边界第一个交点
double xhit_sf, yhit_sf;  //轨迹与算法选择线交点
double xrandom, yrandom;  //随机确定点

double x3 = 120, y3 = 64;   //桌面大小
double x1 = 20, y11 = 16, x2 = 60, y2 = 48;    //击球区左边界，下边界，右边界,上边界
double y4 = 24, y5 = 40;      //球门范围


/*计算冰球轨迹与桌面壁的碰撞点*/
void root(double xa, double ya, double va, double a){
	double t1, t2; 
	if (sin(a) > 0) {
		t1 = (y3 - ya) / (va*sin(a));
		yy1 = y3;
	}
	else {
		t1 = ya / (va*sin(a));                    //确定第一次碰撞的事件t1
		yy1 = 0;                              
	}
	xx1 = xa + t1*va*cos(a);                          //与上下壁碰撞点（xx1,yy1）
	t2 = xx1 / (va*cos(a));
	yy2 = yy1 + t2*va*sin(a);
	xx2 = 0;                                          //与左壁碰撞点（xx2,yy2）
}

/*计算与击球区域边界第一个交点*/
void point1() {
	double t;
	double x22, y22;
	t = (x2 - xa) / (va*cos(a));
	y22 = ya + t*va*sin(a);
	if (y22<y2&&y22>y11) {
		xhit1 = x2;
		yhit1 = y22;
	}
	else {
		t = y11 / (va*sin(a));
		x22 = xx1 + t*va*cos(a);
		xhit1 = x22;
		if (sin(a) > 0)
			yhit1 = y2;
		else
			yhit1 = y11;
	}
}

/*判断是否会进球*/

bool hitin() {
	if (yy2<y5&&yy2>y4)
		return true;
	else
		return false;
}

/*计算与算法选择线交点(xhit_sf,yhit_sf)*/
void point2() {
	double t;
	t = (x1 - xx1) / (va*cos(a));
	xhit_sf = x1;
	if (sin(a) > 0)
		yhit_sf = y3 - va*sin(a)*t;
	else
		yhit_sf = va*sin(a)*t;                          //(xhit_sf,yhit_sf)
}

/*产生随机击球点*/

void pointrandom(){ 
	double num;
    double t;
    srand(time(NULL));
    num = (rand() % 100)*0.01;
	xrandom = x1 + num*(x2 - x1);
	t = (xrandom - x1) / va*cos(a);
	yrandom = yy2 - va*sin(a)*t;
}

/*判断击球点*/
void hitpoint() {
	bool sf;
	point1();
	if (xhit1 > x1) {
		xhit = xhit1;
		yhit = yhit1;
	}
	else {
		sf = hitin();
		if (sf) {
			point2();
			xhit = xhit_sf;
			yhit = yhit_sf;
		}
		else {
			pointrandom();
			xhit = xrandom;
			yhit = yrandom;
		}
	}
}

void main() {
	cout << "请输入冰球当前信息" << endl;
	cout << "x=";
	cin >> xa;
	cout << "y=";
	cin >> ya;
	cout << "va=";
	cin >> va;
	cout << "a=";
	cin >> a;
	root(xa, ya, va, a);
	hitpoint();
	//cout << xx1 << yy1 << xx2 << yy2;
	cout << "(xhit,yhit)=(" << xhit << ","<<yhit << ")";
}
double get_angel(xhit,yhit)//获取击球手速度方向
{
	double alpha, beta;
	double a, b, c, complex_term;

	alpha = (y3 - yhit) / (x3 - xhit);
	complex_term = pow(2 * m1*v1 / ((m2 - m1)(alpha*vx2 - vy2)), 2);
	a = 1 + complex_term;
	b = -2 * complex_term*alpha;
	c = pow(alpha, 2)*complex_term - 1;
	beta = get_result(a, b, c);
	
	return beta;
}
double get_hitx(beta)//获取击球手击球时的x方向坐标
{
	double sin, x1;

	cos = 1 / sqrt(1 + pow(beta, 2));
	x1 = xhit - (r1 + r2)*cos;
	return x1;
}
double get_hity(beta)//获取击球手击球时的y方向坐标
{
	double cos, y1;

	sin = beta / sqrt(1 + pow(beta, 2));
	if (yhit >= y3)//判断冰球相对球门中心的y轴相对位置
	{
		y1 = yhit + (r1 + r2)*sin;
	}
	else{
		y1 = yhit - (r1 + r2)*sin;
	}
		return y1;
}

double get_result(a,b,c)//求根函数
{
	double a, b, c;
	t = b*b - 4 * a*c;
	if (t == 0) //由根的判别式来决定执行哪条分支
	{
		
	  f2(a, b);
    
	}
	else if (t < 0)
	{
		f1();
	}
	else
	{
		f3(a, b, c);
	}
}
void f1()
{
	cout << "此方程无根！" << endl;
}
double f2(double a, double b)
{
	x = -b / (2 * a);
	return x;
}
void f3(double a, double b, double c)
{
	x1 = ((-b + (sqrt(t))) / (2 * a));
	xhit = ((-b - (sqrt(t))) / (2 * a));
	return max(x1, xhit);
}

