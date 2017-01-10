/*
击球手：
坐标（x1，y1），质量m1，半径r1，碰撞前速度（vx1，vy1）合速度v1，碰撞后速度（vx1_，vy1_）
冰球：
坐标（x2，y2）质量m2，半径r1，碰撞前速度（vx2，vy2）,碰撞后速度（vx2_，vy2_）
球门中心坐标（x3，y3）

输入：击球时冰球中心坐标。输出：击球手击球速度大小方向和击球手坐标。
*/
#include <iostream>
#include <algorithm>
#include <cmath>

double get_angel(x2,y2)//获取击球手速度方向
{
	double alpha, beta;
	double a, b, c, complex_term;

	alpha = (y3 - y2) / (x3 - x2);
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
	x1 = x2 - (r1 + r2)*cos;
	return x1;
}
double get_hity(beta)//获取击球手击球时的y方向坐标
{
	double cos, y1;

	sin = beta / sqrt(1 + pow(beta, 2));
	if (y2 >= y3)//判断冰球相对球门中心的y轴相对位置
	{
		y1 = y2 + (r1 + r2)*sin;
	}
	else{
		y1 = y2 - (r1 + r2)*sin;
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
	x2 = ((-b - (sqrt(t))) / (2 * a));
	return max(x1, x2);
}
