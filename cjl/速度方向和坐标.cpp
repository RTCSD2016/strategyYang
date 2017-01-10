/*
�����֣�
���꣨x1��y1��������m1���뾶r1����ײǰ�ٶȣ�vx1��vy1�����ٶ�v1����ײ���ٶȣ�vx1_��vy1_��
����
���꣨x2��y2������m2���뾶r1����ײǰ�ٶȣ�vx2��vy2��,��ײ���ٶȣ�vx2_��vy2_��
�����������꣨x3��y3��

���룺����ʱ�����������ꡣ����������ֻ����ٶȴ�С����ͻ��������ꡣ
*/
#include <iostream>
#include <algorithm>
#include <cmath>

double get_angel(x2,y2)//��ȡ�������ٶȷ���
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
double get_hitx(beta)//��ȡ�����ֻ���ʱ��x��������
{
	double sin, x1;

	cos = 1 / sqrt(1 + pow(beta, 2));
	x1 = x2 - (r1 + r2)*cos;
	return x1;
}
double get_hity(beta)//��ȡ�����ֻ���ʱ��y��������
{
	double cos, y1;

	sin = beta / sqrt(1 + pow(beta, 2));
	if (y2 >= y3)//�жϱ�������������ĵ�y�����λ��
	{
		y1 = y2 + (r1 + r2)*sin;
	}
	else{
		y1 = y2 - (r1 + r2)*sin;
	}
		return y1;
}

double get_result(a,b,c)//�������
{
	double a, b, c;
	t = b*b - 4 * a*c;
	if (t == 0) //�ɸ����б�ʽ������ִ��������֧
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
	cout << "�˷����޸���" << endl;
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
