#include"baseClass.h"
#include<math.h>
#include<random>
const double pi = 3.1415926;
const double Epsilon = 1e-2;//误差阈值
const double Alaser;
const double LaserL;
const double c;//光速
const double T_T;
const double T_R;
const double RouRain;//雨滴的平均光反射系数
const double alpha;//介质散射系数
const std::vector<double> LaserPos{ 1.0,1.0,1.0 };//雷达位置
/*@param diameter 雨滴直径
* @param R 降雨量
*/
double Calculator::rainModel(double diameter, int R)
{
	//N(D)
	//F-L模型
	double N_t = 172 * pow(R, 0.22);
	double D_g = 0.72 * pow(R, 0.23);
	double phi = 1.43 - 3 * 1e-4 * R;
	double result = N_t / (pow(2 * pi, 0.5) * diameter * log(phi));
	result = result * exp(-log(diameter / D_g) / (2 * pow(log(phi), 2)));
	return result;
}

double Calculator::simpson(double a, double b,int R)
{
	double I2n = 0, h = b - a, T2n=h*(rainModel(a,R)+rainModel(b,R))/2, In=T2n,Tn;
	for (int n = 1; abs(I2n - In) > Epsilon; n += n, h /= 2.0)
	{
		In = I2n;
		Tn = T2n;
		double sigma = 0;
		for (int k = 0; k < n; k++)
		{
			sigma += rainModel(a + (k + 0.5) * h,R);
		}
		T2n = (Tn + h * sigma) / 2;
		I2n = (4 * T2n - Tn) / 3;
	}
	return I2n;
}

double Calculator::reflection(const rainDrop& rd, CloudPoint& cp)
{
	//计算反射的强度并返回
	double betaR = 0;
	double v = abs(rd.x - LaserPos[0]) * abs(rd.x - LaserPos[0]) + abs(rd.y - LaserPos[1]) * abs(rd.y - LaserPos[1]) + abs(rd.z - LaserPos[2]) * abs(rd.z - LaserPos[2]);
	v = pow(v, 0.5);
	betaR = cp.intensity * pi * c * RouRain * rd.diameter * rd.diameter * T_R * T_T * exp(alpha * v) / (8 * v * v);
	return betaR;
}

void Calculator::Decay(const rainDrop& rd, CloudPoint& cp)
{
	//衰减
	

}

double Calculator::calculAvgRainNum(int R)
{
	//模型处于0.5mm-6mm
	//使用simpson方法进行积分
	auto avgRainNum=simpson(0.5, 6, R);
	//考虑激光束为底面积为Alaser，高为LaserL的圆柱体
	double avgNum = avgRainNum * Alaser * LaserL;//这里的LaserL要改成到物体距离，如果没有碰到物体的话就是侦测的最远距离
	return avgNum;
}

double Calculator::Pn(int R, int n)
{
	//泊松分布
	//n是雨滴的数量
	//pn<1就当作路径上没有雨滴
	double mu=calculAvgRainNum(R);
	double nJieCheng = 1;
	for (int i = 1; i <= n; i++)
	{
		nJieCheng *= i;
	}
	double result = exp(-mu) * (pow(mu, n)) / nJieCheng;
	return result;
}

void Calculator::ReflectionAndDecay(const std::vector<rainDrop>& rd, CloudPoint& cp)
{
	//输入一个雨滴和一个点云，计算雨滴对于点云的衰减，
	//同时在该点云的反射出来的强度和原来的强度中选一个最大的作为该点云的新的强度
	//计算反射和衰减的时候需要粒子的位置
	//衰减要计算多次，反射应该只用计算一次就好，因为衰减要计算多次，所以n的设计是没问题的
	auto reflectIntensity = Calculator::reflection(rd[0], cp);
	for (auto& curr : rd)
	{
		//对于每一个雨滴计算衰减和反射
		Calculator::Decay(curr, cp);
	}
	cp.intensity = std::max(cp.intensity, reflectIntensity);
}

double Calculator::randomDoubleGen(double lower_bound, double upper_bound)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	// 定义分布
	std::uniform_real_distribution< > dis(lower_bound, upper_bound);
	// 生成随机数
	double random_num = dis(gen);
	return random_num;
}

std::vector<rainDrop> rainDrop::RainDropGenerate(int n,CloudPoint cp)
{
	//按照雨滴谱分布随机生成一个雨滴
	//要求是在激光雷达的线上，可以不用整个雨滴都在里面，只要有一部分在里面就行，这个可以用连线表示
	//随机在0.5-6的直径中生成一个直径，取得对应的N(D)，和积分后的结果做比较，得到一个比值x，在[0-1]之间随机生成一个值y，如果x>y那么就将这个直径作为随机生成的雨滴的直径
	std::vector<rainDrop> rains(n);
	auto rainNumSum = Calculator::simpson(0.5, 6, Calculator::RainIntensity);
	for (int i = 0; i < n; i++)
	{
		bool flag = true;
		double dia;
		while (flag)
		{
			double rainDiameter = Calculator::randomDoubleGen(0.5, 6.0);
			auto Nd = Calculator::rainModel(rainDiameter, Calculator::RainIntensity);
			double x = Calculator::randomDoubleGen(0, 1);
			if (Nd / rainNumSum > x)
			{
				dia = rainDiameter;
				flag = false;
			}
		}
		//随机生成距离，最后按照距离排序
		double dis = abs(cp.x - LaserPos[0]) * abs(cp.x - LaserPos[0]) + abs(cp.y - LaserPos[1]) * abs(cp.y - LaserPos[1]) + abs(cp.z - LaserPos[2]) * abs(cp.z - LaserPos[2]);
		dis = pow(dis, 0.5);
		double rainDis = Calculator::randomDoubleGen(0.1, dis);
		rainDrop tempRd;
		tempRd.distance = rainDis;
		//有角度和距离，可以计算出坐标
		
	}

	return rains;
}

