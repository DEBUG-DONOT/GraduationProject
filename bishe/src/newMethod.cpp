/*
* 新的设计
* 输入：原来的点云
* for each line:
*	如果没有碰到物体
*		直接执行noisy Model
*	如果碰到物体了
*		记录碰到物体返回的信号强度
*		根据分布函数判断路径上是否有雨滴
*		如果有
*			执行noisy Model
*			取雨滴反射的信号强度和物体反射的信号强度的最大值
*			和阈值比较
*		如果没有
*			返回物体返回的信号强度
*/
#include"baseClass.h"
#include<math.h>
#include<vector>
int threshold;//阈值

int main()
{
	std::vector<CloudPoint> pointCloud;//输入的点云
	int R=1;//降雨量
	Calculator::RainIntensity = R;
	for (auto& currPoint : pointCloud)
	{
		//无论是激光束有没有碰到点云，计算的东西是一样的
		///int finalIntensity=0;
		double pnResult;
		int n = 0;
		//光束没有碰到物体
		//计算使得pn(n)>1的最大的n,对于雨滴进行递归的计算，但是我不想要递归，先计算出n然后循环计算要好一点
		while (Calculator::Pn(R, n + 1) > 1)
		{
			n++;
		}
		//计算雨滴的衰减和反射
		Calculator::ReflectionAndDecay(rainDrop::RainDropGenerate(n,currPoint), currPoint);
		if (currPoint.intensity < threshold)currPoint.intensity = -1;
	}

}

