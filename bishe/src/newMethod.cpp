/*
* �µ����
* ���룺ԭ���ĵ���
* for each line:
*	���û����������
*		ֱ��ִ��noisy Model
*	�������������
*		��¼�������巵�ص��ź�ǿ��
*		���ݷֲ������ж�·�����Ƿ������
*		�����
*			ִ��noisy Model
*			ȡ��η�����ź�ǿ�Ⱥ����巴����ź�ǿ�ȵ����ֵ
*			����ֵ�Ƚ�
*		���û��
*			�������巵�ص��ź�ǿ��
*/
#include"baseClass.h"
#include<math.h>
#include<vector>
int threshold;//��ֵ

int main()
{
	std::vector<CloudPoint> pointCloud;//����ĵ���
	int R=1;//������
	Calculator::RainIntensity = R;
	for (auto& currPoint : pointCloud)
	{
		//�����Ǽ�������û���������ƣ�����Ķ�����һ����
		///int finalIntensity=0;
		double pnResult;
		int n = 0;
		//����û����������
		//����ʹ��pn(n)>1������n,������ν��еݹ�ļ��㣬�����Ҳ���Ҫ�ݹ飬�ȼ����nȻ��ѭ������Ҫ��һ��
		while (Calculator::Pn(R, n + 1) > 1)
		{
			n++;
		}
		//������ε�˥���ͷ���
		Calculator::ReflectionAndDecay(rainDrop::RainDropGenerate(n,currPoint), currPoint);
		if (currPoint.intensity < threshold)currPoint.intensity = -1;
	}

}

