#pragma once
#include<vector>
class rainDrop
{
public:
	double x, y, z;
	double diameter;
	double distance;//和激光雷达的距离
	static std::vector<rainDrop> RainDropGenerate(int n, CloudPoint cp);
private:
};

class CloudPoint
{
public:
	double x, y, z;
	double intensity;
	
};

class Calculator
{
public:
	static double RainIntensity;
	static double Pn(int R, int n);
	static double simpson(double a,double b,int R);
	static void ReflectionAndDecay(const std::vector<rainDrop>& rd, CloudPoint& cp);
	static double randomDoubleGen(double lower_bound, double upper_bound);
	static double rainModel(double diameter, int R);
private:
	static double calculAvgRainNum(int R);
	static double reflection(const rainDrop& rd, CloudPoint& cp);
	static void Decay(const rainDrop& rd, CloudPoint& cp);
};






