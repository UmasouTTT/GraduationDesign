#pragma once

#include "perm.h"

using namespace std;
using namespace perm_struct;

class three_perm
{
public:
	three_perm();
	~three_perm();

	static const int predict_worest_energy = -45;//预计最低值

private:
	//长度为n的构型的历史平均权重(需要初始化)
	double *average_weights;
	//长度为n的构型的数量(需要初始化)
	double *weights_numbers;
	//各分支具体构型
	point *configurations_point_three;
	char *configurations_class_three;
	//当前构型能量
	int present_energy;
	//最低能量
	int lowest_energy = 0;
	//最低能量构型
	point *lowest_configurations_point;
	char *lowest_configurations_class;

	//perm最低能量
	int perm_lowest_energy = 0;
	//perm最低能量构型
	point *perm_lowest_configurations_point;
	char *perm_lowest_configurations_class;
	//初次迭代获取的最低能量值
	int lowest_energy_first;

private:
	//求小值
	template <typename T>
	T Min(T num1, T num2) {
		if (num1 < num2) {
			return num1;
		}
		return num2;
	}
	//数组赋值
	template <typename T>
	void ArrayAssignment(T number1[], T number2[], int length) {
		for (size_t i = 0; i < length; i++) {
			number1[i] = number2[i];
		}
	}
	//判断该坐标是否已经被使用
	bool IsThisPositionAlreadyOccupied(point p, int n);
	//计算合法的动作数
	int LegalActions(point p, int n);
	//计算两个点之间的距离
	float DistenceBetweenPoints(point point1, point point2);
	//重构计算合法动作函数，提高计算速率
	int LegalActions(point p, point legal_actions[4], int n);
	//测试运算结果是否正确
	bool TestResultIsSatisfied(int target_energy, int length);
	//初始化（初始化变元，前两个值为定值）
	void InitConfig(string &input, point &p, double &weight);	
	//迭代过程
	void CircleCalculate(int n, int whole_length, point p_before, double weight, string input);
	//迭代过程(只考虑权重)
	void CircleCalculateByWeight(int n, int whole_length, point p_before, double weight, string input);
	//根据选择的更新全局变量
	void  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase);
	//更新Cn,Zn
	void UpdateAverageWeight(double w, int n);
	//计算权重
	double CalculateWeight(double w, int energy_increase);
	//计算能量增量
	int EnergyIncrease(point p, char type, point p_before, int n);
	//更新临时参数
	void UpdateTempVariables(point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase);
	//计算当前构型能量值
	int CalculatePresentConfigEnergy(point _configurations_point[], char _configurations_class[], int length);
	//计算初期构型的权重
	double CalculateStartConfigWeight(point _configurations_point[], char _configurations_class[], int length);
public:
	//算法
	void StartCalculate(string input);
	//改进算法α
	void StartCalculateImproveFirst(string input);
	//迭代改进算法
	void Branch_choose_improve_1(string input, int time);
	//迭代改进算法2
	void StartCalculateByWeight(string input);
	//循环处理迭代改进算法2
	void CircleAlgripham2(string input, int times);
};

