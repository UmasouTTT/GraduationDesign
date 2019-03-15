#pragma once

#include "perm.h"

using namespace std;
using namespace perm_struct;

class three_perm
{
public:
	three_perm();
	~three_perm();



private:
	double average_weights[perm::max_size_of_input];
	//长度为n的构型的数量(需要初始化)
	double weights_numbers[perm::max_size_of_input];
	//各分支具体构型
	point configurations_point[perm::max_size_of_input];
	char configurations_class[perm::max_size_of_input];
	//当前构型能量
	int present_energy;
	//最低能量
	int lowest_energy = 0;
	//最低能量构型
	point lowest_configurations_point[perm::max_size_of_input];
	char lowest_configurations_class[perm::max_size_of_input];

	//perm最低能量
	int perm_lowest_energy = 0;
	//perm最低能量构型
	point perm_lowest_configurations_point[perm::max_size_of_input];
	char perm_lowest_configurations_class[perm::max_size_of_input];

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
	int LegalActions(point p, vector<point> &legal_action, int n);
	//测试运算结果是否正确
	bool TestResultIsSatisfied(int target_energy, int length);
	//初始化（初始化变元，前两个值为定值）
	void InitConfig(string &input, point &p, double &weight);	
	//迭代过程
	void CircleCalculate(int n, int whole_length, point p_before, double weight, string input);
	//根据选择的更新全局变量
	void  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, double _average_weights[], double _weights_numbers[]);
	//更新Cn,Zn
	void UpdateAverageWeightByThree(double w, int n, double _average_weights[], double _weights_numbers[]);
	//计算权重
	double CalculateWeight(double w, int energy_increase);
	//计算能量增量
	int EnergyIncrease(point p, char type, point p_before, int n);
	//更新临时参数
	void UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase);
public:
	//算法
	void StartCalculate(string input);
};

