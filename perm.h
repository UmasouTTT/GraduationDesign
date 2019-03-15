#pragma once
#include "stdafx.h"
#include<stdlib.h>
#include<map>
#include<unordered_map>
#include<vector>
#include<math.h>
#include<cmath>
#include<string>
#include<time.h>
#include<iostream>

namespace perm_struct {
	struct point {
		int x;
		int y;
		bool operator < (const point &p) const {
			return x < p.x || (x == p.x && y < p.y);
		}
	};
}


#define random(a,b) (((double)rand()/RAND_MAX)*(b - a) + a)
#define judge_is_zero 0.000001

using namespace std;
using namespace perm_struct;

class perm
{
public:
	perm();
	~perm();

	static const int max_size_of_input = 100;//最长链长
	static const int max_size_of_legal_input = 4;//每个点的最多合法方向
	static const int max_size_of_possibleConditions = 24;//
	static const double T;
	//static const double C0;
	static const int Z0;
	static const int C;
	static const double MAX_DOUBLE;

private:
	//用于控制人口
	double C0;
	//当前最大分支标识号
	int max_tag = 0;
	//权重算术平均值(需要初始化)
	//double *average_weights = new double[max_size_of_input];
	double average_weights[max_size_of_input];
	//长度为n的构型的数量(需要初始化)
	//double *weights_numbers = new double[max_size_of_input];
	double weights_numbers[max_size_of_input];
	//各分支具体构型
	point configurations_point[max_size_of_input];
	char configurations_class[max_size_of_input];
	//各分支当前构型能量
	int present_energy;
	//vector<int>present_energy;
	//最低能量
	int lowest_energy = 0;
	//最低能量构型
	point lowest_configurations_point[max_size_of_input];
	char lowest_configurations_class[max_size_of_input];
	//用于获取合法组合集合
	vector<int> input_numbers;
	vector<int> combination_one;
	vector<vector<int>> combination_result;
	//最低能量构型数量
	int num_of_lowestConfigurations;

public:
	//获取能量值
	int GetEnergy() { return lowest_energy; }
	//获取最低能量构型点坐标
	void GetPointPosition(point points[max_size_of_input]);
	//获取最低能量构型点类型
	void GetPoint(char points[max_size_of_input]);
	//设置能量值
	void SetEnergy(int energy) { present_energy = energy; }
	//设置最低能量构型点坐标
	void SetPointPosition(point points[max_size_of_input]);
	//设置最低能量构型点类型
	void SetPoint(char points[max_size_of_input]);
	//设置权重算术平均值
	void SetAverageWeight(double _average_weights[max_size_of_input]);
	//设置长度为n的构型的数量
	void SetThisWeightNumber(double _weights_numbers[max_size_of_input]);
	//获取权重算术平均值
	void GetAverageWeight(double _average_weights[max_size_of_input]);
	//获取长度为n的构型的数量
	void GetThisWeightNumber(double _weights_numbers[max_size_of_input]);
	//获取最低构型数
	int GetNumOfLowestConfigurations() { return num_of_lowestConfigurations; }
	//设置人口繁衍情况
	void SetPopulation(double _C0) { C0 = _C0; }

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
	//计算两个点之间的距离
	float DistenceBetweenPoints(point point1, point point2);
	//计算能量增量
	int EnergyIncrease(point p, char type, point p_before, int n);
	//判断该坐标是否已经被使用
	bool IsThisPositionAlreadyOccupied(point p, int n);
	//计算合法的动作数
	int LegalActions(point p, int n);
	//重构计算合法动作函数，提高计算速率
	int LegalActions(point p, vector<pair<int, point>> &legal_actions, int n);
	//**************计算好度*************
	double CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n);
	//**************计算权重*************
	double CalculateWeight(double w, point p, char type, int energy_increase, point p_before);
	//计算生长比例系数
	double CalculatingLengthCoefficient(int n, int length);
	//**************计算预计权重及各个动作的好度（避免重复计算）*************(由于内容较多，分两步进行)
	double CalculatePredictWeightMid(double w, point p_before, char type, vector<double> &good_degrees, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase, int n);
	double CalculatePredictWeight(double w, point p_before, char type, vector<double> &good_degrees, int n, int length, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase);
	//*************更新Cn,Zn***************
	void UpdateAverageWeight(double w, int n);
	//***************计算上门限***********
	double CalculateUpperThreshold(int n);
	//**************计算下门限***********
	double CalculateLowerThreshold(double upper_threshold);
	//*****************根据选择的更新全局变量***************
	int  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, point point_before[], char type_before[]);
	int  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase);
	//**********************按照概率生成随机动作********************************
	point GetNextActionByGoodDegrees(point p_before, vector<double> &good_degrees);
	//递归计算排列组合
	void CalculationCombinations(int offset, int k);
	//获取可能组合数
	vector<vector<int>>GetCombinations(vector<int> &legal_actions, int num);
	//根据数值获取相应的动作
	vector<point> GetActionsByNum(vector<int> &numbers, point p_before);
	//************按照好度概率随机选择动作集合*****************
	vector<point> ChooseActionsGroupByGoodDegrees(int k, vector<double> &good_degrees, point p_before);
	//测试运算结果是否正确
	bool TestResultIsSatisfied(int target_energy, int length);	
	//初始化（初始化变元，前两个值为定值）
	void InitConfig(string &input, point &p, double &weight, int tag);
	//计算最大长度
	void CalculateMaxSize(int length);
	//初始化全局变元
	void InitGlobalVariable(string input);
	//生成初始权重
	void InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, double _average_weights[max_size_of_input], int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _weights_numbers[max_size_of_input]);

public:
	//迭代算法
	void StartCalculate(string input, int num_of_circle);
	//迭代计算各分支情况
	void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input);
	//多次计算各分支情况
	void CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, double _average_weights[max_size_of_input], int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _weights_numbers[max_size_of_input]);
};

