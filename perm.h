#pragma once
#include<stdlib.h>
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
	perm(int predict_energy);
	perm(int predict_energy, bool isSoultionSaved);
	~perm();

	static const int max_size_of_input = 100;//最长链长
	static const int max_size_of_legal_input = 4;//每个点的最多合法方向
	static const int max_size_of_possibleConditions = 24;//
	static const double T;
	static const double C0;
	static const int Z0;
	static const int C;
	static const float part_config_for_save;//初始迭代所保存构型比例
	static const int target_lowest_energy;//目标最低构型

private:
	//当前最大分支标识号
	int max_tag = 0;
	//权重算术平均值(需要初始化)
	//double *average_weights = new double[max_size_of_input];
	double *average_weights; 
	//长度为n的构型的数量(需要初始化)
	//double *weights_numbers = new double[max_size_of_input];
	double *weights_numbers;
	//各分支具体构型
	point *configurations_point;
	char *configurations_class;
	//各分支当前构型能量
	int present_energy;
	//vector<int>present_energy;
	//最低能量
	int lowest_energy = 0;
	//最低能量构型
	point *lowest_configurations_point;
	char *lowest_configurations_class;
	//用于获取合法组合集合

	int **combinations_result;
	point *choose_actions_3;
	//最低能量构型数量
	int num_of_lowestConfigurations;
	//默认最差构型能量，由第一次随机遍历产生
	int worest_energy;
	//定义Perm迭代出构型长度
	int choose_config_length = 0;
	//记录目前出现的最优构型
	point **best_config_ever;
	//记录目前最优构型数量
	int best_config_num;
	//生成初始权重
	bool InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, int _energy, point points[max_size_of_input], char _points[max_size_of_input]);
	//用于记录开始进行剪枝时的链长
	int beginPuningTheBranch;
	//用于判定是否已经开始剪枝
	bool isPuneBegin;
	//用于判定是否需要记录最优解分支
	bool isSolutionNeedToBeSaved;

public:
	//获取平均权重
	void GetAverageWeight(double _average_weights[max_size_of_input]);
	//获取长度为n的构型的数量
	void GetWeightNumber(double _weights_numbers[max_size_of_input]);
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
	//获取最低构型数
	int GetNumOfLowestConfigurations() { return num_of_lowestConfigurations; }
	//获取最低可用构型数量
	int GetNumOfLowestConfigAndCanUse() { return best_config_num; }
	//获取当前最优局部构型
	void GetCurrentOptimalLocalConfiguration(point _best_config_ever[100][perm::max_size_of_input], int num_of_best_config, int length);
	//获取指定链长平均权重(for test)
	double GetTargetLengthWeight(int _index) { return average_weights[_index - 1]; }
	//获取开始剪枝时的链长
	int GetTheLengthStartPuneing() { return beginPuningTheBranch; }

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
	int LegalActions(point p, point legal_actions_p[4], int legal_actions_t[4], int n);
	//**************计算好度*************
	double CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n);
	//**************计算权重*************
	double CalculateWeight(double w, point p, char type, int energy_increase, point p_before);
	//计算生长比例系数
	double CalculatingLengthCoefficient(int n, int length);
	//**************计算预计权重及各个动作的好度（避免重复计算）*************(由于内容较多，分两步进行)
	double CalculatePredictWeightMid(double w, point p_before, char type, double good_degrees[4], int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4], int n);
	double CalculatePredictWeight(double w, point p_before, char type, double good_degrees[4], int n, int length, int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4]);
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
	point GetNextActionByGoodDegrees(point p_before, double good_degrees[4]);
	//递归计算排列组合
	void CalculationCombinations(int offset, int k);
	//获取可能组合数
	void GetCombinations(int(&legal_actions)[4], int num, int num_of_legal_actions, int &num_of_result);
	//根据数值获取相应的动作
	void GetActionsByNum(int index, point p_before, int length_of_result);
	//************按照好度概率随机选择动作集合*****************
    void ChooseActionsGroupByGoodDegrees(int k, double good_degrees[4], point p_before);
	//测试运算结果是否正确
	bool TestResultIsSatisfied(int target_energy, int length);	
	//初始化（初始化变元，前两个值为定值,需要生成初始权重）
	void InitConfig(string &input, point &p, double &weight);
	//初始化（初始化变元，前两个值为定值，无需生成初始权重）
	void InitConfigWithoutInitWeight(string &input, point &p, double &weight);
	//计算最大长度
	void CalculateMaxSize(int length);
	//初始化全局变元
	void InitGlobalVariable(string input);
	//获取动作序号
	int GetActionsNum(point actions[4], point action);

public:
	//迭代算法
	void StartCalculate(string input, int num_of_circle);
	//迭代算法(只考虑权重)
	void StartCalculateOnlyWeight(string input, int num_of_circle);
	//迭代计算各分支情况（记录分支数量）
	void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input);
	//迭代计算各分支情况（不记录分支数量）
	void CalculationProcessOnlyConsiderWeight(int n, int whole_length, int tag, point p_before, double weight, string input);
	//多次计算各分支情况
	void CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]);
	//多次计算各分支情况(只考虑权重)
	void CircleCalculateProcessOnlyConsiderWeight(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]);
	//获取指定链长平均权重值
	double GetAverageWeightWithLengthN(int n) { return average_weights[n-1]; }

	//非算法部分功能
private:
	//向历史最优构型中添加新构型
	void AddNewConfigToBestConfigEver(point points[max_size_of_input], int length);
	//判断两种构型是否可以认为一致（目前认为两种构型保持在所选链长一致即认为两种构型一致）
	bool IsTwoConfigTheSame(point points[max_size_of_input], point _points[max_size_of_input], int start, int end);
};

