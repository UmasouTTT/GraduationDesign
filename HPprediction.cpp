// HPprediction.cpp: 定义控制台应用程序的入口点。
//算法

#include "stdafx.h"
#include<stdlib.h>
#include<map>
#include<vector>
#include<math.h>
#include<cmath>
#include<string>
#include<time.h>

#define random(a,b) (((double)rand()/RAND_MAX)*(b - a) + a)

using namespace std;

struct point {
	int x;
	int y;
	bool operator < (const point &p) const {
		return x < p.x ;
	}
};

const double T = 0.25;
const int C0 = 1;
const int Z0 = 1;
const int C = 1;
//当前最大分支标识号
int max_tag = 0;
//权重算术平均值(需要初始化)
vector<double>average_weights;
//长度为n的构型的数量(需要初始化)
vector<int>weights_numbers;
//各分支具体构型
vector<map<point, char>>configurations;
//各分支当前构型能量
vector<int>present_energy;
//最低能量
int lowest_energy = 0;
//最低能量构型
vector<point>lowest_energy_points;

//求小值
template <typename T> 
T Min(T num1, T num2) {
	if (num1 < num2) {
		return num1;
	}
	return num2;
}


//计算两个点之间的距离
float DistenceBetweenPoints(point point1, point point2) {
	float result = (float)(point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
	return sqrtf(result);
}

//计算能量增量
int EnergyIncrease(point p, char type, const map<point, char> &points) {
	int result = 0;
	if (type == 'P') {
		return 0;
	}
	//遍历所有节点，判断距离
	map<point, char>::const_iterator iter = points.begin();
	map<point, char>::const_iterator endIter = points.end();
	for (; iter != endIter; iter++) {
		point _point = iter->first;
		char c = iter->second;
		if (c == 'H' && DistenceBetweenPoints(p, _point) == 1) {
			result -= 1;
		}
	}
	return result;
}

//判断该坐标是否已经被使用
bool IsThisPositionAlreadyOccupied(point p, const map<point, char> &points) {
	if (points.find(p) == points.end()) {
		return false;
	}
	return true;
}

//计算合法的动作数
int LegalActions(point p, const map<point, char> &points) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, points)) {
		result += 1;
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, points)) {
		result += 1;
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, points)) {
		result += 1;
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, points)) {
		result += 1;
	}
	return result;
}


//**************计算好度*************
double CalculateGoodResults(point p, char type, const map<point, char> &points) {
	double result = 0.0;
	int actions_later = LegalActions(p, points);
	result += ((double)actions_later + 0.5) * exp(-EnergyIncrease(p, type, points) / T);
	return result;	
}

//**************计算权重*************
double CalculateWeight(double w, point p, char type, const map<point, char> &points, int &energy_increase) {
	energy_increase = EnergyIncrease(p, type, points);
	double result = w * exp(-energy_increase / T);
	return result;
}

//**************计算预计权重及各个动作的好度（避免重复计算）*************
double CalculatePredictWeight(double w, point p_before, char type, const map<point, char> &points, vector<double> &good_degrees) {
	double result = 0;
	int legal_action_numbers = 0;
	//n+1步为上端放置
	point p1(p_before);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, points)) {
		good_degrees.push_back(CalculateGoodResults(p1, type, points));
		result += EnergyIncrease(p1, type, points);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为右端放置
	point p2(p_before);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, points)) {
		good_degrees.push_back(CalculateGoodResults(p2, type, points));
		result += EnergyIncrease(p2, type, points);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为下端放置
	point p3(p_before);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, points)) {
		good_degrees.push_back(CalculateGoodResults(p3, type, points));
		result += EnergyIncrease(p3, type, points);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为左端放置
	point p4(p_before);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, points)) {
		good_degrees.push_back(CalculateGoodResults(p4, type, points));
		result += EnergyIncrease(p4, type, points);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//与前一步权重求积
	double energy_increase_average = result / legal_action_numbers;
	result = w * exp(-energy_increase_average / T);
	return result;
}

//*************更新Cn,Zn***************
void UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n] * weights_numbers[n];
	++weights_numbers[n];
	average_weights[n] = (average_weight_before + w) / weights_numbers[n];
}

//***************计算上门限***********
double CalculateUpperThreshold(int n) {
	double result = C * (average_weights[n] / Z0) * (weights_numbers[n] / C0) * (weights_numbers[n] / C0);
	return result;
}

//**************计算下门限***********
double CalculateLowerThreshold(double upper_threshold) {
	double result = 0.2 * upper_threshold;
	return result;
}

//*******************创建新的分支*****************************
void CreateNewBranch(int father_tag) {
	//分支构型
	map<point, char> before_construction = configurations[father_tag];
	configurations.push_back(before_construction);
	//分支能量
	int before_energy = present_energy[father_tag];
	present_energy.push_back(before_energy);
	//分支标识
	++max_tag;
}

//*****************根据选择的更新全局变量***************
void UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//更新权重算术平均值及该种构型长度的数量
	UpdateAverageWeight(weight, n);
	//更新各分支具体构型
	configurations[tag].insert(make_pair(p, type));
	//更新各分支当前构型能量
	present_energy[tag] += energy_increase;
}

//**********************按照概率生成随机动作********************************
point GetNextActionByGoodDegrees(point p_before, vector<double> &good_degrees) {
	double whole_good_degrees = 0;
	double present_goodD_sum = good_degrees[0];
	for (size_t i = 0; i < good_degrees.size(); i++){
		whole_good_degrees += good_degrees[i];
	}
	double result = random(0, whole_good_degrees);
	if (result >= 0 && result < present_goodD_sum) {
		p_before.y = p_before.y + 1;
		return p_before;
	}
	if (result >= present_goodD_sum && result < (present_goodD_sum + good_degrees[1])) {
		p_before.x = p_before.x + 1;
		return p_before;
	}
	present_goodD_sum += good_degrees[1];
	if (result >= present_goodD_sum && result < (present_goodD_sum + good_degrees[2])) {
		p_before.y = p_before.y - 1;
		return p_before;
	}
	p_before.x = p_before.x - 1;
	return p_before;
}

//迭代计算各分支情况
void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		if (present_energy[tag] < lowest_energy) {
			lowest_energy = present_energy[tag];
			lowest_energy_points.clear();
			map<point, char> points = configurations[tag];
			map<point, char>::iterator iter = points.begin();
			map<point, char>::iterator endIter = points.end();
			for (; iter != endIter; iter++) {
				lowest_energy_points.push_back(iter->first);
			}
		}
		return;
	}
	int k_free = LegalActions(p_before, configurations[tag]);
	if (k_free == 0) {
		return;
	}
	//各个动作好度
	vector<double>good_degrees;
	//计算各个动作的好度与权重预测值
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n], configurations[tag], good_degrees);
	//计算上下门限
	double upper_threshold = CalculateUpperThreshold(n);
	double lower_threshold = CalculateLowerThreshold(upper_threshold);
	//根据预测值与上下门限的数值关系分类讨论
	if (predict_wigtht >= lower_threshold && predict_wigtht <= upper_threshold) {
		//根据好度概率选择下一动作
		point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
		//能量增益
		int energy_increase = 0;
		//计算做完该动作的权重
		double present_weight = CalculateWeight(weight, next_action, input[n], configurations[tag], energy_increase);
		//更新
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n], energy_increase);
		//进入分支
		CalculationProcess(n+1, whole_length, tag, next_action, present_weight, input);
	}
	else if (predict_wigtht < lower_threshold) {
		//按照1/2的概率丢弃该分支
		double rand_result = random(0, 1);
		if (rand_result < 0.5) {
			return;
		}
		else {
			//根据好度概率选择下一动作
			point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
			//能量增益
			int energy_increase = 0;
			//计算做完该动作的权重
			double present_weight = CalculateWeight(weight, next_action, input[n], configurations[tag], energy_increase);
			//更新
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n], energy_increase);
			//进入分支
			CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
		}
	}
	else {
		
	}
}


int main()
{
	srand((int)time(0));
	string input = "PHPHPHPHPHPHPH";

	return 0;
}

