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
#include<iostream>

#define random(a,b) (((double)rand()/RAND_MAX)*(b - a) + a)
#define judge_is_zero 0.000001

using namespace std;

struct point {
	int x;
	int y;
	bool operator < (const point &p) const {
		return x < p.x || (x == p.x && y < p.y);
	}
};

const double T = 0.25;
const double C0 = 10000;
const int Z0 = 1;
const int C = 1;
//当前最大分支标识号
int max_tag = 0;
//权重算术平均值(需要初始化)
vector<double>average_weights;
//长度为n的构型的数量(需要初始化)
vector<double>weights_numbers;
//各分支具体构型
vector<map<point, char>>configurations;
//各分支当前构型能量
vector<int>present_energy;
//最低能量
int lowest_energy = 0;
//最低能量构型
map<point, char>lowest_energy_points;
//用于获取合法组合集合
vector<int> input_numbers;
vector<int> combination_one;
vector<vector<int>> combination_result;

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
int EnergyIncrease(point p, char type, const map<point, char> &points, point p_before) {
	int result = 0;
	if (type == 'P') {
		return 0;
	}
	//遍历所有节点，判断距离
	map<point, char>::const_iterator iter = points.begin();
	map<point, char>::const_iterator endIter = points.end();
	for (; iter != endIter; iter++) {
		point _point = iter->first;
		//在链上相邻不影响能量
		if (_point.x == p_before.x && _point.y == p_before.y) {
			continue;
		}
		char c = iter->second;
		if (c == 'H' && DistenceBetweenPoints(p, _point) == 1) {
			result -= 1;
		}
	}
	return result;
}

//判断该坐标是否已经被使用
bool IsThisPositionAlreadyOccupied(point p, const map<point, char> &points) {
	map<point, char>::const_iterator iter = points.begin();
	map<point, char>::const_iterator endIter = points.end();
	for (; iter != endIter; iter++) {
		if (iter->first.x == p.x && iter->first.y == p.y) {
			return true;
		}
	}
	return false;
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
double CalculateGoodResults(point p, char type, const map<point, char> &points, point p_before) {
	double result = 0.0;
	int actions_later = LegalActions(p, points);
	result += ((double)actions_later + 0.5) * exp(-EnergyIncrease(p, type, points, p_before) / T);
	return result;	
}

//**************计算权重*************
double CalculateWeight(double w, point p, char type, const map<point, char> &points, int &energy_increase, point p_before) {
	energy_increase = EnergyIncrease(p, type, points, p_before);
	double result = w * exp(-energy_increase / T);
	return result;
}

//计算生长比例系数
double CalculatingLengthCoefficient(int n, int length) {
	if (n <= length * 0.3) {
		return 1;
	}
	if (n > length * 0.3 && n < length * 0.75) {
		return 30;
	}
	return 5;
}

//**************计算预计权重及各个动作的好度（避免重复计算）*************(由于内容较多，分两步进行)
double CalculatePredictWeightMid(double w, point p_before, char type, const map<point, char> &points, vector<double> &good_degrees, int k_free) {
	double result = 0;
	int legal_action_numbers = 0;
	//n+1步为上端放置
	point p1(p_before);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, points)) {
		good_degrees.push_back(CalculateGoodResults(p1, type, points, p_before));
		result += EnergyIncrease(p1, type, points, p_before);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为右端放置
	point p2(p_before);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, points)) {
		good_degrees.push_back(CalculateGoodResults(p2, type, points, p_before));
		result += EnergyIncrease(p2, type, points, p_before);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为下端放置
	point p3(p_before);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, points)) {
		good_degrees.push_back(CalculateGoodResults(p3, type, points, p_before));
		result += EnergyIncrease(p3, type, points, p_before);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//n+1步为左端放置
	point p4(p_before);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, points)) {
		good_degrees.push_back(CalculateGoodResults(p4, type, points, p_before));
		result += EnergyIncrease(p4, type, points, p_before);
		++legal_action_numbers;
	}
	else {
		good_degrees.push_back(0);
	}
	//与前一步权重求积
	if (type == 'P') {
		return w * k_free;
	}
	double energy_increase_average = result / legal_action_numbers;
	result = w * exp(-energy_increase_average / T);
	return result;
}
double CalculatePredictWeight(double w, point p_before, char type, const map<point, char> &points, vector<double> &good_degrees, int n, int length, int k_free) {
	double temp_result = CalculatePredictWeightMid(w, p_before, type, points, good_degrees, k_free);
	temp_result *= CalculatingLengthCoefficient(n, length);
	return temp_result;
}

//*************更新Cn,Zn***************
void UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n - 1] * weights_numbers[n - 1];
	++weights_numbers[n - 1];
	average_weights[n - 1] = (average_weight_before + w) / weights_numbers[n - 1];
}

//***************计算上门限***********
double CalculateUpperThreshold(int n) {
	double result = C * (average_weights[n - 1] / Z0) * (weights_numbers[n - 1] / C0) * (weights_numbers[n - 1] / C0);
	return result;
}

//**************计算下门限***********
double CalculateLowerThreshold(double upper_threshold) {
	double result = 0.2 * upper_threshold;
	return result;
}

//*******************创建新的分支*****************************
void CreateNewBranch(const map<point, char> &config_before, int energy_before) {
	//分支构型
	map<point, char> before_construction = config_before;
	configurations.push_back(before_construction);
	//分支能量
	int before_energy = energy_before;
	present_energy.push_back(before_energy);
	//分支标识
	++max_tag;
}

//*****************根据选择的更新全局变量***************
pair<map<point, char>, int> UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//记录之前的能量和构型
	pair<map<point, char>, int> result_before;
	map<point, char>config_before = configurations[tag];
	int energy_before = present_energy[tag];
	result_before = make_pair(config_before, energy_before);
	//更新权重算术平均值及该种构型长度的数量
	UpdateAverageWeight(weight, n);
	//更新各分支具体构型
	configurations[tag].insert(make_pair(p, type));
	//更新各分支当前构型能量
	present_energy[tag] += energy_increase;
	return result_before;
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

//递归计算排列组合
void CalculationCombinations(int offset, int k) {
	if (k == 0) {
		combination_result.push_back(combination_one);
		return;
	}
	//每次递归结束后，要考虑i是不是i <= people.size() - k，如果没有继续i++，如果i大于这个，返回上一次递归
	for (rsize_t i = offset; i <= input_numbers.size() - k; ++i) {
		combination_one.push_back(input_numbers[i]);
		CalculationCombinations(i + 1, k - 1);
		combination_one.pop_back();//删除combination最后一个元素
	}
}

//获取可能组合数
vector<vector<int>>GetCombinations(vector<int> &legal_actions, int num) {
	input_numbers.clear();
	combination_one.clear();
	combination_result.clear();
	//初始化输入数据
	for (size_t i = 0; i < legal_actions.size(); i++){
		input_numbers.push_back(legal_actions[i]);
	}
	//迭代计算组合数
	CalculationCombinations(0, num);
	return combination_result;
}

//根据数值获取相应的动作
vector<point> GetActionsByNum(vector<int> &numbers, point p_before) {
	vector<point>result;
	for (size_t i = 0; i < numbers.size(); i++){
		point temp_p(p_before);
		if (numbers[i] == 0) {
			temp_p.y = temp_p.y + 1;
		}
		else if(numbers[i] == 1){
			temp_p.x = temp_p.x + 1;
		}
		else if (numbers[i] == 2) {
			temp_p.y = temp_p.y - 1;
		}
		else if (numbers[i] == 3) {
			temp_p.x = temp_p.x - 1;
		}
		result.push_back(temp_p);
	}
	return result;
}

//************按照好度概率随机选择动作集合*****************
vector<point> ChooseActionsGroupByGoodDegrees(int k, vector<double> &good_degrees, point p_before) {
	//合法动作集合
	vector<int>legal_actions;
	for (size_t i = 0; i < good_degrees.size(); i++){
		if (good_degrees[i] - 0.0 < judge_is_zero && good_degrees[i] - 0.0 > -judge_is_zero) {
			continue;
		}
		legal_actions.push_back(i);
	}
	//计算可行组合数
	vector<vector<int>>combination_actions = GetCombinations(legal_actions, k);
	//计算好度总和
	vector<pair<double, double>>combinations_sum_section;
	//当前区间下限
	double present_good_degrees_sum = 0;
	for (size_t i = 0; i < combination_actions.size(); i++){
		double temp_good_degree_sum = 0;
		//计算该种组合的好度和
		for (size_t j = 0; j < combination_actions[i].size(); j++) {
			temp_good_degree_sum += good_degrees[combination_actions[i][j]];
		}
		combinations_sum_section.push_back(make_pair(present_good_degrees_sum, present_good_degrees_sum + temp_good_degree_sum));
		//更新区间下限
		present_good_degrees_sum += temp_good_degree_sum;
	}
	//随机选取
	double random_result = random(0, present_good_degrees_sum);
	//找到选取的集合
	int choose_com;
	for (size_t i = 0; i < combinations_sum_section.size(); i++){
		if (random_result >= combinations_sum_section[i].first && random_result <= combinations_sum_section[i].second) {
			choose_com = i;
			break;
		}
	}
	vector<point> choose_actions = GetActionsByNum(combination_actions[choose_com], p_before);
	return choose_actions;
}

//迭代计算各分支情况
void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		if (present_energy[tag] <= lowest_energy) {
			if (present_energy[tag] < lowest_energy) {
				cout << present_energy[tag] << endl;
			}
			lowest_energy = present_energy[tag];			
			lowest_energy_points.clear();
			lowest_energy_points = configurations[tag];
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
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n - 1], configurations[tag], good_degrees, n, whole_length, k_free);
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
		double present_weight = CalculateWeight(weight, next_action, input[n - 1], configurations[tag], energy_increase, p_before);
		//更新
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase);
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
			double present_weight = CalculateWeight(weight, next_action, input[n - 1], configurations[tag], energy_increase, p_before);
			//更新
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase);
			//进入分支
			CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
		}
	}
	else {
		//计算分支数量
		int k = Min((double)k_free, (double)(predict_wigtht / upper_threshold));
		//根据好度概率选择下一动作集合
		vector<point>choose_actions = ChooseActionsGroupByGoodDegrees(k, good_degrees, p_before);
		//记录更新前的值
		map<point, char> config_before;
		int energy_before;
		//根据各动作生成新的分支
		for (size_t i = 0; i < choose_actions.size(); i++){
			point next_action = choose_actions[i];
			//能量增益
			int energy_increase = 0;			
			if (i == 0) {//无需新建分支
				//计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], configurations[tag], energy_increase, p_before);
				//更新
				pair<map<point, char>, int>result_before = UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase);
				config_before = result_before.first;
				energy_before = result_before.second;
				//进入分支
				CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
			}
			else {//新建分支
				//计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], config_before, energy_increase, p_before);
				//建立新分支
				CreateNewBranch(config_before, energy_before);
				//更新
				UpdateGlobalVariables(present_weight, n, next_action, max_tag, input[n - 1], energy_increase);
				//进入分支
				CalculationProcess(n + 1, whole_length, max_tag, next_action, present_weight, input);
			}
		}
	}
}

//初始化（初始化变元，前两个值为定值）
void InitConfig(string &input, point &p, double &weight) {
	//清空数据
	average_weights.clear();
	weights_numbers.clear();
	configurations.clear();
	present_energy.clear();
	max_tag = 0;
	weight = 1;
	//权重算术平均值(需要初始化)
	for (size_t i = 0; i < input.length(); i++){
		average_weights.push_back(1);
	}
	//长度为n的构型的数量(需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers.push_back(1);
	}
	//各分支具体构型
	map<point, char>first_config;
	point p1;
	p1.x = 0;
	p1.y = 0;
	first_config.insert(make_pair(p1, input[0]));
	p1.x = p1.x + 1;
	p = p1;
	first_config.insert(make_pair(p1, input[1]));
	configurations.push_back(first_config);
	//各分支当前构型能量
	present_energy.push_back(0);
}

int main()
{
	srand((int)time(0));
	string input_string = "PPHPPHHPPPPHHPPPPHHPPPPHH";
	string input_string1 = "PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPPHHHHHHPPHHHHHHPPPPPPPPPHPHHPHHHHHHHHHHHPPHHHPHHPHPPHPHHHPPPPPPHHH";
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	map<point, char> result_energy_low_config;

	/*time_t tt = time(NULL);
	tm* t = localtime(&tt);
	printf("%d-%02d-%02d %02d:%02d:%02d\n", t->tm_year + 1900,
		t->tm_mon + 1,
		t->tm_mday,
		t->tm_hour,
		t->tm_min,
		t->tm_sec);*/

	int tag_i = 0;
	while (tag_i < 5) {
		InitConfig(input_string, p_second, start_weigtht);
		CalculationProcess(3, input_string.length(), 0, p_second, start_weigtht, input_string);
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			result_energy_low_config = lowest_energy_points;
			cout << result_energy_low << endl;
		}
		++tag_i;
	}

	tag_i = 0;
	while (tag_i < 5) {
		InitConfig(input_string1, p_second, start_weigtht);
		CalculationProcess(3, input_string1.length(), 0, p_second, start_weigtht, input_string1);
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			result_energy_low_config = lowest_energy_points;
			cout << result_energy_low << endl;
		}
		++tag_i;
	}
	
	
	return 0;
}

