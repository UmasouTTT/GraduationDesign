#include "stdafx.h"
#include "perm.h"


const int perm::max_size_of_input;
const int perm::max_size_of_legal_input;
const int perm::max_size_of_possibleConditions;

const double perm::T = 0.25;
const double perm::C0 = 1000;
const int perm::Z0 = 1;
const int perm::C = 1;

perm::perm()
{
}


perm::~perm()
{
}


//获取最低能量构型点坐标
void perm::GetPointPosition(point points[max_size_of_input]) {
	ArrayAssignment(points, lowest_configurations_point, max_size_of_input);
}
//获取最低能量构型点类型
void perm::GetPoint(char points[max_size_of_input]) {
	ArrayAssignment(points, lowest_configurations_class, max_size_of_input);
}

//设置最低能量构型点坐标
void perm::SetPointPosition(point points[max_size_of_input]) {
	ArrayAssignment(configurations_point, points, max_size_of_input);
}
//设置最低能量构型点类型
void perm::SetPoint(char points[max_size_of_input]) {
	ArrayAssignment(configurations_class, points, max_size_of_input);
}

//设置权重算术平均值
void perm::SetAverageWeight(double _average_weights[max_size_of_input]) {
	ArrayAssignment(average_weights, _average_weights, max_size_of_input);
}
//设置长度为n的构型的数量
void perm::SetThisWeightNumber(double _weights_numbers[max_size_of_input]) {
	ArrayAssignment(weights_numbers, _weights_numbers, max_size_of_input);
}



//算法内容
//计算两个点之间的距离
float perm::DistenceBetweenPoints(point point1, point point2) {
	float result = (float)(point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
	return sqrtf(result);
}

//计算能量增量
int perm::EnergyIncrease(point p, char type, point p_before, int n) {
	int result = 0;
	if (type == 'P') {
		return 0;
	}
	//遍历所有节点，判断距离
	for (size_t i = 0; i < n - 1; i++) {
		point _point = configurations_point[i];
		//在链上相邻不影响能量
		if (_point.x == p_before.x && _point.y == p_before.y) {
			continue;
		}
		char c = configurations_class[i];
		if (c == 'H' && DistenceBetweenPoints(p, _point) == 1) {
			result -= 1;
		}
	}
	return result;
}

//判断该坐标是否已经被使用
bool perm::IsThisPositionAlreadyOccupied(point p, int n) {
	for (size_t i = 0; i < n - 1; i++) {
		point _p = configurations_point[i];
		if (p.x == _p.x && p.y == _p.y) {
			return true;
		}
	}
	return false;
}

//计算合法的动作数
int perm::LegalActions(point p, int n) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
	}
	return result;
}
//重构计算合法动作函数，提高计算速率
int perm::LegalActions(point p, vector<pair<int, point>> &legal_actions, int n) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p1));
	}
	else {
		legal_actions.push_back(make_pair(0, p1));
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p2));
	}
	else {
		legal_actions.push_back(make_pair(0, p2));
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p3));
	}
	else {
		legal_actions.push_back(make_pair(0, p3));
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p4));
	}
	else {
		legal_actions.push_back(make_pair(0, p4));
	}
	return result;
}


//**************计算好度*************
double perm::CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n) {
	double result = 0.0;
	int actions_later = LegalActions(p, n);
	result += ((double)actions_later + 0.5) * exp(-energy_increase / T);
	return result;
}

//**************计算权重*************
double perm::CalculateWeight(double w, point p, char type, int energy_increase, point p_before) {
	double result = w * exp(-energy_increase / T);
	return result;
}

//计算生长比例系数
double perm::CalculatingLengthCoefficient(int n, int length) {
	if (n <= length * 0.3) {
		return 1;
	}
	if (n > length * 0.3 && n < length * 0.75) {
		return random(30, 45);
	}
	return random(5, 10);
}

//**************计算预计权重及各个动作的好度（避免重复计算）*************(由于内容较多，分两步进行)
double perm::CalculatePredictWeightMid(double w, point p_before, char type, vector<double> &good_degrees, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase, int n) {
	double result = 0;
	int legal_action_numbers = 0;
	for (size_t i = 0; i < legal_actions.size(); i++) {
		point p = legal_actions[i].second;
		if (legal_actions[i].first == 0) {
			good_degrees.push_back(0);
		}
		else {
			int e_increase = EnergyIncrease(p, type, p_before, n);
			good_degrees.push_back(CalculateGoodResults(p, type, p_before, e_increase, n));
			result += e_increase;
			energy_increase.insert(make_pair(p, e_increase));
			++legal_action_numbers;
		}
	}
	//与前一步权重求积
	if (type == 'P') {
		return w * k_free;
	}
	double energy_increase_average = result / legal_action_numbers;
	result = w * exp(-energy_increase_average / T);
	return result;
}
double perm::CalculatePredictWeight(double w, point p_before, char type, vector<double> &good_degrees, int n, int length, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase) {
	double temp_result = CalculatePredictWeightMid(w, p_before, type, good_degrees, k_free, legal_actions, energy_increase, n);
	temp_result *= CalculatingLengthCoefficient(n, length);
	return temp_result;
}

//*************更新Cn,Zn***************
void perm::UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n - 1] * weights_numbers[n - 1];
	++weights_numbers[n - 1];
	average_weights[n - 1] = (average_weight_before + w) / weights_numbers[n - 1];
}

//***************计算上门限***********
double perm::CalculateUpperThreshold(int n) {
	double result = C * (average_weights[n - 1] / Z0) * (weights_numbers[n - 1] / C0) * (weights_numbers[n - 1] / C0);
	return result;
}

//**************计算下门限***********
double perm::CalculateLowerThreshold(double upper_threshold) {
	double result = 0.2 * upper_threshold;
	return result;
}

//*******************创建新的分支*****************************
/*void CreateNewBranch(const map<point, char> &config_before, int energy_before) {
//分支构型
configurations.push_back(config_before);
//分支能量
//present_energy = (int *)realloc(present_energy, (max_tag + 2) * sizeof(int));
//present_energy[max_tag + 1] = energy_before;
//分支标识
++max_tag;
}*/

//*****************根据选择的更新全局变量***************
int  perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, point point_before[], char type_before[]) {
	//记录之前的能量和构型
	int energy_before = present_energy;
	ArrayAssignment(point_before, configurations_point, n);
	ArrayAssignment(type_before, configurations_class, n);
	//更新权重算术平均值及该种构型长度的数量
	UpdateAverageWeight(weight, n);
	//更新各分支具体构型
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
	//更新各分支当前构型能量
	present_energy += energy_increase;
	return energy_before;
}

int  perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//记录之前的能量和构型
	int energy_before = present_energy;
	//更新权重算术平均值及该种构型长度的数量
	UpdateAverageWeight(weight, n);
	//更新各分支具体构型
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
	//更新各分支当前构型能量
	present_energy += energy_increase;
	return energy_before;
}

//**********************按照概率生成随机动作********************************
point perm::GetNextActionByGoodDegrees(point p_before, vector<double> &good_degrees) {
	double whole_good_degrees = 0;
	double present_goodD_sum = good_degrees[0];
	for (size_t i = 0; i < good_degrees.size(); i++) {
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
void perm::CalculationCombinations(int offset, int k) {
	if (k == 0) {
		combination_result.push_back(combination_one);
		return;
	}
	//每次递归结束后，要考虑i是不是i <= people.size() - k，如果没有继续i++，如果i大于这个，返回上一次递归
	for (size_t i = offset; i <= input_numbers.size() - k; ++i) {
		combination_one.push_back(input_numbers[i]);
		CalculationCombinations(i + 1, k - 1);
		combination_one.pop_back();//删除combination最后一个元素
	}
}

//获取可能组合数
vector<vector<int>> perm::GetCombinations(vector<int> &legal_actions, int num) {
	input_numbers.clear();
	combination_one.clear();
	combination_result.clear();
	//初始化输入数据
	for (size_t i = 0; i < legal_actions.size(); i++) {
		input_numbers.push_back(legal_actions[i]);
	}
	//迭代计算组合数
	CalculationCombinations(0, num);
	return combination_result;
}
//根据数值获取相应的动作
vector<point> perm::GetActionsByNum(vector<int> &numbers, point p_before) {
	vector<point>result;
	for (size_t i = 0; i < numbers.size(); i++) {
		point temp_p(p_before);
		if (numbers[i] == 0) {
			temp_p.y = temp_p.y + 1;
		}
		else if (numbers[i] == 1) {
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
vector<point> perm::ChooseActionsGroupByGoodDegrees(int k, vector<double> &good_degrees, point p_before) {
	//合法动作集合
	vector<int>legal_actions;
	for (size_t i = 0; i < good_degrees.size(); i++) {
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
	for (size_t i = 0; i < combination_actions.size(); i++) {
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
	for (size_t i = 0; i < combinations_sum_section.size(); i++) {
		if (random_result >= combinations_sum_section[i].first && random_result <= combinations_sum_section[i].second) {
			choose_com = i;
			break;
		}
	}
	vector<point> choose_actions = GetActionsByNum(combination_actions[choose_com], p_before);
	return choose_actions;
}
//测试运算结果是否正确
bool perm::TestResultIsSatisfied(int target_energy, int length) {
	int result = 0;
	for (size_t i = 0; i < length; i++) {
		point p = lowest_configurations_point[i];
		char type = lowest_configurations_class[i];
		for (size_t j = i + 2; j < length; j++) {
			point _p = lowest_configurations_point[j];
			char _type = lowest_configurations_class[j];
			if (type == _type && type == 'H') {
				float _result = DistenceBetweenPoints(p, _p);
				if (DistenceBetweenPoints(p, _p) == 1) {
					result -= 1;
				}
			}
		}
	}
	if (result == target_energy) {
		return true;
	}
	return false;
}

//迭代计算各分支情况
void perm::CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		if (present_energy <= lowest_energy) {
			cout << present_energy << endl;
			lowest_energy = present_energy;
			ArrayAssignment(lowest_configurations_point, configurations_point, max_size_of_input);
			ArrayAssignment(lowest_configurations_class, configurations_class, max_size_of_input);
		}
		return;
	}
	vector<pair<int, point>>legal_actions;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	int k_free = LegalActions(p_before, legal_actions, n);
	if (k_free == 0) {
		return;
	}
	//各个动作好度
	vector<double>good_degrees;
	//各个动作导致的能量增益
	map<point, int>energy_increase;
	//计算各个动作的好度与权重预测值
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n - 1], good_degrees, n, whole_length, k_free, legal_actions, energy_increase);
	//计算上下门限
	double upper_threshold = CalculateUpperThreshold(n);
	double lower_threshold = CalculateLowerThreshold(upper_threshold);
	if (upper_threshold <= lower_threshold) {
		int not_ok = 1;
	}
	//根据预测值与上下门限的数值关系分类讨论
	if (predict_wigtht >= lower_threshold && predict_wigtht <= upper_threshold) {
		//根据好度概率选择下一动作
		point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
		//计算做完该动作的权重
		double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
		//更新
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action]);
		//进入分支
		CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
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
			//计算做完该动作的权重
			double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
			//更新
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action]);
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
		int energy_before;
		//根据各动作生成新的分支
		for (size_t i = 0; i < choose_actions.size(); i++) {
			point next_action = choose_actions[i];
			if (i == 0) {//无需新建分支
						 //计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
				//更新
				energy_before = UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action], point_before, type_before);
				//进入分支
				CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
			}
			else {//新建分支
				  //计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
				//建立新分支
				present_energy = energy_before;
				ArrayAssignment(configurations_point, point_before, n);
				ArrayAssignment(configurations_class, type_before, n);
				//更新
				UpdateGlobalVariables(present_weight, n, next_action, max_tag, input[n - 1], energy_increase[next_action]);
				//进入分支
				CalculationProcess(n + 1, whole_length, max_tag, next_action, present_weight, input);
			}
		}
	}
}

//初始化（初始化变元，前两个值为定值）
void perm::InitConfig(string &input, point &p, double &weight) {
	//清空数据
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
	max_tag = 0;
	weight = 1;
	//权重算术平均值(需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		average_weights[i] = 1;
	}
	//长度为n的构型的数量(需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 1;
	}
	//各分支具体构型
	point p1;
	p1.x = 0;
	p1.y = 0;
	configurations_point[0] = p1;
	configurations_class[0] = input[0];
	p1.x = p1.x + 1;
	p = p1;
	configurations_point[1] = p1;
	configurations_class[1] = input[1];
	//各分支当前构型能量
	present_energy = 0;
}


void perm::CalculateMaxSize(int length) {
	double result = 1;
	for (size_t i = 0; i < length; i++) {
		result *= 3;
	}
	int max_size_of_input = length;
}

void perm::InitGlobalVariable(string input) {
	//max_size_of_input = input.length();
	double average_weights[max_size_of_input];
	double weights_numbers[max_size_of_input];
}



void perm::StartCalculate(string input, int num_of_circle) {
	//InitGlobalVariable(input);
	//初始化变元
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	//获取时间
	/*struct tm t;   //tm结构指针
	time_t now;  //声明time_t类型变量
	time(&now);      //获取系统日期和时间
	char ch1[64] = { 0 };
	strftime(ch1, sizeof(ch1) - 1, "%Y-%m-%d %H:%M:%S", localtime_s(&t, &now));*/
	//cout << ch1 << endl;
	while (tag_i < num_of_circle) {
		InitConfig(input, p_second, start_weigtht);
		CalculationProcess(3, input.length(), 0, p_second, start_weigtht, input);
		if (TestResultIsSatisfied(lowest_energy, input.length())) {
			cout << "test satisfied!" << endl;
		}
		else {
			cout << "something wrongQAQ~" << endl;
		}
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			cout << "lowest energy: " << result_energy_low << "length of config : " << input.length() << endl;
		}
		/*t = time(NULL);
		char ch[64] = { 0 };
		strftime(ch, sizeof(ch) - 1, "%Y-%m-%d %H:%M:%S", localtime(&t));*/
		//cout << ch << endl;
		++tag_i;
	}
}
