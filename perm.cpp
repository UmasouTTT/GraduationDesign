#include "perm.h"


const int perm::max_size_of_input;
const int perm::max_size_of_legal_input;
const int perm::max_size_of_possibleConditions;

const double perm::T = 0.25;
const double perm::C0 = 10000;
const int perm::Z0 = 1;
const int perm::C = 1;

const float perm::part_config_for_save = 0.3;

const int perm::target_lowest_energy = -48;//目标最低构型

perm::perm()
{
	choose_config_length = 0;
	best_config_num = 0;
	best_config_ever = new point*[max_size_of_input];
	average_weights = new double[max_size_of_input];
	weights_numbers = new double[max_size_of_input];
	configurations_class = new char[max_size_of_input];
	configurations_point = new point[max_size_of_input];
	for (size_t i = 0; i < max_size_of_input; i++) {
		best_config_ever[i] = new point[max_size_of_input];
	}
	lowest_configurations_point = new point[max_size_of_input];
	lowest_configurations_class = new char[max_size_of_input];
	isPuneBegin = true;
	isSolutionNeedToBeSaved = true;
	beginPuningTheBranch = max_size_of_input;

	choose_actions_3 = new point[4];
	combinations_result = new int*[6];
	for (size_t i = 0; i < 6; i++) {
		combinations_result[i] = new int[4];
	}
}

perm::perm(int predict_energy) {
	choose_config_length = 0;
	best_config_num = 0;
	worest_energy = predict_energy;
	best_config_ever = new point*[max_size_of_input];
	average_weights = new double[max_size_of_input];
	weights_numbers = new double[max_size_of_input];
	configurations_class = new char[max_size_of_input];
	configurations_point = new point[max_size_of_input];
	for (size_t i = 0; i < max_size_of_input; i++){
		best_config_ever[i] = new point[max_size_of_input];
	}
	lowest_configurations_point = new point[max_size_of_input];
	lowest_configurations_class = new char[max_size_of_input];
	isPuneBegin = true;
	isSolutionNeedToBeSaved = true;
	beginPuningTheBranch = max_size_of_input;

	choose_actions_3 = new point[4];
	combinations_result = new int*[6];
	for (size_t i = 0; i < 6; i++) {
		combinations_result[i] = new int[4];
	}
}

perm::perm(int predict_energy, bool isSoultionSaved)
{
	worest_energy = predict_energy;
	choose_config_length = 0;
	best_config_num = 0;
	best_config_ever = new point*[max_size_of_input];
	average_weights = new double[max_size_of_input];
	weights_numbers = new double[max_size_of_input];
	configurations_class = new char[max_size_of_input];
	configurations_point = new point[max_size_of_input];
	for (size_t i = 0; i < max_size_of_input; i++) {
		best_config_ever[i] = new point[max_size_of_input];
	}
	lowest_configurations_point = new point[max_size_of_input];
	lowest_configurations_class = new char[max_size_of_input];
	isPuneBegin = true;
	isSolutionNeedToBeSaved = isSoultionSaved;
	beginPuningTheBranch = max_size_of_input;

	choose_actions_3 = new point[4];
	combinations_result = new int*[6];
	for (size_t i = 0; i < 6; i++) {
		combinations_result[i] = new int[4];
	}
}


perm::~perm()
{
	for (size_t i = 0; i < max_size_of_input; i++){
		delete[]best_config_ever[i];
	}
	delete best_config_ever;
	delete average_weights;
	delete weights_numbers;
	delete configurations_point;
	delete configurations_class;
	delete lowest_configurations_point;
	delete lowest_configurations_class;

	for (size_t i = 0; i < 6; i++) {
		delete[]combinations_result[i];
	}
	delete combinations_result;
	delete choose_actions_3;

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

//获取当前最优局部构型
void perm::GetCurrentOptimalLocalConfiguration(point _best_config_ever[100][perm::max_size_of_input], int num_of_best_config, int length) {
	for (size_t i = 0; i < num_of_best_config; i++){
		ArrayAssignment(_best_config_ever[i], best_config_ever[i], length);
	}
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
int perm::LegalActions(point p, point legal_actions_p[4], int legal_actions_t[4], int n) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
		legal_actions_p[0] = p1;
		legal_actions_t[0] = 1;
	}
	else {
		legal_actions_p[0] = p1;
		legal_actions_t[0] = 0;
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
		legal_actions_p[1] = p2;
		legal_actions_t[1] = 1;
	}
	else {
		legal_actions_p[1] = p2;
		legal_actions_t[1] = 0;
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
		legal_actions_p[2] = p3;
		legal_actions_t[2] = 1;
	}
	else {
		legal_actions_p[2] = p3;
		legal_actions_t[2] = 0;
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
		legal_actions_p[3] = p4;
		legal_actions_t[3] = 1;
	}
	else {
		legal_actions_p[3] = p4;
		legal_actions_t[3] = 0;
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
double perm::CalculatePredictWeightMid(double w, point p_before, char type, double good_degrees[4], int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4], int n) {
	double result = 0;
	int legal_action_numbers = 0;
	for (size_t i = 0; i < 4; i++) {
		point p = legal_actions_p[i];
		if (legal_actions_t[i] == 0) {
			good_degrees[i] = 0;
			energy_increase[i] = 0;
		}
		else {
			int e_increase = EnergyIncrease(p, type, p_before, n);
			good_degrees[i] = CalculateGoodResults(p, type, p_before, e_increase, n);
			result += e_increase;
			energy_increase[i] = e_increase;
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
double perm::CalculatePredictWeight(double w, point p_before, char type, double good_degrees[4], int n, int length, int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4]) {
	double temp_result = CalculatePredictWeightMid(w, p_before, type, good_degrees, k_free, legal_actions_p, legal_actions_t, energy_increase, n);
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
point perm::GetNextActionByGoodDegrees(point p_before, double good_degrees[4]) {
	double whole_good_degrees = 0;
	double present_goodD_sum = good_degrees[0];
	for (size_t i = 0; i < 4; i++) {
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
/*void perm::CalculationCombinations(int offset, int k) {
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
}*/

//获取可能组合数
void perm::GetCombinations(int(&legal_actions)[4], int num, int num_of_legal_actions, int &num_of_result) {
	num_of_result = 0;
	//迭代计算组合数
	if (num == 1) {
		for (size_t i = 0; i < num_of_legal_actions; i++) {
			combinations_result[i][0] = legal_actions[i];
			++num_of_result;
		}
	}
	else if (num == 2) {
		for (size_t i = 0; i < num_of_legal_actions; i++) {
			for (size_t j = i + 1; j < num_of_legal_actions; j++) {
				combinations_result[num_of_result][0] = legal_actions[i];
				combinations_result[num_of_result][1] = legal_actions[j];
				++num_of_result;
			}
		}
	}
	else if (num == 3) {
		for (size_t i = 0; i < num_of_legal_actions; i++) {
			for (size_t j = i + 1; j < num_of_legal_actions; j++) {
				for (size_t k = j + 1; k < num_of_legal_actions; k++) {
					combinations_result[num_of_result][0] = legal_actions[i];
					combinations_result[num_of_result][1] = legal_actions[j];
					combinations_result[num_of_result][2] = legal_actions[k];
					++num_of_result;
				}
			}
		}
	}
	else if (num == 4) {
		combinations_result[0][0] = legal_actions[0];
		combinations_result[0][1] = legal_actions[1];
		combinations_result[0][2] = legal_actions[2];
		combinations_result[0][3] = legal_actions[3];
		num_of_result = 1;
	}


	//CalculationCombinations(0, num);
	//return combination_result;
}
//根据数值获取相应的动作
void perm::GetActionsByNum(int index, point p_before, int length_of_result) {
	for (size_t i = 0; i < length_of_result; i++) {
		point temp_p(p_before);
		if (combinations_result[index][i] == 0) {
			temp_p.y = temp_p.y + 1;
		}
		else if (combinations_result[index][i] == 1) {
			temp_p.x = temp_p.x + 1;
		}
		else if (combinations_result[index][i] == 2) {
			temp_p.y = temp_p.y - 1;
		}
		else if (combinations_result[index][i] == 3) {
			temp_p.x = temp_p.x - 1;
		}
		choose_actions_3[i] = temp_p;
	}
}

//************按照好度概率随机选择动作集合*****************
void perm::ChooseActionsGroupByGoodDegrees(int k, double good_degrees[4], point p_before) {
	//合法动作集合
	int legal_actions[4];
	int j = 0;
	for (size_t i = 0; i < 4; i++) {
		if (good_degrees[i] - 0.0 < judge_is_zero && good_degrees[i] - 0.0 > -judge_is_zero) {
			continue;
		}
		legal_actions[j] = i;
		++j;
	}
	//计算可行组合数
	int num_of_result;
	GetCombinations(legal_actions, k, j, num_of_result);
	//计算好度总和
	double combinations_sum_section_1[6];
	double combinations_sum_section_2[6];
	//当前区间下限
	double present_good_degrees_sum = 0;
	for (size_t i = 0; i < num_of_result; i++) {
		double temp_good_degree_sum = 0;
		//计算该种组合的好度和
		for (size_t j = 0; j < k; j++) {
			temp_good_degree_sum += good_degrees[combinations_result[i][j]];
		}
		combinations_sum_section_1[i] = present_good_degrees_sum;
		combinations_sum_section_2[i] = present_good_degrees_sum + temp_good_degree_sum;
		//更新区间下限
		present_good_degrees_sum += temp_good_degree_sum;
	}
	//随机选取
	double random_result = random(0, present_good_degrees_sum);
	//找到选取的集合
	int choose_com;
	for (size_t i = 0; i < num_of_result; i++) {
		if (random_result >= combinations_sum_section_1[i] && random_result <= combinations_sum_section_2[i]) {
			choose_com = i;
			break;
		}
	}

	GetActionsByNum(choose_com, p_before, k);
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
		if (present_energy < lowest_energy) {
			num_of_lowestConfigurations = 1;
			cout << "find lower energy configuration, present energy is:";
			cout << present_energy << endl;
			lowest_energy = present_energy;
			best_config_num = 0;
			num_of_lowestConfigurations = 0;
			if (present_energy <= worest_energy) {
				AddNewConfigToBestConfigEver(configurations_point, max_size_of_input);
			}			
			ArrayAssignment(lowest_configurations_point, configurations_point, max_size_of_input);
			ArrayAssignment(lowest_configurations_class, configurations_class, max_size_of_input);
			if (present_energy == perm::target_lowest_energy) {
				cout << "find target config!" << endl;
				struct tm t;   //tm结构指针
				time_t now;  //声明time_t类型变量
				time(&now);      //获取系统日期和时间
				localtime_s(&t, &now);   //获取当地日期和时间
				string present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
				cout << present_time << endl;
				cout << "end" << endl;
			}
		}
		else if (present_energy == lowest_energy) {
			if (present_energy <= worest_energy) {
				AddNewConfigToBestConfigEver(configurations_point, max_size_of_input);
			}
			cout << "find new configuration :";
			cout << num_of_lowestConfigurations;
			cout << "  present energy is :";
			cout << present_energy << endl;
		}
		return;
	}
	//vector<pair<int, point>>legal_actions;
	point legal_actions_p[4];
	int legal_actions_t[4];
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	int k_free = LegalActions(p_before, legal_actions_p, legal_actions_t, n);
	if (k_free == 0) {
		return;
	}
	//各个动作好度
	double good_degrees[4];
	//各个动作导致的能量增益
	//map<point, int>energy_increase;
	int energy_increase[4];
	//计算各个动作的好度与权重预测值
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n - 1], good_degrees, n, whole_length, k_free, legal_actions_p, legal_actions_t, energy_increase);
	//计算上下门限
	double upper_threshold;
	double lower_threshold;
	//double upper_threshold = CalculateUpperThreshold(n);
	//double lower_threshold = CalculateLowerThreshold(upper_threshold);
	if (weights_numbers[n - 1] == 0) {
		upper_threshold = predict_wigtht + 1;
		lower_threshold = 0;
	}
	else {
		upper_threshold = CalculateUpperThreshold(n);
		lower_threshold = CalculateLowerThreshold(upper_threshold);
	}
	if (upper_threshold <= lower_threshold) {
		int not_ok = 1;
	}
	//根据预测值与上下门限的数值关系分类讨论
	if (predict_wigtht >= lower_threshold && predict_wigtht <= upper_threshold) {
		//获取剪枝开始链长
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//根据好度概率选择下一动作
		point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
		//获取动作标号
		int action_tag = GetActionsNum(legal_actions_p, next_action);
		//计算做完该动作的权重
		double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
		//更新
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag]);
		//进入分支
		CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
	}
	else if (predict_wigtht < lower_threshold) {
		//获取剪枝开始链长
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//按照1/2的概率丢弃该分支
		double rand_result = random(0, 1);
		if (rand_result < 0.5) {
			return;
		}
		else {
			//根据好度概率选择下一动作
			point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
			//获取动作标号
			int action_tag = GetActionsNum(legal_actions_p, next_action);
			//计算做完该动作的权重
			double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
			//更新
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag]);
			//进入分支
			CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
		}
	}
	else {
		//计算分支数量
		int k = Min((double)k_free, (double)(predict_wigtht / upper_threshold));
		//根据好度概率选择下一动作集合
		ChooseActionsGroupByGoodDegrees(k, good_degrees, p_before);
		//记录更新前的值
		int energy_before;
		//记录合法动作集
		point _choose_actions[4];
		for (size_t i = 0; i < k; i++){
			_choose_actions[i] = choose_actions_3[i];
		}
		//根据各动作生成新的分支
		for (size_t i = 0; i < k; i++) {
			point next_action = _choose_actions[i];
			if (i == 0) {//无需新建分支
						 //获取动作标号
				int action_tag = GetActionsNum(legal_actions_p, next_action);
						 //计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
				//更新
				energy_before = UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag], point_before, type_before);
				//进入分支
				CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
			}
			else {//新建分支
				  //获取动作标号
				int action_tag = GetActionsNum(legal_actions_p, next_action);
				  //计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
				//建立新分支
				present_energy = energy_before;
				ArrayAssignment(configurations_point, point_before, n);
				ArrayAssignment(configurations_class, type_before, n);
				//更新
				UpdateGlobalVariables(present_weight, n, next_action, max_tag, input[n - 1], energy_increase[action_tag]);
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
		average_weights[i] = 0;
	}
	//长度为n的构型的数量(需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 0;
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
	//初始化权重
	CalculationProcess(3, input.length(), 0, p, weight, input);
	if (TestResultIsSatisfied(lowest_energy, input.length())) {
		cout << "test satisfied!" << endl;
	}
	else {
		cout << "something wrongQAQ~" << endl;
	}
	//重置其他数据
	max_tag = 0;
	weight = 1;
	//各分支当前构型能量
	present_energy = 0;
	//若没有制定最低能量，则由第一次迭代结构决定
	if (worest_energy == 0) {
		worest_energy == lowest_energy;
	}
}


//初始化（初始化变元，前两个值为定值，无需生成初始权重）
void perm::InitConfigWithoutInitWeight(string &input, point &p, double &weight) {
	//清空数据
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
	max_tag = 0;
	weight = 1;
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
	choose_config_length = input.length() * part_config_for_save;
	//初始化变元
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	//获取最低能量的最多构型数
	num_of_lowestConfigurations = 0;
	int temp_lowestConfig = 0;
	//获取时间
	/*struct tm t;   //tm结构指针
	time_t now;  //声明time_t类型变量
	time(&now);      //获取系统日期和时间
	char ch1[64] = { 0 };
	strftime(ch1, sizeof(ch1) - 1, "%Y-%m-%d %H:%M:%S", localtime_s(&t, &now));*/
	//cout << ch1 << endl;
	while (tag_i < num_of_circle) {
		if (tag_i == 1) {
			isPuneBegin = true;
		}
		if (tag_i == 0) {
			InitConfig(input, p_second, start_weigtht);
			isPuneBegin = false;
		}
		else {
			InitConfigWithoutInitWeight(input, p_second, start_weigtht);
		}		
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
			//temp_lowestConfig = num_of_lowestConfigurations;
			//num_of_lowestConfigurations = 1;
		}
		/*t = time(NULL);
		char ch[64] = { 0 };
		strftime(ch, sizeof(ch) - 1, "%Y-%m-%d %H:%M:%S", localtime(&t));*/
		//cout << ch << endl;
		//根据拟人策略将权重数设置为1
		for (size_t i = 0; i < max_size_of_input; i++) {
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}


//多次计算各分支情况
void perm::CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]) {
	//初始化局部变量
	int tag_i = 0;
	int result_energy_low = 0;	
	num_of_lowestConfigurations = 0;
	int temp_lowestConfig = 0;
	//循环体
	while (tag_i < cirlce_times) {
		//设置perm参数
		if (tag_i == 0) {
			SetAverageWeight(_average_weight);
			SetThisWeightNumber(_weight_number);
		}
		//SetAverageWeight(_average_weights);
		SetEnergy(_energy);
		SetPointPosition(points);
		SetPoint(_points);
		//SetThisWeightNumber(_weights_numbers);
		CalculationProcess(n, input.length(), 0, p_before, weight, input);
		if (TestResultIsSatisfied(lowest_energy, input.length())) {
			cout << "test satisfied!" << endl;
		}
		else if(lowest_energy != 0){
			cout << "something wrongQAQ~" << endl;
		}
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			cout << "lowest energy: " << result_energy_low << "length of config : " << input.length() << endl;
			//temp_lowestConfig = num_of_lowestConfigurations;
			//num_of_lowestConfigurations = 1;
		}
		//选择最多能量构型数
	
		/*t = time(NULL);
		char ch[64] = { 0 };
		strftime(ch, sizeof(ch) - 1, "%Y-%m-%d %H:%M:%S", localtime(&t));*/
		//cout << ch << endl;

		//根据拟人策略将权重数设置为1
		for (size_t i = 0; i < max_size_of_input; i++){
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}

//初始化起始权重
bool perm::InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, int _energy, point points[max_size_of_input], char _points[max_size_of_input]) {
	max_tag = 0;
	weight = 1;
	//权重算术平均值(第一次迭代时获取初始值需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		average_weights[i] = 0;
	}
	//长度为n的构型的数量(需要初始化)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 0;
	}
	//各分支具体构型
	point p1;
	p1.x = 0;
	p1.y = 0;
	configurations_point[0] = p1;
	configurations_class[0] = input[0];
	p1.x = p1.x + 1;
	p_before = p1;
	configurations_point[1] = p1;
	configurations_class[1] = input[1];
	//各分支当前构型能量
	present_energy = 0;
	//初始化权重
	CalculationProcess(3, input.length(), 0, p_before, weight, input);
	if (TestResultIsSatisfied(lowest_energy, input.length())) {
		cout << "finish start weight data!" << endl;
	}
	else if(lowest_energy == 0){
		cout << "can not grow" << endl;
		return false;
	}
	else {
		cout << "something wrongQAQ~" << endl;
	}
	//重置其他数据
	max_tag = 0;
	return true;
}


//判断两种构型是否可以认为一致（目前认为两种构型保持在所选链长一致即认为两种构型一致）
bool perm::IsTwoConfigTheSame(point points[max_size_of_input], point _points[max_size_of_input], int start, int end) {
	for (size_t i = start; i < end; i++){
		if (points[i].x == _points[i].x && points[i].y == _points[i].y) {
			continue;
		}
		else {
			return false;
		}
	}
	return true;
}


//向历史最优构型中添加新构型
void perm::AddNewConfigToBestConfigEver(point points[max_size_of_input], int length) {
	//对于获取部分最优构型的perm算法执行条件一，对于迭代选择分支的Perm算法执行条件二
	if (isSolutionNeedToBeSaved) {
		if (best_config_num < 100) {
			for (size_t i = 0; i < best_config_num; i++) {
				if (IsTwoConfigTheSame(best_config_ever[i], points, 0, choose_config_length)) {
					return;
				}
			}
			ArrayAssignment(best_config_ever[best_config_num], points, choose_config_length);
			++best_config_num;
		}		
	}
	else if(num_of_lowestConfigurations < 100){
		for (size_t i = 0; i < num_of_lowestConfigurations; i++) {
			if (IsTwoConfigTheSame(best_config_ever[i], points, 0, length)) {
				return;
			}
		}
		ArrayAssignment(best_config_ever[num_of_lowestConfigurations], points, length);
		++num_of_lowestConfigurations;
	}
}


//获取动作序号
int perm::GetActionsNum(point actions[4], point action) {
	for (size_t i = 0; i < 4; i++){
		if ((actions[i].x == action.x) && (actions[i].y == action.y)) {
			return i;
		}
	}
	//错误
	return -1;
}

void perm::GetAverageWeight(double _average_weights[max_size_of_input]) {
	ArrayAssignment(_average_weights, average_weights, max_size_of_input);
}

void perm::GetWeightNumber(double _weights_numbers[max_size_of_input]) {
	ArrayAssignment(_weights_numbers, weights_numbers, max_size_of_input);
}




//迭代计算各分支情况(只考虑均权重)
void perm::CalculationProcessOnlyConsiderWeight(int n, int whole_length, int tag, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		if (present_energy < lowest_energy) {
			num_of_lowestConfigurations = 1;
			//cout << "find lower energy configuration, present energy is:";
			//cout << present_energy << endl;
			lowest_energy = present_energy;
			ArrayAssignment(lowest_configurations_point, configurations_point, max_size_of_input);
			ArrayAssignment(lowest_configurations_class, configurations_class, max_size_of_input);
			if (present_energy == perm::target_lowest_energy) {
				cout << "find target config!" << endl;
				struct tm t;   //tm结构指针
				time_t now;  //声明time_t类型变量
				time(&now);      //获取系统日期和时间
				localtime_s(&t, &now);   //获取当地日期和时间
				string present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
				cout << present_time << endl;
				cout << "end" << endl;
			}
		}
		else if (present_energy == lowest_energy) {
			/*cout << "find new configuration :";
			cout << num_of_lowestConfigurations;
			cout << "  present energy is :";
			cout << present_energy << endl;*/
		}
		return;
	}
	point legal_actions_p[4];
	int legal_actions_t[4];
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	int k_free = LegalActions(p_before, legal_actions_p, legal_actions_t, n);
	if (k_free == 0) {
		return;
	}
	//各个动作好度
	double good_degrees[4];
	//各个动作导致的能量增益
	//map<point, int>energy_increase;
	int energy_increase[4];
	//计算各个动作的好度与权重预测值
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n - 1], good_degrees, n, whole_length, k_free, legal_actions_p, legal_actions_t, energy_increase);
	//计算上下门限
	double upper_threshold;
	double lower_threshold;
	if (weights_numbers[n - 1] == 0) {
		upper_threshold = predict_wigtht + 1;
		lower_threshold = 0;
	}
	else {
		upper_threshold = CalculateUpperThreshold(n);
		lower_threshold = CalculateLowerThreshold(upper_threshold);
	}
	if (upper_threshold <= lower_threshold) {
		int not_ok = 1;
	}
	//根据预测值与上下门限的数值关系分类讨论
	if (predict_wigtht >= lower_threshold && predict_wigtht <= upper_threshold) {
		//获取剪枝开始链长
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//根据好度概率选择下一动作
		point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
		//获取动作标号
		int action_tag = GetActionsNum(legal_actions_p, next_action);
		//计算做完该动作的权重
		double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
		//更新
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag]);
		//进入分支
		CalculationProcessOnlyConsiderWeight(n + 1, whole_length, tag, next_action, present_weight, input);
	}
	else if (predict_wigtht < lower_threshold) {
		//获取剪枝开始链长
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//按照1/2的概率丢弃该分支
		double rand_result = random(0, 1);
		if (rand_result < 0.5) {
			return;
		}
		else {
			//根据好度概率选择下一动作
			point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
			//获取动作标号
			int action_tag = GetActionsNum(legal_actions_p, next_action);
			//计算做完该动作的权重
			double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
			//更新
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag]);
			//进入分支
			CalculationProcessOnlyConsiderWeight(n + 1, whole_length, tag, next_action, present_weight, input);
		}
	}
	else {
		//计算分支数量
		int k = Min((double)k_free, (double)(predict_wigtht / upper_threshold));
		//根据好度概率选择下一动作集合
		ChooseActionsGroupByGoodDegrees(k, good_degrees, p_before);
		//记录更新前的值
		int energy_before;
		//记录合法动作集
		point _choose_actions[4];
		for (size_t i = 0; i < k; i++) {
			_choose_actions[i] = choose_actions_3[i];
		}
		//根据各动作生成新的分支
		for (size_t i = 0; i < k; i++) {
			point next_action = _choose_actions[i];
			if (i == 0) {//无需新建分支
						 //获取动作标号
				int action_tag = GetActionsNum(legal_actions_p, next_action);
				//计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
				//更新
				energy_before = UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[action_tag], point_before, type_before);
				//进入分支
				CalculationProcessOnlyConsiderWeight(n + 1, whole_length, tag, next_action, present_weight, input);
			}
			else {//新建分支
				  //获取动作标号
				int action_tag = GetActionsNum(legal_actions_p, next_action);
				//计算做完该动作的权重
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[action_tag], p_before);
				//建立新分支
				present_energy = energy_before;
				ArrayAssignment(configurations_point, point_before, n);
				ArrayAssignment(configurations_class, type_before, n);
				//更新
				UpdateGlobalVariables(present_weight, n, next_action, max_tag, input[n - 1], energy_increase[action_tag]);
				//进入分支
				CalculationProcessOnlyConsiderWeight(n + 1, whole_length, max_tag, next_action, present_weight, input);
			}
		}
	}
}


//多次计算各分支情况
void perm::CircleCalculateProcessOnlyConsiderWeight(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]) {
	//初始化局部变量
	int tag_i = 0;
	int result_energy_low = 0;
	int temp_lowestConfig = 0;
	//循环体
	while (tag_i < cirlce_times) {
		//设置perm参数
		if (tag_i == 0) {
			SetAverageWeight(_average_weight);
			SetThisWeightNumber(_weight_number);
		}
		SetEnergy(_energy);
		SetPointPosition(points);
		SetPoint(_points);
		CalculationProcessOnlyConsiderWeight(n, input.length(), 0, p_before, weight, input);
		if (TestResultIsSatisfied(lowest_energy, input.length())) {
			//cout << "test satisfied!" << endl;
		}
		else if (lowest_energy != 0) {
			cout << "something wrongQAQ~" << endl;
		}
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			//cout << "lowest energy: " << result_energy_low << "length of config : " << input.length() << endl;
		}

		//根据拟人策略将权重数设置为1
		for (size_t i = 0; i < max_size_of_input; i++) {
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}




void perm::StartCalculateOnlyWeight(string input, int num_of_circle) {
	//初始化变元
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	while (tag_i < num_of_circle) {
		if (tag_i == 0) {
			InitConfig(input, p_second, start_weigtht);
		}
		else {
			InitConfigWithoutInitWeight(input, p_second, start_weigtht);
		}
		CalculationProcessOnlyConsiderWeight(3, input.length(), 0, p_second, start_weigtht, input);
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
		//根据拟人策略将权重数设置为1
		for (size_t i = 0; i < max_size_of_input; i++) {
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}