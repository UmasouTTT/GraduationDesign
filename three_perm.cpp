#include "three_perm.h"


three_perm::three_perm()
{
}


three_perm::~three_perm()
{
}

//判断该坐标是否已经被使用
bool three_perm::IsThisPositionAlreadyOccupied(point p, int n) {
	for (size_t i = 0; i < n - 1; i++) {
		point _p = configurations_point_three[i];
		if (p.x == _p.x && p.y == _p.y) {
			return true;
		}
	}
	return false;
}

//计算合法的动作数
int three_perm::LegalActions(point p, int n) {
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
int three_perm::LegalActions(point p, vector<point> &legal_actions, int n) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
		legal_actions.push_back(p1);
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
		legal_actions.push_back(p2);
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
		legal_actions.push_back(p3);
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
		legal_actions.push_back(p4);
	}
	return result;
}


//测试运算结果是否正确
bool three_perm::TestResultIsSatisfied(int target_energy, int length) {
	int result = 0;
	for (size_t i = 0; i < length; i++) {
		point p = configurations_point_three[i];
		char type = configurations_class_three[i];
		for (size_t j = i + 2; j < length; j++) {
			point _p = configurations_point_three[j];
			char _type = configurations_class_three[j];
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

//计算两个点之间的距离
float three_perm::DistenceBetweenPoints(point point1, point point2) {
	float result = (float)(point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
	return sqrtf(result);
}


//初始化（初始化变元，前两个值为定值）
void three_perm::InitConfig(string &input, point &p, double &weight) {
	//清空数据
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
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
	configurations_point_three[0] = p1;
	configurations_class_three[0] = input[0];
	p1.x = p1.x + 1;
	p = p1;
	configurations_point_three[1] = p1;
	configurations_class_three[1] = input[1];
	//各分支当前构型能量
	present_energy = 0;
}

//算法
void three_perm::StartCalculate(string input) {
	//初始化
	perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	CircleCalculate(3, input.size(), p_second, start_weigtht, input);
	if (TestResultIsSatisfied(present_energy, input.length())) {
		cout << "test satisfied!" << endl;
	}
	else {
		cout << "error TAT" << endl;
	}
}

//迭代过程
void three_perm::CircleCalculate(int n, int whole_length, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		cout << "circle end, present energy :";
		cout << present_energy << endl;
		cout << "circle end, history energy :";
		cout << perm_lowest_energy << endl;
		return;
	}
	//获取当前状态可行的动作
	vector<point>legal_actions;
	int k_free = LegalActions(p_before, legal_actions, n);
	if (k_free == 0) {
		//理论上不会发生
		return;
	}
	int min_energy = 1;
	point best_point[perm::max_size_of_input];
	char best_type[perm::max_size_of_input];
	int best_index = -1;
	int _energy_increase;
	double _present_weight;
	//考虑两条支线找到的最小值相等的情况
	int num_of_lowestConfig = 0;
	//分别计算出各可行动作在perm算法下的值，并且取能量最低的
	for (size_t i = 0; i < k_free; i++) {
		perm _perm(-30, false);
		//计算在选择该点后的情况
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//临时参数
		double temp_average_weights[perm::max_size_of_input];
		//长度为n的构型的数量(需要初始化)
		double temp_weights_numbers[perm::max_size_of_input];
		//各分支具体构型
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//更新
		UpdateTempVariables(temp_average_weights, temp_weights_numbers, temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		/*//设置perm参数
		_perm.SetAverageWeight(temp_average_weights);
		_perm.SetEnergy(present_energy + energy_increase);
		_perm.SetPointPosition(temp_configurations_point);
		_perm.SetPoint(temp_configurations_class);
		_perm.SetThisWeightNumber(temp_weights_numbers);*/
		//_perm.CalculationProcess(n + 1, whole_length, 0, legal_actions[i], present_weight, input);
		//根据进度选择迭代次数，初始时迭代次数多，越往后越少
		int circle_times = 5;
		/*if (n / whole_length < 0.3) {
			circle_times = 30;
		}
		else if (n / whole_length < 0.7) {
			circle_times = 20;
		}
		else {
			circle_times = 10;
		}*/
		_perm.CircleCalculateProcess(n + 1, whole_length, legal_actions[i], present_weight, input, circle_times, temp_average_weights, present_energy + energy_increase, temp_configurations_point, temp_configurations_class, temp_weights_numbers);
		//如果当前分支所找到的最小能量比之前的都小，则将该分支选为最优分支
		if (_perm.GetEnergy() < min_energy) {
			best_index = i;
			num_of_lowestConfig = _perm.GetNumOfLowestConfigurations();
			_energy_increase = energy_increase;
			_present_weight = present_weight;
			min_energy = _perm.GetEnergy();
			if (_perm.GetEnergy() < perm_lowest_energy) {
				perm_lowest_energy = _perm.GetEnergy();
				_perm.GetPoint(perm_lowest_configurations_class);
				_perm.GetPointPosition(perm_lowest_configurations_point);
			}						
		}
		//如果两分支所找到的最小能量一致，则将最小构型数量多的设置为最优分支，若数量一致，则随机选
		else if (_perm.GetEnergy() == min_energy) {
			int present_num_of_config = _perm.GetNumOfLowestConfigurations();
			if (present_num_of_config > num_of_lowestConfig) {
				best_index = i;
				num_of_lowestConfig = present_num_of_config;
				_energy_increase = energy_increase;
				_present_weight = present_weight;
				min_energy = _perm.GetEnergy();
			}
			else if (present_num_of_config == num_of_lowestConfig) {
				int _randnum = random(0, 10);
				if (_randnum < 5) {
					best_index = i;
					num_of_lowestConfig = present_num_of_config;
					_energy_increase = energy_increase;
					_present_weight = present_weight;
					min_energy = _perm.GetEnergy();
				}
			}
		}
	}	
	//根据选择的情况进行更新
	UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
	CircleCalculate(n + 1, whole_length, legal_actions[best_index], _present_weight, input);

}


void  three_perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//更新权重算术平均值及该种构型长度的数量
	UpdateAverageWeight(weight, n);
	//更新各分支具体构型
	configurations_point_three[n - 1] = p;
	configurations_class_three[n - 1] = type;
	//更新各分支当前构型能量
	present_energy += energy_increase;
}

//更新Cn,Zn
void three_perm::UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n - 1] * weights_numbers[n - 1];
	++weights_numbers[n - 1];
	average_weights[n - 1] = (average_weight_before + w) / weights_numbers[n - 1];
}


//计算权重
double three_perm::CalculateWeight(double w, int energy_increase) {
	double result = w * exp(-energy_increase / perm::T);
	return result;
}


//计算能量增量
int three_perm::EnergyIncrease(point p, char type, point p_before, int n) {
	int result = 0;
	if (type == 'P') {
		return 0;
	}
	//遍历所有节点，判断距离
	for (size_t i = 0; i < n - 1; i++) {
		point _point = configurations_point_three[i];
		//在链上相邻不影响能量
		if (_point.x == p_before.x && _point.y == p_before.y) {
			continue;
		}
		char c = configurations_class_three[i];
		if (c == 'H' && DistenceBetweenPoints(p, _point) == 1) {
			result -= 1;
		}
	}
	return result;
}

//更新临时参数
void three_perm::UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//赋值
	ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point_three, perm::max_size_of_input);
	ArrayAssignment(_configurations_class, configurations_class_three, perm::max_size_of_input);
	//更新
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	++_weights_numbers[n - 1];
	_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}




//改进算法α
void three_perm::StartCalculateImproveFirst(string input) {
	//初始化
	perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	//首先迭代perm算法，在开始阶段尽可能的进行展开
	perm _perm(-35);
	_perm.StartCalculate(input, 5);
	//获取当前调用perm算法获取的最优部分构型数量
	int num_of_best_config = _perm.GetNumOfLowestConfigAndCanUse();
	//获取已知最有构型
	point best_config_ever[100][perm::max_size_of_input];
	char best_config_ever_points[perm::max_size_of_input];
	_perm.GetPoint(best_config_ever_points);
	//所使用构型长度
	int length_of_part_config = _perm.GetTheLengthStartPuneing();
	_perm.GetCurrentOptimalLocalConfiguration(best_config_ever, num_of_best_config, length_of_part_config);
	//对于已知的最优部分构型，逐一调用改进算法进行优化
	for (size_t i = 0; i < num_of_best_config; i++){
		//计算当前构型能量
		present_energy = CalculatePresentConfigEnergy(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//获取当前构型平均权重
		//double _weigtht = _perm.GetTargetLengthWeight(length_of_part_config);
		start_weigtht = 1;
		//获取当前点坐标
		p_second = best_config_ever[i][length_of_part_config - 1];
		//读取当前构型
		ArrayAssignment(configurations_point_three, best_config_ever[i], length_of_part_config);
		ArrayAssignment(configurations_class_three, best_config_ever_points, length_of_part_config);
		CircleCalculate(length_of_part_config + 1, input.size(), p_second, start_weigtht, input);
		if (TestResultIsSatisfied(present_energy, input.length())) {
			cout << "test satisfied!" << endl;
		}
		else {
			cout << "error TAT" << endl;
		}
	}	
}

//计算当前构型能量值
int three_perm::CalculatePresentConfigEnergy(point _configurations_point[], char _configurations_class[], int length) {
	int result = 0;
	for (size_t i = 0; i < length; i++) {
		point p = _configurations_point[i];
		char type = _configurations_class[i];
		for (size_t j = i + 2; j < length; j++) {
			point _p = _configurations_point[j];
			char _type = _configurations_class[j];
			if (type == _type && type == 'H') {
				float _temp = (float)(p.x - _p.x) * (p.x - _p.x) + (p.y - _p.y) * (p.y - _p.y);
				float _result = sqrtf(_temp);
				if (_result == 1) {
					result -= 1;
				}
			}
		}
	}
	return result;
}