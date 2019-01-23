#include "stdafx.h"
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
		point _p = configurations_point[i];
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
		point p = configurations_point[i];
		char type = configurations_class[i];
		for (size_t j = i + 2; j < length; j++) {
			point _p = configurations_point[j];
			char _type = configurations_class[j];
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
		cout << present_energy << endl;
		cout << present_energy << endl;
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
	//分别计算出各可行动作在perm算法下的值，并且取能量最低的
	for (size_t i = 0; i < k_free; i++) {
		perm _perm;
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
		//设置perm参数
		_perm.SetAverageWeight(temp_average_weights);
		_perm.SetEnergy(present_energy + energy_increase);
		_perm.SetPointPosition(temp_configurations_point);
		_perm.SetPoint(temp_configurations_class);
		_perm.SetThisWeightNumber(temp_weights_numbers);
		_perm.CalculationProcess(n + 1, whole_length, 0, legal_actions[i], present_weight, input);
		if (_perm.GetEnergy() < min_energy) {
			best_index = i;
			_energy_increase = energy_increase;
			_present_weight = present_weight;
			min_energy = _perm.GetEnergy();
			if (_perm.GetEnergy() < perm_lowest_energy) {
				perm_lowest_energy = _perm.GetEnergy();
				_perm.GetPoint(perm_lowest_configurations_class);
				_perm.GetPointPosition(perm_lowest_configurations_point);
			}			
			//min_energy = _perm.GetEnergy();
			//_perm.GetPointPosition(best_point);
			//_perm.GetPoint(best_type);
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
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
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

//更新临时参数
void three_perm::UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//赋值
	ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point, perm::max_size_of_input);
	ArrayAssignment(_configurations_class, configurations_class, perm::max_size_of_input);
	//更新
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	++_weights_numbers[n - 1];
	_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}