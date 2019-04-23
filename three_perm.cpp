#include "three_perm.h"
const int three_perm::predict_worest_energy;


three_perm::three_perm()
{
	lowest_energy_first = 0;
	average_weights = new double[perm::max_size_of_input];
	weights_numbers = new double[perm::max_size_of_input];
	configurations_point_three = new point[perm::max_size_of_input];
	configurations_class_three = new char[perm::max_size_of_input];
	lowest_configurations_point = new point[perm::max_size_of_input];
	lowest_configurations_class = new char[perm::max_size_of_input];
	perm_lowest_configurations_point = new point[perm::max_size_of_input];
	perm_lowest_configurations_class = new char[perm::max_size_of_input];
}


three_perm::~three_perm()
{
	delete average_weights;
	delete weights_numbers;
	delete configurations_point_three;
	delete configurations_class_three;
	delete lowest_configurations_point;
	delete lowest_configurations_class;
	delete perm_lowest_configurations_point;
	delete perm_lowest_configurations_class;
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
int three_perm::LegalActions(point p, point legal_actions[4], int n) {
	int result = 0;
	//n+1步为上端放置
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		legal_actions[result] = p1;
		result += 1;		
	}
	//n+1步为右端放置
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		legal_actions[result] = p2;
		result += 1;
	}
	//n+1步为下端放置
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		legal_actions[result] = p3;
		result += 1;
	}
	//n+1步为左端放置
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		legal_actions[result] = p4;
		result += 1;
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
	perm _perm(predict_worest_energy);
	_perm.StartCalculateOnlyWeight(input, 1);
	//获取当前平均能量及历史搜索构型数
	_perm.GetAverageWeight(average_weights);
	_perm.GetWeightNumber(weights_numbers);
	//获取初次迭代最优能量值
	lowest_energy_first = _perm.GetEnergy();
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
		return;
	}
	//获取当前状态可行的动作
	point legal_actions[4];
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
	//用于判断是否找到比初次迭代结果更好的分枝
	bool _if_find_better_than_ever = false;
	//用于记录该分支当前构型能量值
	int _energy_present = present_energy;
	//分别计算出各可行动作在perm算法下的值，并且取能量最低的
	for (size_t i = 0; i < k_free; i++) {
		perm _perm(-43, false);
		//计算在选择该点后的情况
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//各分支具体构型
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//更新
		UpdateTempVariables(temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		//根据进度选择迭代次数，初始时迭代次数多，越往后越少
		int circle_times = 3;
		_perm.CircleCalculateProcess(n + 1, whole_length, legal_actions[i], present_weight, input, circle_times, _energy_present + energy_increase, temp_configurations_point, temp_configurations_class, average_weights, weights_numbers);
		//记录下该次迭代之后的权重值
		_perm.GetAverageWeight(average_weights);
		_perm.GetWeightNumber(weights_numbers);
		//如果找到比初次迭代更好的结果，直接进入该分支
		if (_perm.GetEnergy() < lowest_energy_first) {
			_if_find_better_than_ever = true;
			if (_perm.GetEnergy() < min_energy) {
				min_energy = _perm.GetEnergy();
				if (_perm.GetEnergy() < perm_lowest_energy) {
					perm_lowest_energy = _perm.GetEnergy();
					_perm.GetPoint(perm_lowest_configurations_class);
					_perm.GetPointPosition(perm_lowest_configurations_point);
				}
			}
			//获取当前构型能量值
			present_energy = _energy_present;
			//根据选择的情况进行更新
			UpdateGlobalVariables(present_weight, n, legal_actions[i], 0, input[n - 1], energy_increase);
			//进入新的分支
			CircleCalculate(n + 1, whole_length, legal_actions[i], present_weight, input);
		}
		else if(!_if_find_better_than_ever){
			//如果当前分支所找到的最小能量比之前的都小，则将该分支选为最优分支
			if (_perm.GetEnergy() < min_energy) {
				best_index = i;
				//获取本次迭代最少构型数
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
			//如果两分支所找到的最小能量一致，按照这两种分支所得到的最优构型数按概率随机选择分支，若数量一致，则随机选
			else if (_perm.GetEnergy() == min_energy) {
				int present_num_of_config = _perm.GetNumOfLowestConfigurations();
				int _randnum = random(0, present_num_of_config + num_of_lowestConfig);
				if (_randnum <= present_num_of_config) {
					best_index = i;
					num_of_lowestConfig = present_num_of_config;
					_energy_increase = energy_increase;
					_present_weight = present_weight;
				}
			}
		}		
	}	
	//如果没有出现最佳构型，根据选择的情况进行更新
	if (!_if_find_better_than_ever) {
		//获取当前构型能量值
		present_energy = _energy_present;
		UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
		CircleCalculate(n + 1, whole_length, legal_actions[best_index], _present_weight, input);
	}	
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
	average_weights[n - 1] = (average_weights[n - 1] + w) / 2;
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
void three_perm::UpdateTempVariables(point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//赋值
	//ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	//ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point_three, n-1);
	ArrayAssignment(_configurations_class, configurations_class_three, n-1);
	//更新
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	//double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	//++_weights_numbers[n - 1];
	//_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}




//改进算法α
void three_perm::StartCalculateImproveFirst(string input) {
	//初始化
	//perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	//首先迭代perm算法，在开始阶段尽可能的进行展开
	perm _perm(predict_worest_energy);
	_perm.StartCalculate(input, 5);
	//获取当前平均能量及历史搜索构型数
	_perm.GetAverageWeight(average_weights);
	_perm.GetWeightNumber(weights_numbers);
	//获取初次迭代最优能量值
	lowest_energy_first = _perm.GetEnergy();
	//获取当前调用perm算法获取的最优部分构型数量
	int num_of_best_config = _perm.GetNumOfLowestConfigAndCanUse();
	//获取已知最有构型
	point best_config_ever[100][perm::max_size_of_input];
	char best_config_ever_points[perm::max_size_of_input];
	_perm.GetPoint(best_config_ever_points);
	ArrayAssignment(configurations_class_three, best_config_ever_points, perm::max_size_of_input);
	//所使用构型长度
	int length_of_part_config = _perm.GetTheLengthStartPuneing();
	_perm.GetCurrentOptimalLocalConfiguration(best_config_ever, num_of_best_config, length_of_part_config);
	//对于已知的最优部分构型，逐一调用改进算法进行优化
	for (size_t i = 0; i < num_of_best_config; i++){
		//计算当前构型能量
		present_energy = CalculatePresentConfigEnergy(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//获取当前构型平均权重

		//double _weigtht = _perm.GetTargetLengthWeight(length_of_part_config);
		start_weigtht = CalculateStartConfigWeight(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//获取当前点坐标
		p_second = best_config_ever[i][length_of_part_config - 1];
		//读取当前构型
		ArrayAssignment(configurations_point_three, best_config_ever[i], length_of_part_config);
		
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




//计算初期构型权重
double three_perm::CalculateStartConfigWeight(point _configurations_point[], char _configurations_class[], int length) {
	double result = 1;
	for (size_t i = 2; i < length; i++){
		//计算能量增量
		int _part_energy_increase = 0;
		if (_configurations_class[i] == 'H') {
			for (size_t j = 0; j < i - 1; j++){
				if (_configurations_class[j] == 'H' && DistenceBetweenPoints(_configurations_point[i], _configurations_point[j]) == 1) {
					_part_energy_increase -= 1;
				}
			}
		}
		result *= exp(-_part_energy_increase / perm::T);
	}
	return result;
}



void three_perm::Branch_choose_improve_1(string input, int time) {
	perm_lowest_energy = 0;
	for (size_t i = 0; i < time; i++){
		StartCalculateImproveFirst(input);
	}
}





//算法(权重)
void three_perm::StartCalculateByWeight(string input) {
	//初始化
	perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	CircleCalculateByWeight(3, input.size(), p_second, start_weigtht, input);
	if (TestResultIsSatisfied(present_energy, input.length())) {
		cout << "test satisfied!" << endl;
	}
	else {
		cout << "error TAT" << endl;
	}
}


//迭代过程（权重）
void three_perm::CircleCalculateByWeight(int n, int whole_length, point p_before, double weight, string input) {
	//结束条件判断
	if (n > whole_length) {
		cout << "circle end, present energy :";
		cout << present_energy << endl;
		cout << "circle end, history energy :";
		cout << perm_lowest_energy << endl;
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
		return;
	}
	//获取当前状态可行的动作
	point legal_actions[4];
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
	//用于判断是否找到比初次迭代结果更好的分枝
	bool _if_find_better_than_ever = false;
	//用于记录该分支当前构型能量值
	int _energy_present = present_energy;
	//当前最佳权重
	double best_weight = 0;
	//当前最低构型能量
	int _lowest_energy = 0;
	//记录最佳平均权重值
	double best_average_weight[perm::max_size_of_input];
	//分别计算出各可行动作在perm算法下的值，并且取能量最低的
	cout << "各分支最低能量和最长链长权重" << endl;
	for (size_t i = 0; i < k_free; i++) {
		perm _perm;
		//计算在选择该点后的情况
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//各分支具体构型
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//更新
		UpdateTempVariables(temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		//根据进度选择迭代次数，初始时迭代次数多，越往后越少
		int circle_times = 1;
		_perm.CircleCalculateProcessOnlyConsiderWeight(n + 1, whole_length, legal_actions[i], present_weight, input, circle_times, _energy_present + energy_increase, temp_configurations_point, temp_configurations_class, average_weights, weights_numbers);
		cout << _perm.GetEnergy();
		cout << "     ";
		cout << _perm.GetAverageWeightWithLengthN(whole_length - 1);
		cout<<""<<endl;
		if (i == 0) {
			//记录最佳权重
			best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
			//记录当前最佳分支
			best_index = i;
			//记录最低能量
			_lowest_energy = _perm.GetEnergy();
			//获取均权重
			_perm.GetAverageWeight(best_average_weight);
			//判断是否找到更低能量
			if (_perm.GetEnergy() < min_energy) {
				min_energy = _perm.GetEnergy();
				if (_perm.GetEnergy() < perm_lowest_energy) {
					perm_lowest_energy = _perm.GetEnergy();
					_perm.GetPoint(perm_lowest_configurations_class);
					_perm.GetPointPosition(perm_lowest_configurations_point);
				}
			}
			//获取当前分支能量增幅和权重
			_energy_increase = energy_increase;
			_present_weight = present_weight;
		}
		//首先，对于比已知最低能量高超过1的分支不考虑
		else if (_perm.GetEnergy() - _lowest_energy < 1) {
			//如果当前出现的最低能量比之前分支最低能量低超过1，直接选择
			if (_lowest_energy - _perm.GetEnergy() > 1) {
				//记录最佳权重
				best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
				//记录当前最佳分支
				best_index = i;
				//记录最低能量
				_lowest_energy = _perm.GetEnergy();
				//获取均权重
				_perm.GetAverageWeight(best_average_weight);
				//判断是否找到更低能量
				if (_perm.GetEnergy() < min_energy) {
					min_energy = _perm.GetEnergy();
					if (_perm.GetEnergy() < perm_lowest_energy) {
						perm_lowest_energy = _perm.GetEnergy();
						_perm.GetPoint(perm_lowest_configurations_class);
						_perm.GetPointPosition(perm_lowest_configurations_point);
					}
				}
				//获取当前分支能量增幅和权重
				_energy_increase = energy_increase;
				_present_weight = present_weight;
			}
			else {
				//然后，按照比例进行选取
				int first_log = log(best_weight);
				int second_log = log(_perm.GetAverageWeightWithLengthN(whole_length - 1));
				//由于权重取log值，所以对于高权重分支进行补偿
				if (first_log > second_log) {
					first_log = (first_log - second_log) * second_log * n;
				}
				else {
					second_log = (second_log - first_log) * first_log * n;
				}
				int random_num = random(0, first_log + second_log);
				if (random_num > first_log) {
					//记录最佳权重
					best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
					//记录当前最佳分支
					best_index = i;
					//记录最低能量
					_lowest_energy = _perm.GetEnergy();
					//获取均权重
					_perm.GetAverageWeight(best_average_weight);
					//判断是否找到更低能量
					if (_perm.GetEnergy() < min_energy) {
						min_energy = _perm.GetEnergy();
						if (_perm.GetEnergy() < perm_lowest_energy) {
							perm_lowest_energy = _perm.GetEnergy();
							_perm.GetPoint(perm_lowest_configurations_class);
							_perm.GetPointPosition(perm_lowest_configurations_point);
						}
					}
					//获取当前分支能量增幅和权重
					_energy_increase = energy_increase;
					_present_weight = present_weight;
				}
			}			
		}
	}
	//进入分支
	cout << "选择分支序号:";
	cout << best_index << endl;
	ArrayAssignment(average_weights, best_average_weight, whole_length);
	UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
	CircleCalculateByWeight(n + 1, whole_length, legal_actions[best_index], _present_weight, input);
}

//循环改进算法2
void three_perm::CircleAlgripham2(string input, int time) {
	perm_lowest_energy = 0;
	for (size_t i = 0; i < time; i++) {
		StartCalculateByWeight(input);
	}
}
