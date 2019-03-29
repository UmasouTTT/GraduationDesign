#include "three_perm.h"


three_perm::three_perm()
{
}


three_perm::~three_perm()
{
}

//�жϸ������Ƿ��Ѿ���ʹ��
bool three_perm::IsThisPositionAlreadyOccupied(point p, int n) {
	for (size_t i = 0; i < n - 1; i++) {
		point _p = configurations_point_three[i];
		if (p.x == _p.x && p.y == _p.y) {
			return true;
		}
	}
	return false;
}

//����Ϸ��Ķ�����
int three_perm::LegalActions(point p, int n) {
	int result = 0;
	//n+1��Ϊ�϶˷���
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
	}
	//n+1��Ϊ�Ҷ˷���
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
	}
	//n+1��Ϊ�¶˷���
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
	}
	//n+1��Ϊ��˷���
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
	}
	return result;
}
//�ع�����Ϸ�������������߼�������
int three_perm::LegalActions(point p, vector<point> &legal_actions, int n) {
	int result = 0;
	//n+1��Ϊ�϶˷���
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
		legal_actions.push_back(p1);
	}
	//n+1��Ϊ�Ҷ˷���
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
		legal_actions.push_back(p2);
	}
	//n+1��Ϊ�¶˷���
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
		legal_actions.push_back(p3);
	}
	//n+1��Ϊ��˷���
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		result += 1;
		legal_actions.push_back(p4);
	}
	return result;
}


//�����������Ƿ���ȷ
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

//����������֮��ľ���
float three_perm::DistenceBetweenPoints(point point1, point point2) {
	float result = (float)(point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
	return sqrtf(result);
}


//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ��
void three_perm::InitConfig(string &input, point &p, double &weight) {
	//�������
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
	weight = 1;
	//Ȩ������ƽ��ֵ(��Ҫ��ʼ��)
	for (size_t i = 0; i < input.length(); i++) {
		average_weights[i] = 0;
	}
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 0;
	}
	//����֧���幹��
	point p1;
	p1.x = 0;
	p1.y = 0;
	configurations_point_three[0] = p1;
	configurations_class_three[0] = input[0];
	p1.x = p1.x + 1;
	p = p1;
	configurations_point_three[1] = p1;
	configurations_class_three[1] = input[1];
	//����֧��ǰ��������
	present_energy = 0;
}

//�㷨
void three_perm::StartCalculate(string input) {
	//��ʼ��
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

//��������
void three_perm::CircleCalculate(int n, int whole_length, point p_before, double weight, string input) {
	//���������ж�
	if (n > whole_length) {
		cout << "circle end, present energy :";
		cout << present_energy << endl;
		cout << "circle end, history energy :";
		cout << perm_lowest_energy << endl;
		return;
	}
	//��ȡ��ǰ״̬���еĶ���
	vector<point>legal_actions;
	int k_free = LegalActions(p_before, legal_actions, n);
	if (k_free == 0) {
		//�����ϲ��ᷢ��
		return;
	}
	int min_energy = 1;
	point best_point[perm::max_size_of_input];
	char best_type[perm::max_size_of_input];
	int best_index = -1;
	int _energy_increase;
	double _present_weight;
	//��������֧���ҵ�����Сֵ��ȵ����
	int num_of_lowestConfig = 0;
	//�ֱ����������ж�����perm�㷨�µ�ֵ������ȡ������͵�
	for (size_t i = 0; i < k_free; i++) {
		perm _perm(-30, false);
		//������ѡ��õ������
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//��ʱ����
		double temp_average_weights[perm::max_size_of_input];
		//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
		double temp_weights_numbers[perm::max_size_of_input];
		//����֧���幹��
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//����
		UpdateTempVariables(temp_average_weights, temp_weights_numbers, temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		/*//����perm����
		_perm.SetAverageWeight(temp_average_weights);
		_perm.SetEnergy(present_energy + energy_increase);
		_perm.SetPointPosition(temp_configurations_point);
		_perm.SetPoint(temp_configurations_class);
		_perm.SetThisWeightNumber(temp_weights_numbers);*/
		//_perm.CalculationProcess(n + 1, whole_length, 0, legal_actions[i], present_weight, input);
		//���ݽ���ѡ�������������ʼʱ���������࣬Խ����Խ��
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
		//�����ǰ��֧���ҵ�����С������֮ǰ�Ķ�С���򽫸÷�֧ѡΪ���ŷ�֧
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
		//�������֧���ҵ�����С����һ�£�����С���������������Ϊ���ŷ�֧��������һ�£������ѡ
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
	//����ѡ���������и���
	UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
	CircleCalculate(n + 1, whole_length, legal_actions[best_index], _present_weight, input);

}


void  three_perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//����Ȩ������ƽ��ֵ�����ֹ��ͳ��ȵ�����
	UpdateAverageWeight(weight, n);
	//���¸���֧���幹��
	configurations_point_three[n - 1] = p;
	configurations_class_three[n - 1] = type;
	//���¸���֧��ǰ��������
	present_energy += energy_increase;
}

//����Cn,Zn
void three_perm::UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n - 1] * weights_numbers[n - 1];
	++weights_numbers[n - 1];
	average_weights[n - 1] = (average_weight_before + w) / weights_numbers[n - 1];
}


//����Ȩ��
double three_perm::CalculateWeight(double w, int energy_increase) {
	double result = w * exp(-energy_increase / perm::T);
	return result;
}


//������������
int three_perm::EnergyIncrease(point p, char type, point p_before, int n) {
	int result = 0;
	if (type == 'P') {
		return 0;
	}
	//�������нڵ㣬�жϾ���
	for (size_t i = 0; i < n - 1; i++) {
		point _point = configurations_point_three[i];
		//���������ڲ�Ӱ������
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

//������ʱ����
void three_perm::UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//��ֵ
	ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point_three, perm::max_size_of_input);
	ArrayAssignment(_configurations_class, configurations_class_three, perm::max_size_of_input);
	//����
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	++_weights_numbers[n - 1];
	_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}




//�Ľ��㷨��
void three_perm::StartCalculateImproveFirst(string input) {
	//��ʼ��
	perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	//���ȵ���perm�㷨���ڿ�ʼ�׶ξ����ܵĽ���չ��
	perm _perm(-35);
	_perm.StartCalculate(input, 5);
	//��ȡ��ǰ����perm�㷨��ȡ�����Ų��ֹ�������
	int num_of_best_config = _perm.GetNumOfLowestConfigAndCanUse();
	//��ȡ��֪���й���
	point best_config_ever[100][perm::max_size_of_input];
	char best_config_ever_points[perm::max_size_of_input];
	_perm.GetPoint(best_config_ever_points);
	//��ʹ�ù��ͳ���
	int length_of_part_config = _perm.GetTheLengthStartPuneing();
	_perm.GetCurrentOptimalLocalConfiguration(best_config_ever, num_of_best_config, length_of_part_config);
	//������֪�����Ų��ֹ��ͣ���һ���øĽ��㷨�����Ż�
	for (size_t i = 0; i < num_of_best_config; i++){
		//���㵱ǰ��������
		present_energy = CalculatePresentConfigEnergy(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//��ȡ��ǰ����ƽ��Ȩ��
		//double _weigtht = _perm.GetTargetLengthWeight(length_of_part_config);
		start_weigtht = 1;
		//��ȡ��ǰ������
		p_second = best_config_ever[i][length_of_part_config - 1];
		//��ȡ��ǰ����
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

//���㵱ǰ��������ֵ
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