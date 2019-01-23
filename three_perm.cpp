#include "stdafx.h"
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
		point _p = configurations_point[i];
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
		average_weights[i] = 1;
	}
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 1;
	}
	//����֧���幹��
	point p1;
	p1.x = 0;
	p1.y = 0;
	configurations_point[0] = p1;
	configurations_class[0] = input[0];
	p1.x = p1.x + 1;
	p = p1;
	configurations_point[1] = p1;
	configurations_class[1] = input[1];
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
		cout << present_energy << endl;
		cout << present_energy << endl;
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
	//�ֱ����������ж�����perm�㷨�µ�ֵ������ȡ������͵�
	for (size_t i = 0; i < k_free; i++) {
		perm _perm;
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
		//����perm����
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
	//����ѡ���������и���
	UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
	CircleCalculate(n + 1, whole_length, legal_actions[best_index], _present_weight, input);

}


void  three_perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//����Ȩ������ƽ��ֵ�����ֹ��ͳ��ȵ�����
	UpdateAverageWeight(weight, n);
	//���¸���֧���幹��
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
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
		point _point = configurations_point[i];
		//���������ڲ�Ӱ������
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

//������ʱ����
void three_perm::UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//��ֵ
	ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point, perm::max_size_of_input);
	ArrayAssignment(_configurations_class, configurations_class, perm::max_size_of_input);
	//����
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	++_weights_numbers[n - 1];
	_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}