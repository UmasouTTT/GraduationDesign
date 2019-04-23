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
int three_perm::LegalActions(point p, point legal_actions[4], int n) {
	int result = 0;
	//n+1��Ϊ�϶˷���
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		legal_actions[result] = p1;
		result += 1;		
	}
	//n+1��Ϊ�Ҷ˷���
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		legal_actions[result] = p2;
		result += 1;
	}
	//n+1��Ϊ�¶˷���
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		legal_actions[result] = p3;
		result += 1;
	}
	//n+1��Ϊ��˷���
	point p4(p);
	p4.x = p4.x - 1;
	if (!IsThisPositionAlreadyOccupied(p4, n)) {
		legal_actions[result] = p4;
		result += 1;
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
	perm _perm(predict_worest_energy);
	_perm.StartCalculateOnlyWeight(input, 1);
	//��ȡ��ǰƽ����������ʷ����������
	_perm.GetAverageWeight(average_weights);
	_perm.GetWeightNumber(weights_numbers);
	//��ȡ���ε�����������ֵ
	lowest_energy_first = _perm.GetEnergy();
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
		if (present_energy == perm::target_lowest_energy) {
			cout << "find target config!" << endl;
			struct tm t;   //tm�ṹָ��
			time_t now;  //����time_t���ͱ���
			time(&now);      //��ȡϵͳ���ں�ʱ��
			localtime_s(&t, &now);   //��ȡ�������ں�ʱ��
			string present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
			cout << present_time << endl;
			cout << "end" << endl;
		}
		return;
	}
	//��ȡ��ǰ״̬���еĶ���
	point legal_actions[4];
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
	//�����ж��Ƿ��ҵ��ȳ��ε���������õķ�֦
	bool _if_find_better_than_ever = false;
	//���ڼ�¼�÷�֧��ǰ��������ֵ
	int _energy_present = present_energy;
	//�ֱ����������ж�����perm�㷨�µ�ֵ������ȡ������͵�
	for (size_t i = 0; i < k_free; i++) {
		perm _perm(-43, false);
		//������ѡ��õ������
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//����֧���幹��
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//����
		UpdateTempVariables(temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		//���ݽ���ѡ�������������ʼʱ���������࣬Խ����Խ��
		int circle_times = 3;
		_perm.CircleCalculateProcess(n + 1, whole_length, legal_actions[i], present_weight, input, circle_times, _energy_present + energy_increase, temp_configurations_point, temp_configurations_class, average_weights, weights_numbers);
		//��¼�¸ôε���֮���Ȩ��ֵ
		_perm.GetAverageWeight(average_weights);
		_perm.GetWeightNumber(weights_numbers);
		//����ҵ��ȳ��ε������õĽ����ֱ�ӽ���÷�֧
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
			//��ȡ��ǰ��������ֵ
			present_energy = _energy_present;
			//����ѡ���������и���
			UpdateGlobalVariables(present_weight, n, legal_actions[i], 0, input[n - 1], energy_increase);
			//�����µķ�֧
			CircleCalculate(n + 1, whole_length, legal_actions[i], present_weight, input);
		}
		else if(!_if_find_better_than_ever){
			//�����ǰ��֧���ҵ�����С������֮ǰ�Ķ�С���򽫸÷�֧ѡΪ���ŷ�֧
			if (_perm.GetEnergy() < min_energy) {
				best_index = i;
				//��ȡ���ε������ٹ�����
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
			//�������֧���ҵ�����С����һ�£����������ַ�֧���õ������Ź��������������ѡ���֧��������һ�£������ѡ
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
	//���û�г�����ѹ��ͣ�����ѡ���������и���
	if (!_if_find_better_than_ever) {
		//��ȡ��ǰ��������ֵ
		present_energy = _energy_present;
		UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
		CircleCalculate(n + 1, whole_length, legal_actions[best_index], _present_weight, input);
	}	
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
	average_weights[n - 1] = (average_weights[n - 1] + w) / 2;
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
void three_perm::UpdateTempVariables(point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase) {
	//��ֵ
	//ArrayAssignment(_average_weights, average_weights, perm::max_size_of_input);
	//ArrayAssignment(_weights_numbers, weights_numbers, perm::max_size_of_input);
	ArrayAssignment(_configurations_point, configurations_point_three, n-1);
	ArrayAssignment(_configurations_class, configurations_class_three, n-1);
	//����
	_configurations_point[n - 1] = p;
	_configurations_class[n - 1] = type;
	//double average_weight_before = _average_weights[n - 1] * _weights_numbers[n - 1];
	//++_weights_numbers[n - 1];
	//_average_weights[n - 1] = (average_weight_before + weight) / _weights_numbers[n - 1];
}




//�Ľ��㷨��
void three_perm::StartCalculateImproveFirst(string input) {
	//��ʼ��
	//perm_lowest_energy = 0;
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[perm::max_size_of_input];
	char type_before[perm::max_size_of_input];
	InitConfig(input, p_second, start_weigtht);
	//���ȵ���perm�㷨���ڿ�ʼ�׶ξ����ܵĽ���չ��
	perm _perm(predict_worest_energy);
	_perm.StartCalculate(input, 5);
	//��ȡ��ǰƽ����������ʷ����������
	_perm.GetAverageWeight(average_weights);
	_perm.GetWeightNumber(weights_numbers);
	//��ȡ���ε�����������ֵ
	lowest_energy_first = _perm.GetEnergy();
	//��ȡ��ǰ����perm�㷨��ȡ�����Ų��ֹ�������
	int num_of_best_config = _perm.GetNumOfLowestConfigAndCanUse();
	//��ȡ��֪���й���
	point best_config_ever[100][perm::max_size_of_input];
	char best_config_ever_points[perm::max_size_of_input];
	_perm.GetPoint(best_config_ever_points);
	ArrayAssignment(configurations_class_three, best_config_ever_points, perm::max_size_of_input);
	//��ʹ�ù��ͳ���
	int length_of_part_config = _perm.GetTheLengthStartPuneing();
	_perm.GetCurrentOptimalLocalConfiguration(best_config_ever, num_of_best_config, length_of_part_config);
	//������֪�����Ų��ֹ��ͣ���һ���øĽ��㷨�����Ż�
	for (size_t i = 0; i < num_of_best_config; i++){
		//���㵱ǰ��������
		present_energy = CalculatePresentConfigEnergy(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//��ȡ��ǰ����ƽ��Ȩ��

		//double _weigtht = _perm.GetTargetLengthWeight(length_of_part_config);
		start_weigtht = CalculateStartConfigWeight(best_config_ever[i], best_config_ever_points, length_of_part_config);
		//��ȡ��ǰ������
		p_second = best_config_ever[i][length_of_part_config - 1];
		//��ȡ��ǰ����
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




//������ڹ���Ȩ��
double three_perm::CalculateStartConfigWeight(point _configurations_point[], char _configurations_class[], int length) {
	double result = 1;
	for (size_t i = 2; i < length; i++){
		//������������
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





//�㷨(Ȩ��)
void three_perm::StartCalculateByWeight(string input) {
	//��ʼ��
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


//�������̣�Ȩ�أ�
void three_perm::CircleCalculateByWeight(int n, int whole_length, point p_before, double weight, string input) {
	//���������ж�
	if (n > whole_length) {
		cout << "circle end, present energy :";
		cout << present_energy << endl;
		cout << "circle end, history energy :";
		cout << perm_lowest_energy << endl;
		if (present_energy == perm::target_lowest_energy) {
			cout << "find target config!" << endl;
			struct tm t;   //tm�ṹָ��
			time_t now;  //����time_t���ͱ���
			time(&now);      //��ȡϵͳ���ں�ʱ��
			localtime_s(&t, &now);   //��ȡ�������ں�ʱ��
			string present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
			cout << present_time << endl;
			cout << "end" << endl;
		}
		return;
	}
	//��ȡ��ǰ״̬���еĶ���
	point legal_actions[4];
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
	//�����ж��Ƿ��ҵ��ȳ��ε���������õķ�֦
	bool _if_find_better_than_ever = false;
	//���ڼ�¼�÷�֧��ǰ��������ֵ
	int _energy_present = present_energy;
	//��ǰ���Ȩ��
	double best_weight = 0;
	//��ǰ��͹�������
	int _lowest_energy = 0;
	//��¼���ƽ��Ȩ��ֵ
	double best_average_weight[perm::max_size_of_input];
	//�ֱ����������ж�����perm�㷨�µ�ֵ������ȡ������͵�
	cout << "����֧��������������Ȩ��" << endl;
	for (size_t i = 0; i < k_free; i++) {
		perm _perm;
		//������ѡ��õ������
		int energy_increase = EnergyIncrease(legal_actions[i], input[n - 1], p_before, n);
		double present_weight = CalculateWeight(weight, energy_increase);
		//����֧���幹��
		point temp_configurations_point[perm::max_size_of_input];
		char temp_configurations_class[perm::max_size_of_input];
		//����
		UpdateTempVariables(temp_configurations_point, temp_configurations_class, present_weight, n, legal_actions[i], input[n - 1], energy_increase);
		//���ݽ���ѡ�������������ʼʱ���������࣬Խ����Խ��
		int circle_times = 1;
		_perm.CircleCalculateProcessOnlyConsiderWeight(n + 1, whole_length, legal_actions[i], present_weight, input, circle_times, _energy_present + energy_increase, temp_configurations_point, temp_configurations_class, average_weights, weights_numbers);
		cout << _perm.GetEnergy();
		cout << "     ";
		cout << _perm.GetAverageWeightWithLengthN(whole_length - 1);
		cout<<""<<endl;
		if (i == 0) {
			//��¼���Ȩ��
			best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
			//��¼��ǰ��ѷ�֧
			best_index = i;
			//��¼�������
			_lowest_energy = _perm.GetEnergy();
			//��ȡ��Ȩ��
			_perm.GetAverageWeight(best_average_weight);
			//�ж��Ƿ��ҵ���������
			if (_perm.GetEnergy() < min_energy) {
				min_energy = _perm.GetEnergy();
				if (_perm.GetEnergy() < perm_lowest_energy) {
					perm_lowest_energy = _perm.GetEnergy();
					_perm.GetPoint(perm_lowest_configurations_class);
					_perm.GetPointPosition(perm_lowest_configurations_point);
				}
			}
			//��ȡ��ǰ��֧����������Ȩ��
			_energy_increase = energy_increase;
			_present_weight = present_weight;
		}
		//���ȣ����ڱ���֪��������߳���1�ķ�֧������
		else if (_perm.GetEnergy() - _lowest_energy < 1) {
			//�����ǰ���ֵ����������֮ǰ��֧��������ͳ���1��ֱ��ѡ��
			if (_lowest_energy - _perm.GetEnergy() > 1) {
				//��¼���Ȩ��
				best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
				//��¼��ǰ��ѷ�֧
				best_index = i;
				//��¼�������
				_lowest_energy = _perm.GetEnergy();
				//��ȡ��Ȩ��
				_perm.GetAverageWeight(best_average_weight);
				//�ж��Ƿ��ҵ���������
				if (_perm.GetEnergy() < min_energy) {
					min_energy = _perm.GetEnergy();
					if (_perm.GetEnergy() < perm_lowest_energy) {
						perm_lowest_energy = _perm.GetEnergy();
						_perm.GetPoint(perm_lowest_configurations_class);
						_perm.GetPointPosition(perm_lowest_configurations_point);
					}
				}
				//��ȡ��ǰ��֧����������Ȩ��
				_energy_increase = energy_increase;
				_present_weight = present_weight;
			}
			else {
				//Ȼ�󣬰��ձ�������ѡȡ
				int first_log = log(best_weight);
				int second_log = log(_perm.GetAverageWeightWithLengthN(whole_length - 1));
				//����Ȩ��ȡlogֵ�����Զ��ڸ�Ȩ�ط�֧���в���
				if (first_log > second_log) {
					first_log = (first_log - second_log) * second_log * n;
				}
				else {
					second_log = (second_log - first_log) * first_log * n;
				}
				int random_num = random(0, first_log + second_log);
				if (random_num > first_log) {
					//��¼���Ȩ��
					best_weight = _perm.GetAverageWeightWithLengthN(whole_length - 1);
					//��¼��ǰ��ѷ�֧
					best_index = i;
					//��¼�������
					_lowest_energy = _perm.GetEnergy();
					//��ȡ��Ȩ��
					_perm.GetAverageWeight(best_average_weight);
					//�ж��Ƿ��ҵ���������
					if (_perm.GetEnergy() < min_energy) {
						min_energy = _perm.GetEnergy();
						if (_perm.GetEnergy() < perm_lowest_energy) {
							perm_lowest_energy = _perm.GetEnergy();
							_perm.GetPoint(perm_lowest_configurations_class);
							_perm.GetPointPosition(perm_lowest_configurations_point);
						}
					}
					//��ȡ��ǰ��֧����������Ȩ��
					_energy_increase = energy_increase;
					_present_weight = present_weight;
				}
			}			
		}
	}
	//�����֧
	cout << "ѡ���֧���:";
	cout << best_index << endl;
	ArrayAssignment(average_weights, best_average_weight, whole_length);
	UpdateGlobalVariables(_present_weight, n, legal_actions[best_index], 0, input[n - 1], _energy_increase);
	CircleCalculateByWeight(n + 1, whole_length, legal_actions[best_index], _present_weight, input);
}

//ѭ���Ľ��㷨2
void three_perm::CircleAlgripham2(string input, int time) {
	perm_lowest_energy = 0;
	for (size_t i = 0; i < time; i++) {
		StartCalculateByWeight(input);
	}
}
