#include "perm.h"


const int perm::max_size_of_input;
const int perm::max_size_of_legal_input;
const int perm::max_size_of_possibleConditions;

const double perm::T = 0.35;
const double perm::C0 = 100000;
const int perm::Z0 = 1;
const int perm::C = 1;

const float perm::part_config_for_save = 0.3;

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
}


perm::~perm()
{
	input_numbers.clear();
	combination_one.clear();
	combination_result.clear();
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

}


//��ȡ����������͵�����
void perm::GetPointPosition(point points[max_size_of_input]) {
	ArrayAssignment(points, lowest_configurations_point, max_size_of_input);
}
//��ȡ����������͵�����
void perm::GetPoint(char points[max_size_of_input]) {
	ArrayAssignment(points, lowest_configurations_class, max_size_of_input);
}

//��������������͵�����
void perm::SetPointPosition(point points[max_size_of_input]) {
	ArrayAssignment(configurations_point, points, max_size_of_input);
}
//��������������͵�����
void perm::SetPoint(char points[max_size_of_input]) {
	ArrayAssignment(configurations_class, points, max_size_of_input);
}

//����Ȩ������ƽ��ֵ
void perm::SetAverageWeight(double _average_weights[max_size_of_input]) {
	ArrayAssignment(average_weights, _average_weights, max_size_of_input);
}

//���ó���Ϊn�Ĺ��͵�����
void perm::SetThisWeightNumber(double _weights_numbers[max_size_of_input]) {
	ArrayAssignment(weights_numbers, _weights_numbers, max_size_of_input);
}

//��ȡ��ǰ���žֲ�����
void perm::GetCurrentOptimalLocalConfiguration(point _best_config_ever[100][perm::max_size_of_input], int num_of_best_config, int length) {
	for (size_t i = 0; i < num_of_best_config; i++){
		ArrayAssignment(_best_config_ever[i], best_config_ever[i], length);
	}
}



//�㷨����
//����������֮��ľ���
float perm::DistenceBetweenPoints(point point1, point point2) {
	float result = (float)(point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
	return sqrtf(result);
}

//������������
int perm::EnergyIncrease(point p, char type, point p_before, int n) {
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

//�жϸ������Ƿ��Ѿ���ʹ��
bool perm::IsThisPositionAlreadyOccupied(point p, int n) {
	for (size_t i = 0; i < n - 1; i++) {
		point _p = configurations_point[i];
		if (p.x == _p.x && p.y == _p.y) {
			return true;
		}
	}
	return false;
}

//����Ϸ��Ķ�����
int perm::LegalActions(point p, int n) {
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
int perm::LegalActions(point p, vector<pair<int, point>> &legal_actions, int n) {
	int result = 0;
	//n+1��Ϊ�϶˷���
	point p1(p);
	p1.y = p1.y + 1;
	if (!IsThisPositionAlreadyOccupied(p1, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p1));
	}
	else {
		legal_actions.push_back(make_pair(0, p1));
	}
	//n+1��Ϊ�Ҷ˷���
	point p2(p);
	p2.x = p2.x + 1;
	if (!IsThisPositionAlreadyOccupied(p2, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p2));
	}
	else {
		legal_actions.push_back(make_pair(0, p2));
	}
	//n+1��Ϊ�¶˷���
	point p3(p);
	p3.y = p3.y - 1;
	if (!IsThisPositionAlreadyOccupied(p3, n)) {
		result += 1;
		legal_actions.push_back(make_pair(1, p3));
	}
	else {
		legal_actions.push_back(make_pair(0, p3));
	}
	//n+1��Ϊ��˷���
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


//**************����ö�*************
double perm::CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n) {
	double result = 0.0;
	int actions_later = LegalActions(p, n);
	result += ((double)actions_later + 0.5) * exp(-energy_increase / T);
	return result;
}

//**************����Ȩ��*************
double perm::CalculateWeight(double w, point p, char type, int energy_increase, point p_before) {
	double result = w * exp(-energy_increase / T);
	return result;
}

//������������ϵ��
double perm::CalculatingLengthCoefficient(int n, int length) {
	if (n <= length * 0.3) {
		return 1;
	}
	if (n > length * 0.3 && n < length * 0.75) {
		return random(30, 45);
	}
	return random(5, 10);
}

//**************����Ԥ��Ȩ�ؼ����������ĺöȣ������ظ����㣩*************(�������ݽ϶࣬����������)
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
	//��ǰһ��Ȩ�����
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

//*************����Cn,Zn***************
void perm::UpdateAverageWeight(double w, int n) {
	double average_weight_before = average_weights[n - 1] * weights_numbers[n - 1];
	++weights_numbers[n - 1];
	average_weights[n - 1] = (average_weight_before + w) / weights_numbers[n - 1];
}

//***************����������***********
double perm::CalculateUpperThreshold(int n) {
	double result = C * (average_weights[n - 1] / Z0) * (weights_numbers[n - 1] / C0) * (weights_numbers[n - 1] / C0);
	return result;
}

//**************����������***********
double perm::CalculateLowerThreshold(double upper_threshold) {
	double result = 0.2 * upper_threshold;
	return result;
}

//*******************�����µķ�֧*****************************
/*void CreateNewBranch(const map<point, char> &config_before, int energy_before) {
//��֧����
configurations.push_back(config_before);
//��֧����
//present_energy = (int *)realloc(present_energy, (max_tag + 2) * sizeof(int));
//present_energy[max_tag + 1] = energy_before;
//��֧��ʶ
++max_tag;
}*/

//*****************����ѡ��ĸ���ȫ�ֱ���***************
int  perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, point point_before[], char type_before[]) {
	//��¼֮ǰ�������͹���
	int energy_before = present_energy;
	ArrayAssignment(point_before, configurations_point, n);
	ArrayAssignment(type_before, configurations_class, n);
	//����Ȩ������ƽ��ֵ�����ֹ��ͳ��ȵ�����
	UpdateAverageWeight(weight, n);
	//���¸���֧���幹��
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
	//���¸���֧��ǰ��������
	present_energy += energy_increase;
	return energy_before;
}

int  perm::UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase) {
	//��¼֮ǰ�������͹���
	int energy_before = present_energy;
	//����Ȩ������ƽ��ֵ�����ֹ��ͳ��ȵ�����
	UpdateAverageWeight(weight, n);
	//���¸���֧���幹��
	configurations_point[n - 1] = p;
	configurations_class[n - 1] = type;
	//���¸���֧��ǰ��������
	present_energy += energy_increase;
	return energy_before;
}

//**********************���ո��������������********************************
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

//�ݹ�����������
void perm::CalculationCombinations(int offset, int k) {
	if (k == 0) {
		combination_result.push_back(combination_one);
		return;
	}
	//ÿ�εݹ������Ҫ����i�ǲ���i <= people.size() - k�����û�м���i++�����i���������������һ�εݹ�
	for (size_t i = offset; i <= input_numbers.size() - k; ++i) {
		combination_one.push_back(input_numbers[i]);
		CalculationCombinations(i + 1, k - 1);
		combination_one.pop_back();//ɾ��combination���һ��Ԫ��
	}
}

//��ȡ���������
vector<vector<int>> perm::GetCombinations(vector<int> &legal_actions, int num) {
	input_numbers.clear();
	combination_one.clear();
	combination_result.clear();
	//��ʼ����������
	for (size_t i = 0; i < legal_actions.size(); i++) {
		input_numbers.push_back(legal_actions[i]);
	}
	//�������������
	CalculationCombinations(0, num);
	return combination_result;
}
//������ֵ��ȡ��Ӧ�Ķ���
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

//************���պöȸ������ѡ��������*****************
vector<point> perm::ChooseActionsGroupByGoodDegrees(int k, vector<double> &good_degrees, point p_before) {
	//�Ϸ���������
	vector<int>legal_actions;
	for (size_t i = 0; i < good_degrees.size(); i++) {
		if (good_degrees[i] - 0.0 < judge_is_zero && good_degrees[i] - 0.0 > -judge_is_zero) {
			continue;
		}
		legal_actions.push_back(i);
	}
	//������������
	vector<vector<int>>combination_actions = GetCombinations(legal_actions, k);
	//����ö��ܺ�
	vector<pair<double, double>>combinations_sum_section;
	//��ǰ��������
	double present_good_degrees_sum = 0;
	for (size_t i = 0; i < combination_actions.size(); i++) {
		double temp_good_degree_sum = 0;
		//���������ϵĺöȺ�
		for (size_t j = 0; j < combination_actions[i].size(); j++) {
			temp_good_degree_sum += good_degrees[combination_actions[i][j]];
		}
		combinations_sum_section.push_back(make_pair(present_good_degrees_sum, present_good_degrees_sum + temp_good_degree_sum));
		//������������
		present_good_degrees_sum += temp_good_degree_sum;
	}
	//���ѡȡ
	double random_result = random(0, present_good_degrees_sum);
	//�ҵ�ѡȡ�ļ���
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
//�����������Ƿ���ȷ
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

//�����������֧���
void perm::CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input) {
	//���������ж�
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
	vector<pair<int, point>>legal_actions;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	int k_free = LegalActions(p_before, legal_actions, n);
	if (k_free == 0) {
		return;
	}
	//���������ö�
	vector<double>good_degrees;
	//�����������µ���������
	map<point, int>energy_increase;
	//������������ĺö���Ȩ��Ԥ��ֵ
	double predict_wigtht = CalculatePredictWeight(weight, p_before, input[n - 1], good_degrees, n, whole_length, k_free, legal_actions, energy_increase);
	//������������
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
	//����Ԥ��ֵ���������޵���ֵ��ϵ��������
	if (predict_wigtht >= lower_threshold && predict_wigtht <= upper_threshold) {
		//��ȡ��֦��ʼ����
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//���ݺöȸ���ѡ����һ����
		point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
		//��������ö�����Ȩ��
		double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
		//����
		UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action]);
		//�����֧
		CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
	}
	else if (predict_wigtht < lower_threshold) {
		//��ȡ��֦��ʼ����
		if (!isPuneBegin && n < beginPuningTheBranch) {
			beginPuningTheBranch = n;
			choose_config_length = beginPuningTheBranch + 3;
		}
		//����1/2�ĸ��ʶ����÷�֧
		double rand_result = random(0, 1);
		if (rand_result < 0.5) {
			return;
		}
		else {
			//���ݺöȸ���ѡ����һ����
			point next_action = GetNextActionByGoodDegrees(p_before, good_degrees);
			//��������ö�����Ȩ��
			double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
			//����
			UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action]);
			//�����֧
			CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
		}
	}
	else {
		//�����֧����
		int k = Min((double)k_free, (double)(predict_wigtht / upper_threshold));
		//���ݺöȸ���ѡ����һ��������
		vector<point>choose_actions = ChooseActionsGroupByGoodDegrees(k, good_degrees, p_before);
		//��¼����ǰ��ֵ
		int energy_before;
		//���ݸ����������µķ�֧
		for (size_t i = 0; i < choose_actions.size(); i++) {
			point next_action = choose_actions[i];
			if (i == 0) {//�����½���֧
						 //��������ö�����Ȩ��
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
				//����
				energy_before = UpdateGlobalVariables(present_weight, n, next_action, tag, input[n - 1], energy_increase[next_action], point_before, type_before);
				//�����֧
				CalculationProcess(n + 1, whole_length, tag, next_action, present_weight, input);
			}
			else {//�½���֧
				  //��������ö�����Ȩ��
				double present_weight = CalculateWeight(weight, next_action, input[n - 1], energy_increase[next_action], p_before);
				//�����·�֧
				present_energy = energy_before;
				ArrayAssignment(configurations_point, point_before, n);
				ArrayAssignment(configurations_class, type_before, n);
				//����
				UpdateGlobalVariables(present_weight, n, next_action, max_tag, input[n - 1], energy_increase[next_action]);
				//�����֧
				CalculationProcess(n + 1, whole_length, max_tag, next_action, present_weight, input);
			}
		}
	}
}

//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ��
void perm::InitConfig(string &input, point &p, double &weight) {
	//�������
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
	max_tag = 0;
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
	configurations_point[0] = p1;
	configurations_class[0] = input[0];
	p1.x = p1.x + 1;
	p = p1;
	configurations_point[1] = p1;
	configurations_class[1] = input[1];
	//����֧��ǰ��������
	present_energy = 0;
	//��ʼ��Ȩ��
	CalculationProcess(3, input.length(), 0, p, weight, input);
	if (TestResultIsSatisfied(lowest_energy, input.length())) {
		cout << "test satisfied!" << endl;
	}
	else {
		cout << "something wrongQAQ~" << endl;
	}
	//������������
	max_tag = 0;
	weight = 1;
	//����֧��ǰ��������
	present_energy = 0;
	//��û���ƶ�������������ɵ�һ�ε����ṹ����
	if (worest_energy == 0) {
		worest_energy == lowest_energy;
	}
}


//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ���������ɳ�ʼȨ�أ�
void perm::InitConfigWithoutInitWeight(string &input, point &p, double &weight) {
	//�������
	//free(present_energy);
	//present_energy = (int *)malloc(sizeof(int));
	max_tag = 0;
	weight = 1;
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
	//��ʼ����Ԫ
	int tag_i = 0;
	point p_second;
	double start_weigtht;
	int result_energy_low = 0;
	point point_before[max_size_of_input];
	char type_before[max_size_of_input];
	//��ȡ�����������๹����
	num_of_lowestConfigurations = 0;
	int temp_lowestConfig = 0;
	//��ȡʱ��
	/*struct tm t;   //tm�ṹָ��
	time_t now;  //����time_t���ͱ���
	time(&now);      //��ȡϵͳ���ں�ʱ��
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
		//�������˲��Խ�Ȩ��������Ϊ1
		for (size_t i = 0; i < max_size_of_input; i++) {
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}


//��μ������֧���
void perm::CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, double _average_weights[max_size_of_input], int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _weights_numbers[max_size_of_input]){
	//��ʼ���ֲ�����
	int tag_i = 0;
	int result_energy_low = 0;	
	num_of_lowestConfigurations = 0;
	int temp_lowestConfig = 0;
	//ѭ����
	while (tag_i < cirlce_times) {
		//����perm����
		if (tag_i == 0) {
			InitStartAverageWeight(n, whole_length, p_before, weight, input, _energy, points, _points);
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
		else {
			cout << "something wrongQAQ~" << endl;
		}
		if (lowest_energy < result_energy_low) {
			result_energy_low = lowest_energy;
			cout << "lowest energy: " << result_energy_low << "length of config : " << input.length() << endl;
			//temp_lowestConfig = num_of_lowestConfigurations;
			//num_of_lowestConfigurations = 1;
		}
		//ѡ���������������
	
		/*t = time(NULL);
		char ch[64] = { 0 };
		strftime(ch, sizeof(ch) - 1, "%Y-%m-%d %H:%M:%S", localtime(&t));*/
		//cout << ch << endl;

		//�������˲��Խ�Ȩ��������Ϊ1
		for (size_t i = 0; i < max_size_of_input; i++){
			weights_numbers[i] = 1;
		}
		++tag_i;
	}
	//num_of_lowestConfigurations = temp_lowestConfig;
}

//��ʼ����ʼȨ��
void perm::InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, int _energy, point points[max_size_of_input], char _points[max_size_of_input]) {
	max_tag = 0;
	weight = 1;
	//Ȩ������ƽ��ֵ(��һ�ε���ʱ��ȡ��ʼֵ��Ҫ��ʼ��)
	for (size_t i = 0; i < input.length(); i++) {
		average_weights[i] = 0;
	}
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	for (size_t i = 0; i < input.length(); i++) {
		weights_numbers[i] = 0;
	}
	//����֧���幹��
	SetPointPosition(points);
	SetPoint(_points);
	//����֧��ǰ��������
	SetEnergy(_energy);
	//��ʼ��Ȩ��
	CalculationProcess(n, input.length(), 0, p_before, weight, input);
	if (TestResultIsSatisfied(lowest_energy, input.length())) {
		cout << "finish start weight data!" << endl;
	}
	else {
		cout << "something wrongQAQ~" << endl;
	}
	//������������
	max_tag = 0;
}


//�ж����ֹ����Ƿ������Ϊһ�£�Ŀǰ��Ϊ���ֹ��ͱ�������ѡ����һ�¼���Ϊ���ֹ���һ�£�
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


//����ʷ���Ź���������¹���
void perm::AddNewConfigToBestConfigEver(point points[max_size_of_input], int length) {
	//���ڻ�ȡ�������Ź��͵�perm�㷨ִ������һ�����ڵ���ѡ���֧��Perm�㷨ִ��������
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