#pragma once
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<string>
#include<time.h>
#include<iostream>

namespace perm_struct {
	struct point {
		int x;
		int y;
		bool operator < (const point &p) const {
			return x < p.x || (x == p.x && y < p.y);
		}
	};
}


#define random(a,b) (((double)rand()/RAND_MAX)*(b - a) + a)
#define judge_is_zero 0.000001

using namespace std;
using namespace perm_struct;

class perm
{
public:
	perm();
	perm(int predict_energy);
	perm(int predict_energy, bool isSoultionSaved);
	~perm();

	static const int max_size_of_input = 100;//�����
	static const int max_size_of_legal_input = 4;//ÿ��������Ϸ�����
	static const int max_size_of_possibleConditions = 24;//
	static const double T;
	static const double C0;
	static const int Z0;
	static const int C;
	static const float part_config_for_save;//��ʼ���������湹�ͱ���
	static const int target_lowest_energy;//Ŀ����͹���

private:
	//��ǰ����֧��ʶ��
	int max_tag = 0;
	//Ȩ������ƽ��ֵ(��Ҫ��ʼ��)
	//double *average_weights = new double[max_size_of_input];
	double *average_weights; 
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	//double *weights_numbers = new double[max_size_of_input];
	double *weights_numbers;
	//����֧���幹��
	point *configurations_point;
	char *configurations_class;
	//����֧��ǰ��������
	int present_energy;
	//vector<int>present_energy;
	//�������
	int lowest_energy = 0;
	//�����������
	point *lowest_configurations_point;
	char *lowest_configurations_class;
	//���ڻ�ȡ�Ϸ���ϼ���

	int **combinations_result;
	point *choose_actions_3;
	//���������������
	int num_of_lowestConfigurations;
	//Ĭ�������������ɵ�һ�������������
	int worest_energy;
	//����Perm���������ͳ���
	int choose_config_length = 0;
	//��¼Ŀǰ���ֵ����Ź���
	point **best_config_ever;
	//��¼Ŀǰ���Ź�������
	int best_config_num;
	//���ɳ�ʼȨ��
	bool InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, int _energy, point points[max_size_of_input], char _points[max_size_of_input]);
	//���ڼ�¼��ʼ���м�֦ʱ������
	int beginPuningTheBranch;
	//�����ж��Ƿ��Ѿ���ʼ��֦
	bool isPuneBegin;
	//�����ж��Ƿ���Ҫ��¼���Ž��֧
	bool isSolutionNeedToBeSaved;

public:
	//��ȡƽ��Ȩ��
	void GetAverageWeight(double _average_weights[max_size_of_input]);
	//��ȡ����Ϊn�Ĺ��͵�����
	void GetWeightNumber(double _weights_numbers[max_size_of_input]);
	//��ȡ����ֵ
	int GetEnergy() { return lowest_energy; }
	//��ȡ����������͵�����
	void GetPointPosition(point points[max_size_of_input]);
	//��ȡ����������͵�����
	void GetPoint(char points[max_size_of_input]);
	//��������ֵ
	void SetEnergy(int energy) { present_energy = energy; }
	//��������������͵�����
	void SetPointPosition(point points[max_size_of_input]);
	//��������������͵�����
	void SetPoint(char points[max_size_of_input]);
	//����Ȩ������ƽ��ֵ
	void SetAverageWeight(double _average_weights[max_size_of_input]);
	//���ó���Ϊn�Ĺ��͵�����
	void SetThisWeightNumber(double _weights_numbers[max_size_of_input]);
	//��ȡ��͹�����
	int GetNumOfLowestConfigurations() { return num_of_lowestConfigurations; }
	//��ȡ��Ϳ��ù�������
	int GetNumOfLowestConfigAndCanUse() { return best_config_num; }
	//��ȡ��ǰ���žֲ�����
	void GetCurrentOptimalLocalConfiguration(point _best_config_ever[100][perm::max_size_of_input], int num_of_best_config, int length);
	//��ȡָ������ƽ��Ȩ��(for test)
	double GetTargetLengthWeight(int _index) { return average_weights[_index - 1]; }
	//��ȡ��ʼ��֦ʱ������
	int GetTheLengthStartPuneing() { return beginPuningTheBranch; }

private:
	//��Сֵ
	template <typename T>
	T Min(T num1, T num2) {
		if (num1 < num2) {
			return num1;
		}
		return num2;
	}
	//���鸳ֵ
	template <typename T>
	void ArrayAssignment(T number1[], T number2[], int length) {
		for (size_t i = 0; i < length; i++) {
			number1[i] = number2[i];
		}
	}
	//����������֮��ľ���
	float DistenceBetweenPoints(point point1, point point2);
	//������������
	int EnergyIncrease(point p, char type, point p_before, int n);
	//�жϸ������Ƿ��Ѿ���ʹ��
	bool IsThisPositionAlreadyOccupied(point p, int n);
	//����Ϸ��Ķ�����
	int LegalActions(point p, int n);
	//�ع�����Ϸ�������������߼�������
	int LegalActions(point p, point legal_actions_p[4], int legal_actions_t[4], int n);
	//**************����ö�*************
	double CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n);
	//**************����Ȩ��*************
	double CalculateWeight(double w, point p, char type, int energy_increase, point p_before);
	//������������ϵ��
	double CalculatingLengthCoefficient(int n, int length);
	//**************����Ԥ��Ȩ�ؼ����������ĺöȣ������ظ����㣩*************(�������ݽ϶࣬����������)
	double CalculatePredictWeightMid(double w, point p_before, char type, double good_degrees[4], int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4], int n);
	double CalculatePredictWeight(double w, point p_before, char type, double good_degrees[4], int n, int length, int k_free, point legal_actions_p[4], int legal_actions_t[4], int energy_increase[4]);
	//*************����Cn,Zn***************
	void UpdateAverageWeight(double w, int n);
	//***************����������***********
	double CalculateUpperThreshold(int n);
	//**************����������***********
	double CalculateLowerThreshold(double upper_threshold);
	//*****************����ѡ��ĸ���ȫ�ֱ���***************
	int  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, point point_before[], char type_before[]);
	int  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase);
	//**********************���ո��������������********************************
	point GetNextActionByGoodDegrees(point p_before, double good_degrees[4]);
	//�ݹ�����������
	void CalculationCombinations(int offset, int k);
	//��ȡ���������
	void GetCombinations(int(&legal_actions)[4], int num, int num_of_legal_actions, int &num_of_result);
	//������ֵ��ȡ��Ӧ�Ķ���
	void GetActionsByNum(int index, point p_before, int length_of_result);
	//************���պöȸ������ѡ��������*****************
    void ChooseActionsGroupByGoodDegrees(int k, double good_degrees[4], point p_before);
	//�����������Ƿ���ȷ
	bool TestResultIsSatisfied(int target_energy, int length);	
	//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ,��Ҫ���ɳ�ʼȨ�أ�
	void InitConfig(string &input, point &p, double &weight);
	//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ���������ɳ�ʼȨ�أ�
	void InitConfigWithoutInitWeight(string &input, point &p, double &weight);
	//������󳤶�
	void CalculateMaxSize(int length);
	//��ʼ��ȫ�ֱ�Ԫ
	void InitGlobalVariable(string input);
	//��ȡ�������
	int GetActionsNum(point actions[4], point action);

public:
	//�����㷨
	void StartCalculate(string input, int num_of_circle);
	//�����㷨(ֻ����Ȩ��)
	void StartCalculateOnlyWeight(string input, int num_of_circle);
	//�����������֧�������¼��֧������
	void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input);
	//�����������֧���������¼��֧������
	void CalculationProcessOnlyConsiderWeight(int n, int whole_length, int tag, point p_before, double weight, string input);
	//��μ������֧���
	void CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]);
	//��μ������֧���(ֻ����Ȩ��)
	void CircleCalculateProcessOnlyConsiderWeight(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _average_weight[max_size_of_input], double _weight_number[max_size_of_input]);
	//��ȡָ������ƽ��Ȩ��ֵ
	double GetAverageWeightWithLengthN(int n) { return average_weights[n-1]; }

	//���㷨���ֹ���
private:
	//����ʷ���Ź���������¹���
	void AddNewConfigToBestConfigEver(point points[max_size_of_input], int length);
	//�ж����ֹ����Ƿ������Ϊһ�£�Ŀǰ��Ϊ���ֹ��ͱ�������ѡ����һ�¼���Ϊ���ֹ���һ�£�
	bool IsTwoConfigTheSame(point points[max_size_of_input], point _points[max_size_of_input], int start, int end);
};

