#pragma once
#include "stdafx.h"
#include<stdlib.h>
#include<map>
#include<unordered_map>
#include<vector>
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
	~perm();

	static const int max_size_of_input = 100;//�����
	static const int max_size_of_legal_input = 4;//ÿ��������Ϸ�����
	static const int max_size_of_possibleConditions = 24;//
	static const double T;
	//static const double C0;
	static const int Z0;
	static const int C;
	static const double MAX_DOUBLE;

private:
	//���ڿ����˿�
	double C0;
	//��ǰ����֧��ʶ��
	int max_tag = 0;
	//Ȩ������ƽ��ֵ(��Ҫ��ʼ��)
	//double *average_weights = new double[max_size_of_input];
	double average_weights[max_size_of_input];
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	//double *weights_numbers = new double[max_size_of_input];
	double weights_numbers[max_size_of_input];
	//����֧���幹��
	point configurations_point[max_size_of_input];
	char configurations_class[max_size_of_input];
	//����֧��ǰ��������
	int present_energy;
	//vector<int>present_energy;
	//�������
	int lowest_energy = 0;
	//�����������
	point lowest_configurations_point[max_size_of_input];
	char lowest_configurations_class[max_size_of_input];
	//���ڻ�ȡ�Ϸ���ϼ���
	vector<int> input_numbers;
	vector<int> combination_one;
	vector<vector<int>> combination_result;
	//���������������
	int num_of_lowestConfigurations;

public:
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
	//��ȡȨ������ƽ��ֵ
	void GetAverageWeight(double _average_weights[max_size_of_input]);
	//��ȡ����Ϊn�Ĺ��͵�����
	void GetThisWeightNumber(double _weights_numbers[max_size_of_input]);
	//��ȡ��͹�����
	int GetNumOfLowestConfigurations() { return num_of_lowestConfigurations; }
	//�����˿ڷ������
	void SetPopulation(double _C0) { C0 = _C0; }

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
	int LegalActions(point p, vector<pair<int, point>> &legal_actions, int n);
	//**************����ö�*************
	double CalculateGoodResults(point p, char type, point p_before, int energy_increase, int n);
	//**************����Ȩ��*************
	double CalculateWeight(double w, point p, char type, int energy_increase, point p_before);
	//������������ϵ��
	double CalculatingLengthCoefficient(int n, int length);
	//**************����Ԥ��Ȩ�ؼ����������ĺöȣ������ظ����㣩*************(�������ݽ϶࣬����������)
	double CalculatePredictWeightMid(double w, point p_before, char type, vector<double> &good_degrees, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase, int n);
	double CalculatePredictWeight(double w, point p_before, char type, vector<double> &good_degrees, int n, int length, int k_free, const vector<pair<int, point>> &legal_actions, map<point, int> &energy_increase);
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
	point GetNextActionByGoodDegrees(point p_before, vector<double> &good_degrees);
	//�ݹ�����������
	void CalculationCombinations(int offset, int k);
	//��ȡ���������
	vector<vector<int>>GetCombinations(vector<int> &legal_actions, int num);
	//������ֵ��ȡ��Ӧ�Ķ���
	vector<point> GetActionsByNum(vector<int> &numbers, point p_before);
	//************���պöȸ������ѡ��������*****************
	vector<point> ChooseActionsGroupByGoodDegrees(int k, vector<double> &good_degrees, point p_before);
	//�����������Ƿ���ȷ
	bool TestResultIsSatisfied(int target_energy, int length);	
	//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ��
	void InitConfig(string &input, point &p, double &weight, int tag);
	//������󳤶�
	void CalculateMaxSize(int length);
	//��ʼ��ȫ�ֱ�Ԫ
	void InitGlobalVariable(string input);
	//���ɳ�ʼȨ��
	void InitStartAverageWeight(int n, int whole_length, point p_before, double weight, string input, double _average_weights[max_size_of_input], int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _weights_numbers[max_size_of_input]);

public:
	//�����㷨
	void StartCalculate(string input, int num_of_circle);
	//�����������֧���
	void CalculationProcess(int n, int whole_length, int tag, point p_before, double weight, string input);
	//��μ������֧���
	void CircleCalculateProcess(int n, int whole_length, point p_before, double weight, string input, int cirlce_times, double _average_weights[max_size_of_input], int _energy, point points[max_size_of_input], char _points[max_size_of_input], double _weights_numbers[max_size_of_input]);
};

