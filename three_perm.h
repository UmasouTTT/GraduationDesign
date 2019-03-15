#pragma once

#include "perm.h"

using namespace std;
using namespace perm_struct;

class three_perm
{
public:
	three_perm();
	~three_perm();



private:
	double average_weights[perm::max_size_of_input];
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	double weights_numbers[perm::max_size_of_input];
	//����֧���幹��
	point configurations_point[perm::max_size_of_input];
	char configurations_class[perm::max_size_of_input];
	//��ǰ��������
	int present_energy;
	//�������
	int lowest_energy = 0;
	//�����������
	point lowest_configurations_point[perm::max_size_of_input];
	char lowest_configurations_class[perm::max_size_of_input];

	//perm�������
	int perm_lowest_energy = 0;
	//perm�����������
	point perm_lowest_configurations_point[perm::max_size_of_input];
	char perm_lowest_configurations_class[perm::max_size_of_input];

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
	//�жϸ������Ƿ��Ѿ���ʹ��
	bool IsThisPositionAlreadyOccupied(point p, int n);
	//����Ϸ��Ķ�����
	int LegalActions(point p, int n);
	//����������֮��ľ���
	float DistenceBetweenPoints(point point1, point point2);
	//�ع�����Ϸ�������������߼�������
	int LegalActions(point p, vector<point> &legal_action, int n);
	//�����������Ƿ���ȷ
	bool TestResultIsSatisfied(int target_energy, int length);
	//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ��
	void InitConfig(string &input, point &p, double &weight);	
	//��������
	void CircleCalculate(int n, int whole_length, point p_before, double weight, string input);
	//����ѡ��ĸ���ȫ�ֱ���
	void  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase, double _average_weights[], double _weights_numbers[]);
	//����Cn,Zn
	void UpdateAverageWeightByThree(double w, int n, double _average_weights[], double _weights_numbers[]);
	//����Ȩ��
	double CalculateWeight(double w, int energy_increase);
	//������������
	int EnergyIncrease(point p, char type, point p_before, int n);
	//������ʱ����
	void UpdateTempVariables(double _average_weights[], double _weights_numbers[], point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase);
public:
	//�㷨
	void StartCalculate(string input);
};

