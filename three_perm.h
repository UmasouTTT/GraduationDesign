#pragma once

#include "perm.h"

using namespace std;
using namespace perm_struct;

class three_perm
{
public:
	three_perm();
	~three_perm();

	static const int predict_worest_energy = -45;//Ԥ�����ֵ

private:
	//����Ϊn�Ĺ��͵���ʷƽ��Ȩ��(��Ҫ��ʼ��)
	double *average_weights;
	//����Ϊn�Ĺ��͵�����(��Ҫ��ʼ��)
	double *weights_numbers;
	//����֧���幹��
	point *configurations_point_three;
	char *configurations_class_three;
	//��ǰ��������
	int present_energy;
	//�������
	int lowest_energy = 0;
	//�����������
	point *lowest_configurations_point;
	char *lowest_configurations_class;

	//perm�������
	int perm_lowest_energy = 0;
	//perm�����������
	point *perm_lowest_configurations_point;
	char *perm_lowest_configurations_class;
	//���ε�����ȡ���������ֵ
	int lowest_energy_first;

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
	int LegalActions(point p, point legal_actions[4], int n);
	//�����������Ƿ���ȷ
	bool TestResultIsSatisfied(int target_energy, int length);
	//��ʼ������ʼ����Ԫ��ǰ����ֵΪ��ֵ��
	void InitConfig(string &input, point &p, double &weight);	
	//��������
	void CircleCalculate(int n, int whole_length, point p_before, double weight, string input);
	//��������(ֻ����Ȩ��)
	void CircleCalculateByWeight(int n, int whole_length, point p_before, double weight, string input);
	//����ѡ��ĸ���ȫ�ֱ���
	void  UpdateGlobalVariables(double weight, int n, point p, int tag, char type, int energy_increase);
	//����Cn,Zn
	void UpdateAverageWeight(double w, int n);
	//����Ȩ��
	double CalculateWeight(double w, int energy_increase);
	//������������
	int EnergyIncrease(point p, char type, point p_before, int n);
	//������ʱ����
	void UpdateTempVariables(point _configurations_point[], char _configurations_class[], double weight, int n, point p, char type, int energy_increase);
	//���㵱ǰ��������ֵ
	int CalculatePresentConfigEnergy(point _configurations_point[], char _configurations_class[], int length);
	//������ڹ��͵�Ȩ��
	double CalculateStartConfigWeight(point _configurations_point[], char _configurations_class[], int length);
public:
	//�㷨
	void StartCalculate(string input);
	//�Ľ��㷨��
	void StartCalculateImproveFirst(string input);
	//�����Ľ��㷨
	void Branch_choose_improve_1(string input, int time);
	//�����Ľ��㷨2
	void StartCalculateByWeight(string input);
	//ѭ����������Ľ��㷨2
	void CircleAlgripham2(string input, int times);
};

