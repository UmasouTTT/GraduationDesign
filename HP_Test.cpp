// HP_Test.cpp: 定义控制台应用程序的入口点。
//


//#include "stdafx.h"
#include "three_perm.h"

using namespace std;
using namespace perm_struct;

string ReverseString(string input) {
	string result;
	for (size_t i = input.length(); i > 0; i--){
		result += input[i - 1];
	}
	return result;
}

int main()
{
	string input_string = "PPHPPHHPPPPHHPPPPHHPPPPHH";
	string input_string1 = "PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPPHHHHHHPPHHHHHHPPPPPPPPPHPHHPHHHHHHHHHHHPPHHHPHHPHPPHPHHHPPPPPPHHH";
	string input_string2 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP";
	string input_64 = "HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH";
	string input_48 = "PPPPPPHPHHPPPPPHHHPHHHHHPHHPPPPHHPPHHPHHHHHPHHHHHHHHHHPHHPHHHHHHHPPPPPPPPPPPHHHHHHHPPHPHHHPPPPPPHPHH";
	srand((int)time(0));
	/*perm _perm;
	
	
	//StartCalculate(input_string, 10);
	_perm.StartCalculate(input_string, 1);
	//StartCalculate(input_string1, 100);
	//获取能量值
	int energy = _perm.GetEnergy();
	//获取最低能量构型点坐标
	point points_position[100];
	_perm.GetPointPosition(points_position);
	//获取最低能量构型点类型
	char points[100];
	_perm.GetPoint(points);*/

	struct tm t;   //tm结构指针
	time_t now;  //声明time_t类型变量
	time(&now);      //获取系统日期和时间
	localtime_s(&t, &now);   //获取当地日期和时间
	string present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
	cout << present_time << endl;



	three_perm _3_perm;
//	_3_perm.StartCalculate(input_string);
	string re_string1 = ReverseString(input_48);
	//_3_perm.StartCalculate(input_64);
	//_3_perm.Branch_choose_improve_1(re_string1, 100);
	_3_perm.CircleAlgripham2(re_string1, 100);


	//perm _perm(-37);
	//_perm.StartCalculate(input_64, 2);
	
	time(&now);      //获取系统日期和时间
	localtime_s(&t, &now);   //获取当地日期和时间
	present_time = to_string(t.tm_hour) + 'h' + to_string(t.tm_min) + 'm' + to_string(t.tm_sec) + 's';
	cout << present_time << endl;



		return 0;
}