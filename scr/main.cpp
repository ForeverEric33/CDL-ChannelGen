#include "D:\my_programm\channel_modeling\header\Channel.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
int main()
{
    vector<BSStation> BS;
    vector<UEStation> UE;

    // 创建并设置基站参数
    BSStation BS1;
    BS1.setFrequency(6e9);                                  // 设置基站频率
    BS1.setPosition({0, 0, 5});                             // 设置空间坐标
    BS1.setAntenna({0, 0, 0}, 1, 90, 0, {1, 1, 0.5, 0.64}); // 设置基站天线，(旋转角度，极化方式，极化角度，定向增益，面板排列)
    BS.push_back(BS1);

    // BSStation BS1;
    // BS1.setFrequency(3e9);//设置基站频率
    // BS1.setPosition({10, 0, 0});//设置空间坐标
    // BS1.setAntenna({0, 0, 0}, 2, 45, 5, {1, 1, 0.5, 0.64});//设置基站天线，(旋转角度，极化方式，极化角度，定向增益，面板排列)
    // BS.push_back(BS1);`

    // 创建并设置终端参数
    UEStation UE1;
    UE1.setSpeed({0, 10, 0});                               // 设置终端速度
    UE1.setPosition({10, 10, 2});                           // 设置空间坐标
    UE1.setAntenna({0, 0, 0}, 1, 90, 0, {1, 1, 0.5, 0.64}); // 设置终端天线，(旋转角度，极化方式，极化角度，定向增益，面板排列)
    UE.push_back(UE1);

    // UEStation UE2;
    // UE2.setSpeed({10, 0, 0});//设置终端速度
    // UE2.setPosition({10, 0, 0});//设置空间坐标
    // UE2.setAntenna({20, 20, 20}, 2, 90, 3, {1, 1, 0.5, 0.64});//设置终端天线，(旋转角度，极化方式，极化角度，定向增益，面板排列)
    // UE.push_back(UE2);

    // 设置信道参数
    Channel Ch1;
    Ch1.setChannel(1, "D:/my_programm/channel_modeling/channel_model/901_CDL_A.txt", 100, 10, 4, 1, false, false);
    // 设置仿真场景，信道模型，时延扩展（手动输入且Channel.h中有根据场景功率计算时延扩展的代码），最大移动速度(多终端场景下取最大值)，过采样因子，仿真时间,是否计算路损和阴影
    // 1：UMa    2:UMi   3:RMa   4：InF  5:InF-SL    6:InF-DL    7:InF-SH     8:InF-DH   9:InH

    Ch1.RUN(BS, UE);

    Ch1.writeCIR(BS, UE, "D:/my_programm/channel_modeling/result/my_result2.txt"); // 输出文件名和保存位置

    return 0;
}
