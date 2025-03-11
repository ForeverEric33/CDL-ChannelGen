#pragma once
#include "Antenna.h"
// 缺少发送功率
class BSStation
{
public:
    struct Position // 空间坐标结构体
    {
        double pos_x, pos_y, pos_z;
    };

    BSStation(); // 构造函数

    void setAntenna(Antenna::transformAngle tr_angle, int P, double zeta, double dBiGain, Antenna::panelarray panel); // 设置基站天线

    void setPosition(Position pos); // 设置基站空间坐标

    void setFrequency(double freq); // 设置基站工作频率

    Antenna m_antenna;       // 基站天线
    Position m_position;     // 基站空间坐标
    double DL_freq, UP_freq; // 基站上下行频率
};