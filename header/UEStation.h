#pragma once
#include "Antenna.h"

class UEStation
{
public:
    struct Position // 空间坐标结构体
    {
        double pos_x, pos_y, pos_z;
    };

    struct moveDirection // 运动方向结构体
    {
        double dir_x, dir_y, dir_z;
    };

    UEStation();

    void setAntenna(Antenna::transformAngle tr_angle, int P, double zeta, double dBiGain, Antenna::panelarray panel); // 设置终端天线

    void setPosition(Position pos); // 设置终端空间坐标

    void setSpeed(moveDirection dir); // 设置终端运动方向

    Antenna m_antenna;   // 终端天线
    Position m_position; // 终端空间坐标
    moveDirection m_dir; // 终端运动方向
    double m_speed;      // 终端速度大小
};