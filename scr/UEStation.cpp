#include "D:\my_programm\channel_modeling\header\UEStation.h"

UEStation::UEStation()
{
}

// 设置终端天线
void UEStation::setAntenna(Antenna::transformAngle tr_angle, int P, double zeta, double dBiGain, Antenna::panelarray panel)
{
    Antenna ant(tr_angle, P, zeta, dBiGain, panel);
    this->m_antenna = ant;
}

// 设置终端空间坐标
void UEStation::setPosition(Position pos)
{
    this->m_position = pos;
}

// 设置终端运动方向
void UEStation::setSpeed(moveDirection dir)
{
    this->m_dir = dir;
    this->m_speed = sqrt(pow(m_dir.dir_x, 2) + pow(m_dir.dir_y, 2) + pow(m_dir.dir_z, 2));
}