#include "D:\my_programm\channel_modeling\header\BSStation.h"

BSStation::BSStation()
{
}

// 设置基站天线
void BSStation::setAntenna(Antenna::transformAngle tr_angle, int P, double zeta, double dBiGain, Antenna::panelarray panel)
{
    Antenna ant(tr_angle, P, zeta, dBiGain, panel);
    this->m_antenna = ant;
}

// 设置基站空间坐标
void BSStation::setPosition(Position pos)
{
    this->m_position = pos;
}

// 设置基站工作频率
void BSStation::setFrequency(double freq)
{
    this->DL_freq = freq;
}