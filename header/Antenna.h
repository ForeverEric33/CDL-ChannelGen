#pragma once
#include <iostream>
#include <vector>
#include <math.h>

#define PI acos(-1)
#define LIGHT_SPEED 3e8

using namespace std;

class Antenna
{
public:
    struct transformAngle // 定义坐标旋转角度
    {
        double alpha, beta, gamma;
    };
    struct panelarray // 定义天线面板排列
    {
        int M, N;
        double dH, dV;
    };
    struct Field // 生成天线方向图
    {
        double Field_theta, Field_phi;
    };

    struct paneldis // 定义天线间距分布
    {
        double dis_x, dis_y, dis_z;
    };

    Antenna(); // 构造函数

    Antenna(transformAngle tr_angle, int P, double zeta, double dBiGain, panelarray panel); // 设置天线参数

    void ULAantenna(transformAngle tr_angle, int P, double zeta, double dBiGain, double linedis, int num);

    void UCAantenna(transformAngle tr_angle, int P, double zeta, double dBiGain, double radius, int num);

    double deg2rad(double ange_deg); // 角度单位转换

    double getpsi(transformAngle tr_angle, double theta, double phi); // 901 7.1节公式7.1-15,计算从LCS到GCS方向图转化的角度

    double gettheta(transformAngle tr_angle, double theta, double phi); // 901 7.1节公式7.1-7,计算LCS的theta

    double getphi(transformAngle tr_anggle, double theta, double phi); // 901 7.1节7.1-8，计算LCS的phi

    double getAEGain(double theta, double phi, double dBiGain); // 计算天线方向图，901 7.3节表格7.3-1

    vector<double> getzeta(int P, double zeta, int Num); // 计算天线极化角zeta

    Field getField(double theta, double phi, int n); // 计算天线极化功率 901 7.1节公式7.1-21、7.1-22

    vector<paneldis> getpaneldis(double lambda); // 计算天线阵元间距 901 7.3 Figure7.3-1

    int AntennaNum; // 天线总数

    vector<paneldis> m_paneldis; // 天线阵列间隔

private:
    // 输入
    transformAngle m_trangle; // LCS-GCS坐标旋转角度
    int m_P;                  // 极化方式
    double m_zeta;            // 极化角度
    double m_dBiGain;         // 定向增益
    panelarray m_panel;       // 天线面板，排列方式，排列间隔

    // 中间参数
    double psi;
    double theta_lcs;
    double phi_lcs;
    double AEGain;
    vector<double> panel_zeta;
};
