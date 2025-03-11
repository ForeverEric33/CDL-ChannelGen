#include "D:\my_programm\channel_modeling\header\Antenna.h"
// 设置天线（旋转角度，极化方式，极化角度，定向增益，面板排列）
Antenna::Antenna()
{
}
Antenna::Antenna(transformAngle trangle, int P, double zeta, double dBiGain, panelarray panel)
{
    this->m_trangle.alpha = deg2rad(trangle.alpha);
    this->m_trangle.beta = deg2rad(trangle.beta);
    this->m_trangle.gamma = deg2rad(trangle.gamma);
    this->m_P = P;
    this->m_zeta = deg2rad(zeta);
    this->m_dBiGain = dBiGain;
    this->m_panel = panel;

    // 如果是单极化
    if (P == 1)
    {
        this->AntennaNum = panel.M * panel.N;
    }
    // 如果是双极化
    else
    {
        this->AntennaNum = panel.M * panel.N * 2;
    }

    this->panel_zeta = getzeta(m_P, m_zeta, AntennaNum);
}

// 角度转换公式
double Antenna::deg2rad(double angle_deg)
{
    return angle_deg * (PI / 180);
}

// 901 7.1节公式7.1-15
double Antenna::getpsi(transformAngle trangle, double theta, double phi) // 公式7.1-15
{
    theta = theta * PI / 180;
    phi = phi * PI / 180;
    double real = sin(trangle.gamma) * cos(theta) * sin(phi - trangle.alpha) + cos(trangle.gamma) * (cos(trangle.beta) * sin(theta) - sin(trangle.beta) * cos(theta) * cos(theta - trangle.alpha));
    double imag = sin(trangle.gamma) * cos(phi - trangle.alpha) + sin(trangle.beta) * sin(trangle.gamma) * sin(phi - trangle.alpha);
    return atan2(imag, real);
}

// 901 7.1节公式7.1-7
double Antenna::gettheta(transformAngle trangle, double theta, double phi) // 公式7.1-7
{
    theta = theta * PI / 180;
    phi = phi * PI / 180;
    double real = cos(trangle.beta) * cos(trangle.gamma) * cos(theta) + (sin(trangle.beta) * cos(trangle.gamma) * cos(phi - trangle.alpha) - sin(trangle.gamma) * sin(phi - trangle.alpha)) * sin(theta);
    return acos(real);
}

// 901 7.1节7.1-8
double Antenna::getphi(transformAngle trangle, double theta, double phi) // 公式7.1-8
{
    theta = theta * PI / 180;
    phi = phi * PI / 180;
    double real = cos(trangle.beta) * sin(theta) * cos(phi - trangle.alpha) - sin(trangle.beta) * cos(theta);
    double imag = cos(trangle.beta) * sin(trangle.gamma) * cos(theta) + (sin(trangle.beta) * sin(trangle.gamma) * cos(phi - trangle.alpha) + cos(trangle.gamma) * sin(phi - trangle.alpha)) * sin(theta);
    return atan2(imag, real);
}

// 计算天线方向图，901 7.3节表格7.3-1
double Antenna::getAEGain(double theta, double phi, double dBiGain) // Table7.3-1
{
    theta = theta * 180 / PI;
    phi = phi * 180 / PI;
    double gain_theta = -min(12 * pow((theta - 90) / 65, 2), 30.0);
    double gain_phi = -min(12 * pow((phi / 65), 2), 30.0);
    return pow(10.0, 0.1 * (dBiGain - min(-(gain_theta + gain_phi), 30.0)));
}

// 计算极化角度
vector<double> Antenna::getzeta(int P, double zeta, int Num) //
{
    double zeta1, zeta2;
    zeta1 = zeta;
    if (zeta > 0)
    {
        zeta2 = zeta - PI / 2;
    }
    else
    {
        zeta2 = zeta + PI / 2;
    }

    vector<double> ant_zeta;
    if (P == 1)
    {
        ant_zeta.resize(Num);
        for (int n = 0; n < Num; n++)
        {
            ant_zeta[n] = zeta1;
        }
    }
    else
    {
        ant_zeta.resize(Num * 2);
        for (int n = 0; n < Num; n++)
        {
            ant_zeta[2 * n] = zeta1;
            ant_zeta[2 * n + 1] = zeta2;
        }
    }
    return ant_zeta;
}

// 计算与参考天线间隔
vector<Antenna::paneldis> Antenna::getpaneldis(double lambda)
{
    vector<paneldis> paneldis;
    paneldis.resize(AntennaNum);
    for (int n = 0; n < AntennaNum; n++)
    {
        paneldis[n].dis_x = 0;
        if (m_P == 1)
        {
            paneldis[n].dis_y = n % m_panel.N * m_panel.dH * lambda;
            paneldis[n].dis_z = n / m_panel.N * m_panel.dV * lambda;
        }
        else
        {
            paneldis[n].dis_y = (n / 2) % m_panel.N * m_panel.dH * lambda;
            paneldis[n].dis_z = (n / 2) / m_panel.N * m_panel.dV * lambda;
        }
    }
    return paneldis;
}

// 计算天线极化功率 901 7.1节公式7.1-21、7.1-22
Antenna::Field Antenna::getField(double theta, double phi, int n) // 公式77.1-21和公式7.1-22
{
    Field m_Field;
    psi = getpsi(m_trangle, theta, phi);
    theta_lcs = gettheta(m_trangle, theta, phi);
    phi_lcs = getphi(m_trangle, theta, phi);
    AEGain = getAEGain(theta_lcs, phi_lcs, m_dBiGain);

    m_Field.Field_theta = sqrt(AEGain) * (cos(panel_zeta[n]) * cos(psi) - sin(panel_zeta[n]) * sin(psi));
    m_Field.Field_phi = sqrt(AEGain) * (cos(panel_zeta[n]) * sin(psi) + sin(panel_zeta[n]) * cos(psi));

    return m_Field;
}

// 定义均匀线性阵列
void Antenna::ULAantenna(transformAngle trangle, int P, double zeta, double dBiGain, double linedis, int num)
{
    this->m_trangle.alpha = deg2rad(trangle.alpha);
    this->m_trangle.beta = deg2rad(trangle.beta);
    this->m_trangle.gamma = deg2rad(trangle.gamma);
    this->m_P = P;
    this->m_zeta = deg2rad(zeta);
    this->m_dBiGain = dBiGain;
    this->panel_zeta = getzeta(m_P, m_zeta, AntennaNum);
}

// 定义均匀圆形阵列天线
void Antenna::UCAantenna(transformAngle trangle, int P, double zeta, double dBiGain, double radius, int num)
{
    this->m_trangle.alpha = deg2rad(trangle.alpha);
    this->m_trangle.beta = deg2rad(trangle.beta);
    this->m_trangle.gamma = deg2rad(trangle.gamma);
    this->m_P = P;
    this->m_zeta = deg2rad(zeta);
    this->m_dBiGain = dBiGain;
    this->panel_zeta = getzeta(m_P, m_zeta, AntennaNum);
}