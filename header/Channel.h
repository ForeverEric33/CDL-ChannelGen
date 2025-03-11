#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <complex>
#include "BSStation.h"
#include "UEStation.h"
using namespace std;

class Channel
{
public:
    struct Cluster // 簇归一化时延，功率，角度
    {
        double delay, power, aod, aoa, zod, zoa;
    };

    struct ClusterPar // 簇内角度扩展，交叉极化功率比，是否为LoS
    {
        double casd, casa, czsd, czsa;
        int XPR_dB, losflag;
    };

    struct ADangle // 生成水平角和垂直角
    {
        vector<vector<double>> phi_aod, phi_aoa, theta_zoa, theta_zod;
    };

    struct Phase // 初始随即相位
    {
        vector<vector<double>> phase11, phase12, phase21, phase22;
    };

    struct Jitter // 相位扰动
    {
        complex<double> Jittervv, Jittervh, Jitterhv, Jitterhh;
    };

    struct R_rxtx // 极坐标向量转(X,Y,Z)
    {
        double R1, R2, R3;
    };

    struct SSPCongig // 小尺度参数，根据大尺度参数表所得
    {
        double sigma_kf, mu_kf, sigma_asa, mu_asa, sigma_asd, mu_asd, sigma_zsa, mu_zsa, sigma_zsd, mu_zsd, delay_spread;
    };

    struct AS_desired // 角度扩展
    {
        double ASD, ASA, ZSD, ZSA;
    };

    void setChannel(int scenario, string filename, double ds, double maxSpeed, int SD, double time, bool isPathloss, bool isShadowfading); // 设置信道参数，最大终端移动速度、过采样因子，仿真时间，是否启用路损或者阴影

    void writeCIR(vector<BSStation> BS, vector<UEStation> UE, string filename); // 写入信道系数至文件

    void RUN(vector<BSStation> BS, vector<UEStation> UE); // 运行信道仿真功能

    double m_DS; // 时延扩展

    int m_randomseed; // 随机撒点，暂时没用到

protected:
    vector<Cluster> m_cluster;         // 簇
    ClusterPar m_pars;                 // 簇内参数
    vector<double> ray_alpha;          // 子径偏转角
    ADangle m_angle;                   // 水平角和垂直角
    vector<vector<double>> m_kappa;    // 交叉极化功率比(w)
    Phase m_phase;                     // 随即初始相位
    double m_lambda, m_TS;             // 波长，采样间隔
    vector<vector<R_rxtx>> r_rx, r_tx; // 极坐标向量转(X,Y,Z)
    vector<vector<Jitter>> m_jitter;   // 相位干扰
    string m_model_name;               // 模型名称

private:
    void getAngle(vector<Cluster> cluster, ClusterPar pars); // 生成角度

    void getKappa(int clusternum, int raynum, int XPR); // 生成交叉极化功率比

    void getPhase(int clusternun, int raynum); // 生成相位

    void getJitter(vector<vector<double>> kappa, Phase phase); // 生成相位扰动

    void getCIR(BSStation BS, UEStation UE); // 结成链路生成信道系数

    void getLSP(BSStation BS, UEStation UE, int scenerios); // 计算大尺度衰落

    void getSSP(int scenerios, int losflag, double freq, double dis_2D, double h_ut, double h_bs); // 计算小尺度参数

    double m_maxspeed;                                           // 最大终端移动速度
    int m_SD;                                                    // 采样密度
    double running_time;                                         // 仿真时间
    double coeffUpdataRate;                                      // 采样率
    int sampleNum;                                               // 采样点数
    vector<vector<vector<vector<vector<complex<double>>>>>> m_H; // 信道系数，为什么嵌套这么多层
    bool isPL, isSF;                                             // 是否计算路径损耗和阴影衰落
    double m_pathloss = 1;                                       // 路损默认为1
    double m_shadowfading = 1;                                   // 阴影衰落默认为1
    int m_scenario;                                              // 仿真场景

    SSPCongig m_SSP; // 小尺度参数
    double K_factor; // LoS场景下K因子
    AS_desired m_as;
    bool israndom = true;
    int randomseed = 0; // 随机散点
};