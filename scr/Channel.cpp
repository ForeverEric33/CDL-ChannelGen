#include "D:\my_programm\channel_modeling\header\Channel.h"
#include <ctime>
#include <string>
vector<double> cale_angular_spreads(vector<Channel::Cluster> cluster);

// 设置信道参数，并读取模型参数
void Channel::setChannel(int scenario, string filename, double ds, double maxSpeed, int SD, double time, bool isPathloss, bool isShadowfading)
{
    this->m_maxspeed = maxSpeed;
    this->m_SD = SD;
    this->running_time = time;
    this->isPL = isPathloss;
    this->isSF = isShadowfading;
    this->m_scenario = scenario;
    this->m_model_name = filename;
    this->m_DS = ds;

    ifstream ifs;
    ifs.open(filename, ios::in);
    if (!ifs.is_open())
    {
        cout << "簇读取失败！" << endl;
        ifs.close();
        return;
    }
    int clusternum, raynum;
    ifs >> clusternum &&ifs >> raynum;

    Cluster clu;

    for (int i = 0; i < clusternum; i++)
    {
        ifs >> clu.delay &&ifs >> clu.power &&ifs >> clu.aod &&ifs >> clu.aoa &&ifs >> clu.zod &&ifs >> clu.zoa;
        clu.power = pow(10, 0.1 * clu.power);
        this->m_cluster.push_back(clu);
    }
    double alpha;
    for (int i = 0; i < raynum; i++)
    {
        ifs >> alpha;
        ray_alpha.push_back(alpha);
    }
    ifs >> m_pars.casd;
    ifs >> m_pars.casa;
    ifs >> m_pars.czsd;
    ifs >> m_pars.czsa;
    ifs >> m_pars.XPR_dB;
    ifs >> m_pars.losflag;

    ifs.close();
}

// 生成水平角和垂直角，901 7.7.1-stpe1 and step2
void Channel::getAngle(vector<Cluster> clusters, ClusterPar par)
{
    // 设置角度的vector容量和簇数一致，有几个簇，有几个角度
    m_angle.phi_aoa.resize(clusters.size());
    m_angle.phi_aod.resize(clusters.size());
    m_angle.theta_zoa.resize(clusters.size());
    m_angle.theta_zod.resize(clusters.size());
    r_rx.resize(clusters.size());
    r_tx.resize(clusters.size());
    // 对于簇内的每个射线，进行同样的初始化操作，ray_alpha为簇内子径偏转角
    for (int i = 0; i < clusters.size(); i++)
    {
        m_angle.phi_aoa[i].resize(ray_alpha.size());
        m_angle.phi_aod[i].resize(ray_alpha.size());
        m_angle.theta_zoa[i].resize(ray_alpha.size());
        m_angle.theta_zod[i].resize(ray_alpha.size());
        r_rx[i].resize(ray_alpha.size());
        r_tx[i].resize(ray_alpha.size());
    }
    // 簇角变为二维数组，两次循环进行运算
    for (int i = 0; i < clusters.size(); i++)
    {
        for (int j = 0; j < ray_alpha.size(); j++)
        {
            m_angle.phi_aoa[i][j] = clusters[i].aoa + par.casa * ray_alpha[j];
            m_angle.phi_aod[i][j] = clusters[i].aod + par.casd * ray_alpha[j];
            m_angle.theta_zoa[i][j] = clusters[i].zoa + par.czsa * ray_alpha[j];
            m_angle.theta_zod[i][j] = clusters[i].zod + par.czsd * ray_alpha[j];
        }
        // 随机排列元素
        random_shuffle(m_angle.phi_aoa[i].begin(), m_angle.phi_aoa[i].end());
        random_shuffle(m_angle.phi_aod[i].begin(), m_angle.phi_aod[i].end());
        random_shuffle(m_angle.theta_zoa[i].begin(), m_angle.theta_zoa[i].end());
        random_shuffle(m_angle.theta_zod[i].begin(), m_angle.theta_zod[i].end());
        for (int j = 0; j < ray_alpha.size(); j++)
        {
            // step11中公式7.5-23、7.5-24，并且将角度制转化为弧度制
            r_rx[i][j].R1 = sin(m_angle.theta_zoa[i][j] * PI / 180) * cos(m_angle.phi_aoa[i][j] * PI / 180);
            r_rx[i][j].R2 = sin(m_angle.theta_zoa[i][j] * PI / 180) * sin(m_angle.phi_aoa[i][j] * PI / 180);
            r_rx[i][j].R3 = cos(m_angle.theta_zoa[i][j] * PI / 180);

            r_tx[i][j].R1 = sin(m_angle.theta_zod[i][j] * PI / 180) * cos(m_angle.phi_aod[i][j] * PI / 180);
            r_tx[i][j].R2 = sin(m_angle.theta_zod[i][j] * PI / 180) * sin(m_angle.phi_aod[i][j] * PI / 180);
            r_tx[i][j].R3 = cos(m_angle.theta_zod[i][j] * PI / 180);
        }
    }
}

// 计算交叉极化功率比(dB-W)，901 7.7.1-step3
void Channel::getKappa(int clusternum, int raynum, int XPR)
{
    m_kappa.resize(clusternum);
    for (int i = 0; i < clusternum; i++)
    {
        m_kappa[i].resize(raynum);
    }

    for (int i = 0; i < clusternum; i++)
    {
        for (int j = 0; j < raynum; j++)
        {
            m_kappa[i][j] = pow(10, 0.1 * XPR);
        }
    }
}

// 生成随即初始相位，在(-pi,pi)服从均匀随即分布,901 7.5-step10
void Channel::getPhase(int clusternum, int raynum)
{
    // 先变成二维数组
    m_phase.phase11.resize(clusternum);
    m_phase.phase12.resize(clusternum);
    m_phase.phase21.resize(clusternum);
    m_phase.phase22.resize(clusternum);
    for (int i = 0; i < clusternum; i++)
    {
        m_phase.phase11[i].resize(raynum);
        m_phase.phase12[i].resize(raynum);
        m_phase.phase21[i].resize(raynum);
        m_phase.phase22[i].resize(raynum);
    }
    // 随机数生成器
    default_random_engine gen(time(NULL));
    // 随机生成-pi到pi的随机数，服从均匀分布
    uniform_real_distribution<double> dis(-PI, PI);
    for (int i = 0; i < clusternum; i++)
    {
        for (int j = 0; j < raynum; j++)
        {
            m_phase.phase11[i][j] = dis(gen);
            m_phase.phase12[i][j] = dis(gen);
            m_phase.phase21[i][j] = dis(gen);
            m_phase.phase22[i][j] = dis(gen);
        }
    }
}

// 计算相位扰动，901 7.5-step11收发端天线方向图中间的矩阵
void Channel::getJitter(vector<vector<double>> kappa, Phase phase)
{
    int clusternum = kappa.size();
    int raynum = kappa[0].size();
    m_jitter.resize(clusternum);
    for (int i = 0; i < clusternum; i++)
    {
        m_jitter[i].resize(raynum);
    }

    for (int i = 0; i < clusternum; i++)
    {
        for (int j = 0; j < raynum; j++)
        {
            // complex(a,b)=a+bj,原本好像没有对kappa乘-1次方，修改一下
            m_jitter[i][j].Jittervv = exp(complex<double>(0, phase.phase11[i][j]));
            m_jitter[i][j].Jittervh = sqrt(pow((kappa[i][j]), -1)) * exp(complex<double>(0, phase.phase12[i][j]));
            m_jitter[i][j].Jitterhv = sqrt(pow((kappa[i][j]), -1)) * exp(complex<double>(0, phase.phase21[i][j]));
            m_jitter[i][j].Jitterhh = exp(complex<double>(0, phase.phase22[i][j]));
        }
    }
}

// 基站终端结成链路，生成信道冲击响应
void Channel::getCIR(BSStation BS, UEStation UE)
{
    vector<vector<vector<vector<complex<double>>>>> H;                   // 第n个基站到第m个终端在第k次采样中的第x根簇的信道冲击响应
    m_lambda = LIGHT_SPEED / BS.DL_freq;                                 // 波长
    coeffUpdataRate = 2 * m_SD * max(m_maxspeed, UE.m_speed) / m_lambda; // 计算采样率
    m_TS = 1 / coeffUpdataRate;                                          // 采样间隔
    sampleNum = running_time / m_TS;                                     // 采样点数

    // 基站终端3D间距
    double dis_3d = sqrt(pow((BS.m_position.pos_x - UE.m_position.pos_x), 2) + pow((BS.m_position.pos_y - UE.m_position.pos_y), 2) + pow((BS.m_position.pos_z - UE.m_position.pos_z), 2));

    vector<Antenna::paneldis> BS_paneldis = BS.m_antenna.getpaneldis(m_lambda); // 基站天线面板排列间隔
    vector<Antenna::paneldis> UE_paneldis = UE.m_antenna.getpaneldis(m_lambda); // 终端天线面板排列间隔

    H.resize(BS.m_antenna.AntennaNum);
    for (int s = 0; s < BS.m_antenna.AntennaNum; s++) // 第一层循环，遍历所有基站
    {
        H[s].resize(UE.m_antenna.AntennaNum);
        for (int u = 0; u < UE.m_antenna.AntennaNum; u++) // 第二层循环，便利所有终端
        {
            H[s][u].resize(sampleNum);

            for (int t = 0; t < sampleNum; t++) // 第三层循环，遍历所有采样点
            {
                H[s][u][t].resize(m_cluster.size());
                for (int i = 0; i < m_cluster.size(); i++) // 第四层循环，遍历所有簇
                {
                    double power = m_cluster[i].power / ray_alpha.size(); // 将簇功率平均分给所有子径
                    complex<double> Hsunt(0, 0);
                    for (int j = 0; j < ray_alpha.size(); j++) // 第五层循环，遍历所有子径
                    {
                        Antenna::Field BS_field = BS.m_antenna.getField(m_angle.theta_zod[i][j], m_angle.phi_aod[i][j], s); // 计算终端天线方向图
                        Antenna::Field UE_field = UE.m_antenna.getField(m_angle.theta_zoa[i][j], m_angle.phi_aoa[i][j], u); // 计算基站天线方向图
                        complex<double> Hsunmt_NLoS(0, 0), Hsunmt_LoS(0, 0);
                        if (m_pars.losflag == 0) // NLoS场景
                        {
                            // 901 7.5-step11 公式7.5-28
                            Hsunmt_NLoS = sqrt(power) * (BS_field.Field_theta * (UE_field.Field_theta * m_jitter[i][j].Jittervv + UE_field.Field_phi * m_jitter[i][j].Jitterhv) + BS_field.Field_phi * (UE_field.Field_theta * m_jitter[i][j].Jittervh + UE_field.Field_phi * m_jitter[i][j].Jitterhh)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE_paneldis[u].dis_x + r_rx[i][j].R2 * UE_paneldis[u].dis_y + r_rx[i][j].R3 * UE_paneldis[u].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_tx[i][j].R1 * BS_paneldis[s].dis_x + r_tx[i][j].R2 * BS_paneldis[s].dis_y + r_tx[i][j].R3 * BS_paneldis[s].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE.m_dir.dir_x + r_rx[i][j].R2 * UE.m_dir.dir_y + r_rx[i][j].R3 * UE.m_dir.dir_z) * m_TS * t / m_lambda));
                            Hsunt += Hsunmt_NLoS;
                        }
                        else // LoS场景
                        {
                            if (i == 0) // LoS径
                            {
                                // 901 7.5-step11 公式7.5-29
                                Hsunmt_LoS = (1 / 20) * (BS_field.Field_theta * (UE_field.Field_theta * 1 + UE_field.Field_phi * 0) + BS_field.Field_phi * (UE_field.Field_theta * 0 + UE_field.Field_phi * -1)) * exp(complex<double>(0, -2 * PI * dis_3d / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE_paneldis[u].dis_x + r_rx[i][j].R2 * UE_paneldis[u].dis_y + r_rx[i][j].R3 * UE_paneldis[u].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_tx[i][j].R1 * BS_paneldis[s].dis_x + r_tx[i][j].R2 * BS_paneldis[s].dis_y + r_tx[i][j].R3 * BS_paneldis[s].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE.m_dir.dir_x + r_rx[i][j].R2 * UE.m_dir.dir_y + r_rx[i][j].R3 * UE.m_dir.dir_z) * m_TS * t / m_lambda));
                            }
                            else // NLoS径
                            {
                                // 901 7.5-step11 公式7.5-28
                                Hsunmt_NLoS = sqrt(power) * (BS_field.Field_theta * (UE_field.Field_theta * m_jitter[i][j].Jittervv + UE_field.Field_phi * m_jitter[i][j].Jitterhv) + BS_field.Field_phi * (UE_field.Field_theta * m_jitter[i][j].Jittervh + UE_field.Field_phi * m_jitter[i][j].Jitterhh)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE_paneldis[u].dis_x + r_rx[i][j].R2 * UE_paneldis[u].dis_y + r_rx[i][j].R3 * UE_paneldis[u].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_tx[i][j].R1 * BS_paneldis[s].dis_x + r_tx[i][j].R2 * BS_paneldis[s].dis_y + r_tx[i][j].R3 * BS_paneldis[s].dis_z) / m_lambda)) * exp(complex<double>(0, 2 * PI * (r_rx[i][j].R1 * UE.m_dir.dir_x + r_rx[i][j].R2 * UE.m_dir.dir_y + r_rx[i][j].R3 * UE.m_dir.dir_z) * m_TS * t / m_lambda));
                            }
                            Hsunt += (sqrt(1 / (K_factor + 1)) * Hsunmt_NLoS + sqrt(K_factor / (K_factor + 1)) * Hsunmt_LoS);
                        }
                    }
                    H[s][u][t][i] = sqrt(m_pathloss * m_shadowfading) * Hsunt; // 901 7.5-step12
                }
            }
        }
    }
    m_H.push_back(H);
}

void Channel::writeCIR(vector<BSStation> BS, vector<UEStation> UE, string filename)
{

    ofstream ofs;
    ofs.open(filename, ios::out);
    int bsNum = BS.size();
    ofs << "基站数量为：" << bsNum << " ";  // 基站数量
    int utNum = UE.size();                  //
    ofs << "终端数量为：" << utNum << endl; // 终端数量
    for (int i = 0; i < bsNum; i++)
    {
        ofs << "第" << i + 1 << "个基站天线数量为：" << BS[i].m_antenna.AntennaNum << " "; // 基站天线数
    }
    for (int i = 0; i < utNum; i++)
    {
        ofs << "第" << i + 1 << "个终端天线数量为：" << UE[i].m_antenna.AntennaNum << endl; // 终端天线数
    }

    int LinkNum = bsNum * utNum;
    ofs << "链路数：" << LinkNum << endl; // 链路数
    ofs << "采样密度：" << m_SD << endl;  // 采样密度
    int linkflag = 0;
    for (int i = 0; i < BS.size(); i++)
    {
        for (int j = 0; j < UE.size(); j++)
        {

            ofs << "基站" << i + 1 << " ";  // 基站编号
            ofs << "终端" << j + 1 << endl; // 终端编号

            ofs << "中心频率为：" << BS[i].DL_freq << endl; // 中心频率
            ofs << "采样率为：" << coeffUpdataRate << endl; // 系统更新率
            int clusterNum = m_cluster.size();
            ofs << "簇的数量：" << clusterNum << endl;
            ofs << "采样次数：" << sampleNum << endl;
            // 取实部
            double IData = 0;
            double QData = 0;
            double Delay = 0;
            for (int s = 0; s < BS[i].m_antenna.AntennaNum; s++)
            {
                for (int u = 0; u < UE[j].m_antenna.AntennaNum; u++)
                {
                    // for (int n = 0; n < m_cluster.size(); n++)
                    for (int t = 0; t < sampleNum; t++)
                    {
                        // for (int t = 0; t < sampleNum; t++)
                        for (int n = 0; n < m_cluster.size(); n++)
                        {
                            IData = m_H[linkflag][s][u][t][n].real();
                            QData = m_H[linkflag][s][u][t][n].real();
                            Delay = m_cluster[n].delay;
                            ofs << "第" << s + 1 << "个基站天线"
                                << "到第" << u + 1 << "个终端天线"
                                << "在第" << t + 1 << "次"
                                << "采样"
                                << "中的第" << n + 1 << "个簇的信道冲击响应为：" << IData << "+" << QData << "j"
                                << " ";
                            ofs << "时延为：" << Delay << endl;
                        }
                    }
                }
            }
            linkflag++;
            cout << " 第 " << linkflag << " 条链路写入完毕！" << endl;
        }
    }
    ofs.close();
}

// 运行函数，导入基站终端容器（包含所有参与仿真的基站终端）
void Channel::RUN(vector<BSStation> BS, vector<UEStation> UE)
{
    int linkflag = 1;
    for (int i = 0; i < BS.size(); i++) // 第一层循环，遍历所有基站
    {
        for (int j = 0; j < UE.size(); j++) // 第二层循环，遍历所有终端
        {
            // 基站终端2D距离
            double dis_2D = sqrt(pow((BS[i].m_position.pos_x - UE[j].m_position.pos_x), 2) + pow((BS[i].m_position.pos_y - UE[j].m_position.pos_y), 2));

            this->getLSP(BS[i], UE[j], m_scenario);                                                                                // 计算大尺度衰减，包括路径损耗和阴影衰落
            this->getSSP(m_scenario, m_pars.losflag, BS[i].DL_freq / 1e9, dis_2D, UE[j].m_position.pos_z, BS[i].m_position.pos_z); // 生成小尺度参数，并进行功率缩放，时延扩展，角度扩展
            this->getAngle(m_cluster, m_pars);                                                                                     // 生成水平角和垂直角
            this->getKappa(m_cluster.size(), ray_alpha.size(), m_pars.XPR_dB);                                                     // 计算交叉极化功率比
            this->getPhase(m_cluster.size(), ray_alpha.size());                                                                    // 生成随即初始相位
            this->getJitter(m_kappa, m_phase);                                                                                     // 生成相位扰动
            this->getCIR(BS[i], UE[j]);                                                                                            // 生成信道冲击响应
            cout << " 第 " << linkflag << " 条链路仿真完毕！" << endl;
            linkflag++;
        }
    }
}

// 计算大尺度衰减，包括路径损耗和阴影衰落
void Channel::getLSP(BSStation BS, UEStation UE, int scenarios)
{
    // pow（x，y）对x求y次幂操作
    double dis_3d = sqrt(pow((BS.m_position.pos_x - UE.m_position.pos_x), 2) + pow((BS.m_position.pos_y - UE.m_position.pos_y), 2) + pow((BS.m_position.pos_z - UE.m_position.pos_z), 2));
    double dis_2d = sqrt(pow((BS.m_position.pos_x - UE.m_position.pos_x), 2) + pow((BS.m_position.pos_y - UE.m_position.pos_y), 2));
    double dis_bp1 = 2 * PI * BS.m_position.pos_z * UE.m_position.pos_z * BS.DL_freq / 3e8;
    double dis_bp2 = 4 * (BS.m_position.pos_z - 1) * (UE.m_position.pos_z - 1) * BS.DL_freq / 3e8;
    double pathloss1, pathloss2, pathloss3;
    double pathloss;
    double sigma;
    switch (scenarios)
    {
    // 先判断是否为los径，然后选取场景，计算pathloss，获取sigma
    case 1: // UMA
        if (m_pars.losflag == 1)
        {
            if (dis_2d >= 10 && dis_2d <= dis_bp2)
            {
                pathloss = 28.0 * 22 * log10(dis_3d) + 20 * log10(BS.DL_freq);
            }
            else
            {
                pathloss = 28.0 * 40 * log10(dis_3d) + 20 * log10(BS.DL_freq) - 9 * log10(pow(dis_bp2, 2) + pow(BS.m_position.pos_z - UE.m_position.pos_z, 2));
            }
            sigma = 4;
        }
        else
        {
            pathloss = 32.4 + 20 * log10(BS.DL_freq) + 30 * log10(dis_3d);
            sigma = 7.8;
        }
        break;
    case 2: // UMI
        if (m_pars.losflag == 1)
        {
            if (dis_2d >= 10 && dis_2d <= dis_bp2)
            {
                pathloss = 32.4 + 21 * log10(dis_3d) + 20 * log10(BS.DL_freq);
            }
            else
            {
                pathloss = 32.4 + 40 * log10(dis_3d) + 20 * log10(BS.DL_freq) - 9.5 * log10(pow(dis_bp2, 2) + pow(BS.m_position.pos_z - UE.m_position.pos_z, 2));
            }
            sigma = 4;
        }
        else
        {
            pathloss = 32.4 + 20 * log10(BS.DL_freq) + 32.9 * log10(dis_3d);
            sigma = 7.82;
        }
        break;
    case 3: // RMA
        if (m_pars.losflag == 1)
        {
            if (dis_2d >= 10 && dis_2d <= dis_bp1)
            {
                pathloss = 20 * log10(40 * PI * dis_3d * BS.DL_freq / 3) + min(0.03 * pow(5, 1.72), 10.0) * log10(dis_3d) - min(0.044 * pow(5, 1.72), 14.77) + 0.002 * log10(5) * dis_3d;
                sigma = 4;
            }
            else
            {
                pathloss = 20 * log10(40 * PI * dis_bp1 * BS.DL_freq / 3) + min(0.03 * pow(5, 1.72), 10.0) * log10(dis_bp1) - min(0.044 * pow(5, 1.72), 14.77) + 0.002 * log10(5) * dis_bp1 + 40 * log10(dis_3d / dis_bp1);
                sigma = 6;
            }
        }
        else
        {

            if (dis_2d >= 10 && dis_2d <= dis_bp1)
            {
                pathloss1 = 20 * log10(40 * PI * dis_3d * BS.DL_freq / 3) + min(0.03 * pow(5, 1.72), 10.0) * log10(dis_3d) - min(0.044 * pow(5, 1.72), 14.77) + 0.002 * log10(5) * dis_3d;
            }
            else
            {
                pathloss1 = 20 * log10(40 * PI * dis_bp1 * BS.DL_freq / 3) + min(0.03 * pow(5, 1.72), 10.0) * log10(dis_bp1) - min(0.044 * pow(5, 1.72), 14.77) + 0.002 * log10(5) * dis_bp1 + 40 * log10(dis_3d / dis_bp1);
            }
            pathloss2 = 161.04 - 7.1 * log10(20) + 7.5 * log10(5) - (24.37 - 3.7 * pow(5 / UE.m_position.pos_z, 2)) * log10(BS.m_position.pos_z) + (43.42 - 3.1 * log10(BS.m_position.pos_z)) * (log10(dis_3d) - 3) + 20 * log10(BS.DL_freq) - (3.2 * pow((log10(11.75 * UE.m_position.pos_z)), 2) - 4.97);
            pathloss = max(pathloss1, pathloss2);
            sigma = 8;
        }
        break;
    case 4: // INF
        pathloss = 31.84 + 21.50 * log10(dis_3d) + 19.00 * log10(BS.DL_freq);
        sigma = 4;
        break;
    case 5: // INF_SL
        pathloss1 = 31.84 + 21.50 * log10(dis_3d) + 19.00 * log10(BS.DL_freq);
        pathloss2 = 33 + 25.5 * log10(dis_3d) + 20 * log10(BS.DL_freq);
        pathloss = max(pathloss1, pathloss2);
        sigma = 5.7;
        break;
    case 6: // INF_DL
        pathloss1 = 31.84 + 21.50 * log10(dis_3d) + 19.00 * log10(BS.DL_freq);
        pathloss2 = 33 + 25.5 * log10(dis_3d) + 20 * log10(BS.DL_freq);
        pathloss3 = 18.6 + 35.7 * log10(dis_3d) + 20 * log(BS.DL_freq);
        // pathloss = max(pathloss1, pathloss2, pathloss3);
        pathloss = max(pathloss1, pathloss2);
        sigma = 7.2;
        break;
    case 7: // INF_SH
        pathloss1 = 31.84 + 21.50 * log10(dis_3d) + 19.00 * log10(BS.DL_freq);
        pathloss2 = 32.4 + 23.0 * log10(dis_3d) + 20 * log10(BS.DL_freq);
        pathloss = max(pathloss1, pathloss2);
        sigma = 5.9;
        break;
    case 8: // INF_DH
        pathloss1 = 31.84 + 21.50 * log10(dis_3d) + 19.00 * log10(BS.DL_freq);
        pathloss2 = 33.63 + 21.9 * log10(dis_3d) + 20 * log10(BS.DL_freq);
        pathloss = max(pathloss1, pathloss2);
        sigma = 4.0;
    case 9: // INH
        if (m_pars.losflag == 1)
        {
            pathloss = 32.4 + 17.3 * log10(dis_3d) + 20 * log10(BS.DL_freq);
            sigma = 3;
        }
        else
        {
            pathloss = 32.4 + 20 * log10(BS.DL_freq) + 31.9 * log10(dis_3d);
            sigma = 8.29;
        }
        break;
    default:
        pathloss = 32.4 + 20 * log10(BS.DL_freq) + 20 * log10(dis_3d);
        sigma = 4;
        break;
    }
    if (isPL)
    {
        m_pathloss = pathloss;
    }
    // 为什么还要做一个对sigma乘（0,1）之间随机数的操作
    if (isSF)
    {
        default_random_engine gen(time(NULL));
        normal_distribution<double> dis(0, 1);
        m_shadowfading = sigma * dis(gen);
    }
}

// 时延扩展计算函数，RMS
double cale_delay_spread(vector<Channel::Cluster> cluster)
{
    double total_power1 = 0.0;
    for (int i = 0; i < cluster.size(); i++)
    {
        total_power1 += cluster[i].power;
    }
    double mean_delay = 0.0;
    for (int i = 0; i < cluster.size(); i++)
    {
        cluster[i].power = cluster[i].power / total_power1; // 每一个簇的功率占总功率的比例
        mean_delay += cluster[i].power * cluster[i].delay;  //
    }
    vector<double> tmp;
    tmp.resize(cluster.size());
    for (int i = 0; i < cluster.size(); i++)
    {
        tmp[i] = cluster[i].delay - mean_delay;
    }
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < cluster.size(); i++)
    {
        sum1 += cluster[i].power * pow(tmp[i], 2);
        sum2 += cluster[i].power * tmp[i];
    }
    return sqrt(sum1 - pow(sum2, 2));
}

// 角度扩展计算函数，RMS DX=sqrt（E(X^2)-(EX)^2）
vector<double> cale_angular_spreads(vector<Channel::Cluster> cluster)
{
    double pt = 0.0; // 总功率

    // 计算总功率
    for (int i = 0; i < cluster.size(); i++)
    {
        pt += cluster[i].power;
    }

    // 归一化每一个簇的功率
    for (int i = 0; i < cluster.size(); i++)
    {
        cluster[i].power = cluster[i].power / pt;
    }

    // 计算复功率
    complex<double> pow_aod(0, 0), pow_aoa(0, 0), pow_eod(0, 0), pow_eoa(0, 0);
    for (int i = 0; i < cluster.size(); i++)
    {
        pow_aod += cluster[i].power * exp(complex<double>(0, cluster[i].aod));
        pow_aoa += cluster[i].power * exp(complex<double>(0, cluster[i].aoa));
        pow_eod += cluster[i].power * exp(complex<double>(0, cluster[i].zod));
        pow_eoa += cluster[i].power * exp(complex<double>(0, cluster[i].zoa));
    }

    double mean_aod, mean_aoa, mean_eod, mean_eoa;

    // 计算平均角度
    mean_aod = atan2(pow_aod.imag(), pow_aod.real());
    mean_aoa = atan2(pow_aoa.imag(), pow_aoa.real());
    mean_eod = atan2(pow_eod.imag(), pow_eod.real());
    mean_eoa = atan2(pow_eoa.imag(), pow_eoa.real());

    // 调整角度
    for (int i = 0; i < cluster.size(); i++)
    {
        cluster[i].aod = cluster[i].aod - mean_aod;
        cluster[i].aoa = cluster[i].aoa - mean_aoa;
        cluster[i].zod = cluster[i].zod - mean_eod;
        cluster[i].zoa = cluster[i].zoa - mean_eoa;

        // 将角度调整到[-π, π]范围内
        cluster[i].aod = atan2(exp(complex<double>(0, cluster[i].aod)).imag(), exp(complex<double>(0, cluster[i].aod)).real());
        cluster[i].aoa = atan2(exp(complex<double>(0, cluster[i].aoa)).imag(), exp(complex<double>(0, cluster[i].aoa)).real());
        cluster[i].zod = atan2(exp(complex<double>(0, cluster[i].zod)).imag(), exp(complex<double>(0, cluster[i].zod)).real());
        cluster[i].zoa = atan2(exp(complex<double>(0, cluster[i].zoa)).imag(), exp(complex<double>(0, cluster[i].zoa)).real());
    }

    double sum1_aod = 0, sum1_aoa = 0, sum1_eod = 0, sum1_eoa = 0, sum2_aod = 0, sum2_aoa = 0, sum2_eod = 0, sum2_eoa = 0, ds_aod = 0, ds_aoa = 0, ds_eod = 0, ds_eoa = 0;
    // 计算各种和值
    for (int i = 0; i < cluster.size(); i++)
    {
        sum1_aod += cluster[i].power * pow(cluster[i].aod, 2);
        sum1_aoa += cluster[i].power * pow(cluster[i].aoa, 2);
        sum1_eod += cluster[i].power * pow(cluster[i].zod, 2);
        sum1_eoa += cluster[i].power * pow(cluster[i].zoa, 2);
        sum2_aod += pow(cluster[i].power * cluster[i].aod, 2);
        sum2_aoa += pow(cluster[i].power * cluster[i].aoa, 2);
        sum2_eod += pow(cluster[i].power * cluster[i].zod, 2);
        sum2_eoa += pow(cluster[i].power * cluster[i].zoa, 2);
    }
    // 计算标准差
    ds_aod = sqrt(sum1_aod - sum2_aod);
    ds_aoa = sqrt(sum1_aoa - sum2_aoa);
    ds_eod = sqrt(sum1_eod - sum2_eod);
    ds_eoa = sqrt(sum1_eoa - sum2_eoa);
    vector<double> ans;
    ans.push_back(mean_aod);
    ans.push_back(mean_aoa);
    ans.push_back(mean_eod);
    ans.push_back(mean_eoa);
    ans.push_back(ds_aod);
    ans.push_back(ds_aoa);
    ans.push_back(ds_eod);
    ans.push_back(ds_eoa);
    return ans;
}

// 计算小尺度参数，对功率缩放，时延扩展，角度扩展
void Channel::getSSP(int scenario, int losflag, double freq, double dis_2D, double h_ut, double h_bs)
{
    // uma
    if (scenario == 1)
    {
        if (losflag == 1)
        {
            m_SSP.mu_asd = 1.06 + 0.1114 * log10(freq);
            m_SSP.sigma_asd = 0.28;
            m_SSP.mu_asa = 1.81;
            m_SSP.sigma_asa = 0.2;
            m_SSP.mu_zsa = 0.95;
            m_SSP.sigma_zsa = 0.16;
            m_SSP.mu_zsd = max(-0.5, -2.1 * (dis_2D / 100) - 0.01 * (h_ut - 1.5) + 0.75);
            m_SSP.sigma_zsd = 0.40;
            m_SSP.mu_kf = 9;
            m_SSP.sigma_kf = 3.5;
        }
        else
        {
            m_SSP.mu_asd = 1.5 - 0.1144 * log10(freq);
            m_SSP.sigma_asd = 0.28;
            m_SSP.mu_asa = 2.08 - 0.27 * log10(freq);
            m_SSP.sigma_asa = 0.11;
            m_SSP.mu_zsa = -0.3236 * log10(freq) + 1.512;
            m_SSP.sigma_zsa = 0.16;
            m_SSP.mu_zsd = max(-0.5, -2.1 * (dis_2D / 1000) - 0.01 * (h_ut - 1.5) + 0.9);
            m_SSP.sigma_zsd = 0.49;
        }

        if (freq < 4)
        {
            m_SSP.delay_spread = 363;
        }
        else if (freq < 11 && freq > 4)
        {
            m_SSP.delay_spread = 363;
        }
        else if (freq < 22 && freq >= 11)
        {
            m_SSP.delay_spread = 302;
        }
        else if (freq < 34 && freq >= 22)
        {
            m_SSP.delay_spread = 266;
        }
        else if (freq < 50 && freq >= 34)
        {
            m_SSP.delay_spread = 249;
        }
        else if (freq < 65 && freq >= 50)
        {
            m_SSP.delay_spread = 228;
        }
        else
        {
            m_SSP.delay_spread = 221;
        }
    }
    else if (scenario == 2)
    {
        if (losflag == 1)
        {
            m_SSP.mu_asd = -0.05 * log10(1 + freq) + 1.21;
            m_SSP.sigma_asd = 0.41;
            m_SSP.mu_asa = -0.08 * log10(1 + freq) + 1.73;
            m_SSP.sigma_asa = 0.014 * log10(1 + freq) + 0.28;
            m_SSP.mu_zsa = -0.1 * log10(1 + freq) + 0.73;
            m_SSP.sigma_zsa = -0.04 * log10(1 + freq) + 0.34;
            m_SSP.mu_zsd = max(-0.21, -14.8 * (dis_2D / 1000) + 0.01 * (h_ut - h_bs) + 0.83);
            m_SSP.sigma_zsd = 0.35;
            m_SSP.mu_kf = 9;
            m_SSP.sigma_kf = 5;
        }
        else
        {
            m_SSP.mu_asd = -0.23 * log10(1 + freq) + 1.53;
            m_SSP.sigma_asd = 0.11 * log10(1 + freq) + 0.33;
            m_SSP.mu_asa = -0.08 * log10(1 + freq) + 1.81;
            m_SSP.sigma_asa = 0.05 * log10(1 + freq) + 0.3;
            m_SSP.mu_zsa = -0.04 * log10(1 + freq) + 0.92;
            m_SSP.sigma_zsa = -0.07 * log10(1 + freq) + 0.41;
            m_SSP.mu_zsd = max(-0.5, -3.1 * (dis_2D / 1000) + 0.01 * max(h_ut - h_bs, 0.0) + 0.2);
            m_SSP.sigma_zsd = 0.35;
        }

        if (freq < 4e9)
        {
            m_SSP.delay_spread = 129;
        }
        else if (freq < 11e9 && freq > 4e9)
        {
            m_SSP.delay_spread = 93;
        }
        else if (freq < 22e9 && freq >= 11e9)
        {
            m_SSP.delay_spread = 76;
        }
        else if (freq < 34e9 && freq >= 22e9)
        {
            m_SSP.delay_spread = 66;
        }
        else if (freq < 50e9 && freq >= 34e9)
        {
            m_SSP.delay_spread = 61;
        }
        else if (freq < 65 && freq >= 50e9)
        {
            m_SSP.delay_spread = 55;
        }
        else
        {
            m_SSP.delay_spread = 53;
        }
    }
    else if (scenario == 3)
    {
        if (losflag == 1)
        {
            m_SSP.mu_asd = 0.90;
            m_SSP.sigma_asd = 0.38;
            m_SSP.mu_asa = 1.52;
            m_SSP.sigma_asa = 0.24;
            m_SSP.mu_zsa = 0.47;
            m_SSP.sigma_zsa = 0.40;
            m_SSP.mu_zsd = max(-1.0, -0.17 * (dis_2D / 1000) - 0.01 * (h_ut - 1.5) + 0.22);
            m_SSP.sigma_zsd = 0.34;
            m_SSP.mu_kf = 7;
            m_SSP.sigma_kf = 4;
        }
        else
        {
            m_SSP.mu_asd = 0.95;
            m_SSP.sigma_asd = 0.45;
            m_SSP.mu_asa = 1.52;
            m_SSP.sigma_asa = 0.13;
            m_SSP.mu_zsa = 0.58;
            m_SSP.sigma_zsa = 0.37;
            m_SSP.mu_zsd = max(-1.0, -0.19 * (dis_2D / 1000) - 0.01 * (h_ut - 1.5) + 0.28);
            m_SSP.sigma_zsd = 0.30;
        }

        m_SSP.delay_spread = 37;
    }
    else if (scenario == 4 || scenario == 5 || scenario == 6 || scenario || 7 || scenario == 8)
    {
        if (losflag == 1)
        {
            m_SSP.mu_asd = 1.56;
            m_SSP.sigma_asd = 0.25;
            m_SSP.mu_asa = -0.18 * log10(1 + freq) + 1.78;
            m_SSP.sigma_asa = 0.12 * log10(1 + freq) + 0.2;
            m_SSP.mu_zsa = -0.2 * log10(1 + freq) + 1.5;
            m_SSP.sigma_zsa = 0.35;
            m_SSP.mu_zsd = 1.35;
            m_SSP.sigma_zsd = 0.35;
            m_SSP.mu_kf = 7;
            m_SSP.sigma_kf = 8;
        }
        else
        {
            m_SSP.mu_asd = 1.57;
            m_SSP.sigma_asd = 0.2;
            m_SSP.mu_asa = 1.72;
            m_SSP.sigma_asa = 0.3;
            m_SSP.mu_zsa = -0.13 * log10(1 + freq) + 1.45;
            m_SSP.sigma_zsa = 0.45;
            m_SSP.mu_zsd = 1.2;
            m_SSP.sigma_zsd = 0.55;
        }

        if (freq < 4e9)
        {
            m_SSP.delay_spread = 39;
        }
        else if (freq < 11e9 && freq > 4e9)
        {
            m_SSP.delay_spread = 30;
        }
        else if (freq < 22e9 && freq >= 11e9)
        {
            m_SSP.delay_spread = 24;
        }
        else if (freq < 34e9 && freq >= 22e9)
        {
            m_SSP.delay_spread = 20;
        }
        else if (freq < 50e9 && freq >= 34e9)
        {
            m_SSP.delay_spread = 18;
        }
        else
        {
            m_SSP.delay_spread = 16;
        }
    }
    else if (scenario == 9)
    {
        if (losflag == 1)
        {
            m_SSP.mu_asd = 1.60;
            m_SSP.sigma_asd = 0.18;
            m_SSP.mu_asa = -0.19 * log10(1 + freq) + 1.781;
            m_SSP.sigma_asa = 0.12 * log10(1 + freq) + 0.119;
            m_SSP.mu_zsa = -0.26 * log10(1 + freq) + 1.44;
            m_SSP.sigma_zsa = -0.04 * log10(1 + freq) + 0.264;
            m_SSP.mu_zsd = -1.43 * log10(1 + freq) + 2.228;
            m_SSP.sigma_zsd = 0.13 * log10(1 + freq) + 0.30;
            m_SSP.mu_kf = 7;
            m_SSP.sigma_kf = 4;
        }
        else
        {
            m_SSP.mu_asd = 1.62;
            m_SSP.sigma_asd = 0.25;
            m_SSP.mu_asa = -0.11 * log10(1 + freq) + 1.863;
            m_SSP.sigma_asa = 0.12 * log10(1 + freq) + 0.059;
            m_SSP.mu_zsa = -0.15 * log10(1 + freq) + 1.387;
            m_SSP.sigma_zsa = -0.09 * log10(1 + freq) + 0.746;
            m_SSP.mu_zsd = 1.08;
            m_SSP.sigma_zsd = 0.36;
        }

        if (freq < 4e9)
        {
            m_SSP.delay_spread = 39;
        }
        else if (freq < 11e9 && freq > 4e9)
        {
            m_SSP.delay_spread = 30;
        }
        else if (freq < 22e9 && freq >= 11e9)
        {
            m_SSP.delay_spread = 24;
        }
        else if (freq < 34e9 && freq >= 22e9)
        {
            m_SSP.delay_spread = 20;
        }
        else if (freq < 50e9 && freq >= 34e9)
        {
            m_SSP.delay_spread = 18;
        }
        else
        {
            m_SSP.delay_spread = 16;
        }
    }
    if (m_DS == -1)
    {
        this->m_DS = m_SSP.delay_spread;
    }
    double ASD_desried_in, ASA_desried_in, ZSD_desried_in, ZSA_desried_in;
    double ASD_desried, ASA_desried, ZSD_desried, ZSA_desried;
    double ASD_desried_out, ASA_desried_out, ZSD_desried_out, ZSA_desried_out;
    // 判断角度是否在范围内
    default_random_engine gen(time(NULL));
    normal_distribution<double> dis(0, 1);
    ASD_desried = min(pow(10, m_SSP.sigma_asd * dis(gen) + m_SSP.mu_asd), 104.0);
    ASA_desried = min(pow(10, m_SSP.sigma_asa * dis(gen) + m_SSP.mu_asa), 104.0);
    ZSD_desried = min(pow(10, m_SSP.sigma_zsd * dis(gen) + m_SSP.mu_zsd), 52.0);
    ZSA_desried = min(pow(10, m_SSP.sigma_zsa * dis(gen) + m_SSP.mu_zsa), 52.0);
    // 转换为弧度
    ASD_desried = ASD_desried * PI / 180;
    ASA_desried = ASA_desried * PI / 180;
    ZSD_desried = ZSD_desried * PI / 180;
    ZSA_desried = ZSA_desried * PI / 180;

    if (m_pars.losflag == 1)
    {

        default_random_engine gen(time(NULL));
        normal_distribution<double> dis(0, 1);
        double K_factor = m_SSP.sigma_kf * dis(gen) + m_SSP.mu_kf;

        double total_power1 = 0.0, total_power2 = 0.0;
        for (int i = 0; i < m_cluster.size(); i++)
        {
            total_power1 += m_cluster[i].power;
        }
        m_cluster[0].power = (total_power1 - m_cluster[0].power) * pow(10.0, 0.1 * K_factor);
        for (int i = 0; i < m_cluster.size(); i++)
        {
            total_power2 += m_cluster[i].power;
        }

        for (int i = 0; i < m_cluster.size(); i++)
        {
            m_cluster[i].power = m_cluster[i].power * (total_power1 / total_power2);
        }
        //////delay_spread for QuaDriGa
        for (int i = 0; i < m_cluster.size(); i++)
        {
            m_cluster[i].delay = m_cluster[i].delay / cale_delay_spread(m_cluster);
        }
    }
    double ds111 = cale_delay_spread(m_cluster);
    for (int i = 0; i < m_cluster.size(); i++)
    {
        m_cluster[i].delay = m_DS * m_cluster[i].delay / ds111;
    }

    for (int i = 0; i < m_cluster.size(); i++)
    {
        m_cluster[i].aod = m_cluster[i].aod * (PI / 180);
        m_cluster[i].aoa = m_cluster[i].aoa * (PI / 180);
        m_cluster[i].zod = (90 - m_cluster[i].zod) * (PI / 180);
        m_cluster[i].zoa = (90 - m_cluster[i].zoa) * (PI / 180);
    }
    vector<double> spread_angle = cale_angular_spreads(m_cluster); // 前四个是mean_angle.后四个是as_model

    vector<double> tmp_aod, tmp_aoa, tmp_zod, tmp_zoa;
    tmp_aod.resize(m_cluster.size());
    tmp_aoa.resize(m_cluster.size());
    tmp_zod.resize(m_cluster.size());
    tmp_zoa.resize(m_cluster.size());
    for (int i = 0; i < m_cluster.size(); i++)
    {
        tmp_aod[i] = atan2(exp(complex<double>(0.0, m_cluster[i].aod - spread_angle[0])).imag(), exp(complex<double>(0.0, m_cluster[i].aod - spread_angle[0])).real()) * (ASD_desried / spread_angle[4]);
        tmp_aoa[i] = atan2(exp(complex<double>(0.0, m_cluster[i].aoa - spread_angle[1])).imag(), exp(complex<double>(0.0, m_cluster[i].aoa - spread_angle[1])).real()) * (ASA_desried / spread_angle[5]);
        tmp_zod[i] = atan2(exp(complex<double>(0.0, m_cluster[i].zod - spread_angle[2])).imag(), exp(complex<double>(0.0, m_cluster[i].zod - spread_angle[2])).real()) * (ZSD_desried / spread_angle[6]);
        tmp_zoa[i] = atan2(exp(complex<double>(0.0, m_cluster[i].zoa - spread_angle[3])).imag(), exp(complex<double>(0.0, m_cluster[i].zoa - spread_angle[3])).real()) * (ZSA_desried / spread_angle[7]);
    }
    if (m_pars.losflag == 1)
    {
        double rot_angle_aod, rot_angle_aoa, rot_angle_zod, rot_angle_zoa;
        rot_angle_aod = m_cluster[0].aod - tmp_aod[0];
        rot_angle_aoa = m_cluster[0].aoa - tmp_aoa[0];
        rot_angle_zod = m_cluster[0].zod - tmp_zod[0];
        rot_angle_zoa = m_cluster[0].zoa - tmp_zoa[0];
        for (int i = 0; i < m_cluster.size(); i++)
        {
            tmp_aod[i] = atan2(exp(complex<double>(0, tmp_aod[i] + rot_angle_aod)).imag(), exp(complex<double>(0, tmp_aod[i] + rot_angle_aod)).real());
            tmp_aoa[i] = atan2(exp(complex<double>(0, tmp_aoa[i] + rot_angle_aoa)).imag(), exp(complex<double>(0, tmp_aoa[i] + rot_angle_aoa)).real());
            tmp_zod[i] = atan2(exp(complex<double>(0, tmp_zod[i] + rot_angle_zod)).imag(), exp(complex<double>(0, tmp_zod[i] + rot_angle_zod)).real());
            tmp_zoa[i] = atan2(exp(complex<double>(0, tmp_zoa[i] + rot_angle_zoa)).imag(), exp(complex<double>(0, tmp_zoa[i] + rot_angle_zoa)).real());
        }
    }
    else
    {
        for (int i = 0; i < m_cluster.size(); i++)
        {
            tmp_aod[i] = atan2(exp(complex<double>(0, tmp_aod[i] + spread_angle[0])).imag(), exp(complex<double>(0, tmp_aod[i] + spread_angle[0])).real());
            tmp_aoa[i] = atan2(exp(complex<double>(0, tmp_aoa[i] + spread_angle[1])).imag(), exp(complex<double>(0, tmp_aoa[i] + spread_angle[1])).real());
            tmp_zod[i] = atan2(exp(complex<double>(0, tmp_zod[i] + spread_angle[2])).imag(), exp(complex<double>(0, tmp_zod[i] + spread_angle[2])).real());
            tmp_zoa[i] = atan2(exp(complex<double>(0, tmp_zoa[i] + spread_angle[3])).imag(), exp(complex<double>(0, tmp_zoa[i] + spread_angle[3])).real());
        }
    }
    ///////////////////////////验证角度扩展//////////////////////////
    for (int i = 0; i < m_cluster.size(); i++)
    {
        m_cluster[i].aod = tmp_aod[i];
        m_cluster[i].aoa = tmp_aoa[i];
        m_cluster[i].zod = tmp_zod[i];
        m_cluster[i].zoa = tmp_zoa[i];
    }
    spread_angle = cale_angular_spreads(m_cluster);
    ASD_desried = spread_angle[4] * 180 / PI;
    ASA_desried = spread_angle[5] * 180 / PI;
    ZSD_desried = spread_angle[6] * 180 / PI;
    ZSA_desried = spread_angle[7] * 180 / PI;

    for (int i = 0; i < m_cluster.size(); i++)
    {
        m_cluster[i].aod = tmp_aod[i] * 180 / PI;
        m_cluster[i].aoa = tmp_aoa[i] * 180 / PI;
        m_cluster[i].zod = 90 - tmp_zod[i] * 180 / PI;
        m_cluster[i].zoa = 90 - tmp_zoa[i] * 180 / PI;
    }
}
