/*
 * @Description: 测试HA Mij = Mj^-1 * Mi
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-26 13:59:32
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 16:21:52
 */
#include "HAPipeline.hpp"
#include "Est_overlop.hpp"
#include "third_party/cmdLine/cmdLine.h"

bool PCRwithOursInit(const string & overlopfile, const string & RelMsfile, const string & Truthfile, 
                    const double & param_overlop, const double & param_threshold,const int & ScansNumber,
                    const string &BunnyPLYfile)
{
    vector<EOP::Pairs_overlop> Vec_Overlop;
    if(!EOP::Read_Overlop(overlopfile, Vec_Overlop, param_overlop)) // 读取overlop文件
    {   
        std::cerr << "初始化失败" << endl;
        return false;
    }
    Relativemotions AllRelMs; 
    Matrix_Read::Read_RelMatrix(RelMsfile, AllRelMs, false); // 读取总的RelMs文件

    Relativemotions lopRelMsij; 
    EOP::getRelMSfromlop(AllRelMs, Vec_Overlop, lopRelMsij); // 得到对应于Vec_Overlop的RelMSij

    PointRegistrationPipeline::HAInit_ByMSE HAInit;
    int number = ScansNumber;
    const string truthfile = Truthfile;
    clock_t Initstart = clock();
    HAInit.Init(lopRelMsij, truthfile, number);
    HAInit.Initnodes(0.0005); // 初始化第一轮顶点
    HAInit.expandNodes(); // 扩展顶点
    clock_t Initend = clock();

    cout << "Our初始化时间为" <<(double) (Initend - Initstart)/CLOCKS_PER_SEC*1000 << "ms" << endl;
    cout << "初值误差为:" << endl;
    HAInit.comparewithTruth(); // 比较结果
    
    vector<ShowPointCloud::PCD_model> Bunny;
    ShowPointCloud::Get_PCDmodel(HAInit.GlobalInitMotion, BunnyPLYfile,Bunny);
    ShowPointCloud::visualAllCloud(Bunny);
}


int main(int argc, char **argv)
{    
    
    std::string overlopfile = argv[1]; // 存放每个匹配对对应的重叠面积
    std::string RelMsfile = argv[2]; // 每个相对运动的结果
    std::string Truthfile = argv[3]; // 真值
    std::string ModelFile = argv[4]; // 存放model的路径 用于最后的结果展示
    double param_overlop = atof(argv[5]); // 选取重叠面积大于该参数的匹配对所为算法的输入 一般0.3
    double param_threshold = atof(argv[6]); // 设定边的过滤的参数  0.12
    int ScansNumber = atoi(argv[7]); // 全局运动的个数

    PCRwithOursInit(overlopfile,RelMsfile, Truthfile, param_overlop, param_threshold, ScansNumber,ModelFile);
    // PCRwithOursInit("../imagedata/BunnyData/Bunnyoverlop.txt","../imagedata/BunnyData/Bunny_ICPAns.txt",
    //                 "../imagedata/BunnyData/Bunny_Groundtruth.txt", 0.2, 0.12, 10,
    //                 "../imagedata/Bunny_PCD");
    
}

