/*
 * @Description: 该程序用于读和写相对运动矩阵文件
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-02-23 19:22:45
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-03 10:45:02
 */

#ifndef ReadWrite_RelMatrixfile_hpp
#define ReadWrite_RelMatrixfile_hpp

#include <fstream>
#include <ostream>
#include <string>
#include <iostream>
#include "Reg_definition.hpp"
#include <pcl/console/print.h>



using namespace std;

namespace Matrix_Read
{
    /**
     * @brief: 基于文件，读入相对运动
     * @param[in] Inputfile 文件的路径
     * @param[out] RelMs 相对运动
     * 注意： 
     * 文件保存形式为
             0  ,  3 // 该行为索引 索引中逗号隔开
            0.401595 -0.000928777    -0.915817   0.00139429
        -0.00890152     0.999948   -0.0049175  0.000206707
            0.915774     0.010127     0.401566 -0.000156769
                0            0            0            1
        如果需要读取MSE, 则第一行为0 , 1 0.1555
     * @author: hyg
     */    
    bool Read_RelMatrix(const string & InputRelMotions, Relativemotions & RelMs, bool isReadMSE = true);



    /**
     * @brief: 基于文件，读入全局运动
     * @param[in] Inputfile 文件的路径
     * @param[in] GlobalNumber 全局点云帧数
     * @param[out] GlobalMats 全局运动 std<vector>
     * @author: hyg
     */    
    bool Read_GlobalMatrix(const string & InputGlobalfile, Matrix4x4Arr & GlobalMats, int  GlobalNumber, bool format = true);


}

namespace Matrix_Write
{   
    /**
     * @brief: 输出相对运动矩阵
     * @param[out] OutputTxt
     * @param[in] RelMs
     * @author: hyg
     */    
    bool Write_Relmatrix(const string & OutputTxt, Relativemotions & RelMs, bool iswirteMSE = true);


    /// @brief 保存全局矩阵
    /// @param OutputTxt 输入:全局矩阵保存文件地址 @param vec_Global 输入:存储全局矩阵的vector  
    bool Write_GlobalMatrix(const string & OutputTxt, Matrix4x4Arr & vec_Global);
}

#endif