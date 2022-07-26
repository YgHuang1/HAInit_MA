/*
 * @Description: 数值分析
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-02-24 14:14:05
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 13:49:36
 */

#ifndef Analysis_numeric_hpp
#define Analysis_numeric_hpp

#include <iostream>
#include "Reg_definition.hpp"
#include <pcl/console/print.h>
#include <sophus/se3.hpp>
#include <numeric>
#include <iterator>
#include <Eigen/Eigen>


using namespace std;

namespace Analysis_numeric
{   
    /**
     * @brief:  
     * @param[in] begin 统计数据vector的初始迭代器
     * @param[in] end 统计数据vector的末尾迭代器
     * @param[out] min max mean median
     * @author: hyg
     */    
    template <typename Type, typename DataInputIterator>
    bool minMaxMeanMedian( DataInputIterator begin, DataInputIterator end,
                       Type & min, Type & max, Type & mean, Type & median )
    {
        if(std::distance(begin,end) < 1){
            return false;
        }

        std::vector<Type> vec_val(begin,end); // 将[begin，end]中的元素拷贝给vec_value
        // Get the median value:
        const auto middle = vec_val.begin() + vec_val.size() / 2;
        std::nth_element(vec_val.begin(), middle, vec_val.end());
        median = *middle;
        min = *std::min_element(vec_val.begin(), middle);
        max = *std::max_element(middle, vec_val.end());
        mean = std::accumulate( vec_val.cbegin(), vec_val.cend(), Type( 0 ) )
            / static_cast<Type>( vec_val.size() );
        return true;
    }



    /**
     * @brief: 统计匹配结果与真实结果的误差比较结果
     * @param[in] Reg_Ans 匹配结果
     * @param[in] GroundTruth 真实结果
     * @param[out] Analysis_Rotans 旋转误差统计结果
     * @param[out] Analysis_TransAns 移动误差统计结果
     * @author: hyg
     */    
    bool AnsAnalysis_betweenTruth(
        const Matrix4x4Arr & Reg_Ans, const Matrix4x4Arr & GroundTruth, 
        std::vector<double> & Analysis_Rotans, std::vector<double> & Analysis_TransAns
    )
    {
        if(Reg_Ans.size() != GroundTruth.size())
        {
            cerr << "匹配得到的结果与真实值阵数不一致" << endl;
            return false;
        }

        std::vector<double> error_rots, error_trans;
        error_rots.reserve(Reg_Ans.size()); error_trans.reserve(Reg_Ans.size());

        for(int i = 0; i < Reg_Ans.size(); i++){
            Trans Truthi = GroundTruth[i].inverse();
            Trans Reg_Ansi = Reg_Ans[i];
            Trans errorMatrix = Reg_Ansi * Truthi;
            Sophus::SO3d e(errorMatrix.block(0,0,3,3));
            double temp_error_rotation = e.log().norm();
            // double temp_error_rotation = (GroundTruth[i] - Reg_Ans[i]).block(0,0,3,3).norm();

            error_rots.push_back(temp_error_rotation); // 旋转误差
            cout << temp_error_rotation << endl;
            double temp1 = Reg_Ans[i](0,3) - GroundTruth[i](0,3);
            double temp2 = Reg_Ans[i](1,3) - GroundTruth[i](1,3);
            double temp3 = Reg_Ans[i](2,3) - GroundTruth[i](2,3);
            Vec3 temp_error_trans(temp1,temp2,temp3);
            error_trans.push_back(temp_error_trans.norm()); // 移动误差
        }

        float Rmin,Rmax,Rmean,Rmedian; 
        minMaxMeanMedian(error_rots.begin(),error_rots.end(),Rmin,Rmax,Rmean,Rmedian);
        Analysis_Rotans.clear(); Analysis_Rotans.reserve(4);
        Analysis_Rotans.push_back(Rmin); Analysis_Rotans.push_back(Rmax); 
        Analysis_Rotans.push_back(Rmean); Analysis_Rotans.push_back(Rmedian); 


        float tmin,tmax,tmean,tmedian; 
        minMaxMeanMedian(error_trans.begin(), error_trans.end(),tmin,tmax,tmean,tmedian);
        Analysis_TransAns.clear(); Analysis_TransAns.reserve(4);
        Analysis_TransAns.push_back(tmin); Analysis_TransAns.push_back(tmax);
        Analysis_TransAns.push_back(tmean);Analysis_TransAns.push_back(tmedian);

        return true;
    }
    

    /**
     * @brief: 统计相对运动集的误差结果 注意 Mij = Mi * Mj^(-1)
     * @param[in] GroundTruth 真实值
     * @param[in] RelMs 相对运动集
     * @param[out] Analysis_Rotans 旋转误差统计结果
     * @param[out] Analysis_TransAns 移动误差统计结果
     * @author: hyg
     */    
    bool RelAnalysis_betweenTruth(
        const Matrix4x4Arr & GroundTruth, const Relativemotions & RelMs,
        std::vector<double> & Analysis_Rotans, std::vector<double> & Analysis_TransAns
    )
    {   
        const int PairSize = RelMs.size();
        std::vector<double> error_rots, error_trans;
        error_rots.reserve(PairSize), error_trans.reserve(PairSize);

        for(int p = 0; p < PairSize; p++){
            const Relativemotion & relM = RelMs[p];
            const Trans & Mi = GroundTruth[relM.i];
            const Trans & Mj = GroundTruth[relM.j];
            const Trans & Mij = relM.Tij;
            const Trans eMij = Mj * Mij* Mi.inverse();

            Sophus::SO3d e(eMij.block(0,0,3,3));
            error_rots.push_back(e.log().norm());

            double temp1 = eMij(0,3), temp2 = eMij(1,3), temp3 = eMij(2,3);
            Vec3 temp_error_trans(temp1,temp2,temp3);
            error_trans.push_back(temp_error_trans.norm());
        }

        float Rmin,Rmax,Rmean,Rmedian; 
        minMaxMeanMedian(error_rots.begin(),error_rots.end(),Rmin,Rmax,Rmean,Rmedian);
        Analysis_Rotans.clear(); Analysis_Rotans.reserve(4);
        Analysis_Rotans.push_back(Rmin); Analysis_Rotans.push_back(Rmax); 
        Analysis_Rotans.push_back(Rmean); Analysis_Rotans.push_back(Rmedian); 


        float tmin,tmax,tmean,tmedian; 
        minMaxMeanMedian(error_trans.begin(), error_trans.end(),tmin,tmax,tmean,tmedian);
        Analysis_TransAns.clear(); Analysis_TransAns.reserve(4);
        Analysis_TransAns.push_back(tmin); Analysis_TransAns.push_back(tmax);
        Analysis_TransAns.push_back(tmean);Analysis_TransAns.push_back(tmedian);

        return true;
    }

    /**
     * @brief: 
     * @param[in] error_rots 旋转误差vector
     * @param[in] error_trans 移动误差vector
     * @author: hyg
     */    
    void cout_RotandTrans(const std::vector<double> &error_rots,const std::vector<double> &error_trans)
    {
        cout << "旋转误差: \n" << "min: " << error_rots[0] << "  max:" << error_rots[1] << "  mean:" << error_rots[2];
        cout << " median:" << error_rots[3] << endl;

        cout << "移动误差: \n" << "min: " << error_trans[0] << "  max:" << error_trans[1] << "  mean:" << error_trans[2];
        cout << " median:" << error_trans[3] << endl;
    }
    

    /**
     * @brief: 统计相对运动集的误差结果 注意 Mij = Mj^(-1) * Mi
     * @param[in] GroundTruth 真实值
     * @param[in] RelMs 相对运动集
     * @param[out] Analysis_ans 误差矩阵的李代数的norm
     * @author: hyg
     */ 
    bool RelAnalysis_betweenTruth(const Matrix4x4Arr & GroundTruth, const Relativemotions & RelMs,
        std::vector<double> & Analysis_ans)
    {
        const int PairSize = RelMs.size();
        Analysis_ans.clear();
        Analysis_ans.reserve(PairSize);
        for(int p = 0; p < PairSize; p++){
            const Relativemotion & relM = RelMs[p];
            const Trans & Mi = GroundTruth[relM.i];
            const Trans & Mj = GroundTruth[relM.j];
            const Trans & Mij = relM.Tij;
            Trans eMij  = Mj * Mij * Mi.inverse();
            Sophus::SE3d e(eMij);
            Analysis_ans.push_back(e.log().norm());
        }
        return true;
    }

} // namespace Analysis_numeric

#endif