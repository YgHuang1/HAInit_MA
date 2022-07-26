/*
 * @Description: 过滤边
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-02 19:05:23
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 15:35:07
 */
#ifndef Edgefitting_hpp
#define Edgefitting_hpp

#include <iostream>
#include "Reg_definition.hpp"
#include <sophus/se3.hpp>


namespace Edgefilting
{   
    /// 基于初始值对原始的RelMs运动集进行判断 返回判断后的RelMs
    /// @param OriginRelMs 输入:初始RelMs @param RelMsthreshold 输入:阈值, RelMsthreshold = SE(Mj*Mij*Mi^(-1)).log().norm(),
    /// @param Initmotion 输入:初值
    Relativemotions edgefilt(const Relativemotions &OriginRelMs, double RelMsthreshold, const Matrix4x4Arr & Initmotion)
    {
        Relativemotions filtedRelMs;
        filtedRelMs.reserve(OriginRelMs.size());
        for(Relativemotions::const_iterator iter = OriginRelMs.begin(); iter != OriginRelMs.end(); iter++)
        {
            Trans eMij = Initmotion[iter->j] * iter->Tij * Initmotion[iter->i].inverse();
            if(Sophus::SE3d(eMij).log().norm() <= RelMsthreshold)
            {
                filtedRelMs.push_back(*iter);
            }
        }
        return filtedRelMs;
    }
}

#endif