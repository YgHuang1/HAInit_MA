/*
 * @Description:  读取overlop
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-28 21:36:11
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 15:35:27
 */

#ifndef Est_overlop_hpp
#define Est_overlop_hpp

#include <iostream>
#include "Reg_definition.hpp"
#include "ReadWrite_RelMatrixfile.hpp"
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/transforms.h>
#include <pcl/registration/icp.h>
#include "show_model.hpp"
#include <pcl/io/pcd_io.h>
#include <string>
#include <numeric>


namespace EOP
{
    struct Pairs_overlop
    {
    int i,j;
    double overlopij;
    Pairs_overlop(
        int i_ = 0,
        int j_ = 0,
        double overlopij_ = 0.0
    ) : i(i_), j (j_), overlopij(overlopij_)
    {

    }
    
    };

    
    ///@brief  读取overlop文件 文件格式 0 1 0.785 分别表示 datashapeID modelshapeID Overlop
    bool Read_Overlop(const string & overlopfile, vector<EOP::Pairs_overlop> & Scanslop, const double & threshold)
    {
        ifstream fin(overlopfile, ios::in);
        if (!fin.is_open()){
            cout << "overlop读取失败" << endl;
            return false;
        }
        cout << "正在读取overlop文件" << endl;
        Scanslop.clear(); Scanslop.reserve(200);

        string line_info,input_result;
        for(int i = 0; getline(fin, line_info); i++ ){
            EOP::Pairs_overlop templop;
            stringstream input(line_info);
            for(int j = 0; input >> input_result; ++j)
            {
                string::size_type size;
                if(j == 0)
                {
                    templop.i = stoi(input_result, &size);
                }
                if(j == 1)
                {
                    templop.j = stoi(input_result, &size);
                }
                if(j == 2)
                {
                    templop.overlopij = stod(input_result, &size);
                }
            }
            if(templop.overlopij >= threshold)
            {
                Scanslop.push_back(templop);
            }
        }
        fin.close();
        return true;
    }

    ///@brief 基于从lop文件中读取的Vec_overlop 将大于一定重叠面积的匹配对提出出来
    void getRelMSfromlop(const Relativemotions & AllRelMs, const  vector<EOP::Pairs_overlop> & Vec_Overlop, 
                     Relativemotions & lopRelMs, bool ishalf = true)
    {
        lopRelMs.clear(); lopRelMs.reserve(AllRelMs.size());
        Relativemotions_map tempMap = getmap(AllRelMs);
        Pair_set tempset;
        for(auto it : Vec_Overlop)
        {   
            if(it.i > it.j && ishalf) // 如(3,1) 
            {
                // if(tempset.find({it.j, it.i}) != tempset.end()) // 如果找到了(1,3)
                // {
                //     continue;
                // } else
                // {
                //     Relativemotion tempmotion;
                //     tempmotion.i = it.j, tempmotion.j = it.i;
                //     tempmotion.Tij = tempMap.at({it.i,it.j}).Tij.inverse();
                //     lopRelMs.push_back(tempmotion);
                //     tempset.insert({it.j, it.i});
                // }
                continue;
            } else
            {
                lopRelMs.push_back(tempMap.at({it.i, it.j}));
                tempset.insert({it.i,it.j});
            }
            
        }
    }

    bool writr_overlopfile(const string & filename,const vector<EOP::Pairs_overlop> & Scanslop)
    {   
        std::ofstream fout(filename, std::ios::out);
        for(auto it : Scanslop)
        {
            fout << it.i << " " << it.j << " " << it.overlopij << endl;
        }
        fout.flush();
    }
};

#endif