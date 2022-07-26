/*
 * @Description: 用于进行顶点的初始化
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-02 16:50:03
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-11 16:07:58
 */
#ifndef Initmotion_hpp
#define Initmotion_hpp

#include "Reg_definition.hpp"
#include "Graph_Theory.hpp"

namespace InitRelMs
{
    /**
     * @brief: 
     * @param[in] RelMs 相对运动集的RelMs
     * @param[in] MSTpair 最小生成树
     * @param[in] is_ZeroMain 是否将第0视角帧设为参考帧
     * @param[out] initmotion 初始化的运动集
     * @author: hyg
     */    
    bool InitMotion_fromMST(
        const Relativemotions &RelMs,
        const std::set<Pair> &MSTpair,
        Matrix4x4Arr &initmotion,
        const int nmainviewer = 0,
        bool is_ZeroMain = true
    )
    {    
        cout << "根据MST进行BFS搜索并进行初始值的确认" << endl;
        Relativemotions_map map_RelMs = getmap(RelMs);
        vector<Pair> pairs;
        Graph_Theory::BFS(MSTpair, pairs, nmainviewer);
        if(pairs.empty())
        {
            std::cerr << "生成树初始化失败" << endl;
        }
        for(auto iter : pairs)
        {
            cout << iter.first << " , " << iter.second << endl;
        }
        cout << " -------------" << endl;

    
        for(auto iter : pairs)
        {   
            Pair temppair(iter.first, iter.second);
            Trans Mij = Trans::Identity();
            if(MSTpair.find(temppair) != MSTpair.end())
            {   
                initmotion[iter.second] = initmotion[iter.first] *  map_RelMs.find(temppair)->second.Tij.inverse();               
            } else
            {   
                temppair.first = iter.second, temppair.second = iter.first;
                initmotion[iter.second] = initmotion[iter.first] *  map_RelMs.find(temppair)->second.Tij;
            }
        }

        Trans nMainvierinv = initmotion[0].inverse();
        for(int i = 0; i < initmotion.size(); i++)
        {
            initmotion[i] = initmotion[i] * nMainvierinv;
        }
        
        
    }
}


#endif