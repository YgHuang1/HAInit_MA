/*
 * @Description: 将Relativemotions中全局下的索引转换为局部下的索引
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-11 21:40:08
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-27 15:04:49
 */

#ifndef Global2clus_hpp
#define Global2clus_hpp

#include "Reg_definition.hpp"


namespace Global2clus
{
    void convertRelMs(const set<int> & nodes, Relativemotions & RelMsneedtoconnvert)
    {   
        if(RelMsneedtoconnvert.empty())
        {
            std::cerr << "需要进行转换的RelMs为空" << endl;
        }
        Relativemotions convertdRelms;
        convertdRelms.reserve(RelMsneedtoconnvert.size());

        for(Relativemotions::const_iterator iter = RelMsneedtoconnvert.begin(); 
            iter != RelMsneedtoconnvert.end(); iter++ )
        {
            IndexT indexi = iter->i, indexj = iter->j;
            Relativemotion singRelMs = *iter;
            IndexT ii = std::distance(nodes.begin(), nodes.find(indexi));
            IndexT jj = std::distance(nodes.begin(), nodes.find(indexj));
            singRelMs.i = ii, singRelMs.j = jj;
            convertdRelms.push_back(singRelMs);
        }

        RelMsneedtoconnvert = std::move(convertdRelms);
        convertdRelms.clear();
    }

    void convertRelMs(const std::vector<int> & nodes, Relativemotions & RelMsneedtoconnvert)
    {   
        if(RelMsneedtoconnvert.empty())
        {
            std::cerr << "需要进行转换的RelMs为空" << endl;
        }
        Relativemotions convertdRelms;
        convertdRelms.reserve(RelMsneedtoconnvert.size());

        for(Relativemotions::const_iterator iter = RelMsneedtoconnvert.begin(); 
            iter != RelMsneedtoconnvert.end(); iter++ )
        {
            IndexT indexi = iter->i, indexj = iter->j;
            Relativemotion singRelMs = *iter;
            IndexT ii = std::distance(nodes.begin(), find(nodes.begin(), nodes.end(), indexi));
            IndexT jj = std::distance(nodes.begin(), find(nodes.begin(), nodes.end(), indexj));
            singRelMs.i = ii, singRelMs.j = jj;
            convertdRelms.push_back(singRelMs);
        }

        RelMsneedtoconnvert = std::move(convertdRelms);
        convertdRelms.clear();
    }

}


#endif