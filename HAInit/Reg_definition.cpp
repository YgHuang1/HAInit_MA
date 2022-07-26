/*
 * @Description:  对应于Reg_definition.hpp
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-02-23 16:07:44
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-27 20:15:42
 */


#include "Reg_definition.hpp"
#include <iostream>
#include <sophus/se3.hpp>

using namespace std;

bool Tripletcompute(
    std::vector<mygraph::Triplet> &vec_triplets,
    Relativemotions & relativemotions,
    std::vector<Motion_Triplets> & motion_triplets
)
{
    std::vector<Motion_Triplets> CyclesAllTri;
    CyclesAllTri.reserve(vec_triplets.size());

    Relativemotions_map map_relatives = getmap(relativemotions);
    Relativemotions_map map_relatives_validated;

    for ( size_t r = 0; r < vec_triplets.size(); ++r){
        const mygraph::Triplet & triplet = vec_triplets[r];
        const IndexT I = triplet.i, J = triplet.j, K = triplet.k;
        // Find the three relative motions
   
        const Pair ij(I,J), ji(J,I);
        const Trans Mij = (map_relatives.count(ij)) ?
        map_relatives.at(ij).Tij : Trans (map_relatives.at(ji).Tij.inverse());

        const Pair jk(J,K), kj(K,J);
        const Trans Mjk = (map_relatives.count(jk)) ?
        map_relatives.at(jk).Tij : Trans (map_relatives.at(kj).Tij.inverse());

        const Pair ki(K,I), ik(I,K);
        const Trans Mki = (map_relatives.count(ki)) ?
        map_relatives.at(ki).Tij : Trans (map_relatives.at(ik).Tij.inverse());

        Trans Mot2Indentity = Mij*Mki*Mjk;
  
        Sophus::SE3d SE3Trans(Mot2Indentity); 
        // cout << i << endl;
        // cout << SE3Trans.log().transpose() << endl;
        // 其中Sophus：：se(3)中前三项为移动，后三项为旋转so(3)
        Vec6 se3_Mot2Indentity = SE3Trans.log().matrix();
        
        Motion_Triplets tempTriplet(I,J,K,se3_Mot2Indentity.norm());
        CyclesAllTri.push_back(tempTriplet); 
    }
 
    map_relatives = std::move(map_relatives_validated);
    motion_triplets.clear();
    motion_triplets = std::move(CyclesAllTri);

    return(!CyclesAllTri.empty()); // 为空表示不存在Triplets

}



bool Tripletcompute_directed(
    const  std::vector<mygraph::Triplet> &vec_triplets,
    const Relativemotions & relativemotions,
    std::vector<NewMotion_Triplets> & Newmotion_triplets
)
{
    std::vector<NewMotion_Triplets> CyclesAllTri; // 用于存储所有的Triplets 以及对应的边 norm
    CyclesAllTri.reserve(vec_triplets.size());

    Relativemotions_map map_relatives = getmap(relativemotions); // 得到relativemotions的map 
    Relativemotions_map map_relatives_validated; // 

    Trans Bigerror;
    Bigerror << 1.0,0.0,0.0,10000,
                0.0,1.0,0.0,10000,
                0.0,0.0,1.0,10000,
                0.0,0.0,0.0,1.0;
    
    for ( size_t r = 0; r < vec_triplets.size(); ++r){
        const mygraph::Triplet & triplet = vec_triplets[r];
        const IndexT I = triplet.i, J = triplet.j, K = triplet.k;
        // 共可以形成6对
        const Pair ij(I,J), ji(J,I); // ij
        const Pair jk(J,K), kj(K,J); // jk
        const Pair ki(K,I), ik(I,K); // ki  Mij*Mki*Mjk = I
        
        // 因此可以形成8对triplets
        Trans Mij = Bigerror, Mji = Bigerror; Trans Mjk = Bigerror, Mkj = Bigerror; Trans Mki = Bigerror, Mik = Bigerror; // 存在6个RelMs
        Mij = (map_relatives.count(ij)) ? map_relatives.at(ij).Tij : Trans (Bigerror);
        Mji = (map_relatives.count(ji)) ? map_relatives.at(ji).Tij : Trans (Bigerror);
        Mjk = (map_relatives.count(jk)) ? map_relatives.at(jk).Tij : Trans (Bigerror);
        Mkj = (map_relatives.count(kj)) ? map_relatives.at(kj).Tij : Trans (Bigerror);
        Mki = (map_relatives.count(ki)) ? map_relatives.at(ki).Tij : Trans (Bigerror);
        Mik = (map_relatives.count(ik)) ? map_relatives.at(ik).Tij : Trans (Bigerror);
        
        std::vector<double> cycleerror; cycleerror.reserve(8); // 存放8个误差
        
        cycleerror.push_back(Sophus::SE3d(Mij * Mki * Mjk).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mij * Mki * Mkj.inverse()).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mij * Mik.inverse() * Mjk).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mij * Mik.inverse() * Mkj.inverse()).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mji.inverse() * Mki * Mjk).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mji.inverse() * Mki * Mkj.inverse()).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mji.inverse() * Mik.inverse() * Mjk).log().norm());
        cycleerror.push_back(Sophus::SE3d(Mji.inverse() * Mik.inverse() * Mkj.inverse()).log().norm());

        auto smallest = std::min_element(cycleerror.begin(), cycleerror.end());
        int Count = std::distance(cycleerror.begin(), smallest);
        
        switch (Count)
        {
            case 0:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ij, jk, ki, Sophus::SE3d(Mij * Mki * Mjk).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 1:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ij, kj, ki, Sophus::SE3d(Mij * Mki * Mkj.inverse()).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 2:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ij, jk, ik, Sophus::SE3d(Mij * Mik.inverse() * Mjk).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 3:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ij, kj, ik, Sophus::SE3d(Mij * Mik.inverse() * Mkj.inverse()).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 4:
                {
                     NewMotion_Triplets temptriplets(I, J, K, ji, jk, ki, Sophus::SE3d(Mji.inverse() * Mki * Mjk).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 5:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ji, kj, ki, Sophus::SE3d(Mji.inverse() * Mki * Mkj.inverse()).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 6:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ji, jk, ik, Sophus::SE3d(Mji.inverse() * Mik.inverse() * Mjk).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            case 7:
                {
                    NewMotion_Triplets temptriplets(I, J, K, ji, kj, ik, Sophus::SE3d(Mji.inverse() * Mik.inverse() * Mkj.inverse()).log().norm());
                    CyclesAllTri.push_back(temptriplets);
                    break;
                }
            default:
                break;
        }

        // if(I == 4 && J == 5 && K == 6)
        // {
        //     cout << Count << endl;
        //     cout << Sophus::SE3d(Mij * Mik.inverse() * Mjk).log().norm() << endl;
        //     cout << Mij << endl;
        //     cout << Mik << endl;
        //     cout << Mjk << endl;
        // }
        // Mji = (map_relatives.count)
    }
    Newmotion_triplets.clear();
    Newmotion_triplets = std::move(CyclesAllTri);
}