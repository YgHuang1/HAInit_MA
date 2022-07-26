/*
 * @Description: L1-IRLS 部分代码参考OpenMVG
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-04-07 19:54:53
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 15:36:10
 */

#ifndef L1_IRLS_hpp
#define L1_IRLS_hpp

#include "Reg_definition.hpp"
#include <lemon/lemon/adaptors.h>
#include <lemon/lemon/dfs.h>
#include <lemon/lemon/kruskal.h>
#include <lemon/lemon/list_graph.h>
#include <lemon/lemon/path.h>
#include <queue>

#include <sophus/se3.hpp>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <openMVG/numeric/l1_solver_admm.hpp>

namespace sloveL1_IRLS
{
    struct Node
    {
        using InternalType = IndexArr; // IndexArr -> std::vector<uint32_t>;
        InternalType edges; //  array of vertex indices
    };

    using NodeArr = std::vector<Node>;

     struct Link {
      uint32_t ID; // node index
      uint32_t parentID; // parent link
      inline Link(uint32_t ID_=0, uint32_t parentID_=0) : ID(ID_), parentID(parentID_) {}
    };
    
    using LinkQue = std::queue<Link>; // queue 先进先出的数据结构

    using graph_t = lemon::ListGraph;
    using map_EdgeMap = graph_t::EdgeMap<double>;

    using MapEdgeIJ2M = std::map<std::pair<uint32_t,uint32_t>, Matrix4x4>;

    // 根据相对运动的图寻找最小生成树
    uint32_t FindMaximumSpanningTree(const Relativemotions& RelRs, graph_t& g, MapEdgeIJ2M& mapIJ2R, NodeArr& minGraph)
    {
        if(RelRs.empty()){
            return -1;
        };

        // A- Compute the number of node we need
        std::set<uint32_t> setNodes;
        for(size_t p = 0;p < RelRs.size(); ++p){
            const Relativemotion & relR = RelRs[p];
            setNodes.insert(relR.i);
            setNodes.insert(relR.j);
        }
        
        
        // B- Create a node graph for each element of the set
        using map_NodeT = std::map<uint32_t,graph_t::Node>;
        map_NodeT map_index_to_node;
        for(const auto & iter : setNodes){
            map_index_to_node[iter] = g.addNode();
        }

        // C- create a graph from RelRs with weighted edges
        map_EdgeMap map_edgeMap(g);
        for (size_t p = 0; p < RelRs.size(); ++p) {
            const Relativemotion& relR = RelRs[p];
            pair<uint32_t,uint32_t> IJ(relR.i,relR.j);
            pair<uint32_t,uint32_t> JI(relR.j,relR.i);
            mapIJ2R.emplace(IJ,relR.Tij);
            Matrix4x4 Tji = relR.Tij;
            mapIJ2R.emplace(JI,Tji);
        
            // add edge to the graph
            graph_t::Edge edge =  g.addEdge(map_index_to_node[relR.i], map_index_to_node[relR.j]);
            map_edgeMap[ edge ] = -1.0f;
        }

        // D- Computer the MST of the graph
        std::vector<graph_t::Edge> tree_edge_vec;
        lemon::kruskal(g, map_edgeMap, std::back_inserter(tree_edge_vec));

        const size_t nViews = lemon::countNodes(g);
        minGraph.resize(nViews);
        for (size_t i= 0; i < tree_edge_vec.size(); i++)
        {
            minGraph[g.id(g.u(tree_edge_vec[i]))].edges.push_back(g.id(g.v(tree_edge_vec[i])));
            minGraph[g.id(g.v(tree_edge_vec[i]))].edges.push_back(g.id(g.u(tree_edge_vec[i])));
        }
        return tree_edge_vec.size();

    }

 
    /// @param[in] RelRs -> Relativemotions的集合 
    /// @return Rs -> 输出的初值
    /// @brief 根据Relativemotions初始化MA的初值, 该函数不需要进行DFS确定最大连通块，但不是BFS方法
    bool InitMotionsMAT(
      const Relativemotions & RelRs,
      Matrix4x4Arr & Rs,
      const uint32_t nMainViewID
    )
    {
        if(Rs.empty())
            return false;
        graph_t g;
        MapEdgeIJ2M mapIJ2M;
        NodeArr minGraph;
        FindMaximumSpanningTree(RelRs, g, mapIJ2M, minGraph);
        g.clear();

        for(int i = 0;i < minGraph.size(); i++){
            auto temp = minGraph[i];
            cout << "************" << endl;
            for(int j = 0; j < temp.edges.size();j++){
            cout << temp.edges[j] << endl;
            }
        }
        cout << minGraph.size() << endl;

        LinkQue stack;
        stack.push(Link(nMainViewID,uint32_t(0)));
        // Rs[nMainViewID] = Matrix4x4::Identity();
        do{
            const Link &link = stack.front();
            const Node &node = minGraph[link.ID];
            
            for(Node::InternalType::const_iterator pEdge = node.edges.begin();
            pEdge != node.edges.end(); ++ pEdge ){
                const size_t edge = *pEdge;
                if (edge == link.parentID){
                    assert(mapIJ2M.find({link.parentID,link.ID}) != mapIJ2M.end());
                    std::pair<uint32_t,uint32_t> temppair(link.parentID,link.ID);
                    Matrix4x4 Mij = mapIJ2M[temppair];
                    if(link.ID > link.parentID){
                        Trans temp = Mij.inverse() * Rs[link.parentID];
                        Rs[link.ID] << temp(0,0), temp(0,1), temp(0,2), temp(0,3),
                                    temp(1,0), temp(1,1), temp(1,2), temp(1,3),
                                    temp(2,0), temp(2,1), temp(2,2), temp(2,3),
                                    temp(3,0), temp(3,1), temp(3,2), temp(3,3);
                    }else{
                        Trans temp = Mij*Rs[link.parentID];
                        Rs[link.ID] << temp(0,0), temp(0,1), temp(0,2), temp(0,3),
                                    temp(1,0), temp(1,1), temp(1,2), temp(1,3),
                                    temp(2,0), temp(2,1), temp(2,2), temp(2,3),
                                    temp(3,0), temp(3,1), temp(3,2), temp(3,3);
                    }

                cout << link.ID << "  ,  " << link.parentID << endl; // 显示节点

            } else{
                stack.push(Link(edge,link.ID));
            }
            //  cout << link.parentID << "  ,  " << link.ID << endl;
            }
            stack.pop();
        } while(!stack.empty());
    };


    bool BFSMotionMAT(
        const Relativemotions & RelMs,
        const vector<Pair> &SrcPair_set, 
        Matrix4x4Arr & Rs,
        bool is_ZeroMain = true)
    {   
        Relativemotions_map map_RelMs = getmap(RelMs);
        for(vector<Pair>::const_iterator iter = SrcPair_set.begin(); iter != SrcPair_set.end(); iter++){
            if(iter->first < iter->second){ 
                assert(map_RelMs.find({iter->first, iter->second}) != map_RelMs.end());
                Pair tempPair(iter->first, iter->second);
                Trans Mij = map_RelMs[tempPair].Tij;
                Trans temp = Mij.inverse()*Rs[iter->first];
                Rs[iter->second] = temp;
            } else{
                assert(map_RelMs.find({iter->second, iter->first}) != map_RelMs.end());
                Pair tempPair(iter->second, iter->first);
                Trans Mij = map_RelMs[tempPair].Tij;
                Trans temp = Mij*Rs[iter->first];
                Rs[iter->second] = temp;
            }
        }

        if(is_ZeroMain)
        {
            Trans ZeroM_inv = Rs[0].inverse();
            for(size_t r = 0; r < Rs.size(); r++)
            {
                Trans &Ri = Rs[r];
                Ri = Ri * ZeroM_inv;
            }
        }

        return true;
    }


    /**
     * @Description: apply correction to global Motions
     * @param x：x为误差迭代累乘项
     * @param nMainViewID：主视图
     * @param Rs:全局运动
     * @return {*}
     * @author: hyg
     */    
    void CorrectMatrix(
        const Eigen::MatrixXd & x,
        const uint32_t nMainViewID,
        Matrix4x4Arr& Rs
    )
    {
        for(size_t r = 0; r < Rs.size(); ++r){
            if (r == nMainViewID)
            continue;
            Matrix4x4 &Ri = Rs[r];
            const uint32_t i = (r<nMainViewID ? r : r-1);
            const Eigen::Matrix<double,6,1> eRid = Eigen::Matrix<double,6,1>(x.block<6,1>(6*i,0));
            const Trans eRi;
            Trans updated_Ri = Sophus::SE3d::exp(eRid).matrix();
            Ri = updated_Ri * Ri;
            // const Eigen::Matrix<double,
        }
    }


    void FillErrorMatrix(
        const Relativemotions& RelRs,
        const Matrix4x4Arr& Rs,
        Eigen::VectorXd &b
    ){
        for(size_t r = 0; r < RelRs.size(); ++r){
            const Relativemotion& relR = RelRs[r];
            const Trans& Ri = Rs[relR.i];
            const Trans& Rj = Rs[relR.j];
            const Trans& Rij = relR.Tij;
            const Trans eRij = Rj * Rij * Ri.inverse();// 计算误差矩阵
            Eigen::Matrix<double,6,1> erij;
            
            erij = Sophus::SE3d(eRij).log();
            b.block<6,1>(6*r,0) = erij; // 将相对误差矩阵转换为李代数的误差形式保存在向量b中
        }
    }

    void FillMappingMatrix(const Relativemotions& RelRs,
        const uint32_t nMainViewID,
        Eigen::SparseMatrix<double> &A)
    {
        A.reserve(A.rows()*2);
        Eigen::SparseMatrix<double>::Index i = 0,j = 0;
        for (size_t r = 0; r < RelRs.size(); ++r){
            Relativemotion relR = RelRs[r];
            if(relR.i != nMainViewID){
            j = 6*(relR.i<nMainViewID ? relR.i : relR.i -1);
            A.insert(i+0,j+0) = -1.0;
            A.insert(i+1,j+1) = -1.0;
            A.insert(i+2,j+2) = -1.0;
            A.insert(i+3,j+3) = -1.0;
            A.insert(i+4,j+4) = -1.0;
            A.insert(i+5,j+5) = -1.0;
            } 
            if (relR.j != nMainViewID){
            j = 6*(relR.j<nMainViewID ? relR.j : relR.j-1);
            A.insert(i+0,j+0) = 1.0;
            A.insert(i+1,j+1) = 1.0;
            A.insert(i+2,j+2) = 1.0;
            A.insert(i+3,j+3) = 1.0;
            A.insert(i+4,j+4) = 1.0;
            A.insert(i+5,j+5) = 1.0;
            }
            i+=6;
        }
        A.makeCompressed();
    }

    bool SolveL1RA(
      const Relativemotions& RelRs,
      Matrix4x4Arr& Rs,
      const Eigen::SparseMatrix<double>& A,
      const unsigned int nMainViewID
    ){
        const unsigned nObss = (unsigned)RelRs.size();
        const unsigned nVars = (unsigned)Rs.size() -1;
        const unsigned m = nObss*6;
        const unsigned n = nVars*6;

        const unsigned b_m = nObss*6;

        Eigen::VectorXd x(Eigen::VectorXd::Zero(n)), b(b_m);
        
        // Current error and the previous one
        double e = std::numeric_limits<double>::max(), ep;
        unsigned iter = 0;
        // L1RA iterate optimation till the desired precision is reached 
        using namespace openMVG;
        do{
            FillErrorMatrix(RelRs,Rs,b);
            cout << "第" << iter << "轮迭代误差" << endl;
            // cout << b << endl;
            L1Solver<Eigen::SparseMatrix<double>>::Options options;
            L1Solver<Eigen::SparseMatrix<double>> l1_solver(options,A);
            l1_solver.Solve(b,&x);

            ep = e; e= x.norm();
            if (ep<e)
            break;
            CorrectMatrix(x, nMainViewID, Rs);
            // cout << iter << endl;
        
        } while (++iter <32 && e > 1e-5 &&(ep-e)/e>1e-2);

        std::cout << "L1RA Converged in " << iter << " iterations." << std::endl;
        return true;
    }



}



#endif
