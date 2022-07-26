/*
 * @Description: 该程序用于定义多视角点云配准所需要的定义
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-02-23 15:26:22
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 14:08:50
 */

#ifndef REG_DEFINITION_HPP
#define REG_DEFINITION_HPP

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <unordered_map>
#include <sophus/se3.hpp>

using namespace std;

using IndexT = uint32_t;       // 索引
using Trans = Eigen::Matrix4d; // 变换矩阵
using Pair = std::pair<IndexT, IndexT>;
using Pair_set = std::set<Pair>;
using Matrix4x4 = Trans;
using Matrix4x4Arr = std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>;
using IndexArr = std::vector<uint32_t>;
using Nodesset = std::set<IndexT>; // 用来存放节点

typedef Eigen::Matrix<double, 6, 1> Vec6; // 李代数
typedef Eigen::Vector3d Vec3;

struct Relativemotion // 定义相对运动 i j motionij
{
    IndexT i, j;
    Trans Tij;
    double TrimmedMSE;

    Relativemotion(
        IndexT i_ = 0,
        IndexT j_ = 0,
        const Trans &Tij_ = Trans::Identity()
        ) : i(i_), j(j_), Tij(Tij_)
    {
    }

    void Print()
    {
        std::cout << "Index: " << i << " , " << j << "   ";
        std::cout << "TrimmedMSE: " << TrimmedMSE << std::endl;
        std::cout << "Transformation: " << Tij << std::endl;
    }
};

class CompareRelativeMotion
{
    public:
        bool operator()(const Relativemotion & M1, const Relativemotion & M2)
        {
            return M1.TrimmedMSE <= M2.TrimmedMSE;
        }
};

struct InitializedMotion // 定义用于初始化的全局运动
{
    IndexT i;
    Trans Ti;

    InitializedMotion(
        IndexT i_ = 0,
        const Trans & Ti_ = Trans::Identity()
    ): i(i_), Ti(Ti_){
        
    }
};

using InitializedMotions = std::vector<InitializedMotion>;

using Relativemotions = std::vector<Relativemotion>;
using Relativemotions_MSEset = std::set<Relativemotion, CompareRelativeMotion>;
using Relativemotions_map = std::map<Pair, Relativemotion>; // 可以根据Pair 获得Relativemotion

struct Motion_Triplets
{
    IndexT i, j, k;
    double CyclesNorm;
    Motion_Triplets(
        IndexT i_ = 0,
        IndexT j_ = 0,
        IndexT k_ = 0,
        double CyclesNorm_ = 0) : i(i_), j(j_), k(k_), CyclesNorm(CyclesNorm_)
    {
    }

    void Print()
    {
        cout << i << " , " << j << " , " << k << endl;
        cout << CyclesNorm << endl;
    }
};

struct NewMotion_Triplets
{
    IndexT i, j, k;
    Pair IJ, JK, KI; 
    double CyclesNorm;
    NewMotion_Triplets(
        IndexT i_,
        IndexT j_,
        IndexT k_,
        Pair ij_, 
        Pair jk_, 
        Pair ki_,
        double CyclesNorm_ ) : i(i_), j(j_), k(k_), IJ(ij_), JK(jk_), KI(ki_), CyclesNorm(CyclesNorm_)
    {
    }

    NewMotion_Triplets()
    {

    };
    
};

/// @brief NewMotion_Triplets排序
class CompareNewMotion_Triplets
{
    public:
        bool operator()(const NewMotion_Triplets & M1, const NewMotion_Triplets & M2)
        {
            return M1.CyclesNorm < M2.CyclesNorm;
        }
};

class CompareMotion_Triplets
{
    public:
        bool operator()(const Motion_Triplets & M1, const Motion_Triplets & M2)
        {
            return M1.CyclesNorm < M2.CyclesNorm;
        }
};


/**
 * @brief  根据Relativemotions得到Pair_set
 * @param relMots Relativemotions
 * @return Pair_set
 * @author: hyg
 */
inline Pair_set getPairs(const Relativemotions & relMots){
    Pair_set pairs;
    for (const auto & cur_motion : relMots)
        pairs.insert({cur_motion.i,cur_motion.j});
    return pairs;
}

/**
 * @brief 根据Relativemotions得到Relativemotions_map
 * @param retMots Relativemotions
 * @return Relativemotions_map
 * @author: hyg
 */
inline Relativemotions_map getmap(const Relativemotions & retMots)
{
    Relativemotions_map map_mots;
    cout << "正在生成相对运动集的map结构体" << endl;
    for (auto iter = retMots.begin(); iter != retMots.end(); iter ++ ){
        Pair temp_pair(iter->i, iter->j);
        map_mots.emplace(temp_pair, *iter);
    }

    return map_mots;
}

// 创建set 其中key值为Pair(i,j) value值为MSE 同时对MSE进行升序排列
inline Relativemotions_MSEset getMSEset(const Relativemotions & retMots)
{   
    Relativemotions_MSEset setMSE_mots;
    cout << "正在生成相对运动集的set结构体" << endl;
    for(auto iter = retMots.begin(); iter != retMots.end(); iter++)
    {
        setMSE_mots.insert(*iter);
        // cout << iter->i << " " << iter->j << endl;
    }
    // cout << setMSE_mots.size() << endl;
    return setMSE_mots;

}

namespace mygraph
{
/**
* @brief Simple container for a tuple of three value
* @note It is used to store the node id of triplets of a graph.
*/
struct Triplet
{

    /**
     * @brief Constructor
     * @param ii First element of the triplet
     * @param jj Second element of the triplet
     * @param kk Third element of the triplet
     */
    Triplet( IndexT ii, IndexT jj, IndexT kk )
        : i( ii ), j( jj ), k( kk )
    { }

    /**
     * @brief Indicate if an edge contains one of the element of the triplet
     * @param edge Edge to test
     * @retval true if edge contains at least one of index of the triplet
     * @retval false if edge contains none of the index of the triplet
     */
    bool contain( const std::pair<IndexT, IndexT> & edge ) const
    {
        const IndexT It = edge.first;
        const IndexT Jt = edge.second;
        return ( ( It == i || It == j || It == k ) &&
                ( Jt == i || Jt == j || Jt == k ) && It != Jt );
    }

    /// the three triplet index id
    IndexT i, j, k;
};


/**
* @brief Return triplets contained in the graph build from IterablePairs
* @param[in] pairs A list of pairs
* @param[out] triplets List of triplet found in graph
* @return boolean return true if some triplet are found
**/
template <typename IterablePairs, class TTripletContainer>
bool ListTriplets
(
  const IterablePairs & pairs,
  TTripletContainer & triplets
)
{
  triplets.clear();

  // Build an adjacency list corresponding to the edge list
  std::unordered_map<IndexT, std::set<IndexT>> adjacency_list;
  for (const auto & edge : pairs)
  {
    adjacency_list[edge.first].insert(edge.second);
    adjacency_list[edge.second].insert(edge.first);
  }

  std::vector<IndexT> node_candidate_for_triplet;

  // List the pair and find all triplets thanks to the adjacency list
  for (const auto & edge_it : pairs)
  {
    // Find any targeting edge that contains the first and the second node index
    const auto & node1_edges = adjacency_list.find(edge_it.first)->second;
    const auto & node2_edges = adjacency_list.find(edge_it.second)->second;

    // Compute the intersection between the two adjacency lists to find
    //  triplets (it will list the nodes that are connected to the first and
    //  second)
    node_candidate_for_triplet.clear();
    std::set_intersection(node1_edges.cbegin(),
                          node1_edges.cend(),
                          node2_edges.cbegin(),
                          node2_edges.cend(),
                          std::back_inserter(node_candidate_for_triplet));
    // Add a triplet
    for (const auto & node_index_it : node_candidate_for_triplet)
    {
      std::array<IndexT, 3> triplet_indexes {{
        static_cast<IndexT>(edge_it.first),
        static_cast<IndexT>(edge_it.second),
        node_index_it}};
      // sort the triplet indexes as i<j<k (monotonic ascending sorting)
      std::sort(triplet_indexes.begin(), triplet_indexes.end());
      triplets.emplace_back(triplet_indexes[0],
                            triplet_indexes[1],
                            triplet_indexes[2]);
    }

    // Since we have already listed all the triplets than contain this edge.
    // We can now remove the edge from the adjacency graph to reduce the
    // node array size for the next iterations.
    adjacency_list[edge_it.first].erase(edge_it.second);
    adjacency_list[edge_it.second].erase(edge_it.first);
  }
  return ( !triplets.empty() );
}

/**
* @brief Return triplets contained in the graph build from IterablePairs
* @param pairs Graph pairs
* @return List of triplet found in graph
*/
template <typename IterablePairs>
static std::vector<mygraph::Triplet> TripletListing
(
  const IterablePairs & pairs
)
{
  std::vector<mygraph::Triplet> triplets;
  mygraph::ListTriplets( pairs, triplets );
  return triplets;
}

}

/**
 * @Description:  根据vec_triplets 和 relativemotions 得到 motion_triplets
 * @param[in] vec_triplets
 * @param[in] relativemotions
 * @param[out] motion_triplets
 */

bool Tripletcompute(
    std::vector<mygraph::Triplet> &vec_triplets,
    Relativemotions & relativemotions,
    std::vector<Motion_Triplets> & motion_triplets
);


bool Tripletcompute_directed(
    const  std::vector<mygraph::Triplet> &vec_triplets,
    const Relativemotions & relativemotions,
    std::vector<NewMotion_Triplets> & Newmotion_triplets
);


#endif