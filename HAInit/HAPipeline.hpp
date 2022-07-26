/* 
 * @Description: HA Init 类
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-26 13:59:32
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 14:42:55
 */
#ifndef HAPipeline_hpp
#define HAPipeline_hpp

#include "Reg_definition.hpp"
#include "Graph_Theory.hpp"
#include "Initmotion.hpp"
#include "ReadWrite_RelMatrixfile.hpp"
#include "Analysis_numeric.hpp"
#include "Edgefitting.hpp"
#include "Global2clus.hpp"
#include "show_model.hpp"
#include "NewL1_IRLS.hpp"
#include "L1_IRLS.hpp"

namespace PointRegistrationPipeline
{   
       ///  @brief: 根据RelMs MSE进行初始化类
    class Initmotion_ByMSEMST
    {   
        public:
            std::set<Pair> MSTset; // 存放初始化图的边
            std::vector<Pair> pairs; // 存放用于基于MSTset进行BFS后的边 用于顺序扩展
            std::vector<double> Rots, trans; // 用于存放初值与真值的比较结果
            bool isSaveInitmotion = false; // 用于判断是否保存初始化后的矩阵

        public:
            Relativemotions RelMs; // 读入的相对矩阵
            Matrix4x4Arr initmotion; // 初始化后的配准初值
            Relativemotions_map RelMsmap; // RelMs的map结构体 用来根据Pair<IndexT,IndexT> 寻找RelM
            Relativemotions_MSEset RelMsMSEset;  // RelMs 的set结构体
            Matrix4x4Arr Groundtruth; // 真实值
            
            string RelMsfile; // 相对矩阵存放地址
            string initmotionfile; // 初始化后的初值存放地址
            int PointNumber; // 点云数量
            

        public:
            ///  读入相对运动集 和 点云帧数
            /// @param RelMsMSEfile 输入的相对运动集文件路径 @param number 点云帧数
            bool Init(const string &RelMsMSEfile,int number)
            {
                this->RelMsfile = RelMsMSEfile;
                this->PointNumber = number;
                this->initmotion.assign(this->PointNumber, Trans::Identity());
                if(Matrix_Read::Read_RelMatrix(this->RelMsfile, this->RelMs, true) == false)
                {
                    std::cerr << "读取相对运动集失败" << endl; return false;
                }
                this->RelMsmap = getmap(this->RelMs);
                this->RelMsMSEset = getMSEset(this->RelMs);
                return true;
            }

            bool Init(const Relativemotions & originRelMs,int number)
            {
                this->RelMsfile = " ";
                this->PointNumber = number;
                this->RelMs = originRelMs;
                this->initmotion.assign(this->PointNumber, Trans::Identity());
                this->RelMsmap = getmap(this->RelMs);
                this->RelMsMSEset = getMSEset(this->RelMs);
                cout << this->RelMsMSEset.size() << endl;
                return true;
            }


            // 基于krushkal 算法得到MST 
            void getMSTBykrushkal(){
                // Graph_Theory::KruskalMST GraphMST;
                // GraphMST.Greatealledge(this->RelMsMSEset);
                // this->MSTset = GraphMST.krushkal(this->PointNumber);
                this->MSTset.clear();
                for(auto it = this->RelMsMSEset.begin(); it != this->RelMsMSEset.end(); it++)
                {
                    this->MSTset.insert({it->i, it->j});
                    cout << it->i << " " << it->j << endl;
                }
                // cout << "MSTset:" << endl;
                // for(auto it : this->MSTset)
                // {
                //     cout << it.first << " " << it.second << endl;
                // }
            }

            /// 根据MST计算得到初值
            /// @param isSaveInitmotion 是否保存初值 @param initmotionfile 保存初始化配准初值的文件路径
            void InitMotion_fromMST(const string &initmotionfile = "", bool isSaveInitmotion = false)
            {
                InitRelMs::InitMotion_fromMST(this->RelMs, this->MSTset, this->initmotion, 0);
                if(isSaveInitmotion){
                    this->initmotionfile = initmotionfile;
                    this->isSaveInitmotion = isSaveInitmotion;
                    Matrix_Write::Write_GlobalMatrix(this->initmotionfile, this->initmotion);
                }   
            }

            // 进行初值与真实的对比
            /// @param truthfile 真实值的文件路径
            void AnalysisComparedTruth(const string &truthfile)
            {
                Matrix_Read::Read_GlobalMatrix(truthfile, this->Groundtruth, this->PointNumber, true);
                Analysis_numeric::AnsAnalysis_betweenTruth(this->initmotion, this->Groundtruth, this->Rots, this->trans);
                Analysis_numeric::cout_RotandTrans(this->Rots, this->trans);
            }

            void AnalysisComparedTruth(const Matrix4x4Arr & truthmat)
            {
                this->Groundtruth = truthmat;
                Analysis_numeric::AnsAnalysis_betweenTruth(this->initmotion, this->Groundtruth, this->Rots, this->trans);
                Analysis_numeric::cout_RotandTrans(this->Rots, this->trans);
            }
    };

    /// @brief 
    class TripletPipeline
    {
        public:
            Relativemotions filtedRelMs; // 经过RelMsthreshold过滤后的RelMs
            Pair_set pairs; // filtedRelMs中包含的边
            std::vector<mygraph::Triplet> vec_triplets; //得到所有的triplets
            std::vector<NewMotion_Triplets> Motiontri; // 计算所有triplets的误差
            std::vector<NewMotion_Triplets> filtedMotiontri; // 存储所有的经过Tripletthreshold过滤的Motiontri
            Pair_set filtedpairs; // 存储filtedMotiontri中的边 即进行最后clu的边
            Relativemotions filtedMotiontri_RelMs; // 存储filtedpairs对应的RelMS
        public:
            bool Init(const Relativemotions &OrigiaRelMs, double RelMsthreshold, const Matrix4x4Arr & Groundtruth)
            {   
                // 得到过滤后的RelMs
                this->filtedRelMs = Edgefilting::edgefilt(OrigiaRelMs, RelMsthreshold, Groundtruth);
                this->pairs = getPairs(this->filtedRelMs); // 根据RelMs得到其中包含的边
                this->vec_triplets = mygraph::TripletListing(this->pairs); // 计算得到存在的所有的边长为3的回环
                // 计算所有边长为3的回环的error
                Tripletcompute_directed(this->vec_triplets, this->filtedRelMs, this->Motiontri); 
                // sort对this->Motiontri排升序
                sort(this->Motiontri.begin(), this->Motiontri.end(), CompareNewMotion_Triplets());       
            }

            /// @brief 根据Tripletthreshold对Triplet进行过滤 得到 filtedpairs 和 filtedMotiontri
            void Tripletfilted(const double & Tripletthreshold)
            {
                this->filtedMotiontri.reserve(this->Motiontri.size());
                for(std::vector<NewMotion_Triplets>::const_iterator iter = Motiontri.begin();
                    iter != Motiontri.end(); iter++ )
                {
                    if(iter->CyclesNorm <= Tripletthreshold)
                    {
                        this->filtedMotiontri.push_back(*iter);
                        this->filtedpairs.insert(iter->IJ); // 将Motiontri中的边存入到filtedpairs中
                        this->filtedpairs.insert(iter->JK);
                        this->filtedpairs.insert(iter->KI);
                    }
                }

                this->filtedMotiontri_RelMs.reserve(filtedpairs.size());
                Relativemotions_map tempmap = getmap(this->filtedRelMs);
                for(auto iter : this->filtedpairs)
                {
                    this->filtedMotiontri_RelMs.push_back(tempmap.at(iter)); // 得到经过Tripletthreshold的RelMs
                }
            }
    };

    /// @brief 根据全局的边 以及 部分全局的顶点 得到这些顶点对应的边
    void getGlobaledgesfromcom(const Pair_set & allpairs,const set<int> nodes, Pair_set & filtedPairsByCom) // 根据得到的最大连通分量得到对应的边
    {
        filtedPairsByCom.clear();
        for(auto iter : allpairs)
        {
            IndexT firstnode = iter.first, secondnode = iter.second;
            if(nodes.find(firstnode) != nodes.end()
                && nodes.find(secondnode) != nodes.end())
            {
                filtedPairsByCom.insert({firstnode, secondnode});
            }
        }
    }
    /// @brief 根据全局的边 以及 部分全局的顶点的集合 得到这些顶点对应的边的集合 用于扩展
    void getAllGlobaledgesfromcoms(const Pair_set & allpairs, const vector<set<int>> & nodess, vector<Pair_set>  & filtedPairsByComs)
    {
        filtedPairsByComs.clear(); filtedPairsByComs.reserve(nodess.size()); // filtedPairsByComs 的个数与nodess的个数一样
        for(vector<set<int>>::const_iterator nodes_iter = nodess.begin(); nodes_iter != nodess.end(); nodes_iter++)
        {
            Pair_set newpairset;
            getGlobaledgesfromcom(allpairs, *nodes_iter, newpairset);
            filtedPairsByComs.push_back(newpairset);
        }
    }
    
    /// @brief 根据已经确定的全局顶点 以及 扩展过程中需要确定的顶点的集合(包含了部分全局顶点 即基顶点) 分离得到需要确定的顶点 和父顶点
    /// @param[in] allnodes_hadinit 所有已经确定了初值的顶点 即基顶点的集合
    /// @param[in] allnodes_needdepart 所有需要分离顶点的vector
    /// @param[out] vec_unfixdnodes 输出: 得到所有需要确定初值的顶点 不包含基顶点
    /// @param[out] vec_fixednodes  输出: 对应于 nodes_needinit 得到其对应的基顶点
    template<typename T>
    void getNodes_needInit(const set<T> & allnodes_hadinit, const vector<set<int>> & allnodes_needdepart, vector<set<int>> & vec_unfixdnodes, 
                            vector<set<int>> & vec_fixednodes )
    {
        vec_unfixdnodes.clear(); vec_fixednodes.clear();
        for(vector<set<int>>::const_iterator iter = allnodes_needdepart.begin(); iter != allnodes_needdepart.end(); iter++)
        {
            set<int> fixednodes, unfixednodes;
            for(auto it : *iter)
            {
                if(allnodes_hadinit.find(it) != allnodes_hadinit.end()) // 说明it已经初始化了
                {
                    fixednodes.insert(it);} 
                else{
                    unfixednodes.insert(it);}}
            vec_unfixdnodes.push_back(unfixednodes), vec_fixednodes.push_back(fixednodes);
        }
    }
    

    class TripletPipelineForHA // 用于HA Init 计算所有TriPlet的结果
    {
        public:
            Relativemotions filtedRelMs; // 经过RelMsthreshold过滤后的RelMs
            Pair_set pairs; // filtedRelMs中包含的边
            std::vector<mygraph::Triplet> vec_triplets; //得到所有的triplets
            std::vector<Motion_Triplets> Motiontri; // 计算所有triplets的误差

        public:
            bool Init(const Relativemotions &OrigiaRelMs)
            {   
                // 得到过滤后的RelMs
                this->filtedRelMs = OrigiaRelMs;
                this->pairs = getPairs(this->filtedRelMs); // 根据RelMs得到其中包含的边
                this->vec_triplets = mygraph::TripletListing(this->pairs); // 计算得到存在的所有的边长为3的回环
                // 计算所有边长为3的回环的error
                Tripletcompute(this->vec_triplets, this->filtedRelMs, this->Motiontri);
                // Tripletcompute_directed(this->vec_triplets, this->filtedRelMs, this->Motiontri); 
                // sort对this->Motiontri排升序
                sort(this->Motiontri.begin(), this->Motiontri.end(), CompareMotion_Triplets());       
            }
    };


    /// @brief 根据MSE进行HA初始化
    class HAInit_ByMSE // HA Init 类
    {
        /// @brief 用于存放第一轮的顶点和边 初始化的结构体 
        struct HAstruct 
        {
            std::set<Pair> MSTset; // 存放用于初始化的边
            std::vector<Pair> pairs; // 存放用于基于MSTset进行BFS后的边 用于顺序扩展
            std::set<int> NodesInit; // 存放初始化的顶点 局部
            std::set<int> GlobalNodesInit; // 存放对应于全局的顶点 set会自动排序

            Relativemotions allRelMs; //  相对矩阵全局 没有经过最大连通分量
            Pair_set allpairs; // 没有经过最大连通分量前的边

            Pair_set filtedPairsByCom; // 经过最大连通分量判断后的边 先全局 后局部
            Relativemotions filtedRelMs;  // 相对矩阵 先全局 后局部 经过最大连通分量

            Matrix4x4Arr initmotion; // 初始化后的配准初值 局部

            public:
                void getGlobaledgesBycom() // 根据得到的最大连通分量得到对应的边
                {
                    this->filtedPairsByCom.clear();
                    for(auto iter : this->allpairs)
                    {
                        IndexT firstnode = iter.first, secondnode = iter.second;
                        if(this->GlobalNodesInit.find(firstnode) != this->GlobalNodesInit.end()
                         && this->GlobalNodesInit.find(secondnode) != this->GlobalNodesInit.end())
                        {
                            this->filtedPairsByCom.insert({firstnode, secondnode});
                        }
                    }
                }

                void InitByGraph()
                {
                    Graph_Theory::UndirectedGraph_Find_MaxComp MaxComp;
                    cout << "寻找最大连通分量" << endl;
                    MaxComp.Find_Max_components(this->allpairs, this->GlobalNodesInit); //  得到最大连通分量 
                    getGlobaledgesBycom(); // 得到最大连通分量对应的边
                    this->filtedRelMs.reserve(this->filtedPairsByCom.size());
                    
                    Relativemotions_map tempmap = getmap(this->allRelMs);
                    for(auto iter : this->filtedPairsByCom)
                    {
                        this->filtedRelMs.push_back(tempmap.at(iter)); // 得到了对应最大连通分量边的RelMS
                    }
                    // 更新nodes filtedRelMs  全局->局部
                    Global2clus::convertRelMs(this->GlobalNodesInit, this->filtedRelMs); // 更改filtedRelMs
                    for(int i = 0; i < this->GlobalNodesInit.size(); i++)
                    {
                        this->NodesInit.insert(i); // 更新nodes 
                    }

                    for(auto iter : this->filtedRelMs)
                    {
                        cout << iter.i << " " << iter.j << endl;
                    }

                    // 进行初始化 
                    Matrix4x4Arr tempinitmotion;
                    PointRegistrationPipeline::Initmotion_ByMSEMST tempInit; // 创建进行初始化初值的类
                    
                    tempInit.Init(this->filtedRelMs, this->GlobalNodesInit.size());
                    tempInit.getMSTBykrushkal();
                    tempInit.InitMotion_fromMST();
                    tempinitmotion = tempInit.initmotion; // 将初始化后的顶点输入tempinitmotion

                    
                    const unsigned nobss = (unsigned)this->filtedRelMs.size();
                    const unsigned nVars = (unsigned)this->GlobalNodesInit.size() - 1;
                    const unsigned m = nobss*6;
                    const unsigned n = nVars*6;
                    Eigen::SparseMatrix<double> A(m,n);
                    sloveL1_IRLS::FillMappingMatrix(this->filtedRelMs, (uint32_t)0, A);
                    sloveL1_IRLS::SolveL1RA(this->filtedRelMs, tempinitmotion, A, 0); // 进行L1优化
                    this->initmotion = tempinitmotion;
                }
        };
        
        /// @brief 用于进行顶点扩展的结构体
        struct expandstruct
        {
            set<int> unfixednode_Globalset; // 需要初始化的子顶点 全局
            set<int> fixednode_Globalset; // 基顶点 全局
            set<int> unfixednode_localset; // 需要初始化的子顶点 局部
            set<int> fixednode_localset; // 基顶点 局部

            map<int, Trans> fiexdNodes_map; // 所有的基顶点 及其对应的map结构体 全局
            
            Relativemotions RelMs; // 存放所有的edges 全局 -> 局部
            
            Matrix4x4Arr unfixednode_initmotion; // 需要进行初始化的顶点对应初始化后的初值
            Matrix4x4Arr fixedNode_initmotion; // 已经进行了初始化的顶点

            set<int> Allnodes; // 包含所有的顶点 全局
            Pair_set Edges_L1; // 包含所有进行L1的边 全局
            public:
                void localtoGlobal()
                {
                    vector<int> nodes_fatherwithson;
                    nodes_fatherwithson.reserve(this->unfixednode_Globalset.size() + this->fixednode_Globalset.size());
                    for(set<int>::const_iterator iter = unfixednode_Globalset.begin(); iter != unfixednode_Globalset.end();
                        iter++ )
                    {
                        nodes_fatherwithson.push_back(*iter);
                    }
                    for(set<int>::const_iterator iter = fixednode_Globalset.begin(); iter != fixednode_Globalset.end();
                        iter++ )
                    {
                        nodes_fatherwithson.push_back(*iter);
                    }

                    cout << "正在将RelMs从全局转到局部" << endl;
                    Global2clus::convertRelMs(nodes_fatherwithson, this->RelMs);
                    // for(auto iter : this->RelMs)
                    // {
                    //     cout << iter.i << " " << iter.j << endl;
                    //     cout << iter.Tij << endl;
                    // }

                    /// 将fixedNode_initmotion中的值赋值
                    this->fixedNode_initmotion.reserve(this->fixednode_Globalset.size());
                    for(auto iter : this->fixednode_Globalset)
                    {
                        this->fixedNode_initmotion.push_back(this->fiexdNodes_map.at(iter));
                    }
                
                    
                    for(int i = 0; i < this->unfixednode_Globalset.size(); i++) // 将全局的子顶点转换到局部的子顶点下
                    {
                        this->unfixednode_localset.insert(i);
                    }
                    for(int i = 0; i < this->fixednode_Globalset.size(); i++) // 将全局的基顶点转换到局部的基顶点下
                    {
                        this->fixednode_localset.insert(i + this->unfixednode_Globalset.size());
                    }
                    this->unfixednode_initmotion.reserve(this->unfixednode_Globalset.size());

                    /// 输出fixednode_localset unfixednode_localset
                    // cout << "fixednode_localset:" << endl;
                    // for(auto it : this->fixednode_localset)
                    // {
                    //     cout << it << endl;
                    // } 
                    // cout << "unfixednode_localset:" << endl;
                    // for(auto it : this->unfixednode_localset)
                    // {
                    //     cout << it << endl;
                    // }
                }

                

                void sloveL1() // 求解之前需要将所有的nodes转换到局部坐标下
                {
                    /// 首先对初值进行初始化
               
                    Matrix4x4Arr tempMotionArr;
                    tempMotionArr.assign(this->Allnodes.size(), Trans::Identity());
                    // sloveL1_IRLS::InitMotionsMAT(this->RelMs, tempMotionArr, 0); 
                    PointRegistrationPipeline::Initmotion_ByMSEMST tempInit; // 创建进行初始化初值的类
                    tempInit.Init(this->RelMs, this->Allnodes.size());
                    tempInit.getMSTBykrushkal();
                    tempInit.InitMotion_fromMST();
                    tempMotionArr = tempInit.initmotion; // 将初始化后的顶点输入至tempMotionArr
                    
                    // for(auto it : tempMotionArr)
                    // {
                    //     cout << it << endl;
                    // }

                    Trans endTrans = (*(tempMotionArr.end() - 1)).inverse(); // 最后一个元素求逆
                    Trans fatherTrans = fiexdNodes_map.find(*(fixednode_Globalset.rbegin()))->second; // 基顶点的Trans
                    cout << "基顶点为" << *fixednode_Globalset.rbegin() << endl;
                    for(int i = 0; i < tempMotionArr.size(); i++)
                    {
                        tempMotionArr[i] = fatherTrans* endTrans * tempMotionArr[i] ; // 全部转到基顶点
                        cout << tempMotionArr[i] << endl;
                    }
                    this->unfixednode_initmotion.reserve(this->unfixednode_Globalset.size());
                    for(int i = 0; i < this->unfixednode_Globalset.size(); i++) // 获得经过最小生成树遍历的初值
                    {
                        this->unfixednode_initmotion.push_back(tempMotionArr.at(i));
                    }

                    const unsigned nObss = (unsigned)this->RelMs.size();
                    const unsigned nVars = (unsigned)this->unfixednode_localset.size();
                    const unsigned m = nObss*6;
                    const unsigned n = nVars * 6;
                    Eigen::SparseMatrix<double> J(m,n);
                    sloveL1_IRLSNew::FillMappingMatrix(this->unfixednode_localset, this->fixednode_localset,
                                                        this->RelMs, J);
                    sloveL1_IRLSNew::SolveL1MA(this->RelMs, this->fixednode_localset, this->unfixednode_localset,
                                               this->unfixednode_initmotion, this->fixedNode_initmotion, J);
                    
                }
            
        };

        public:
            std::vector<HAstruct> HAstructs; // 用于存放每轮扩展的边等信息
            std::set<int> AllNodesInit; // 用于存放已经确定初值的顶点
            Relativemotions filtedRelMs; // 经过RelMsthreshold过滤后的RelMs
            Relativemotions_map filtedRelMs_map; // filtedRelMs 对应的map 用于查找RelMS
            Matrix4x4Arr GlobalInitMotion; // 存放所有的已经初始化的Global Motion 开始全部初始化为I
            int Number;                    //  点云帧数
            TripletPipelineForHA HATriplet; // 类 计算并保存所有的Triplet 以及 cyclenorm
            double threshold, thresholdend; // 每轮扩展增加的阈值 每轮扩展时最大的阈值
            Matrix4x4Arr truthmat; // 真实值 用于比较
        
        public:
            void Init(const string &RelMsMSEfile,const string &truthfile,int number) // 加载所有的RelMS truth 已经总的点云帧数
            {
                if(Matrix_Read::Read_RelMatrix(RelMsMSEfile, this->filtedRelMs, false) == false)
                {
                    std::cerr << "读取相对运动集失败" << endl; 
                }

                this->Number = number;
                this->GlobalInitMotion.assign(this->Number, Trans::Identity()); // 假定所有的初值为单位阵
                this->HATriplet.Init(this->filtedRelMs);
                // this->filtedRelMs = HATriplet.filtedRelMs;
                this->filtedRelMs_map = getmap(this->filtedRelMs);
                this->HAstructs.reserve(10); //假定10轮完成初始化
                Matrix_Read::Read_GlobalMatrix(truthfile, this->truthmat, this->Number, true);
            }

            void Init(const Relativemotions & RelMS,const string &truthfile,int number) // 加载所有的RelMS truth 已经总的点云帧数
            {
                this->filtedRelMs = RelMS;
                this->Number = number;
                this->GlobalInitMotion.assign(this->Number, Trans::Identity()); // 假定所有的初值为单位阵
                this->HATriplet.Init(this->filtedRelMs);
                // this->filtedRelMs = HATriplet.filtedRelMs;
                this->filtedRelMs_map = getmap(this->filtedRelMs);
                this->HAstructs.reserve(10); //假定10轮完成初始化
                Matrix_Read::Read_GlobalMatrix(truthfile, this->truthmat, this->Number, true);
            }
            
            /// @brief 进行第一轮的初始化
            void Initnodes(const double & threshold = 0.01)
            {
                // 每次cycle的增量为0.0003
                // 初始化
                this->threshold = threshold; // 每轮增加阈值
                this->thresholdend = threshold; // 开始阈值
                cout << "每轮扩展的阈值为" << this->threshold << endl;
                cout << "正在进行顶点的初始化" << endl;
                HAstruct InitHAstruct;
                { // 进行顶点的初始化
                    InitHAstruct.initmotion.reserve(this->Number);
                    InitHAstruct.allRelMs.reserve(this->HATriplet.filtedRelMs.size());
                    Graph_Theory::UndirectedGraph_Find_MaxComp MaxComp; set<int> MaxTemp;
                    do
                    {
                       for(std::vector<Motion_Triplets>::const_iterator iter = this->HATriplet.Motiontri.begin();
                            iter != this->HATriplet.Motiontri.end(); iter++)
                        {
                            if(iter->CyclesNorm < this->thresholdend && iter->CyclesNorm >=(this->thresholdend-this->threshold))
                            {
                                InitHAstruct.allpairs.insert({iter->i, iter->j});
                                InitHAstruct.allpairs.insert({iter->j, iter->k});
                                InitHAstruct.allpairs.insert({iter->i, iter->k});
                                cout << iter->i << " " << iter->j << " " << iter->k << " " << iter->CyclesNorm << endl;
                            }
                        }
                        this->thresholdend =  this->thresholdend + this->threshold; // 更新区间结束阈值
                        // MaxComp.Find_Max_components(InitHAstruct.allpairs, MaxTemp);
                    } while ( InitHAstruct.allpairs.size() <= 2); // 直到不为0

                    InitHAstruct.allRelMs.reserve(InitHAstruct.allpairs.size());
                    for(auto iter : InitHAstruct.allpairs)
                    {
                        InitHAstruct.allRelMs.push_back(this->filtedRelMs_map.at(iter)); // 得到对应allpairs的allRelMs
                    }
                }              
                InitHAstruct.InitByGraph();      
                this->HAstructs.push_back(InitHAstruct);    
                for(auto iter : InitHAstruct.GlobalNodesInit)    
                {
                    this->AllNodesInit.insert(iter);  // 得到了所有已经初始化的顶点
                } 
                for(int i = 0; i < InitHAstruct.initmotion.size(); i++) // 将局部处值重新赋予到全局中
                {
                    std::set<int>::iterator iter = this->AllNodesInit.begin();
                    std::advance(iter, i);
                    cout << *iter << endl;
                    this->GlobalInitMotion.at(*iter) = InitHAstruct.initmotion.at(i);
                }

                cout << "已经初始化的顶点为:" << endl;
                for(auto it : this->AllNodesInit)
                {
                    cout << it << " ";
                }
                cout << endl;
                // for(auto it : this->GlobalInitMotion)
                // {
                //     cout << it << endl;
                // }

            }

            /// @brief 顶点的扩展
            void expandNodes()
            {   
            
                while (this->AllNodesInit.size() < this->Number) 
                {   
                    // cout << "正在进行顶点的扩展" << endl;
                    Pair_set newpair;
                    do
                    {
                        for(std::vector<Motion_Triplets>::const_iterator iter = this->HATriplet.Motiontri.begin();
                                iter != this->HATriplet.Motiontri.end(); iter++)
                        {
                            if(iter->CyclesNorm < this->thresholdend ) // 阈值区间搜索
                            {   
                                // 获得该阈值区间内 并且没有进行L1优化的边
                                if(!(this->AllNodesInit.find(iter->i) != this->AllNodesInit.end() && 
                                this->AllNodesInit.find(iter->j) != this->AllNodesInit.end())) // 找到了ij说明 顶点已经初始化
                                {
                                    newpair.insert({iter->i, iter->j});
                                    cout << "扩展的边为:" << iter->i << " " << iter->j << endl;
                                }
                                if(!(this->AllNodesInit.find(iter->i) != this->AllNodesInit.end() && 
                                this->AllNodesInit.find(iter->k) != this->AllNodesInit.end())) // 找到了ij说明 顶点已经初始化
                                {
                                    newpair.insert({iter->i, iter->k});
                                    cout << "扩展的边为:" << iter->i << " " << iter->k << endl;
                                }
                                if(!(this->AllNodesInit.find(iter->k) != this->AllNodesInit.end() && 
                                this->AllNodesInit.find(iter->j) != this->AllNodesInit.end())) // 找到了ij说明 顶点已经初始化
                                {
                                    newpair.insert({iter->j, iter->k});
                                    cout << "扩展的边为:" << iter->j << " " << iter->k << endl;
                                }
                                // cout << iter->i << " " << iter->j << " " << iter->k << " " << iter->CyclesNorm << endl;
                            }
                        } 
                        this->thresholdend = this->thresholdend + this->threshold; // 更新阈值区间
                    } while (newpair.size() <= 1); // 当newpair中大于1时 退出阈值区间搜索
                
                    for(auto iter : newpair)
                    {
                        cout << iter.first << " " << iter.second << endl;
                    }

                    std::vector<set<int>> nodesss; // **** 保存着多个连通分量对应的顶点 注意: 该顶点的集合中包含着已经初始化的顶点 后续需要将这些顶点进行分离
                    Graph_Theory::UndirectedGraph_Find_MaxComp MaxComp;
                    MaxComp.Find_All_components(newpair, nodesss); // 找到所有的连通分量

                    std::vector<Pair_set> pairss; // **** 保存着多个连通分量对应的pairs
                    getAllGlobaledgesfromcoms(newpair, nodesss, pairss); // 得到对应nodesss的所有的pairs
                    
                    vector<set<int>> vec_unfixednodes, vec_fixednodes; // ****
                    getNodes_needInit(this->AllNodesInit, nodesss, vec_unfixednodes, vec_fixednodes); // 得到基顶点 和 子顶点
                    if(vec_fixednodes.empty()) // 如果为空 说明没有基顶点
                    {
                        continue;
                    }

                    for(int i = 0; i <vec_unfixednodes.size(); i++)
                    {
                        expandstruct newexpand;
                        newexpand.unfixednode_Globalset = vec_unfixednodes.at(i); // 得到所有的基顶点
                        cout << "扩展的顶点为" << endl;
                        for(auto iter : newexpand.unfixednode_Globalset)
                        {
                            cout << iter << " ";
                        }
                        cout << endl;
                        newexpand.Edges_L1 = pairss.at(i);  // 得到所有要进行L1的边
                        newexpand.Allnodes = nodesss.at(i); // 得到所有的顶点
                        newexpand.fixednode_Globalset = vec_fixednodes.at(i); // 得到所有的父顶点
                        if(newexpand.fixednode_Globalset.empty()) // 如果为空 说明当前不存在父顶点
                        {
                            continue;
                        }
                        newexpand.RelMs.reserve(newexpand.Edges_L1.size()); 
                        for(auto iter : newexpand.Edges_L1)
                        {
                            newexpand.RelMs.push_back(this->filtedRelMs_map.at(iter)); // 得到Edges_L1中边对应的Relativemotions
                        }

                        for(auto iter : newexpand.RelMs)
                        {
                            cout << iter.i << " " << iter.j << endl;
                        }
                        cout  <<  "基顶点" << endl;
                        for(auto iter : newexpand.fixednode_Globalset)
                        {
                            cout << iter << endl;
                        }
                        cout << "子顶点" << endl;
                        for(auto iter : newexpand.unfixednode_Globalset)
                        {
                            cout << iter << endl;
                        }

                        for(auto iter : newexpand.Edges_L1)
                        {
                            cout << iter.first << " " << iter.second << endl;
                        }

                        for(auto iter : vec_fixednodes.at(i))
                        {
                            newexpand.fiexdNodes_map.emplace(iter, this->GlobalInitMotion.at(iter));
                            // cout << "已经固定的顶点为:" << iter << endl << this->GlobalInitMotion.at(iter) << endl;
                        }
                        newexpand.localtoGlobal();
                        newexpand.sloveL1(); // 求解sloveL1
                        // 更新全局顶点 求交集
                        auto itend = set_union(newexpand.unfixednode_Globalset.begin(), newexpand.unfixednode_Globalset.end(),
                                                    this->AllNodesInit.begin(), this->AllNodesInit.end(), 
                                                    inserter(this->AllNodesInit, this->AllNodesInit.begin())); 
                        // 更新全局motion 
                        for(int i = 0; i < newexpand.unfixednode_Globalset.size(); i++)
                        {
                            set<int>::const_iterator setbegin = newexpand.unfixednode_Globalset.begin();
                            advance(setbegin, i); // 向前相加
                            this->GlobalInitMotion.at(*setbegin) = newexpand.unfixednode_initmotion.at(i);
                        }
                    }
                }
                cout << "初始化后的初值为:" << endl;
                Trans ZeroMain = this->GlobalInitMotion[0].inverse();
                for(int i = 0; i < this->GlobalInitMotion.size(); i++)
                {
                    this->GlobalInitMotion[i] = this->GlobalInitMotion[i] * ZeroMain;
                    cout << this->GlobalInitMotion[i] << endl;
                }
                
            }

            void comparewithTruth()
            {   
                vector<double> Rots, trans;
                Analysis_numeric::AnsAnalysis_betweenTruth(this->GlobalInitMotion, this->truthmat, Rots, trans);
                Analysis_numeric::cout_RotandTrans(Rots, trans);
            }
    };
};

#endif
