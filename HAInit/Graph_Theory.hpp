/*
 * @Description: 图论相关运算
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-04-06 21:47:50
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-30 09:31:29
 */
#ifndef Graph_Theory_hpp
#define Graph_Theory_hpp

#include "Reg_definition.hpp"
#include <vector>
#include <set>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/smart_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <queue>
#include <math.h>
#include <iostream>
#include <memory.h>
#include <iostream>
#include <algorithm>


using namespace std;

namespace Graph_Theory
{    
    class UndirectedGraph_Find_MaxComp{

        private:
            vector<int> G[100]; // 
            int seen[100]; // 判断节点是否被访问
            int Count_componts; // 连通块个数
        
        public:

            UndirectedGraph_Find_MaxComp(){
                for(int i = 0; i < sizeof(seen)/ sizeof(seen[0]); i++){
                    seen[i] = 0;
                }
            }


            // 深度搜索
            void dfs(int u, set<int> &nodes){
                seen[u] = 1;
                for(int i = 0; i < G[u].size(); i++){
                    if(!seen[G[u][i]]){
                        // cout << G[u][i] << "\t";
                        nodes.insert(G[u][i]-1);
                        dfs(G[u][i], nodes);
                    }
                }
            }
            // 找到边集合中的最大连通块， 将最大连通块中的顶点输出到Pair_set
            void Find_Max_components(const set<Pair> & SrcPair_set, set<int> & Pair_set){
            
                set<int> All_nodes;
                std::vector<set<int>> nodess;

                for(auto iter : SrcPair_set){
                    G[iter.first+1].push_back(iter.second+1);
                    G[iter.second+1].push_back(iter.first+1);
                    All_nodes.insert(iter.first+1);
                    All_nodes.insert(iter.second+1);
                }      
                
                auto iter = All_nodes.end();
                iter--;
                int n = *iter;
               
                cout <<" nodes : " << n << endl;
                // int m = SrcPair_set.size(); // 边的个数 即SrcPair_set的大小
                
                for(int i = 1; i <= n; i++){
                    set<int> nodes;
                    nodes.insert(i-1);

                    if(!seen[i])
                    {  
                        //  cout << i << endl;
                        dfs(i, nodes);
                        if(nodes.size() > 1){
                            nodess.push_back(nodes);
                        }
                    }
                }

                int maxnodes = nodess[0].size();
                set<int> finalNodes;
                for(int i = 0; i < nodess.size(); i++)
                {
                    if(nodess[i].size() >= maxnodes){
                        maxnodes = nodess[i].size();
                        finalNodes = nodess[i];
                    }
                }

                Pair_set.clear();
                Pair_set.swap(finalNodes);
            }

            // 找到所有的连通分量
            void Find_All_components(const set<Pair> & SrcPair_set, std::vector<set<int>> & allnodess){
            
                set<int> All_nodes;
                std::vector<set<int>> nodess;

                for(auto iter : SrcPair_set){
                    G[iter.first+1].push_back(iter.second+1);
                    G[iter.second+1].push_back(iter.first+1);
                    All_nodes.insert(iter.first+1);
                    All_nodes.insert(iter.second+1);
                }      
                
                auto iter = All_nodes.end();
                iter--;
                int n = *iter;
               
                cout <<" nodes : " << n << endl;
                // int m = SrcPair_set.size(); // 边的个数 即SrcPair_set的大小
                
                for(int i = 1; i <= n; i++){
                    set<int> nodes;
                    nodes.insert(i-1);

                    if(!seen[i])
                    {  
                        //  cout << i << endl;
                        dfs(i, nodes);
                        if(nodes.size() > 1){
                            nodess.push_back(nodes);
                        }
                    }
                }
                allnodess.clear();
                allnodess = std::move(nodess);
                nodess.clear();  
            }
    };


    // Pair_set 中包含了最大连通块中的顶点 srcPair_set 包含了所有pair 即边  edgesPair 得到最大连通块的所有边
    void GetEdges(const std::set<Pair> &srcPair_set, const std::set<int> &Pair_set, std::set<Pair> &edgesPair)
    {
        edgesPair.clear();
        for(int i = 0; i < Pair_set.size(); i++){
            set<int>::iterator iter = Pair_set.begin();
            std::advance(iter, i);
            int firstNode = *iter;
            for(int j = i +1; j < Pair_set.size(); j++)
            {
                std::advance(iter, 1);
                int secondNode = *iter;
                if(srcPair_set.find({firstNode, secondNode}) != srcPair_set.end())
                {
                    edgesPair.insert({firstNode, secondNode});
                }
            }

        }
    }

    /**
     * @brief: 根据输入的边的集合与边的集合中包含的顶点 输出最大度数点
     * @param[in] EdgesPair 边的set集合 
     * @param[in] Nodes 顶点的set集合 
     * @author: hyg
     */    
    int GetMainviewr(const std::set<Pair> &EdgesPair, const std::set<int> &Nodes)
    {   
        int Max_nodes = *Nodes.rbegin() + 1;
        cout << "Max_nodes : " <<  Max_nodes << endl;

        Eigen::SparseMatrix<double, Eigen::ColMajor> adjMat(Max_nodes, Max_nodes);
        for(auto iter : EdgesPair){
            adjMat.insert(iter.first, iter.second) = 1;
        }

        Eigen::SparseMatrix<double, Eigen::ColMajor> adjMatTranspose = adjMat.transpose();
        Eigen::MatrixXd AdjMat = adjMat + adjMatTranspose;
        Eigen::MatrixXd CoDeg = (AdjMat * AdjMat).array() * AdjMat.array();

        int temp = 0, Mainviewer;

        for(int j = 0; j < Max_nodes; j++){
            int degree = 0;
            for(int i = 0; i < Max_nodes; i++){
                if(AdjMat(i, j) == 1)
                {
                    degree ++;
                }
            }
            
            if(degree >= temp){
                temp = degree;
                Mainviewer = j;
            }

            cout << "顶点" << j << "的度数为" << degree << endl;
        }

        return Mainviewer;
    }

    ///  进行BFS 如果失败返回false
    ///  @param[in] SrcPair_set -> 边的集合 @param[out] Pairs -> 输出的BFS顺序对
    ///  @param[in] Mainviewer ->主视角
    bool BFS(const set<Pair> &SrcPair_set, vector<Pair> &Pairs, const int &Mainviewer){
        
        if(SrcPair_set.empty())
        {   
            std::cerr << "输入为空" << endl;
            return false; 
        }

        int graph[100][100], color[100], dist[100];
        memset(graph, 0, sizeof(graph));
        memset(color, 0, sizeof(color));
        memset(dist, 0, sizeof(dist));
        const int WHITE = 0, GRAY = 1, BLACK = 2;
        
        std::set<int> Allnodes;
        Allnodes.clear();

        for(auto iter : SrcPair_set){
            graph[iter.first][iter.second] = 1;
            graph[iter.second][iter.first] = 1;
            Allnodes.insert(iter.first);
            Allnodes.insert(iter.second);
        }

        auto iter = Allnodes.end();
        iter--;
        int nodes = *iter;
        // cout << nodes << endl;
        Pairs.clear();
        Pairs.reserve(nodes + 1);
        
        // run BFS
        std::queue<int> q;
        q.push(Mainviewer);
        dist[Mainviewer] = 1;

        do{
            int u = q.front();
            q.pop();

            for (int i = 0; i < nodes + 1 ; i++){
                if((graph[u][i]==1) && (color[i] == WHITE )){
                    q.push(i);
                    color[i] = GRAY;
                    dist[i] = dist[u] + 1;
                    Pairs.push_back({u, i });
                }
            }
            color[u] = BLACK;
        }while(!q.empty());

        if(Pairs.empty())
        {
            return false;
        }
        return true;
    }


    class FindMST_directedGraph
    {   
       
        public:
            void dijkstra(const Relativemotions_MSEset & RelMSE_set, int number,
                            int startnode, Pair_set &Edges)
            {
                double graph[number][number];
                double distances[number];
                int father[number];
                bool visit[number];

                memset(graph, 0, sizeof(graph));
                memset(distances, 0, sizeof(distances));
                for(auto iter : RelMSE_set)
                {
                    graph[iter.i][iter.j] = iter.TrimmedMSE;
                };
                for(int i  = 0; i < number; i++)
                {
                    for(int j = 0; j < number; j++)
                    {
                        cout << graph[i][j] << " ";
                    }
                    cout << endl;
                }

                priority_queue<pair<double, int> > queue;
                pair <double, int> nodotmp;
                int i, j;
                for (int i = 0; i < number; i++) {
                    distances[i] = 100;
                    father[i] = -1;
                    visit[i] = false;
                }
                distances[startnode] = 0.0;
                queue.push(pair <double, int>(distances[startnode], startnode));

                while (!queue.empty()) {
                    nodotmp = queue.top();
                    queue.pop();
                    i = nodotmp.second;
                    if (!visit[i]) {
                        visit[i] = true;
                        for (j = 0; j < number; j++)
                            if (!visit[j] && graph[i][j] > 0 && distances[i] + graph[i][j] < distances[j]) {
                                distances[j] = distances[i] + graph[i][j];
                                father[j] = i;
                                cout << i << " , " << j << endl;
                                queue.push(pair <double, int>(-distances[j], j));
                            }
                    }
                }

                Edges.clear();
                for (int i = 0; i < number; i++)
                {   
                    if (father[i] == -1)
                    {
                        continue;
                    }
                    Edges.insert(make_pair(father[i], i));
                }
            };
    };

    class KruskalMST
    {
    public:
        struct Edge {
            // the two endpoints and the weight
            int u, v;
            double w;
            // a comparator that stores by least wight
            bool operator<(const Edge& p) const
            {
                return w < p.w;
            }
        };
        int pr[100];
        // int n, e; // nodes edges
        vector<Edge> edges;
        std::set<Pair> MST;
        double ShortestTree;
        
    public:
        
        KruskalMST()
        {
            memset(this->pr, 0, sizeof(pr));
            this->edges.reserve(100);
            // this->n = 0;
        }

        void Greatoneedge(int startnode, int endnode, double weight) // 生成一条边
        {
            Edge edge;
            edge.u = startnode, edge.v = endnode, edge.w = weight;
            this->edges.push_back(edge);
            // edge.v = startnode, edge.u = endnode, edge.w = weight;
            // this->edges.push_back(edge);
        }

        void Greatealledge(const Relativemotions_MSEset & RelMsSet) // 根据Set生成所有的边
        {   
            cout << "正在根据RelMsSet生成图的边" << endl;
            for(auto iter : RelMsSet)
            {
                Greatoneedge(iter.i, iter.j, iter.TrimmedMSE);
            }
        }

        int findset(int r)
        {
            if (this->pr[r] == r) return r;
            return findset(this->pr[r]);
        }

        void makeset(int n)
        {
            for (int i = 0; i < n; i++) this->pr[i] = i;
        }

        std::set<Pair> krushkal(int n) // 进行krushkal 返回MST
        {   
            cout << "正在根据Krushkal算法计算图的最小生成树" << endl;
            // sort the edges
            sort(edges.begin(), edges.end());

            makeset(n); // create set

            int count = 0;
            double sum = 0;
            for (int i = 0; i < (int)edges.size(); i++) {
                int u = findset(edges[i].u);
                int v = findset(edges[i].v);
                if (u != v)
                {
                    pr[u] = v;
                    count++;
                    sum += edges[i].w;
                    cout << edges[i].u << "," << edges[i].v << " , " << edges[i].w << endl;
                    this->MST.insert(make_pair(edges[i].u, edges[i].v));
                    if (count == n - 1) break;
                }
            }
            this->ShortestTree = sum;
            return this->MST;
        }

    };
    
}


#endif