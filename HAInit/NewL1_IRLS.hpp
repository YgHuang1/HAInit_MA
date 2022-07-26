/*
 * @Description: 新的 L1_IRLS 用于顶点扩展
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-05-30 13:45:16
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-30 15:05:20
 */
#ifndef NewL1_IRLS_hpp
#define NewL1_IRLS_hpp

#include "Reg_definition.hpp"
#include "sophus/se3.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <openMVG/numeric/l1_solver_admm.hpp>

namespace sloveL1_IRLSNew // 用来求解L1 IRLS 
{   
    /// @brief 根据已经固定的顶点 和 未固定的顶点得到J
    void FillMappingMatrix(const set<int> & unfixedNodes,
                           const set<int> & fixedNodes,
                           const Relativemotions & RelMS,
                           Eigen::SparseMatrix<double> &J)
    {
        J.reserve(J.rows()*2);
        Eigen::SparseMatrix<double>::Index i = 0, j = 0;
        for(size_t r = 0; r < RelMS.size(); ++r)
        {
            Relativemotion relM = RelMS[r];
            if(unfixedNodes.find(relM.i) != unfixedNodes.end()) // 如果在没有固定的点中找到了i
            {
                j = 6 * relM.i;
                J.insert(i+0,j+0) = -1.0;
                J.insert(i+1,j+1) = -1.0;
                J.insert(i+2,j+2) = -1.0;
                J.insert(i+3,j+3) = -1.0;
                J.insert(i+4,j+4) = -1.0;
                J.insert(i+5,j+5) = -1.0;
            }
            if(unfixedNodes.find(relM.j) != unfixedNodes.end()) // 如果在没有固定的点中找到了j
            {
                j = 6* relM.j;
                J.insert(i+0,j+0) = 1.0;
                J.insert(i+1,j+1) = 1.0;
                J.insert(i+2,j+2) = 1.0;
                J.insert(i+3,j+3) = 1.0;
                J.insert(i+4,j+4) = 1.0;
                J.insert(i+5,j+5) = 1.0;
            }
            i += 6;
        }
        J.makeCompressed();
    }

    void CorrectMatrix(
        const Eigen::MatrixXd &x,
        Matrix4x4Arr & unfixedMs
    ){
        for(size_t r = 0; r < unfixedMs.size(); ++r)
        {
            Matrix4x4 &Mi = unfixedMs[r];
            const uint32_t i = r;
            const Eigen::Matrix<double, 6,1> eMid = Eigen::Matrix<double,6,1>(x.block<6,1>(6*i,0));
            const Trans eMi;
            Trans updated_Mi = Sophus::SE3d::exp(eMid).matrix();
            Mi = updated_Mi * Mi;
        }
    }

    void FillErrorMatrix(
        const set<int> & unfixedNodes,
        const set<int> & fixedNodes,
        const Relativemotions & RelMS,
        const Matrix4x4Arr& unfixedMs,
        const Matrix4x4Arr& fixedMs,
        Eigen::VectorXd &b
    )
    {
        for(size_t r = 0; r < RelMS.size(); r++)
        {
            const Relativemotion & relM = RelMS[r];
            Trans Mi, Mj;
            if(relM.i < *(fixedNodes.begin()))  /// 说明i未固定 
            {
                 Mi = unfixedMs[relM.i];
            } else // 说明i已经固定
            {
                const int & relMi = distance(fixedNodes.begin(), fixedNodes.find(relM.i)); // 返回其与begin之间的距离
                 Mi = fixedMs[relMi];
            }

            if(relM.j < *(fixedNodes.begin()))  /// 说明j未固定 
            {
                 Mj = unfixedMs[relM.j];
            } else // 说明j已经固定
            {
                const int & relMj = distance(fixedNodes.begin(), fixedNodes.find(relM.j)); // 返回其与begin之间的距离
                 Mj = fixedMs[relMj];
            }
            const Trans eRij = Mj * relM.Tij * Mi.inverse(); // 得到误差矩阵 Mj*Mij*Mi^-1

            Eigen::Matrix<double, 6,1> erij; 
            erij = Sophus::SE3d(eRij).log(); // 得到误差矩阵对应的李代数
            b.block<6,1>(6*r,0) = erij;  
        }
    }

    bool SolveL1MA(
        const Relativemotions & RelMS,
        const set<int> & fixedNodes,
        const set<int> & unfixdeNodes,
        Matrix4x4Arr & unfixedMs,
        const Matrix4x4Arr & fixedMs,
        const Eigen::SparseMatrix<double> &J
    )
    {
        const unsigned nObss = (unsigned)RelMS.size();
        const unsigned nVars = (unsigned)unfixdeNodes.size();
        const unsigned m = nObss*6;
        const unsigned n = nVars * 6;
        const unsigned b_m = nObss*6;
        Eigen::VectorXd x(Eigen::VectorXd::Zero(m)), b(b_m);
        // Current error and the previous one
        double e = std::numeric_limits<double>::max(), ep;
        unsigned iter = 0;
        using namespace openMVG;
        do 
        {
            FillErrorMatrix(unfixdeNodes, fixedNodes, RelMS, unfixedMs, fixedMs, b);
            cout << "第" << iter << "轮迭代误差" << endl;
            L1Solver<Eigen::SparseMatrix<double>>::Options options;
            L1Solver<Eigen::SparseMatrix<double>> l1_solver(options,J);
            l1_solver.Solve(b,&x);

            ep = e; e= x.norm();
            if (ep<e)
            break;
            CorrectMatrix(x, unfixedMs);
        } while (++iter <32 && e > 1e-5 &&(ep-e)/e>1e-2);

        std::cout << "L1RA Converged in " << iter << " iterations." << std::endl;
        return true;
    }
}

#endif