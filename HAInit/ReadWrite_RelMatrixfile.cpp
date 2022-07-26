/*
 * @Description: 该程序用于读和写相对运动矩阵文件
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-02-23 19:22:45
 * @LastEditors: hyg
 * @LastEditTime: 2022-05-30 09:17:23
 */

#include "ReadWrite_RelMatrixfile.hpp"

using namespace pcl::console;
using namespace std;

namespace Matrix_Read
{
     
    bool Read_RelMatrix(const string & InputRelMotions, Relativemotions & RelMs, bool isReadMSE)
    {
        ifstream fin(InputRelMotions, ios::in); // 判断是否可以打开
        if (!fin.is_open()){
            cout << "矩阵读取失败" << endl;
            return false;
        }

        // 将相对运动的索引和相对运动矩阵读入到RelRs中
        string line_info,input_result;

        RelMs.clear();
        RelMs.reserve(1000);
        cout << "正在读取相对运动矩阵" << endl;
        {
            int mm = 0, nn =0;
            double TrimmedMSEij = 0.0;
            if(fin){
            for(int i = 0; getline(fin, line_info); i++ ){
                stringstream input(line_info);
                Eigen::Matrix4d temp_motion;
                for (int j = 0; input >> input_result; ++j){
                if ( i%5 == 0 ){
                    string::size_type size;
                    if(j == 0){
                        mm = stoi(input_result, &size);
                        continue;
                    }
                    if (j == 1 ){
                        string douhao;
                        douhao = input_result;
                        continue;
                    }  // 将索引输入到pair中  
                    if(j == 2 ){
                        nn = stoi(input_result, &size);
                        continue;
                    }
                    if(isReadMSE)
                    {
                        if(j == 3)
                        {
                            TrimmedMSEij = stod(input_result, &size);
                            continue;
                        }
                    } 
                } else {
                    int index = i%5 - 1;
                    string::size_type size1;
                    temp_motion(index,j) = stod(input_result,&size1);
                }
                } 
                if ( (i%5 -1) == 3 ){
                    // cout << temp_motion << endl;
                    Relativemotion tempMotion(mm,nn,temp_motion);
                    cout << mm << " , " << nn << endl;
                    if(isReadMSE)
                    {
                        tempMotion.TrimmedMSE = TrimmedMSEij;
                    } else
                    {
                        tempMotion.TrimmedMSE = 1.0;
                    }
                    // RelRs.push_back(tempMotion);
                    RelMs.push_back(tempMotion);
                }
            }
            } else {
            cout << "no such file" << endl;
            }
            fin.close();
        }

        cout << "Total " << RelMs.size() << " Relative pair-wise" << endl;
        return true;
    }

    // typedef Read_RelMatrix Read_initialRelmotion;

    bool Read_GlobalMatrix(const string & InputGlobalfile, Matrix4x4Arr & GlobalMats, int GlobalNumber,  bool format)
    {
        ifstream finTruth(InputGlobalfile, ios::in);
        if (!finTruth.is_open()){
            cout << "真实矩阵读取失败" << endl;
            return false;
        }

        Eigen::MatrixXd GroundTruth;
        GroundTruth.resize(GlobalNumber*4,4);
        
        GroundTruth.setZero();
        string line_infoT,input_resultT;
        {
            if(finTruth){
                for(int i = 0; getline(finTruth,line_infoT);i++){
                    stringstream input(line_infoT);
                    for(int j = 0; input >> input_resultT;++j){
                    string::size_type size;
                    GroundTruth(i,j) = stod(input_resultT,&size);
                    }
                }
            }
        }
        // finTruth.close();

        GlobalMats.clear();
        GlobalMats.reserve(GlobalNumber);
        for(int i = 0; i < GlobalNumber; i++){
            Trans temp = GroundTruth.block(i*4,0,4,4);
            GlobalMats.push_back(temp);
        }

        if(format){
            GlobalMats.clear();
            GlobalMats.reserve(GlobalNumber);
            Trans Mviewer0_inv = GroundTruth.block(0,0,4,4).inverse();
            for(int i = 0; i < GlobalNumber; i++ ){
                Trans temp = GroundTruth.block(i*4,0,4,4) * Mviewer0_inv;
                GlobalMats.push_back(temp);
            }
        }

        // if(format){
        //     Trans Mviewer0_inv = GlobalMats[0].inverse();
        //     for(int i = 0; i < GlobalMats.size(); i++){
        //         GlobalMats[i] = GlobalMats[i] * Mviewer0_inv;
        //     }
        // }
        return true;

    }


}

namespace Matrix_Write
{
    bool Write_Relmatrix(const string & OutputTxt, Relativemotions & RelMs, bool iswirteMSE)
    {
        std::ofstream fout(OutputTxt, std::ios::out);
        if(!fout.is_open()){
            cerr << "输出失败" << endl;
            return false;
        } else{
            print_highlight("Saving RelMs to"); print_value("%s", OutputTxt.c_str());
            cout << endl;
        }
        for(auto iter : RelMs){
            if(iswirteMSE)
            {
                fout << iter.i << "  ,  " << iter.j << "  " << iter.TrimmedMSE << endl;
                fout << iter.Tij << endl;
            }
            else
            {
                fout << iter.i << "  ,  " << iter.j << endl;
                fout << iter.Tij << endl;
            }
            
        }

        fout.flush();
        
        return true;
    }
    /// @brief 保存全局矩阵
    /// @param OutputTxt 输入:全局矩阵保存文件地址 @param vec_Global 输入:存储全局矩阵的vector
    bool Write_GlobalMatrix(const string & OutputTxt, Matrix4x4Arr & vec_Global)
    {
        std::ofstream fout(OutputTxt, std::ios::out);
        if(!fout.is_open()){
            cerr << "输出失败" << endl;
            return false;
        } else {
            print_highlight("Saving vec_Global to"); print_value("%s", OutputTxt.c_str());
            cout << endl;
        }
        for(auto iter : vec_Global){
            fout << iter << endl;
        }

        fout.flush();

        return true;
    }   
}



