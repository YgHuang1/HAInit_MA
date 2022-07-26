/*
 * @Description: 点云相关处理
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-04-06 10:06:05
 * @LastEditors: hyg
 * @LastEditTime: 2022-07-26 15:36:23
 */

#ifndef PointCloud_Processing_hpp
#define PointCloud_Processing_hpp

#include "show_model.hpp"
#include <pcl/filters/extract_indices.h>
#include "Analysis_numeric.hpp"
#include "ReadWrite_RelMatrixfile.hpp"

using namespace pcl;

namespace PointCloud_processing
{   
    /// @brief 得到整体模型
    void GetScanSets(const string & sPCDDir, 
                    pcl::PointCloud<pcl::PointXYZ> &TotalCloud, 
                    const string & InitialMotionfile)
    {
        std::vector<std::string> Vec_PCD;
        {
        if( !stlplus::folder_exists(sPCDDir))
        {
            std::cerr << "\n the input directory doesn't exist" << std::endl;
        }
        std::vector<std::string> vec_pcd = stlplus::folder_files(sPCDDir);
        std::sort(vec_pcd.begin(),vec_pcd.end());
        for (auto iter = vec_pcd.begin(); iter != vec_pcd.end(); iter++){
            const std::string sPCDFileName = stlplus::create_filespec(sPCDDir,*iter);

            if(!stlplus::extension_part(sPCDFileName).compare("pcd") == 0){
                print_value("%s", sPCDFileName.c_str());
                cout << "unkown file format" << endl;
                continue;
            }
            Vec_PCD.push_back(sPCDFileName);
            cout << sPCDFileName << endl;
        }
        }
        int ScansNumber = Vec_PCD.size();
        Matrix4x4Arr GlobalMaotion;
        Matrix_Read::Read_GlobalMatrix(InitialMotionfile, GlobalMaotion, ScansNumber);
        vector<pcl::PointCloud<pcl::PointXYZ>> Allmodel;
        Allmodel.reserve(Vec_PCD.size());
        for(size_t i = 0; i < Vec_PCD.size(); i++)
        {
            pcl::PointCloud<pcl::PointXYZ> newcloud;
            pcl::io::loadPCDFile(Vec_PCD.at(i), newcloud);
            cout << "Loading " << Vec_PCD.at(i) << endl;
            cout << "Total " << newcloud.points.size() << " points" << endl;
            pcl::PointCloud<pcl::PointXYZ> cloud_Trans;
            pcl::transformPointCloud(newcloud, cloud_Trans, GlobalMaotion.at(i));
            Allmodel.push_back(cloud_Trans);
        }
        TotalCloud.points.clear();
        for(auto it : Allmodel)
        {
            TotalCloud = TotalCloud + it;
        }
    }
    
    void Cross_model(const std::shared_ptr<const pcl::PointCloud<pcl::PointXYZ>> &cloud,
                     std::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud_Cross, 
                     char axis_cross = 'y', double distance = 0.001f)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud_downsampled(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::VoxelGrid<pcl::PointXYZ> down_sampled;
        down_sampled.setInputCloud(cloud);
        down_sampled.setLeafSize(0.005f, 0.005f, 0.005f);
        down_sampled.filter(*Cloud_downsampled);

        std::vector<double> PointDis;
        PointDis.reserve(Cloud_downsampled->points.size());

        std::vector<int> indexs;
        double min, max, mean, median;
        int j = 0;

        switch (axis_cross)
        {
        case 'x':
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.x);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].x < (mean + distance/4) &&
                    cloud->points[i].x > (mean - distance/4) )
                {
                    indexs.push_back(j);
                }
                j++;
            }
            cout << mean - distance/4  << " " << mean + distance/4 << endl;
            break;

        case 'z':
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.z);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].z < (mean + distance/4) &&
                    cloud->points[i].z > (mean - distance/4) )
                {
                    indexs.push_back(j);
                }
                j++;
            }


            break;
        default:
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.y);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].y < (mean+ distance/4) &&
                    cloud->points[i].y > (mean - distance/4) )
                {
                    indexs.push_back(j);
                }
                j++;
            }
            
            break;
        }

        pcl::PointIndices index_ptr;
        index_ptr.indices = indexs;
        pcl::PointIndices::ConstPtr index_ptr2(new pcl::PointIndices(index_ptr));
        pcl::ExtractIndices<pcl::PointXYZ> extract;
        extract.setInputCloud(cloud);
        extract.setIndices(index_ptr2);
        extract.setNegative(false);
        extract.filter(*cloud_Cross);

    }


    void Cross_model(const std::shared_ptr<const pcl::PointCloud<pcl::PointXYZ>> &cloud,
                     std::shared_ptr<pcl::PointCloud<pcl::PointXYZ>> &cloud_Cross, 
                     double distanceMax , double distancemin , char axis_cross = 'y')
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud_downsampled(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::VoxelGrid<pcl::PointXYZ> down_sampled;
        down_sampled.setInputCloud(cloud);
        down_sampled.setLeafSize(0.005f, 0.005f, 0.005f);
        down_sampled.filter(*Cloud_downsampled);

        std::vector<double> PointDis;
        PointDis.reserve(Cloud_downsampled->points.size());

        std::vector<int> indexs;
        double min, max, mean, median;
        int j = 0;

        switch (axis_cross)
        {
        case 'x':
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.x);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].x < distanceMax &&
                    cloud->points[i].x > distancemin )
                {
                    indexs.push_back(j);
                }
                j++;
            }

            break;

        case 'z':
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.z);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].z < distanceMax &&
                    cloud->points[i].z > distancemin )
                {
                    indexs.push_back(j);
                }
                j++;
            }


            break;
        default:
            PointDis.clear();
            for(auto iter : Cloud_downsampled->points){
                PointDis.push_back(iter.y);
            }

            Analysis_numeric::minMaxMeanMedian(PointDis.begin(), PointDis.end(),
                                min, max, mean, median);
            cout << "min:" << min << " max:" << max << "mean:" << mean << "median:" << median << endl;
            
            j = 0;
            for(int i = 0; i <cloud->points.size(); i++){
                if(cloud->points[i].y < distanceMax &&
                    cloud->points[i].y > distancemin )
                {
                    indexs.push_back(j);
                }
                j++;
            }


            break;
        }

        pcl::PointIndices index_ptr;
        index_ptr.indices = indexs;
        pcl::PointIndices::ConstPtr index_ptr2(new pcl::PointIndices(index_ptr));
        pcl::ExtractIndices<pcl::PointXYZ> extract;
        extract.setInputCloud(cloud);
        extract.setIndices(index_ptr2);
        extract.setNegative(false);
        extract.filter(*cloud_Cross);

    }
    
    /**
     * @brief: 将PCD文件中的无效点去除后保存新的pcd文件
     * @param[in] PCDfilename 输入的pcd文件路径
     * @param[out] PCDsavedfilename 输出的pcd文件路径
     * @author: hyg
     */    
    bool removeNaNPoints(const string & PCDfilename, const string & PCDsavedfilename){ // 去除文件中pcd中包含的无效点
        pcl::PointCloud<pcl::PointXYZ> cloud;
        if(pcl::io::loadPCDFile(PCDfilename, cloud) <0)
        {
            std::cerr << "pcd文件读取失败" << endl;
            return false;
        }
        std::vector<int> indices;
        cout << "去除文件";
        pcl::console::print_value("%s", PCDfilename.c_str());
        cout << "中的无效点" << endl;
        cout << cloud.points.size() << " to ";
        pcl::removeNaNFromPointCloud(cloud, cloud, indices);
        cout << cloud.points.size() << endl;
        if(pcl::io::savePCDFile(PCDsavedfilename, cloud) < 0){
            std::cerr << "pcd文件保存失败" << endl;
            return false;
        }
        return true;
    }


    /**
     * @brief: 将整个目录中PCD文件中的无效点去除后保存新的pcd文件至新的目录
     * @param[in] origin_sPCDDir 输入的pcd文件目录
     * @param[out] New_sPCDDir 输出的pcd文件目录
     * @author: hyg
     */  
    bool removeNanPoints_Totalfile(const std::string &origin_sPCDDir, const std::string & New_sPCDDir){
          // 判断输入路径是否存在
        if( !stlplus::folder_exists(origin_sPCDDir))
        {
            std::cerr << "\n the input directory doesn't exist" << std::endl;
            return -1;
        }
        // 判断输出路径是否为空
        if (New_sPCDDir.empty()){
            std::cerr << "\n Invaild output directory" << std::endl;
            return -1;
        }
        // 判断输出文件夹是否可以生成
        if (!stlplus::folder_exists(New_sPCDDir))
        {
            if (! stlplus::folder_create(New_sPCDDir))
            {
                std::cerr<< "\n cannot create the output directory" << std::endl;
                return -1;
            }
        }

        // vec_pcd -> bun000.pcd
        std::vector<std::string> vec_pcd = stlplus::folder_files(origin_sPCDDir);
        std::sort(vec_pcd.begin(),vec_pcd.end());

        for(auto iter = vec_pcd.begin(); iter != vec_pcd.end(); iter++){
            const std::string sPCDfilename = stlplus::create_filespec(origin_sPCDDir,*iter);
            
            const std::string sPCDfilenamePart = stlplus::extension_part(sPCDfilename);
            if(!(sPCDfilenamePart.compare("pcd") == 0))
            {
                print_value("%s", sPCDfilenamePart.c_str());
                cout << "unkown file format" << endl;
                continue;
            }

            const std::string NewPCDfile = stlplus::create_filespec(New_sPCDDir, *iter);
            if(!removeNaNPoints(sPCDfilename, NewPCDfile)){
                continue;
            }

        }

        return true;

    }

    bool All2OnePointCloud(const std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &Allcloud,
        const std::vector<Trans> &GlobalMatrix,
        pcl::PointCloud<pcl::PointXYZ> &TotalPointCloud)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr TotalCloud(new pcl::PointCloud<pcl::PointXYZ>);
        for(int i = 0; i < Allcloud.size(); i++)
        {
            pcl::PointCloud<pcl::PointXYZ> cloud_Trans;
            pcl::transformPointCloud(*Allcloud[i], cloud_Trans, GlobalMatrix[i]);
            *TotalCloud = *TotalCloud + cloud_Trans;
        }
        
        TotalPointCloud = *TotalCloud;
        cout << Allcloud.size() << " Point Clouds in " << TotalPointCloud.points.size() << " Points " << endl;
        return (!Allcloud.empty()) && (!GlobalMatrix.empty()) && (GlobalMatrix.size() == Allcloud.size());
    }


}


#endif