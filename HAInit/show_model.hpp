/*
 * @Description: 显示模型
 * @Version: 1.0
 * @Author: hyg
 * @Date: 2022-04-05 16:50:20
 * @LastEditors: hyg
 * @LastEditTime: 2022-06-03 11:16:24
 */
#ifndef show_model_hpp
#define show_model_hpp

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <iostream>
#include <pcl/recognition/trimmed_icp.h>
#include <pcl/recognition/ransac_based/auxiliary.h>
#include <time.h>
#include <boost/thread/thread.hpp>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/filter.h>
#include <regex>

#include <pcl/io/ply_io.h>
#include <pcl/console/print.h>
#include <pcl/console/time.h>
#include <pcl/console/parse.h>
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "Reg_definition.hpp"

using namespace std;
using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

namespace ShowPointCloud
{
    struct PCD_model
    {
        std::string filename; // 存放的pcd文件地址
        Eigen::Matrix4d ClobalTrans; // 存放pcd的全局motion
        pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud; //  初始pcd中的点云转换到全局坐标下的点云
        PCD_model(): Cloud(new pcl::PointCloud<pcl::PointXYZ>){};
    };

    bool SaveTotalModel(const std::vector<PCD_model> & alldata, const string & outfilename ){
        pcl::PointCloud<pcl::PointXYZ>::Ptr AllModel(new pcl::PointCloud<pcl::PointXYZ>);
        
        for(int i = 0; i < alldata.size(); i++){
            *AllModel = *AllModel + *alldata[i].Cloud;
        }
        if(pcl::io::savePCDFile(outfilename, *AllModel) < 0){
            std::cerr << "保存失败" << std::endl;\
            return false;
        }
        return true;
        
    }



    void visualAllCloud(const std::vector<PCD_model> & alldata){
        boost::shared_ptr<pcl::visualization::PCLVisualizer>  viewer(new pcl::visualization::PCLVisualizer("1"));
        
        const int MAX = 240, MIN = 30;

        std::vector<pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>> src_model;
        for(int i = 0; i < alldata.size(); i ++){
            unsigned long R = rand() % (MAX - MIN + 1) + MIN;
            unsigned long G = rand() % (MAX - MIN + 1) + MIN;
            unsigned long B = rand() % (MAX - MIN + 1) + MIN;
            pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h1(alldata[i].Cloud, R, G, B );
            src_model.push_back(src_h1);
        }

        viewer->setBackgroundColor(255,255,255);
        for (int i = 0; i < alldata.size(); i++){
            std::string s = std::to_string(i);
            viewer->addPointCloud(alldata[i].Cloud, src_model[i], s);
            
        }

        pcl::console::TicToc time;
        time.tic();
        while (!viewer->wasStopped())
        {
            viewer->spinOnce(100);
            boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            // if(time.toc()/1000 > 5){
            //     viewer.close();  
            //     break;  
            // }
        }
    }

    void visualAllCloud(const std::vector<pcl::PointCloud<pcl::PointXYZ>> & alldata){
        boost::shared_ptr<pcl::visualization::PCLVisualizer>  viewer(new pcl::visualization::PCLVisualizer("1"));
        
        const int MAX = 255, MIN = 0;

        std::vector<pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>> src_model;
        for(int i = 0; i < alldata.size(); i ++){
            unsigned long R = rand() % (MAX - MIN + 1) + MIN;
            unsigned long G = rand() % (MAX - MIN + 1) + MIN;
            unsigned long B = rand() % (MAX - MIN + 1) + MIN;
            pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h1(alldata[i].makeShared(), R, G, B );
            src_model.push_back(src_h1);
        }

        viewer->setBackgroundColor(255,255,255);
        for (int i = 0; i < alldata.size(); i++){
            std::string s = std::to_string(i);
            viewer->addPointCloud(alldata[i].makeShared(), src_model[i], s);
            
        }

        pcl::console::TicToc time;
        time.tic();
        while (!viewer->wasStopped())
        {
            viewer->spinOnce(100);
            boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            // if(time.toc()/1000 > 5){
            //     viewer.close();  
            //     break;  
            // }
        }
    }
    
    
    void show_CompareModel_intwoWindows(const pcl::PointCloud<pcl::PointXYZ>::Ptr &first_cloud, 
                            const pcl::PointCloud<pcl::PointXYZ>::Ptr &second_cloud)
    {
        boost::shared_ptr<pcl::visualization::PCLVisualizer>  Mviewer(new pcl::visualization::PCLVisualizer("model1"));

        int v1 = 0, v2 = 1;
        Mviewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
        Mviewer->setBackgroundColor(0, 0, 0, v1);
        Mviewer->addText("first cloud", 10, 10, "v1_text", v1);
        
        
        Mviewer->createViewPort(0.5, 0.0, 1, 1.0,v2);
        Mviewer->setBackgroundColor(0.2, 0.2, 0.2, v2);
        Mviewer->addText("second cloud", 10, 10, "v2_text", v2);

        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h1(first_cloud, 0 , 255 ,0 );
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h2(second_cloud, 90 , 255 ,15 );
        Mviewer->addPointCloud(first_cloud, src_h1, "cloud1", v1);
        Mviewer->addPointCloud(second_cloud, src_h2, "cloud2", v2);

        pcl::console::TicToc time;
        time.tic();
        while (!Mviewer->wasStopped())
        {
            Mviewer->spinOnce(100);
            boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            // if(time.toc()/1000 > 5){
            //     viewer.close();  
            //     break;  
            // }
        }

    }

    void show_SingleModel(const pcl::PointCloud<pcl::PointXYZ>::Ptr &SrcPointCloud)
    {
        boost::shared_ptr<pcl::visualization::PCLVisualizer>  Mviewer(new pcl::visualization::PCLVisualizer("1"));
        int v1 = 0;
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h1(SrcPointCloud, 0 , 255 ,0 );
        Mviewer->addPointCloud<pcl::PointXYZ>(SrcPointCloud, src_h1, "Cloud1", v1);
        Mviewer->setBackgroundColor( 0, 0, 0);
        pcl::console::TicToc time;
        time.tic();
        while (!Mviewer->wasStopped())
        {
            Mviewer->spinOnce(100);
            boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            // if(time.toc()/1000 > 5){
            //     viewer.close();  
            //     break;  
            // }
        }
    }


    void show_Comparemodels_inOneWindow(const pcl::PointCloud<pcl::PointXYZ>::Ptr &first_cloud, 
                            const pcl::PointCloud<pcl::PointXYZ>::Ptr &second_cloud)
    {
        boost::shared_ptr<pcl::visualization::PCLVisualizer>  viewer(new pcl::visualization::PCLVisualizer("model"));
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>  src_h1(first_cloud, 0, 255, 0);
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>  src_h2(second_cloud, 0, 0, 255);
    
        viewer->setBackgroundColor(255, 255, 255);
        viewer->addPointCloud(first_cloud, src_h1, "cloud1", 0);
        viewer->addPointCloud(second_cloud, src_h2,"cloud2", 0);
        // viewer.addPointCloud(first_cloud, src_h1, 0);
        // viewer.addPointCloud(second_cloud, src_h2, 0);

        pcl::console::TicToc time;
        time.tic();
        while (!viewer->wasStopped())
        {
            viewer->spinOnce(100);
            boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            // if(time.toc()/1000 > 5){
            //     viewer.close();  
            //     break;  
            // }
        }
        
    }



    int Target2SourceShow(const string &sourcePCDfile, const string &targetPCDfile, const Trans &TransMat){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_source(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_target(new pcl::PointCloud<pcl::PointXYZ>);

        if(pcl::io::loadPCDFile(sourcePCDfile, *cloud_source) < 0){
            std::cerr << "第一个点云路径读取失败" << endl;
            return -1;
        } else 
        {
            pcl::console::print_info("loading cloud_source points ");
            cout << cloud_source->points.size() << endl;
        }

        if(pcl::io::loadPCDFile(targetPCDfile, *cloud_target) < 0){
            std::cerr << "第二个点云路径读取失败" << endl;
            return -1;
        } else
        {
            pcl::console::print_info("loading cloud_target points ");
            cout << cloud_target->points.size() << endl;
        }
        {
            pcl::console::print_highlight("Processing:");
            pcl::console::print_value("%s", sourcePCDfile.c_str());
            cout << " and ";
            pcl::console::print_value("%s", targetPCDfile.c_str());
            cout << endl;
        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_end(new pcl::PointCloud<pcl::PointXYZ>); 
        pcl::transformPointCloud(*cloud_source, *cloud_end, TransMat);
        show_Comparemodels_inOneWindow(cloud_target, cloud_end);

        return 1;
        
    };

    int Target2SourceShow(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_source, const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_target, const Trans &TransMat){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_end(new pcl::PointCloud<pcl::PointXYZ>); 
        pcl::transformPointCloud(*cloud_source, *cloud_end, TransMat);
        show_Comparemodels_inOneWindow(cloud_target, cloud_end);

        return 1;
        
    };

    bool Get_PCDmodel(const Matrix4x4Arr & GlobalMatrix, 
                                    const string & sPCDDir,
                                std::vector<PCD_model> & PCDModels, bool isgetOriginshape = false)
    {
        
        PCDModels.clear();
        PCDModels.reserve(GlobalMatrix.size());

        std::vector<std::string> Vec_PCD;
        Vec_PCD.reserve(GlobalMatrix.size());

        if( !stlplus::folder_exists(sPCDDir))
        {
            std::cerr << "\n the input directory doesn't exist" << std::endl;
            return false;
        }
        std::vector<std::string> vec_pcd = stlplus::folder_files(sPCDDir);
        std::sort(vec_pcd.begin(),vec_pcd.end());
        for (auto iter = vec_pcd.begin(); iter != vec_pcd.end(); iter++){
            const std::string sPCDFileName = stlplus::create_filespec(sPCDDir,*iter);
            if(! stlplus::extension_part(sPCDFileName).compare("pcd") == 0)
            {
                continue;
            }
            Vec_PCD.push_back(sPCDFileName);
            std::cout << sPCDFileName << std::endl;
        }
        
        for(int i = 0; i < Vec_PCD.size(); i++){
            PCD_model temp_pcd;
            temp_pcd.filename = Vec_PCD[i];
            temp_pcd.ClobalTrans = GlobalMatrix[i];
            if(isgetOriginshape)
            {
                temp_pcd.ClobalTrans = Trans::Identity();
            }
            cout << temp_pcd.filename << endl;
        
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Trans(new pcl::PointCloud<pcl::PointXYZ>);

            if(pcl::io::loadPCDFile(temp_pcd.filename, *cloud_Trans) < 0){
                continue;
            };
            pcl::transformPointCloud(*cloud_Trans, *temp_pcd.Cloud, temp_pcd.ClobalTrans);
            PCDModels.push_back(temp_pcd);
            
        }
        return true;
    }


    bool Get_PCDmodel(const Matrix4x4Arr & GlobalMatrix, 
                                    const vector<string> & PCDfiles,
                                std::vector<PCD_model> & PCDModels)
    {
        
        PCDModels.clear();
        PCDModels.reserve(GlobalMatrix.size());

        std::vector<std::string> Vec_PCD;
        Vec_PCD.assign(PCDfiles.begin(), PCDfiles.end());
        for(int i = 0; i < Vec_PCD.size(); i++){
            PCD_model temp_pcd;
            temp_pcd.filename = Vec_PCD[i];
            temp_pcd.ClobalTrans = GlobalMatrix[i];
            cout << temp_pcd.filename << endl;
        
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Trans(new pcl::PointCloud<pcl::PointXYZ>);

            if(pcl::io::loadPCDFile(temp_pcd.filename, *cloud_Trans) < 0){
                continue;
            };
            pcl::transformPointCloud(*cloud_Trans, *temp_pcd.Cloud, temp_pcd.ClobalTrans); // 
            PCDModels.push_back(temp_pcd);
            
        }
        cout << endl;
        return true;
    }
}



#endif