
 // -------------------------------------------------------------------------------- 
 // Copyright (c) 2017 CNR-ITIA <iras@itia.cnr.it>
 // All rights reserved.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are met:
 //
 // 1. Redistributions of source code must retain the above copyright notice,
 // this list of conditions and the following disclaimer.
 // 2. Redistributions in binary form must reproduce the above copyright
 // notice, this list of conditions and the following disclaimer in the
 // documentation and/or other materials provided with the distribution.
 // 3. Neither the name of mosquitto nor the names of its
 // contributors may be used to endorse or promote products derived from
 // this software without specific prior written permission.
 //
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 // POSSIBILITY OF SUCH DAMAGE.
 // -------------------------------------------------------------------------------- 

#include <ros/ros.h>

#include <kdl_parser/kdl_parser.hpp>
#include <sensor_msgs/JointState.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <eigen_conversions/eigen_msg.h>
#include <dynamic_reconfigure/server.h>


#include <itia_clik/itia_clik.h>
#include <itia_rutils/itia_rutils.h>
#include <itia_rutils/itia_dynpar.h>


int main(int argc, char **argv){
  ros::init(argc, argv, "itia_clik_test");
  ros::NodeHandle nh;
  
  ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME,ros::console::levels::Debug);
  
  ROS_DEBUG("TEST OF ITIA_CLICK");
  
  std::string urdf_name;
  std::vector<std::string> joint_names;
  joint_names.push_back("shoulder_pan_joint");
  joint_names.push_back("shoulder_lift_joint");
  joint_names.push_back("elbow_joint");
  joint_names.push_back("wrist_1_joint");
  joint_names.push_back("wrist_3_joint");
  joint_names.push_back("wrist_2_joint");
  
  
  std::string urdf_initial_link_name = "base";
  std::string urdf_final_link_name   = "flange";
  int linear_flag=false;
  ROS_DEBUG_STREAM("QUI1");
  
//  ENRICO PATCH 2017-01-20 path for new itia_rutils
  
//   if (!itia::rutils::loadParam(nh,"clik_test/linear",&linear_flag))
//     return -1;
//   
//   if (!itia::rutils::loadParam(nh,"clik_test/robot_name",&urdf_name))
//     return -1;
//   
//   ROS_DEBUG_STREAM("robot name: " << urdf_name);
//   if (!itia::rutils::loadParam(nh,"clik_test/base_link",&urdf_initial_link_name))
//     return -1;
//   
//   if (!itia::rutils::loadParam(nh,"clik_test/end_link",&urdf_final_link_name))
//     return -1;

  if (!nh.getParam("clik_test/linear",linear_flag))
  {
    ROS_ERROR("Impossible to find %s\n", std::string("clik_test/linear").c_str() );
    return -1;
  }
  
  if (!nh.getParam("clik_test/robot_name",urdf_name))
  {
    ROS_ERROR("Impossible to find %s\n", std::string("clik_test/robot_name").c_str() );
    return -1;
  }
  
  ROS_DEBUG_STREAM("robot name: " << urdf_name);
  if (!nh.getParam("clik_test/base_link",urdf_initial_link_name))
  {
    ROS_ERROR("Impossible to find %s\n", std::string("clik_test/base_link").c_str() );
    return -1;
  }
  
  if (!nh.getParam("clik_test/end_link",urdf_final_link_name))
  {
    ROS_ERROR("Impossible to find %s\n", std::string("clik_test/end_link").c_str() );
    return -1;
  }
  
//  ENRICO PATCH 2017-01-20
  

  KDL::Tree tree;
  kdl_parser::treeFromParam(urdf_name,tree);
  
  
  int njoint=6;
  
  Eigen::VectorXd q0           = Eigen::VectorXd::Constant(njoint,1, 0.1);
  Eigen::VectorXd q            = Eigen::VectorXd::Constant(njoint,1, 0.1);
  Eigen::VectorXd meas_q       = Eigen::VectorXd::Constant(njoint,1, 0.1);
  Eigen::VectorXd dq           = Eigen::VectorXd::Constant(njoint,1, 0);
  Eigen::VectorXd twist_a_in_b = Eigen::VectorXd::Constant(6,1, 0);
  
  Eigen::VectorXd dq_min = Eigen::VectorXd::Constant(njoint,1,-M_PI);
  Eigen::VectorXd dq_max = Eigen::VectorXd::Constant(njoint,1, M_PI);
  Eigen::VectorXd q_min  = Eigen::VectorXd::Constant(njoint,1,-M_PI);
  Eigen::VectorXd q_max  = Eigen::VectorXd::Constant(njoint,1, M_PI);
  Eigen::VectorXd v_min  = Eigen::VectorXd::Constant(3,1,-0.1);
  Eigen::VectorXd v_max  = Eigen::VectorXd::Constant(3,1, 0.1);
  Eigen::VectorXd w_min  = Eigen::VectorXd::Constant(3,1,-M_PI);
  Eigen::VectorXd w_max  = Eigen::VectorXd::Constant(3,1, M_PI);
  
  ROS_DEBUG_STREAM("clik initialization");
  q0 << 0.25487217702390336, 0.12538685683415196, -1.535588355638758, 1.4102015051685415, 1.0402703393644015, -3.1415926491743593;
  
  itia::clik::Clik clik(joint_names);
  ROS_DEBUG("Clik constructed!");
  if (!clik.changeChain(urdf_name,urdf_initial_link_name,urdf_final_link_name))
    return -1;
  ROS_DEBUG("Clik chain load!");
//   clik.setLimits(dq_min,dq_max,q_min,q_max,v_min,v_max,w_min,w_max);
  clik.setInitialState(q0);
  ROS_DEBUG("Clik initial conditions set!");

  
  clik.setGains(50,0);
  
  
  Eigen::Affine3d T_bt(Eigen::Affine3d::Identity());
  Eigen::Vector3d x_t_in_b;
  x_t_in_b << 0.7,0.2,0.5;
  
  Eigen::Vector3d axis;
  axis << 0,0,1;
  
  T_bt.translate(x_t_in_b);
  T_bt.rotate(Eigen::AngleAxisd(M_PI*0.25,axis));
  Eigen::VectorXd v_t_in_b=Eigen::VectorXd::Constant(6,1, 0.0);
  
  ROS_DEBUG_STREAM("clik update");
  
  Eigen::Affine3d T_ba;
  double st=0.01;
  
  ros::Rate loop_rate(1.0/st);
  ros::Time t0=ros::Time::now();
  double t=0;
  
  ros::Publisher joint_pub = nh.advertise<sensor_msgs::JointState>("joint_states",1);
  ros::Publisher pose_pub  = nh.advertise<geometry_msgs::PoseStamped>("pose_state",1);
  ros::Publisher twist_pub = nh.advertise<geometry_msgs::TwistStamped>("twist_state",1);
  sensor_msgs::JointState     joint_msg;
  geometry_msgs::PoseStamped  pose_msg;
  geometry_msgs::TwistStamped twist_msg;
  
  
  itia::rutils::MsgReceiver<geometry_msgs::Pose> pose_rec;
  ros::Subscriber pose_sub=nh.subscribe<geometry_msgs::Pose>("target_pose",
                                                             1,
                                                             &itia::rutils::MsgReceiver<geometry_msgs::Pose>::callback,
                                                             &pose_rec);
  
  itia::rutils::MsgReceiver<sensor_msgs::JointState> joint_rec;
  ros::Subscriber joint_sub=nh.subscribe<sensor_msgs::JointState>(  "target_joint_state",
                                                                    1,
                                                                    &itia::rutils::MsgReceiver<sensor_msgs::JointState>::callback,
                                                                    &joint_rec);
  
  
  joint_msg.position.resize( njoint );
  joint_msg.velocity.resize( njoint );
  joint_msg.effort.resize(   njoint );
  joint_msg.name=joint_names;
  
  Eigen::VectorXd target_q(6);
  target_q << -1.0396842072150925, 0.014592895101750105, -2.1023522136154345, 2.091736009930085, 0.5276911685268104, -4.187455525141445;
  
  Eigen::VectorXd target_dq(6);
  target_dq << 0,0,0,0,0,0;
  
  ROS_INFO_STREAM("QUI");
  while (ros::ok())
  {
    if (linear_flag)
    {
      if (pose_rec.isANewDataAvailable())
        tf::poseMsgToEigen(pose_rec.getData(),T_bt);
      
      clik.update(st,T_bt,v_t_in_b,meas_q,itia::clik::MotionType::Linear,&T_ba,&twist_a_in_b,&q,&dq);
    }
    else 
    {
      if (joint_rec.isANewDataAvailable())
        for (int idx=0;idx<6;idx++)
          target_q(idx)=joint_rec.getData().position.at(idx);
      
      clik.update(st,target_q,target_dq,meas_q,itia::clik::MotionType::Linear,&T_ba,&twist_a_in_b,&q,&dq);
    }
    meas_q=q;
    
    t=(ros::Time::now()-t0).toSec();
    
    
    joint_msg.header.stamp=ros::Time::now();
    pose_msg.header.stamp=joint_msg.header.stamp;
    twist_msg.header.stamp=joint_msg.header.stamp;
    
    for (int idx=0;idx<njoint;idx++)
    {
      joint_msg.position.at(idx) =  q(idx);
      joint_msg.velocity.at(idx) = dq(idx);
    }
    
    tf::twistEigenToMsg(twist_a_in_b,twist_msg.twist);
    tf::poseEigenToMsg(T_ba,pose_msg.pose);
    
    joint_pub.publish(joint_msg);
    pose_pub.publish(pose_msg);
    twist_pub.publish(twist_msg);
    
    loop_rate.sleep();
    ros::spinOnce();
  }
  
  ROS_DEBUG_STREAM("end");
  return 0;  
}