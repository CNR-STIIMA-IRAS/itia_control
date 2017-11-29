
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

#ifndef __cImpendence__
#define __cImpendence__

#include <iostream>
#include <fstream>
#include <kdl/chain.hpp>
#include <kdl/jntarrayacc.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <ros/ros.h>
#include <eigen_conversions/eigen_kdl.h>
#include <geometry_msgs/PoseStamped.h>

namespace itia {
namespace controller {

class c_impedance_control {

private:

    ros::NodeHandle nh;
    Eigen::Matrix<double,6,6> a;
    Eigen::Matrix<double,6,6>Kp;
    Eigen::Matrix<double,6,6>Kd;
    Eigen::Matrix<double,6,6>M;
    Eigen::Matrix<double,6,6>Minv;

    int cont1;
    int cont2;
    int cont3;
    
    double m_domega_min_psp_deg_s2;
    
    KDL::Vector phi_cd, phi_cd_old;
    
    KDL::Vector delta_w_psp, delta_iw_psp, pos_t_sp, pos_r_sp, pos_t_msr, pos_r_msr, deltapos_t, deltapos_r, vel_t_sp, vel_r_sp, vel_t_msr, deltavel_r, deltavel_t;
    KDL::Vector start_pos_t, start_pos_r, phi_cd_start, phi_cd_cmd, phi_ed, dphi_cd_cmd, dphi_ed, ddphi_cd_cmd, ddphi_ed, delta_r_psp, delta_w_psp_old;
    bool firstEntry_psp;
    int imp_mode_r, imp_mode_t;
    KDL::Wrench force_start;
    
    KDL::Rotation R_start;
    
    KDL::Vector xDev_psp;
    KDL::Vector xdDev_psp;
    KDL::Vector xddDev_psp;
    
    Eigen::Quaterniond q_psp;
    KDL::Vector w_psp;
    KDL::Vector dw_psp;
    
    KDL::Vector xDev;
    KDL::Vector xdDev;
    KDL::Vector xddDev;
    Eigen::Matrix<double,3,1> force_filt;

    Eigen::Quaterniond q;
    KDL::Vector w;
    KDL::Vector dw;
    Eigen::Matrix<double,3,1> torque_filt;

    double forceFiltTime;
    double min_mass;
    double maxDevLinVel;
    double maxDevAngVel;
    double maxDevLinPos;
    double maxDevAngPos;
    double forceDeadBand;
    double torqueDeadBand;
    bool saturation;
    bool firstEntry;

    KDL::FrameAcc frameCom;
    KDL::FrameAcc frameDev;
    ros::Time initTime;
    
    geometry_msgs::PoseStamped impPosesetpoint;
    
    ros::Publisher  impPoseSp_pub;

    void integration ( const double& timeInterval );
    
    void integration_psp ( const double& timeInterval );

public:
    enum AngleConvection { EULER_ZYX    = 0, //default
                           RPY          = 1
                         } angleConvection;

    c_impedance_control ( const std::string& name );

    KDL::FrameAcc get_control_action ( const KDL::FrameAcc& frameSp,
                                       const KDL::Wrench& force_sp,
                                       const KDL::Wrench& force,
                                       const double& timeInterval );
    
    KDL::FrameAcc get_control_action_position_sp ( const KDL::FrameAcc cartFrame_Sp,
                                                   const KDL::FrameVel cartFrame_msr,
                                                   const KDL::Wrench& force, //tool frame
                                                   const double& timeInterval,
                                                   const bool& setfirstEntry_psp,
                                                   const double& add_damp
                                                 ) ;

    void set_parameters ( const std::vector<double>& mass,
                          const std::vector<double>& stiff,
                          const std::vector<double>& damp );
};

};
}

#endif