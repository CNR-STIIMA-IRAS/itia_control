
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
#include <ros/time.h>
#include <control_models/ss.h>

int main(int argc,char** argv){

  
    ros::init(argc,&*argv,"controller_models_test");
    ros::NodeHandle nh;
    ros::Rate loop_rate(20);

    int order=4;
    itia::model::Integrator integrator(order);
    
    Eigen::VectorXd x0(order);
    x0.setZero();
    
    integrator.setState(x0);
    Eigen::VectorXd u(1);
    u(0)=1;
    ROS_INFO(" ");
    for (double t=0;t<=0.2;)
    {
      t+=0.1;
      integrator.simulate(0.1,u);
      Eigen::VectorXd x=integrator.getFutureState();
      std::cout << "t= " << t << ", u = " << u << ", x = " << x.transpose() << std::endl;
      integrator.setState(x);
    }
    
    order=2;
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(order,order);
    Eigen::MatrixXd B=Eigen::MatrixXd::Zero(order,1);
    Eigen::MatrixXd C=Eigen::MatrixXd::Zero(1,order);
    Eigen::MatrixXd D=Eigen::MatrixXd::Zero(1,1);
    
    A=-Eigen::MatrixXd::Identity(order,order);
    B.setOnes();
    
    x0.resize(order);
    x0.setZero();
    itia::model::ContinuosStateSpace sys(A,B,C,D);
    sys.setState(x0);
    for (double t=0;t<=0.2;)
    {
      t+=0.1;
      sys.simulate(0.1,u);
      Eigen::VectorXd x=sys.getFutureState();
      std::cout << "t= " << t << ", u = " << u << ", sys.x = " << x.transpose() << std::endl;
      sys.setState(x);
    }
    
    
    
    return 0;
}