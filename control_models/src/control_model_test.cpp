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