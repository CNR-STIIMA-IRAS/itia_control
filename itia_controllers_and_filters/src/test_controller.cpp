#include <ros/ros.h>
#include <ros/time.h>
#include <itia_controllers_and_filters/pid.h>

#include <itia_controllers_and_filters/discrete_state_space.h>

int main(int argc,char** argv)
{
	// ------ Init ROS ------
	ros::init(argc,&*argv,"pid_controller");
	ros::NodeHandle nh;
	ros::Rate loop_rate(20);

//     itia::controller::pid _pid;
//     double sp_value=1;
//     double output=0;
//     while (ros::ok()){
//         ros::spinOnce();
//         loop_rate.sleep();
//         _pid.step(sp_value,output,0.0);
//         std::cout << _pid.control_action << std::endl;
//     }
		
	return 0;	
}