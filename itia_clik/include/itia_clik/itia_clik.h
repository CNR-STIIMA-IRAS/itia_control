#ifndef __ITIA_MOTION_CLIK
#define __ITIA_MOTION_CLIK

#include <Eigen/Dense>
#include <ros/ros.h>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_kdl.h>
#include <kdl_parser/kdl_parser.hpp>
#include <urdf/model.h>



// #include <itia_clik/clikConfig.h>
#include <itia_rutils/itia_rutils.h>
# include <itia_mutils/frame_distance.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif 


# include <std_msgs/Float64MultiArray.h>

namespace itia {
namespace clik {
  
  enum MotionType { Joint=0, Linear=1, FiveAxisJoint = 2, FiveAxisLinear=3};  
  enum State {CONSTRUCTED, INITIALIZED, ERROR};
  
  class Clik
  {
  protected:
    // b = fix frame
    // a = actual tool frame
    // t = target tool frame
    
    unsigned int m_nr_jnt;
    KDL::ChainJntToJacSolver* m_jac_sol; //jacobian solver
    KDL::ChainFkSolverPos_recursive* m_fk_sol; //forward kinematics solver
    KDL::Tree m_tree;
    KDL::Chain m_chain;
    KDL::JntArray   m_jnt_array;                            //last observer  position
    KDL::Jacobian   m_jac;
    Eigen::VectorXd m_dq; //last joint velocity
    Eigen::Affine3d m_T_ba; //last joint position
    Eigen::VectorXd m_dq_min; //min joint velocity
    Eigen::VectorXd m_dq_max; //max joint velocity
    Eigen::VectorXd m_q_min; //min joint position
    Eigen::VectorXd m_q_max; //max joint position
    Eigen::VectorXd m_twist_min; //min linear twist
    Eigen::VectorXd m_twist_max;                            //max linear twist
    
    
    Eigen::VectorXd m_cartesian_error;
    Eigen::MatrixXd m_error_jacobian;
    
    double m_gain;                                          //controller gain
    double m_obs_gain;                                      //observer gain
    
    
    MotionType m_motion_type;
    State m_state;
    
    std::string m_robot_description_name;
    std::string m_base_link_name;
    std::string m_ee_link_name;
    std::vector<std::string> m_joint_names;
    std::vector<std::string> m_kdl_names;
    Eigen::MatrixXd m_kdl_to_clik_perm;
    Eigen::MatrixXd m_clik_to_kdl_perm;
    bool getJointProperties();
    ros::Publisher m_err_pub;
    ros::NodeHandle m_nh;
    std_msgs::Float64MultiArray m_msg;
    
  public:
    
    Clik(const std::vector<std::string>& joint_names);
    
    ~Clik();
    
    bool changeChain( const std::string& robot_description_name,
                      const std::string& base_link_name,
                      const std::string& ee_link_name);
    
    
    std::vector<std::string> getJointNames(){return m_joint_names;};
    
    void setGains(const double&  gain,const double& obs_gain){m_gain=gain;m_obs_gain=obs_gain;};
    
    void setInitialState(const Eigen::VectorXd& q0)
    {
      if (m_state != INITIALIZED)
      {
        ROS_ERROR("Try to set initial state to not initialed Clik object!");
        return;
      }
      m_jnt_array.data=m_clik_to_kdl_perm * q0;
    };
    
    void setLimits( const Eigen::VectorXd& dq_min,
                    const Eigen::VectorXd& dq_max,
                    const Eigen::VectorXd& q_min,
                    const Eigen::VectorXd& q_max,
                    const Eigen::VectorXd& v_min,
                    const Eigen::VectorXd& v_max,
                    const Eigen::VectorXd& w_min,
                    const Eigen::VectorXd& w_max);

    
    void update(  const double& delta_time, const Eigen::Affine3d& T_bt, const Eigen::VectorXd& V_t_in_b, const Eigen::VectorXd meas_joint, const itia::clik::MotionType& motion_type, Eigen::Affine3d* T_ba, Eigen::VectorXd* V_a_in_b, Eigen::VectorXd* q_obsv, Eigen::VectorXd* dq);
    void update(  const double& delta_time, Eigen::VectorXd target_q, Eigen::VectorXd target_dq, Eigen::VectorXd meas_joint, const itia::clik::MotionType& motion_type, Eigen::Affine3d* T_ba, Eigen::VectorXd* V_a_in_b, Eigen::VectorXd* q_obsv, Eigen::VectorXd* dq);
    
    bool ikine(const Eigen::Affine3d& T_bt, const Eigen::VectorXd q0, const MotionType& motion_type, Eigen::VectorXd* q);
    
  };
  
//   void callback(Clik* clik, itia_clik::clikConfig& config, uint32_t level);
  
}
}

#endif