
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

#include <itia_clik/itia_clik.h>


namespace itia
{
namespace clik
{


Clik::Clik(const std::vector<std::string>& joint_names)
:
m_joint_names(joint_names)
{
  
  m_state=CONSTRUCTED;
  m_robot_description_name="";
  m_base_link_name="";
  m_ee_link_name="";
  m_jac_sol=NULL;
  m_fk_sol=NULL;
  
  m_nr_jnt= m_joint_names.size();
  m_jnt_array.resize(m_nr_jnt);
  m_jnt_array.data.setZero();
  m_jac.resize(m_nr_jnt);
  m_dq.resize(m_nr_jnt);
  m_dq.setZero();
  
  m_gain=0;
  m_dq.resize(m_nr_jnt);
  m_dq.setZero();
  
  m_twist_min.resize(6);
  m_twist_max.resize(6);
  m_twist_min.block(0,0,3,1)=Eigen::VectorXd::Constant(3,1,-1);
  m_twist_max.block(0,0,3,1)=Eigen::VectorXd::Constant(3,1, 1);
  m_twist_min.block(3,0,3,1)=Eigen::VectorXd::Constant(3,1,-0.5);
  m_twist_max.block(3,0,3,1)=Eigen::VectorXd::Constant(3,1, 0.5);
  m_q_max  = Eigen::VectorXd::Constant(6,1,  M_PI);
  m_q_min  = Eigen::VectorXd::Constant(6,1, -M_PI);
  m_dq_max = Eigen::VectorXd::Constant(6,1,  0.5);
  m_dq_min = Eigen::VectorXd::Constant(6,1, -0.5);
  
  
  m_err_pub = m_nh.advertise<std_msgs::Float64MultiArray>("error", 1);
  m_msg.data.resize(6);
 
  m_cartesian_error.resize(6);
  m_cartesian_error.setZero();
  m_error_jacobian.resize(6, 6);
  m_error_jacobian.setIdentity();
}

  


Clik::~Clik()
{
  if (m_jac_sol != NULL)
    delete m_jac_sol;
  if (m_fk_sol != NULL)
    delete m_fk_sol;  
}


void Clik::update(  const double& delta_time, 
                    const Eigen::Affine3d& T_bt, 
                    const Eigen::VectorXd& V_t_in_b, 
                    Eigen::VectorXd meas_joint,
                    const itia::clik::MotionType& motion_type, 
                    Eigen::Affine3d* T_ba, 
                    Eigen::VectorXd* V_a_in_b, 
                    Eigen::VectorXd* q_obsv, 
                    Eigen::VectorXd* dq
                 )
{
  
  ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME,ros::console::levels::Info);
  if (m_state != INITIALIZED)
  {
    ROS_ERROR("NOT INITIALIZED YET!");
    return;
  }
  
  if (motion_type==MotionType::Joint)
  {
    ROS_WARN("Joint interpolation is not allowed in cartesian move. Set equal to Linear interpolation");
    m_motion_type=MotionType::Linear;
  }
  if (motion_type==MotionType::FiveAxisJoint)
  {
    ROS_WARN("FiveAxisJoint interpolation is not allowed in cartesian move. Set equal to FiveAxisLinear interpolation");
    m_motion_type=MotionType::FiveAxisLinear;
  }
  if (motion_type==MotionType::FiveAxisLinear)
  {
    ROS_WARN("FiveAxisLinear interpolation is not implemented yet. Set equal to Linear interpolation");
    m_motion_type=MotionType::Linear;
  }
  
  // Sort names
  meas_joint = m_clik_to_kdl_perm * meas_joint;
  
  
  /* NOTE
   * A posteriori update of observer state
   * q_obsv(k|k)   a-posteriori estimate
   * q_obsv(k|k-1) a-priori estimate
   * q(k)          measure output
   * q_obsv(k|k)= inv(1+m_obs_gain*dt)*( q_obsb(k|k-1) + m_obs_gain*dt* q(k) )
   */
  m_jnt_array.data = ( m_jnt_array.data + m_obs_gain*delta_time*meas_joint ) / (1+m_obs_gain*delta_time);
  
  *q_obsv   = m_kdl_to_clik_perm * m_jnt_array.data;
  // compute position
  KDL::Frame frame_ba;
  m_fk_sol->JntToCart(m_jnt_array,frame_ba);
  tf::transformKDLToEigen(frame_ba,m_T_ba);

  // compute jacobian
  m_jac_sol->JntToJac(m_jnt_array,m_jac);

  
  // compute error and its jacobian
  itia::mutils::getFrameDistanceJac(m_T_ba, T_bt, &m_cartesian_error, &m_error_jacobian);
//   ROS_WARN_STREAM("m_T_ba = \n" << m_T_ba.matrix());
//   ROS_WARN_STREAM("T_bt = \n" << T_bt.matrix());
//   ROS_WARN_STREAM("error = " << m_cartesian_error.transpose());
  
  for (int idx=0;idx<m_cartesian_error.size();idx++)
  if (m_cartesian_error(idx)>0.005)
    m_cartesian_error(idx)=0.005;
  else if (m_cartesian_error(idx)<-0.005)
    m_cartesian_error(idx)=-0.005;
  
  Eigen::VectorXd velocity_in_b;
  Eigen::VectorXd next_error;
  
  /* NOTE
   * Controller update 
   * dq(k) = (J+(K*Je*J*st))\((v_t+K*e ) + K*Je*v_t*dt)
   */
  bool  flag = true;
  double ratio = 1;
  double ratio_cart = 1;
  double ratio_joint = 1;
  int ntrial = 5;
  int itrial = 0;
  
  Eigen::FullPivLU<Eigen::MatrixXd> jac_decomp(m_jac.data);
  jac_decomp.setThreshold(1e-2);

  while (flag && itrial<ntrial)
  {
    Eigen::MatrixXd den_control_matrix = Eigen::MatrixXd::Identity(6, 6) + m_error_jacobian*m_gain* delta_time/ratio;
    Eigen::FullPivLU<Eigen::MatrixXd> den_control_matrix_decomp(den_control_matrix);
    den_control_matrix_decomp.setThreshold(1e-2);
    Eigen::VectorXd num_control_vector = m_cartesian_error+ m_error_jacobian*V_t_in_b*delta_time/ratio;
    next_error = den_control_matrix_decomp.solve(num_control_vector);
    if (den_control_matrix_decomp.rank()<6)
    {
      Eigen::MatrixXd null=den_control_matrix_decomp.kernel();
      ROS_WARN_THROTTLE(0.1,"Singolarity point!");
      
      for (int iC=0;iC<null.cols();iC++)
      {
        Eigen::VectorXd null_versor=null.col(iC);
        null_versor.normalize();
        next_error=next_error-(next_error.dot(null_versor))*null_versor;
      }
    }
    for (int idx=0;idx<next_error.size();idx++)
      if (next_error(idx)>0.005)
        next_error(idx)=0.01;
      else if (next_error(idx)<-0.005)
        next_error(idx)=-0.005;
      
    velocity_in_b = (V_t_in_b + m_gain * next_error)/ratio;
    m_dq = jac_decomp.solve(velocity_in_b);
    if (jac_decomp.rank()<m_dq.rows())
    {
      Eigen::MatrixXd null=jac_decomp.kernel();
      ROS_WARN_THROTTLE(0.1,"Singolarity point!");
      
      for (int iC=0;iC<null.cols();iC++)
      {
        Eigen::VectorXd null_versor=null.col(iC);
        null_versor.normalize();
        m_dq=m_dq-(m_dq.dot(null_versor))*null_versor;
      }
    }
    
    velocity_in_b = m_jac.data*m_dq;
    
    ratio_cart = std::max((velocity_in_b.cwiseQuotient(m_twist_max)).maxCoeff(),
                    (velocity_in_b.cwiseQuotient(m_twist_min)).maxCoeff());
    ratio_joint = std::max((m_dq.cwiseQuotient(m_dq_max)).maxCoeff(),
                    (m_dq.cwiseQuotient(m_dq_min)).maxCoeff());
    
    double ratio_trial = std::max(ratio_cart, ratio_joint);
    if (ratio_trial >= (1+1e-3))
    {
      ratio *= 1.02*ratio_trial; 
      ROS_DEBUG_STREAM("qui, itrial = " <<  itrial << ", ratio = " << ratio << ",  ratio_trial = " << ratio_trial);
      itrial++;
    } 
    else
      flag = false;
  }
  
  bool saturated=false;
  double tol=1e-1;
  for (int idx=0;idx<m_jnt_array.data.rows();idx++)
  {
    if (m_jnt_array.data(idx)>(m_q_max(idx)-tol))
    {
      m_jnt_array.data(idx)=m_q_max(idx);
      ROS_DEBUG_STREAM_THROTTLE(15,"Joint(" << m_joint_names.at(idx) << ") limit reached! value = " << m_jnt_array.data(idx) << ", limit = " << m_q_max(idx));
      saturated=true;
    }
    else if (m_jnt_array.data(idx)<(m_q_min(idx)+tol))
    {
      m_jnt_array.data(idx)=m_q_min(idx);
      ROS_DEBUG_STREAM_THROTTLE(15,"Joint(" << m_joint_names.at(idx) << ") limit reached! value = " << m_jnt_array.data(idx) << ", limit = " << m_q_min(idx));
      saturated=true;
    }
  }
  if (saturated)
  {
    m_dq= Eigen::VectorXd::Constant(m_nr_jnt,1,0);
    velocity_in_b=Eigen::VectorXd::Constant(6,1,0);
  }
  
  /* NOTE
   * A priori update of observer state
   * q_obsv(k|k)   a-posteriori estimate
   * q_obsv(k+1|k) a-priori estimate
   * dq(k)         update from controller
   * q_obsv(k+1|k) = q_obsv(k|k) + dq(k)*st
   */
  m_jnt_array.data = m_jnt_array.data + m_dq*delta_time;

  for (int idx = 0;idx<6;idx++)
    m_msg.data.at(idx) = m_cartesian_error(idx);
  
  m_err_pub.publish(m_msg);

  *T_ba     = m_T_ba;
  *V_a_in_b = velocity_in_b;
  *dq       = m_kdl_to_clik_perm * m_dq;
  
  
};

void Clik::update(  const double& delta_time, Eigen::VectorXd target_q, Eigen::VectorXd target_dq, Eigen::VectorXd meas_joint, const itia::clik::MotionType& motion_type, Eigen::Affine3d* T_ba, Eigen::VectorXd* V_a_in_b, Eigen::VectorXd* q, Eigen::VectorXd* dq)
{
  if (m_state != INITIALIZED)
  {
    ROS_ERROR("NOT INITIALIZED YET!");
    return;
  }
  
  
  if (motion_type==MotionType::FiveAxisLinear)
  {
    ROS_DEBUG_THROTTLE(30,"FiveAxisLinear interpolation is not allowed in cartesian move. Set equal to FiveAxisJoint interpolation");
    m_motion_type=MotionType::FiveAxisJoint;
  }
  
  if (motion_type==MotionType::FiveAxisJoint)
  {
    ROS_DEBUG_THROTTLE(30,"FiveAxisJoint interpolation is not implemented yet. Set equal to Joint interpolation");
    m_motion_type=MotionType::Joint;
  }
  
  if (motion_type==MotionType::Linear)
  {
    ROS_DEBUG_THROTTLE(30,"Linear interpolation is not allowed in cartesian move. Set equal to Joint interpolation");
    m_motion_type=MotionType::Joint;
  }

  
  // Sort names
  target_q   = m_clik_to_kdl_perm * target_q;
  target_dq  = m_clik_to_kdl_perm * target_dq;
  meas_joint = m_clik_to_kdl_perm * meas_joint;
  
  /* NOTE
   * A posteriori update of observer state
   * q_obsv(k|k)   a-posteriori estimate
   * q_obsv(k|k-1) a-priori estimate
   * q(k)          measure output
   * q_obsv(k|k)= inv(1+m_obs_gain*dt)*( q_obsb(k|k-1) + m_obs_gain*dt* q(k) )
   */
  m_jnt_array.data = ( m_jnt_array.data + m_obs_gain*delta_time*meas_joint ) / (1+m_obs_gain*delta_time);
  *q        = m_kdl_to_clik_perm * m_jnt_array.data;
  
  // compute jacobian
  m_jac_sol->JntToJac(m_jnt_array,m_jac);
  
  /* NOTE
   * A priori update of observer state
   * q_obsv(k|k)   a-posteriori estimate
   * q_obsv(k+1|k) a-priori estimate
   * dq(k)         update from controller
   * q_obsv(k+1|k) = ( q_obsv(k|k) + (dqt(t) + K*q(k))*st)/(1+K*st)
   */
  KDL::JntArray q_priori;
  q_priori.data = ( m_jnt_array.data + (target_dq+ m_gain*target_q) *delta_time ) / (1+m_gain*delta_time);
  /* NOTE
   * note control law compute with k+1 estimate (backward intergration
   */
  Eigen::VectorXd error=target_q-q_priori.data;
  for (int idx=0;idx<m_nr_jnt;idx++)
    if (std::abs(error(idx))<3e-4)
      error(idx)=0;
    
  m_dq =  target_dq + m_gain * (target_q-q_priori.data);
  
  // compute jacobian with a priori information
  m_jac_sol->JntToJac(q_priori,m_jac);
  Eigen::VectorXd velocity_in_b = m_jac.data*m_dq;
  
  double ratio_jnt = std::max((m_dq.cwiseQuotient(m_dq_max)).maxCoeff(),
                          (m_dq.cwiseQuotient(m_dq_min)).maxCoeff());
  
  double ratio_cart = std::max((velocity_in_b.cwiseQuotient(m_twist_max)).maxCoeff(),
                   (velocity_in_b.cwiseQuotient(m_twist_min)).maxCoeff());
  
  double ratio = std::max(ratio_cart,ratio_jnt);
  if (ratio>1)
  {
    m_dq=m_dq/ratio;
    velocity_in_b/=ratio;
  }
  else
    ratio = 1;
  
  bool saturated=false;
  double tol=1e-1;
  for (int idx=0;idx<m_jnt_array.data.rows();idx++)
  {
    if (m_jnt_array.data(idx)>(m_q_max(idx)-tol))
    {
      m_jnt_array.data(idx)=m_q_max(idx);
      ROS_DEBUG_STREAM_THROTTLE(15,"Joint(" << m_joint_names.at(idx) << ") limit reached! value = " << m_jnt_array.data(idx) << ", limit = " << m_q_max(idx));
      saturated=true;
    }
    else if (m_jnt_array.data(idx)<(m_q_min(idx)+tol))
    {
      m_jnt_array.data(idx)=m_q_min(idx);
      ROS_DEBUG_STREAM_THROTTLE(15,"Joint(" << m_joint_names.at(idx) << ") limit reached! value = " << m_jnt_array.data(idx) << ", limit = " << m_q_min(idx));
      saturated=true;
    }
  }
  if (saturated)
  {
    velocity_in_b = Eigen::VectorXd::Constant(6,1,0);
    m_dq=Eigen::VectorXd::Constant(m_nr_jnt,1,0);
  }
  else
  {
    /* NOTE
   * A priori update of observer state with scaled velocity
   */
    m_jnt_array.data = ( m_jnt_array.data + (target_dq+ m_gain*target_q) *delta_time/ratio ) / (1+m_gain*delta_time/ratio);
  }
  
  KDL::Frame frame_ba;
  m_fk_sol->JntToCart(m_jnt_array,frame_ba);
  tf::transformKDLToEigen(frame_ba,m_T_ba);
  
  *V_a_in_b =velocity_in_b;
  *T_ba     = m_T_ba;
  *dq       = m_kdl_to_clik_perm * m_dq;
  
}


bool Clik::changeChain(const std::string& robot_description_name, const std::string& base_link_name, const std::string& ee_link_name)
{
  bool tree_changed=false;
  
  if (robot_description_name.compare(m_robot_description_name))
  {
    tree_changed=true;
    m_robot_description_name=robot_description_name;
    if (!kdl_parser::treeFromParam(m_robot_description_name,m_tree))
    {
      ROS_ERROR("Error on extracting tree from param");
      ROS_ERROR_STREAM("robot = " << m_robot_description_name);
      m_state = ERROR;
      return false;
    }
  }
  
  bool chain_changed=false;
  if (ee_link_name.compare(m_ee_link_name))
  {
    chain_changed=true;
    m_ee_link_name=ee_link_name;
    ROS_DEBUG_STREAM("Tool changing. robot = " << m_robot_description_name << ", base link = " << m_base_link_name << ", end effector link = " << ee_link_name);
  }
  
  if (base_link_name.compare(m_base_link_name))
  {
    chain_changed=true;
    m_base_link_name=base_link_name;
    ROS_DEBUG_STREAM("Base changing. robot = " << m_robot_description_name << ", base link = " << m_base_link_name << ", end effector link = " << ee_link_name);
  }
  
  if (chain_changed || tree_changed )
  {
    if (!m_tree.getChain(m_base_link_name,ee_link_name,m_chain))
    {
      ROS_ERROR("Error on extracting chain from tree");
      ROS_ERROR_STREAM("robot = " << m_robot_description_name << ", base link = " << m_base_link_name << ", end effector link = " << m_ee_link_name);
      m_state = ERROR;
      return false;
    }
    
    if (!getJointProperties())
    {
      m_state=ERROR;
      return false;
    }
    m_dq_max = Eigen::VectorXd::Constant(m_nr_jnt,1,  1);
    m_dq_min = Eigen::VectorXd::Constant(m_nr_jnt,1, -1);

    if (m_jac_sol != NULL)
      delete m_jac_sol;
    if (m_fk_sol != NULL)
      delete m_fk_sol;
    
    m_jac_sol = new KDL::ChainJntToJacSolver(m_chain);
    m_fk_sol  = new KDL::ChainFkSolverPos_recursive(m_chain);
    
    KDL::Frame frame_ba;
    m_fk_sol->JntToCart(m_jnt_array,frame_ba);
    tf::transformKDLToEigen(frame_ba,m_T_ba);
  }
  
  m_state = INITIALIZED;
  return true;
  
}

bool Clik::getJointProperties()
{
  urdf::Model model;
  model.initParam(m_robot_description_name);
  
  int nJnt=m_chain.getNrOfJoints();
  if (nJnt!=m_nr_jnt)
  {
    ROS_ERROR("desired KDL tree has a different number of joints!");
    return false;
  }
  m_kdl_names.resize(nJnt);
  m_q_max.resize(nJnt);
  m_q_min.resize(nJnt);
  m_dq_max.resize(nJnt);
  m_dq_min.resize(nJnt);
  
  int iJnt=0;
  for (int idx=0;idx<m_chain.getNrOfSegments();idx++)
  {
    ROS_DEBUG_STREAM("m_chain.getSegment(idx).getJoint().getName() = " << m_chain.getSegment(idx).getJoint().getName());
    ROS_DEBUG_STREAM("m_chain.getSegment(idx).getJoint().getType() = " << m_chain.getSegment(idx).getJoint().getTypeName());
      
    if (m_chain.getSegment(idx).getJoint().getType()!=KDL::Joint::JointType::None)
    {
      auto joint=model.getJoint(m_chain.getSegment(idx).getJoint().getName());
      m_kdl_names.at(iJnt) = m_chain.getSegment(idx).getJoint().getName();
      ROS_DEBUG_STREAM("m_joint_names.at(iJnt) = " << m_joint_names.at(iJnt));
      m_q_max(iJnt)  = joint->limits->upper;
      m_q_min(iJnt)  = joint->limits->lower;
      m_dq_max(iJnt) = joint->limits->velocity;
      m_dq_min(iJnt) = -joint->limits->velocity;
      iJnt++;
    }
  }
  
  ROS_DEBUG_STREAM("chain parsed!");
  
  if (itia::rutils::permutationName(m_joint_names,m_kdl_names,&m_kdl_to_clik_perm,&m_clik_to_kdl_perm))
  {
    ROS_DEBUG_STREAM("joints loaded!");
    m_state = INITIALIZED;
    return true;
  } 
  else
  {
    ROS_ERROR_STREAM("Joint names not found in KDL tree");
    return false;
  }
  
  return true;
}

void Clik::setLimits(const Eigen::VectorXd& dq_min, const Eigen::VectorXd& dq_max, const Eigen::VectorXd& q_min, const Eigen::VectorXd& q_max, const Eigen::VectorXd& v_min, const Eigen::VectorXd& v_max, const Eigen::VectorXd& w_min, const Eigen::VectorXd& w_max)
{
  m_twist_min.block(0,0,3,1)=v_min;
  m_twist_max.block(0,0,3,1)=v_max;
  m_twist_min.block(3,0,3,1)=w_min;
  m_twist_max.block(3,0,3,1)=w_max;
  m_q_max  = m_clik_to_kdl_perm * q_max;
  m_dq_max = m_clik_to_kdl_perm * dq_max;
  m_q_min  = m_clik_to_kdl_perm * q_min;
  m_dq_min = m_clik_to_kdl_perm * dq_min;
}


bool Clik::ikine(const Eigen::Affine3d& T_bt, const Eigen::VectorXd q0, const MotionType& motion_type, Eigen::VectorXd* q)
{
  if (m_state != INITIALIZED)
  {
    ROS_ERROR("NOT INITIALIZED YET!");
    return false;
  }
  
  Eigen::VectorXd q0_kdl_order = m_clik_to_kdl_perm * q0;
  Eigen::Affine3d T_ba;
  bool flag = true;
  double st = 1e-3;
  Eigen::VectorXd v_a_in_b;
  Eigen::VectorXd dq;
  
  int max_iter = 1e6;
  int iter = 0;
  double toll = 1e-3; 
  Eigen::VectorXd distance = Eigen::VectorXd::Constant(6, 1, 0);
  while (flag && iter<max_iter)
  {
    update(st, T_bt, Eigen::VectorXd::Constant(6, 1, 0.0), q0_kdl_order, motion_type,  &T_ba, &v_a_in_b, q, &dq);
    itia::mutils::getFrameDistance(T_bt, T_ba, &distance);
    iter++;
    if (distance.norm()>toll);
      return true;
  }
  
  return false;
};



};
};