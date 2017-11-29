#ifndef __ITIA_CONTROL_DISCRETE_STATE_SPACE__
#define __ITIA_CONTROL_DISCRETE_STATE_SPACE__

#include <Eigen/Dense>
# include <exception>
# include <stdexcept>
# include <boost/shared_ptr.hpp>
# include <itia_rutils/itia_rutils.h>
# include <ros/ros.h>
namespace itia
{
namespace control
{

inline bool loadControlParam(  ros::NodeHandle& nh, 
                              std::string& par_namespace, 
                              Eigen::MatrixXd& A, 
                              Eigen::MatrixXd& B, 
                              Eigen::MatrixXd& Baw, 
                              Eigen::MatrixXd& C, 
                              Eigen::MatrixXd& D, 
                              Eigen::VectorXd& max_output, 
                              Eigen::VectorXd& min_output, 
                              Eigen::VectorXd& initial_state)
{
  std::string par_name = par_namespace+"/A";
  if (!itia::rutils::getParam(nh, par_name, A))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  
  par_name = par_namespace+"/B";
  if (!itia::rutils::getParam(nh, par_name, B))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  
  par_name = par_namespace+"/Baw";
  if (!itia::rutils::getParam(nh, par_name, Baw))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  
  par_name = par_namespace+"/C";
  if (!itia::rutils::getParam(nh, par_name, C))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  
  par_name = par_namespace+"/D";
  if (!itia::rutils::getParam(nh, par_name, D))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  
  par_name = par_namespace+"/max_output";
  Eigen::MatrixXd matrix;
  if (!itia::rutils::getParam(nh, par_name, matrix))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  max_output = matrix.col(0);
  
  par_name = par_namespace+"/min_output";
  if (!itia::rutils::getParam(nh, par_name, matrix))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }    
  min_output = matrix.col(0);
  
  
  par_name = par_namespace+"/initial_state";
  if (!itia::rutils::getParam(nh, par_name, matrix))
  {
    ROS_ERROR("[DiscreteStateSpaceController] cannot find '%s' parameter!",par_name.c_str());
    return false;
  }
  initial_state = matrix.col(0);
  return true;
}


class DiscreteStateSpace
{
private:
  unsigned int m_order;
  unsigned int m_nin;
  unsigned int m_nout;
  
  Eigen::VectorXd m_state;
  Eigen::VectorXd m_output;
  Eigen::VectorXd m_output_nosat;
  Eigen::VectorXd m_max_output;
  Eigen::VectorXd m_min_output;
  
  Eigen::MatrixXd m_A;
  Eigen::MatrixXd m_B;
  Eigen::MatrixXd m_Baw;                                    // for anti-windup
  Eigen::MatrixXd m_C;
  Eigen::MatrixXd m_D;
  Eigen::MatrixXd m_initialization_matrix;

public:
  DiscreteStateSpace( const Eigen::Ref<Eigen::MatrixXd> A, 
                      const Eigen::Ref<Eigen::MatrixXd> B, 
                      const Eigen::Ref<Eigen::MatrixXd> Baw, 
                      const Eigen::Ref<Eigen::MatrixXd> C, 
                      const Eigen::Ref<Eigen::MatrixXd> D, 
                      const Eigen::Ref<Eigen::VectorXd> max_output, 
                      const Eigen::Ref<Eigen::VectorXd> min_output,
                      const Eigen::Ref<Eigen::VectorXd> initial_state)
  {
    m_A = A;
    m_B = B;
    m_Baw = Baw;
    m_C = C;
    m_D = D;
    m_max_output = max_output;
    m_min_output = min_output;

    m_order = m_A.rows();
    if (m_A.cols() != m_order)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix A is not square");
    
    m_nin = m_B.cols();
    if (m_B.rows() != m_order)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix B has wrong rows dimension");
    
    
    m_nout = m_C.rows();
    if (m_C.cols() != m_order)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix C has wrong cols dimension");
    
    if (m_Baw.rows() != m_order)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix Baw has wrong rows dimension");

    if (m_Baw.cols() != m_nout)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix Baw has wrong cols dimension");
    
    if (m_D.rows() != m_nout)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix D has wrong rows dimension");

    if (m_D.cols() != m_nin)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix D has wrong cols dimension");
    
    if (m_max_output.rows() != m_nout)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix max_output has wrong rows dimension");
    
    if (m_min_output.rows() != m_nout)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix max_output has wrong rows dimension");
    
    if (initial_state.rows() == m_order)
      m_state = initial_state;
    else
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix initial_state has wrong rows dimension");
      
    m_output_nosat.resize(m_nout);
    m_output_nosat.setZero();
    m_output.resize(m_nout);
    m_output.setZero();
    
  };
  
	void setInitializationMatrix(const Eigen::Ref<Eigen::MatrixXd> initialization_matrix){m_initialization_matrix=initialization_matrix;};
  Eigen::MatrixXd getInitializationMatrix(){return m_initialization_matrix;};
  Eigen::MatrixXd getAMatrix(){return m_A;};
  Eigen::MatrixXd getBMatrix(){return m_B;};
  Eigen::MatrixXd getCMatrix(){return m_C;};
  Eigen::MatrixXd getDMatrix(){return m_D;};
  
  Eigen::VectorXd getOutput(){return m_output;};
	
	void setState(const Eigen::Ref<Eigen::VectorXd> state){m_state=state;};
	
  
  Eigen::VectorXd& update(const Eigen::Ref<Eigen::VectorXd> input)
  {
    if (input.rows() != m_nin)
      throw std::invalid_argument("[ DiscreteStateSpace ] Matrix input has wrong rows dimension");
   
    m_output_nosat = m_C*m_state + m_D*input;
    m_state = m_A*m_state+m_B*input;
    for (unsigned int iO = 0;iO<m_nout;iO++)
    {
      if (m_output_nosat(iO) > m_max_output(iO))
        m_output(iO) = m_max_output(iO);
      else if (m_output_nosat(iO) < m_min_output(iO))
        m_output(iO) = m_min_output(iO);
      else
        m_output(iO) = m_output_nosat(iO);
      
    }
    // ANTIWINDUP
    m_state += m_Baw*(m_output-m_output_nosat);
    return m_output;
  }
  
  
  
};


typedef boost::shared_ptr<DiscreteStateSpace> DiscreteStateSpacePtr;
inline DiscreteStateSpacePtr createDiscreteStateSpareFromParam(ros::NodeHandle& nh, 
                                         std::string& par_namespace, 
                                         const  Eigen::Ref<Eigen::VectorXd> initialization_vector)
{
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd Baw; 
  Eigen::MatrixXd C; 
  Eigen::MatrixXd D;
  Eigen::VectorXd max_output; 
  Eigen::VectorXd min_output;
  Eigen::VectorXd initial_state;
  Eigen::MatrixXd initialization_matrix;
  
  if (!itia::control::loadControlParam(nh, par_namespace, A, B, Baw, C, D, max_output, min_output, initial_state))
  {
    ROS_ERROR("No valid parameters in %s/%s",nh.getNamespace().c_str(),par_namespace.c_str());
    throw std::invalid_argument("no valid parameters");
  }
  
  bool flag=false;
  if (!itia::rutils::getParam(nh, par_namespace+"/initialization_matrix", initialization_matrix))
    ROS_DEBUG("cannot find 'position/initialization_matrix' parameter, using default initial state!");
  else
  {
    flag=true;
    initial_state=initialization_matrix*initialization_vector;
  }
  DiscreteStateSpacePtr ptr;
  ptr.reset(new  itia::control::DiscreteStateSpace(A, B, Baw, C, D, max_output.col(0), min_output.col(0), initial_state.col(0)));
  if (flag)
    ptr->setInitializationMatrix(initialization_matrix);
  return ptr;
}


}
}

# endif
