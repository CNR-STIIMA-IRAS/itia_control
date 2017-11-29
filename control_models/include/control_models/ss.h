#ifndef __ITIA_STATE_SPACE__
#define __ITIA_STATE_SPACE__

#include <Eigen/Core>
#include <ros/console.h>
#include <unsupported/Eigen/MatrixFunctions>

namespace itia
{
namespace model
{
  class ContinuosStateSpace
  {
  protected:
    Eigen::MatrixXd m_A;
    Eigen::MatrixXd m_Ainv;
    Eigen::MatrixXd m_B;
    Eigen::MatrixXd m_C;
    Eigen::MatrixXd m_D;
    Eigen::MatrixXd m_I;
    Eigen::VectorXd m_x;
    Eigen::VectorXd m_x_new;
    
    unsigned int m_nr_state;
    unsigned int m_nr_input;
    unsigned int m_nr_output;
    
  public:
    ContinuosStateSpace(const Eigen::MatrixXd A,
                        const Eigen::MatrixXd B,
                        const Eigen::MatrixXd C,
                        const Eigen::MatrixXd D);
    
    
    void setState(const Eigen::VectorXd& x);
    
    Eigen::VectorXd getState() const; 
    Eigen::VectorXd getFutureState() const; 
    Eigen::VectorXd simulate(const double& t, const Eigen::VectorXd& u);
    Eigen::VectorXd simulate(const double& t, const Eigen::VectorXd& u, const Eigen::VectorXd& x_in);
    Eigen::VectorXd simulate(const double& t, const Eigen::VectorXd& u, const Eigen::VectorXd& x_in, Eigen::VectorXd* x_out);
    virtual Eigen::MatrixXd getFreeResponseState(const double& t) const;
    virtual Eigen::MatrixXd getForcedResponseState(const double& t) const;
    Eigen::MatrixXd getFreeResponseOutput(const double& t) const;
    Eigen::MatrixXd getForcedResponseOutput(const double& t) const;
    Eigen::MatrixXd getA()const {return m_A;};
    Eigen::MatrixXd getB()const {return m_B;};
    Eigen::MatrixXd getC()const {return m_C;};
    Eigen::MatrixXd getD()const {return m_D;};
    
    
  };
  
  class Integrator: public ContinuosStateSpace
  {
  public:
    Integrator(const int& order);
    Eigen::VectorXd getFreeResponseOutput(const double& t) const;
    Eigen::MatrixXd getFreeResponseState(const double& t) const;
    Eigen::MatrixXd getForcedResponseState(const double& t) const;
  };
  
}
}



#endif