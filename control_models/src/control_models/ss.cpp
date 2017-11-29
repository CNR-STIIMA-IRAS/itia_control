
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

#include <control_models/ss.h>
#include <boost/math/special_functions/factorials.hpp>

namespace itia{
namespace model{

ContinuosStateSpace::ContinuosStateSpace(const Eigen::MatrixXd A, const Eigen::MatrixXd B, const Eigen::MatrixXd C, const Eigen::MatrixXd D)
: m_A(A),
  m_B(B),
  m_C(C),
  m_D(D)
{
  m_nr_input=m_B.cols();
  m_nr_output=m_C.rows();
  m_nr_state=m_A.cols();
  
  if (m_A.rows()!=m_nr_state)
  {
    ROS_ERROR("matrix A has to be square");
    throw "matrix A has to be square";
  }
  if (m_B.rows()!=m_nr_state)
  {
    ROS_ERROR("B has not the same number of rows of A");
    throw "B has not the same number of rows of A";
  }
  if (m_C.cols()!=m_nr_state)
  {
    ROS_ERROR("C has not the same number of columns of A");
    throw "C has not the same number of columns of A";
  }
  if (m_D.cols()!=m_nr_input)
  {
    ROS_ERROR("D has not the same number of columns of B");
    throw "D has not the same number of columns of B";
  }
  if (m_D.cols()!=m_nr_output)
  {
    ROS_ERROR("D has not the same number of rows of C");
    throw "D has not the same number of rows of C";
  }
  m_x.resize(m_nr_state,1);
  m_x.setZero();
  m_x_new.resize(m_nr_state,1);
  m_x_new.setZero();
  
  if (std::abs(A.determinant())<1e-4)
  {
    //ROS_ERROR("zero eigen values not supperted yet");
    //throw "zero eigen values not supperted yet";
  } 
  else 
    m_Ainv=A.inverse();
  m_I.resize(m_nr_state,m_nr_state);
  m_I.setIdentity();
}

void ContinuosStateSpace::setState(const Eigen::VectorXd& x)
{
  if (x.rows()!=m_nr_state)
  {
    ROS_ERROR("x has not the same number of rows of A");
    throw "x has not the same number of rows of A";
  }
  m_x=x;
}

Eigen::VectorXd ContinuosStateSpace::simulate(const double& t, const Eigen::VectorXd& u)
{
  m_x_new=this->getFreeResponseState(t)*m_x+this->getForcedResponseState(t)*u;
  return m_C*m_x_new+m_D*u;
}

Eigen::VectorXd ContinuosStateSpace::simulate(const double& t, const Eigen::VectorXd& u, const Eigen::VectorXd& x)
{
  setState(x);
  simulate(t,u);
}

Eigen::VectorXd ContinuosStateSpace::simulate(const double& t, const Eigen::VectorXd& u, const Eigen::VectorXd& x_in, Eigen::VectorXd* x_out)
{
  setState(x_in);
  simulate(t,u);
  *x_out=m_x_new;
}

Eigen::VectorXd ContinuosStateSpace::getState() const
{
  return m_x;
}

Eigen::VectorXd ContinuosStateSpace::getFutureState() const
{
  return m_x_new;
}




Eigen::MatrixXd ContinuosStateSpace::getFreeResponseOutput(const double& t) const
{
  return m_C*getFreeResponseState(t);
}

Eigen::MatrixXd ContinuosStateSpace::getForcedResponseOutput(const double& t) const
{
  return m_C*getForcedResponseState(t)+m_D;
}

Eigen::MatrixXd ContinuosStateSpace::getFreeResponseState(const double& t) const
{
  ROS_ERROR("eigen MatrixExponential does not work!!!!!");
  Eigen::MatrixXd expA(m_nr_state,m_nr_state);
  expA.setZero();
  return expA;
}

Eigen::MatrixXd ContinuosStateSpace::getForcedResponseState(const double& t) const
{
  ROS_ERROR("eigen MatrixExponential does not work!!!!!");
  Eigen::VectorXd IexpAB(m_nr_state);
  IexpAB.setZero();
  return IexpAB;
}



Integrator::Integrator(const int& order): 
ContinuosStateSpace(Eigen::MatrixXd::Zero(order,order),
		    Eigen::MatrixXd::Zero(order,1),
		    Eigen::MatrixXd::Zero(1,order),
		    Eigen::MatrixXd::Zero(1,1))
{
  m_B(0,0)=1.0;
  m_C(0,m_nr_state-1)=1.0;
  m_A.block(1,0,order-1,order-1)=Eigen::MatrixXd::Identity(order-1,order-1);
}

Eigen::VectorXd Integrator::getFreeResponseOutput(const double& t) const
{
  Eigen::VectorXd CexpA=Eigen::VectorXd::Zero(m_nr_state);

  CexpA(0)=1;
  for (int idx_r=1;idx_r<m_nr_state;idx_r++)
  {
    CexpA(idx_r)=CexpA(idx_r-1)*t/(idx_r);
  }
  return CexpA;
}

Eigen::MatrixXd Integrator::getFreeResponseState(const double& t) const
{
  Eigen::MatrixXd expA=Eigen::MatrixXd::Zero(m_nr_state,m_nr_state);

  expA(0,0)=1;
  for (int idx_r=1;idx_r<m_nr_state;idx_r++)
  {
    expA(idx_r,0)=expA(idx_r-1,0)*t/(idx_r);
  }
  for (int idx_c=1;idx_c<m_nr_state;idx_c++)
    expA.block(idx_c,idx_c,m_nr_state-idx_c,1)=expA.block(0,0,m_nr_state-idx_c,1);
    
      
  
  return expA;
}

Eigen::MatrixXd Integrator::getForcedResponseState(const double& t) const
{
  Eigen::VectorXd IexpAB(m_nr_state);
  
  IexpAB(0)=t;
  for (int idx_r=1;idx_r<m_nr_state;idx_r++)
    IexpAB(idx_r)=IexpAB(idx_r-1)*t/(idx_r+1);
//   for (int idx_r=0;idx_r<m_nr_state;idx_r++)
//     IexpAB(idx_r)=pow(t,idx_r+1)/boost::math::factorial<double>(idx_r+1);
  
  return IexpAB;
}




}
}