
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