#ifndef __ITIA_TF__
#define __ITIA_TF__

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>

//--------------------------------------------------------------------------
//------------------ generic transfer function -----------------------------
//--------------------------------------------------------------------------

namespace itia{
    namespace filter{

class filter {
public:
  
  filter();
  
  void first_ord_LP_filt( const double time_constant
                          , const double gain
                          , const double sample_period
                          , const double initial_output
                          , const double initial_input  );
  
  void first_ord_HP_filt( const double time_constant
                          , const double gain
                          , const double sample_period
                          , const double initial_output
                          , const double initial_input  );
  
  void set_st(const double sample_period);
  
  void setCoeffs(const std::vector<double>& coeff_num, const std::vector<double>& coeff_den)
  {
    m_coeff_num = coeff_num;
    m_coeff_den = coeff_den;
    order = m_coeff_den.size()-1;
    m_out_state.resize(order);
    m_in_state.resize(order+1);
    std::fill(m_out_state.begin(), m_out_state.end(), 0.0);
    std::fill(m_in_state.begin(), m_in_state.end(), 0.0);
  }
  
  void setInitConditions(const std::vector<double>& out_state, const std::vector<double>& in_state)
  {
    m_out_state = out_state;
    m_in_state = in_state;
  }
  
protected:
  double st;
  int order;
  std::vector<double> m_out_state;    //order past output
  std::vector<double> m_in_state;         //order+1 past inputs
  std::vector<double> m_coeff_num;    //order+1 numerator coefficients
  std::vector<double> m_coeff_den;    //order denominator coefficients (of a monic polynomial)
  
};

class tf: public filter {
public: 
    
    // costructor
    tf();
        
    // data
    double tf_out;
  
    // method
    void step(const double tf_in);  
    
};


inline filter::filter(){
  set_st(1.0);
  order=0;
  m_coeff_num.reserve(order+1);
  m_coeff_num[0]=1.0;
};


inline void filter::set_st(const double sample_period){
  st=sample_period;
};


inline void filter::first_ord_LP_filt(const double time_constant, 
                                      const double gain, 
                                      const double sample_period, 
                                      const double initial_output=0.0, 
                                      const double initial_input=0.0){
  double a;
  
  set_st(sample_period);
  if (time_constant<=0)
    a=0;
  else
    a=exp(-st/time_constant);
  
  order = 1;
  // tf=(1-a)z^-1/(1-az^-1);
  m_coeff_den.resize(order);
  m_coeff_num.resize(order+1);
  m_out_state.resize(order);
  m_in_state.resize(order+1);
  
  m_coeff_num.at(0)   =   0;
  m_coeff_num.at(1)   =   (1-a)*gain;
  m_coeff_den.at(0)   =   -a;
  m_in_state.at(0)    =   initial_input;
  m_in_state.at(1)    =   initial_input;
  m_out_state.at(0)   =   initial_output;

};


inline void filter::first_ord_HP_filt(const double time_constant, 
                                      const double gain, 
                                      const double sample_period, 
                                      const double initial_output=0.0, 
                                      const double initial_input=0.0){
  double a;
  
  set_st(sample_period);
  if (time_constant<=0)
    a=0;
  else
    a=exp(-st/time_constant);
  
  order = 1;
  // tf=(1-a)z^-1/(1-az^-1);
  m_coeff_den.resize(order);
  m_coeff_num.resize(order+1);
  m_out_state.resize(order);
  m_in_state.resize(order+1);
  
  m_coeff_num.at(0)   =   (1-a)*gain/st;
  m_coeff_num.at(1)   =   -(1-a)*gain/st;
  m_coeff_den.at(0)   =   -a;
  m_in_state.at(0)    =   initial_input;
  m_in_state.at(1)    =   initial_input;
  m_out_state.at(0)   =   initial_output;
    
};



inline tf::tf(){

};


inline void tf::step(const double tf_in) {
  if (order==0){
    tf_out = m_coeff_num[0] * (tf_in);
  } else {
    int idx;
    
    tf_out = m_coeff_num.at(0) *tf_in;
    
    for (idx=1;idx <= order;idx++)
    {
//       printf("-m_coeff_den[idx] * m_out_state[idx-1]=%e\n",-m_coeff_den[idx] * m_out_state[idx-1]); 
      tf_out += m_coeff_num[idx] * m_in_state[idx-1] - m_coeff_den[idx] * m_out_state[idx-1];
    }
    
    std::rotate(m_out_state.rbegin(), m_out_state.rbegin() + 1, m_out_state.rend());
    std::rotate(m_in_state.rbegin(), m_in_state.rbegin() + 1, m_in_state.rend());
    m_out_state.at(0) = tf_out;  
    m_in_state.at(0)  = tf_in;  
    
//     printf("b(k)=%e, in(k)=%e\n",m_coeff_num.at(0), tf_in); 
//     for (idx = 0;idx<order;idx++)
//     {
//       printf("b(%d)=%e, in(k-%d)=%e, a(%d)=%e,out(k-%d)=%e, \n", idx,m_coeff_num.at(idx+1), idx, m_in_state.at(idx), idx,m_coeff_den.at(idx+1),idx, m_out_state.at(idx));
//     }
    
//     ROS_INFO("tf_out=%e", tf_out);
//     std::string str;
//     std::getline(std::cin, str);
  }
};


}
}
#endif
