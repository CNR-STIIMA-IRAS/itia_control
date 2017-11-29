#ifndef __PID__
#define __PID__

#include <itia_controllers_and_filters/tf.h>
#include <itia_futils/itia_futils.h>
#include <itia_msgs/send_float64.h>
#include <itia_msgs/send_string.h>
#include <ros/console.h>
#include <ros/time.h>
#include <ros/ros.h>

namespace itia
{
namespace controller
{

  // piddf
    class pid 
    {

    public: 

      struct pidPar 
      {
        // PID Parameters (defined for each PID loop)
        double Kp;
        double Ki;
        double Kd;
        double Kaw;
        double Td;
        double Tf;
        double st;
        
        // Control action limits 
        double Umax;
        double Umin;
        
        // Control action slope limits 
        double DUmax;
        double DUmin;
      };
        
      pid(ros::NodeHandle* ptr_nh) : m_ptr_nh(ptr_nh), m_setpoint(0.0) 
      {
        Istate              = 0.0;
        control_action      = 0.0;
        control_action_past = 0.0;
        
        m_internal_pid_param.st     = 1.0;
        m_internal_pid_param.Kp     = 1.0;
        m_internal_pid_param.Ki     = 0.0;
        m_internal_pid_param.Kd     = 0.0;
        m_internal_pid_param.Kaw    = 0.0;
        m_internal_pid_param.Td     = 0.0;
        m_internal_pid_param.Tf     = 0.0;

        m_internal_pid_param.Umax   =  1e2;
        m_internal_pid_param.Umin   = -1e2;
        m_internal_pid_param.DUmax  =  1e5;
        m_internal_pid_param.DUmin  = -1e5;
        
        last_time=ros::Time::now();
      }

      void update_param(const pidPar& in_pid_par, const double initial_control_action)
      {
          
        m_internal_pid_param.st   = in_pid_par.st;
        m_internal_pid_param.Kp   = in_pid_par.Kp;
        m_internal_pid_param.Ki   = in_pid_par.Ki;
        m_internal_pid_param.Kd   = in_pid_par.Kd;
        m_internal_pid_param.Kaw  = in_pid_par.Kaw;
        m_internal_pid_param.Td   = in_pid_par.Td;
        m_internal_pid_param.Tf   = in_pid_par.Tf;
    
        vel_error_filter.first_ord_LP_filt( m_internal_pid_param.Td, m_internal_pid_param.Kd, m_internal_pid_param.st, 0.0, 0.0 );
        control_action_filter.first_ord_LP_filt( m_internal_pid_param.Tf,1.0, m_internal_pid_param.st, 0.0, 0.0 );
        
        Istate = initial_control_action;
        control_action_past = initial_control_action;
          
      }

      void step(const double output, const double ff_action,const double tracking_signal = 0.0 )
      {
    
        double pos_err;
        double control_action_nofilt;
        double control_action_nosat;
        
        pos_err = m_setpoint - output;

        
        ros::Time act_time = ros::Time::now();
        ros::Duration delta_time = act_time - last_time;
                
        double time_interval = delta_time.toSec();
        if ( time_interval > m_internal_pid_param.st*10 )
           time_interval = m_internal_pid_param.st*10;  
        
        if ( m_internal_pid_param.Ki != 0)
        {
          double Ti = m_internal_pid_param.Kp / m_internal_pid_param.Ki;
          Istate = Istate + m_internal_pid_param.Kp * exp( -time_interval / Ti ) * pos_err; //PIDPLUS 
        }
        
        vel_error_filter.step(pos_err);
        last_time = act_time;
        
        control_action_nofilt = m_internal_pid_param.Kp * pos_err + Istate ;//+ vel_error_filter.tf_out;
              
        control_action_filter.step(control_action_nofilt);
        control_action_nosat = control_action_filter.tf_out + ff_action;
        
        double tmp_control_action;
        // slew-rate constraints
        if ( (control_action_nosat-control_action_past) > m_internal_pid_param.DUmax )
          tmp_control_action = control_action_past + m_internal_pid_param.DUmax;
        else if ( (control_action_nosat-control_action_past) < m_internal_pid_param.DUmin )
          tmp_control_action = control_action_past + m_internal_pid_param.DUmin;
        else
          tmp_control_action = control_action_nosat;
        
        saturate( m_internal_pid_param.Umax, m_internal_pid_param.Umin, &tmp_control_action );
        control_action = tmp_control_action;
        control_action_past = control_action;
               
        // anti-windup action
        if ( m_internal_pid_param.Kaw != 0 )
        {
          double Taw = m_internal_pid_param.Kp / m_internal_pid_param.Kaw;
          Istate = Istate + m_internal_pid_param.Kp * exp( -time_interval/Taw ) * ( tracking_signal + control_action - control_action_nosat ); //PIDPLUS 
        }
      }

      double get_control_action()
      { 
        return control_action;
      }
      
      double get_setpoint()
      {
        return m_setpoint;
      }
      
      void set_setpoint(const double setpoint)
      {
        m_setpoint = setpoint;
      }

      void saturate( const double& uBound, const double& lBound,double *sat_output )
      {
        // amplitude constraints
        if (*sat_output  > uBound )
            *sat_output  = uBound;
        else if ( *sat_output  < lBound )
            *sat_output  = lBound;
          
      }

      bool cb_receive_setpoint( itia_msgs::send_float64::Request& req,
                                itia_msgs::send_float64::Response& res)
      {
        ROS_INFO( " [ %s%s:%d%s ]\t receive new setpoint value: %f\n",GREEN, __FUNCFILE__, __LINE__, RESET,req.value);
        m_setpoint = req.value;
        return res.res = (true);
      }
      
      bool cb_update_param(  itia_msgs::send_string::Request&  req,
                             itia_msgs::send_string::Response& res )
      {
          
        std::string tmp_str=req.name.data+"/pid_param";
        if(!m_ptr_nh->getParam(tmp_str+"/st",m_internal_pid_param.st))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param st=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.st);

        if(!m_ptr_nh->getParam(tmp_str+"/Kp",m_internal_pid_param.Kp))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Kp=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Kp);
        
        
        if(!m_ptr_nh->getParam(tmp_str+"/Ki",m_internal_pid_param.Ki))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Ki=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Ki);

        
        if(!m_ptr_nh->getParam(tmp_str+"/Kd",m_internal_pid_param.Kd))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Kd=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Kd);

        if(!m_ptr_nh->getParam(tmp_str+"/Kaw",m_internal_pid_param.Kaw))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Kaw=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Kaw);

        if(!m_ptr_nh->getParam(tmp_str+"/Td",m_internal_pid_param.Td))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Td=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Td);

        if(!m_ptr_nh->getParam(tmp_str+"/Tf",m_internal_pid_param.Tf))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        } 
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Tf=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Tf);

        if(!m_ptr_nh->getParam(tmp_str+"/Umin",m_internal_pid_param.Umin))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Umin=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Umin);

        if(!m_ptr_nh->getParam(tmp_str+"/Umax",m_internal_pid_param.Umax))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param Umax=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.Umax);

        if(!m_ptr_nh->getParam(tmp_str+"/DUmin",m_internal_pid_param.DUmin))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param DUmin=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.DUmin);

        if(!m_ptr_nh->getParam(tmp_str+"/DUmax",m_internal_pid_param.DUmax))
        {
          ROS_WARN( " [ %s%s:%d%s ]\t call /pid_update_param FAIL!!",GREEN, __FUNCFILE__, __LINE__, RESET);
          return res.res = (false);
        }
        else
          ROS_INFO( " [ %s%s:%d%s ]\t %s/pid_param DUmax=%5.4f!!",GREEN, __FUNCFILE__, __LINE__, RESET,tmp_str.c_str(),m_internal_pid_param.DUmax);
        
        vel_error_filter.first_ord_LP_filt(m_internal_pid_param.Td, m_internal_pid_param.Kd, m_internal_pid_param.st, 0.0, 0.0);
        control_action_filter.first_ord_LP_filt(m_internal_pid_param.Tf,1.0, m_internal_pid_param.st, 0.0, 0.0);
        return res.res = (true);
      }

    private: 

        ros::NodeHandle*  m_ptr_nh;
        double            m_setpoint;
        
        pidPar            m_internal_pid_param;

        double control_action_past;
        double Istate; //integrator state
        double control_action;

        itia::filter::tf  vel_error_filter;
        itia::filter::tf  control_action_filter;
        ros::Time         last_time;
        
};    

}
}
#endif




