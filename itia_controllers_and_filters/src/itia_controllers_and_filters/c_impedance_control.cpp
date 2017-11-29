#include <itia_controllers_and_filters/c_impedance_control.h>
#include <boost/concept_check.hpp>

namespace itia {
namespace controller {
c_impedance_control::c_impedance_control ( const std::string& name ) {

    cont1 = 50;
    cont2 = 50;
    cont3 = 0;
    
    imp_mode_t = 1;
    imp_mode_r = 1;
    
    for (int ii=0; ii<3; ii++)
    {
        phi_cd(ii) = 0.; 
        phi_cd_old(ii) = 0.;
        force_start(ii) = 0.;
        force_start(ii+3) = 0.;
        xDev_psp(ii) = 0.;
        xdDev_psp(ii) = 0.;
        xddDev_psp(ii) = 0.;
        start_pos_t(ii) = 0.;
        start_pos_r(ii) = 0.;
        delta_w_psp_old(ii) = 0.;
        delta_r_psp(ii) = 0.;
        phi_cd_cmd(ii) = 0.;
        dphi_cd_cmd(ii) = 0.; 
        ddphi_cd_cmd(ii) = 0.; 
        delta_w_psp(ii) = 0.;
        pos_t_sp(ii) = 0.;
        pos_r_sp(ii) = 0.;
        pos_t_msr(ii) = 0.;
        pos_r_msr(ii) = 0.;
        deltapos_t(ii) = 0.;
        vel_t_sp(ii) = 0.;
        vel_t_msr(ii) = 0.;
        deltavel_t(ii) = 0.;
        torque_filt(ii) = 0.;
        force_filt(ii) = 0.;
    }
    
    initTime=ros::Time::now();
    maxDevLinVel=0.1;
    maxDevAngVel=0.1;
    maxDevLinPos=0.02;
    maxDevAngPos=0.02;
    forceDeadBand=0;
    saturation=false;
    min_mass=0.01;
    forceFiltTime=0.1;
    force_filt.setZero();
    firstEntry=true;
    firstEntry_psp=true;

    q.setIdentity();
    q_psp.setIdentity();
    bool welldefined=true;
    int tmp;
    if ( !nh.getParam ( name+"/angle_convection",tmp ) )
        angleConvection=EULER_ZYX;
    else
        angleConvection= ( AngleConvection ) tmp;


    std::vector<double> stiff;
    if ( !nh.getParam ( name+"/stiff",stiff ) )
        ROS_ERROR ( "impedance control parameters '%s/stiff' not found\n",name.c_str() );
    if ( stiff.size() !=6 ) {
        ROS_ERROR ( "Stiffness parameters have wrong dimensions\n" );
        welldefined=false;
    }

    std::vector<double> damp;
    if ( !nh.getParam ( name+"/damp",damp ) )
        ROS_ERROR ( "impedance control parameters not found\n" );
    if ( damp.size() !=6 ) {
        ROS_ERROR ( "damping parameters have wrong dimensions\n" );
        welldefined=false;
    }

    std::vector<double> mass;
    if ( !nh.getParam ( name+"/mass",mass ) )
        ROS_ERROR ( "impedance control parameters not found\n" );
    if ( mass.size() !=6 ) {
        ROS_ERROR ( "mass parameters have wrong dimensions\n" );
        welldefined=false;
    }

    Kp.setZero();
    Kd.setZero();
    M.setZero();
    Minv.setZero();

    double dparam;
    if ( nh.getParam ( name+"/maxDevLinPos",dparam ) )
        maxDevLinPos=dparam;
    if ( nh.getParam ( name+"/maxDevAngPos",dparam ) )
        maxDevAngPos=dparam;
    if ( nh.getParam ( name+"/maxDevLinVel",dparam ) )
        maxDevLinVel=dparam;
    if ( nh.getParam ( name+"/maxDevAngVel",dparam ) )
        maxDevAngVel=dparam;
    if ( nh.getParam ( name+"/forceFiltTime",dparam ) )
        forceFiltTime=dparam;
    if ( nh.getParam ( name+"/forceDeadBand",dparam ) )
        if ( dparam>0 )
            forceDeadBand=dparam;
        else
            forceDeadBand=0;
    if ( nh.getParam ( name+"/torqueDeadBand",dparam ) )
        if ( dparam>0 )
            torqueDeadBand=dparam;
        else
            forceDeadBand=0;

    if ( welldefined )
        set_parameters ( mass,stiff,damp );

    impPoseSp_pub=nh.advertise<geometry_msgs::PoseStamped> ( "imp_setpoint",1);
};

void c_impedance_control::set_parameters ( const std::vector<double>& mass,
        const std::vector<double>& stiff,
        const std::vector<double>& damp ) {
    for ( int idx=0; idx<6; idx++ ) {
      
        if ( stiff.at ( idx ) <0.0 )
            Kp ( idx,idx ) =0.0;
        else
            Kp ( idx,idx ) =stiff.at ( idx );

        if ( mass.at ( idx ) <min_mass )
            M ( idx,idx ) =min_mass;
        else
            M ( idx,idx ) =mass.at ( idx );
        
        if ( damp.at ( idx ) <0.0 )
            Kd ( idx,idx ) =0.0;
        else
            Kd ( idx,idx ) =2. * damp.at ( idx ) * M ( idx,idx ) * sqrt( Kp ( idx,idx ) / M ( idx,idx ) ) ;

        std::cout << "iAx = "<<idx<<", mass = "<<M ( idx,idx ) <<", spring = "<<Kp ( idx,idx ) <<", damping = "<<Kd ( idx,idx ) <<std::endl;
        
    }
    
    Minv=M.inverse();
    
    std::cout << "Minv: " << Minv << std::endl;
}

void c_impedance_control::integration ( const double& timeInterval ) {
    //tra
    xDev+=timeInterval* ( xdDev+ xddDev*timeInterval/2.0 );
    xdDev+=xddDev*timeInterval;

    //rot
    Eigen::Quaterniond omega ( 0,w ( 0 ),w ( 1 ),w ( 2 ) );
    Eigen::Quaterniond dOmega ( 0,dw ( 0 ),dw ( 1 ),dw ( 2 ) );
    Eigen::Quaterniond dq, ddq, dddq, ddddq;
    dq.coeffs()    =0.5* ( omega*q ).coeffs();
    ddq.coeffs()   =0.5* ( ( dOmega*q ).coeffs()  + ( omega*dq ).coeffs() );
    dddq.coeffs()  =0.5* ( 2* ( dOmega*dq ).coeffs() + ( omega*ddq ).coeffs() );
    ddddq.coeffs() =0.5* ( 3* ( dOmega*ddq ).coeffs() + ( omega*dddq ).coeffs() );
    q.coeffs() +=dq.coeffs() *timeInterval+
                 ddq.coeffs() *pow ( timeInterval,2.0 ) /2.0+
                 dddq.coeffs() *pow ( timeInterval,3.0 ) /6.0+
                 ddddq.coeffs() *pow ( timeInterval,4.0 ) /24.0;
    q.normalize();
    w+=dw*timeInterval;

}

KDL::FrameAcc c_impedance_control::get_control_action ( const KDL::FrameAcc& frameSp,
        const KDL::Wrench& forceSp, //tool frame
        const KDL::Wrench& force, //tool frame
        const double& timeInterval ) {

    double digfilt=0.0;
    if ( forceFiltTime>0 )
        digfilt=exp ( -timeInterval/forceFiltTime );

    if ( firstEntry ) {
        firstEntry=false;
        frameCom=frameSp;
        q.setIdentity();
    }


    // la forza e gli sostamenti sono in terna frameCom del passo prima
    KDL::Wrench deltaForce= ( force-forceSp );

//                      double sinu=10*sin(2*M_PI/10*(ros::Time::now()-initTime).toSec());
//                      KDL::Wrench deltaForceBase=KDL::Wrench(KDL::Vector( sinu,0,0 ),
//                                                                                                                                                                               KDL::Vector(0,0,10));
//                      KDL::Wrench deltaForce=frameCom.M.R.Inverse()*deltaForceBase;


    for ( int idx=0; idx<3; idx++ ) {
        if ( deltaForce.force ( idx ) >forceDeadBand )
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( deltaForce.force ( idx )-forceDeadBand );
        else if ( deltaForce.force ( idx ) <-forceDeadBand )
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( deltaForce.force ( idx ) +forceDeadBand );
        else
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt );

        if ( deltaForce.torque ( idx ) >torqueDeadBand )
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( deltaForce.torque ( idx )-torqueDeadBand );
        else if ( deltaForce.torque ( idx ) <-torqueDeadBand )
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( deltaForce.torque ( idx ) +torqueDeadBand );
        else
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt );
    }

    Eigen::Matrix<double,6,1> acc,vel,pos,genForce;
    genForce.block ( 0,0,3,1 ) =force_filt;
    genForce.block ( 3,0,3,1 ) =torque_filt;

    for ( int idx=0; idx<3; idx++ ) {
        vel ( idx ) =xdDev ( idx );
        pos ( idx ) =xDev ( idx );
        vel ( idx+3 ) =w ( idx );
    }

    KDL::Rotation R_sp, R_msr;

    R_sp = frameSp.M.R;

    KDL::Rotation rotDev=KDL::Rotation::Quaternion ( q.x(),q.y(),q.z(),q.w() );
    if ( angleConvection==c_impedance_control::EULER_ZYX )
        rotDev.GetEulerZYX ( pos ( 5 ),pos ( 4 ),pos ( 3 ) );
    else if ( angleConvection==c_impedance_control::RPY )
        rotDev.GetRPY ( pos ( 3 ),pos ( 4 ),pos ( 5 ) );
    else
        rotDev.GetEulerZYX ( pos ( 5 ),pos ( 4 ),pos ( 3 ) );

    acc=Minv* ( ( genForce ) - ( Kp*pos + Kd*vel ) );

    for ( int idx=0; idx<3; idx++ ) {
        xddDev ( idx ) =acc ( idx );
        dw ( idx ) =acc ( idx+3 );
    }
//                      std::cout.precision(4);
//                      std::cout << "genForce " <<std::setw(8)<< genForce.transpose()<<std::endl;
//                      std::cout << "acc      " <<std::setw(8)<< acc.transpose()<<std::endl;
//                      std::cout << "vel      " <<std::setw(8)<< vel.transpose()<<std::endl;
//                      std::cout << "pos      " <<std::setw(8)<< pos.transpose()<<std::endl;
    integration ( timeInterval );


    frameDev.M.R=KDL::Rotation::Quaternion ( q.x(),q.y(),q.z(),q.w() );
    frameDev.M.w=w;
    frameDev.M.dw=dw;

    frameDev.p.p=xDev;
    frameDev.p.v=xdDev;
    frameDev.p.dv=xddDev;

//                      frameCom.M.w =cImpSp.M.w+  frameCom.M.R*frameDev.M.w;
//                      frameCom.M.dw=cImpSp.M.dw+ frameCom.M.R*frameDev.M.dw;
//
//                      frameCom.p.v =cImpSp.p.v+  frameCom.M.R*frameDev.p.v;
//                      frameCom.p.dv=cImpSp.p.dv+ frameCom.M.R*frameDev.p.dv;
//
//                      KDL::Rotation devKtoDevKpp=((cImpSp.M.R*frameDev.M.R).Inverse()*frameCom.M.R).Inverse(); //rotazione che manda le deviazioni al passo k al sistema di riferimento del passo k+1
//                      frameCom.M.R =cImpSp.M.R*frameDev.M.R;
//                      frameCom.p.p =cImpSp.p.p*frameDev.p.p;

    frameCom=frameSp*frameDev;
//                      KDL::Rotation devKtoDevKpp=frameDev.M.R.Inverse(); //rotazione che manda le deviazioni al passo k al sistema di riferimento del passo k+1
// //                   KDL::Rotation devKtoDevKpp=KDL::Rotation::Identity(); //rotazione che manda le deviazioni al passo k al sistema di riferimento del passo k+1
//                      // le deviazioni espresse in nella nuova terna
//                      xDev=devKtoDevKpp*xDev;
//                      xdDev=devKtoDevKpp*xdDev;
//                      xddDev=devKtoDevKpp*xddDev;
//                      w=devKtoDevKpp*w;
//                      dw=devKtoDevKpp*dw;
//                      //che fare delle rotazioni?????

    KDL::Rotation R_tmp;

    KDL::Vector p_tmp;

    KDL::Rotation R_dev;

    KDL::Vector p_dev;

    p_dev = frameDev.p.p;

    R_dev = frameDev.M.R;

    R_tmp = frameCom.M.R;

    p_tmp = frameCom.p.p;

//     if ( cont1++ >= 50 ) {
// 
//         cont1 = 0;
// 
//         std::cout << "1' method" << std::endl;
// 
//         std::cout << "timeInterval: " << timeInterval << std::endl;
// 
//         std::cout << "q: " << q.x() << "; " << q.y() << "; " << q.z() << "; " << q.w() << "; " << std::endl;
// 
//         std::cout << "P_cmd: " << p_tmp ( 0 ) << "; " << p_tmp ( 1 ) << "; " << p_tmp ( 2 ) << std::endl;
// 
//         std::cout << "p_dev: " << p_dev ( 0 ) << "; " << p_dev ( 1 ) << "; " << p_dev ( 2 ) << std::endl;
// 
//         std::cout << "R_cmd: " << R_tmp ( 0,0 ) << "; " << R_tmp ( 0,1 ) << "; " << R_tmp ( 0,2 ) << "; " << R_tmp ( 1,0 ) << "; "
//                   << R_tmp ( 1,1 ) << "; " << R_tmp ( 1,2 ) << "; " << R_tmp ( 2,0 ) << "; " << R_tmp ( 2,1 ) << "; " << R_tmp ( 2,2 ) << std::endl;
// 
//         std::cout << "R_dev: " << R_dev ( 0,0 ) << "; " << R_dev ( 0,1 ) << "; " << R_dev ( 0,2 ) << "; " << R_dev ( 1,0 ) << "; "
//                   << R_dev ( 1,1 ) << "; " << R_dev ( 1,2 ) << "; " << R_dev ( 2,0 ) << "; " << R_dev ( 2,1 ) << "; " << R_dev ( 2,2 ) << std::endl;
// 
//         /*std::cout << "R_sp: " << R_sp(0,0) << "; " << R_sp(0,1) << "; " << R_sp(0,2) << "; " << R_sp(1,0) << "; "
//                   << R_sp(1,1) << "; " << R_sp(1,2) << "; " << R_sp(2,0) << "; " << R_sp(2,1) << "; " << R_sp(2,2) << std::endl;
// 
//         std::cout << "R_msr: " << R_msr(0,0) << "; " << R_msr(0,1) << "; " << R_msr(0,2) << "; " << R_msr(1,0) << "; "
//                   << R_msr(1,1) << "; " << R_msr(1,2) << "; " << R_msr(2,0) << "; " << R_msr(2,1) << "; " << R_msr(2,2) << std::endl;
//            */
// 
//         std::cout << "-------------------------------" << std::endl;
// 
// //   std::cout << "cmdV " << frameCom.p.v.x()<< " "<< frameCom.p.v.y()<< " "<< frameCom.p.v.z()<< " "<<std::endl;
// 
//     }

    return frameCom;

}

KDL::FrameAcc c_impedance_control::get_control_action_position_sp ( const KDL::FrameAcc cartFrame_Sp,
        const KDL::FrameVel cartFrame_msr,
        const KDL::Wrench& force, //end-effector frame
        const double& timeInterval,
        const bool& setfirstEntry_psp,
        const double& add_damp ) {
    
    KDL::Frame TooltoSetPoint;
    Eigen::Matrix<double,6,1> acc, Force_act, deltavel;
    KDL::Rotation R_compliant_to_sp, R_compliant_to_sp_start, T_dangles_to_w_cd, R_msr, R_sp, T_dangles_to_w_delta_r_psp;
    KDL::Vector delta_r_psp_acc, delta_r_psp_vel, new_r;
    KDL::Wrench msrForce_tool, wrench_tmp;
    Eigen::VectorXd msrForce ( 3 );
    Eigen::VectorXd msrTorque_tool ( 3 );
    Eigen::VectorXd msrForceTorque_ImpedanceControl ( 6 ); // forces in base frame and torques in setpoint frame
    Eigen::VectorXd msrforce_Baseframe ( 3 );
    Eigen::Matrix<double,3,3> HoleToBaseRotationMatrix;
    Eigen::Matrix<double,3,3> ToolToBaseRotationMatrix;
    Eigen::Matrix<double,3,3> M_r_inv, M_t_inv, Kp_t, Kp_r, Kd_t, Kd_r, T_dangles_to_w_cd_eigen;
    Eigen::Matrix<double,3,1> ddphi_cd_eigen, acc_t, Forces, Torques, deltapos_t_eigen, deltavel_t_eigen, phi_cd_eigen, dphi_cd_eigen;
    KDL::Rotation R_new, R_cmd;
    KDL::Vector p_tmp, p_dev;
    
    acc_t.setZero();
    
    firstEntry_psp = setfirstEntry_psp;

    double digfilt=0.0;
    if ( forceFiltTime>0 )
        digfilt=exp ( -timeInterval/forceFiltTime );
    
    if ( firstEntry_psp ) {
      
        if ( !nh.getParam ( "ur_cimp_psp/domega_min_psp_deg_s2",m_domega_min_psp_deg_s2 ) ) {
          ROS_WARN ( "angular acc min position setpoint topic name not found!\n" );
          m_domega_min_psp_deg_s2=1.;
        }
        
        for (int ii=0; ii<3; ii++)
        {
            phi_cd(ii) = 0.; 
            phi_cd_old(ii) = 0.;
            force_start(ii) = 0.;
            force_start(ii+3) = 0.;
            xDev_psp(ii) = 0.;
            xdDev_psp(ii) = 0.;
            xddDev_psp(ii) = 0.;
            start_pos_t(ii) = 0.;
            start_pos_r(ii) = 0.;
            delta_w_psp_old(ii) = 0.;
            delta_r_psp(ii) = 0.;
            phi_cd_cmd(ii) = 0.;
            dphi_cd_cmd(ii) = 0.; 
            ddphi_cd_cmd(ii) = 0.; 
            delta_w_psp(ii) = 0.;
            pos_t_sp(ii) = 0.;
            pos_r_sp(ii) = 0.;
            pos_t_msr(ii) = 0.;
            pos_r_msr(ii) = 0.;
            deltapos_t(ii) = 0.;
            vel_t_sp(ii) = 0.;
            vel_t_msr(ii) = 0.;
            deltavel_t(ii) = 0.;
        }
        
        force_start(0) = force(0);
        force_start(1) = force(1);
        force_start(2) = force(2);
        force_start(3) = force(3);
        force_start(4) = force(4);
        force_start(5) = force(5);
        
        force_filt.setZero();

        start_pos_t = cartFrame_msr.p.p;
        
        R_start = cartFrame_msr.M.R;
        
        if ( angleConvection==c_impedance_control::EULER_ZYX )
            R_start.GetEulerZYX ( start_pos_r ( 2 ),start_pos_r ( 1 ),start_pos_r ( 0 ) );
        else if ( angleConvection==c_impedance_control::RPY )
            R_start.GetRPY ( start_pos_r ( 0 ),start_pos_r ( 1 ),start_pos_r ( 2 ) ); //alpha,beta,gamma
        else
            R_start.GetEulerZYX ( start_pos_r ( 2 ),start_pos_r ( 1 ),start_pos_r ( 0 ) );
        
        frameCom=cartFrame_Sp;
        q_psp.setIdentity();
    }

    pos_t_sp = cartFrame_Sp.p.p;
    pos_t_msr = cartFrame_msr.p.p;
    
    vel_t_sp = cartFrame_Sp.p.v;
    vel_r_sp = cartFrame_Sp.M.w;

    vel_t_msr = cartFrame_msr.p.v;
    
    R_sp = cartFrame_Sp.M.R;
    R_msr = cartFrame_msr.M.R;
    
    R_compliant_to_sp = R_sp.Inverse()*R_msr;
    
    R_compliant_to_sp_start = R_sp.Inverse()*R_start;
    
    TooltoSetPoint.p = cartFrame_Sp.p.p - cartFrame_Sp.p.p;  
//     TooltoSetPoint.M = R_compliant_to_sp;
    TooltoSetPoint.M = R_msr;
//     TooltoSetPoint.p = cartFrame_Sp.p.p - cartFrame_msr.p.p;
//     TooltoSetPoint.p = cartFrame_msr.p.p;
    
    msrForce_tool = force - force_start;
    
    msrForce ( 0 ) = msrForce_tool.force ( 0 );
    msrForce ( 1 ) = msrForce_tool.force ( 1 );
    msrForce ( 2 ) = msrForce_tool.force ( 2 );
    
    msrTorque_tool ( 0 ) = msrForce_tool.torque ( 0 );
    msrTorque_tool ( 1 ) = msrForce_tool.torque ( 1 );
    msrTorque_tool ( 2 ) = msrForce_tool.torque ( 2 );
    
    int s_1,s_2;
    
    s_1 = 500;
    s_2 = 3000;
    
    cont2++;
    cont3++;
    
    bool try_torque = false;
    
    if (try_torque == true)
    { 
        if (cont3<s_1)
        {
            msrForce_tool.torque ( 2 ) = 0.;
        }
        else if (cont3<s_2)
        {
            msrForce_tool.torque ( 2 ) = 5.;
        }
        else if (cont3==s_2)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "-------------------------" << std::endl;
            std::cout << "-------------------------" << std::endl;
            std::cout << "stop applying torque/force" << std::endl;
            std::cout << "-------------------------" << std::endl;
            std::cout << "-------------------------" << std::endl;
            std::cout << "-------------------------" << std::endl;
        }
    } 
    
    for ( int ii = 0; ii<3; ii++ ) {
        for ( int jj = 0; jj<3; jj++ ) {
            ToolToBaseRotationMatrix ( ii,jj ) = R_msr ( ii,jj );
        }
    }
    
    msrforce_Baseframe = ToolToBaseRotationMatrix * msrForce;
    
    wrench_tmp.force ( 0 ) = msrForce_tool.force ( 0 );
    wrench_tmp.force ( 1 ) = msrForce_tool.force ( 1 );
    wrench_tmp.force ( 2 ) = msrForce_tool.force ( 2 );
    wrench_tmp.torque ( 0 ) = msrForce_tool.torque ( 0 );
    wrench_tmp.torque ( 1 ) = msrForce_tool.torque ( 1 );
    wrench_tmp.torque ( 2 ) = msrForce_tool.torque ( 2 );
    
    wrench_tmp=wrench_tmp.RefPoint(TooltoSetPoint.p);
    wrench_tmp=TooltoSetPoint.M*wrench_tmp;
    
    msrForceTorque_ImpedanceControl ( 0 ) = msrforce_Baseframe ( 0 );
    msrForceTorque_ImpedanceControl ( 1 ) = msrforce_Baseframe ( 1 );
    msrForceTorque_ImpedanceControl ( 2 ) = msrforce_Baseframe ( 2 );
    msrForceTorque_ImpedanceControl ( 3 ) = wrench_tmp.torque ( 0 );
    msrForceTorque_ImpedanceControl ( 4 ) = wrench_tmp.torque ( 1 );
    msrForceTorque_ImpedanceControl ( 5 ) = wrench_tmp.torque ( 2 );

    double force_max = 100;
    double torque_max = 20;
    
    for ( int idx=0; idx<3; idx++ ) {
        if ( msrForceTorque_ImpedanceControl ( idx ) >forceDeadBand )
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( msrForceTorque_ImpedanceControl ( idx )-forceDeadBand );
        else if ( msrForceTorque_ImpedanceControl ( idx ) <-forceDeadBand )
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( msrForceTorque_ImpedanceControl ( idx ) +forceDeadBand );
        if ( msrForceTorque_ImpedanceControl ( idx ) >force_max )
            force_filt ( idx ) =force_max;
        if ( msrForceTorque_ImpedanceControl ( idx ) <-force_max )
            force_filt ( idx ) =-force_max;
        else if ( msrForceTorque_ImpedanceControl ( idx ) > torque_max )
            force_filt ( idx ) =torque_max;
        else if ( msrForceTorque_ImpedanceControl ( idx ) < -torque_max )
            force_filt ( idx ) =-torque_max; 
        else
            force_filt ( idx ) =force_filt ( idx ) * ( digfilt );

        if ( msrForceTorque_ImpedanceControl ( idx+3 ) >torqueDeadBand )
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( msrForceTorque_ImpedanceControl ( idx+3 )-torqueDeadBand );
        else if ( msrForceTorque_ImpedanceControl ( idx+3 ) <-torqueDeadBand )
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt ) + ( 1-digfilt ) * ( msrForceTorque_ImpedanceControl ( idx+3 ) +torqueDeadBand );
        else
            torque_filt ( idx ) =torque_filt ( idx ) * ( digfilt );
    }
        
    Force_act.block ( 0,0,3,1 ) = force_filt;
    Force_act.block ( 3,0,3,1 ) = torque_filt;

    if ( angleConvection==c_impedance_control::EULER_ZYX ) {
        
        R_sp.GetEulerZYX ( pos_r_sp ( 2 ),pos_r_sp ( 1 ),pos_r_sp ( 0 ) );
        R_msr.GetEulerZYX ( pos_r_msr ( 2 ),pos_r_msr ( 1 ),pos_r_msr ( 0 ) );
        R_compliant_to_sp.GetEulerZYX ( phi_cd ( 2 ),phi_cd ( 1 ),phi_cd ( 0 ) ); //alpha, beta, gamma
        
        if (imp_mode_r == 1)
        {
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd_cmd(2) + phi_cd_start(2));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd_cmd(2) + phi_cd_start(2))*cos(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd_cmd(2) + phi_cd_start(2));
            T_dangles_to_w_cd(2,2) = cos(phi_cd_cmd(2) + phi_cd_start(2))*cos(phi_cd_cmd(1) + phi_cd_start(1));
        }
        else
        {
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd(2));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd(2))*cos(phi_cd(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd(2));
            T_dangles_to_w_cd(2,2) = cos(phi_cd(2))*cos(phi_cd(1));
        }
        
    } else if ( angleConvection==c_impedance_control::RPY ) {
        
        R_sp.GetRPY ( pos_r_sp ( 0 ),pos_r_sp ( 1 ),pos_r_sp ( 2 ) );
        R_msr.GetRPY ( pos_r_msr ( 0 ),pos_r_msr ( 1 ),pos_r_msr ( 2 ) );
        R_compliant_to_sp.GetRPY ( phi_cd ( 0 ),phi_cd ( 1 ),phi_cd ( 2 ) ); //alpha, beta, gamma
        
        if (imp_mode_r == 1)
        {
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd_cmd(0) + phi_cd_start(0));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd_cmd(0) + phi_cd_start(0))*cos(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd_cmd(0) + phi_cd_start(0));
            T_dangles_to_w_cd(2,2) = cos(phi_cd_cmd(0) + phi_cd_start(0))*cos(phi_cd_cmd(1) + phi_cd_start(1));
        }
        else
        {
            
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd(0));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd(0))*cos(phi_cd(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd(0));
            T_dangles_to_w_cd(2,2) = cos(phi_cd(0))*cos(phi_cd(1));
            
        }
        
    } else {
        
        R_sp.GetEulerZYX ( pos_r_sp ( 2 ),pos_r_sp ( 1 ),pos_r_sp ( 0 ) );
        R_msr.GetEulerZYX ( pos_r_msr ( 2 ),pos_r_msr ( 1 ),pos_r_msr ( 0 ) );
        R_compliant_to_sp.GetEulerZYX ( phi_cd ( 2 ),phi_cd ( 1 ),phi_cd ( 0 ) ); //alpha, beta, gamma
        
       if (imp_mode_r == 1)
        {
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd_cmd(2) + phi_cd_start(2));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd_cmd(2) + phi_cd_start(2))*cos(phi_cd_cmd(1) + phi_cd_start(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd_cmd(2) + phi_cd_start(2));
            T_dangles_to_w_cd(2,2) = cos(phi_cd_cmd(2) + phi_cd_start(2))*cos(phi_cd_cmd(1) + phi_cd_start(1));
        }
        else
        {
            T_dangles_to_w_cd(0,0) = 1;
            T_dangles_to_w_cd(0,1) = 0.;
            T_dangles_to_w_cd(0,2) = sin(phi_cd(1));
            T_dangles_to_w_cd(1,0) = 0.;
            T_dangles_to_w_cd(1,1) = cos(phi_cd(2));
            T_dangles_to_w_cd(1,2) = -sin(phi_cd(2))*cos(phi_cd(1));
            T_dangles_to_w_cd(2,0) = 0.;
            T_dangles_to_w_cd(2,1) = sin(phi_cd(2));
            T_dangles_to_w_cd(2,2) = cos(phi_cd(2))*cos(phi_cd(1));
        }
        
    }
    
    //     if ( angleConvection==c_impedance_control::EULER_ZYX )
    //         R_compliant_to_sp_start.GetEulerZYX ( phi_cd_start ( 2 ),phi_cd_start ( 1 ),phi_cd_start ( 0 ) );
    //     else if ( angleConvection==c_impedance_control::RPY )
    //         R_compliant_to_sp_start.GetRPY ( phi_cd_start ( 0 ),phi_cd_start ( 1 ),phi_cd_start ( 2 ) ); //alpha,beta,gamma
    //     else
    //         R_compliant_to_sp_start.GetEulerZYX ( phi_cd_start ( 2 ),phi_cd_start ( 1 ),phi_cd_start ( 0 ) );
    
    phi_cd_start = - pos_r_sp + start_pos_r;
    
    for ( int ii=0; ii<3; ii++ ) {
    
        if (imp_mode_t == 1)
        {
            deltavel_t(ii) = vel_t_sp(ii) - xdDev_psp(ii);
            deltapos_t(ii) = pos_t_sp(ii) - (xDev_psp(ii)+start_pos_t(ii));
        }
        else if (imp_mode_t == 2)
        {
            deltavel_t(ii) = vel_t_sp(ii) - vel_t_msr(ii);
            deltapos_t(ii) = pos_t_sp(ii) - pos_t_msr(ii);
        }
        
        if (imp_mode_r == 1)
        {
            phi_cd_eigen(ii) = - ( phi_cd_cmd(ii) + phi_cd_start(ii) ); //
            dphi_cd_eigen(ii) = - ( dphi_cd_cmd(ii) ); //
        }
        else if (imp_mode_r == 2)
        {
            dphi_cd_eigen(ii) = dphi_cd_cmd ( ii );//( phi_cd ( ii ) - phi_cd_old(ii) ) / timeInterval; // da verificare
            phi_cd_eigen(ii) = phi_cd ( ii ); // da verificare
            phi_cd_old(ii) = phi_cd ( ii );
        }
        
    }
    
    for (int ii=0; ii<3; ii++)
    {
        
        for (int jj=0; jj<3; jj++)
        {
            
            M_r_inv(ii,jj) = Minv(3+ii,3+jj);
            M_t_inv(ii,jj) = Minv(ii,jj);
            
            Kp_t(ii,jj) = Kp(ii,jj);
            Kp_r(ii,jj) = Kp(ii+3,jj+3);
            
            Kd_t(ii,jj) = Kd(ii,jj) + 2. * add_damp * M ( ii,jj ) * sqrt( Kp_t ( ii,jj ) * M_t_inv( ii,jj ) ) ;
            Kd_r(ii,jj) = Kd(ii+3,jj+3);
            
            T_dangles_to_w_cd_eigen(ii,jj) = T_dangles_to_w_cd(ii,jj);
            
        }
        
        deltapos_t_eigen(ii) = deltapos_t(ii);
        deltavel_t_eigen(ii) = deltavel_t(ii);
        
        Forces(ii) = Force_act(ii);
        Torques(ii) = Force_act(3+ii);
        
    }
    
    acc_t = M_t_inv * ( Forces + ( Kp_t*deltapos_t_eigen + Kd_t*deltavel_t_eigen ) );
    
    Eigen::Matrix<double,3,1> T_plot;
    T_plot = T_dangles_to_w_cd_eigen.transpose() * (Torques);
    
    if (imp_mode_r == 1)
    {
        ddphi_cd_eigen = M_r_inv * ( T_dangles_to_w_cd_eigen.transpose() * (Torques) + ( Kp_r*phi_cd_eigen + Kd_r*dphi_cd_eigen ) );
    }
    else
    {
        ddphi_cd_eigen = M_r_inv * ( T_dangles_to_w_cd_eigen.transpose() * (Torques) - ( Kp_r*phi_cd_eigen + Kd_r*dphi_cd_eigen ) );
    }
        
    for ( int idx=0; idx<3; idx++ ) {
        xddDev_psp(idx) = acc_t(idx);
        ddphi_cd_cmd(idx) = ddphi_cd_eigen(idx);
    }

    if (imp_mode_t == 1)
    {
        for (int ii=0; ii<3; ii++)
        {
            xdDev_psp(ii)+=xddDev_psp(ii)*timeInterval;
            xDev_psp(ii)+=timeInterval* ( xdDev_psp(ii) + xddDev_psp(ii)*timeInterval/2.0 );
        }
    }
    else if (imp_mode_t == 2)
    {
        for (int ii=0; ii<3; ii++)
        {
            xdDev_psp(ii)=xddDev_psp(ii)*timeInterval;
            xDev_psp(ii)=timeInterval* ( xdDev_psp(ii) + xddDev_psp(ii)*timeInterval/2.0 );
        }
    }
    
    if (imp_mode_r == 1)
    {
        for (int ii=0; ii<3; ii++)
        {
            dphi_cd_cmd(ii) += ddphi_cd_cmd(ii)*timeInterval;
            phi_cd_cmd(ii) += timeInterval* ( dphi_cd_cmd(ii) + ddphi_cd_cmd(ii)*timeInterval/2.0 );
            
            delta_r_psp(ii) = phi_cd_cmd(ii) + phi_cd_start(ii);
            delta_r_psp_vel(ii) = dphi_cd_cmd(ii);
        }
    }
    else if (imp_mode_r == 2)
    {
        for (int ii=0; ii<3; ii++)
        {
            dphi_cd_cmd(ii) = ddphi_cd_cmd(ii)*timeInterval;
            phi_cd_cmd(ii) = timeInterval* ( dphi_cd_cmd(ii) + ddphi_cd_cmd(ii)*timeInterval/2.0 );
            
            delta_r_psp(ii) = phi_cd_cmd(ii);
            delta_r_psp_vel(ii) = dphi_cd_cmd(ii);
        }
    }
    
    if ( angleConvection==c_impedance_control::EULER_ZYX ) {
        
        for (int ii=0; ii<3; ii++)
        {
            delta_r_psp_acc(ii) = ddphi_cd_cmd(ii);
        }
        
        T_dangles_to_w_delta_r_psp(0,0) = 1;
        T_dangles_to_w_delta_r_psp(0,1) = 0.;
        T_dangles_to_w_delta_r_psp(0,2) = sin(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(1,0) = 0.;
        T_dangles_to_w_delta_r_psp(1,1) = cos(delta_r_psp ( 2 ));
        T_dangles_to_w_delta_r_psp(1,2) = -sin(delta_r_psp ( 2 ))*cos(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(2,0) = 0.;
        T_dangles_to_w_delta_r_psp(2,1) = sin(delta_r_psp ( 2 ));
        T_dangles_to_w_delta_r_psp(2,2) = cos(delta_r_psp ( 2 ))*cos(delta_r_psp ( 1 ));
        
        delta_w_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_vel;
        
//         dw_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_acc;
        
        for (int ii=0; ii<3; ii++)
        {
            dw_psp(ii) = (delta_w_psp(ii) - delta_w_psp_old(ii)) / timeInterval;
            delta_w_psp_old(ii) = delta_w_psp(ii);
        }
        
    } else if ( angleConvection==c_impedance_control::RPY ) {
        
        for (int ii=0; ii<3; ii++)
        {
            delta_r_psp_acc(ii) = ddphi_cd_cmd(ii);
        }
        
        T_dangles_to_w_delta_r_psp(0,0) = 1;
        T_dangles_to_w_delta_r_psp(0,1) = 0.;
        T_dangles_to_w_delta_r_psp(0,2) = sin(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(1,0) = 0.;
        T_dangles_to_w_delta_r_psp(1,1) = cos(delta_r_psp ( 0 ));
        T_dangles_to_w_delta_r_psp(1,2) = -sin(delta_r_psp ( 0 ))*cos(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(2,0) = 0.;
        T_dangles_to_w_delta_r_psp(2,1) = sin(delta_r_psp ( 0 ));
        T_dangles_to_w_delta_r_psp(2,2) = cos(delta_r_psp ( 0 ))*cos(delta_r_psp ( 1 ));
        
        delta_w_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_vel;
        
//         dw_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_acc;
        
        for (int ii=0; ii<3; ii++)
        {
            dw_psp(ii) = (delta_w_psp(ii) - delta_w_psp_old(ii)) / timeInterval;
            delta_w_psp_old(ii) = delta_w_psp(ii);
        }
        
    } else {
        
        for (int ii=0; ii<3; ii++)
        {
            delta_r_psp_acc(ii) = ddphi_cd_cmd(ii);
        }
        
        T_dangles_to_w_delta_r_psp(0,0) = 1;
        T_dangles_to_w_delta_r_psp(0,1) = 0.;
        T_dangles_to_w_delta_r_psp(0,2) = sin(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(1,0) = 0.;
        T_dangles_to_w_delta_r_psp(1,1) = cos(delta_r_psp ( 2 ));
        T_dangles_to_w_delta_r_psp(1,2) = -sin(delta_r_psp ( 2 ))*cos(delta_r_psp ( 1 ));
        T_dangles_to_w_delta_r_psp(2,0) = 0.;
        T_dangles_to_w_delta_r_psp(2,1) = sin(delta_r_psp ( 2 ));
        T_dangles_to_w_delta_r_psp(2,2) = cos(delta_r_psp ( 2 ))*cos(delta_r_psp ( 1 ));
        
        delta_w_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_vel;
        
//         dw_psp = T_dangles_to_w_delta_r_psp * delta_r_psp_acc;
        
        for (int ii=0; ii<3; ii++)
        {
            dw_psp(ii) = (delta_w_psp(ii) - delta_w_psp_old(ii)) / timeInterval;
            delta_w_psp_old(ii) = delta_w_psp(ii);
        }
        
    }
    
    if (imp_mode_t == 1)
    {
        
        impPosesetpoint.pose.position.x = xDev_psp(0)+start_pos_t(0);
        impPosesetpoint.pose.position.y = xDev_psp(1)+start_pos_t(1);
        impPosesetpoint.pose.position.z = xDev_psp(2)+start_pos_t(2);
        
        frameDev.p.p=xDev_psp+start_pos_t;
        
        frameDev.p.v=xdDev_psp;
    }
    else if (imp_mode_t == 2)
    {
        impPosesetpoint.pose.position.x = xDev_psp(0)+pos_t_msr(0);
        impPosesetpoint.pose.position.y = xDev_psp(1)+pos_t_msr(1);
        impPosesetpoint.pose.position.z = xDev_psp(2)+pos_t_msr(2);
        
        frameDev.p.p=xDev_psp+pos_t_msr;
        
        frameDev.p.v=xdDev_psp+vel_t_msr;
    }
    
    if (imp_mode_r == 1)
    {
        new_r = start_pos_r + phi_cd_cmd;
    }
    else if (imp_mode_r == 2)
    {
        new_r = pos_r_msr + delta_r_psp;
    }
    
    impPosesetpoint.header.stamp=ros::Time::now();
    
    impPoseSp_pub.publish ( impPosesetpoint );

    if ( angleConvection==c_impedance_control::EULER_ZYX ) {
        R_new = KDL::Rotation::EulerZYX( new_r ( 2 ),new_r ( 1 ),new_r ( 0 ) );
    } else if ( angleConvection==c_impedance_control::RPY ) {
        R_new = KDL::Rotation::RPY( new_r ( 0 ),new_r ( 1 ),new_r ( 2 ) );
    } else {
        R_new = KDL::Rotation::EulerZYX( new_r ( 2 ),new_r ( 1 ),new_r ( 0 ) );
    }

    frameDev.M.R = R_new;
    frameDev.M.w=delta_w_psp+vel_r_sp;
    frameDev.M.dw=dw_psp;
    frameDev.p.dv=xddDev_psp;
        
    frameCom=frameDev;

    p_dev = frameDev.p.p;

    R_cmd = frameDev.M.R;

    p_tmp = frameCom.p.p;

    if ( cont2 < -50) {

        cont2 = 0;

//         if ( setfirstEntry_psp == true )
//             std::cout << "true" << std::endl;
//         else
//             std::cout << "false" << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "K: " << Kp_t ( 0,0 ) << "; " << Kp_t ( 1,1 ) << "; " << Kp_t ( 2,2 ) << "; " << std::endl;
//         std::cout << "D: " << Kd_t ( 0,0 ) << "; " << Kd_t ( 1,1 ) << "; " << Kd_t ( 2,2 ) << "; " << std::endl;
        
//         std::cout << "2' method" << std::endl;
//         std::cout << " " << std::endl;
        std::cout << "xDev_psp: " << xDev_psp ( 0 ) << "; " << xDev_psp ( 1 ) << "; " << xDev_psp ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "start_pos_t: " << start_pos_t ( 0 ) << "; " << start_pos_t ( 1 ) << "; " << start_pos_t ( 2 ) << "; " << std::endl;
//         std::cout << "xddDev_psp: " << xddDev_psp ( 0 ) << "; " << xddDev_psp ( 1 ) << "; " << xddDev_psp ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "acc_t: " << acc_t ( 0 ) << "; " << acc_t ( 1 ) << "; " << acc_t ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "Torques: " << Torques ( 0 ) << "; " << Torques ( 1 ) << "; " << Torques ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "pos_t_sp: " << pos_t_sp ( 0 ) << "; " << pos_t_sp ( 1 ) << "; " << pos_t_sp ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "phi_cd_start: " << phi_cd_start ( 0 ) << "; " << phi_cd_start ( 1 ) << "; " << phi_cd_start ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "start_pos_r: " << start_pos_r ( 0 ) << "; " << start_pos_r ( 1 ) << "; " << start_pos_r ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "pos_r_msr: " << pos_r_msr ( 0 ) << "; " << pos_r_msr ( 1 ) << "; " << pos_r_msr ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "pos_r_sp: " << pos_r_sp ( 0 ) << "; " << pos_r_sp ( 1 ) << "; " << pos_r_sp ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "phi_cd_eigen: " << phi_cd_eigen ( 0 ) << "; " << phi_cd_eigen ( 1 ) << "; " << phi_cd_eigen ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "dphi_cd_eigen: " << dphi_cd_eigen ( 0 ) << "; " << dphi_cd_eigen ( 1 ) << "; " << dphi_cd_eigen ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "phi_cd_cmd: " << phi_cd_cmd ( 0 ) << "; " << phi_cd_cmd ( 1 ) << "; " << phi_cd_cmd ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "dphi_cd_cmd: " << dphi_cd_cmd ( 0 ) << "; " << dphi_cd_cmd ( 1 ) << "; " << dphi_cd_cmd ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "ddphi_cd_cmd: " << ddphi_cd_cmd ( 0 ) << "; " << ddphi_cd_cmd ( 1 ) << "; " << ddphi_cd_cmd ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "Force_act: " << Force_act(0)  << "; " << Force_act ( 1 ) << "; " << Force_act ( 2 ) << "; " << Force_act ( 3 )  << "; " << Force_act ( 4 )  << "; " << Force_act ( 5 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "timeInterval: " << timeInterval << std::endl;
        std::cout << " " << std::endl;
        std::cout << "T_plot: " << T_plot ( 0 ) << "; " << T_plot( 1 ) << "; " << T_plot ( 2 ) << "; " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "msrForceTorque_ImpedanceControl: " << msrForceTorque_ImpedanceControl ( 0 ) << "; " << msrForceTorque_ImpedanceControl ( 1 ) << "; " << msrForceTorque_ImpedanceControl ( 2 ) << "; " << msrForceTorque_ImpedanceControl ( 3 ) << "; " << msrForceTorque_ImpedanceControl ( 4 ) << "; " << msrForceTorque_ImpedanceControl ( 5 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "delta_w_psp: " << delta_w_psp ( 0 ) << "; " << delta_w_psp ( 1 ) << "; " << delta_w_psp ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "P_cmd: " << p_tmp ( 0 ) << "; " << p_tmp ( 1 ) << "; " << p_tmp ( 2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "new_r: " << new_r ( 0 ) << "; " << new_r ( 1 ) << "; " << new_r ( 2 ) << "; " << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "delta_r_psp: " << delta_r_psp ( 0 ) << "; " << delta_r_psp ( 1 ) << "; " << delta_r_psp ( 2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "R_cmd: " << R_cmd ( 0,0 ) << "; " << R_cmd ( 0,1 ) << "; " << R_cmd ( 0,2 ) << "; " << R_cmd ( 1,0 ) << "; "
//                   << R_cmd ( 1,1 ) << "; " << R_cmd ( 1,2 ) << "; " << R_cmd ( 2,0 ) << "; " << R_cmd ( 2,1 ) << "; " << R_cmd ( 2,2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "R_msr: " << R_msr ( 0,0 ) << "; " << R_msr ( 0,1 ) << "; " << R_msr ( 0,2 ) << "; " << R_msr ( 1,0 ) << "; "
//                   << R_msr ( 1,1 ) << "; " << R_msr ( 1,2 ) << "; " << R_msr ( 2,0 ) << "; " << R_msr ( 2,1 ) << "; " << R_msr ( 2,2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "R_sp: " << R_sp ( 0,0 ) << "; " << R_sp ( 0,1 ) << "; " << R_sp ( 0,2 ) << "; " << R_sp ( 1,0 ) << "; "
//                   << R_sp ( 1,1 ) << "; " << R_sp ( 1,2 ) << "; " << R_sp ( 2,0 ) << "; " << R_sp ( 2,1 ) << "; " << R_sp ( 2,2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "R_compliant_to_sp_start: " << R_compliant_to_sp_start ( 0,0 ) << "; " << R_compliant_to_sp_start ( 0,1 ) << "; " << R_compliant_to_sp_start ( 0,2 ) << "; " << R_compliant_to_sp_start ( 1,0 ) << "; "
//                   << R_compliant_to_sp_start ( 1,1 ) << "; " << R_compliant_to_sp_start ( 1,2 ) << "; " << R_compliant_to_sp_start ( 2,0 ) << "; " << R_compliant_to_sp_start ( 2,1 ) << "; " << R_compliant_to_sp_start ( 2,2 ) << std::endl;
//         std::cout << " " << std::endl;
//         std::cout << "R_start: " << R_start ( 0,0 ) << "; " << R_start ( 0,1 ) << "; " << R_start ( 0,2 ) << "; " << R_start ( 1,0 ) << "; "
//                   << R_start ( 1,1 ) << "; " << R_start ( 1,2 ) << "; " << R_start ( 2,0 ) << "; " << R_start ( 2,1 ) << "; " << R_start ( 2,2 ) << std::endl;
//         std::cout << " " << std::endl;  
//         std::cout << "T_dangles_to_w_cd: " << T_dangles_to_w_cd(0,0) << "; " << T_dangles_to_w_cd(0,1) << "; " << T_dangles_to_w_cd(0,2) << "; " << T_dangles_to_w_cd(1,0) << "; "
//                   << T_dangles_to_w_cd(1,1) << "; " << T_dangles_to_w_cd(1,2) << "; " << T_dangles_to_w_cd(2,0) << "; " << T_dangles_to_w_cd(2,1) << "; " << T_dangles_to_w_cd(2,2) << std::endl;
//         std::cout << " " << std::endl;  
    std::cout << "-------------------------------" << std::endl;

    }

    return frameCom;

}



}
}