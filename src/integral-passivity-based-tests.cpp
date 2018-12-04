
#include <yaml-cpp/yaml.h>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/dynamics.hpp>


#include <state-observation/tools/definitions.hpp>

#include <control-prototyping/integral-passivity-based-control.hpp>



int main()
{

  using namespace stateObservation;

  YAML::Node config = YAML::LoadFile("config.yaml");

  bool verbose = false;
  const std::string filename = config["filename_model"].as<std::string>();

  se3::Model model;
  se3::urdf::buildModel (filename, model, verbose);
  se3::Data data(model);

  double dt = config["dt"].as<double>();

  int iterations = int(config["Duration"].as<double>() / dt); //10seconds

  ///////////// CREATION OF THE DESIRED TRAJECTORY ////////////

  Vector q_des_current = Vector::Random(model.nq);
  Vector vq_des_current = Vector::Random(model.nv);
  Vector aq_des_current = Vector::Random(model.nv);
  Vector jq_des_current = Vector::Random(model.nv); ///jerk

  IndexedVectorArray q_des;
  IndexedVectorArray vq_des;
  IndexedVectorArray aq_des;
  

  for (int i=0; i< iterations; ++i)
  {
    q_des.pushBack(q_des_current);
    vq_des.pushBack(vq_des_current);
    aq_des.pushBack(aq_des_current);

    q_des_current += vq_des_current*dt + aq_des_current*dt*dt/2;
    vq_des_current += aq_des_current*dt;
    aq_des_current += jq_des_current*dt;
    jq_des_current = 10*Vector::Random(model.nv)-q_des_current- vq_des_current- aq_des_current;
        
  }


  ////////////////////// Loading the "real" robot for simulation /////////////
  const std::string filenameReal = config["filename_real"].as<std::string>();

  se3::Model modelReal;
  se3::urdf::buildModel (filenameReal, modelReal, verbose);
  se3::Data dataReal(modelReal);

  /////////////// initialization of the real robot ///////
  IndexedVectorArray q ;
  IndexedVectorArray vq ;
  IndexedVectorArray aq ;
  IndexedVectorArray aq_ref ;
  IndexedVectorArray power ;
  IndexedVectorArray energy ;

  controlPrototyping::IntegralPassivityBasedControl intglPassCntrl(model.nv);
  intglPassCntrl.setSamplingTime(dt);

  IndexedVectorArray q_error;
  IndexedVectorArray vq_error;

  Vector q_current = Vector::Random(modelReal.nq);
  Vector vq_current = Vector::Random(modelReal.nv);
  Vector aq_current = Vector::Random(modelReal.nv);

  Vector aq0 = Vector::Zero(modelReal.nv);/// just a zero acceleration

  Vector q_error_current(model.nv+1);
  Vector vq_error_current(model.nv+1);

  Vector power_current(1);
  Vector energy_current(1);
  energy_current.setZero();

 
 
  double mu =  config["Passivity-basedPID_mu"].as<double>();
  double lambda = config["Passivity-basedPID_lambda"].as<double>();
  double gamma =  config["Passivity-basedPID_gamma"].as<double>();
  double L_multiplier =  config["Passivity-basedPID_L"].as<double>();
  double Ka_multiplier = config["Passivity-basedPID_Ka"].as<double>();

  bool L_mass  =  config["Passivity-basedPID_L_Mass"].as<bool>();
  bool Ka_mass  =  config["Passivity-basedPID_Ka_Mass"].as<bool>();

  Matrix L(modelReal.nv,modelReal.nv),Ka(modelReal.nv,modelReal.nv);

  if (!L_mass)
  {
    L.setIdentity();        
    L*=L_multiplier;
  }

  if (!Ka_mass)
  {
    Ka.setIdentity();        
    Ka*=Ka_multiplier;
  }

  enum mode
  {
    kinematicFeedback,
    integralPassivity
  } torqueMode;

  std::string filePrefix;

  switch (config["Torquemode"].as<int>())
  {
  case 0: 
    torqueMode = kinematicFeedback;
    filePrefix = "kf";
    break;
  case 1:
    torqueMode = integralPassivity;
    filePrefix = "ip";
  }

  ///////////////////Tracking of the desired trajectory ////////////
  for (int i=0; i< iterations; ++i)
  {

    q.pushBack(q_current);
    vq.pushBack(vq_current);
    aq.pushBack(aq_current);

    

    q_error_current.head(model.nv)=q_des[i]-q_current;
    vq_error_current.head(model.nv)=vq_des[i]-vq_current;

    q_error_current(model.nv)= q_error_current.head(model.nv).norm();   
    vq_error_current(model.nv)= vq_error_current.head(model.nv).norm();   

    q_error.pushBack(q_error_current);
    vq_error.pushBack(vq_error_current);

    ////////////////////getting the torque///////////////////////
    Vector tau;

    if (torqueMode == kinematicFeedback)
    {
      //////////////////////getting the reference acceleration /////////////////
      double gainp = config["Kinematic_GainP"].as<double>(); 
      double gainv = config["Kinematic_GainV"].as<double>();
      
      Vector aq_ref_current  = aq_des[i] + gainv * (vq_des[i]-vq_current) + gainp * (q_des[i]-q_current);
      aq_ref.pushBack(aq_ref_current);      
      
      computeAllTerms(model, data, q_current, vq_current);
      data.M.triangularView<Eigen::StrictlyLower>() = data.M.transpose().triangularView<Eigen::StrictlyLower>();
      
      tau = data.M * aq_ref_current + data.nle;

    }
    else if (torqueMode == integralPassivity)
    {

      computeAllTerms(model, data, q_current, vq_current);
      data.M.triangularView<Eigen::StrictlyLower>() = data.M.transpose().triangularView<Eigen::StrictlyLower>();

      Matrix C = se3::computeCoriolisMatrix(model,data,q_current,vq_current );
      
      Vector tau_des = data.M * aq_des[i] + data.nle;

      if (L_mass)
      {
        L=L_multiplier*data.M;        
      }

      if (Ka_mass)
      {
        Ka=Ka_multiplier*data.M;        
      }
      


      intglPassCntrl.setGains(Ka,L,lambda,mu,gamma);

      tau = intglPassCntrl.getTorque(tau_des,vq_des[i]-vq_current,q_des[i]-q_current,data.M,C);

     // std::cout << tau.transpose() << std::endl;
    }

    /////////////////////getting real acceleration//////////
    aq_current=se3::forwardDynamics(modelReal, dataReal, q_current, vq_current, 
                                           tau,Matrix::Zero(0,modelReal.nv),Vector::Zero(0));

    power_current(0)=tau.dot(vq_current);
    energy_current(0)+=fabs(power_current(0));
    power.pushBack(power_current);
    energy.pushBack(energy_current);

    q_current += vq_current*dt + aq_current*dt*dt/2;
    vq_current += aq_current*dt;



  }
  q.writeInFile(filePrefix+"-q");
  vq.writeInFile(filePrefix+"-vq");
  aq.writeInFile(filePrefix+"-aq");
  q_des.writeInFile(filePrefix+"-q_des");
  vq_des.writeInFile(filePrefix+"-vq_des");
  aq_des.writeInFile(filePrefix+"-aq_des");
  aq_ref.writeInFile(filePrefix+"-aq_ref");
  q_error.writeInFile(filePrefix+"-q_error");
  vq_error.writeInFile(filePrefix+"-vq_error");
  power.writeInFile(filePrefix+"-power");
  energy.writeInFile(filePrefix+"-energy");

 
  return 0;


  
}