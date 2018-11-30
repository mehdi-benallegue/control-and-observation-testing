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
  bool verbose = false;
  const std::string filename = "/opt/openrobots/share/ur_description/urdf/ur5_gripper.urdf";

  

  se3::Model model;
  se3::urdf::buildModel (filename, model, verbose);
  se3::Data data(model);

  double dt = 0.0001;

  int iterations = int(10. / dt); //10seconds

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
  const std::string filenameReal = "/opt/openrobots/share/ur_description/urdf/ur5_gripper.urdf";

  se3::Model modelReal;
  se3::urdf::buildModel (filenameReal, modelReal, verbose);
  se3::Data dataReal(modelReal);

  /////////////// initialization of the real robot ///////
  IndexedVectorArray q ;
  IndexedVectorArray vq ;
  IndexedVectorArray aq ;
  IndexedVectorArray aq_ref ;

  controlPrototyping::IntegralPassivityBasedControl intglPassCntrl(model.nv);
  intglPassCntrl.setSamplingTime(dt);

  IndexedVectorArray q_error;
  IndexedVectorArray vq_error;

  Vector q_current = Vector::Random(modelReal.nq);
  Vector vq_current = Vector::Random(modelReal.nv);
  Vector aq_current = Vector::Random(modelReal.nv);

  Vector aq0 = Vector::Zero(modelReal.nv);/// just a zero acceleration



  ///////////////////Tracking of the desired trajectory ////////////
  for (int i=0; i< iterations; ++i)
  {

    q.pushBack(q_current);
    vq.pushBack(vq_current);
    aq.pushBack(aq_current);

    q_error.pushBack(q_des[i]-q_current);
    vq_error.pushBack(vq_des[i]-vq_current);

    ////////////////////getting the torque///////////////////////
    Vector tau;

    enum mode
    {
      kinematicFeedback,
      integralPassivity
    } torqueMode;

    torqueMode = integralPassivity;

    if (torqueMode == kinematicFeedback)
    {
      //////////////////////getting the reference acceleration /////////////////
      double gainp = 50 ; 
      double gainv = 2*sqrt(gainp);
      
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
      

      tau = intglPassCntrl.getTorque(tau_des,vq_des[i]-vq_current,q_des[i]-q_current,data.M,C);

      // std::cout << tau.transpose() << std::endl;
    }

    /////////////////////getting real acceleration//////////
    aq_current=se3::forwardDynamics(modelReal, dataReal, q_current, vq_current, 
                                           tau,Matrix::Zero(0,modelReal.nv),Vector::Zero(0));

    q_current += vq_current*dt + aq_current*dt*dt/2;
    vq_current += aq_current*dt;



  }


  



  q.writeInFile("TMP-q");
  vq.writeInFile("TMP-vq");
  aq.writeInFile("TMP-aq");
  q_des.writeInFile("TMP-q_des");
  vq_des.writeInFile("TMP-vq_des");
  aq_des.writeInFile("TMP-aq_des");
  aq_ref.writeInFile("TMP-aq_ref");
  q_error.writeInFile("TMP-q_error");
  vq_error.writeInFile("TMP-vq_error");

  return 0;


  
}