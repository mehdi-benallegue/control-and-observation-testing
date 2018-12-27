#include <iostream>
#include <fstream>

#include <stdlib.h>

#include <yaml-cpp/yaml.h>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/examples/offline-model-base-flex-estimation.hpp>

#include <time.h>


using namespace stateObservation;

int main (int argc, char *argv[])
{
  if (argc!=6 && argc !=2)
  {
    std::cout << "Error: Misused application"<< std::endl;
    std::cout << "Use either by parameters or by configuration file" <<std::endl;
    std::cout << "measurements inputs numberOfContacts dt m"<<std::endl;
    std::cout << "Parameters : "<<std::endl;
    std::cout << "Leaving"<<std::endl;

  }
  else
  {
    IndexedVectorArray y,u,numberOfContacts;


    bool withForces =true;

    IndexedMatrixArray Q, R;

    int stateSize = flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU::staticGetStateSize();

    Matrix3 Kfe,Kfv,Kte,Ktv;

    Kfe.setZero();
    Kfv.setZero();
    Kte.setZero();
    Ktv.setZero();

    
    if (argc!=2)
    {
      std::cout << "Read Sensors filename: "<<argv[1]<< std::endl;
      y.readVectorsFromFile(argv[1]);
      std::cout << "Sensors loaded, size:"<< y.size() << std::endl;
      std::cout << "Read Inputs, filename: "<<argv[2]<< std::endl;
      u.readVectorsFromFile(argv[2]);
      std::cout << "Inputs loaded, size:"<< u.size() << std::endl;

      std::cout << "Read contact number, filename: "<<argv[3]<< std::endl;
      numberOfContacts.readVectorsFromFile(argv[3]);
      std::cout << "Contact numbers loaded, size:"<< numberOfContacts.size() << std::endl;
    }
    else
    {
      std::cout << "Read Configuration file : " << argv[1] << std::endl;
      YAML::Node config = YAML::LoadFile(argv[1]);


      std::cout << "Read Sensors filename: " << config["measurementFile"].as<std::string>() << std::endl;
      y.readVectorsFromFile(config["measurementFile"].as<std::string>());
      std::cout << "Sensors loaded, size:"<< y.size() << std::endl;
      std::cout << "Read Inputs, filename: " << config["inputFile"].as<std::string>() << std::endl;
      u.readVectorsFromFile(config["inputFile"].as<std::string>());
      std::cout << "Inputs loaded, size:"<< u.size() << std::endl;


      bool constantNumberOfContact=true;


      if (config["contactsNumberFile"].IsDefined())
      {

        std::cout << "Read contact number, filename: " << config["contactsNumberFile"].as<std::string>() << std::endl;
        
          
        numberOfContacts.readVectorsFromFile(config["contactsNumberFile"].as<std::string>());
        std::cout << "Contact numbers loaded, size:"<< numberOfContacts.size() << std::endl;

        unsigned i = 1;

        while ( i<numberOfContacts.size() && constantNumberOfContact )
        {
          if (numberOfContacts[i](0) < numberOfContacts[i-1](0))
          {
            constantNumberOfContact=false;
          }
          ++i;
        }
      }
      else
      {
        std::cout << "Read contact number " << config["contactsNumber"].as<std::string>() << std::endl;
        
        for(size_t i = 0; i < u.size(); i++)
        {
          
        }
        
        numberOfContacts.pushBack(Vector1::Constant(config["contactsNumber"].as<int>()) );
      }

      if (config["withForces"].IsDefined())
      {
        withForces = config["withForces"].as<bool>();
      }
      
      if (config["Qfile"].IsDefined())
      {
        std::cout << "Read Q from file " << config["Qfile"] << std::endl;
        Q.readFromFile(config["Qfile"].as<std::string>(),
              stateSize,
              stateSize);
      }
      else if (config["Q"].IsDefined())
      {
        std::cout << "use constant Q" << std::endl;
        Q.pushBack(tools::stringToMatrix(config["Q"].as<std::string>(), stateSize,  stateSize));
      }
      else if (config["Qdiagonal"].IsDefined())
      {
        std::cout << "use diagonal Q" << std::endl;
        Matrix Qmat(stateSize,    stateSize);
        Qmat.setZero();
        Qmat.diagonal() = tools::stringToVector(config["Qdiagonal"].as<std::string>());
        Q.pushBack(Qmat);
      }
      
      if (config["RFile"].IsDefined())
      {
        std::cout << "Read R from file " << config["RFile"].as<std::string>() <<std::endl;
        if (constantNumberOfContact && withForces)
        {
          if (withForces)
          {
            R.readFromFile(config["RFile"].as<std::string>(), 6 + size_t(numberOfContacts[0](0))*6 );
          }
          else
          {
            R.readFromFile(config["RFile"].as<std::string>(), 6);
          }
          
        } 
        else
        {
          throw std::runtime_error(" Impossible to set R matrix from a file with non constant number of contacts and force sensors");
        }
      }
      else if (config["RIMUDiag"].IsDefined())
      {
        std::cout << "Read R IMU diagonal " << std::endl;
        Vector6 RimuDiag  = tools::stringToVector(config["RIMUDiag"].as<std::string>());
        Matrix6 Rimu;
        Rimu.setZero();
        Rimu.diagonal()=RimuDiag;

        Vector6 RftDiag;

        if (withForces)
        {
          std::cout << "Read R Force Torque diagonal " << std::endl;
          RftDiag = tools::stringToVector(config["RFTDiag"].as<std::string>());
        }

        Matrix Ri;
        Ri = Rimu;

        for(size_t i = 0; i < u.size(); i++)
        {
          if (withForces)
          {
            unsigned contacts = unsigned(numberOfContacts[i](0));
            Ri.conservativeResize(6+contacts*6,6+contacts*6);
            Ri.block(0,6,6,contacts*6).setZero();
            Ri.block(6,0,contacts*6,contacts*6+6).setZero();
            for(size_t j = 0; j < contacts; j++)
            {
              Ri.block(6+j,6+j,6,6).diagonal()=RftDiag;
            }            
          }          
        }

        R.pushBack(Ri);
      }

      if (config["Kfe"].IsDefined())
      {
        std::cout << "Reading linear stiffness matrix "<< std::endl;
        Kfe = tools::stringToMatrix(config["Kfe"].as<std::string>(),3,3);
      }
      else if (config["kfeDiag"].IsDefined()) 
      {
        std::cout << "Reading linear stiffness matrix diagonal"<< std::endl;
        Kfe.diagonal() = tools::stringToVector(config["KfeDiag"].as<std::string>());
      }

      if (config["Kfv"].IsDefined())
      {
        std::cout << "Reading linear damping matrix "<< std::endl;
        Kfv = tools::stringToMatrix(config["Kfv"].as<std::string>(),3,3);
      }
      else if (config["kfvDiag"].IsDefined()) 
      {
        std::cout << "Reading linear damping matrix diagonal"<< std::endl;
        Kfv.diagonal() = tools::stringToVector(config["KfvDiag"].as<std::string>());
      }

      if (config["Kte"].IsDefined())
      {
        std::cout << "Reading angular stiffness matrix "<< std::endl;
        Kte = tools::stringToMatrix(config["Kte"].as<std::string>(),3,3);
      }
      else if (config["kteDiag"].IsDefined()) 
      {
        std::cout << "Readinf angular damping matrix diagonal"<< std::endl;
        Kfe.diagonal() = tools::stringToVector(config["KfeDiag"].as<std::string>());
      }

      if (config["Ktv"].IsDefined())
      {
        std::cout << "Reading angular damping matrix "<< std::endl;
        Ktv = tools::stringToMatrix(config["Ktv"].as<std::string>(),3,3);
      }
      else if (config["ktvDiag"].IsDefined()) 
      {
        std::cout << "Reading angular damping matrix diagonal"<< std::endl;
        Kfv.diagonal() = tools::stringToVector(config["KtvDiag"].as<std::string>());
      }
      
    }

    

   // std::cout <<"numberOfContacts " <<numberOfContacts.getLastIndex()<< " "
   //           << numberOfContacts[numberOfContacts.getLastIndex()].size() << " "
   //           << numberOfContacts[numberOfContacts.getLastIndex()] << " "
   //           << std::endl;//<< " "<< numberOfContacts[2]<< " "<< numberOfContacts[3]<< " ";

    Matrix xh0 = Vector::Zero(flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::size,1);

    double dt=atof(argv[4]);
    double mass=atof(argv[5]);

    std::cout << "Rebuiding state" << std::endl;
    IndexedVectorArray xhat=
      examples::offlineModelBaseFlexEstimation( y, u, xh0, numberOfContacts, dt, mass, withForces,
                                               IndexedMatrixArray(), IndexedMatrixArray(),Kfe, Kfv, Kte, Ktv,
                                                0x0,0x0,0x0,0x0,1);
    std::cout << "State rebuilt" << std::endl;

    xhat.writeInFile("xhat");
  }
  return 0;

}





