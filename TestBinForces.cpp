#include <iostream>
#include <fstream>
#include <stdio.h>
#include "OpenMM.h"
#include "rapidjson/document.h"
#include "gtest/gtest.h"


TEST(TestBinForce,DISABLED_CheckControlTools)
{
    std::vector<std::string> ctools(1);
    ctools[0] = "MeasureBinProperties";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::Vec3 binForces[2];
    binForces[0] = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    binForces[1] = OpenMM::Vec3(0.0f,0.0f,0.0f);
    OpenMM::ControlTools& controls = *new OpenMM::ControlTools(ctools,(double)292,0.001,startPoint,endPoint,2,1);
    controls.setBinForces(binForces);
    OpenMM::Vec3* temp = controls.getBinForces();
//    printf("Bin forces %f, %f, %f\n",temp[0][0],temp[0][1],temp[0][2]);
//    printf("Bin forces %f, %f, %f\n",temp[1][0],temp[1][1],temp[1][2]);
    for(int j=0;j<2;j++){
        EXPECT_EQ(binForces[j],temp[j]);
    }
}

TEST(TestBinForce,DISABLED_CheckKernel)
{
    std::vector<std::string> ctools(1);
    ctools[0] = "ControlBinForces";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::Vec3* forces;
    forces = new OpenMM::Vec3[1];
    forces[0] = OpenMM::Vec3(0.0f,0.0f,0.0f);
    OpenMM::ControlTools tools(ctools,(double)292,0.001,startPoint,endPoint);
    tools.setBinForces(forces);
    OpenMM::Vec3* temp1 = tools.getBinForces();
    OpenMM::Platform::loadPluginsFromDirectory(
           OpenMM::Platform::getDefaultPluginsDirectory());
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    HarmonicBondForce* forceField = new HarmonicBondForce();
	forceField->addBond(0, 1, 1.5, 1);
	system.addForce(forceField);
    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("OpenCL");
    OpenMM::VelocityVerletIntegrator integrator(0.01);
    OpenMM::Context context(system,integrator,platform,tools);
    string Platformname = context.getPlatform().getName();
    ASSERT_STREQ("OpenCL",Platformname.c_str());
    OpenMM::Vec3* temp2 = context.getControls().getBinForces();
  	
    EXPECT_EQ(forces[0],temp2[0]);
    EXPECT_EQ(1,context.getControls().getNBins());
}

TEST(TestBinForce,DISABLED_CheckKernelWithNBins)
{
    std::vector<std::string> ctools(1);
    ctools[0] = "ControlBinForces";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::Vec3* forces;
    forces = new OpenMM::Vec3[10];
    
    for(int j=0;j<10;j++)
        forces[j] = OpenMM::Vec3(0.0f,0.0f,0.0f);
    OpenMM::ControlTools tools(ctools,(double)292,0.001,startPoint,endPoint,0.1,10,1);
    tools.setBinForces(forces);
    
    OpenMM::Platform::loadPluginsFromDirectory(
           OpenMM::Platform::getDefaultPluginsDirectory());
    
    System system;
    system.addParticle(2.0);
    system.addParticle(2.0);
    HarmonicBondForce* forceField = new HarmonicBondForce();
	forceField->addBond(0, 1, 1.5, 1);
	system.addForce(forceField);
    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("OpenCL");
    OpenMM::VelocityVerletIntegrator integrator(0.01);
    OpenMM::Context context(system,integrator,platform,tools);
    string Platformname = context.getPlatform().getName();
    ASSERT_STREQ("OpenCL",Platformname.c_str());
    OpenMM::Vec3* temp = context.getControls().getBinForces();
    
    for(int j=0;j<10;j++)
        EXPECT_EQ(forces[j],temp[j]);
    EXPECT_EQ(10,context.getControls().getNBins());
    
}

TEST(TestBinForce,ParseJson){
	std::ifstream file("../foam-1728.json");
	if(!file.is_open()){
		std::cout<<"File not found"<<std::endl;
		FAIL();
	}
	std::string jsonstr;
	std::string line;
	while (!file.eof()) {
		getline(file, line);
		jsonstr += line;
	}
	file.close();
    
	rapidjson::Document document;
	if(document.Parse<0>(jsonstr.c_str()).HasParseError()){
		FAIL();
	}
	else
	{
		SUCCEED();
	}
}

TEST(TestBinForce,SingleBinTest)
{
    std::ifstream file("../foam-1728.json");
    if(!file.is_open()){
        std::cout<<"File not found"<<std::endl;
        exit(0);
    }
    
    std::string jsonstr;
    std::string line;
    while (!file.eof()) {
        getline(file, line);
        jsonstr += line;
    }
    file.close();
    
    rapidjson::Document document;
    if(document.Parse<0>(jsonstr.c_str()).HasParseError()){
        std::cout<<"Parsing error "<<std::endl;
        exit(0);
    }
    
    const double stepSizeInFs = document["StepSizeInFs"].GetDouble();
    const int numParticles = document["NumberAtoms"].GetInt();
    std::string equationStr = document["Equation"].GetString();
    const double rCut = document["rCut"].GetDouble();
    const rapidjson::Value& b = document["Boxsize"];
    double bx = b[(rapidjson::SizeType)0].GetDouble();
    double by = b[(rapidjson::SizeType)1].GetDouble();
    double bz = b[(rapidjson::SizeType)2].GetDouble();
    int numSpecies = document["NumberSpecies"].GetInt();
    std::vector<OpenMM::Vec3> posInNm;
    std::vector<OpenMM::Vec3> velInNm;
    
    OpenMM::Platform::loadPluginsFromDirectory(
                                               OpenMM::Platform::getDefaultPluginsDirectory());
    
    OpenMM::System system;
    //add particles to the system by reading from json file
    system.setDefaultPeriodicBoxVectors(OpenMM::Vec3(bx,0,0),OpenMM::Vec3(0,by,0),OpenMM::Vec3(0,0,bz));
    const rapidjson::Value& mass = document["masses"];
    const rapidjson::Value& pos = document["Positions"];
    const rapidjson::Value& vel = document["Velocities"];
    
    posInNm.clear();
    velInNm.clear();
    for(rapidjson::SizeType m=0;m<mass.Size();m++){
        double tm = mass[(rapidjson::SizeType) m].GetDouble();
        system.addParticle(tm);
        posInNm.push_back(OpenMM::Vec3(pos[(rapidjson::SizeType) m]["0"].GetDouble(),
                                       pos[(rapidjson::SizeType) m]["1"].GetDouble(),
                                       pos[(rapidjson::SizeType) m]["2"].GetDouble()
                                       ));
        velInNm.push_back(OpenMM::Vec3(vel[(rapidjson::SizeType) m]["0"].GetDouble(),
                                       vel[(rapidjson::SizeType) m]["1"].GetDouble(),
                                       vel[(rapidjson::SizeType) m]["2"].GetDouble()
                                       ));
    }
    std::cout<<"Initialized System with "<<system.getNumParticles()<<" Particles."<<std::endl;
    
    OpenMM::CustomNonbondedForce* nonbonded = new OpenMM::CustomNonbondedForce(equationStr);
    nonbonded->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(rCut);// * OpenMM::NmPerAngstrom
    nonbonded->addPerParticleParameter("Ar");//later on collect it from json string based on OF
    for(int p=0;p<numParticles;p++){
        std::vector<double> params(numSpecies);
        params[0] = (double) 1;
        nonbonded->addParticle(params);
    }
    system.addForce(nonbonded);
    std::vector<std::string> ctools(1);
    ctools[0] = "ControlBinForces";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::Vec3* forces;
    forces = new OpenMM::Vec3[1];
    
    for(int j=0;j<1;j++)
        forces[j] = OpenMM::Vec3(0.0f,0.0f,0.0f);
    OpenMM::ControlTools tools(ctools,(double)292,0.001,startPoint,endPoint);
    tools.setBinForces(forces);
    int num_plat = OpenMM::Platform::getNumPlatforms();
    std::cout<<"Number of registered platforms: "<< num_plat <<std::endl;
    for(int i=0;i<num_plat;i++)
    {
	    OpenMM::Platform& tempPlatform = OpenMM::Platform::getPlatform(i);
	    std::string tempPlatformName = tempPlatform.getName();
	    std::cout<<"Platform "<<(i+1)<<" is "<<tempPlatformName.c_str()<<std::endl;
    }
    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("OpenCL");
    OpenMM::VelocityVerletIntegrator integrator(stepSizeInFs * OpenMM::PsPerFs);
    OpenMM::Context context(system,integrator,platform,tools);
    string Platformname = context.getPlatform().getName();
    ASSERT_STREQ("OpenCL",Platformname.c_str());
    
    context.setPositions(posInNm);
    context.setVelocities(velInNm);
    
    int counter=0;
    std::cout<<"starting the loop\n";
    for(int frame=1;frame<400;++frame){
        counter++;
        integrator.step(1);
    }
    std::cout<<"Simulation completed with counter value "<<counter<<std::endl;
}
int main(int argc, char **argv){
	testing::InitGoogleTest(&argc,argv); 
	return RUN_ALL_TESTS();
}
