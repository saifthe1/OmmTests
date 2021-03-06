#include<iostream>
#include<fstream>
#include<stdio.h>
#include "OpenMM.h"
#include "rapidjson/document.h"
#include "gtest/gtest.h"

TEST(TestBinProperties,ReadJson){
	std::ifstream file("../foam.json");
	if(!file.is_open()){
		std::cout<<"File not found"<<std::endl;
		FAIL();
	}
	else{
		SUCCEED();
	}
}


TEST(TestBinProperties,ParseJson){
	std::ifstream file("../foam.json");
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

TEST(TestBinProperties,MeasurementTools){
    std::vector<std::string> mtools(1);
    mtools[0] = "MeasureBinProperties";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::MeasurementTools measurements(mtools,startPoint,endPoint,1,1);
    EXPECT_EQ(1,measurements.getNBins());
    EXPECT_EQ(1,measurements.getWriteInterval());
}

TEST(TestBinProperties,SingleBin){
    std::ifstream file("../foam-1728.json");
    if(!file.is_open()){
        std::cout<<"File not found"<<std::endl;
        FAIL();
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
        FAIL();
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
    std::vector<std::string> mtools(1);
    mtools[0] = "MeasureBinProperties";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::MeasurementTools measurements(mtools,startPoint,endPoint,1,1);
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
    OpenMM::Context context(system,integrator,platform,measurements);
    string Platformname = context.getPlatform().getName();
    std::cout<<"Using OpenMM "<<Platformname.c_str()<<" Platform"<<std::endl;
    std::vector<OpenMM::Vec3> temppos;
    context.setPositions(posInNm);
    context.setVelocities(velInNm);
    
    
    int nbs = measurements.getNBins();
    int writeinterval = measurements.getWriteInterval();
    //start the simulation loop
    int counter=0;
    std::cout<<"starting the loop\n";
    for(int frame=1;frame<400;++frame){
        integrator.step(1);
        int* mola = context.getMeasurements().getMols();
        OpenMM::Vec3* mom = context.getMeasurements().getBinMom();
        double gpuKE = 0.0;
        int gpuMol = 0.0;
        for(int k=0;k<nbs;k++){
            gpuMol += mola[k];
        }
        EXPECT_EQ(numParticles,gpuMol);
    }
    std::cout<<"Simulation completed with counter value "<<counter<<std::endl;
}

TEST(TestBinProperties,DISABLED_TwoBins){
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
    std::vector<std::string> mtools(1);
    mtools[0] = "MeasureBinProperties";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::MeasurementTools measurements(mtools,startPoint,endPoint,2,100);
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
    OpenMM::Context context(system,integrator,platform,measurements);
    string Platformname = context.getPlatform().getName();
    std::cout<<"Using OpenMM "<<Platformname.c_str()<<" Platform"<<std::endl;
    std::vector<OpenMM::Vec3> temppos;
    context.setPositions(posInNm);
    context.setVelocities(velInNm);
    
    
    int nbs = measurements.getNBins();
    int writeinterval = measurements.getWriteInterval();
    //start the simulation loop
    integrator.step(1);
    int counter=0;
    std::cout<<"starting the loop\n";
    for(int frame=1;frame<400;++frame){
            int* mola = context.getMeasurements().getMols();
            OpenMM::Vec3* mom = context.getMeasurements().getBinMom();
            double gpuKE = 0.0;
            int gpuMol = 0.0;
            for(int k=0;k<nbs;k++){
                gpuMol += mola[k];
            }
            EXPECT_EQ(numParticles,gpuMol);
        integrator.step(100);
    }
    std::cout<<"Simulation completed with counter value "<<counter<<std::endl;
}

TEST(TestBinProperties,DISABLED_LargeSize){
    std::ifstream file("../foam.json");
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
    std::vector<std::string> mtools(1);
    mtools[0] = "MeasureBinProperties";
    OpenMM::Vec3 startPoint = OpenMM::Vec3((10*0.34),(0*0.34),(10*0.34));
    OpenMM::Vec3 endPoint = OpenMM::Vec3((10*0.34),(20*0.34),(10*0.34));
    OpenMM::MeasurementTools measurements(mtools,startPoint,endPoint,10,100);
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
    OpenMM::Context context(system,integrator,platform,measurements);
    string Platformname = context.getPlatform().getName();
    std::cout<<"Using OpenMM "<<Platformname.c_str()<<" Platform"<<std::endl;
    std::vector<OpenMM::Vec3> temppos;
    context.setPositions(posInNm);
    context.setVelocities(velInNm);
    
    
    int nbs = measurements.getNBins();
    int writeinterval = measurements.getWriteInterval();
    //start the simulation loop
    integrator.step(1);
    int counter=0;
    std::cout<<"starting the loop\n";
    for(int frame=1;frame<400;++frame){
        counter++;
        OpenMM::State state = context.getState(OpenMM::State::Energy);
        const double time = state.getTime();
        if(counter==writeinterval){
            int* mola = context.getMeasurements().getMols();
            OpenMM::Vec3* mom = context.getMeasurements().getBinMom();
            double gpuKE = 0.0;
            int gpuMol = 0.0;
            for(int k=0;k<nbs;k++){
                gpuMol += mola[k];
            }
            EXPECT_EQ(numParticles*counter,gpuMol);
            counter=1;
        }
        integrator.step(1);
    }
    std::cout<<"Simulation completed with counter value "<<counter<<std::endl;
}

int main(int argc, char **argv){
	testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}
