#include<cstdlib>
#include<cstdio>

#include<Common.h>
#include<Sinkhorn.h>

#include <chrono>
#include<iostream>
#include<vector>
int Interpolation();
int main(){
  Interpolation();
  return 0;

}
double  Multiscale_solver(int image_size, std::vector<double>& muXdat, std::vector<double>& muYdat);
int main2(int argc, char** argv) {

    if (argc <= 1)
        return -1;
    int image_size = atoi(argv[1]);
    std::string filex, filey;


    int cnt = 0;

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= 1; ++i) {
        for (int j = 2; j <= 2; j++) {

            cnt++;
            printf("i=%d\n", i);
            filex = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(image_size) + \
                "_100" + std::to_string(i) + ".dat";
            filey = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(image_size) + \
                "_100" + std::to_string(j) + ".dat";
            std::vector<double> muXdat = readFile<double>(filex.c_str());
            std::vector<double> muYdat = readFile<double>(filey.c_str());

            std::chrono::high_resolution_clock::time_point periter = std::chrono::high_resolution_clock::now();
            double res = Multiscale_solver(image_size, muXdat, muYdat);
            std::chrono::high_resolution_clock::time_point periter_endTime = std::chrono::high_resolution_clock::now();

            std::cout << "time: " << (double)std::chrono::duration_cast<std::chrono::microseconds>(periter_endTime - periter).count() / (1000000.) << std::endl;
            std::cout << res << std::endl;
        }
    }

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();

    std::cout << "time: " << (double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / (1000000. * (double)cnt) << std::endl;
    return 0;

}

 

int W2Grid() {

    // this example is the Sinkhorn analogon to the CPLEX W2Grid() example.

    double muXdat[]={6.101275e-03, 7.204110e-03, 8.535062e-03, 9.910280e-03, 1.104128e-02, 1.166015e-02, 1.167389e-02, 1.123175e-02, 6.704433e-03, 8.128397e-03, 9.854252e-03, 1.163492e-02, 1.309322e-02, 1.388739e-02, 1.390795e-02, 1.335771e-02, 7.292185e-03, 9.029484e-03, 1.114963e-02, 1.334184e-02, 1.513976e-02, 1.613150e-02, 1.620053e-02, 1.561280e-02, 7.908921e-03, 9.950044e-03, 1.245024e-02, 1.503625e-02, 1.715816e-02, 1.834562e-02, 1.848901e-02, 1.792387e-02, 8.567936e-03, 1.091527e-02, 1.379087e-02, 1.675992e-02, 1.919501e-02, 2.057699e-02, 2.080940e-02, 2.028992e-02, 9.247266e-03, 1.191367e-02, 1.517757e-02, 1.854434e-02, 2.131180e-02, 2.291030e-02, 2.324255e-02, 2.273922e-02, 9.914861e-03, 1.291900e-02, 1.659968e-02, 2.040437e-02, 2.355167e-02, 2.540296e-02, 2.582486e-02, 2.523923e-02, 1.055807e-02, 1.392742e-02, 1.806860e-02, 2.236741e-02, 2.594533e-02, 2.805923e-02, 2.849642e-02, 2.764292e-02};

    double muYdat[]={3.908070e-03, 4.687770e-03, 5.569132e-03, 6.429051e-03, 7.125903e-03, 7.561235e-03, 7.732784e-03, 7.742170e-03, 5.244827e-03, 5.949688e-03, 6.792079e-03, 7.684573e-03, 8.490201e-03, 9.083672e-03, 9.421825e-03, 9.576329e-03, 7.285414e-03, 7.776036e-03, 8.404373e-03, 9.162126e-03, 9.953290e-03, 1.064454e-02, 1.115002e-02, 1.149372e-02, 1.034929e-02, 1.044027e-02, 1.058721e-02, 1.093244e-02, 1.148517e-02, 1.214090e-02, 1.277317e-02, 1.333126e-02, 1.502886e-02, 1.456200e-02, 1.391201e-02, 1.346945e-02, 1.344750e-02, 1.383385e-02, 1.447324e-02, 1.520699e-02, 2.202719e-02, 2.099936e-02, 1.930211e-02, 1.767088e-02, 1.664565e-02, 1.640100e-02, 1.678885e-02, 1.752175e-02, 3.157802e-02, 3.026055e-02, 2.748538e-02, 2.437437e-02, 2.190360e-02, 2.055179e-02, 2.025589e-02, 2.062372e-02, 4.272636e-02, 4.166590e-02, 3.810900e-02, 3.349480e-02, 2.927364e-02, 2.633617e-02, 2.482542e-02, 2.433618e-02};

    // dimensions of grid
    int muXdim[]={8, 8};
    int muYdim[]={8, 8};


    // put raw marginal data into small container structs
    TDoubleMatrix muX,muY;

    muX.data=muXdat;
    muX.dimensions=muXdim;
    muX.depth=2;

    muY.data=muYdat;
    muY.dimensions=muYdim;
    muY.depth=2;


    // fundamental parameters
    int depth=3; // hierarchy depth
    int layerFinest=depth; // what is the finest layer we want to solve?
    int layerCoarsest=1; // coarsest layer to solve on. sometimes skip a few layers at the top

    int msg; // store return codes from functions

    ///////////////////////////////////////////////
    // problem setup
    TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }



    TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    /////////////////////////////////////////////////////////////////////////////////////////
    // setup a cost function provider
    THierarchicalCostFunctionProvider_SquaredEuclidean costProvider(
        MultiScaleSetupX.posH, MultiScaleSetupY.posH,
        MultiScaleSetupX.radii, MultiScaleSetupY.radii,
        muX.depth, layerFinest,
        true,
        MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
        1.);

    /////////////////////////////////////////////////////////////////////////////////////////
    // epsScaling
    double epsStart=1E2;
    double epsTarget=1E-1;
    int epsSteps=10;

    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalesFromBox(muXdim[0],2,MultiScaleSetupX.nLayers); // eps scales for each layer
    msg=epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists for each layer
        // according to scales in line above
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    TSinkhornSolverBase::TSinkhornSolverParameters cfg={
        1E-4, // maxError
        100000, // maxIterations
        100, // innerIterations
        20, // maxAbsorptionLoops
        1E2, // absorption_scalingBound
        1E2, // absorption_scalingLowerBound
        1E-10, // truncation_thresh
        true // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };

    /////////////////////////////////////////////////////////////////////////////////////////
    // create solver object
    TSinkhornSolverStandard SinkhornSolver(
            MultiScaleSetupX.nLayers,
            epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
            layerCoarsest, layerFinest,
            cfg,
            MultiScaleSetupX.HP, MultiScaleSetupY.HP,
            MultiScaleSetupX.muH, MultiScaleSetupY.muH,
            MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
            &costProvider
            );


    SinkhornSolver.initialize();
    msg=SinkhornSolver.solve();
    printf("return code: %d\n",msg);
    if(msg!=0) { return msg; }

    // recompute kernel one last time
    SinkhornSolver.generateKernel();


    // unregularized transport cost
    double transportCost=SinkhornSolver.scorePrimalUnreg();
    printf("unregularized transport cost: %e\n",transportCost);



    // output coupling data
    TSparsePosContainer couplingData=SinkhornKernelGetPosData(SinkhornSolver.kernel);
    // hint: this contains all non-zero elements of truncated kernel,
    // so it will contain a large number of very small entries
    printf("==========================\n");
    printf("coupling data:\n");
    for(int i=0;i<couplingData.nParticles;i++) {
        printf("%d\t%d\t%e\n",couplingData.posStart[i],couplingData.posEnd[i],couplingData.mass[i]);
    }
    printf("==========================\n");



    return 0;
}

using namespace std;
int Interpolation() {
    // computes optimal transport
    // then generates a snapshot of the displacement interpolation
    // and raterizes this onto a grid

    // marginals are simply two Dirac masses each, located on a 10x10 grid
    char filenameMuX[]="../data/muX.dat";
    char filenameMuY[]="../data/muY.dat";
    char xdat[]="../data/Xdat.dat";
    char ydat[]="../data/Ydat.dat";
    std::vector<double> muXdat=readFile<double>(filenameMuX);
    std::vector<double> muYdat=readFile<double>(filenameMuY);
    std::vector<double> posXdat=readFile<double>(xdat);
    std::vector<double> posYdat=readFile<double>(ydat);    
    // masses
    
    //vector<double> muXdat = {0.5, 0.5};
    //vector<double> muYdat = {0.5, 0.5};
 
    // positions (two consecutive entries are xy-coordinates of one Dirac mass)
    //vector<double> posXdat ={1., 2., 1, 8., 2., 3};
    //vector<double> posYdat ={2., 8., 5, 5., 7.,8 };

    // first nr: nr of points, second nr: nr of spatial dimensions
    vector<int> posXdim = {2, 3};
    posXdim[0] = posXdat.size();
    
    vector<int> posYdim = {2, 3};
    posYdim[0] = posYdat.size();double min = 1e100;
    for(int i =0;i < muYdat.size();++i){
    if(min>muYdat[i]){
     std::cout << muYdat[i] << "  ";
     min=muYdat[i];
     }
    }
    //std::cout << std::endl;
    //std::cout << muXdat.size();
    //return 0;

    // plug marginal positions into simple matrix container structs
    TDoubleMatrix posX,posY;

    posX.data=posXdat.data();
    posX.dimensions=posXdim.data();
    posX.depth=3;

    posY.data=posYdat.data();
    posY.dimensions=posYdim.data();
    posY.depth=3;

    // fundamental parameters
    int depth=0; // hierarchy depth
        // (not really needed here: just one trivial layer with full data)
    int layerFinest=depth; // what is the finest layer we want to solve?
    int layerCoarsest=0; // coarsest layer to solve on. sometimes skip a few layers at the top
    int dim=posX.dimensions[1]; // spatial dimension

    int msg; // store return codes from functions


    /////////////////////////////////////////////////////////////////////////////////////////
    // basic hierarchical setup
    TMultiScaleSetupSingleBase MultiScaleSetupX(&posX,muXdat.data(),depth);
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    TMultiScaleSetupSingleBase MultiScaleSetupY(&posY,muYdat.data(),depth);
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    /////////////////////////////////////////////////////////////////////////////////////////
    // setup a cost function provider
    THierarchicalCostFunctionProvider_SquaredEuclidean costProvider(
        MultiScaleSetupX.posH, MultiScaleSetupY.posH,
        MultiScaleSetupX.radii, MultiScaleSetupY.radii,
        dim, layerFinest,
        true,
        MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
        1.);


    /////////////////////////////////////////////////////////////////////////////////////////
    // epsScaling
    double epsStart=1E2;
    double epsTarget=1E-1;
    int epsSteps=15;

    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalingAllFinest(layerFinest);


    TSinkhornSolverBase::TSinkhornSolverParameters cfg={
        1E-4, // maxError
        100000, // maxIterations
        100, // innerIterations
        20, // maxAbsorptionLoops
        1E2, // absorption_scalingBound
        1E2, // absorption_scalingLowerBound
        1E-10, // truncation_thresh
        true // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };

    /////////////////////////////////////////////////////////////////////////////////////////
    // create solver object
    TSinkhornSolverStandard SinkhornSolver(MultiScaleSetupX.nLayers,
            epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
            layerCoarsest, layerFinest,
            cfg,
            MultiScaleSetupX.HP, MultiScaleSetupY.HP,
            MultiScaleSetupX.muH, MultiScaleSetupY.muH,
            MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
            &costProvider
            );


    /////////////////////////////////////////////////////////////////////////////////////////
    // initialize
    SinkhornSolver.initialize();
    // solve
    msg=SinkhornSolver.solve();
    printf("return code: %d\n",msg);
    if(msg!=0) { return msg; }

    // recompute kernel one last time
    SinkhornSolver.generateKernel();


    // unregularized transport cost
    double transportCost=SinkhornSolver.scorePrimalUnreg();
    printf("unregularized transport cost: %e\n",transportCost);

    // extract particles for interpolation
    TSparsePosContainer couplingData=SinkhornKernelGetPosData(SinkhornSolver.kernel);

    // output particle data
    printf("==========================\n");
    printf("coupling data:\n");
    for(int i=0;i<couplingData.nParticles;i++) {
        printf("%d\t%d\t%e\n",couplingData.posStart[i],couplingData.posEnd[i],couplingData.mass[i]);
    }
    printf("==========================\n");
    return 0;

    // generate interpolation

    // characterize base space geometry
    TGeometry_Euclidean geometry;

    int nTime=5;
    for(int tDisc=0;tDisc<=nTime;tDisc++) {
        double t=(double) tDisc/nTime; // at what time do we want to generate it?
        printf("t=%f\n",t);
        TParticleContainer particles=ModelOT_Interpolate<TGeometry_Euclidean>(couplingData,
                posX.data, posY.data,
                dim, t, geometry);

//		printf("==========================\n");
//		printf("particle data:\n");
//		for(int i=0;i<particles.nParticles;i++) {
//			printf("\t%e\t%e\t%e\n",particles.mass[i],particles.pos[2*i],particles.pos[2*i+1]);
//		}
//		printf("==========================\n");

        // rasterize interpolation
        double img[10*10]; // allocate 10x10 pixel grid which will hold rasterized image
        int res[2]={10,10};
        // set to zero
        for(int x=0;x<res[0];x++) {
            for(int y=0;y<res[1];y++) {
                img[res[1]*x+y]=0;
            }
        }
        // bi-linear rasterization of particle data to 2d Cartesian grid
        ProjectInterpolation<2>(particles.pos.data(), particles.mass.data(),
            img, particles.nParticles, res);

        // print rasterization
        printf("==========================\n");
        for(int x=0;x<res[0];x++) {
            for(int y=0;y<res[1];y++) {
                printf("%f  ",img[res[1]*x+y]);
            }
            printf("\n");
        }
        printf("==========================\n");
    }

    return 0;




}

int WpSphere() {
    // this example is the Sinkhorn analogon to the CPLEX WpSphere() example.


    int msg;

    // read raw data from file
    char filenamePoints[]="data/sphere/3081/points.dat";
    char filenameMuX[]="data/sphere/3081/m0.dat";
    char filenameMuY[]="data/sphere/3081/m1.dat";

    std::vector<double> posDat=readFile<double>(filenamePoints);
    std::vector<double> muXdat=readFile<double>(filenameMuX);
    std::vector<double> muYdat=readFile<double>(filenameMuY);


    // fundamental parameters
    int dim=3;
    int res=muXdat.size();
    int depth=4;
    int layerCoarsest=2;
    int layerFinest=depth;
    double p=2;

    ///////////////////////////////////////////////

    // setup problem data
    int posXdim[]={res, dim};

    TDoubleMatrix posX;

    posX.data=posDat.data();
    posX.dimensions=posXdim;
    posX.depth=2;


    ///////////////////////////////////////////////

    TMultiScaleSetupSingleSphere MultiScaleSetupX(&posX,muXdat.data(),depth);
    MultiScaleSetupX.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    // project points back to sphere
    MultiScaleSetupX.SetupProjectPoints();
    // compute sphere radii
    msg=MultiScaleSetupX.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    // allocate hierarchical dual variables
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }

    TMultiScaleSetupSingleSphere MultiScaleSetupY(&posX,muYdat.data(),depth);
    MultiScaleSetupY.HierarchyBuilderChildMode=THierarchyBuilder::CM_Tree;
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    // project points back to sphere
    MultiScaleSetupY.SetupProjectPoints();
    // compute sphere radii
    msg=MultiScaleSetupY.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    // allocate hierarchical dual variables
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    printf("hierarchical cardinalities:\n");
    for(int layer=0;layer<MultiScaleSetupX.nLayers;layer++) {
        printf("%d\t%d\n",layer,MultiScaleSetupX.HP->layers[layer]->nCells);
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup a cost function provider
    THierarchicalCostFunctionProvider_Sphere costProvider(
        MultiScaleSetupX.posH, MultiScaleSetupY.posH,
        MultiScaleSetupX.radii, MultiScaleSetupY.radii,
        dim, 0,
        true,
        MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
        p);


    /////////////////////////////////////////////////////////////////////////////////////////
    // epsScaling
    double epsStart=1E1;
    double epsTarget=.2E-4;
    int epsSteps=25;
    double epsBoxScale=.7;

    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalesFromBox(epsBoxScale,p,MultiScaleSetupX.nLayers); // eps scales for each layer
    msg=epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    /////////////////////////////////////////////////////////////////////////////////////////
    // other parameters
    TSinkhornSolverBase::TSinkhornSolverParameters cfg={
        1E-4, // maxError
        100000, // maxIterations
        100, // innerIterations
        20, // maxAbsorptionLoops
        1E2, // absorption_scalingBound
        1E2, // absorption_scalingLowerBound
        1E-10, // truncation_thresh
        true // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };

    /////////////////////////////////////////////////////////////////////////////////////////
    // create solver object
    TSinkhornSolverStandard SinkhornSolver(MultiScaleSetupX.nLayers, epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
            layerCoarsest, layerFinest,
            cfg,
            MultiScaleSetupX.HP, MultiScaleSetupY.HP,
            MultiScaleSetupX.muH, MultiScaleSetupY.muH,
            MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
            &costProvider
            );


    SinkhornSolver.initialize();
    msg=SinkhornSolver.solve();
    printf("return code: %d\n",msg);
    if(msg!=0) { return msg; }

    // recompute kernel one last time
    SinkhornSolver.generateKernel();
    // return primal objective value
    double primalScore=SinkhornSolver.scorePrimalUnreg();
    printf("primal score: %e\n",primalScore);


    return 0;
}

int WFR() {
    // Wasserstein--Fisher--Rao distance

    double muXdat[]={6.101275e-03, 7.204110e-03, 8.535062e-03, 9.910280e-03, 1.104128e-02, 1.166015e-02, 1.167389e-02, 1.123175e-02, 6.704433e-03, 8.128397e-03, 9.854252e-03, 1.163492e-02, 1.309322e-02, 1.388739e-02, 1.390795e-02, 1.335771e-02, 7.292185e-03, 9.029484e-03, 1.114963e-02, 1.334184e-02, 1.513976e-02, 1.613150e-02, 1.620053e-02, 1.561280e-02, 7.908921e-03, 9.950044e-03, 1.245024e-02, 1.503625e-02, 1.715816e-02, 1.834562e-02, 1.848901e-02, 1.792387e-02, 8.567936e-03, 1.091527e-02, 1.379087e-02, 1.675992e-02, 1.919501e-02, 2.057699e-02, 2.080940e-02, 2.028992e-02, 9.247266e-03, 1.191367e-02, 1.517757e-02, 1.854434e-02, 2.131180e-02, 2.291030e-02, 2.324255e-02, 2.273922e-02, 9.914861e-03, 1.291900e-02, 1.659968e-02, 2.040437e-02, 2.355167e-02, 2.540296e-02, 2.582486e-02, 2.523923e-02, 1.055807e-02, 1.392742e-02, 1.806860e-02, 2.236741e-02, 2.594533e-02, 2.805923e-02, 2.849642e-02, 2.764292e-02};

    double muYdat[]={3.908070e-03, 4.687770e-03, 5.569132e-03, 6.429051e-03, 7.125903e-03, 7.561235e-03, 7.732784e-03, 7.742170e-03, 5.244827e-03, 5.949688e-03, 6.792079e-03, 7.684573e-03, 8.490201e-03, 9.083672e-03, 9.421825e-03, 9.576329e-03, 7.285414e-03, 7.776036e-03, 8.404373e-03, 9.162126e-03, 9.953290e-03, 1.064454e-02, 1.115002e-02, 1.149372e-02, 1.034929e-02, 1.044027e-02, 1.058721e-02, 1.093244e-02, 1.148517e-02, 1.214090e-02, 1.277317e-02, 1.333126e-02, 1.502886e-02, 1.456200e-02, 1.391201e-02, 1.346945e-02, 1.344750e-02, 1.383385e-02, 1.447324e-02, 1.520699e-02, 2.202719e-02, 2.099936e-02, 1.930211e-02, 1.767088e-02, 1.664565e-02, 1.640100e-02, 1.678885e-02, 1.752175e-02, 3.157802e-02, 3.026055e-02, 2.748538e-02, 2.437437e-02, 2.190360e-02, 2.055179e-02, 2.025589e-02, 2.062372e-02, 4.272636e-02, 4.166590e-02, 3.810900e-02, 3.349480e-02, 2.927364e-02, 2.633617e-02, 2.482542e-02, 2.433618e-02};

    // dimensions of grid
    int muXdim[]={8, 8};
    int muYdim[]={8, 8};


    // put raw marginal data into small container structs
    TDoubleMatrix muX,muY;

    muX.data=muXdat;
    muX.dimensions=muXdim;
    muX.depth=2;

    muY.data=muYdat;
    muY.dimensions=muYdim;
    muY.depth=2;


    // fundamental parameters
    int depth=3; // hierarchy depth
    int layerFinest=depth; // what is the finest layer we want to solve?
    int layerCoarsest=1; // coarsest layer to solve on. sometimes skip a few layers at the top
    int dim=2;

    // parameters for WFR distance
    double WFRlenScale=10;
    double WFRKLweight=4*std::pow(WFRlenScale,2);

    int msg; // store return codes from functions


    ///////////////////////////////////////////////
    // problem setup
    TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }

    TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup a cost function provider
    THierarchicalCostFunctionProvider_SquaredEuclidean costProvider(
        MultiScaleSetupX.posH, MultiScaleSetupY.posH,
        MultiScaleSetupX.radii, MultiScaleSetupY.radii,
        dim, layerFinest,
        true,
        MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
        1.,
        true, WFRlenScale);


    /////////////////////////////////////////////////////////////////////////////////////////
    // epsScaling
    double epsStart=1E2;
    double epsTarget=1E-1;
    int epsSteps=10;

    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalesFromBox(muXdim[0],2,MultiScaleSetupX.nLayers); // eps scales for each layer
    msg=epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists for each layer
        // according to scales in line above
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    TSinkhornSolverBase::TSinkhornSolverParameters cfg={
        1E-4, // maxError
        100000, // maxIterations
        100, // innerIterations
        20, // maxAbsorptionLoops
        1E2, // absorption_scalingBound
        1E2, // absorption_scalingLowerBound
        1E-10, // truncation_thresh
        true // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };

    /////////////////////////////////////////////////////////////////////////////////////////
    // create solver object
    TSinkhornSolverKLMarginals SinkhornSolver(MultiScaleSetupX.nLayers,
            epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
            layerCoarsest, layerFinest,
            cfg,
            MultiScaleSetupX.HP, MultiScaleSetupY.HP,
            MultiScaleSetupX.muH, MultiScaleSetupY.muH,
            MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
            &costProvider,
            WFRKLweight
            );


    SinkhornSolver.initialize();
    msg=SinkhornSolver.solve();
    printf("return code: %d\n",msg);
    if(msg!=0) { return msg; }

    // recompute kernel one last time
    SinkhornSolver.generateKernel();


    // unregularized transport cost
    double transportCost=SinkhornSolver.scorePrimalUnreg();
    printf("len scale: %e\n",WFRlenScale);
    printf("unregularized unbalanced transport cost: %e\n",transportCost);


    return 0;
}


double Multiscale_solver(int image_size, std::vector<double>& muXdat, std::vector<double>& muYdat) {

 
    int msg;



    // dimensions of grid
    std::vector<int> muXdim = { image_size, image_size };
    std::vector<int> muYdim = { image_size, image_size };

    // put raw marginal data into small container structs
    TDoubleMatrix muX,muY;

    muX.data=muXdat.data();
    muX.dimensions=muXdim.data();
    muX.depth=2;

    muY.data=muYdat.data();
    muY.dimensions=muYdim.data();
    muY.depth=2;


    // fundamental parameters
    int depth=4; // hierarchy depth
    int layerFinest=depth; // what is the finest layer we want to solve?
    int layerCoarsest=1; // coarsest layer to solve on. sometimes skip a few layers at the top

 

    ///////////////////////////////////////////////
    // problem setup
    TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,depth);
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }


    TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,depth);
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }

    /////////////////////////////////////////////////////////////////////////////////////////
    // setup a cost function provider
    THierarchicalCostFunctionProvider_SquaredEuclidean costProvider(
        MultiScaleSetupX.posH, MultiScaleSetupY.posH,
        MultiScaleSetupX.radii, MultiScaleSetupY.radii,
        muX.depth, layerFinest,
        true,
        MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
        0.5);
 


    /////////////////////////////////////////////////////////////////////////////////////////
    // epsScaling
    double epsStart=1E5;
    double epsTarget=1E-2;
    int epsSteps=20;

    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalesFromBox(muXdim[0],2,MultiScaleSetupX.nLayers); // eps scales for each layer
    msg=epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists for each layer
        // according to scales in line above
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
 

    TSinkhornSolverBase::TSinkhornSolverParameters cfg={
        1E-4, // maxError
        100000, // maxIterations
        100, // innerIterations
        20, // maxAbsorptionLoops
        1E2, // absorption_scalingBound
        1E2, // absorption_scalingLowerBound
        1E-10, // truncation_thresh
        true // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };

    /////////////////////////////////////////////////////////////////////////////////////////
    // create solver object
    TSinkhornSolverStandard SinkhornSolver(
            MultiScaleSetupX.nLayers,
            epsScalingHandler.nEpsLists, epsScalingHandler.epsLists,
            layerCoarsest, layerFinest,
            cfg,
            MultiScaleSetupX.HP, MultiScaleSetupY.HP,
            MultiScaleSetupX.muH, MultiScaleSetupY.muH,
            MultiScaleSetupX.alphaH, MultiScaleSetupY.alphaH,
            &costProvider
            );


    SinkhornSolver.initialize();
    msg=SinkhornSolver.solve();
    printf("return code: %d\n",msg);
    if(msg!=0) { return msg; }

    // recompute kernel one last time
    SinkhornSolver.generateKernel();


    // unregularized transport cost
    double transportCost=SinkhornSolver.scorePrimalUnreg();




    // output coupling data
//	TSparsePosContainer couplingData=SinkhornKernelGetPosData(SinkhornSolver.kernel);
    // hint: this contains all non-zero elements of truncated kernel,
    // so it will contain a large number of very small entries
//	printf("==========================\n");
//	printf("coupling data:\n");
//	for(int i=0;i<couplingData.nParticles;i++) {
//		printf("%d\t%d\t%e\n",couplingData.posStart[i],couplingData.posEnd[i],couplingData.mass[i]);
//	}
//	printf("==========================\n");

        printf("unregularized transport cost: %e\n",transportCost);

    return 0;
}

