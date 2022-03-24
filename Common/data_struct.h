#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H
#include"basic_settings.h"

struct CuSinkhornSolverParameters {
        double maxError; // error accuracy to be achieved
        int maxIterations; // maximal iterations per singleSolve
        int innerIterations; // nr of iterations between absorption and error checks
        int maxAbsorptionLoops; // nr of maximal loops of absorption until algorithm must stabilize
        double absorption_scalingBound; // maximal value of a scaling variable
        double absorption_scalingLowerBound; // value above which we do a safety absorption // NOT YET IMPLEMENTED!
        double truncation_thresh; // truncation threshold in kernel sparsification
        bool refineKernel; // whether initial kernel should be generated via refinement on subsequent layers
};
struct CuCost{
        CuCost(double weight_,
               double WFlenscale_,double WFprefactor_,
               int posDim_,int layerBottom_):
        weight(weight_),
        WFlenscale(WFlenscale_),
        WFprefactor(WFprefactor_),
          posDim(posDim_),
        layerBottom(layerBottom_)
        {}
        double weight; // global rescaling parameter for Euclidean distance
        double WFlenscale; // max transport distance in WF mode
        double WFprefactor; // weight prefactor to avoid repeated computation
        int posDim; // dimensionality of coordinates in xPos and yPos arrays
        int layerBottom; // number of finest layer (0 is coarsest)

};
struct TSinkhornSolverParameters {
        double maxError; // error accuracy to be achieved
        int maxIterations; // maximal iterations per singleSolve
        int innerIterations; // nr of iterations between absorption and error checks
        int maxAbsorptionLoops; // nr of maximal loops of absorption until algorithm must stabilize
        double absorption_scalingBound; // maximal value of a scaling variable
        double absorption_scalingLowerBound; // value above which we do a safety absorption // NOT YET IMPLEMENTED!
        double truncation_thresh; // truncation threshold in kernel sparsification
        bool refineKernel; // whether initial kernel should be generated via refinement on subsequent layers
};
struct CuMultiscal {
        double **xPos, **yPos; // pointers to coordinates of points:
        // hierarchical: xPos is list of pointers to coordinates at each hierarchy level
        double **xRadii, **yRadii; // likewise: radii of each hierarchical cell, used to compute lower bounds
        double **alpha, **beta; // pointers to hierarchical dual variables, to compute effective costs if required
        //THierarchicalPartition *HP; // hierarchical partition classes
        double **muH; // hierarchical masses
};






struct XData
{
    int  ind[MAXLAYER + 1];
    double data[MAXSIZE * MAXSIZE * 3];
};



template <int MAXELEMENT>
 class COO_Sparese_Matrix
{
public:
    unsigned int nonZeros;
    int row_ptrl[MAXELEMENT];
    int col_ptrl[MAXELEMENT];
    double val[MAXELEMENT];
    int crs_row_ptrl[MAXSIZE*MAXSIZE + 1];
};


typedef  double COO_Buff;
typedef  int  Ind_Ptr_;
#endif // DATA_STRUCT_H
