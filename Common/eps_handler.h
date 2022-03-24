#ifndef EPS_HANDLER_H_
#define EPS_HANDLER_H_


#include<TEpsScaling.h>
#ifndef MAXEPSLIST
    #define MAXEPSLIST 20
#endif
#ifndef MAXEPS_EACH_LAYER
    #define MAXEPS_EACH_LAYER 30
#endif
class CuEspsScalingHandler_
{
public:
    int nlayers;
    int nEpsLists[MAXEPSLIST]; // number of eps in each layer
    double epsLists[MAXEPSLIST][MAXEPS_EACH_LAYER]; // one eps list per layer
    CuEspsScalingHandler_(){};
    int init_eps(TEpsScalingHandler &,bool refine = true);
};
#endif // EPS_HANDLER_H
