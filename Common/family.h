#ifndef FAMILY_H
#define FAMILY_H
#include"MultiScaleTools.h"
#include<vector>
#include <math.h>
#include<iostream>
#include"basic_settings.h"

inline int total_children(int dimension,int nlayers,int totalnum){

    int mn = dimension * nlayers;
    int n1 = std::pow(2,mn);
    int n2 = std::pow(2,dimension);
    return (n1-n2)/(n2-1) + totalnum;
}
#define MAX(a, b) ((a)>(b)?(a):(b))


class Family
{
public:
    int dimension;
    int nLayers;
    const std::vector<int> dim;
    const std::vector<double> mu;
    unsigned int totalnum;

    std::vector<int> nCells; // how many cell at each layer?
    std::vector<int> children_at_each_layer;

    Family(std::vector<int> &dim_, std::vector<double> &mu_, int nlayers);
    int init();



};

inline unsigned int get_children_ind(int layer){
    return ((1 << (DIM * (layer + 1 ))) - (1 << DIM ))/((1 << DIM) - 1) ;
}
//inline unsigned int get_n_children_ind(int layer){
//    return 1 + ((1 << (DIM * (layer))) - (1 << DIM ))/((1 << DIM) - 1);
//}

#define get_n_children_ind(layer)  (1 + ((1 << (DIM * (layer))) - (1 << DIM ))/((1 << DIM) - 1))

template <int TOTALNUM_,int MAXLAYER_>
class Familys
{
public:
    int nlayers;
    int nCells[MAXLAYER_ + 1]; // [1 4 16 64  ... 4096]

    int nChildren_length;
    int nChildren[2 * TOTALNUM_];

    int nChildren_ind_length;
    int nChildren_ind[2 * TOTALNUM_];

    int children_length;
    int children[2 * TOTALNUM_];

    int check_for_layers(TMultiScaleSetupSingleBase *MultiScaleSetupX);
};

template <int TOTALNUM_,int MAXLAYER_>
Familys<TOTALNUM_,MAXLAYER_> * createFamily
(TMultiScaleSetupSingleBase *MultiScaleSetupX){

    Familys<TOTALNUM_,MAXLAYER_>  *familyX =(Familys<TOTALNUM_,MAXLAYER_> *)
            malloc(sizeof(Familys<TOTALNUM_,MAXLAYER_>));
    familyX->nChildren_ind[0] = 0;
    familyX->nlayers = MultiScaleSetupX->nLayers;

    for (int i = 0 ; i < familyX->nlayers ; ++i){
        familyX->nCells[i] = MultiScaleSetupX->HP->layers[i]->nCells;
        int nchildren_ind = get_n_children_ind(i);
        for (int j = 0; j < MultiScaleSetupX->HP->layers[i]->nCells;++j){
            int num_nChildren = MultiScaleSetupX->HP->layers[i]->nChildren[j];

            familyX->nChildren[nchildren_ind + j] = MultiScaleSetupX->HP->layers[i]->nChildren[j];
            familyX->nChildren_ind[nchildren_ind + j + 1] =
                     familyX->nChildren_ind[nchildren_ind + j] + num_nChildren;
             for (int k = 0 ; k < num_nChildren; ++k){
              familyX->children[familyX->nChildren_ind[nchildren_ind + j] +  k]
                      =  MultiScaleSetupX->HP->layers[i]->children[j][k];
            }

        }

    }
    return familyX;

}

template <int TOTALNUM_,int MAXLAYER_>
int Familys<TOTALNUM_,MAXLAYER_>::check_for_layers(TMultiScaleSetupSingleBase *MultiScaleSetupX){

    for (int i = 0 ; i < this->nlayers ; ++i){
        if ( MultiScaleSetupX->HP->layers[i]->nCells != this->nCells[i]){
            printf("error! for i");exit(0);
        }
       int nchildren_ind = get_n_children_ind(i);
        for (int j = 0; j < MultiScaleSetupX->HP->layers[i]->nCells;++j){
            int num_nChildren = MultiScaleSetupX->HP->layers[i]->nChildren[j];
            if (MultiScaleSetupX->HP->layers[i]->nChildren[j] !=
                    this->nChildren[nchildren_ind + j]){
                printf("error! for j");exit(0);
            }
            for (int k = 0 ; k < num_nChildren; ++k){
               if (this->children[this->nChildren_ind[nchildren_ind + j] +  k] !=
                                  MultiScaleSetupX->HP->layers[i]->children[j][k])
               {
                   printf("error for k!\n");
                     exit(0);

               }

            }
        }
    }
}


#endif // FAMILY_H
