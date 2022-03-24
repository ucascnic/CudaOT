#include"family.h"
#include"MultiScaleTools.h"
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
Family::Family(std::vector<int> &dim_,std::vector<double> &mu_,int nlayers):
    dim(dim_),mu(mu_)
    {
    this->totalnum = 1;
    for (int i = 0 ; i < dim.size() ; ++i){
        this->totalnum *= dim[i];
    }
    this->dimension = dim.size();
    this->nLayers = nlayers;
    this->nCells.resize(nLayers);


}


int Family::init(){

    TDoubleMatrix muX;
    std::vector<double> mu_ = this->mu;
    std::vector<int> dim_ = this->dim;
    muX.data=mu_.data();
    muX.dimensions=dim_.data();
    muX.depth=dim.size();

    TMultiScaleSetupSingleGrid  MultiScaleSetup(&muX,this->nLayers);
    int msg;
    msg=MultiScaleSetup.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetup.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetup.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }



    int nLayers = MultiScaleSetup.nLayers;
    printf("total layers %d\n",nLayers);
    int num_nChildren = 0;
    int num_parent = 0;
    for (int i = 0 ; i < nLayers ; ++i){
        this->nCells[i] =  MultiScaleSetup.HP->layers[i]->nCells;
        printf("i = %d\t,cellnum = %d\n",i,this->nCells[i]);

        for (int j = 0; j < this->nCells[i] ;++j){
            num_nChildren += MultiScaleSetup.HP->layers[i]->nChildren[j];
            num_parent += MultiScaleSetup.HP->layers[i]->parent[j];

//            printf("parent = %d\t",parent);
//            layersY->layer[i].nChildren[j] = MultiScaleSetupY.HP->layers[i]->nChildren[j];
//            layersY->layer[i].nLeaves[j] = MultiScaleSetupY.HP->layers[i]->nLeaves[j];
//            for (int k = 0 ; k < num_nChildren; ++k){
//              layersY->layer[i].children[j][k] =MultiScaleSetupY.HP->layers[i]->children[j][k];
//              layersY->layer[i].leaves[j][k] =MultiScaleSetupY.HP->layers[i]->leaves[j][k];
//            }
        }
        printf("\n");
    }
    printf("number num_nChildren = %d\t",num_nChildren);
    printf("number num_parent = %d\t",num_parent);
    int t = total_children(this->dimension,this->nLayers,this->totalnum);
    printf("t = %d\n",t);
    return 0;
}

//create a family based on 1dimension vector

