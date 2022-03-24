#include"eps_handler.h"
#include"basic_settings.h"
#include<vector>
int CuEspsScalingHandler_::init_eps(TEpsScalingHandler &epsScalingHandler,bool refine){
    this->nlayers = epsScalingHandler.nLayers;
    if ( this->nlayers <= 4)
        return -1;
    for (int i = 0 ; i< this->nlayers; ++i){
        this->nEpsLists[i] = epsScalingHandler.nEpsLists[i];
        for (int j = 0; j < epsScalingHandler.nEpsLists[i]; ++j){
                this->epsLists[i][j] = epsScalingHandler.epsLists[i][j];
        }
    }

    if (refine){
        std::vector<double> eps_final_three_layer ;
        std::vector<double> eps_final_three_layer_start ;
        std::vector<int> neps_final_three_layer ;
        int cnt = 0;
        if (p_index < 1.4) {
            cnt = 4;
             neps_final_three_layer = std::vector<int>{ 5 ,10,10,8 };
            if (this->nlayers == 10) {
                cnt = 5;
                neps_final_three_layer = std::vector<int>{5, 10 ,10,10,8 };
                epsScalingHandler.epsTarget = 1e-3;
                eps_final_three_layer = std::vector<double>{   2  ,  0.1,   0.1,  0.00001, epsScalingHandler.epsTarget };
                // for the large problem, the eps plays an important role during iteration, 
                // Ithis eps is for solving OT between two images with the cost |x-y|^p
                // for point clouds, I suggest to scale the point's coordinates to [0, log2(n)];
                eps_final_three_layer_start = std::vector<double>{100,   100 ,  2  ,    0.1  ,     0.1 };
                 

            }
            if (this->nlayers == 9) {

                epsScalingHandler.epsTarget = 1e-2;
                eps_final_three_layer = std::vector<double>{  2, 0.1, 0.001, epsScalingHandler.epsTarget };
                 eps_final_three_layer_start = std::vector<double>{ 100, 2,  0.1 ,1 };

            }
            if (this->nlayers == 8) {
                epsScalingHandler.epsTarget = 1e-3;
                eps_final_three_layer = std::vector<double>{ 2,  1,    0.01, epsScalingHandler.epsTarget };
               eps_final_three_layer_start = std::vector<double>{ 100, 2, 0.1 , 0.01 };
                //eps_final_three_layer_start = std::vector<double>{ 20, 2, 1 , 1 };
            }
            if (this->nlayers <= 7) {
                epsScalingHandler.epsTarget = 1e-4;
                eps_final_three_layer = std::vector<double>{     2,    1,   0.1, epsScalingHandler.epsTarget };
                eps_final_three_layer_start = std::vector<double>{ 100, 20,  0.1 , 1 };
            }
        }else{ 
            cnt = 3;
            epsScalingHandler.epsTarget = 1e-1;
                // this eps solves the p  >= 1.5
            neps_final_three_layer = std::vector<int>{  10,10,8 };
            
            eps_final_three_layer = std::vector<double>{   5,  2 , epsScalingHandler.epsTarget };
            eps_final_three_layer_start = std::vector<double>{ 20, 6, 20 };
    }

    for (int ind = 0 ; ind < cnt; ++ind)
    {
        int i = this->nlayers - (cnt - ind);
        if(i > 2){
        this->nEpsLists[i] = neps_final_three_layer[ind];
        if(p_index <  1.5 || p_index > 2.2)
            this->epsLists[i][0] = eps_final_three_layer_start[ind];
        for (int j = 1; j < this->nEpsLists[i]; ++j){
           this->epsLists[i][j] = (this->epsLists[i][0]-eps_final_three_layer[ind])*exp(((double)-j))
                    + eps_final_three_layer[ind];
        }
        }
    }
    }

}
