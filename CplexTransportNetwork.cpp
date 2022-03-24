
#include<LP_CPLEX.h>
#include<Tools.h>
#include<cstdlib>
#include <vector>
#include<iostream>
#include<Common.h>


template<typename T>
double sqrdistance2d(T* a, T* b) {
        return sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]));
}


template<int XRow,int Ycol>
double CplexSolver(std::vector<double> &, std::vector<double> &);


double TransportNetworkSolver(std::vector<double> &muXdat, std::vector<double> &muYdat,int image_size){
        double res = 0.0;
        switch (image_size) {
        case 32:
            res =  CplexSolver<32,32>( muXdat,  muYdat);
 
            break;

        case 64:
            res =  CplexSolver<64,64>( muXdat,  muYdat);
 
            break;
        case 128:
            res =  CplexSolver<128,128>( muXdat,  muYdat);
            break;
        default:
            break;
        }
        return res;
  
}

int main(int argc, char* argv[]){
 
    if (argc <= 1)
        return -1;
    int image_size = atoi(argv[1]);
    std::string filex,filey;

 

    clock_t begin = clock();
    int cnt = 0;
    for (int i = 1;i <= 9;++i){
        for (int j = 2; j <= 2; j++){
             
        cnt++;
        printf("i=%d\n",i);
        filex = "../data/ot_benchmark/CauchyDensity/data"+ std::to_string(image_size)+ \
                "_100"+std::to_string(i)+".dat";
        filey = "../data/ot_benchmark/CauchyDensity/data"+ std::to_string(image_size)+ \
                "_100"+std::to_string(j)+".dat";
        std::vector<double> muXdat=readFile<double>(filex.c_str());
        std::vector<double> muYdat=readFile<double>(filey.c_str());
        double res =  TransportNetworkSolver(muXdat,muYdat,image_size);
        std::cout << res << std::endl;
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "\n" << elapsed_secs/(double)cnt << "s" << std::endl;



}


template<int XRow,int Ycol>
double CplexSolver(std::vector<double> &muXdat, std::vector<double> &muYdat ){


    int n1 = muXdat.size();
    int n2 = muYdat.size();
    std::vector<double> coordsXY(n1 + n2);
    int cnt = 0;
    for (int i = 0; i < XRow ; i++) {
        for (int j = 0 ; j < Ycol; j++){
            coordsXY[cnt] = i;
            coordsXY[cnt+1] = j;
            cnt += 2;
        }
     }
    double *c = (double *) malloc(sizeof(double) * n1 * n2);
    int idarc = 0;
     for (int i = 0; i < n1; i++) {
             for (int j = 0; j < n2; j++) {

                     c[idarc] = sqrdistance2d(&coordsXY[i * 2], &coordsXY[j * 2]);
                     idarc++;
             }
     }
    double * mu = (double *) malloc(sizeof(double) * n1 * n2);
    TCouplingHandlerDensePrototype<double>  coupling(n1,n2,c, mu);

    TCPLEXNetSolver<TCouplingHandlerDensePrototype<double>>
            solver(&coupling,muXdat.data(),muYdat.data());


    solver.setup();

    solver.solve();

    double res =  solver.getObjective();
    free(mu);
    free(c);
    return res;


}
