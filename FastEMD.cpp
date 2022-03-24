//// Fast Network Simplex for optimal mass transport -- test unit
//// by Nicolas Bonneel (Nov. 2013)
///
//// Typical runtime: 10k source nodes + 10k destination nodes, OpenMP enabled, 6 cores, 12 threads :
///  -  double precision: 7.8 seconds using 1.6GB of RAM (or 7.2 seconds with 2.4 GB of RAM without #define SPARSE_FLOW)
///  -  single precision: 3.5 seconds using 1.2GB of RAM
///  -  int : 5 seconds using 1.2 GB of RAM



#include <iostream>
#include <vector>

#include "FastTransport/network_simplex_simple.h"
#include <chrono>

#include <stdio.h>


using namespace lemon;

// all types should be signed
typedef int64_t arc_id_type; // {short, int, int64_t} ; Should be able to handle (n1*n2+n1+n2) with n1 and n2 the number of nodes (INT_MAX = 46340^2, I64_MAX = 3037000500^2)
typedef double supply_type; // {float, double, int, int64_t} ; Should be able to handle the sum of supplies and *should be signed* (a demand is a negative supply)
typedef double cost_type;  // {float, double, int, int64_t} ; Should be able to handle (number of arcs * maximum cost) and *should be signed*

struct TsFlow {
        int from, to;
        supply_type amount;
};

template<typename T>
cost_type sqrdistance2d(T* a, T* b) {
        return sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]));
}
#include<string>
#include<vector>
#include<Tools.h>

double  FastEMD_solver(int image_size,std::vector<double> &muXdat, std::vector<double> &muYdat);
int main(int argc, char** argv) {

    if (argc <= 1)
        return -1;
    int image_size = atoi(argv[1]);
    std::string filex, filey;


    int cnt = 0;

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= 1; ++i) {
        for (int j = 1; j <= 9; j++) {

            cnt++;
            printf("i=%d\n", i);
            filex = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(image_size) + \
                "_100" + std::to_string(i) + ".dat";
            filey = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(image_size) + \
                "_100" + std::to_string(j) + ".dat";
            std::vector<double> muXdat = readFile<double>(filex.c_str());
            std::vector<double> muYdat = readFile<double>(filey.c_str());

            std::chrono::high_resolution_clock::time_point periter = std::chrono::high_resolution_clock::now();
            double res = FastEMD_solver(image_size, muXdat, muYdat);
            std::chrono::high_resolution_clock::time_point periter_endTime = std::chrono::high_resolution_clock::now();

            std::cout << "time: " << (double)std::chrono::duration_cast<std::chrono::microseconds>(periter_endTime - periter ).count() / (1000000.) << std::endl;
            std::cout << res << std::endl;
        }
    }

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();

    std::cout << "time: " << (double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / (1000000. * (double)cnt) << std::endl;
    return 0;

}


double  FastEMD_solver(int image_size,std::vector<double> &muXdat, std::vector<double> &muYdat) {

        typedef FullBipartiteDigraph Digraph;
        DIGRAPH_TYPEDEFS(FullBipartiteDigraph);



        arc_id_type n1 = image_size*image_size, n2 = n1; // just demo for the case n1=n2 ; adapt otherwise
        std::vector<double> coordsXY(n1 * 2);
        std::vector<supply_type> weights1(n1), weights2(n2); // works faster with integer weights though
        for (int i = 0 ; i < n1 ; ++i){
            weights1[i] = muXdat[i];
            weights2[i] = -muYdat[i];
        }

        Digraph di(n1, n2);
        NetworkSimplexSimple<Digraph, supply_type, cost_type, arc_id_type> net(di, true, n1 + n2, n1*n2);

        // random point coordinates in 2D ; the same for both distributions here
        int cnt = 0;
        for (int i = 0; i < image_size ; i++) {
            for (int j = 0 ; j < image_size; j++){
                coordsXY[cnt] = i;
                coordsXY[cnt+1] = j;
                cnt += 2;
            }
         }

        arc_id_type idarc = 0;
        for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++) {

                        cost_type d = sqrdistance2d(&coordsXY[i * 2], &coordsXY[j * 2]);

                        net.setCost(di.arcFromId(idarc), d);
                        idarc++;
                }
        }


        net.supplyMap(&weights1[0], n1, &weights2[0], n2);

        // We run the network simplex
       
        int ret = net.run();

        double resultdist = net.totalCost();  // resultdist is the EMD


         // here is how to get the flow
        std::vector<TsFlow> flow;
        flow.reserve(n1 + n2 - 1);
        double recomputed_cost = 0;
        for (arc_id_type i = 0; i < n1; i++) {
                for (arc_id_type j = 0; j < n2; j++)
                {
                        TsFlow f;
                        f.amount = net.flow(di.arcFromId(i*n2 + j));
                        f.from = i;
                        f.to = j;

                        recomputed_cost += f.amount*sqrdistance2d(&coordsXY[i * 2], &coordsXY[j * 2]);			 // another way to get the EMD ; same as "resultdist"
                        if (fabs((double)f.amount) > 1E-18)
                                flow.push_back(f);
                }
        }
        //std::cout << "EMD = "<<resultdist<<" == "<<recomputed_cost << std::endl;

        return resultdist;
}
