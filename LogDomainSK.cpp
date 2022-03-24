

#include<LogSinkhorn/log_sinkhorn.h>

#include<time.h>
#include<Tools.h>
int main(){

    std::string filex, filey;
    int N = 64;
    filex = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(N) + "_1001.dat";
    filey = "../data/ot_benchmark/CauchyDensity/data" + std::to_string(N) + "_1002.dat";
    std::vector<double> muXdat = readFile<double>(filex.c_str());
    std::vector<double> muYdat = readFile<double>(filey.c_str());
    switch (N) {
    case 32:
        test_32x32(muXdat,muYdat);
        break;
    case 64:
        test_64x64(muXdat,muYdat);
        break;
    default:
        break;
    }

}

