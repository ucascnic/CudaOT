#include<SparseSinkhorn/image_ot_solver.h>
#include<string>
#include<Tools.h>
#include<iostream>
#include <ctime>
// this function tests the full data set of the DOTbechmark
// the file path of the data should be adjusted before running the program
// this file is the first example I created to test the performance of the
// OT between two images, it shows a good performance compared with the LP
// with regard to the meomory usage and solving time.
// the next importange job I should do is to realize the OT between two color images!



int main2(int argc, char **argv){


    //if (argc <= 1)
     //   return -1;
    //int image_size = atoi(argv[1]);
    int image_size = 64;
    std::string filex,filey;
    std::vector<double> result;
 

    clock_t begin = clock();
    int cnt = 0;
    for (int i = 1;i <= 5; i=i+1){
        for (int j = 1; j <= 165; j = j + 11){
             
        cnt++;
 
        filex = "../data/Yaledat/Yale"+ std::to_string(i)+".dat";
        filey = "../data/Yaledat/Yale"+ std::to_string(j)+".dat";
        std::vector<double> muXdat=readFile<double>(filex.c_str());
        std::vector<double> muYdat=readFile<double>(filey.c_str());
        double res  = image_ot_solver(image_size,muXdat,muYdat);
        result.push_back(res);
        //printf("%.4f\n", res);
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "\n" << elapsed_secs/(double)cnt << "s" << std::endl;
    
    
    for (int i = 0;i < 5 ;++i){
    for (int j = 0;j < 15;++j)
      std::cout << result[i*15+j] << " " ;
    std::cout << " " << std::endl;}
    return 0;

}

int main(int argc, char **argv){


    if (argc <= 1)
        return -1;
    int image_size = atoi(argv[1]);
    std::string filex,filey;

 

    clock_t begin = clock();
    int cnt = 0;
    for (int i = 1;i <= 1;++i){
        for (int j = 7; j <= 7; j++){
             
        cnt++;
 
        filex = "../data/ot_benchmark/ClassicImages/picture"+ std::to_string(image_size)+ \
                "_100"+std::to_string(i)+".dat";
        filey = "../data/ot_benchmark/ClassicImages/picture"+ std::to_string(image_size)+ \
                "_100"+std::to_string(j)+".dat";
        std::vector<double> muXdat=readFile<double>(filex.c_str());
        std::vector<double> muYdat=readFile<double>(filey.c_str());
        double res  = image_ot_solver(image_size,muXdat,muYdat);
        printf("%.4f\n", res);
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "\n" << elapsed_secs/(double)cnt << "s" << std::endl;
    return 0;

}





