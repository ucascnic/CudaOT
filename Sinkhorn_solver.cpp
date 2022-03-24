#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<Common.h>
#include<string>
#include<vector>
#include<omp.h>
#include<ctime>
#include<time.h>
#define n 64
#define XSIZE n*n
#define YSIZE n*n
#define PARALLEL
constexpr int KERNALSIZE =  XSIZE*YSIZE*sizeof(double);
using namespace  std;

constexpr int nlist[20] = {5,5,5,5,5,20,20,20,20,20,20,20,20,20,20,20,20,20,20,100};
typedef double  datatype;
//for (int j = 0;j < YSIZE;++j){
//    int ind = i + XSIZE*j;
//    kernal[ind] = c[ind] - g[j];
//    xmindata = xmindata<kernal[ind]?xmindata:kernal[ind];
//}



int main(int argc, char* argv[]);
int W2Grid();
int Interpolation();
int WpSphere();
int WFR();
int W2Grid_256x256();
int modified();

template<int j>class Calsum{
public:
    static void cal(double *res,  double *  xmindat,
                      double *  kernal,double *_eps,int i){
        Calsum<j-1>::cal(res,xmindat,kernal,_eps,i);
        *res +=  exp((*xmindat - kernal[i + XSIZE * j])/ (*_eps)) ;
    }
};
template<>class Calsum<-1>{
public:
    static void cal(double *res,  double *  xmindat,
                    const  double *  kernal,double *_eps,int i){}
};
template<int inneriter>class Innerinter{
public:
    static void per_sinkhorn(double *_eps,double *kernal,double *c,double *f,double *g,
                             double *logx,double *logy){
        Innerinter<inneriter-1>::per_sinkhorn(_eps,kernal,c,f,g,
                                              logx,logy);


          for (int i = 0 ; i <3 ; ++i)
              cout << c[i];
           cout  << endl;
for (int kk=0;kk<2;++kk)
{
    *_eps = (1e6 - 1e-6) * std::exp(-kk) + 1e-6;
            memcpy(kernal,c,KERNALSIZE);
//#pragma omp parallel for
            for(int i = 0 ; i< XSIZE;++i){
                double xmindata = 1e15;

                for (int j = 0;j < YSIZE;++j){
                    int ind = i * XSIZE +  j;
                    kernal[ind] -=  g[j] ;
                    xmindata = xmindata<kernal[ind]?xmindata:kernal[ind];
                }


                double resi = 0.0;
                for (int j =0 ; j < YSIZE; ++j){
                      resi +=  exp((xmindata - kernal[i * XSIZE+ j])/ *_eps) ;
                }


                f[i] =  xmindata +   (*_eps) *  ( logx[i] -std::log(resi) ) ;

            }

//            memcpy(kernal,c,KERNALSIZE);
//#pragma omp parallel for
            for(int j = 0 ; j< YSIZE; ++j){
                double mindata = 1e15;
                for (int i = 0; i  < XSIZE; ++i){
                    int ind = i * XSIZE+ j;
                    kernal[ind] -=  f[i]  ;
                    mindata = mindata<kernal[ind]?mindata:kernal[ind];
                }

                double resj = 0.0;
                for (int i =0 ; i < XSIZE; ++i){
                      resj +=  exp((mindata - kernal[i * XSIZE+ j])/ *_eps) ;

                }
               g[j] +=  mindata + (*_eps) *  (logy[j] - std::log(resj) ) ;

        }
            for (int i = 0 ; i<3; ++i)
                std::cout << "gg = "  << g[i] << "  ";
            std::cout << std::endl;
            for (int i = 0 ; i<3; ++i)
                std::cout << "ff = " << f[i] << "  ";std::cout << std::endl;
            if (kk==1){ exit(0);}
//            for(int i = 0 ; i < 3;++i)
//                cout << f[i] << endl;
//            exit(0);
//            for(int i = 0 ; i < 3;++i)
//                cout << g[i] << endl;
}

//        double result = 0.0;
//        for(int i = 0 ; i < XSIZE ;++i)
//            for(int j = 0; j < YSIZE; ++j){
//                result +=  std::exp((f[i] + g[j] - c[i+j*YSIZE])/ *_eps);
//                cout << result << endl;
//                if (j > 10)
//                    exit(0);
//            }

//


    }

};

template<>class Innerinter<-1>{
public:
    static void per_sinkhorn(double *_eps,double *kernal,double *c,double *f,double *g,
                             double *logx,double *logy){
    }
};

template<int iter> class Iter
{
public:

    static void  periter(double *_eps,
                         double *kernal,double *c,double *f,double *g,double *logx,double *logy)
    {

        Iter<iter-1>::periter(_eps,kernal,c,f,g, logx, logy);

        *_eps = (1e6 - 1e-6) * std::exp(-iter) + 1e-6;
        cout << iter << endl;
        {
        Innerinter<nlist[iter]-1>::per_sinkhorn(_eps,kernal,c,f,g,
                                                  logx,logy);
        }
    }
};
template<> class Iter<-1>
{
public:
    static void  periter(double *eps,
                         double *kernal,double *c,double *f,double *g,double *logx,double *logy)
    {}

};

int main(int argc, char* argv[]) {
//        W2Grid();
    modified();
}

int modified(){
    string filex,filey;

    filex = "/home/dongdong/data/ot_benchmark/CauchyDensity/data"+ to_string(n)+"_1001.dat";
    filey = "/home/dongdong/data/ot_benchmark/CauchyDensity/data"+ to_string(n)+"_1002.dat";


    std::vector<double> muXdat=readFile<double>(filex.c_str());
    std::vector<double> muYdat=readFile<double>(filey.c_str());

    // dimensions of grid
    int dim = 2;
    std::vector<int> muXdim = {n ,n };
    std::vector<int> muYdim = {n ,n };


    // put raw marginal data into small container structs
    TDoubleMatrix muX,muY;



    muX.data=muXdat.data();
    muX.dimensions=muXdim.data();

    muY.data=muYdat.data();
    muY.dimensions=muYdim.data();


    TDoubleMatrix *posX=GridToolsGetGridMatrix(dim, muX.dimensions);
    TDoubleMatrix *posY=GridToolsGetGridMatrix(dim, muY.dimensions);

    int xres=posX->dimensions[0];
    int yres=posY->dimensions[0];
    TCostFunctionProvider_Dynamic costFunctionProvider(
        &xres, &yres,
        &(posX->data), &(posY->data),
        1, dim);

    int datasize = xres * yres;
    int kernelsize = datasize * (sizeof(datatype));

    datatype *c = (datatype *)malloc(kernelsize);
    datatype *kernal = (datatype *)malloc(kernelsize);

    costFunctionProvider.getCDense(c);



    datatype f[XSIZE]={0.0};
    datatype g[YSIZE]={0.0};
    datatype logx[XSIZE];
    datatype logy[YSIZE];

    for (int i = 0 ;i < xres ; ++i)
        logx[i] = std::log(muXdat[i]);
    for (int i = 0 ;i < yres ; ++i)
        logy[i] = std::log(muYdat[i]);

    double eps = 1e-8;
    int iternum = 20;
//    for (int i = 0 ;i < iternum-1 ; ++i){
//        nlist[i] = 20;
//    }

//    nlist[iternum-1] = 100;
//      omp_set_num_threads(512);


      double _eps;
#ifdef PARALLEL
      double start = omp_get_wtime();

#else

      clock_t start,end;
      start = clock();
#endif

//      cout <<"hello" << omp_get_thread_num() <<  endl;
      Iter<19>::periter(&_eps,kernal,c,f,g,logx,logy);
#ifdef PARALLEL
     double end = omp_get_wtime();
     printf("time=%f\n",(double) (end-start));
#else
     printf("time=%f\n",((double)(clock() - start)/CLOCKS_PER_SEC));
#endif

//    for (int iter = 0 ; iter < iternum; ++iter){
//         _eps = (1e6 - eps) * std::exp(-iter) + eps;
//        clock_t start,end;
//        for(int k=0;k<nlist[iter];k++) {
//            start = clock();
//#pragma omp parallel for
//            for(int i = 0 ; i< XSIZE;++i){
//                double xmindata = 1e10;
//                for (int j = 0;j < YSIZE;++j){
//                    int ind = i + XSIZE*j;
//                    kernal[ind] = c[ind] - g[j];
//                    xmindata = xmindata<kernal[ind]?xmindata:kernal[ind];
//                }
//                double resi = 0.0;
//                for (int j =0 ; j < YSIZE; ++j){
//                    resi +=  exp((xmindata - kernal[i + XSIZE*j])/_eps) ;
//                }
//                f[i] = xmindata - _eps *  (std::log(resi) -   logx[i]) ;

//            }
//            end = clock();
//            printf("time=%f\n",(double) (end-start)/CLOCKS_PER_SEC);
//#pragma omp parallel for
//            for(int j = 0 ; j< YSIZE; ++j){
//                double mindata = 1e10;
//                for (int i = 0; i < XSIZE; ++i){
//                    int ind = i + XSIZE*j;
//                    kernal[ind] = c[ind] - f[i];
//                    mindata = std::min(mindata,kernal[ind]);
//                }
//                double resj = 0.0;

//#pragma omp parallel for reduction(+:resj)
//                for (int i =0 ; i < XSIZE; ++i){
//                    int ind = i + XSIZE*j;
//                    resj +=  exp((mindata - kernal[ind])/_eps) ;
//                }
//                g[j] = mindata - _eps *  (std::log(resj) -   logy[j]) ;


//            }
//        }


//    }
      double result = 0.0;
      for(int i = 0 ; i < XSIZE ;++i)
          for(int j = 0; j < YSIZE; ++j){
              result += std::exp((f[i] + g[j] - c[i+j*YSIZE])/ _eps) * c[i+j*YSIZE]  ;
          }
      cout <<result <<endl;

    free(c);
    free(kernal);
    return 0;
}

