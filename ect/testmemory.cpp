#include<iostream>
struct XData

{
    int  ind[10];
    double data[500000];
};
#include<stdlib.h>

int main(){
     XData *muXH  = (XData*) malloc(sizeof(XData));
     std::cout << sizeof(XData) << std::endl;
     return 0;

}
