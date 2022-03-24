#ifndef BASIC_SETTINGS_H
#define BASIC_SETTINGS_H
#define MODE_MAX_ 1
#define MODE_MIN_ 0
#include <math.h>

// this calculate the promblem size for an optimal transproblem
// for a two-dimension image the problemsize = image size * image size


static constexpr int MSG_ABSORB_REITERATE=1;
static constexpr int MSG_EXCEEDMAXITERATIONS=30101;
static constexpr int MSG_NANSCALING=30102;
static constexpr int MSG_NANINERROR=30103;
static constexpr int MSG_ABSORB_TOOMANYABSORPTIONS=30201;


#define  OT_PROBLEM_SIZE(image_size)  (image_size*image_size)
#define MIN(a,b) ((a)<(b)?(a):(b))

#define GETLOG2(__MIN)  __MIN <= 32 ? 5 : ( __MIN<= 64 ? 6 : ( __MIN<= 128 ? 7 :  \
    ( __MIN<= 256 ? 8 : ( __MIN<= 512 ? 9 : ( __MIN<= 1024 ? 9 : ( 0 ) ) ) ) ) )
    
    
#define ALLOCALTE_MEMORY(__MIN)  __MIN <= 32 ? 19: ( __MIN<= 64 ? 21 : ( __MIN<= 128 ? 23 :  \
    ( __MIN<= 256 ? 25 : ( __MIN<= 512 ? 27 : ( __MIN<= 1024 ? 27 : ( 0 ) ) ) ) ) )

#define p_index 2


#define DIMENSIONS  2


#define DIM 2

#define MAXSIZE 512
#define MAXLAYER  GETLOG2( MAXSIZE )

//#define _RECORD_
#endif // BASIC_SETTINGS_H
