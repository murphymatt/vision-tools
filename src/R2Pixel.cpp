// Source file for the R2Pixel class 



// Include files 

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include "R2Pixel.h"



// Public variables 

R2Pixel R2null_pixel(0.0, 0.0, 0.0, 1.0);
R2Pixel R2black_pixel(0.0, 0.0, 0.0, 1.0);
R2Pixel R2red_pixel(1.0, 0.0, 0.0, 1.0);
R2Pixel R2green_pixel(0.0, 1.0, 0.0, 1.0);
R2Pixel R2blue_pixel(0.0, 0.0, 1.0, 1.0);
R2Pixel R2yellow_pixel(1.0, 1.0, 0.0, 1.0);
R2Pixel R2cyan_pixel(0.0, 1.0, 1.0, 1.0);
R2Pixel R2magenta_pixel(1.0, 0.0, 1.0, 1.0);
R2Pixel R2white_pixel(1.0, 1.0, 1.0, 1.0);


// Pixel Functions
double R2Pixel::
PixelDistance(R2Pixel * p)
{
  return sqrt(pow(Red() - p->Red(), 2) +
	      pow(Green() - p->Green(), 2) +
	      pow(Blue() - p->Blue(), 2));
}


double R2Pixel::
GetPixelWeight(R2Pixel * p)
{
  return 2 / (1 + exp(-1 * PixelDistance(p))) - 1;
}


R2Pixel * R2Pixel::
BlendByWeight(R2Pixel * p, double l, double h)
{
  R2Pixel* ret = new R2Pixel(*this);
  //ret->Clamp((1 - PixelDistance(p)) * h + l);
  *ret = *ret * GetPixelWeight(p);
  return ret;
}
