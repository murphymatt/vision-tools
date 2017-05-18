// Source file for image class

#include <cmath>
#include <limits>
#include <string>
// #include <thread>
#include <utility>
#include <vector>

// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"


struct pointComp
{
  inline bool operator() (std::pair< R2Point*, double > p1,
			  std::pair< R2Point*, double > p2)
  {
    return p1.second > p2.second;
  }
};


struct pixelComp
{
  inline bool operator() (R2Pixel * p1, R2Pixel * p2)
  {
    return p1->Luminance() > p2->Luminance();
  }
};


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}


R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
  // fit a 2D conic to five points
  R2Point p1(1.2,3.5);
  R2Point p2(2.1,2.2);
  R2Point p3(0.2,1.6);
  R2Point p4(0.0,0.5);
  R2Point p5(-0.2,4.2);

  // build the 5x6 matrix of equations
  double** linEquations = dmatrix(1,5,1,6);

  linEquations[1][1] = p1[0]*p1[0];
  linEquations[1][2] = p1[0]*p1[1];
  linEquations[1][3] = p1[1]*p1[1];
  linEquations[1][4] = p1[0];
  linEquations[1][5] = p1[1];
  linEquations[1][6] = 1.0;

  linEquations[2][1] = p2[0]*p2[0];
  linEquations[2][2] = p2[0]*p2[1];
  linEquations[2][3] = p2[1]*p2[1];
  linEquations[2][4] = p2[0];
  linEquations[2][5] = p2[1];
  linEquations[2][6] = 1.0;

  linEquations[3][1] = p3[0]*p3[0];
  linEquations[3][2] = p3[0]*p3[1];
  linEquations[3][3] = p3[1]*p3[1];
  linEquations[3][4] = p3[0];
  linEquations[3][5] = p3[1];
  linEquations[3][6] = 1.0;
  
  linEquations[4][1] = p4[0]*p4[0];
  linEquations[4][2] = p4[0]*p4[1];
  linEquations[4][3] = p4[1]*p4[1];
  linEquations[4][4] = p4[0];
  linEquations[4][5] = p4[1];
  linEquations[4][6] = 1.0;

  linEquations[5][1] = p5[0]*p5[0];
  linEquations[5][2] = p5[0]*p5[1];
  linEquations[5][3] = p5[1]*p5[1];
  linEquations[5][4] = p5[0];
  linEquations[5][5] = p5[1];
  linEquations[5][6] = 1.0;

  printf("\n Fitting a conic to five points:\n");
  printf("Point #1: %f,%f\n",p1[0],p1[1]);
  printf("Point #2: %f,%f\n",p2[0],p2[1]);
  printf("Point #3: %f,%f\n",p3[0],p3[1]);
  printf("Point #4: %f,%f\n",p4[0],p4[1]);
  printf("Point #5: %f,%f\n",p5[0],p5[1]);

  // compute the SVD
  double singularValues[7]; // 1..6
  double** nullspaceMatrix = dmatrix(1,6,1,6);
  svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

  // get the result
  printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

  // make sure the solution is correct:
  printf("Equation #1 result: %f\n",  p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
	 p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
	 p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
	 p1[0]*nullspaceMatrix[4][smallestIndex] + 
	 p1[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  printf("Equation #2 result: %f\n",  p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
	 p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
	 p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
	 p2[0]*nullspaceMatrix[4][smallestIndex] + 
	 p2[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  printf("Equation #3 result: %f\n",  p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
	 p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
	 p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
	 p3[0]*nullspaceMatrix[4][smallestIndex] + 
	 p3[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  printf("Equation #4 result: %f\n",  p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
	 p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
	 p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
	 p4[0]*nullspaceMatrix[4][smallestIndex] + 
	 p4[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  printf("Equation #5 result: %f\n",  p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
	 p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
	 p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
	 p5[0]*nullspaceMatrix[4][smallestIndex] + 
	 p5[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  R2Point test_point(0.34,-2.8);

  printf("A point off the conic: %f\n", test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
	 test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
	 test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
	 test_point[0]*nullspaceMatrix[4][smallestIndex] + 
	 test_point[1]*nullspaceMatrix[5][smallestIndex] + 
	 nullspaceMatrix[6][smallestIndex]);

  return; 
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
  // create copy of original image
  R2Image * tmp = new R2Image(this->width, this->height);

  // Apply the Sobel oprator to the image in X direction
  int H[3][3] = {{1, 0, -1},
		 {2, 0, -2},
		 {1, 0, -1}};
  
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      tmp->Pixel(x,y) = 
	Pixel(x-1,y-1)*H[0][0] + Pixel(x+1,y-1)*H[0][2] +
	Pixel(x-1,y)*H[1][0] + Pixel(x+1,y)*H[1][2] +
	Pixel(x-1,y+1)*H[2][0] + Pixel(x+1,y+1)*H[2][2];

      /*
	double red = tempImage.Pixel(x,y).Red();
	double green = tempImage.Pixel(x,y).Green();
	double blue = tempImage.Pixel(x,y).Blue();
	double alpha = tempImage.Pixel(x,y).Alpha();
      */
      // normalize to grayscale
      //tempImage.Pixel(x,y).SetRed((red + blue + green)/3 + 0.5);
      //tempImage.Pixel(x,y).SetGreen((red + blue + green)/3 + 0.5);
      //tempImage.Pixel(x,y).SetBlue((red + blue + green)/3 + 0.5);
      //tempImage.Pixel(x,y).SetAlpha(alpha + 0.5);
      
      // tempImage.Pixel(x,y).Clamp();
    }
  }
  *this = *tmp;
}

void R2Image::
SobelY(void)
{
  // create copy of original image
  R2Image * tmp = new R2Image(this->width, this->height);
  
  // Apply the Sobel oprator to the image in Y direction
  int H[3][3] = {{-1, -2, -1},
		 {0, 0, 0},
		 {1, 2, 1}};
  
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      tmp->Pixel(x,y) =
	Pixel(x-1,y-1)*H[0][0] + Pixel(x-1,y+1)*H[2][0] + 
	Pixel(x,y-1)*H[0][1] + Pixel(x,y+1)*H[2][1] +
	Pixel(x+1,y-1)*H[0][2] + Pixel(x+1,y+1)*H[2][2];

      /*
	double red = tempImage.Pixel(x,y).Red();
	double green = tempImage.Pixel(x,y).Green();
	double blue = tempImage.Pixel(x,y).Blue();
	double alpha = tempImage.Pixel(x,y).Alpha();
      */ 
      // normalize to grayscale
      // tempImage.Pixel(x,y).SetRed((red + blue + green)/3 + 0.5);
      // tempImage.Pixel(x,y).SetGreen((red + blue + green)/3 + 0.5);
      // tempImage.Pixel(x,y).SetBlue((red + blue + green)/3 + 0.5);
      // tempImage.Pixel(x,y).SetAlpha(alpha + 0.5);
      
      // tempImage.Pixel(x,y).Clamp();
    }
  }
  *this = *tmp;
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  BilateralFilter(4, sigma, 0);
}


void R2Image::
BlurD(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
  const double PI = 3.1415926535;
  const double EULER = 2.7182818285;
  
  // instantiate and resize kernel
  std::vector<double> kernel;
  int s = int(sigma);
  int dim = (s * 6) + 1;
  kernel.resize(dim);

  // compute kernel
  for (int i=0; i<dim; i++) {
    kernel[i] = pow(EULER, -1 * pow(i-(3*s),2) / (2.*pow(sigma,2))) 
      / (sqrt(2. * PI) * sigma);
  }

  R2Image * tmp = new R2Image(this->width, this->height);

  // apply kernel to image over y direction
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      R2Pixel * val = new R2Pixel();
      double weights = 0;
      for (int ly=(-3*s); ly<=(3*s); ly++) {
        if (y+ly>=0 && y+ly<height) {
	  *val += Pixel(x,y+ly)*kernel[ly+3*s];
	  weights += kernel[ly+3*s];
        }
      }
      *val /= weights;
      tmp->Pixel(x,y) = *val;
      // tempImage.Pixel(x,y).Clamp();
    }
  }

  // apply kernel to image over x direction
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      R2Pixel * val = new R2Pixel();
      double weights = 0;
      for (int lx=(-3*s); lx<=(3*s); lx++) {
        if (x+lx>=0 && x+lx<width) {
	  *val += tmp->Pixel(x+lx,y)*kernel[lx+3*s];
	  weights += kernel[lx+3*s];
        }
      }
      *val /= weights;
      Pixel(x,y) = *val;
      // Pixel(x,y).Clamp();
    }
  }
}


void R2Image::
DrawBox(int x, int y, bool good)
{
  int sideLen = 3, low, high;
    
  low  = (x - sideLen >= 0)     ? x-sideLen : 0;
  high = (x + sideLen <= width) ? x+sideLen : width;
  for (int lx = low; lx < high; lx++) {
    if (good) {
      if (y-sideLen >= 0)     Pixel(lx, y-sideLen) = R2green_pixel;
      if (y+sideLen < height) Pixel(lx, y+sideLen) = R2green_pixel;
    } else {
      if (y-sideLen >= 0)     Pixel(lx, y-sideLen) = R2red_pixel;
      if (y+sideLen < height) Pixel(lx, y+sideLen) = R2red_pixel;
    }
  }

  low  = (y - sideLen >= 0)     ? y-sideLen : 0;
  high = (y + sideLen <= width) ? y+sideLen : height;
  for (int ly = low; ly < high; ly++) {
    if (good) {
      if (x-sideLen >= 0)    Pixel(x-sideLen, ly) = R2green_pixel;
      if (x+sideLen < width) Pixel(x+sideLen, ly) = R2green_pixel;  
    } else {
      if (x-sideLen >= 0)    Pixel(x-sideLen, ly) = R2red_pixel;
      if (x+sideLen < width) Pixel(x+sideLen, ly) = R2red_pixel;
    }
  }
}


// helper function to determine harris score of each pixel s.t. score = r + g + b + a
double R2Image::
getPixelMagnitude(int x, int y)
{
  return Pixel(x,y).Red() + Pixel(x,y).Green() + Pixel(x,y).Blue() + Pixel(x,y).Alpha();
}


std::vector< R2Point* > R2Image::
GetBestFeatures(void)
{
  std::vector< std::pair< R2Point*, double > > points;

  for (int y=5; y<height-5; y++) {
    for (int x=5; x<width-5; x++) {
      std::pair <R2Point*, double> pair (new R2Point(x,y), Pixel(x,y).Luminance());
      points.push_back(pair);
    }
  }

  // sort points by harris score
  std::sort(points.begin(), points.end(), pointComp());

  // determine 150 greatest threshold areas with distance at least 10 points
  std::vector< R2Point* > feats;
  feats.push_back(points.at(0).first);
  R2Point* cur = feats.at(0);
  bool isFar;
  int topScoreCount = 1;

  for (int i=1; i<points.size(); i++) {
    isFar = true;
    cur = points.at(i).first;
    for (int j=0; j<feats.size(); j++) {
      if (abs(cur->X() - feats.at(j)->X()) < 10 &&
	  abs(cur->Y() - feats.at(j)->Y()) < 10) {
	isFar = false;
	break;
      }
    }
    if (isFar) {
      feats.push_back(cur);
      topScoreCount++;
    }
    if (topScoreCount >= 150) break;
  }

  return feats;
}


void R2Image::
Harris(double sigma)
{
  // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
  // Output should be 50% grey at flat regions, white at corners and black/dark near edges
  
  // create ix and iy components
  R2Image * ix = new R2Image(this->width, this->height, this->pixels);
  ix->SobelX();
  R2Image * iy = new R2Image(this->width, this->height, this->pixels);
  iy->SobelY();

  // create 3 images from ix and iy components
  R2Image * img1 = new R2Image(this->width, this->height);
  R2Image * img2 = new R2Image(this->width, this->height);
  R2Image * img3 = new R2Image(this->width, this->height);
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      img1->Pixel(x,y) = ix->Pixel(x,y) * ix->Pixel(x,y);
      img2->Pixel(x,y) = iy->Pixel(x,y) * iy->Pixel(x,y);
      img3->Pixel(x,y) = ix->Pixel(x,y) * iy->Pixel(x,y);
    }
  }
 
  delete ix;
  delete iy;

  // blur each temporary image
  img1->Blur(2.0);
  img2->Blur(2.0);
  img3->Blur(2.0);

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      Pixel(x,y) = 
	img1->Pixel(x,y)*img2->Pixel(x,y) - 
	img3->Pixel(x,y)*img3->Pixel(x,y) -
	sigma * ((img1->Pixel(x,y)+img2->Pixel(x,y)) *
		 (img1->Pixel(x,y)+img2->Pixel(x,y)));
        
      // normalize each pixel to grayscale
      double red = Pixel(x,y).Red();
      double green = Pixel(x,y).Green();
      double blue = Pixel(x,y).Blue();
      double alpha = Pixel(x,y).Alpha();
        
      Pixel(x,y).SetRed(red + 0.5);
      Pixel(x,y).SetGreen(green + 0.5);
      Pixel(x,y).SetBlue(blue + 0.5);
      Pixel(x,y).SetAlpha(alpha + 0.5);
      Pixel(x,y).Clamp();
    }
  }

  delete img1;
  delete img2;
  delete img3;

  std::vector< R2Point* > top_scores = GetBestFeatures();
 
  // highlight 150 greatest threshold areas
  R2Point * feat = top_scores.at(0);
  for (int i = 0; i < top_scores.size(); i++) { 
    feat = top_scores.at(i);
    DrawBox(feat->X(), feat->Y(), true);
  }
}


void R2Image::
Sharpen()
{
  // create temporary image
  R2Image * tmp = new R2Image(this->width, this->height);

  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  int H[3][3] = {{0, -1, 0},
		 {-1, 5, -1},
		 {0, -1, 0}};

  for (int y=1; y<height-1; y++) {
    for (int x=1; x<width-1; x++) {
      for (int ly=-1; ly<2; ly++) {
        for (int lx=-1; lx<2; lx++) {
          tmp->Pixel(x,y) += Pixel(x+lx,y+ly)*H[ly+1][lx+1];
        }
      }
      tmp->Pixel(x,y).Clamp();
    }
  }

  *this = *tmp;
}


void R2Image::
HighPassFilter(double sigma)
{
  // create and blur temporary image
  R2Image * tmp = new R2Image(*this);
  tmp->Blur(sigma);

  // subtract pixels of temp image from original
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      Pixel(x,y) *= 2;
      Pixel(x,y) -= tmp->Pixel(x,y);
      Pixel(x,y).Clamp();
    }
  }
}


void R2Image::
HighPassSharpen(double sigma)
{
  R2Image * h = new R2Image(*this);
  h->HighPassFilter(sigma);
  R2Image * l = new R2Image(*this);
  l->Blur(sigma);

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      Pixel(x,y) = 0.5 * h->Pixel(x,y) + 0.5 * l->Pixel(x,y);
    }
  }
}


void R2Image::
MedianFilter(int window)
{
  R2Image * tmp = new R2Image(*this);
  
  // iterate over image
  for (int y=0; y<height; y++) {
    for (int x=0; x<height; x++) {

      int ind=0;
      std::vector<R2Pixel*> pixels;
      
      // iterate over search window, adding pixel intensities
      for (int ly=y-window; ly<y+window; ly++) {
	for (int lx=x-window; lx<x+window; lx++) {
	  if ((ly>=0 && ly<height) && (lx>=0 && lx<width)) {
	    pixels.push_back(new R2Pixel(Pixel(lx,ly)));
	    ind++;
	  }
	}
      }
      
      // sort pixel intensities
      std::sort(pixels.begin(), pixels.end(), pixelComp());
      const R2Pixel * p = pixels.at(ind/2);
      tmp->SetPixel(x, y, *p);
      delete p;
    }
  }

  *this = *tmp;
  delete tmp;
}


void R2Image::
BilateralFilter(double sigma, int window, double threshold)
{
  // Bilateral filter: gausian blur with threshold implemented
  const double PI = 3.1415926535;
  const double EULER = 2.7182818285;
  
  // instantiate and resize kernel
  std::vector<double> kernel;
  int s = int(sigma);
  int dim = (s * 6) + 1;
  kernel.resize(dim);

  // compute kernel
  for (int i=0; i<dim; i++) {
    kernel[i] = pow(EULER, -1 * pow(i-(3*s),2) / (2.*pow(sigma,2))) 
      / (sqrt(2. * PI) * sigma);
  }

  R2Image * tmp = new R2Image(this->width, this->height);

  // apply kernel to image over y direction
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      R2Pixel * val = new R2Pixel();
      double weights = 0;
      for (int ly=(-3*s); ly<=(3*s); ly++) {
        if (y+ly>=0 && y+ly<height) {
	  *val += Pixel(x,y+ly)*kernel[ly+3*s];
	  weights += kernel[ly+3*s];
        }
      }
      *val /= weights;
      tmp->Pixel(x,y) = *val;
      // tempImage.Pixel(x,y).Clamp();
    }
  }

  // apply kernel to image over x direction
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      R2Pixel * val = new R2Pixel();
      double weights = 0;
      for (int lx=(-3*s); lx<=(3*s); lx++) {
        if (x+lx>=0 && x+lx<width) {
	  *val += tmp->Pixel(x+lx,y)*kernel[lx+3*s];
	  weights += kernel[lx+3*s];
        }
      }
      *val /= weights;
      Pixel(x,y) = *val;
      // Pixel(x,y).Clamp();
    }
  }
}


void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1)
    {
      int x=y1;
      y1=y0;
      y0=x;

      x=x1;
      x1=x0;
      x0=x;
    }
  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;
  if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
  // note that this division needs to be done in a way that preserves the fractional part
  int y = y0;
  for(int x=x0;x<=x1;x++)
    {
      Pixel(x,y).Reset(r,g,b,1.0);
      error = error + deltaerr;
      if(error>=0.5)
	{
	  if(deltay>0) y = y + 1;
	  else y = y - 1;

	  error = error - 1.0;
	}
    }
  if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
    {
      for(int x=x0-3;x<=x0+3;x++)
	{
	  for(int y=y0-3;y<=y0+3;y++)
	    {
	      Pixel(x,y).Reset(r,g,b,1.0);
	    }
	}
    }
}


double R2Image::
SSD(double x1, double y1, R2Image * otherImage, double x2, double y2, double dx, double dy) 
{
  double sum = 0.;
  for (int i=-1*dx; i<=dx; i++) {
    for (int j=-1*dy; j<=dy; j++) {
      sum +=
        pow(getPixelMagnitude(x1+i,y1+j) - otherImage->getPixelMagnitude(x2+i,y2+j), 2);
    }
  }
  
  return sum;
}


double 
vecDiffMagnitude(std::pair<R2Point*, R2Point*> pair0,
		 std::pair<R2Point*, R2Point*> pair1)
{
  int dx1, dx2, dy1, dy2;

  // normalize pairs to 0,0 base
  dx1 = pair0.first->X() - pair0.second->X();
  dx2 = pair1.first->X() - pair1.second->X();
  dy1 = pair0.first->Y() - pair0.second->Y();
  dy2 = pair1.first->Y() - pair1.second->Y();

  return sqrt((dx1 - dx2)*(dx1 - dx2) + (dy1 - dy2)*(dy1 - dy2));
}  


std::vector< std::pair<R2Point*, R2Point*> > R2Image::
computeFeaturePairs(R2Image * otherImage) {
  R2Image * tmp = new R2Image(this->width, this->height, this->pixels);
  tmp->Harris(0.04);

  int searchWindowX = width / 10;
  int searchWindowY = height / 10;

  std::vector< std::pair< R2Point*, R2Point* > > pairs;
  std::vector< R2Point* > features = tmp->GetBestFeatures();
  delete tmp;
  
  // iterate over features
  // get surrounding pixels in original image
  // compare surrounding pixels to 20% area within other image
  // match using sum of squared differences of pixel values
  for (int i=0; i<features.size(); i++) {
    R2Point * pt = features.at(i);
    int minX = fmax(pt->X() - searchWindowX, 5),
      maxX = fmin(pt->X() + searchWindowX, width-5),
      minY = fmax(pt->Y() - searchWindowY, 5),
      maxY = fmin(pt->Y() + searchWindowY, height-5);
    double minDiff = std::numeric_limits<double>::infinity(), diff;
    int cX = 0, cY = 0;
  
    for(int ly=minY; ly<maxY; ly++) {
      for (int lx=minX; lx<maxX; lx++) {
        diff = SSD(pt->X(), pt->Y(), otherImage, lx, ly, 5, 5);
        if (diff < minDiff) { 
          minDiff = diff;
          cX = lx; 
          cY = ly;
        }
      }
    }
    R2Point *pt0 = new R2Point(pt->X(), pt->Y()), *pt1 = new R2Point (cX, cY);
    pairs.push_back(std::make_pair(pt0, pt1));
  }

  return pairs;
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
  // into this image with a 50% opacity.
  std::vector< std::pair < R2Point*, R2Point* > > pairs =
    computeFeaturePairs(otherImage);


  // run RANSAC to determine which pairs are good
  int ITERATIONS = 1000;
  int THRESHOLD  = 4;
  int bestVector = 0;
  int mostMatches = 0;

  for (int i=0; i<ITERATIONS; i++) {
    int r = rand() % 150;
    int matchesToThis = 0;
    std::pair<R2Point*, R2Point*>
      p0(pairs.at(r).first, pairs.at(r).second);
    for (int j=0; j<pairs.size(); j++) {
      std::pair<R2Point*, R2Point*>
	p1(pairs.at(j).first, pairs.at(j).second);
      if (vecDiffMagnitude(p0,p1) < THRESHOLD) {
	matchesToThis++;
      }
    }
    if (matchesToThis > mostMatches) {
      mostMatches = matchesToThis;
      bestVector = r;
    }
  }
  
  // draw lines matching features
  std::pair<R2Point*, R2Point*>
    pr0(pairs.at(bestVector).first, pairs.at(bestVector).second);
  for (int i=0; i<pairs.size(); i++) {
    std::pair<R2Point*, R2Point*>
      pr1(pairs.at(i).first, pairs.at(i).second);
    if (vecDiffMagnitude(pr0,pr1) < THRESHOLD) {
      DrawBox(pr0.first->X(), pr0.first->Y(), true);
      line(pr1.first->X(), pr1.second->X(),
	   pr1.first->Y(), pr1.second->Y(), 0, 1, 0);
    } else {
      DrawBox(pr1.first->X(), pr1.first->Y(), false);
      line(pr1.first->X(), pr1.second->X(),
	   pr1.first->Y(), pr1.second->Y(), 1, 0, 0);
    }
  }
  
  return;
}


R2Point* R2Image::
applyTransformationMatrix(R2Point * p, double * H)
{
  double x =
    H[0] * p->X() +
    H[1] * p->Y() +
    H[2] * 1;
  double y =
    H[3] * p->X() +
    H[4] * p->Y() +
    H[5] * 1;

  return new R2Point(x,y);
}


// gets determinant of 2x2 matrix
double R2Image::
twoDeterminant(double m[4])
{
  return m[0]*m[3] - m[1]*m[2];
}


// gets determinant of 3x3 matrix
double R2Image::
threeDeterminant(double m[9])
{
  return
    m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] -
    m[0]*m[5]*m[7] - m[1]*m[3]*m[8] - m[2]*m[4]*m[6];
}


R2Pixel * R2Image::
interpolate(double w, double h)
{
  double xRatio = w - floor(w);
  double yRatio = h - floor(h);
  
  int xMin = (w < 0)            ? 0 :            static_cast<int>(floor(w));
  int xMax = (w > this->width)  ? this->width :  static_cast<int>(ceil(w));
  int yMin = (h < 0)            ? 0 :            static_cast<int>(floor(h));
  int yMax = (h > this->height) ? this->height : static_cast<int>(ceil(h));

  R2Pixel *top, *bottom;
  bottom = new R2Pixel((1 - xRatio) * Pixel(xMin, yMin)
		       + xRatio * Pixel(xMax, yMin));
  top = new R2Pixel((1 - xRatio) * Pixel(xMin, yMax)
		    + xRatio * Pixel(xMax, yMax));
  
  return new R2Pixel((1 - yRatio) * *bottom + yRatio * *top);
}


void R2Image::
MergePixels(int h_s, int h_f, double* H,
	    R2Image * otherImage, R2Image * outputImage) {

  R2Point *t;
  for (int h=h_s; h<h_f; h++) {
    for (int w=0; w<outputImage->Width(); w++) {
      
      t = new R2Point(w,h);
      t = applyTransformationMatrix(t, H);

      if (w < width && h < height) {
	outputImage->Pixel(w, h) += 0.5 * otherImage->Pixel(w,h);
      }

      if (ceil(t->X()) >= 0 && ceil(t->X()) < width &&
	  ceil(t->Y()) >= 0 && ceil(t->Y()) < height) {
	outputImage->Pixel(w, h) +=
	  //0.5 * Pixel(static_cast<int>(ceil(t->X())),
	  //static_cast<int>(ceil(t->Y())));
	  0.5 * *interpolate(t->X(), t->Y());
      }
    }
  }
  delete t;
}


double * R2Image::
BuildH(std::vector< std::pair< R2Point*, R2Point* > > cor)
{
  // for each point, compute the first two rows of A_i
  int r;
  double** A = dmatrix(1,2*cor.size(),1,9);
  for (int i=0; i<cor.size(); i++) {
    r=i*2+1;
    A[r][1] = 0;
    A[r][2] = 0;
    A[r][3] = 0;
    A[r][4] = -1 * cor.at(i).first->X();
    A[r][5] = -1 * cor.at(i).first->Y();
    A[r][6] = -1;
    A[r][7] = cor.at(i).second->Y() * cor.at(i).first->X();
    A[r][8] = cor.at(i).second->Y() * cor.at(i).first->Y();
    A[r][9] = cor.at(i).second->Y();

    r++;
    A[r][1] = cor.at(i).first->X();
    A[r][2] = cor.at(i).first->Y();
    A[r][3] = 1;
    A[r][4] = 0;
    A[r][5] = 0;
    A[r][6] = 0;
    A[r][7] = -1 * cor.at(i).second->X() * cor.at(i).first->X();
    A[r][8] = -1 * cor.at(i).second->X() * cor.at(i).first->Y();
    A[r][9] = -1 * cor.at(i).second->X();
  }
    
  // Obtain the SVD, use the singular vector corresponding to the smallest
  // singular value
  double singularValues[10]; // 1..10
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(A, 2*cor.size(), 9, singularValues, nullspaceMatrix);

  // get the result
  int minIndex=1;
  for (int i=2; i<10; i++) {
    if (singularValues[i]<singularValues[minIndex]) minIndex = i;
  }

  double * H = (double*) malloc(sizeof(double) * 9);
  for(int i=0; i<9; i++) {
    H[i] = (1 / nullspaceMatrix[9][minIndex]) *
      nullspaceMatrix[i+1][minIndex];
  }
  
  return H;
}


double * R2Image::
ComputeHomographyMatrix(R2Image *otherImage)
{
  std::vector< std::pair< R2Point*, R2Point* > > pairs =
    computeFeaturePairs(otherImage);

  // run RANSAC to determine which pairs are good
  int ITERATIONS = 100;
  double * bestH = (double*) malloc(10 * sizeof(double));
  
  int THRESHOLD  = 4;
  int mostMatches = 0;

  double* H;
  for (int j=0; j<ITERATIONS; j++) {
    
    // DLT to compute H matrix
    // randomly select 4 point coordinances
    std::vector< std::pair<R2Point*, R2Point*> > cor;
    for (int i=0; i<4; i++) {
      int r = rand() % pairs.size();
      cor.push_back(pairs.at(r));
    }
        
    H = BuildH(cor);
    
    // apply H to each point and compute difference between points
    int matchesToThis = 0;
    for (int i=0; i<pairs.size(); i++) {
      std::pair<R2Point*, R2Point*>
	p(pairs.at(i).first, pairs.at(i).second);

      double x = H[0] * p.first->X() + H[1] * p.first->Y() + H[2];
      double y = H[3] * p.first->X() + H[4] * p.first->Y() + H[5];

      if (abs(x - p.second->X()) < THRESHOLD &&
	  abs(y - p.second->Y()) < THRESHOLD) {
	matchesToThis++;
      }

      if (matchesToThis > mostMatches) {
	mostMatches = matchesToThis;
	for (int i=0; i < 10; i++) {
	  bestH[i] = H[i];
	}
      }
    }
  }
  return bestH;
}


double * R2Image::
InvertHomographyMatrix(double* H)
{
  double det = threeDeterminant(H);
  double * ret = (double*)malloc(9 * sizeof(double));

  double tmp0[4] = {H[4], H[5], H[7], H[8]};
  ret[0] = (1/ det) * twoDeterminant(tmp0);
  double tmp1[4] = {H[2], H[1], H[8], H[7]};
  ret[1] = (1/ det) * twoDeterminant(tmp1);
  double tmp2[4] = {H[1], H[2], H[4], H[5]};
  ret[2] = (1/ det) * twoDeterminant(tmp2);
  double tmp3[4] = {H[5], H[3], H[8], H[6]};
  ret[3] = (1/ det) * twoDeterminant(tmp3);
  double tmp4[4] = {H[0], H[2], H[6], H[8]};
  ret[4] = (1/ det) * twoDeterminant(tmp4);
  double tmp5[4] = {H[2], H[0], H[5], H[3]};
  ret[5] = (1/ det) * twoDeterminant(tmp5);
  double tmp6[4] = {H[3], H[4], H[6], H[7]};
  ret[6] = (1/ det) * twoDeterminant(tmp6);
  double tmp7[4] = {H[1], H[0], H[7], H[6]};
  ret[7] = (1/ det) * twoDeterminant(tmp7);
  double tmp8[4] = {H[0], H[1], H[3], H[4]};
  ret[8] = (1/ det) * twoDeterminant(tmp8);

  return ret;
}


void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.

  double * H = ComputeHomographyMatrix(otherImage);
  // invert H

  H = InvertHomographyMatrix(H);
  // apply transformation matrix to corners of original image to
  //   determine bounds of transformed image
  R2Point * tl = new R2Point(0,0);
  R2Point * tr = new R2Point(0,this->height-1);
  R2Point * bl = new R2Point(this->width-1,0);
  R2Point * br = new R2Point(this->width-1, this->height-1);
  tl = applyTransformationMatrix(tl, H);
  tr = applyTransformationMatrix(tr, H);
  bl = applyTransformationMatrix(bl, H);
  br = applyTransformationMatrix(br, H);
  
  double xCoords[3] = {tr->X(), bl->X(), br->X()};
  double yCoords[3] = {tr->Y(), bl->Y(), br->Y()};
  double minX = tl->X(), minY = tl->Y();
  double maxX = tl->X(), maxY = tl->Y();

  delete tl;
  delete tr;
  delete bl;
  delete br;
  
  for (int i=0; i<3; i++) {
    if (xCoords[i] < minX) minX = xCoords[i];
    if (xCoords[i] > maxX) maxX = xCoords[i];
    if (yCoords[i] < minY) minY = yCoords[i];
    if (yCoords[i] > maxY) maxY = yCoords[i];
  }

  int nWidth = fmax(width, ceil(maxX - minX)),
    nHeight = fmax(height, ceil(maxY - minY));

  // apply transformation to each pixel in the original image
  R2Image * tmp = new R2Image(nWidth, nHeight);
  MergePixels(0, tmp->Height(), H, otherImage, tmp);
  
  // multithreading variables
  // const int nThreads = 4;
  // std::vector<std::thread> t;

  // // initialize threads
  // for (unsigned int i=0; i<height; i += height/nThreads) {
  //   t.push_back(std::thread(mergePixels, i, i+height/nThreads,
  //          minX, minY, H, otherImage, tmp));
  // }

  // // finish threads
  // for (unsigned int i=0; i<t.size(); i++) {
  //   t[i].join();
  // }
    
  *this = *tmp;
  delete tmp;
  
  return;
}


double R2Image::
greenRatio(double x, double y)
{
  R2Pixel p = Pixel(x,y);
  return p.Green() / (p.Red() + p.Green() + p.Blue());
}


R2Point* R2Image::
Convolve(R2Image * subImage, double x, double y, double dx, double dy, bool t)
{
  //if(t) printf("x = %f, y = %f\n", x, y);

  int sW = subImage->Width(), sH = subImage->Height();
  int lower_x = fmax(x-dx, sW/2), upper_x = fmin(x + dx, Width() - sW/2);
  int lower_y = fmax(y-dy, sH/2), upper_y = fmin(y + dy, Height() - sH/2);

  double fx=-1, fy=-1;
  double minDiff = std::numeric_limits<double>::infinity(), diff;

  // iterate over a search window to determine the center location of the subimage
  for(int ly=lower_y; ly<upper_y; ly++) {
    for(int lx=lower_x; lx<upper_x; lx++) {
      if (t) printf("lx = %d, ly = %d\n", lx, ly);
      diff = SSD(lx, ly, subImage, sW/2, sH/2, sW/2, sH/2);
      if (diff < minDiff) {
	minDiff = diff;
	fx = lx;
	fy = ly;
      }
      if (t) printf("diff = %f\n", diff);
    }
  }

  if (t) printf("upper_x = %d, upper_y = %d\n", upper_x, upper_y);

  if (fx < 0 || fy < 0) return new R2Point(x,y);
  return new R2Point(fx,fy);
}


std::vector< R2Point* > R2Image::
TrackMarkers(R2Image * marker1, R2Image * marker2, R2Image * marker3, R2Image * marker4)
{
  // get coordinates of each marker
  std::vector< R2Point* > markerCoords;
  markerCoords.resize(4);

  int w = Width()/4, h = Height()/4;

  // 2 each marker, convolve over image to determine location
  //assumes markers are originally in their respective quadrants of the image

  /*
  markerCoords.at(0)=Convolve(marker1, w, h, w, h, true);
  markerCoords.at(1)=Convolve(marker2, 3*w, h, w, h, false);
  markerCoords.at(2)=Convolve(marker3, w, 3*h, w, h, false);
  markerCoords.at(3)=Convolve(marker4, 3*w, 3*h, w, h, false);
  */

  markerCoords.at(0) = new R2Point(48, 129);
  markerCoords.at(1) = new R2Point(440, 130);
  markerCoords.at(2) = new R2Point(40, 285);
  markerCoords.at(3) = new R2Point(443, 286);

 
  return markerCoords;
}


std::vector< R2Point* > R2Image::
TrackMarkerMovement(R2Image * marker1, R2Image * marker2,
		    R2Image * marker3, R2Image * marker4,
		    std::vector< R2Point* > markers)
{
  int SEARCHWINDOW = 10;
  std::vector< R2Point* > ret;
  ret.resize(4);

  R2Point *m1 = markers.at(0), *m2 = markers.at(1), *m3 = markers.at(2), *m4 = markers.at(3);

  printf("m1 x = %f, y = %f\n", m1->X(), m1->Y());
  
  ret.at(0)=Convolve(marker1, m1->X(), m1->Y(), SEARCHWINDOW, SEARCHWINDOW, true);
  ret.at(1)=Convolve(marker2, m2->X(), m2->Y(), SEARCHWINDOW, SEARCHWINDOW, false);
  ret.at(2)=Convolve(marker3, m3->X(), m3->Y(), SEARCHWINDOW, SEARCHWINDOW, false);
  ret.at(3)=Convolve(marker4, m4->X(), m4->Y(), SEARCHWINDOW, SEARCHWINDOW, false);

  return ret;
}


R2Image* R2Image::
GetSubImage(R2Point* coord, double w, double h)
{
  R2Image * ret = new R2Image(w,h);
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      ret->Pixel(x,y) = Pixel(coord->X() - w/2 + x, coord->Y() - h/2 + y);
    }
  }
  return ret;
}


void R2Image::
ResizeImage(int w, int h)
{
  R2Image * tmp = new R2Image(w,h);
  int iw = Width(), ih = Height();
  for (int i = iw/2 - w/2; i < iw/2 + w/2; i++) {
    for (int j = ih/2 - h/2; j < ih/2 + h/2; j++) {
      tmp->Pixel(i - (iw/2 - w/2), j - (ih/2 - h/2)) = Pixel(i, j);
    }
  }
  *this = *tmp;
  delete tmp;
}


void R2Image::
LabelPoints(std::vector< R2Point* > points)
{
  R2Point* pt;
  for (int i=0; i<points.size(); i++) {
    pt = points.at(i);
    DrawBox(pt->X(), pt->Y(), false);
  }
}


void R2Image::
ProjectImage(R2Image * otherImage,
	     R2Image * m1, R2Image * m2, R2Image * m3, R2Image * m4)
{
  // normalize each of the marker subimages
  int SUBIMAGE_WIDTH = 32, SUBIMAGE_HEIGHT = 32;
  m1->ResizeImage(SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
  m2->ResizeImage(SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
  m3->ResizeImage(SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
  m4->ResizeImage(SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
  // locate 4 markers in original image
  std::vector< R2Point* > markerCoords = TrackMarkers(m1, m2, m3, m4);
  LabelPoints(markerCoords);
  ProjectPixels(otherImage, markerCoords);
  
  
  // implement feature tracking for the rest of the images
  R2Image *frame;
  std::string in_path = "frames4/frame_",
    out_path = "frames_out/frame_";
  for (int i = 2; i < 40; i++) {
    std::string in = in_path, out = out_path;
    if (i < 10) {
      in.append("0");
      out.append("0");
    }
    in.append(std::to_string(i));
    in.append(".jpeg\0");
    out.append(std::to_string(i));
    out.append(".jpeg\0");

    frame = new R2Image(in.c_str());

    markerCoords = frame->TrackMarkerMovement(m1, m2, m3, m4, markerCoords);

    m1 = frame->GetSubImage(markerCoords.at(0), SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
    m2 = frame->GetSubImage(markerCoords.at(1), SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
    m3 = frame->GetSubImage(markerCoords.at(2), SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);
    m4 = frame->GetSubImage(markerCoords.at(3), SUBIMAGE_WIDTH, SUBIMAGE_HEIGHT);

    // TODO: determine wtf is up with labelpoints method
    frame->LabelPoints(markerCoords);

    frame->ProjectPixels(otherImage, markerCoords);

    // Write output image
    if (!frame->Write(out.c_str())) {
      fprintf(stderr, "Unable to read image from %s\n", out.c_str());
      exit(-1);
    }
  }
}


// applies homogrphy transformation from otherImage onto original image
// markerCoords is a 4-vector
void R2Image::
ProjectPixels(R2Image* otherImage, std::vector< R2Point* > markerCoords)
{
  // compute correlation matrix and project pixels
  std::vector< std::pair< R2Point*, R2Point* > > cor;
  cor.resize(4);

  std::pair< R2Point*, R2Point* > p0 (new R2Point(0,0), markerCoords.at(0));
  std::pair< R2Point*, R2Point* > p1 (new R2Point(Width(),0), markerCoords.at(1));
  std::pair< R2Point*, R2Point* > p2 (new R2Point(0,Height()), markerCoords.at(2));
  std::pair< R2Point*, R2Point* > p3 (new R2Point(Width(),Height()), markerCoords.at(3));

  cor.at(0) = p0;
  cor.at(1) = p1;
  cor.at(2) = p2;
  cor.at(3) = p3;
  
  // compute homography matrix mapping points from full image to points within markers
  double * H = BuildH(cor);

  R2Point * p;
  for (int y = 0; y < otherImage->Height(); y++) {
    for (int x = 0; x < otherImage->Width(); x++) {
      p = applyTransformationMatrix(new R2Point(x,y), H);
      if (p->X() >= 0 & p->X() < width && p->Y() >= 0 && p->Y() < height &&
	greenRatio(p->X(), p->Y()) > 0.35) {
	Pixel(p->X(), p->Y()) = otherImage->Pixel(x,y);
      }
    }
  } 
}


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	double r = (double) getc(fp) / max_value;
	double g = (double) getc(fp) / max_value;
	double b = (double) getc(fp) / max_value;
	R2Pixel pixel(r, g, b, 1);
	SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
	R2Pixel pixel(r, g, b, 1);
	SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
	const R2Pixel& p = (*this)[i][j];
	int r = (int) (255 * p.Red());
	int g = (int) (255 * p.Green());
	int b = (int) (255 * p.Blue());
	fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
	if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
	const R2Pixel& p = (*this)[i][j];
	int r = (int) (255 * p.Red());
	int g = (int) (255 * p.Green());
	int b = (int) (255 * p.Blue());
	fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
};
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
	r = g = b = (double) *(p++) / 255;
	a = 1;
      }
      else if (ncomponents == 1) {
	r = g = b = (double) *(p++) / 255;
	a = 1;
	p++;
      }
      else if (ncomponents == 3) {
	r = (double) *(p++) / 255;
	g = (double) *(p++) / 255;
	b = (double) *(p++) / 255;
	a = 1;
      }
      else if (ncomponents == 4) {
	r = (double) *(p++) / 255;
	g = (double) *(p++) / 255;
	b = (double) *(p++) / 255;
	a = (double) *(p++) / 255;
      }
      else {
	fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
	return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}

  
int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width;  /* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;   /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}

