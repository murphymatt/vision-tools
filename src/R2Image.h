b// Include file for image class
#ifndef R2_IMAGE_INCLUDED
#define R2_IMAGE_INCLUDED

#include <vector>


// Constant definitions

typedef enum {
  R2_IMAGE_RED_CHANNEL,
  R2_IMAGE_GREEN_CHANNEL,
  R2_IMAGE_BLUE_CHANNEL,
  R2_IMAGE_ALPHA_CHANNEL,
  R2_IMAGE_NUM_CHANNELS
} R2ImageChannel;

typedef enum {
  R2_IMAGE_POINT_SAMPLING,
  R2_IMAGE_BILINEAR_SAMPLING,
  R2_IMAGE_GAUSSIAN_SAMPLING,
  R2_IMAGE_NUM_SAMPLING_METHODS
} R2ImageSamplingMethod;

typedef enum {
  R2_IMAGE_OVER_COMPOSITION,
  R2_IMAGE_IN_COMPOSITION,
  R2_IMAGE_OUT_COMPOSITION,
  R2_IMAGE_ATOP_COMPOSITION,
  R2_IMAGE_XOR_COMPOSITION,
} R2ImageCompositeOperation;


// Class definition

class R2Image {
 public:
  // Constructors/destructor
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height);
  R2Image(int width, int height, const R2Pixel *pixels);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Image properties
  int NPixels(void) const;
  int Width(void) const;
  int Height(void) const;

  // Pixel access/update
  R2Pixel& Pixel(int x, int y);
  R2Pixel *Pixels(void);
  R2Pixel *Pixels(int row);
  R2Pixel *operator[](int row);
  const R2Pixel *operator[](int row) const;
  void SetPixel(int x, int y,  const R2Pixel& pixel);

  // Image processing
  R2Image& operator=(const R2Image& image);

  // Per-pixel operations
  void Brighten(double factor);
  void ChangeSaturation(double factor);

  // show how SVD works
  void svdTest();


  // helper functions
  void DrawBox(int x, int y, bool good);
  double getPixelMagnitude(int x, int y);
  std::vector< R2Point* > GetBestFeatures(void);
  void line(int x0, int x1, int y0, int y1, float r, float g, float b);
  double SSD(double x0, double y0, R2Image * otherImage,
	     double x1, double y1, double dx, double dy);
  std::vector< std::pair <R2Point*, R2Point*> >
    computeFeaturePairs(R2Image* otherImage);
  R2Point* applyTransformationMatrix(R2Point* p, double* H);
  double* BuildH(std::vector< std::pair< R2Point*, R2Point* > > cor);
  double* ComputeHomographyMatrix(R2Image *otherImage);
  double* InvertHomographyMatrix(double* H);
  void MergePixels(int h_s, int h_f, double* H,
		   R2Image * otherImage, R2Image * outputImage);
  double twoDeterminant(double m[4]);
  double threeDeterminant(double m[9]);
  R2Pixel* interpolate(double width, double height);
  double greenRatio(double x, double y);
  R2Point* Convolve(R2Image * subImage, double x, double y, double dx, double dy);
  double CompareBlock(R2Image * subImage);
  std::vector< R2Point* >
    TrackMarkers(R2Image * marker1, R2Image * marker2, R2Image * marker3, R2Image * marker4);
  std::vector< R2Point* >
    TrackMarkerMovement(R2Image * marker1, R2Image * marker2,
			R2Image * marker3, R2Image * marker4,
			std::vector< R2Point* > markers);
  R2Image * GetSubImage(R2Point* coord, double w, double h);
  void ResizeImage(int w, int h);
  void ProjectPixels(R2Image* otherImage, std::vector< R2Point* > markerCoords);
  void LabelPoints(std::vector< R2Point* > points);
  

  // Linear filtering operations
  void SobelX();
  void SobelY();
  void LoG();
  void Blur(double sigma);
  void BlurD(double sigma);
  void Harris(double sigma);
  void Sharpen(void);
  void HighPassFilter(double sigma); 
  void HighPassSharpen(double sigma);


  // Non-Linear filtering operations
  void MedianFilter(int window);
  void BilateralFilter(double sigma, int window, double threshold);


  // further operations
  void blendOtherImageTranslated(R2Image * otherImage);
  void blendOtherImageHomography(R2Image * otherImage);
  void ProjectImage(R2Image * otherImage,
		    R2Image * m1, R2Image * m2, R2Image * m3, R2Image * m4);


  // File reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadJPEG(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteJPEG(const char *filename) const;

 private:
  // Utility functions
  void Resize(int width, int height);
  R2Pixel Sample(double u, double v,  int sampling_method);

 private:
  R2Pixel *pixels;
  int npixels;
  int width;
  int height;
};



// Inline functions

inline int R2Image::
NPixels(void) const
{
  // Return total number of pixels
  return npixels;
}



inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline R2Pixel& R2Image::
Pixel(int x, int y)
{
  // Return pixel value at (x,y)
  // (pixels start at lower-left and go in row-major order)
  return pixels[x*height + y];
}



inline R2Pixel *R2Image::
Pixels(void)
{
  // Return pointer to pixels for whole image 
  // (pixels start at lower-left and go in row-major order)
  return pixels;
}



inline R2Pixel *R2Image::
Pixels(int x)
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline R2Pixel *R2Image::
operator[](int x) 
{
  // Return pixels pointer for row at x
  return Pixels(x);
}



inline const R2Pixel *R2Image::
operator[](int x) const
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline void R2Image::
SetPixel(int x, int y, const R2Pixel& pixel)
{
  // Set pixel
  pixels[x*height + y] = pixel;
}



#endif
