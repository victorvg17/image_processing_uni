/**
 * Image Processing and Computer Graphics
 *
 * Exercise 1: Noise, basic operators and filters
 */

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include <math.h> ///log, cos
#include "CMatrix.h"
#include "CTensor.h"
#include <assert.h>

// #define GAUSSIAN_NOISE
// #define IMAGE_DIFF

using namespace std;

typedef vector<double> Array;
typedef vector<Array> Matrix;

Matrix getGaussian(int height, int width, double sigma) {
  Matrix kernel(height, Array(width));
  double sum=0.0;
  int i,j;
  double exp_power;

  double den_r = 2*M_PI*sigma*sigma;

  for (i=0 ; i<height ; i++) {
    for (j=0 ; j<width ; j++) {
      exp_power = (i*i+j*j)/(2*sigma*sigma);
      kernel[i][j] = exp(-exp_power)/(den_r);
      sum += kernel[i][j];
    }
  }

  for (i=0 ; i<height ; i++) {
    for (j=0 ; j<width ; j++) {
      kernel[i][j] /= sum;
    }
  }

  return kernel;
}

CMatrix<float> applyGaussianFilter(CMatrix<float>& sourceImage, Matrix& filter){
  
  int height = sourceImage.ySize();
  int width = sourceImage.xSize();
  int filterHeight = filter.size();
  int filterWidth = filter[0].size();
  int destImageHeight = height - filterHeight + 1;
  int destImageWidth = width - filterWidth + 1;
  int i,j,h,w;

  cout << "filterHeight, filterwidth: " << height << " " << width << endl;

  CMatrix<float> destImage(destImageHeight, destImageWidth);

  
  for (i=1 ; i<destImageHeight ; i++) {
    for (j=1 ; j<destImageWidth ; j++) {
      for (h=i-1 ; h<i+filterHeight ; h++) {
        for (w=j-1 ; w<j+filterWidth ; w++) {
            destImage(i, j) += filter[h-i][w-j]*sourceImage(h, w);
        }
      }
    }
  }

  return destImage;

}

//generates gaussian noise using box-muller method as given in slides
const double box_muller_method(const int sigma_1, const int sigma_2) {
  const double u = (double) rand() / RAND_MAX;
  const double v = (double) rand() / RAND_MAX;

  const double n = sqrt(-2*log(u)) * cos(2*M_PI*v);
  const double m = sqrt(-2*log(u)) * sin(2*M_PI*v);
  const double tot_noise = sigma_1*n + sigma_2*m;
  return tot_noise;
}

int main(int argc, char** args) 
{  
  /// Tell the compiler not to throw warnings for unused variables
  /// Remove these lines if you want to use command line arguments.
  (void)argc;
  (void)args;

  /// Print something so we know the program actually runs
  std::cout << "Hello, World!\n";

  #ifdef GAUSSIAN_NOISE
    /// Define image
    CMatrix<float> aImage;
    /// Read image from a PGM file
    aImage.readFromPGM("lena.pgm");

    const int sigma_1 = 10;
    const int sigma_2 = 20;

    const double psnr_numr = pow(aImage.max() - aImage.min(), 2);
    double psnr_denr = 0;
    /// Add Gaussian noise here
    for (int y = 0; y < aImage.ySize(); ++y) {
      for (int x = 0; x < aImage.xSize(); ++x) {
        const double tot_noise = box_muller_method(sigma_1, sigma_2);
        psnr_denr += pow(tot_noise, 2);
        aImage(x, y) += tot_noise;
      }
    }
    aImage.clip(0, 255);
    /// Write noisy image to PGM file
    aImage.writeToPGM("lenaNoisy.pgm");

    const double psnr = 10*log(psnr_numr/psnr_denr) / log(10);
    std::cout << "PSNR: " << psnr << std::endl;
  #endif

  #ifdef IMAGE_DIFF
    /// Define image
    CTensor<float> aImage_1;
    /// Read image from a PGM file
    aImage_1.readFromPPM("Sidenbladh.ppm");

    CTensor<float> aImage_2;
    /// Read image from a PGM file
    aImage_2.readFromPPM("SidenbladhBG.ppm");

    aImage_1 *= -1;
    aImage_2 += aImage_1;

    aImage_2.writeToPPM("Sidenbladhdiff.ppm");
  #endif

  Matrix filter = getGaussian(5, 5, 1.0);

  /// Define image
  CMatrix<float> sourceImage;
  /// Read image from a PGM file
  sourceImage.readFromPGM("chinaToilet.pgm");

  ///destination image
  CMatrix<float> destImage = applyGaussianFilter(sourceImage, filter);
  destImage.writeToPGM("chinaToiletGaussFilter.pgm");
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 5; ++j) 
        cout << filter[i][j] << "\t"; 
    cout << endl; 
  }




  return EXIT_SUCCESS;
}

