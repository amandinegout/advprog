// Threaded two-dimensional Discrete FFT transform
// Amandine GOUT
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <pthread.h>
#include "Complex.h"
#include "InputImage.h"


// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.


// Global variables

// Each thread needs to know how many threads there are
int nb_threads = 16;
Complex* hin;
// Image 
int im_width, im_height;
//Weights
Complex* weights; 
Complex* weights_inv; 

// The mutex and condition variables allow the main thread to
// know when all helper threads are completed.
pthread_mutex_t countMutex;
pthread_mutex_t exitMutex;
pthread_cond_t exitCond;
bool* localSense; // array of bools, one per thread
bool globalSense;
int startCount;

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v, int N)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init(){ // you will likely need some parameters
	// All mutex and condition variables must be "initialized"
	pthread_mutex_init(&countMutex, 0);
    pthread_mutex_init(&exitMutex,0);
    pthread_cond_init(&exitCond, 0);
    // Holds the exit mutex until waiting for exitCond condition
    pthread_mutex_lock(&exitMutex);
	startCount = nb_threads+1; // main thread
	// Create and initialize the localSense arrar, 1 entry per thread
	localSense = new bool[nb_threads+1];
	for (int i = 0; i < nb_threads+1; ++i) 
		localSense[i] = true;
	// Initialize global sense
	globalSense = true;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(unsigned long myId){ // Again likely need parameters 
	localSense[myId] = !localSense[myId]; // Toggle private sense variable
	pthread_mutex_lock(&countMutex);
	int myCount = startCount;
	startCount--;
	pthread_mutex_unlock(&countMutex);
	
	if (myCount == 1) { // All threads here, reset count and toggle global sense
		startCount = nb_threads+1; // reset of startCount
		globalSense = localSense[myId];
		pthread_mutex_lock(&exitMutex);
		pthread_cond_signal(&exitCond);
		pthread_mutex_unlock(&exitMutex);
	}else{
		while (globalSense != localSense[myId]) { } // Spin
	}
}

void ReverseInitialArray(Complex* h, int N, int start_index)
{
	Complex swapArray[N];
	unsigned n_opp = 0;

	for(int n=0;n<N;n++){
	  n_opp = ReverseBits(n,N);
	  swapArray[n] = *(h+start_index+n_opp);
	}
	for(int n=0;n<N;n++){
	  *(h+start_index+n) = swapArray[n];
	}
}

void ComputeWeights(Complex* W, int N)
{
  int size = N/2;
  for (int n = 0; n < size;++n){
	  double theta = (2*M_PI*n)/N;
	  W[n]= Complex(cos(theta),-1*sin(theta));
  }
}

void ComputeWeightsInv(Complex* W, int N)
{
  int size = N/2;
  for (int n = 0; n < size;++n){
	  double theta = (2*M_PI*n)/N;
	  W[n]= Complex(cos(theta),+1*sin(theta));
  }
}

void transpose(Complex* src, int wd, int he) {
	Complex temp;
	for (int i=0;i<he;++i){
		for (int j=i;j<wd;++j){
			temp = *(src+j*he+i);
			*(src+j*he+i) = *(src+i*wd+j);
			*(src+i*wd+j)=temp;
		}
	}
}

void Transform1D(Complex* h, int N, Complex* W, int start)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  // "W" is the array of weights to compute the FFT
  
  // Reorder row from intial array
  ReverseInitialArray(h,N,start*N);
  
  int nb_iter = log2(N);
  
  Complex temp1;
  Complex temp2;
  int weight_index;
  for (int k =1;k<nb_iter+1;k++){
	  int nb_transform = (int) (N/(pow(2,k)));
	  for (int i=0;i<nb_transform;i++){
		  int transform_start = i*pow(2,k);
		  for (int n=0;n<pow(2,k-1);n++){
			  weight_index = (int)(n*N/pow(2,k));
			  temp1 = h[(int) (start*N+transform_start+n)];
			  temp2 = W[weight_index].operator*(h[(int) (start*N+transform_start+n+pow(2,k-1))]);
			  h[(int) (start*N+transform_start+n)]= temp1.operator+(temp2);
			  h[(int) (start*N+transform_start+n+pow(2,k-1))]= temp1.operator-(temp2);
		  } 
	  }
  }
} 

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  unsigned long myId = (unsigned long)v;
  
  //cout << "Thread number " << myId << " starting 1D transform" << endl;
  // Calculate 1d DFT for assigned rows
  int rowsPerThread = im_height / nb_threads;
  
  int startingRow = myId * rowsPerThread;
  for (int r = 0; r < rowsPerThread; ++r){
	  int thisRow = startingRow + r;
	  Transform1D(hin, im_width, weights, thisRow);
  }
  //cout << "Thread number " << myId << " completed 1D transform" << endl;
  //cout << "Thread number " << myId << " entering barrier" << endl;
  MyBarrier(myId);
  
  return 0;
}

void Transform2D(const char* inputFN) 
{ 
  // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  im_width = image.GetWidth();
  im_height = image.GetHeight();
  // Create image array data
  //cout << "Loading image data" <<endl;
  hin = image.GetImageData();
  
  // Compute weights for FFT
  //cout << "Computing weights for FFT" <<endl;
  weights = new Complex[im_width/2];
  ComputeWeights(weights, im_width);
  
  // Create 16 threads and start 1D FFT
  
  //cout << "Starting 1D transform " << endl;
  // Now start the threads
  for (int i = 0; i < nb_threads; ++i){
	// Now create the thread
	pthread_t pt; // pThread variable (output param from create)
	// Third param is the thread starting function
	// Fourth param is passed to the thread starting function
	pthread_create(&pt, 0, Transform2DTHread, (void*)i);
  }
  // Wait for all threads complete
  MyBarrier(nb_threads);
  //cout << "1D transform complete on rows" << endl;
  
  //cout << "Writing image data for 1D" << endl;
  // Write the transformed data
  char const* filename = "MyAfter1D.txt";
  image.SaveImageData(filename, hin, im_width, im_height);
  
  //cout << "Starting 1D transform on columns" << endl;
  // Do TRANSPOSE
  transpose(hin, im_width, im_height);
  // Now start the threads
  for (int i = 0; i < nb_threads; ++i){
	// Now create the thread
	pthread_t pt; // pThread variable (output param from create)
	pthread_create(&pt, 0, Transform2DTHread, (void*)i);
  }
  // Wait for all threads complete
  MyBarrier(nb_threads);
  // Do TRANSPOSE
  transpose(hin, im_height, im_width);
  //cout << "1D transform complete on columns " << endl;
  //cout << "2D transform complete" << endl;
  
  //cout << "Writing image data for 2D" << endl;
  // Write the transformed data
  char const* filename2 = "MyAfter2D.txt";
  image.SaveImageData(filename2, hin, im_width, im_height);
  
  /********************************************************************/
  // Compute weights for FFT Inverse
  //cout << "Computing weights for FFT Inverse" <<endl;
  ComputeWeightsInv(weights, im_width);
  
  //cout << "Starting 1D inverse transform " << endl;
  // Now start the threads
  for (int i = 0; i < nb_threads; ++i){
	// Now create the thread
	pthread_t pt; // pThread variable (output param from create)
	pthread_create(&pt, 0, Transform2DTHread, (void*)i);
  }
  // Wait for all threads complete
  MyBarrier(nb_threads);
  //cout << "1D inverse transform complete on rows" << endl;
  
  //cout << "Starting 1D inverse transform on columns" << endl;
  // Do TRANSPOSE
  transpose(hin, im_width, im_height);
  // Now start the threads
  for (int i = 0; i < nb_threads; ++i){
	// Now create the thread
	pthread_t pt; // pThread variable (output param from create)
	pthread_create(&pt, 0, Transform2DTHread, (void*)i);
  }
  // Wait for all threads complete
  MyBarrier(nb_threads);
  // Do TRANSPOSE
  transpose(hin, im_height, im_width);
  //cout << "1D inverse transform complete on columns " << endl;
  //cout << "2D inverse transform complete" << endl;
  
  // Normalize inverse transform
  Complex norm = Complex(1/((float)im_height*(float)im_width));
  for (int k=0;k<im_width*im_height;++k){
	hin[k] = norm.operator*(hin[k]);
  }
  
  //cout << "Writing image data for Inverse 2D" << endl;
  // Write the transformed data
  char const* filename3 = "MyAfterInverse.txt";
  image.SaveImageData(filename3, hin, im_width, im_height);
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  MyBarrier_Init();
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
