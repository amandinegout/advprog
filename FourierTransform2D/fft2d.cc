// Distributed two-dimensional Discrete FFT transform
// AMANDINE GOUT
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "Complex.h"
#include "InputImage.h"
#include <math.h> 

using namespace std;

void Transform1D(Complex* h, int w, Complex* H,int rank_nb, int row_point);

void Transform1Dinverse(Complex* H, int w, Complex* h, int in_start, int row_point);

void transpose(Complex* src, Complex* dst, int wd, int he) {
	for (int i=0;i<he;++i){
		for (int j=0;j<wd;++j){
			dst[j*he+i]=src[i*wd+j];
		}
	}
}
   
void Transform2D(const char* inputFN,int argc, char** argv) 
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  // 9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().
  
  /* Step 1 */
  InputImage image(inputFN);  // Create the helper object for reading the image
  int im_width = image.GetWidth();
  int im_height = image.GetHeight();
  
  /* Step 2 */
  int  numtasks, rank, rc; 
  
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks); // Returns total number of available CPUs
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  /* Step 3 */
  Complex* Hout = new Complex[im_width*im_height];
  /* Step 4 */
  Complex* hin = image.GetImageData();
  
  /* Step 5 */
  int nb_rows_cpu = im_width/numtasks;
  int nb_cols_cpu = im_height/numtasks;
  int start_rows = im_width*nb_rows_cpu*rank;
  //int start_cols = im_height*nb_cols_cpu*rank;
  Complex* row_buff = new Complex[im_width*nb_rows_cpu];
  Complex* col_buff = new Complex[im_height*nb_cols_cpu];
  Complex* row_buff_res = new Complex[im_width*nb_rows_cpu];
  Complex* col_buff_res = new Complex[im_height*nb_cols_cpu];
  
  
  if (rank==0){
	cout << "Width is " << im_width << endl;
    cout << "Height is " << im_height << endl;
	printf ("Number of processors = %d My rank= %d\n", numtasks,rank);
	cout << "Each CPU has to compute " << nb_rows_cpu << " rows" << endl;
    
	for (int i=0;i<nb_rows_cpu;++i){
		Transform1D(hin, im_width, row_buff, start_rows, i*im_width);
	}
	/* Step 6 */
	// CPU 0 fills his rows in Hout
	for (int k=0;k<nb_rows_cpu*im_width;++k){
		Hout[k] = row_buff[k];
	}
	// CPU 0 receives from other CPUs and write their rows in Hout
	for(int round = 1;round<numtasks;++round){
		MPI_Status status; // Will be used to indicate the source of the message
		rc = MPI_Recv(row_buff, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " recv failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
		int count = 0;
		MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
		cout << "Rank " << rank 
		<< " received " << count << " bytes from "
		<< status.MPI_SOURCE << endl;
		for (int k=0;k<nb_rows_cpu*im_width;++k){
			Hout[status.MPI_SOURCE*im_width*nb_rows_cpu+k] = row_buff[k];
		} 
	}
	// Save text file with fft1d
	char const* filename = "MyAfter1D.txt";
	image.SaveImageData(filename, Hout, im_width, im_height);
	
	// CPU 0 transposes Hout
	Complex* Hout_trans = new Complex[im_height*im_width];
	transpose(Hout, Hout_trans, im_width, im_height);
	
	/* Step 6a */
	// CPU 0 sends the columns (now ligns) to other CPUs
	for(int rk = 1;rk<numtasks;++rk){
		for(int i = 0;i<im_height*nb_cols_cpu;++i){
			col_buff[i] = Hout_trans[im_height*nb_cols_cpu*rk+i] ;
		}
		rc = MPI_Send(col_buff, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, rk,
                        0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " send failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
	}
	
	// CPU 0 computes the 1d-fft of his columns
	for(int i = 0;i<im_height*nb_cols_cpu;++i){
		col_buff[i] = Hout_trans[i] ;
	}
    for (int i=0;i<nb_cols_cpu;++i){
		Transform1D(col_buff, im_height, col_buff_res, 0,i*im_height);
	}
	// CPU 0 puts his results from 1d_fft columns in Hout_trans
	for (int k=0;k<nb_cols_cpu*im_height;++k){
		Hout_trans[k] = col_buff_res[k];
	}
	// CPU 0 receives from other CPUs and write their results in Hout_trans
	for(int round = 1;round<numtasks;++round){
		MPI_Status status; // Will be used to indicate the source of the message
		rc = MPI_Recv(col_buff_res, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " recv failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
		int count = 0;
		MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
		cout << "Rank " << rank 
		<< " received " << count << " bytes from "
		<< status.MPI_SOURCE << endl;
		for (int k=0;k<nb_cols_cpu*im_height;++k){
			Hout_trans[status.MPI_SOURCE*im_height*nb_cols_cpu+k] = col_buff_res[k];
		}
	}
	// CPU 0 transposes Hout_trans to get Hout
	transpose(Hout_trans, Hout, im_height, im_width);
	
	// Save text file with fft1d
	char const* filename2 = "MyAfter2D.txt";
	image.SaveImageData(filename2, Hout, im_width, im_height);
    
    /***************************************/
    // Graduate students part : FFT inverse
    /***************************************/
    cout << "PERFORMING FFT INVERSE " << endl;
    
    Complex* hin_fftinv = new Complex[im_height*im_width];
    
    // CPU 0 sends the rows of Hout to other CPUs
	for(int rk = 1;rk<numtasks;++rk){
		for(int i = 0;i<im_width*nb_rows_cpu;++i){
			row_buff[i] = Hout[im_width*nb_rows_cpu*rk+i] ;
		}
		rc = MPI_Send(row_buff, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, rk,
                        0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " send failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
	}  
	
	// CPU 0 computes the 1d-fft inverse of his rows
	for(int i = 0;i<im_width*nb_rows_cpu;++i){
		row_buff[i] = Hout[i] ;
	}

    for (int i=0;i<nb_rows_cpu;++i){
		Transform1Dinverse(row_buff, im_width, row_buff_res, 0,i*im_width);
	}

	// CPU 0 fills his rows in hin_fftinv
	for (int k=0;k<nb_rows_cpu*im_width;++k){
		hin_fftinv[k] = row_buff_res[k];
	}
	
	// CPU 0 receives from other CPUs and write their rows in hin_fftinv
	for(int round = 1;round<numtasks;++round){
		MPI_Status status; // Will be used to indicate the source of the message
		rc = MPI_Recv(row_buff_res, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " recv failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
		int count = 0;
		MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
		cout << "Rank " << rank 
		<< " received " << count << " bytes from "
		<< status.MPI_SOURCE << endl;
		for (int k=0;k<nb_rows_cpu*im_width;++k){
			hin_fftinv[status.MPI_SOURCE*im_width*nb_rows_cpu+k] = row_buff_res[k];
		} 
	}
	
	// CPU 0 transposes hin_fftinv
	Complex* hin_fftinv_trans = new Complex[im_height*im_width];
	transpose(hin_fftinv, hin_fftinv_trans, im_width, im_height);
	
	// CPU 0 sends the columns (now rows) to other CPUs
	for(int rk = 1;rk<numtasks;++rk){
		for(int i = 0;i<im_height*nb_cols_cpu;++i){
			col_buff[i] = hin_fftinv_trans[im_height*nb_cols_cpu*rk+i] ;
		}
		rc = MPI_Send(col_buff, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, rk,
                        0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " send failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
	}
	
	// CPU 0 computes the 1d-fft inverse of his columns
	for(int i = 0;i<im_height*nb_cols_cpu;++i){
		col_buff[i] = hin_fftinv_trans[i] ;
	}
    for (int i=0;i<nb_cols_cpu;++i){
		Transform1Dinverse(col_buff, im_height, col_buff_res, 0,i*im_height);
	}
	
	// CPU 0 normalises the results and puts the 1d_fft inverse columns in hin_fftinv_trans
	Complex norm = Complex(1/((float)im_height*(float)im_width));
	for (int k=0;k<nb_cols_cpu*im_height;++k){
		hin_fftinv_trans[k] = norm.operator*(col_buff_res[k]);
	}
	
	// CPU 0 receives from other CPUs and write their results in hin_fftinv_trans
	for(int round = 1;round<numtasks;++round){
		MPI_Status status; // Will be used to indicate the source of the message
		rc = MPI_Recv(col_buff_res, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		if (rc != MPI_SUCCESS){
			cout << "Rank " << rank
			<< " recv failed, rc " << rc << endl;
			MPI_Finalize();
			exit(1);
		}
		int count = 0;
		MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
		cout << "Rank " << rank 
		<< " received " << count << " bytes from "
		<< status.MPI_SOURCE << endl;
		for (int k=0;k<nb_cols_cpu*im_height;++k){
			hin_fftinv_trans[status.MPI_SOURCE*im_height*nb_cols_cpu+k] = norm.operator*(col_buff_res[k]);
		}
	}
	
	// CPU 0 transposes hin_fftinv_trans to get hin_fftinv
	transpose(hin_fftinv_trans, hin_fftinv, im_height, im_width);
	
	// Save text file with fft inverse
	char const* filename3 = "MyAfterInverse.txt";
	image.SaveImageData(filename3, hin_fftinv, im_width, im_height);
    
  } else {
    for (int i=0;i<nb_rows_cpu;++i){
		Transform1D(hin, im_width, row_buff, start_rows,i*im_width);
	}

	rc = MPI_Send(row_buff, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, 0,
                        0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS){
        cout << "Rank " << rank
        << " send failed, rc " << rc << endl;
        MPI_Finalize();
        exit(1);
    }
    
    // Other CPUs receive their columns from CPU 0
	MPI_Status status; // Will be used to indicate the source of the message
	rc = MPI_Recv(col_buff, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, 0,0, MPI_COMM_WORLD, &status);
	if (rc != MPI_SUCCESS){
		cout << "Rank " << rank
		<< " recv failed, rc " << rc << endl;
		MPI_Finalize();
		exit(1);
	}
	int count = 0;
	MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
	cout << "Rank " << rank 
	<< " received " << count << " bytes from "
	<< status.MPI_SOURCE << endl;
	
	// CPUs compute their 1d-FFT for columns
    for (int i=0;i<nb_cols_cpu;++i){
		Transform1D(col_buff, im_height, col_buff_res, 0,i*im_height);
	}
	
	//CPUs send their results to CPU 0
	rc = MPI_Send(col_buff_res, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, 0,
                        0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS){
        cout << "Rank " << rank
        << " send failed, rc " << rc << endl;
        MPI_Finalize();
        exit(1);
    }
    
    /***************************************/
    // Graduate students part : FFT inverse
    /***************************************/
    
    // Other CPUs receive their rows from CPU 0
    
	rc = MPI_Recv(row_buff, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, 0,0, MPI_COMM_WORLD, &status);
	if (rc != MPI_SUCCESS){
		cout << "Rank " << rank
		<< " recv failed, rc " << rc << endl;
		MPI_Finalize();
		exit(1);
	}
	
	MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
	cout << "Rank " << rank 
	<< " received " << count << " bytes from "
	<< status.MPI_SOURCE << endl;
	
	// CPUs compute their 1d-FFT inverse for rows
    for (int i=0;i<nb_rows_cpu;++i){
		Transform1Dinverse(row_buff, im_width, row_buff_res, 0,i*im_width);
	}
    
    //CPUs send their results to CPU 0
	rc = MPI_Send(row_buff_res, sizeof(Complex)*im_width*nb_rows_cpu, MPI_CHAR, 0,
                        0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS){
        cout << "Rank " << rank
        << " send failed, rc " << rc << endl;
        MPI_Finalize();
        exit(1);
    }
    
    // Other CPUs receive their columns from CPU 0

	rc = MPI_Recv(col_buff, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, 0,0, MPI_COMM_WORLD, &status);
	if (rc != MPI_SUCCESS){
		cout << "Rank " << rank
		<< " recv failed, rc " << rc << endl;
		MPI_Finalize();
		exit(1);
	}

	MPI_Get_count(&status, MPI_CHAR, &count); // returns source, tag and number of element of datatype received
	cout << "Rank " << rank 
	<< " received " << count << " bytes from "
	<< status.MPI_SOURCE << endl;
    
    
    // CPUs compute their 1d-FFT inverse for columns
    for (int i=0;i<nb_cols_cpu;++i){
		Transform1Dinverse(col_buff, im_height, col_buff_res, 0,i*im_height);
	}
	
	//CPUs send their results to CPU 0
	rc = MPI_Send(col_buff_res, sizeof(Complex)*im_height*nb_cols_cpu, MPI_CHAR, 0,
                        0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS){
        cout << "Rank " << rank
        << " send failed, rc " << rc << endl;
        MPI_Finalize();
        exit(1);
    }
	
  }
  
    MPI_Finalize();
  
  
}

void Transform1D(Complex* h, int w, Complex* H, int in_start, int row_point)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  for (int n = 0; n < w;++n){
	  H[n+row_point] = Complex(0,0);
	  for (int k = 0; k < w;++k){
		  Complex Htmp = Complex(cos((2*M_PI*n*k)/w),-sin((2*M_PI*n*k)/w)).operator *(h[k+in_start+row_point]);
		  H[n+row_point] = H[n+row_point].operator +(Htmp);
	  }
  }
}

void Transform1Dinverse(Complex* H, int w, Complex* h, int in_start, int row_point)
{
  // Implement a simple 1-d DFT inverse using the double summation equation
  // given in the assignment handout.  h is the time-domain output
  // data, w is the width (N), and H is the FFT input array.
  for (int k = 0; k < w;++k){
	  h[k+row_point] = Complex(0,0);
	  for (int n = 0; n < w;++n){
		  Complex htmp = Complex(cos((2*M_PI*n*k)/w),sin((2*M_PI*n*k)/w)).operator *(H[n+in_start+row_point]);
		  h[k+row_point] = h[k+row_point].operator +(htmp);
	  }
  }
}

int main(int argc, char** argv)
{
  //string fn("Tower.txt"); // default file 
  string fn("Unit_test.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str(),argc,argv); // Perform the transform.
}  
  

  
