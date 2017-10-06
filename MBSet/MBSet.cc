// Calculate and display the Mandelbrot set
// ECE4893/8893 final project, Fall 2011

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <stack>
#include <cmath>  

#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "complex.h"

using namespace std;

// Min and max complex plane values
Complex  minC(-2.0, -1.2);
Complex  maxC( 1.0, 1.8);
// Global stack declaration to maintain a history
stack<double> realCStack;
stack<double> imgCStack;
int      maxIt = 2000;     // Max iterations for the set computations
double ***mbset_map = NULL;
double **points_colors = NULL;
static float updateRate = 50.0;
int nb_threads = 16; //16
int im_dim = 512; //512
double step_real = (maxC.real-minC.real)/(im_dim+1.0);
double step_imag = (maxC.imag-minC.imag)/(im_dim+1.0);
// The mutex and condition variables allow the main thread to
// know when all helper threads are completed.
pthread_mutex_t countMutex;
pthread_mutex_t exitMutex;
pthread_cond_t exitCond;
bool* localSense; // array of bools, one per thread
bool globalSense;
int startCount;
int xpos1 = 0; 
int ypos1 = 0;
int xpos2 = 0;
int ypos2 = 0;
bool selecting = false;

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
	
	//cout << "Thread number " << myId << "entered My Barrier" << endl;
	
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

void generate_random_colors(int color_nb){
  points_colors = new double*[color_nb];
  srand(time(NULL));
  for (int j=0;j<color_nb;j++){
	  points_colors[j] = new double[3];
	  points_colors[j][0]= (double) rand() / (RAND_MAX);
	  points_colors[j][1]= (double) rand() / (RAND_MAX);
      points_colors[j][2]= (double) rand() / (RAND_MAX);
  }
}

void draw_MBSet(){
	for (int i=0;i<im_dim;i++){
		for (int j=0;j<im_dim;j++){
			glColor3f(mbset_map[i][j][0], mbset_map[i][j][1], mbset_map[i][j][2]);
			glBegin(GL_POINTS);
				glVertex2f(i, j);
			glEnd();
		}
	}
}

void ComputeMBSetRow(int row_nb)
{
	double real_part = minC.real+(row_nb+1)*step_real;
	for (int col = 0; col < im_dim; col++){
		double imag_part = minC.imag+(col+1)*(step_imag);
		Complex z_0 = Complex(real_part,imag_part);
		Complex magnitude = z_0.Mag();
		double magnitude2 = z_0.Mag2();
		Complex z_n = Complex(z_0.real,z_0.imag);
		int iter_n = 0;
		while (iter_n < maxIt-1 && sqrt(magnitude2)<=2.0){
			z_n = z_0.operator +(z_n.operator *(z_n));
			magnitude2 = z_n.Mag2();
			iter_n++;
		}
		if (sqrt(magnitude2)<=2.0){
			// mag<=2 : in and black
			mbset_map[row_nb][col][0] = 0.0; mbset_map[row_nb][col][1] = 0.0; mbset_map[row_nb][col][2] =0.0; 
		} else {
			// mag>2 : out and color
			mbset_map[row_nb][col][0] = points_colors[iter_n][0];
			mbset_map[row_nb][col][1] = points_colors[iter_n][1];
			mbset_map[row_nb][col][2] = points_colors[iter_n][2];
			//mbset_map[row_nb][col][0] = 1.0; mbset_map[row_nb][col][1] = 1.0; mbset_map[row_nb][col][2] = 1.0;
		}
	}
} 

void* ComputeMBSetThread(void* v)
{ // This is the thread starting point.  "v" is the thread number
	
  unsigned long myId = (unsigned long)v;
  // Calculate MBSet for assigned rows
  int rowsPerThread = im_dim / nb_threads;
  int startingRow = myId * rowsPerThread;
  for (int r = 0; r < rowsPerThread; ++r){
	  int thisRow = startingRow + r;
	  ComputeMBSetRow(thisRow);
  }
  MyBarrier(myId);
  return 0;
}

void ComputeMBSet(){
  // Computing MBSet map values (color) with nb_threads (16) threads
  for (int i = 0; i < nb_threads; ++i){
	// Now create the thread
	pthread_t pt; // pThread variable (output param from create)
	// Third param is the thread starting function
	// Fourth param is passed to the thread starting function
	pthread_create(&pt, 0, ComputeMBSetThread, (void*)i);
  }
  // Wait for all threads complete
  MyBarrier(nb_threads);
}

void plot_square(){
	glLoadIdentity();
	glColor3f(1.0, 1.0, 1.0);
	glLineWidth(2.0);
	glBegin(GL_LINE_LOOP);
		glVertex2f(xpos1,ypos1); 
		glVertex2f(xpos2,ypos1);
		glVertex2f(xpos2,ypos2);
		glVertex2f(xpos1,ypos2);
	glEnd();
}

void display(void)
{ // Your OpenGL display code here
  static int pass;
  //cout << "Displaying pass " << ++pass << endl;
  // clear all
  glClear(GL_COLOR_BUFFER_BIT);
  // Clear the matrix
  glLoadIdentity();
  // Set the viewing transformation
  gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  
  // Draw MBSet
  draw_MBSet();
  if (selecting = true && ypos2!=0 && xpos2 !=0){
	  plot_square();
  }
  //cout << "End of display" << endl;
  glutSwapBuffers();
}

void init()
{ // Your OpenGL initialization code here
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_FLAT);
  // Initialization of mbset_map and previous_mbset_map
  mbset_map = new double**[im_dim];
  for (int k=0;k<im_dim;k++){
	  mbset_map[k] = new double*[im_dim];
	  for (int j=0;j<im_dim;j++){
		  mbset_map[k][j] = new double[3];
	  }
  }
}

void reshape(int w, int h)
{ // Your OpenGL window reshape code here
  glViewport(0,0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, (GLdouble)w, (GLdouble)0.0, h, (GLdouble)-w, (GLdouble)w);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void mouse(int button, int state, int x, int y)
{ // Your mouse click processing here
  // state == 0 means pressed, state != 0 means released
  // Note that the x and y coordinates passed in are in
  // PIXELS, with y = 0 at the top.
  // int xpos, ypos;
  // glfwGetMousePos(&xpos, &ypos);
  
  if (button==0 && state==0) {
	  xpos1 = x;
	  ypos1 = im_dim-y;
	  selecting = true;
	  cout << "Start of selection" << endl;
  } else if (state !=0 && xpos2 != 0 && ypos2 != 0){
	  selecting = false;
	  cout << "End of selection" << endl;
	  
	  // Store previous C values
	  realCStack.push(minC.real);
	  realCStack.push(maxC.real);
	  imgCStack.push(minC.imag);
	  imgCStack.push(maxC.imag);
	  
	  // Adjusting pos1 and pos2 for pos1 to be the min value (lower left) and pos2 the max value (upper right)
	  int x_temp = xpos2;
	  int y_temp = ypos2;
	  xpos2 = std::max(xpos1,x_temp);
	  xpos1 = std::min(xpos1,x_temp);
	  ypos2 = std::max(ypos1,y_temp);
	  ypos1 = std::min(ypos1,y_temp);
	  
	  // Update new C value 
	  double temp_minr = minC.real;
	  double temp_minim = minC.imag;
	  minC.real = temp_minr+(xpos1+1)*step_real; //   => DEBUG !!
	  maxC.real = temp_minr+(xpos2+1)*step_real;
	  minC.imag = temp_minim+(ypos1+1)*step_imag;
	  maxC.imag = temp_minim+(ypos2+1)*step_imag;
	  
	  step_real = (maxC.real-minC.real)/(im_dim+1.0);
	  step_imag = (maxC.imag-minC.imag)/(im_dim+1.0);
	  
	  // Display
	  ComputeMBSet();
	  glutPostRedisplay();  
	  
	  // initialize xpos2 and ypos2 back to 0 for the next selection
	  xpos2 = 0;
	  ypos2 = 0;
  }
}

void motion(int x, int y)
{ // Your mouse motion here, x and y coordinates are as above
	int y_temp = im_dim-y;
	if (abs(x-xpos1)>=abs(y_temp-ypos1) && x!=xpos1){
		xpos2 = xpos1+(x-xpos1)*abs(y_temp-ypos1)/abs(x-xpos1);
		ypos2 = y_temp;
	} else if (y_temp!=ypos1){
		xpos2 = x;
		ypos2 = ypos1+(y_temp-ypos1)*abs(x-xpos1)/abs(y_temp-ypos1);
	}
}

void keyboard(unsigned char c, int x, int y) // SEE IF NEED TO SAVE IMAGE OR RECOMPUTE IT FROM C VALUES
{ // Your keyboard processing here
	if (c == 'b'){
		if( (!imgCStack.empty()) && (!realCStack.empty())){
			maxC.imag = imgCStack.top();
			imgCStack.pop();
			maxC.real = realCStack.top();
			realCStack.pop();
			minC.imag = imgCStack.top();
			imgCStack.pop();
			minC.real = realCStack.top();
			realCStack.pop();
			step_real = (maxC.real-minC.real)/(im_dim+1.0);
			step_imag = (maxC.imag-minC.imag)/(im_dim+1.0);
			ComputeMBSet();
			glutPostRedisplay();
		} else {
			cout << "ERROR : cannot go back. Initial position." << endl;
		}
	} else {
		cout << "ERROR : invalid key. Should press b if you want to go back to the previous screen." << endl;
	}
}



int main(int argc, char** argv)
{
  // Initialize OpenGL, but only on the "master" thread or process.
  // See the assignment writeup to determine which is "master" 
  // and which is slave.
  MyBarrier_Init();
  generate_random_colors(maxIt);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(im_dim, im_dim);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("MBSet");
  init();
  ComputeMBSet();
  glutIdleFunc(display);
  glutDisplayFunc(display);  
  glutReshapeFunc(reshape);
  glutMotionFunc(motion);
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboard);
  glutMainLoop();
  
  return 0;
}

