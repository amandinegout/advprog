// Draw an Icosahedron
// ECE4893/8893 Project 4
// YOUR NAME HERE

#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

static const float DEG2RAD = M_PI/180;
static float updateRate = 10.0;
GLfloat rotlim = (GLfloat) 360.0;

#define NFACE 20
#define NVERTEX 12

#define X .525731112119133606 
#define Z .850650808352039932

// These are the 12 vertices for the icosahedron
static GLfloat vdata[NVERTEX][3] = {    
   {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
   {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
   {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
};

// These are the 20 faces.  Each of the three entries for each 
// vertex gives the 3 vertices that make the face.
static GLint tindices[NFACE][3] = { 
   {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
   {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
   {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
   {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

int testNumber; // Global variable indicating which test number is desired
int depth; // Global variable indicating which depth is chosen for test 5 and 6
double **triangles_colors = NULL;

void generate_random_colors(int color_nb){
  triangles_colors = new double*[color_nb];
  srand(time(NULL));
  //cout << "number of colours :" << color_nb << endl;
  for (int j=0;j<color_nb;j++){
	  triangles_colors[j] = new double[3];
	  triangles_colors[j][0]= (double) rand() / (RAND_MAX);
	  triangles_colors[j][1]= (double) rand() / (RAND_MAX);
      triangles_colors[j][2]= (double) rand() / (RAND_MAX);
  }
}

void DrawTriangles(GLfloat* v1, GLfloat* v2, GLfloat* v3)
{
  // Drawing triangles
  glBegin(GL_TRIANGLES);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
  glEnd(); 
}

void DrawLines(GLfloat* v1, GLfloat* v2, GLfloat* v3)
{
  // Drawing white lines
  glBegin(GL_LINE_LOOP);
	//glLineWidth(0.5);
    glColor3f(1.0, 1.0, 1.0);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
  glEnd();
}

void normalize(GLfloat* vect){
	GLfloat norm = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
	vect[0] = (GLfloat) vect[0]/norm;
	vect[1] = (GLfloat) vect[1]/norm;
	vect[2] = (GLfloat) vect[2]/norm;
}

void subdivide(GLfloat* v1, GLfloat* v2, GLfloat* v3, int depth_sub, int color_count){
	//cout << "color count: " << color_count << endl;
	if (depth_sub==0){
		glColor3f(0.0, 0.0, 1.0);
		glColor3f(triangles_colors[color_count][0],triangles_colors[color_count][1],triangles_colors[color_count][2]);
		DrawTriangles(v1,v2,v3);
		DrawLines(v1,v2,v3);
	}
	else{
		GLfloat v12[3];
		GLfloat v23[3];
		GLfloat v31[3];
		for(int i = 0; i < 3; i++)
		{
			v12[i] = (GLfloat) (v1[i] + v2[i])/2.0; 
			v23[i] = (GLfloat) (v2[i] + v3[i])/2.0; 
			v31[i] = (GLfloat) (v3[i] + v1[i])/2.0; 
		}
		normalize(v12);
		normalize(v23);
		normalize(v31);
		subdivide(v1,v12,v31,depth_sub-1,color_count*4);
		subdivide(v2,v12,v23,depth_sub-1,color_count*4+1);
		subdivide(v3,v23,v31,depth_sub-1,color_count*4+2);
		subdivide(v12,v23,v31,depth_sub-1,color_count*4+3);
	}
}

// Test cases.  Fill in your code for each test case
void Test1()
{
  for (int i=0;i<NFACE;i++){
      glColor3f(triangles_colors[i][0],triangles_colors[i][1],triangles_colors[i][2]);
      DrawTriangles(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0], &vdata[tindices[i][2]][0]);
	  DrawLines(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0], &vdata[tindices[i][2]][0]);
  }
}

void Test2()
{
  // Rotation
  static GLfloat rotation = 0.0;
  glRotatef(rotation, 1.0, 1.0, 0.0);
  rotation += 1.0;
  if(rotation > rotlim)
    rotation = (GLfloat)0.0;
  //static GLfloat rotX = 0.0;
  //static GLfloat rotY = 0.0;
  //glRotatef(rotX, 1.0, 0.0, 0.0);
  //glRotatef(rotY, 0.0, 1.0, 0.0);
  //rotX += 1.0;
  //rotY += 1.0;
  //if(rotX > rotlim)
    //rotX = (GLfloat)0.0;
  //if(rotY > rotlim)
    //rotY = (GLfloat)0.0;
  Test1();
}

void Test3()
{
	for (int i=0;i<NFACE;i++){
      subdivide(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0], &vdata[tindices[i][2]][0],depth,i);
    }
}

void Test4()
{
  // Rotation
  static GLfloat rotX = 0.0;
  static GLfloat rotY = 0.0;
  glRotatef(rotX, 1.0, 0.0, 0.0);
  glRotatef(rotY, 0.0, 1.0, 0.0);
  rotX += 1.0;
  rotY += 1.0;
  if(rotX > rotlim)
    rotX = (GLfloat)0.0;
  if(rotY > rotlim)
    rotY = (GLfloat)0.0;
	Test3();
}

void Test5(int depth)
{
	for (int i=0;i<NFACE;i++){
      subdivide(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0], &vdata[tindices[i][2]][0],depth,i);
    }
}

void Test6(int depth)
{
  // Rotation
  static GLfloat rotX = 0.0;
  static GLfloat rotY = 0.0;
  glRotatef(rotX, 1.0, 0.0, 0.0);
  glRotatef(rotY, 0.0, 1.0, 0.0);
  rotX += 1.0;
  rotY += 1.0;
  if(rotX > rotlim)
    rotX = (GLfloat)0.0;
  if(rotY > rotlim)
    rotY = (GLfloat)0.0;
  Test5(depth);
}

void display(void){
  static int pass;
  cout << "Displaying pass " << ++pass << endl;
  glEnable(GL_LINE_SMOOTH);
  glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
  // clear all
  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);
  // Depth buffer
  glEnable(GL_DEPTH_TEST);
  // Clear the matrix
  glLoadIdentity();
  // Set the viewing transformation
  gluLookAt(0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  glPushMatrix();
  glTranslatef(250, 250, 0);
  glScalef(250.0, 250.0, 250.0);
  
  // Select which test to run
  switch(testNumber){
	case 1:
	  Test1();
	  break;
	case 2:
	  Test2();
	  break;
	case 3:
	  Test3();
	  break;
	case 4:
	  Test4();
	  break;
	case 5:
	  Test5(depth);
	  break;
	case 6:
	  Test6(depth);
	  break;
	default:
	  cout << "Wrong test number : between 1 and 6" << endl;
	  abort();
  }
  glPopMatrix();
  glutSwapBuffers(); // If double buffering
  //glFlush();
}

void init()
{
  //select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_FLAT);
  glLoadIdentity();
}

void reshape(int w, int h)
{
  glViewport(0,0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, (GLdouble)w, (GLdouble)0.0, h, (GLdouble)w, (GLdouble)-w);
  //glFrustum(0.0, (GLdouble)w, (GLdouble)0.0, h, (GLdouble)1, (GLdouble)40);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void timer(int)
{
  glutPostRedisplay();
  glutTimerFunc(1000.0 / updateRate, timer, 0);
}

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      std::cout << "Usage: icosahedron testnumber" << endl;
      exit(1);
    }
  // Set the global test number
  testNumber = atol(argv[1]);
  if (testNumber==3 || testNumber==4){
	  depth=1;
  }
  // Set the global depth
  if (argc==3){
	  depth = atol(argv[2]);
	  if (depth>5){
		  cout << "Error : Depth should not be more than 5" << endl;
		  abort();
	  }
  }
  generate_random_colors(NFACE*pow(4,depth));
  // Initialize glut  and create your window here
  // Set your glut callbacks here
  // Enter the glut main loop here
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  //glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(500, 500);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Icosahedron");
  init();
  glutDisplayFunc(display);  
  glutReshapeFunc(reshape);
  glutTimerFunc(1000.0 / updateRate, timer, 0);
  glutMainLoop();
  return 0;
}

