// This program is from the OpenGL Programming Guide.  It shows a robot arm
// that you can rotate by pressing the arrow keys.

#include <iostream>
//#include "Interactable.cpp"
#include <cstdlib>
#include <vector>

#include <chrono>
#include <ctime>


#ifdef __APPLE_CC__
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif

#include "atoms.cu"


static Atom* all_objects[NUMBER_OF_ATOMS];

// The robot arm is specified by (1) the angle that the upper arm makes
// relative to the x-axis, called shoulderAngle, and (2) the angle that the
// lower arm makes relative to the upper arm, called elbowAngle.  These angles
// are adjusted in 5 degree increments by a keyboard callback.


// Handles the reshape event by setting the viewport so that it takes up the
// whole visible region, then sets the projection matrix to something reason-
// able that maintains proper aspect ratio.
void reshape(GLint w, GLint h) {

  glViewport(0, 0, w, h);   //This sets up the viewport so that the coordinates (0, 0) are at the top left of the window

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(65.0, GLfloat(w)/GLfloat(h), 1.0, 20.0);
  gluLookAt( 0,0,10, 0,0,0, 0,1,0);

  //Set up the orthographic projection so that coordinates (0, 0) are in the top left
  //and the minimum and maximum depth is -10 and 10. To enable depth just put in
  //glEnable(GL_DEPTH_TEST)

  glutPostRedisplay();

}

void initalize_atoms(){

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){

    Atom* new_atom = (Atom*)malloc(sizeof(Atom));

    glm::vec2 startingPosition;
    startingPosition = glm::vec2( std::fmod((2*i*ATOM_RADIUS*ATOM_SPACING) + ATOM_RADIUS, STARTING_WIDTH / 50 )- (STARTING_WIDTH / 100),
                                   floor(((2*i*ATOM_RADIUS*ATOM_SPACING) + ATOM_RADIUS) / (STARTING_WIDTH / 50 )) * ATOM_SPACING * 2 * ATOM_RADIUS + ATOM_RADIUS - (STARTING_HEIGHT / 100));

    glm::vec2 startingVelocity;
    startingVelocity = glm::vec2(((rand() % 20000) / 1000.00) - 10,((rand() % 20000) / 1000.00) - 10);


    GLfloat color[3] = {(rand() % 100) / (float)100.0, (rand() % 100) / (float)100.0,(rand() % 100) / (float)100.0 };

    new_atom->setPosition(startingPosition);
    new_atom->setVelocity(startingVelocity);
    new_atom->setColor(color);

    new_atom->print();

    //std::cout << "Created Atom - \n";

    all_objects[i] = new_atom;
    /*
    all_objects[i]->setPosition(startingPosition);
    all_objects[i]->setVelocity(startingVelocity);
    all_objects[i]->setColor(color);
    */
  }

}

std::chrono::system_clock::time_point start;
    // Some computation here
std::chrono::system_clock::time_point end;

void timer(int time){

    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "FPS = " <<  1/ elapsed_seconds.count() << "\n";

    glutPostRedisplay();

    glutTimerFunc( (1.0 / 60) * 1000, timer, time);

    start = std::chrono::system_clock::now();
}

void DrawRectWithOpenGL(float x, float y, float w, float h)
{

  
    glPushMatrix();  //Make sure our transformations don't affect any other transformations in other code
    glTranslatef(x, y, 0.0f);  //Translate rectangle to its assigned x and y position
    //Put other transformations here
    glBegin(GL_QUADS);   //We want to draw a quad, i.e. shape with four sides
      glColor3f(1, 0, 0); //Set the colour to red 
      glVertex2f(0, 0);            //Draw the four corners of the rectangle
      glVertex2f(0, h);
      glVertex2f(w, h);
      glVertex2f(w, 0);       
    glEnd();
  glPopMatrix();
}

void special(int key, int, int) {
  //return;
  switch (key) {
    //case GLUT_KEY_RIGHT:  glutPostRedisplay(); break;
    default: return;
  }
}

// Clears the current window and draws a triangle.
void display() {

  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){
    //all_objects[i]->print();
    all_objects[i]->move(1.0 / 60);
  }
  //std::cout << "+++++++++++++++++++++++++++++++++++++++=\n";

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){
    //std::cout << "Collision for ball - " << i << "---------\n";
    all_objects[i]->update_vel = false;
    all_objects[i]->handleCollisions(all_objects);
  }   

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){
    all_objects[i]->updateVelocity();
    all_objects[i]->draw();
  }   

  // Flush drawing command buffer to make drawing happen as soon as possible.
  glFlush();
}

// Initializes GLUT, the display mode, and main window; registers callbacks;
// enters the main event loop.
int main(int argc, char** argv) {

  start = std::chrono::system_clock::now();

  initalize_atoms();

  //std::cout << " x - " << all_objects[0]->coords[0] << ", y - " << all_objects[0]->coords[1] << "\n";

  // Use a single buffered window in RGB mode (as opposed to a double-buffered
  // window or color-index mode).
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

  // Position window at (80,80)-(480,380) and give it a title.
  glutInitWindowPosition(80, 80);
  glutInitWindowSize(STARTING_WIDTH, STARTING_HEIGHT);
  glutCreateWindow("A Simple Triangle");

  // Tell GLUT that whenever the main window needs to be repainted that it
  // should call the function display().
  glutDisplayFunc(display);

  glutReshapeFunc(reshape);

  glutTimerFunc(100, timer, 1);

  glutSpecialFunc(special);

  // Tell GLUT to start reading and processing events.  This function
  // never returns; the program only exits when the user closes the main
  // window or kills the process.
  glutMainLoop();
}