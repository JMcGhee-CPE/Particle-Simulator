#include <iostream>
#include <cstdlib>
#include <vector>

#include <chrono>
#include <ctime>
#include <stdio.h>

#ifdef __APPLE_CC__
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif
#include "atoms.cu"

#include <curand.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 512

__global__
void initalize_atoms(Atom** d_atoms)
{

  int linearIdx = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
  
  curandState_t state;

  /* we have to initialize the state */
  curand_init(0, /* the seed controls the sequence of random values that are produced */
              0, /* the sequence number is only important with multiple cores */
              linearIdx * 5, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
              &state);


  if(linearIdx < NUMBER_OF_ATOMS){
    Atom* new_atom = (Atom*)malloc(sizeof(Atom));

    int x_position_i = ((2*linearIdx*ATOM_RADIUS*ATOM_SPACING) + ATOM_RADIUS) * 10000;

    int y_position_i = floor(((double)x_position_i) / (2*STARTING_WIDTH * 100)) * 2 * ATOM_RADIUS*ATOM_SPACING * 10000;

    x_position_i %= 200*STARTING_WIDTH;

    double x_position = (x_position_i - (STARTING_WIDTH * 100)) / 10000.0;
    double y_position = (y_position_i - (STARTING_HEIGHT * 100)) / 10000.0;

    glm::vec2 startingPosition;
    startingPosition = glm::vec2( x_position, y_position);

    glm::vec2 startingVelocity;
    startingVelocity = glm::vec2(((curand(&state) % 20000) / 1000.00) - 10,((curand(&state) % 20000) / 1000.00) - 10);

    new_atom->coords = startingPosition;
    new_atom->vel = startingVelocity;

    new_atom->color[0] = (curand(&state) % 100) / (float)100.0;
    new_atom->color[1] = (curand(&state) % 100) / (float)100.0;
    new_atom->color[2] = (curand(&state) % 100) / (float)100.0;

    d_atoms[linearIdx] = new_atom;

    //printf("Idx = %d; Coords: { %f, %f }; Vel: { %f, %f }\n", linearIdx, d_atoms[linearIdx]->coords[0], d_atoms[linearIdx]->coords[1], d_atoms[linearIdx]->vel[0], d_atoms[linearIdx]->vel[1]);
   
  }
}

__global__
void move_atoms(Atom** d_atoms)
{
  int linearIdx = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
  if(linearIdx < NUMBER_OF_ATOMS){
   float time = 1.0 / 60;
    if (d_atoms[linearIdx]->coords[1] >= (STARTING_HEIGHT / 100) || d_atoms[linearIdx]->coords[1] <= -(STARTING_HEIGHT / 100)){
      d_atoms[linearIdx]->vel *= glm::vec2(1,-1);
      if(d_atoms[linearIdx]->coords[1] >= (STARTING_HEIGHT / 100)) {
        d_atoms[linearIdx]->coords[1] = (STARTING_HEIGHT / 100);
      } else {
        d_atoms[linearIdx]->coords[1] = -(STARTING_HEIGHT / 100);
      }
    }
    if (d_atoms[linearIdx]->coords[0] >= (STARTING_WIDTH / 100) || d_atoms[linearIdx]->coords[0] <= -(STARTING_WIDTH / 100)){
      d_atoms[linearIdx]->vel *= glm::vec2(-1,1);
      if(d_atoms[linearIdx]->coords[0] >= (STARTING_WIDTH / 100)) {
        d_atoms[linearIdx]->coords[0] = (STARTING_WIDTH / 100);
      } else {
        d_atoms[linearIdx]->coords[0] = -(STARTING_WIDTH / 100);
      }
    }
    glm::vec2 displacement = d_atoms[linearIdx]->vel;
    displacement *= time;

    d_atoms[linearIdx]->coords += displacement;
  }
}

__global__
void collide_atoms(Atom** d_atoms)
{
  int linearIdx = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
  if(linearIdx < NUMBER_OF_ATOMS){

    glm::vec2 new_vector(0,0);
    d_atoms[linearIdx]->new_pos = d_atoms[linearIdx]->coords;
    int num_collisions = 0;
    for(int i = 0; i < NUMBER_OF_ATOMS; i++){
      bool is_current_atom = d_atoms[i] != d_atoms[linearIdx];
      //printf("Current position atom 1 = { %f, %f}; atom 2 = { %f, %f}\n", d_atoms[i]->coords[0], d_atoms[i]->coords[1], d_atoms[linearIdx]->coords[0], d_atoms[linearIdx]->coords[1]);
      float dist = sqrtf( powf(d_atoms[i]->coords[0] - d_atoms[linearIdx]->coords[0], 2) + powf(d_atoms[i]->coords[1] - d_atoms[linearIdx]->coords[1], 2));
      if(is_current_atom && (dist <= (2 * ATOM_RADIUS))){

        // Formula for elastic collisions https://williamecraver.wixsite.com/elastic-equations

        glm::vec2 position_vector = d_atoms[linearIdx]->coords;
        position_vector -= d_atoms[i]->coords;

        glm::vec2 unit_vector = position_vector;
        unit_vector /= glm::length(position_vector);


        glm::vec2 offset_vector = unit_vector;
        offset_vector *= (2 * ATOM_RADIUS * COLLISION_MOVE);

        glm::vec2 velocity_difference = d_atoms[linearIdx]->vel;
        velocity_difference -= d_atoms[i]->vel;

        float dotproduct = glm::dot(velocity_difference, position_vector);

        glm::vec2 contribution = position_vector;
        contribution *= -dotproduct / powf( glm::length(position_vector), 2);

        contribution += d_atoms[linearIdx]->vel;

        new_vector += contribution;

        d_atoms[linearIdx]->new_pos += offset_vector;

        num_collisions++;
      }
    }

    if(num_collisions > 0){
      //std::cout << "Num Collisions = " << num_collisions << " \n";
      new_vector /= num_collisions;
      //std::cout << "New Vector = {" << new_vector[0] << ", " << new_vector[1] << "} \n";
      d_atoms[linearIdx]->new_vel = new_vector;
      //std::cout << "New vector - { " << new_vel[0] << ", " << new_vel[1] << " }\n";
      d_atoms[linearIdx]->update_vel = true;
    }

    // Else no collisions, velocity unchanged

  }
}

__global__
void update_atoms(Atom** d_atoms)
{
  int linearIdx = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
  if(linearIdx < NUMBER_OF_ATOMS){
    if(d_atoms[linearIdx]->update_vel){
      //std::cout << "velocity updated!\n";
      d_atoms[linearIdx]->vel = d_atoms[linearIdx]->new_vel;
      d_atoms[linearIdx]->coords = d_atoms[linearIdx]->new_pos;
      d_atoms[linearIdx]->update_vel = false;
    }
  }
}

__global__
void prepare_draw(Atom** d_atoms, Atom* d_atoms_ref)
{
  int linearIdx = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
  if(linearIdx < NUMBER_OF_ATOMS){
    memcpy(&d_atoms_ref[linearIdx], d_atoms[linearIdx],sizeof(Atom));
  }
}

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


void special(int key, int, int) {
  //return;
  switch (key) {
    //case GLUT_KEY_RIGHT:  glutPostRedisplay(); break;
    default: return;
  }
}

Atom **atoms, **d_atoms, *atom_refs, *d_atom_ref;

// Clears the current window and draws a triangle.
void display() {

  int number_of_blocks = (int)ceil(NUMBER_OF_ATOMS / (float)THREADS_PER_BLOCK);

  move_atoms<<< number_of_blocks, THREADS_PER_BLOCK>>>(d_atoms);

  cudaDeviceSynchronize();

  collide_atoms<<< number_of_blocks, THREADS_PER_BLOCK>>>(d_atoms);

  cudaDeviceSynchronize();

  update_atoms<<< number_of_blocks, THREADS_PER_BLOCK>>>(d_atoms);

  cudaDeviceSynchronize();

  prepare_draw<<< number_of_blocks, THREADS_PER_BLOCK>>>(d_atoms, d_atom_ref);

  cudaMemcpy(atom_refs, d_atom_ref, NUMBER_OF_ATOMS * sizeof(Atom), cudaMemcpyDeviceToHost);

  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){
    //atom_refs[i].print();
    atom_refs[i].draw();
  }   

  // Flush drawing command buffer to make drawing happen as soon as possible.
  glFlush();
}

int main(int argc, char** argv)
{

  atoms = (Atom**)malloc(NUMBER_OF_ATOMS* sizeof(Atom*));
  atom_refs = (Atom*)malloc(NUMBER_OF_ATOMS * sizeof(Atom));

  for(int i = 0; i < NUMBER_OF_ATOMS; i++){
    atoms[i] = (Atom*)malloc(sizeof(Atom));
  }

  cudaMalloc((void**)&d_atoms, NUMBER_OF_ATOMS * sizeof(Atom*));
  cudaMalloc(&d_atom_ref, NUMBER_OF_ATOMS * sizeof(Atom));
  //std::cout << "Blocks  = { " << blocks.x << ", " << blocks.y << ", " << blocks.z << "}\n";
  //std::cout << "thread_num  = { " << thread_num.x << ", " << thread_num.y << ", " << thread_num.z << "}\n"; 

  int number_of_blocks = (int)ceil(NUMBER_OF_ATOMS / (float)THREADS_PER_BLOCK);

  //std::cout << "Number of blocks = " << number_of_blocks << "\n";

  initalize_atoms<<< number_of_blocks, THREADS_PER_BLOCK>>>(d_atoms);

  cudaDeviceSynchronize();

  start = std::chrono::system_clock::now();

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


  cudaFree(d_atoms);
  free(atoms);
}