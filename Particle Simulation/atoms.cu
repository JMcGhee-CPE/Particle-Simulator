#pragma once

#include "glm/glm.hpp"
#include "glm/gtx/projection.hpp"

#include <windows.h>

#include <stdio.h>
#include <cstring>

#define ATOM_RADIUS 0.05
#define ATOM_RESOLUTION 12
#define NUMBER_OF_ATOMS 2048
#define ATOM_SPACING 2
#define STARTING_WIDTH 800
#define STARTING_HEIGHT 600
#define COLLISION_MOVE 0.1

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

void DrawCircle(float cx, float cy, float r, int num_segments, const GLfloat* color)
{

    glPushMatrix();

      glColor3fv(color);

      glBegin(GL_POLYGON);
      for(int i = 0; i < num_segments; i++)
      {
          float theta = 2.0f * 3.1415926f * float(i) / float(num_segments);//get the current angle

          float x = r * cosf(theta);//calculate the x component
          float y = r * sinf(theta);//calculate the y component

          glVertex2f(x + cx, y + cy);//output vertex

      }
      glEnd();

    glPopMatrix();
}


class Atom {


public:

  glm::vec2 coords;
  glm::vec2 vel;
  glm::vec2 new_pos;
  glm::vec2 new_vel;
  bool update_vel = false;

  GLfloat color[3];


  CUDA_CALLABLE_MEMBER
  void setPosition(glm::vec2 position){
    coords = position;
  }

  CUDA_CALLABLE_MEMBER
  void setColor(GLfloat* new_color){
    color[0] = new_color[0];
    color[1] = new_color[1];
    color[2] = new_color[2];
  }

  CUDA_CALLABLE_MEMBER
  void setVelocity(glm::vec2 velocity){
    vel = velocity;
  }

  CUDA_CALLABLE_MEMBER
  glm::vec2 getVelocity(){
    return vel;
  }

  CUDA_CALLABLE_MEMBER
  glm::vec2 getPosition(){
    return coords;
  }

  CUDA_CALLABLE_MEMBER
  void move(float time){
    if (coords[1] >= (STARTING_HEIGHT / 100) || coords[1] <= -(STARTING_HEIGHT / 100)){
      vel *= glm::vec2(1,-1);
      if(coords[1] >= (STARTING_HEIGHT / 100)) {
        coords[1] = (STARTING_HEIGHT / 100);
      } else {
        coords[1] = -(STARTING_HEIGHT / 100);
      }
    }
    if (coords[0] >= (STARTING_WIDTH / 100) || coords[0] <= -(STARTING_WIDTH / 100)){
      vel *= glm::vec2(-1,1);
      if(coords[0] >= (STARTING_WIDTH / 100)) {
        coords[0] = (STARTING_WIDTH / 100);
      } else {
        coords[0] = -(STARTING_WIDTH / 100);
      }
    }
    glm::vec2 displacement = vel;
    displacement *= time;

    coords += displacement;
  }

  CUDA_CALLABLE_MEMBER
  void handleCollisions( Atom** potential_Collisions ){
    glm::vec2 new_vector(0,0);
    new_pos = coords;
    int num_collisions = 0;
    for(int i = 0; i < NUMBER_OF_ATOMS; i++){
      if((potential_Collisions[i] != this) && (abs(glm::distance( potential_Collisions[i]->coords, coords )) <= (2 * ATOM_RADIUS))){

        // Formula for elastic collisions https://williamecraver.wixsite.com/elastic-equations

        glm::vec2 position_vector = coords;
        position_vector -= potential_Collisions[i]->coords;

        glm::vec2 unit_vector = position_vector;
        unit_vector /= glm::length(position_vector);


        glm::vec2 offset_vector = unit_vector;
        offset_vector *= (2 * ATOM_RADIUS * COLLISION_MOVE);

        glm::vec2 velocity_difference = vel;
        velocity_difference -= potential_Collisions[i]->vel;

        float dotproduct = glm::dot(velocity_difference, position_vector);

        glm::vec2 contribution = position_vector;
        contribution *= -dotproduct / powf( glm::length(position_vector), 2);

        contribution += vel;

        new_vector += contribution;

        new_pos += offset_vector;

        num_collisions++;
      }
    }

    if(num_collisions > 0){
      //std::cout << "Num Collisions = " << num_collisions << " \n";
      new_vector /= num_collisions;
      //std::cout << "New Vector = {" << new_vector[0] << ", " << new_vector[1] << "} \n";
      new_vel = new_vector;
      //std::cout << "New vector - { " << new_vel[0] << ", " << new_vel[1] << " }\n";
      update_vel = true;
    }

    // Else no collisions, velocity unchanged

  }

  CUDA_CALLABLE_MEMBER
  void updateVelocity(){
    if(update_vel){
      //std::cout << "velocity updated!\n";
      vel = new_vel;
      coords = new_pos;
      update_vel = false;
    }
  }

  void draw(){

    // Set every pixel in the frame buffer to the current clear color.

    DrawCircle(coords[0], coords[1], ATOM_RADIUS, ATOM_RESOLUTION, color);
  }

  void print(){
    std::cout << "Coords: { " << coords[0] << ", " << coords[1] << " }; Vel: { " << vel[0] << ", " << vel[1] << "}\n";
  }

};