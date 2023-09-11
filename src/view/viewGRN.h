/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#define GL_GLEXT_PROTOTYPES
#include <iostream>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>

struct sCellSphere
{
  int id;
  float position[3];
  float color[3];
  float radius;
};
typedef struct sCellSphere CellSphere;


GLfloat _colormap[64][3];
GLfloat _min_signal = 0.0, _max_signal = 100.0;
GLfloat _OnOffThreshold = 20.0; // Default: 1.0;
FILE * trajFile = NULL;
float _worldsize = 0.0;
CellSphere * _mycells = NULL;
int _popsize = 0;
float _time = 0.0;
int _line = 0;
int _pause = 0;
GLuint _SphereDList;
float _speed = 0.02;  /* default is 0.02 frames per second, i.e. wait 50 ms between two frames */
enum actions { MOVE_EYE, TWIST_EYE, ZOOM, MOVE_NONE };
GLint _action;
GLdouble _xStart = 0.0, _yStart = 0.0;
GLfloat _nearClip, _farClip, _distance, _twistAngle, _incAngle, _azimAngle;

int _width = 500;
int _height = 500;

int _frame = 0;

int _nrChannels = 4;
GLsizei stride = _nrChannels * _width;
GLsizei bufferSize = stride * _height;
std::vector<char> _buffer(bufferSize);

void init_color_map();
void read_min_max_signals();
void read_file_line();
void animate(int value);
void init();
void display();
void reshape (int w, int h);
void resetView();
void polarView( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist);
void keyboard(unsigned char key, int x, int y);
GLvoid mouse( GLint button, GLint state, GLint x, GLint y );
GLvoid motion( GLint x, GLint y );
void print_help(char* prog_name);
void paintGL();
void GenerateVideo();
