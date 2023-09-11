
#define GL_GLEXT_PROTOTYPES
#define GL_SILENCE_DEPRECATION

#include <GLUT/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define DEFAULT_SIGNAL -1
#define SIGNAL_NOT_FOUND -1
#define MAX_SIG_STRLEN 16
#define MAX_NBR_SIGNAL 1024

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
FILE * trajFile = NULL;
float _worldsize = 0.0;
CellSphere * _mycells = NULL;
int _popsize = 0;
float _time = 0.0;
int _line = 0;
int _pause = 0;
GLuint _SphereDList;
float _speed = 0.02;  /* default is 0.02 frames per second, i.e. wait 50 ms between two frames */
int _first_signal_col;
int _signal_index; 
char* _signal_name;
struct list_signals { unsigned int size; char name[MAX_NBR_SIGNAL][128]; };
struct list_signals _list_signals = { .size = 0 };
char _signal_found = 0;
enum actions { MOVE_EYE, TWIST_EYE, ZOOM, MOVE_NONE };
GLint _action;
GLdouble _xStart = 0.0, _yStart = 0.0;
GLfloat _nearClip, _farClip, _distance, _twistAngle, _incAngle, _azimAngle;

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
void next_signal();
void previous_signal();
GLvoid mouse( GLint button, GLint state, GLint x, GLint y );
GLvoid motion( GLint x, GLint y );
void print_help(char* prog_name);






int main(int argc, char** argv)
{
  /* Define allowed command-line options */
  const char * options_list = "hf:s:i:p";
  static struct option long_options_list[] = {
    { "help",     no_argument, NULL, 'h' },
    { "file",     required_argument, NULL, 'f' },
    { "speed",    required_argument, NULL, 's' },
    { "signal",   required_argument, NULL, 'i' },
    { "pause",    no_argument, NULL, 'p' },
    { 0, 0, 0, 0 }
  };

  /* Get actual values of the command-line options */
  int option;
  char * filename;
  char filename_specified = 0;
  char signal_specified = 0;
  while ((option = getopt_long(argc, argv,
                               options_list, long_options_list, NULL)) != -1)
  {
    switch ( option )
    {
      case 'h' :
      {
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
      }
      case 'f' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          fprintf(stderr, "ERROR : Option -f or --file : missing argument.\n" );
          exit( EXIT_FAILURE );
        }

        filename = (char *) malloc((strlen(optarg) + 1)*sizeof(char));
        sprintf( filename, "%s", optarg );
        filename_specified = 1;
        break;
      }
      case 's' :
      {
        if ( strcmp( optarg, "" ) == 0 )
          {
            fprintf(stderr,
                "ERROR : Option -s or --speed : missing argument.\n" );
            exit( EXIT_FAILURE );
          }
        _speed = atof(optarg); /* in frames per millisecond */
        break;
      }
      case 'i' :
        if ( strcmp( optarg, "" ) == 0 )
        {
            fprintf(stderr,
                "ERROR : Option -i or --signal : missing argument.\n" );
            exit( EXIT_FAILURE );
        }
        _signal_name = (char *)malloc((strlen(optarg) + 1)*sizeof(char));
        sprintf( _signal_name, "%s", optarg );
        _signal_index = SIGNAL_NOT_FOUND; /* the index will need to be found in the header file */
        signal_specified = 1;
        break;
      case 'p':
        _pause = 1;
        break;
    }
  }

  if (!filename_specified) trajFile = fopen("trajectory.txt", "r");
  else
  {
  	trajFile = fopen(filename, "r");
  	free(filename);
  }
  if (trajFile == NULL)
  {
    fprintf(stderr, "Error, file trajectory.txt is missing.\n");
    exit(EXIT_FAILURE);
  }

  if (!signal_specified)
  {
    _signal_name = (char *)malloc(sizeof(char));
    _signal_name[0] = 0;
    _signal_index = DEFAULT_SIGNAL; /* by default */
  }

  while (_worldsize == 0.0) read_file_line();
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize (500, 500);
  glutInitWindowPosition (100, 100);
  glutCreateWindow ( "Simuscale - View cells" );
  init ();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);


  glutTimerFunc(0,animate,0);
  glutMainLoop();

  free(_signal_name);
  return 0;
}



void init_color_map()
{
  _colormap[0][0] = 0.000000; _colormap[0][1] = 0.000000; _colormap[0][2] = 0.562500;
  _colormap[1][0] = 0.000000; _colormap[1][1] = 0.000000; _colormap[1][2] = 0.625000;
  _colormap[2][0] = 0.000000; _colormap[2][1] = 0.000000; _colormap[2][2] = 0.687500;
  _colormap[3][0] = 0.000000; _colormap[3][1] = 0.000000; _colormap[3][2] = 0.750000;
  _colormap[4][0] = 0.000000; _colormap[4][1] = 0.000000; _colormap[4][2] = 0.812500;
  _colormap[5][0] = 0.000000; _colormap[5][1] = 0.000000; _colormap[5][2] = 0.875000;
  _colormap[6][0] = 0.000000; _colormap[6][1] = 0.000000; _colormap[6][2] = 0.937500;
  _colormap[7][0] = 0.000000; _colormap[7][1] = 0.000000; _colormap[7][2] = 1.000000;
  _colormap[8][0] = 0.000000; _colormap[8][1] = 0.062500; _colormap[8][2] = 1.000000;
  _colormap[9][0] = 0.000000; _colormap[9][1] = 0.125000; _colormap[9][2] = 1.000000;
  _colormap[10][0] = 0.000000; _colormap[10][1] = 0.187500; _colormap[10][2] = 1.000000;
  _colormap[11][0] = 0.000000; _colormap[11][1] = 0.250000; _colormap[11][2] = 1.000000;
  _colormap[12][0] = 0.000000; _colormap[12][1] = 0.312500; _colormap[12][2] = 1.000000;
  _colormap[13][0] = 0.000000; _colormap[13][1] = 0.375000; _colormap[13][2] = 1.000000;
  _colormap[14][0] = 0.000000; _colormap[14][1] = 0.437500; _colormap[14][2] = 1.000000;
  _colormap[15][0] = 0.000000; _colormap[15][1] = 0.500000; _colormap[15][2] = 1.000000;
  _colormap[16][0] = 0.000000; _colormap[16][1] = 0.562500; _colormap[16][2] = 1.000000;
  _colormap[17][0] = 0.000000; _colormap[17][1] = 0.625000; _colormap[17][2] = 1.000000;
  _colormap[18][0] = 0.000000; _colormap[18][1] = 0.687500; _colormap[18][2] = 1.000000;
  _colormap[19][0] = 0.000000; _colormap[19][1] = 0.750000; _colormap[19][2] = 1.000000;
  _colormap[20][0] = 0.000000; _colormap[20][1] = 0.812500; _colormap[20][2] = 1.000000;
  _colormap[21][0] = 0.000000; _colormap[21][1] = 0.875000; _colormap[21][2] = 1.000000;
  _colormap[22][0] = 0.000000; _colormap[22][1] = 0.937500; _colormap[22][2] = 1.000000;
  _colormap[23][0] = 0.000000; _colormap[23][1] = 1.000000; _colormap[23][2] = 1.000000;
  _colormap[24][0] = 0.062500; _colormap[24][1] = 1.000000; _colormap[24][2] = 0.937500;
  _colormap[25][0] = 0.125000; _colormap[25][1] = 1.000000; _colormap[25][2] = 0.875000;
  _colormap[26][0] = 0.187500; _colormap[26][1] = 1.000000; _colormap[26][2] = 0.812500;
  _colormap[27][0] = 0.250000; _colormap[27][1] = 1.000000; _colormap[27][2] = 0.750000;
  _colormap[28][0] = 0.312500; _colormap[28][1] = 1.000000; _colormap[28][2] = 0.687500;
  _colormap[29][0] = 0.375000; _colormap[29][1] = 1.000000; _colormap[29][2] = 0.625000;
  _colormap[30][0] = 0.437500; _colormap[30][1] = 1.000000; _colormap[30][2] = 0.562500;
  _colormap[31][0] = 0.500000; _colormap[31][1] = 1.000000; _colormap[31][2] = 0.500000;
  _colormap[32][0] = 0.562500; _colormap[32][1] = 1.000000; _colormap[32][2] = 0.437500;
  _colormap[33][0] = 0.625000; _colormap[33][1] = 1.000000; _colormap[33][2] = 0.375000;
  _colormap[34][0] = 0.687500; _colormap[34][1] = 1.000000; _colormap[34][2] = 0.312500;
  _colormap[35][0] = 0.750000; _colormap[35][1] = 1.000000; _colormap[35][2] = 0.250000;
  _colormap[36][0] = 0.812500; _colormap[36][1] = 1.000000; _colormap[36][2] = 0.187500;
  _colormap[37][0] = 0.875000; _colormap[37][1] = 1.000000; _colormap[37][2] = 0.125000;
  _colormap[38][0] = 0.937500; _colormap[38][1] = 1.000000; _colormap[38][2] = 0.062500;
  _colormap[39][0] = 1.000000; _colormap[39][1] = 1.000000; _colormap[39][2] = 0.000000;
  _colormap[40][0] = 1.000000; _colormap[40][1] = 0.937500; _colormap[40][2] = 0.000000;
  _colormap[41][0] = 1.000000; _colormap[41][1] = 0.875000; _colormap[41][2] = 0.000000;
  _colormap[42][0] = 1.000000; _colormap[42][1] = 0.812500; _colormap[42][2] = 0.000000;
  _colormap[43][0] = 1.000000; _colormap[43][1] = 0.750000; _colormap[43][2] = 0.000000;
  _colormap[44][0] = 1.000000; _colormap[44][1] = 0.687500; _colormap[44][2] = 0.000000;
  _colormap[45][0] = 1.000000; _colormap[45][1] = 0.625000; _colormap[45][2] = 0.000000;
  _colormap[46][0] = 1.000000; _colormap[46][1] = 0.562500; _colormap[46][2] = 0.000000;
  _colormap[47][0] = 1.000000; _colormap[47][1] = 0.500000; _colormap[47][2] = 0.000000;
  _colormap[48][0] = 1.000000; _colormap[48][1] = 0.437500; _colormap[48][2] = 0.000000;
  _colormap[49][0] = 1.000000; _colormap[49][1] = 0.375000; _colormap[49][2] = 0.000000;
  _colormap[50][0] = 1.000000; _colormap[50][1] = 0.312500; _colormap[50][2] = 0.000000;
  _colormap[51][0] = 1.000000; _colormap[51][1] = 0.250000; _colormap[51][2] = 0.000000;
  _colormap[52][0] = 1.000000; _colormap[52][1] = 0.187500; _colormap[52][2] = 0.000000;
  _colormap[53][0] = 1.000000; _colormap[53][1] = 0.125000; _colormap[53][2] = 0.000000;
  _colormap[54][0] = 1.000000; _colormap[54][1] = 0.062500; _colormap[54][2] = 0.000000;
  _colormap[55][0] = 1.000000; _colormap[55][1] = 0.000000; _colormap[55][2] = 0.000000;
  _colormap[56][0] = 0.937500; _colormap[56][1] = 0.000000; _colormap[56][2] = 0.000000;
  _colormap[57][0] = 0.875000; _colormap[57][1] = 0.000000; _colormap[57][2] = 0.000000;
  _colormap[58][0] = 0.812500; _colormap[58][1] = 0.000000; _colormap[58][2] = 0.000000;
  _colormap[59][0] = 0.750000; _colormap[59][1] = 0.000000; _colormap[59][2] = 0.000000;
  _colormap[60][0] = 0.687500; _colormap[60][1] = 0.000000; _colormap[60][2] = 0.000000;
  _colormap[61][0] = 0.625000; _colormap[61][1] = 0.000000; _colormap[61][2] = 0.000000;
  _colormap[62][0] = 0.562500; _colormap[62][1] = 0.000000; _colormap[62][2] = 0.000000;
  _colormap[63][0] = 0.500000; _colormap[63][1] = 0.000000; _colormap[63][2] = 0.000000;
}


void read_min_max_signals()
{
  FILE * minmaxfile = fopen("normalization.txt", "r");
  if (minmaxfile == NULL)
  {
    printf("No file called normalization.txt in current directory.\n");
    printf("Using default values for coloring: min=0.0 and max=100.0.\n");
    return;
  }

  int retval;
  char c;
  char str[128];
  retval = fscanf(minmaxfile, "%f", &_max_signal);
  while (retval != 1)
  {
    if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}

    if (fscanf(minmaxfile, "%s", str) == 1)
    {
      if (str[0]=='#' || (strncmp(str, _signal_name, 127) != 0 && strlen(_signal_name)) ) 
      {
        /* ignore the whole line */
        /* printf("comment line :\n"); */
        do {retval = fscanf(minmaxfile, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
        if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}
      }
      else if ( strncmp(str, _signal_name, 127) == 0 || strlen(_signal_name) == 0 )
      {
        /* nothing to do */
      }
      else
      {
        /* line starting with a letter, for example */
        fprintf(stderr, "Unknown file format.\n");
        return;
      }
    }
    retval = fscanf(minmaxfile, "%f", &_max_signal);
  }
  _max_signal += 1e-6; /* to avoid clipping for the max value */

  retval = fscanf(minmaxfile, "%f", &_min_signal);
  while (retval != 1)
  {
    if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}

    if (fscanf(minmaxfile, "%c", &c) == 1)
    {
      if (c=='#')
      {
        /* comment line, ignore the whole line */
        /* printf("comment line :\n"); */
        do {retval = fscanf(minmaxfile, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
        if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}
      }
      else
      {
        /* line strating with a letter, for example */
        fprintf(stderr, "Unknown file format.\n");
        return;
      }
    }
    retval = fscanf(minmaxfile, "%f", &_min_signal);
  }

  fclose(minmaxfile);
}



void read_file_line()
{
  float signal;
  int i, col, ind;
  int colorindex;
  int minclip = 0, maxclip = 0;
  int retval;
  char c;
  char str[128];


  // Read the header...
  // TODO: refactor this !

  /*  printf("reading line %d\n", line); */
  retval = fscanf(trajFile, "%f", &_time);

  while (retval != 1)
  {
    if (retval == EOF) {exit(0);}

    if (fscanf(trajFile, "%c", &c) == 1)
    {
      /* putchar(c); */
      switch(c)
      {
        case '#' :
          /* comment line, ignore the whole line */
          /* printf("comment line :\n"); */
          break; 
        case '!' :
          /* param line specifying _worldsize */
          fscanf(trajFile, "%f", &_worldsize);
          break;
        case '$' :
          fscanf(trajFile, "%d = %s", &ind, str); 
          if ( _list_signals.size == 0 ) 
          {
            _first_signal_col = ind;
            /* printf("_first_signal_col = %d\n", _first_signal_col); */
          }
          _list_signals.size++;
          if ( _list_signals.size < MAX_NBR_SIGNAL )
          {
            strncpy(_list_signals.name[_list_signals.size-1], str, 127);
          }
          /* signal name e.g. $ 8 = SIGNAL_1 */
          if (!_signal_found)
          {
            if ( strncmp(str,_signal_name,127) == 0 )
            {
              _signal_index = ind; 
              _signal_found = 1;
            }
            if ( _signal_index == DEFAULT_SIGNAL )
            {
              _signal_name = (char*)realloc(_signal_name,(strlen(str) + 1)*sizeof(char));
              _signal_index = ind;
              sprintf(_signal_name, "%s", str);
              _signal_found = 1;
            }
          }
          break;
        default :  
          fprintf(stderr, "Unknown file format.\n");
          exit(EXIT_FAILURE);
      }

      /* go to next line or exit if EOF is reached */
      do {retval = fscanf(trajFile, "%c", &c); /* printf("%c", c); */ }  while ((c != '\n') && (retval != EOF));
      if (retval == EOF) {exit(0);}
      _line ++;
    }

    if ( ( retval = fscanf(trajFile, "%f", &_time) ) && _signal_index == SIGNAL_NOT_FOUND )
    {
      fprintf(stderr, "Signal %s not found.\n", _signal_name);
      exit(EXIT_FAILURE);
    }
  }
  // </TODO>


  retval = fscanf(trajFile, "%d", &_popsize);
  /* printf("time %f pop %d \n", _time, _popsize); */

  if (_mycells != NULL) free(_mycells);
  _mycells = (CellSphere *) malloc(_popsize * sizeof(CellSphere));

  for(i = 0; i < _popsize; i++)
  {
    fscanf(trajFile, "%d %f %f %f %f", &_mycells[i].id, &_mycells[i].position[0], &_mycells[i].position[1], &_mycells[i].position[2], &_mycells[i].radius);

    col = 8;
    while ( col < _first_signal_col )
    {
      retval = fscanf(trajFile, "%*s"); /* skip to next entry */
      if ( retval == EOF )
      {
        fprintf(stderr, "File contains only %d column, but needs %d.\n", col, _signal_index);
        exit( EXIT_FAILURE );
      }
      ++col;
    }

    while ( col < _signal_index )
    {
      retval = fscanf(trajFile, "%f", &signal); /* skip to next entry */
      if ( retval == EOF )
      {
        fprintf(stderr, "File contains only %d column, but needs %d.\n", col, _signal_index);
        exit( EXIT_FAILURE );
      }
      ++col;
    }
    fscanf(trajFile, "%f", &signal);

    /* compute color from signal */
    colorindex = (int) 64 * (signal - _min_signal) / (_max_signal - _min_signal);
    if (colorindex < 0) {colorindex = 0; minclip = 1;}
    else if (colorindex > 63) {colorindex = 63; maxclip = 1;}

    _mycells[i].color[0] = _colormap[colorindex][0];
    _mycells[i].color[1] = _colormap[colorindex][1];
    _mycells[i].color[2] = _colormap[colorindex][2];

    do {retval = fscanf(trajFile, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
    if (retval == EOF) {exit(0);}

    // TODO : tmp patch
    if (i < _popsize-1)
      fscanf(trajFile, "\n%f %d", &_time, &_popsize);
    // </TODO>
  }

  if (minclip) printf("Warning, at least one cell had min color clipping at t=%f.\n", _time);
  if (maxclip) printf("Warning, at least one cell had max color clipping at t=%f.\n", _time);

  _line++;
}


void animate(int value)
{
  read_file_line();
  glutPostRedisplay();
  if (_pause == 0) glutTimerFunc(1.0/_speed, animate, 0);
}


void init(void)
{
  init_color_map();
  read_min_max_signals();

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_SMOOTH);

   /* color of the light */
   GLfloat white_light[] = {1.0, 1.0, 1.0, 1.0};
   glLightfv(GL_LIGHT0, GL_AMBIENT, white_light);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
   glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);

   /* position of the light */
   GLfloat light_position[] = {0.0, 0.0, 0.0, 1.0}; /* will move with the viewpoint */
   /* w=0 thus parallel rays like the sun, called directional light source */
   /* The direction (along z axis here) is transformed by the current modelview matrix */
   /* There is no attenuation for a directional light source */
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   /* lighting model */
   GLfloat global_ambient_light[] = {0.3, 0.3, 0.3, 1.0}; /* to see objects even if no light source */
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient_light);
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE); /* infinite viewpoint */
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE); /* back faces are inside the spheres, never seen */

   /* material for the objects */
   GLfloat mat_specular[] = {0.0, 0.0, 0.0, 1.0};
   GLfloat mat_ambient_refl[] = {0.2, 0.2, 0.2, 1.0};
   GLfloat mat_shininess[] = {0.0};
   GLfloat mat_diffuse_color[] = {0.5, 0.5, 1.0};
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient_refl);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_color);
   glColorMaterial(GL_FRONT, GL_DIFFUSE); /* now glColor changes diffusion color */

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_RESCALE_NORMAL);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_COLOR_MATERIAL);

   /* Set up nearClip and farClip so that ( farClip - nearClip ) > maxObjectSize, */
   /* and determine the viewing distance (adjust for zooming) */
   _nearClip = 0.5*_worldsize;
   _farClip = _nearClip + 4.0*_worldsize;
   resetView();

   /* display list to store the goemetry of a sphere */
   _SphereDList = glGenLists(1);
   glNewList(_SphereDList, GL_COMPILE);
   glutSolidSphere(1.0, 10, 10);
   glEndList();
}

void display(void)
{
  int i;
  GLfloat NicheColor[] = {0.5,0.5,0.5,0.1};
  char mystring[128];
  GLubyte raster_rect[16] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /*************** Draw the text and legends ******************/
  sprintf(mystring, "t = %.2f,  nb cells = %d ", _time, _popsize);
  glColor3f (1.0, 1.0, 1.0);
  glWindowPos2i(15, 15); /* also sets current raster color to white */
  for (i = 0; i < strlen(mystring); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);

  glColor3fv(_colormap[0]);
  glWindowPos2i(250, 10);
  for(i = 0; i < 64; i++)
  {
    glColor3fv(_colormap[i]);
    glWindowPos2i(250 + 3*i, 10);
    glBitmap(3, 16, 0.0, 0.0, 3, 0, raster_rect);
  }

  glColor3f (1.0, 1.0, 1.0);
  sprintf(mystring, "%.2f", _min_signal);
  glWindowPos2i(250, 30); /* also sets current raster color to white */
  for (i = 0; i < strlen(mystring); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);

  glColor3f (1.0, 1.0, 1.0);
  sprintf(mystring, "%s", _signal_name);
  glWindowPos2i(280, 30); /* also sets current raster color to white */
  if ( strlen(mystring) > MAX_SIG_STRLEN ) /* print first three chars ... last (MAX_SIG_STRLEN - 6) */
  {
    for (i = 0; i < 3; i++)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);
    for (i = 3; i < 6; i++)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, '.');
    for (i = strlen(mystring) - (MAX_SIG_STRLEN - 6); i < strlen(mystring); i++)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);
  }
  else
  {
    for (i = 0; i < strlen(mystring); i++)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);
  }

  glColor3f (1.0, 1.0, 1.0);
  sprintf(mystring, "%.2f", _max_signal);
  glWindowPos2i(420, 30); /* also sets current raster color to white */
  for (i = 0; i < strlen(mystring); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystring[i]);


  /*********************** Set the view **************************/

  glLoadIdentity ();             /* clear the matrix */
  /* viewing transformation  */
  /* gluLookAt (-1.5*_worldsize, -1.5*_worldsize, 1.5*_worldsize, _worldsize, _worldsize, 0.0, _worldsize, _worldsize, 1.5*_worldsize); */
  polarView( _distance, _azimAngle, _incAngle, _twistAngle );

  /* current matrix is V */

  /*********************** Draw the box **************************/
  glColor3f (1.0, 1.0, 1.0); // White
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Mesh

  float world_size[3] = {_worldsize, _worldsize, _worldsize};
  float world_margins[3] = {2.0, 2.0, 0.0};
  float total_world_size[3] = {world_size[0] + 2 * world_margins[0],
                               world_size[1] + 2 * world_margins[1],
                               world_size[2] + 2 * world_margins[2]};

  // Draw outer world borders
  glBegin(GL_QUAD_STRIP);
  glVertex3d(0, 0, 0);
  glVertex3d(0, total_world_size[1], 0);
  glVertex3d(total_world_size[0], 0, 0);
  glVertex3d(total_world_size[0], total_world_size[1], 0);
  glVertex3d(total_world_size[0], 0, total_world_size[2]);
  glVertex3d(total_world_size[0], total_world_size[1], total_world_size[2]);
  glVertex3d(0, 0, total_world_size[2]);
  glVertex3d(0, total_world_size[1], total_world_size[2]);
  glVertex3d(0, 0, 0);
  glVertex3d(0, total_world_size[1], 0);
  glEnd();

  // Draw outer world borders
  glBegin(GL_QUAD_STRIP);
  glVertex3d(world_margins[0], world_margins[1], world_margins[2]);
  glVertex3d(world_margins[0], world_size[1], world_margins[2]);
  glVertex3d(world_size[0], world_margins[1], world_margins[2]);
  glVertex3d(world_size[0], world_size[1], world_margins[2]);
  glVertex3d(world_size[0], world_margins[1], world_size[2]);
  glVertex3d(world_size[0], world_size[1], world_size[2]);
  glVertex3d(world_margins[0], world_margins[1], world_size[2]);
  glVertex3d(world_margins[0], world_size[1], world_size[2]);
  glVertex3d(world_margins[0], world_margins[1], world_margins[2]);
  glVertex3d(world_margins[0], world_size[1], world_margins[2]);
  glEnd();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Back to fill mode



  /******************* Draw the cells       *********************/
  for (i = 0; i < _popsize; i++)
  {
    glPushMatrix(); /* save V */
    if (_mycells[i].position[2] > 0) {
            glColor3fv(_mycells[i].color);
        } else {
            glColor3fv(NicheColor);
        }
    glTranslatef(_mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2]); /* current matrix is VT  */
    glScalef(_mycells[i].radius, _mycells[i].radius, _mycells[i].radius ); /* current matrix is VTS */
    glCallList(_SphereDList);
    glPopMatrix(); /* current matrix is restored to V */
  }

  glutSwapBuffers();
}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  /*gluPerspective(30.0, (1.0*w)/h, 0.5*_worldsize, 5*_worldsize); */
  gluPerspective( 45.0f, (1.0*w)/h, _nearClip, _farClip );
  glMatrixMode (GL_MODELVIEW);
}


void resetView()
{
  _distance = _nearClip + (_farClip - _nearClip) / 2.0;
  _twistAngle = 0.0;	/* rotation of viewing volume (camera) */
  _incAngle = 85.0;
  _azimAngle = 30.0;
}


void polarView( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist)
{
  /* printf(" incidence %f azimuth %f\n", incidence, azimuth); */
  glTranslatef( -_worldsize/2.0, -_worldsize/2.0, -distance);
  glRotatef( -twist, 0.0f, 0.0f, 1.0);
  glRotatef( -incidence, 1.0f, 0.0f, 0.0);
  glRotatef( -azimuth, 0.0f, 0.0f, 1.0);
}


void keyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'r': /* r or R : reset viewpoint */
      resetView();
      glutPostRedisplay();
      break;
    case 'R':
      resetView();
      glutPostRedisplay();
      break;
    case 32: /* space bar = pause the animation */
      if (_pause == 0)
      {_pause = 1;}
      else
      {
        _pause = 0;
        glutTimerFunc(0,animate,0);
      }
      break;
    case 'j': /* decrease max signal */
      _max_signal -= 0.1*(_max_signal - _min_signal); 
      break;
    case 'k': /* increase max signal */
      _max_signal += 0.1*(_max_signal - _min_signal); 
    case 'h': /* decrease min signal */
      _min_signal -= 0.1*(_max_signal - _min_signal); 
      break;
    case 'l': /* increase min signal */
      _min_signal += 0.1*(_max_signal - _min_signal); 
      break;
    case ']': /* display next signal */
      next_signal();
      read_min_max_signals();
      break;
    case '[':
      previous_signal();
      read_min_max_signals();
      break;
    case 27: /* escape to quit */
      exit(0);
  }
}

void next_signal()
{
      _signal_index++;
      _signal_index = ( ( _signal_index - _first_signal_col ) % _list_signals.size ) + _first_signal_col;
      _signal_name = (char*)realloc(_signal_name,(strlen(_list_signals.name[_signal_index - _first_signal_col]) + 1)*sizeof(char));
      sprintf(_signal_name, "%s", _list_signals.name[_signal_index - _first_signal_col]); 
}

void previous_signal()
{
      _signal_index--;
      _signal_index = ( ( _signal_index - _first_signal_col + _list_signals.size ) % _list_signals.size );
      _signal_index +=  _first_signal_col;
      _signal_name = (char*)realloc(_signal_name,(strlen(_list_signals.name[_signal_index - _first_signal_col]) + 1)*sizeof(char));
      sprintf(_signal_name, "%s", _list_signals.name[_signal_index - _first_signal_col]); 
}

GLvoid mouse( GLint button, GLint state, GLint x, GLint y )
{
  if (state == GLUT_DOWN)
  {
    switch (button)
    {
      case GLUT_LEFT_BUTTON:
        _action = MOVE_EYE;
        break;
      case GLUT_MIDDLE_BUTTON:
        _action = TWIST_EYE;
        break;
      case GLUT_RIGHT_BUTTON:
        _action = ZOOM;
        break;
    }

    /* Update the saved mouse position */
    _xStart = x;
    _yStart = y;
  }
  else
  {
    _action = MOVE_NONE;
  }
}

GLvoid motion( GLint x, GLint y )
{
  switch (_action)
  {
    case MOVE_EYE:
      /* Adjust the eye position based on the mouse position */
      _azimAngle += (GLdouble) (x - _xStart);
      _incAngle -= (GLdouble) (y - _yStart);
      break;
    case TWIST_EYE:
      /* Adjust the eye twist based on the mouse position */
      _twistAngle = fmod(_twistAngle+(x - _xStart), 360.0);
      break;
    case ZOOM:
      /* Adjust the eye distance based on the mouse position */
      _distance -= (GLdouble) (y - _yStart)/10.0;
      break;
  }

  /* Update the stored mouse position for later use */
  _xStart = x;
  _yStart = y;

  glutPostRedisplay();
}



void print_help( char* prog_name )
{
  printf( "\n************* Simuscale - Viewing program ************* \n\n" );
  printf( "\n\
Usage : %s -h\n\
   or : %s [options]\n\n\
\t-h or --help      : Display this screen\n\n\
Options (i : integer, d : double, s : string) :\n\n\
\t-f or --file   s  : View the simulation stored in file named s (default: trajectory.txt)\n\
\t-s or --speed  d  : Set the playing speed to d, in frames per ms (default: 0.02)\n\
\t-i or --signal s  : Set the signal to display (default is column 8 in file)\n\n\
\t-p or --pause     : Launch in pause mode on (resume with space bar)\n\n\
If an option is not set, the program uses the default value for this parameter.\n\n\
While the animation is playing, you can use the following commands:\n\n\
\tleft button of the mouse to change the viewpoint,\n\
\tright button to zoom in or zoom out,\n\
\tR to reset the perspective.\n\
\tspace bar to pause the animation or to resume.\n\
\tj and k keys to increase or decrease max signal intensity.\n\
\th and l keys to increase or decrease min signal intensity.\n\
\t] and [ keys to switch to next or previous signal.\n\
\tPress ESC to quit.\n\n",\
   prog_name+2, prog_name+2 );
}

