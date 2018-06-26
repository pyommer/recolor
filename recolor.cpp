//recolor.cpp
#include <OpenImageIO/imageio.h>
#include <GL/glut.h>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
using namespace std;
OIIO_NAMESPACE_USING

#include "recolor.h"

// main program
int main(int argc, char *argv[])
{
    // verify command line
    if (argc < 2)
    {
        cout << "usage: ./warper <input.img> [<output.img>] [-OPTIONS]" << endl;
        exit(1);
    }

    // determine recoloring mode
    char mode = (argc>3 && argv[3][0]=='-') ? argv[3][1] : PROTAN1;

    // read the input image
    const pixmap_t pixmap = readImage((const char*) argv[1]);

    // perform the recoloring transformation on the image
    pixmap_t result = recolorImage(mode, pixmap);

    // write output image
    if (argc > 2)
        writeImage((const char*) argv[2], (const pixmap_t) result);

    // initialize display structure
    display.pixmap = result;
    display.window[0] = (int) display.pixmap[0].size();
    display.window[1] = (int) display.pixmap.size();

    // destroy no-longer needed data
    destroy(result);
    destroy(pixmap);

    // start up GLUT
    glutInit(&argc,argv);

    // create the graphics window, giving width, height, and title text
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_RGB);
    glutInitWindowSize(display.window[0], display.window[1]);
    glutCreateWindow("Recolor");

    // display result image transformation
    displayImage(display.pixmap);
    glFlush();                          // flush OpenGL pipeline

    // set up the callback routines
    glutDisplayFunc(handleDisplay);     // display update callback
    glutKeyboardFunc(handleKey);        // keyboard key press callback

    // specify window clear (background) color to be opaque black
    glClearColor(0,0,0,1);

    // Enter GLUT's event loop
    glutMainLoop();
    return 0;
}

