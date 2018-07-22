# recolor

Files:

    - recolor.cpp   - the source code file of the recolor program
    - recolor.h     - the library file with constants and functions for recolor
    - recolor.pdf   - the project writeup describing the recolor project
    - recolor.html  - HTML file breifly describing the project with examples
    - test.jpg      - a test image to use with recolor
                    - https://cdn-images-1.medium.com/max/1280/1*c4y7KbWvyQFSz9u41ephVw.jpeg
    - Makefile      - the makefile to use to build the recolor project
    - README.md     - this readme file for the recolor project

# Description:

    Program: recolor.cpp
    Author:  Phillip Yommer
    Course:  CPSC 4040
    Title:   Computer Graphics Images
    Problem: Final Programming Project
    Title:   Recoloring Images for Vision Deficient Viewers
    Date:    12/04/2017

Purpose:

	Color mapping to a limited palette. The main goal is to allow a
	user to select a default color vision deficiency (Protanomaly,
	Deutanomaly, or Tritanomaly) and have the program filter the colors
	of an image to the visible colors for the deficiency, creating an
	output image clearly visible in detail to a viewer with the
	selected color deficiency.

    For more detail see the Final Report at
    https://github.com/pyommer/recolor/blob/master/recolor.pdf.
    Example images are included in the Results section of the Final Report.

Compile:

    make recolor

Usage:

    ./recolor <input.img> [<output.img>] [-OPTIONS]

Arguments:

    "input.img"  - the required filename of the input image to be recolored.
    "output.img" - the optional filename to write the recolored output image to.
    -OPTIONS     - the optional recoloring mode for the initial image recolor.

Interface:

    'a' - recolor the image for Protanomaly
    'b' - recolor the image for Deutanomaly
    'c' - recolor the image for Tritanomaly
    'q' - quit


Assumptions:

    1.  Images are stored in either RGB or RGBA color space.
    2.  The difference between RGB and sRGB space is negligible.
    3.  The difference between XYZ and xyY space is negligible.

Notes:

    1.  Pixels are converted from RGB space to XYZ space.
    2.  XYZ pixels are transformed via some color deficiency-correcting
        transformation matrix.
    3.  XYZ color-corrected pixels are converted back to RGB space.
    4.  Alpha values for RGBA pixels are maintained.
