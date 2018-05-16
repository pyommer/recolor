/**
 * Program: recolor.cpp
 * Author:  Phillip Yommer
 * Email:   pyommer@g.clemson.edu
 * Course:  CPSC 4040
 * Title:   Computer Graphics Images
 * Problem: Final Programming Project
 * Title:   Recoloring Images for Vision Deficient Viewers
 * Date:    12/04/2017
 *
 *
 * Purpose:	Color mapping to a limited palette. The main goal is to allow a
 *          user to select a default color vision deficiency (Protanomaly,
 *          Deutanomaly, or Tritanomaly) and have the program filter the colors
 *          of an image to the visible colors for the deficiency, creating an
 *          output image clearly visible in detail to a viewer with the
 *          selected color deficiency.
 *
 * Usage:   make recolor ; ./recolor <input.img> [<output.img>] [-OPTIONS]
 *
 */

Usage:      ./recolor <input.img> [<output.img>] [-OPTIONS]

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


Bugs:


Notes:

    1.  Pixels are converted from RGB space to XYZ space.
    2.  XYZ pixels are transformed via some color deficiency-correcting
        transformation matrix.
    3.  XYZ color-corrected pixels are converted back to RGB space.
    4.  Alpha values for RGBA pixels are maintained.
