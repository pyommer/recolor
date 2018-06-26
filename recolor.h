//recolor.h
/*-------------------------------PREPROCESSOR---------------------------------*/
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

// pi, rounded to (10^-10)
#ifndef PI
#define PI      3.1415926536
#endif

// rounding factor,  -SIGMA < x < SIGMA  =>  x == 0.0
#ifndef SIGMA
#define SIGMA   0.001
#endif

// recoloring color vision deficiency modes
#define PROTAN1 'a'
#define DEUTAN1 'b'
#define TRITAN1 'c'

/*---------------------------STRUCTURE-DEFINITIONS----------------------------*/

// defined types for storing channel values
typedef vector<double>   vector_t;
typedef vector<vector_t> matrix_t;
typedef vector<matrix_t> pixmap_t;

// display structure to hold OpenGL data
struct Display
{
    int window[2];          // window dimensions (win[0]=width, win[1]=height)
    pixmap_t pixmap;        // the pixmap defining the image to be displayed
} display;              // the global display structure


/*----------------------------FUNCTION-PROTOTYPES-----------------------------*/

// image file i/o functions
matrix_t read_input();
const pixmap_t readImage(const char *name);
void writeImage(const char *name, const pixmap_t pixmap);
void displayImage(const pixmap_t pixmap);

// image recoloring
pixmap_t recolorImage(char mode, const pixmap_t pixmap);
matrix_t recolorPixel(char mode, vector_t pixel);
matrix_t getFilter(char mode);

// color spaces
matrix_t toRGB(vector_t pixel);
matrix_t toRGB(matrix_t pixel);
matrix_t toXYZ(vector_t pixel);
matrix_t toXYZ(matrix_t pixel);

// OpenGL handlers
void handleDisplay();
void handleKey(unsigned char key, int x, int y);

// destructors
void destroy(pixmap_t pixmap);
void destroy(matrix_t matrix);
void destroy(vector_t vector);

// matrix statistics
double normalize(double max, double min, double value);
double normalize(double value, matrix_t matrix);
matrix_t normalize(matrix_t matrix);
double getMax(pixmap_t pixmap);
double getMin(pixmap_t pixmap);
double getMax(matrix_t matrix);
double getMin(matrix_t matrix);
double getMax(vector_t matrix);
double getMin(vector_t matrix);
double getMax(double a, double b);
double getMin(double a, double b);

// printers
void printMatrix(const char *name, matrix_t matrix);
void printMatrix(const char *name, vector_t matrix);
void printMatrix(matrix_t matrix);
void printMatrix(vector_t matrix);

// matrix initialization and equality
matrix_t identity(int n);
matrix_t zero(int n);
matrix_t zero(int m, int n);
matrix_t set(matrix_t matrix);
matrix_t set(vector_t matrix);
matrix_t set(double x, double y, double z);
matrix_t set(double x, double y);
matrix_t set(int x, int y, int z);
matrix_t set(int x, int y);
bool equals(matrix_t a, matrix_t b);

// matrix arithmetic
matrix_t divide(double value, matrix_t matrix);
matrix_t multiply(double value, matrix_t matrix);
matrix_t add(double value, matrix_t matrix);
matrix_t add(matrix_t a, matrix_t b);

// matrix operations
matrix_t cross(matrix_t a, matrix_t b);
matrix_t cross(matrix_t a, vector_t b);
matrix_t transpose(matrix_t matrix);
matrix_t transpose(vector_t matrix);
double det(matrix_t matrix);
matrix_t invert(matrix_t matrix);
matrix_t cofactor(int r, int c, matrix_t matrix);
matrix_t adjoint(matrix_t matrix);
double sum(matrix_t matrix);

/*---------------------------FUNCTION-DEFINITIONS-----------------------------*/

/* -- pixmap recoloring -- */

// returns the recolored pixmap for the specified mode
pixmap_t recolorImage(char mode, const pixmap_t pixmap)
{
    int c = (int) pixmap[0][0].size();
    pixmap_t result(pixmap.size(), matrix_t(pixmap[0].size(), vector_t(pixmap[0][0].size())));
    for (int i=0; i<(int)pixmap.size(); i++)
    {
        for (int j=0; j<(int)pixmap[0].size(); j++)
        {
            // determine the recolored pixel values
            matrix_t temp = recolorPixel(mode, pixmap[i][j]);

            // set the result pixel to the recolored pixel values
            for (int k=0; k<3; k++)
                result[i][j][k] = temp[k][0];

            // determine and set alpha channel
            if (c > 3)
            {
                for (int k=3; k<c; k++)
                    result[i][j][k] = pixmap[i][j][k];
            }
        }
    }
    return result;
}

// returns the recolored pixel for the specified mode
matrix_t recolorPixel(char mode, vector_t pixel)
{
    // verify pixel channels
    if (pixel.size() < 3)
    {
        printf("recolor pixel error: too few pixel channels");
        printf(",\t|pixel| = %d\n", (int) pixel.size());
        exit(1);
    }

    // remove alpha channel
    matrix_t p = zero(3,1);
    for (int i=0; i<3; i++)
        p[i][0] = pixel[i];

    // recolor the pixel
    return toRGB(cross(getFilter(mode), toXYZ(p)));
}

// returns the recoloring transformation matrix for the specified mode
matrix_t getFilter(char mode)
{
    switch (mode)
    {
        case PROTAN1:
        {
            matrix_t protan(3, vector_t(3));
            protan[0][0] = 0.567;
            protan[0][1] = 0.433;
            protan[0][2] = 0.0;
            protan[1][0] = 0.558;
            protan[1][1] = 0.442;
            protan[1][2] = 0.0;
            protan[2][0] = 0.0;
            protan[2][1] = 0.242;
            protan[2][2] = 0.478;
            return set(protan);
        }
        case DEUTAN1:
        {
            matrix_t deutan(3, vector_t(3));
            deutan[0][0] = 0.625;
            deutan[0][1] = 0.375;
            deutan[0][2] = 0.0;
            deutan[1][0] = 0.7;
            deutan[1][1] = 0.3;
            deutan[1][2] = 0.0;
            deutan[2][0] = 0.0;
            deutan[2][1] = 0.3;
            deutan[2][2] = 0.7;
            return set(deutan);
        }
        case TRITAN1:
        {
            matrix_t tritan(3, vector_t(3));
            tritan[0][0] = 0.95;
            tritan[0][1] = 0.05;
            tritan[0][2] = 0.0;
            tritan[1][0] = 0.0;
            tritan[1][1] = 0.433;
            tritan[1][2] = 0.567;
            tritan[2][0] = 0.0;
            tritan[2][1] = 0.475;
            tritan[2][2] = 0.525;
            return set(tritan);
        }
    }
    return identity(3);
}

/* -- color space conversions -- */

// returns the RGB pixel matrix for the CIE XYZ coordinate vector
matrix_t toRGB(vector_t pixel)
{
    return toRGB(transpose(pixel));
}
// returns the RGB pixel matrix for the CIE XYZ coordinate matrix
matrix_t toRGB(matrix_t pixel)
{
    // verify pixel is [3x1]
    if (!(pixel.size()==3 && pixel[0].size()==1) || (pixel.size()==1 && pixel[0].size()==3))
    {
        printf("pixel to RGB error: wrong number of pixel dimensions");
        printf(",\t|pixel| = [%dx%d]\n", (int) pixel.size(), (int) pixel[0].size());
        exit(1);
    }
    if (pixel.size()==1 && pixel[0].size()==3)
        pixel = transpose(pixel);

    // build the XYZ-to-RGB transformation matrix
    matrix_t rgb = zero(3,3);
    rgb[0][0] = 0.4124;
    rgb[0][1] = 0.3576;
    rgb[0][2] = 0.1805;
    rgb[1][0] = 0.2126;
    rgb[1][1] = 0.7152;
    rgb[1][2] = 0.0722;
    rgb[2][0] = 0.0193;
    rgb[2][1] = 0.1192;
    rgb[2][2] = 0.9505;
    return cross(rgb,pixel);
}

// returns the CIE XYZ coordinates for the RGB pixel vector
matrix_t toXYZ(vector_t pixel)
{
    return toXYZ(transpose(pixel));
}
// returns the CIE XYZ coordinates for the RGB pixel matrix
matrix_t toXYZ(matrix_t pixel)
{
    // verify pixel is [3x1]
    if (!(pixel.size()==3 && pixel[0].size()==1) || (pixel.size()==1 && pixel[0].size()==3))
    {
        printf("pixel to XYZ error: wrong number of pixel dimensions");
        printf(",\t|pixel| = [%dx%d]\n", (int) pixel.size(), (int) pixel[0].size());
        exit(1);
    }
    if (pixel.size()==1 && pixel[0].size()==3)
        pixel = transpose(pixel);

    // build the RGB-to-XYZ transformation matrix
    matrix_t xyz = zero(3,3);
    xyz[0][0] = 3.2406;
    xyz[0][1] = -1.5372;
    xyz[0][2] = -0.4986;
    xyz[1][0] = -0.9689;
    xyz[1][1] = 1.8758;
    xyz[1][2] = 0.0415;
    xyz[2][0] = 0.0557;
    xyz[2][1] = -0.0204;
    xyz[2][2] = 1.057;
    return cross(xyz,pixel);
}

/* -- i/o functions -- */

// returns an image matrix of the normalized values pixels from the image file
const pixmap_t readImage(const char *name)
{
    // open input image file
    ImageInput *in = ImageInput::open(name);
    if (!in)
    {
        cout << "Could not open image " << name;
        cerr << ",\terror = " << geterror() << endl;
        exit(1);
    }

    // get image spec data
    int h = in->spec().height;
    int w = in->spec().width;
    int c = in->spec().nchannels;

    // read the pixel data from the input image to the local pixmap
    unsigned char local_pixmap[h*w*c];
    if (!in->read_image(TypeDesc::UINT8,local_pixmap))
    {
        cout << "Could not read image from " << name;
        cerr << ", error = " << geterror() << endl;
        delete in;
        exit(1);
    }

    // initialize the image matrix to normalized values for the input image
    pixmap_t image(h, matrix_t(w, vector_t(c)));
    for (int i=0; i<h; i++)
        for (int j=0; j<w; j++)
            for (int k=0; k<c; k++)
                image[i][j][k] = local_pixmap[i*w*c+j*c+k]/255.0;   // [0.0,1.0]

    // close input image
    if (!in->close())
    {
        cout << "Could not close " << name;
        cerr << ",\terror = " << geterror() << endl;
        delete in;
        exit(1);
    }
    delete in;
    printf("image file %s read success, [hxwxc] = [%dx%dx%d]\n", name, h,w,c);
    return image;
}

// writes the image matrix to the specified image file
void writeImage(const char *name, const pixmap_t pixmap)
{
    // create output image file
    ImageOutput *out = ImageOutput::create(name);
    if (!out)
    {
        cout<<"Could not create output image for " << name;
        cerr << ", error = " << geterror() << endl;
        exit(1);
    }

    // initialize image spec data
    int h = (int) pixmap.size();
    int w = (int) pixmap[0].size();
    int c = (int) pixmap[0][0].size();

    // initialize a local pixmap for the output image matrix
    int i=0, j=0, k=0;
    unsigned char local_pixmap[h*w*c];

    // set the local pixmap values from the output image matrix
    for (int l=0; l<h*w*c; l++)
    {
        local_pixmap[l] = 255*pixmap[i][j][k++];      // [0.0,1.0] -> [0,255]

        // pixel stored locally
        if (k == c)
        {
            j++;    // increment column
            k=0;    // reset channel
        }
        // scanline stored locally
        if (j == w)
        {
            i++;    // increment row
            j=0;    // reset column
        }
        // image stored locally
        if (i == h)
            break;
    }

    // write image spec
    ImageSpec spec(w,h,c,TypeDesc::UINT8);
    if (!out->open(name,spec))
    {
        cout << "Could not open " << name;
        cerr << ", error = " << geterror() << endl;
        delete out;
        exit(1);
    }

    // write the local pixmap to output image file
    if(!out->write_image(TypeDesc::UINT8,local_pixmap))
    {
        cout << "Could not write image to " << name;
        cerr << ", error = " << geterror() << endl;
        delete out;
        exit(1);
    }

    // close output file
    if (!out->close())
    {
        cout << "Could not close " << name;
        cerr << ", error = " << geterror() << endl;
        delete out;
        exit(1);
    }
    delete out;
    printf("image file %s write success, [hxwxc] = [%dx%dx%d]\n", name, h,w,c);
}

// draws the image matrix to the OpenGL viewport
void displayImage(pixmap_t pixmap)
{
    // initialize image spec data
    int h = (int) pixmap.size();
    int w = (int) pixmap[0].size();
    int c = (int) pixmap[0][0].size();
    int pixform;

    // verify pixmap dimensions
    if (!(h>0 && w>0))
    {
        cout << "display error: invalid pixmap dimensions" << endl;
        return;
    }

    // determine the output pixform from the image size
    switch (c)
    {
        case 1:
            pixform = GL_LUMINANCE;
            break;
        case 3:
            pixform = GL_RGB;
            break;
        case 4:
            pixform = GL_RGBA;
            break;
        default:
            cout << "display error: invalid pixel format, nchan = "<< c << endl;
            return;
    }

    // initialize the local pixmap from the input image
    int i = h-1;    // start at last row of image
    int j = 0;
    int k = 0;
    unsigned char local_pixmap[h*w*c];

    // set the local pixmap values from the output image matrix
    for (int l=0; l<h*w*c; l++)
    {
        local_pixmap[l] = 255*pixmap[i][j][k++];

        // pixel stored locally
        if (k == c)
        {
            j++;    // increment column
            k=0;    // reset channel
        }
        // scanline stored locally
        if (j == w)
        {
            i--;    // decrement row
            j=0;    // reset column
        }
        // image stored locally, done when i == -1
        if (i == -1)
        {
            if (l != h*w*c-1)
                cout << "display image error: local pixmap size" << endl;
            break;
        }
    }

    // draw the local pixmap to the viewport
    glPixelZoom(1.0,1.0);
    //glRasterPos2i(0,0);           // start display at the lower left corner
    /* when commented-out, the image seems to be displayed correctly.
     * when raster is enabled, the image looks to be drawn from the origin,
     * (h/2,w/2), instead of the lower-left corner.
     */
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glDrawPixels(w, h, pixform, GL_UNSIGNED_BYTE, local_pixmap);
    printf("output image draw success\n");
}

/* -- OpenGL handlers -- */

// handles displaying an image to the window
void handleDisplay()
{
    glClear(GL_COLOR_BUFFER_BIT);   // clear window to background color
    displayImage(display.pixmap);    // draw the image
    glFlush();                      // flush OpenGL pipeline to the viewport
}

// handles keypresses - display image until 'q' or ESC pressed
void handleKey(unsigned char key, int x, int y)
{
    switch(key)
    {
        case PROTAN1:
        case DEUTAN1:
        case TRITAN1:
            display.pixmap = recolorImage((char) key, display.pixmap);
            break;
        case 'q':       // q or ESC - quit
        case 'Q':
        case 27:
          destroy(display.pixmap);
          cout << "exit." << endl;
          exit(0);
        default:        // not a valid key -- just ignore it
            return;
    }
}


/* -- matrix operations -- */

// returns the cross product matrix: Axb, ([3x3]x[3x1] = [3x1])
matrix_t cross(matrix_t a, vector_t b)
{
    return cross(a, transpose(b));
}
// returns the cross product matrix: AxB, ([3x3]x[3x3] = [3x3])
matrix_t cross(matrix_t a, matrix_t b)
{
    // verify matrix dimensions, a[m1xn] x b[nxm2]
    if (a[0].size() != b.size())
    {
        // transpose if b is vector-matrix
        if (a[0].size() == b[0].size())
            return cross(a,transpose(b));
        printf("cross product error: matrix size mismatch");
        printf(",\t|A[0]| = %d, |B| = %d\n", (int) a[0].size(), (int) b.size());
        printMatrix("A", a);
        printMatrix("B", b);
        exit(1);
    }

    // compute the cross product
    matrix_t result = zero((int) a.size(), (int) b[0].size());
    for (int i=0; i<(int) a.size(); i++)
        for (int j=0; j<(int) b[0].size(); j++)
            for (int k=0; k<(int) a[0].size(); k++)
                result[i][j] += a[i][k]*b[k][j];
    return result;
}

// returns the transpose of the matrix, [mxn] -> [nxm]
matrix_t transpose(matrix_t matrix)
{
    matrix_t result = zero((int) matrix[0].size(), (int) matrix.size());
    for (int i=0; i<(int) matrix[0].size(); i++)
        for (int j=0; j<(int) matrix.size(); j++)
            result[i][j] = matrix[j][i];
    return result;
}
// returns the transpose of the vector, [1xn] -> [nx1]
matrix_t transpose(vector_t matrix)
{
    matrix_t result = zero((int) matrix.size(), 1);
    for (int i=0; i<(int) matrix.size(); i++)
        result[i][0] = matrix[i];
    return result;
}

// returns the determinant of an nxn matrix
double det(matrix_t matrix)
{
    // verify matrix is square
    if (matrix.size()!=matrix[0].size())
    {
        printf("matrix determinant error: non-square matrix");
        printf(",\t|M| = [%dx%d]\n", (int)matrix.size(), (int)matrix[0].size());
        printMatrix("M", matrix);
        exit(1);
    }

    // verify matrix size
    if ((int) matrix.size()<2 || (int) matrix[0].size()<2)
    {

        printf("matrix determinant error: matrix is less than [2x2]");
        printf(",\t|M| = [%dx%d]\n", (int)matrix.size(), (int)matrix[0].size());
        printMatrix("M", matrix);
        exit(1);
    }

    // return [2x2] determinant, |{{a,b},{c,d}}| = (ad-bc)
    if ((int) matrix.size()==2 && (int) matrix[0].size()==2)
        return ((matrix[0][0]*matrix[1][1]) - (matrix[0][1]*matrix[1][0]));

    // recursively call det() for any [nxn], n>2
    double sum = 0.0;
    for (int i=0; i<(int) matrix.size(); i++)
        sum += matrix[0][i]*det(cofactor(0,i,matrix));
    return sum;
}

// returns the inverse of the matrix
matrix_t invert(matrix_t matrix)
{
    // verify matrix is square and det != 0
    if (matrix.size()!=matrix[0].size() || det(matrix)==0.0)
    {
        printf("invert matrix error: non-invertible matrix");
        printf(",\t|M| = [%dx%d]\n", (int)matrix.size(), (int)matrix[0].size());
        printMatrix("M", matrix);
        exit(1);
    }
    return divide(det(matrix), adjoint(matrix));    // M^-1 = adj(M)/det(M)
}

// returns the cofactor matrix for the element at position r,c in the matrix
matrix_t cofactor(int r, int c, matrix_t matrix)
{
    // verify row and column
    if ((r<0 || r>=(int) matrix.size()) || (c<0 || c>=(int) matrix[0].size()))
    {
        printf("matrix cofactor error: row or column out of matrix bounds");
        printf(",\tr = %d, c = %d\n", r, c);
        printMatrix("M", matrix);
        exit(1);
    }

    // determine cofactor matrix
    matrix_t result = zero((int) matrix.size()-1, (int) matrix[0].size()-1);
    for (int i=1; i<(int) matrix.size(); i++)
        for (int j=1; j<(int) matrix.size(); j++)
            result[i-1][j-1] = matrix[(r+i)%matrix.size()][(c+j)%matrix.size()];
    return result;
}

// returns the adjoint of the matrix
matrix_t adjoint(matrix_t matrix)
{
    // verify matrix is square
    if (matrix.size()!=matrix[0].size())
    {
        printf("adjoint matrix error: non-square matrix");
        printf(",\t|M| = [%dx%d]\n", (int)matrix.size(), (int)matrix[0].size());
        printMatrix("M", matrix);
        exit(1);
    }

    // compute the determinant of the cofactor for each element
    matrix_t result = zero((int) matrix.size(), (int) matrix[0].size());
    for (int i=0; i<(int) matrix.size(); i++)
        for (int j=0; j<(int) matrix[0].size(); j++)
            result[i][j] = det(cofactor(i,j,matrix));

    // determine sign multiplier for each element position
    /*
    for (int i=0; i<(int) result.size(); i++)
        for (int j=0; j<(int) result[0].size(); j++)
            result[i][j] *= ((i+j)%2 == 1) ? -1.0 : 1.0;
    */
    return transpose(result);
}

// returns the sum of the elements in the matrix
double sum(matrix_t matrix)
{
    double result = 0.0;
    for (int i=0; i<(int)matrix.size(); i++)
        for (int j=0; j<(int)matrix[0].size(); j++)
            result += matrix[i][j];
    return result;
}

/* -- matrix arithmetic -- */

// returns the matrix divided by the quotient value
matrix_t divide(double value, matrix_t matrix)
{
    // verify the value
    if (value == 0.0)
    {
        printf("divide matrix error: divide by 0\n");
        printMatrix("M", matrix);
        exit(1);
    }

    // divide the matrix by the value
    matrix_t result = set(matrix);
    for (int i=0; i<(int) matrix.size(); i++)
        for (int j=0; j<(int) matrix[0].size(); j++)
            result[i][j] /= value;
    return result;
}

// returns the product matrix matrix*value
matrix_t multiply(double value, matrix_t matrix)
{
    matrix_t result = set(matrix);
    for (int i=0; i<(int) matrix.size(); i++)
        for (int j=0; j<(int) matrix[0].size(); j++)
            result[i][j] *= value;
    return result;
}

// returns the sum of the two matrices
matrix_t add(matrix_t a, matrix_t b)
{
    // verify the value
    if (a.size()!=b.size() || a[0].size()!=b[0].size())
    {
        printf("add matrix error: mismatched matrix sizes");
        printf(",\t|A| = [%dx%d]", (int) a.size(), (int) a[0].size());
        printf(", |B| = [%dx%d]\n", (int) b.size(), (int) b[0].size());
        printMatrix("A", a);
        printMatrix("B", b);
        exit(1);
    }

    // add the matrices
    matrix_t result = set(a);
    for (int i=0; i<(int) a.size(); i++)
        for (int j=0; j<(int) a[0].size(); j++)
            result[i][j] += b[i][j];
    return result;
}
// returns the sum matrix matrix+value
matrix_t add(double value, matrix_t matrix)
{
    matrix_t result = set(matrix);
    for (int i=0; i<(int) matrix.size(); i++)
        for (int j=0; j<(int) matrix[0].size(); j++)
            result[i][j] += value;
    return result;
}

/* -- matrix initialization and equality -- */

// returns true if a[i][j]==b[i][j], for all i,j in a,b
bool equals(matrix_t a, matrix_t b)
{
    // verify dimensions
    if (a.size()!=b.size() || a[0].size()!=b[0].size())
    {
        // transpose if a and b are vector matrices
        if ((a.size()==b[0].size()) && (a[0].size()==b.size()))
            return equals(a,transpose(b));
        return false;
    }

    // check equality for each i,j value in a and b
    for (int i=0; i<(int) a.size(); i++)
        for (int j=0; j<(int) a[0].size(); j++)
            if ((a[i][j]<b[i][j]-SIGMA) || (a[i][j]>b[i][j]+SIGMA))
                return false;
    return true;
}

// returns the [nxn] identity matrix
matrix_t identity(int n)
{
    matrix_t result = zero(n);                  // identity = {{1, 0, 0},
    for (int i=0; i<n; i++)                     //             {0, 1, 0},
        result[i][i] = 1.0;                     //             {0, 0, 1}}
    return result;
}

// returns an [nxn] matrix of elements set to zero
matrix_t zero(int n)
{
    return zero(n,n);
}
// returns an [mxn] matrix of elements set to zero
matrix_t zero(int m, int n)
{
    matrix_t result(m, vector_t(n));
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            result[i][j] = 0.0;
    return result;
}

// returns a matrix set to the values of the elements in the matrix
matrix_t set(matrix_t matrix)
{
    matrix_t result(matrix.size(), vector_t(matrix[0].size()));
    for (int i=0; i<(int)matrix.size(); i++)
        for (int j=0; j<(int)matrix[0].size(); j++)
            result[i][j] = matrix[i][j];
    return result;
}
// returns a matrix set to the values of the elements in the matrix (vector)
matrix_t set(vector_t matrix)
{
    return set(transpose(matrix));
}
// returns a matrix (vector) set to the values, result = {x,y,0}
matrix_t set(double x, double y, double z)
{
    matrix_t result = zero(3,1);
    result[0][0] = x;
    result[1][0] = y;
    result[2][0] = z;
    return result;
}
// returns a matrix (vector) set to the values, result = {x,y,0}
matrix_t set(double x, double y)
{
    return set(x, y, 0.0);
}
// returns a matrix (vector) set to the values, result = {x,y,z}
matrix_t set(int x, int y, int z)
{
    return set(x/1.0, y/1.0, z/1.0);
}
// returns a matrix (vector) set to the values, result = {x,y,0}
matrix_t set(int x, int y)
{
    return set(x,y,0);
}

/* -- destructors -- */

// clears the data in the pixmap matrix
void destroy(pixmap_t pixmap)
{
    for (int i=0; i<(int) pixmap.size(); i++)
        destroy(pixmap[i]);
    pixmap.clear();
}
// clears the data in the matrix
void destroy(matrix_t matrix)
{
    for (int i=0; i<(int) matrix.size(); i++)
        destroy(matrix[i]);
    matrix.clear();
}
// clears the data in the vector matrix
void destroy(vector_t vector)
{
    vector.clear();
}

/* -- matrix statistics -- */

// returns the normalized value to the getMax and getMin values
double normalize(double max, double min, double value)
{
    return (value-min)/(max-min);
}
// returns the normalized value in the matrix
double normalize(double value, matrix_t matrix)
{
    return normalize(getMax(matrix),getMin(matrix),value);
}
// returns the normalized matrix
matrix_t normalize(matrix_t matrix)
{
    matrix_t result = zero((int) matrix.size(), (int) matrix[0].size());
    for (int i=0; i<(int) matrix.size(); i++)
        for (int j=0; j<(int) matrix[0].size(); j++)
            result[i][j] = normalize(matrix[i][j], matrix);
    return result;
}

// returns the maximum value in the pixmap
double getMax(pixmap_t pixmap)
{
    double result = getMax(pixmap[0]);
    for (int i=1; i<(int) pixmap.size(); i++)
        if (getMax(pixmap[i]) > result)
            result = getMax(pixmap[i]);
    return result;
}
// returns the minimum value in the pixmap
double getMin(pixmap_t pixmap)
{
    double result = getMin(pixmap[0]);
    for (int i=1; i<(int) pixmap.size(); i++)
        if (getMin(pixmap[i]) < result)
            result = getMin(pixmap[i]);
    return result;
}
// returns the maximum value in the matrix
double getMax(matrix_t matrix)
{
    double result = getMax(matrix[0]);
    for (int i=1; i<(int) matrix.size(); i++)
        if (getMax(matrix[i]) > result)
            result = getMax(matrix[i]);
    return result;
}
// returns the minimum value in the matrix
double getMin(matrix_t matrix)
{
    double result = getMin(matrix[0]);
    for (int i=1; i<(int) matrix.size(); i++)
        if (getMin(matrix[i]) < result)
            result = getMin(matrix[i]);
    return result;
}
// returns the maximum value in the matrix (vector)
double getMax(vector_t matrix)
{
    double result = matrix[0];
    for (int i=0; i<(int) matrix.size(); i++)
        if (matrix[i] > result)
            result = matrix[i];
    return result;
}
// returns the minimum value in the matrix (vector)
double getMin(vector_t matrix)
{
    double result = matrix[0];
    for (int i=0; i<(int) matrix.size(); i++)
        if (matrix[i] < result)
            result = matrix[i];
    return result;
}

// returns the maximum of the two values
double getMax(double a, double b)
{
    return (a>b) ? a : b;
}
// returns the minimum of the two values
double getMin(double a, double b)
{
    return (a<b) ? a : b;
}

/* -- printers -- */

// prints the matrix, with a specified name
void printMatrix(const char *name, matrix_t matrix)
{
    printf("%s = ", name);
    printMatrix(matrix);
}
// prints the vector matrix, with a specified name
void printMatrix(const char *name, vector_t matrix)
{
    printf("%s = ", name);
    printMatrix(matrix);
}
// prints the matrix
void printMatrix(matrix_t matrix)
{
    for (int i=0; i<(int) matrix.size(); i++)
        printMatrix(matrix[i]);
    printf("\n");
}
// prints the vector matrix
void printMatrix(vector_t matrix)
{
    for (int i=0; i<(int) matrix.size(); i++)
        printf("\t%7.2f", matrix[i]);
    printf("\n");
}

