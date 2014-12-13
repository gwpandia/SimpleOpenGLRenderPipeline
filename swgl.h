/* Copyright (c) 2009 by pandia. All Rights Reserved */

/***************************************************************************
   NAME
     swgl.cpp
   PURPOSE
     OpenGL pipeline simulation.
     BresenhamLine, Triangle plane filling, and Z-Buffer implementation.
     Phong shading
   NOTES
     I hate matrix.
     And I also hate precising problem...
   AUTHOR
     Chih-hung, Liu (pandia@nccu.edu.tw)
   HISTORY
     Chih-hung, Liu - Oct 24, 2009: Created.
     Chih-hung, Liu - Nov 24, 2009: HW2 Final Version, 
                                    with resize event handled
     Chih-hung, Liu - Dec 15, 2009: HW3 RC1 multi-model added.
     Chih-hung, Liu - Dec 15, 2009: Toggle interpolate Normal, 
                                    Light, View Vector.
***************************************************************************/

#ifndef __swgl_h__
#define __swgl_h__
#include "GLee.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//implement the following function call to archive the opengl pipeline

void swTranslated(GLdouble x, GLdouble y, GLdouble z);
void swScaled(GLdouble x, GLdouble y, GLdouble z);
void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z);

void swMatrixMode(GLenum mode);
void swLoadIdentity(void);
void swLoadMatrixd(const GLdouble * m);
void swMatrixMulplication(const GLdouble *first, const GLdouble *m, 
    GLdouble *result, int firstRow, int middle, int lastColumn);
void swMatrixMulplication(const GLdouble *first, const GLfloat *m, 
    GLfloat *result, int fR, int mid, int lC);
void swMultMatrixd(const GLdouble * m);
void swPushMatrix(void);
void swPopMatrix(void);

void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ);
void swFrustum(	GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble nearVal, GLdouble farVal);
void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

void swViewport(GLint x, GLint y, GLsizei width, GLsizei height);

bool swTransformation(const GLdouble h[4], GLdouble w[4]);
void copyMatrix(const GLdouble *source, GLdouble *dest, int size);
GLdouble innerProduct(GLdouble *a, int);

//---------------------------------------------------------------------------
//cghw2
//---------------------------------------------------------------------------
void writepixel(int x, int y, GLdouble r, GLdouble g, GLdouble b);

//Bresenham's algorithm
bool BresenhamLine(int x1, int y1, int x2, int y2, GLdouble r, GLdouble g, GLdouble b);
void ScanLine(int& x1, int& y1, int& x2, int& y2);
void doScanLine(int x1, int x2, int y, int ystep, int deltax, int deltay, int error, bool steep);
bool swTriangle(GLdouble x1, GLdouble y1, GLdouble z1,
			 GLdouble x2, GLdouble y2, GLdouble z2,
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble r, GLdouble g, GLdouble b);
			 
bool isInBound(int test, int min, int max);
bool swInitZbuffer(int width, int height);
bool swClearZbuffer();

//---------------------------------------------------------------------------
//cghw3
//---------------------------------------------------------------------------
void writepixelfast(int x, int y, GLdouble r, GLdouble g, GLdouble b);

bool swZbufferResize(int w, int h);
GLdouble swZbufferValue(int x, int y);

void swClear(GLbitfield mask);

bool swNormalTransformation(const GLdouble h[4], GLdouble w[4]);

void swLightfv(GLenum light, GLenum pname, const GLfloat *params);
void swMaterialfv (GLenum face, GLenum pname, const GLfloat *params);
void swMateriali (GLenum face, GLenum pname, GLint param);
void normalize(GLdouble *w);
GLdouble innerProduct(GLdouble *a, GLdouble *b, int size);
void gouraudShad(GLdouble *L, GLdouble *V, GLdouble *N, GLdouble *C);

void ScanLineG(int &x1, int &y1, int &x2, int &y2, GLdouble &z1, GLdouble &z2, 
    GLdouble &r1, GLdouble &g1, GLdouble &b1, GLdouble &r2, GLdouble &g2, GLdouble &b2,
    GLdouble &nx1, GLdouble &ny1, GLdouble &nz1, GLdouble &nw1,
    GLdouble &nx2, GLdouble &ny2, GLdouble &nz2, GLdouble &nw2,
    GLdouble &vx1, GLdouble &vy1, GLdouble &vz1, GLdouble &vw1,
    GLdouble &vx2, GLdouble &vy2, GLdouble &vz2, GLdouble &vw2,
    GLdouble &lx1, GLdouble &ly1, GLdouble &lz1, GLdouble &lw1,
    GLdouble &lx2, GLdouble &ly2, GLdouble &lz2, GLdouble &lw2);

void doScanLineG(int x1, int x2, int y, int ystep, int deltax, int deltay, int error, bool steep, GLdouble z1, GLdouble z2, 
    GLdouble r1, GLdouble g1, GLdouble b1, GLdouble r2, GLdouble g2, GLdouble b2,
    GLdouble nx1, GLdouble ny1, GLdouble nz1, GLdouble nw1,
    GLdouble nx2, GLdouble ny2, GLdouble nz2, GLdouble nw2,
    GLdouble vx1, GLdouble vy1, GLdouble vz1, GLdouble vw1,
    GLdouble vx2, GLdouble vy2, GLdouble vz2, GLdouble vw2,
    GLdouble lx1, GLdouble ly1, GLdouble lz1, GLdouble lw1,
    GLdouble lx2, GLdouble ly2, GLdouble lz2, GLdouble lw2);
//Gouraud shading
//vertex position:	(x1, y1, z1)   in object space
//vertex normal:	(nx1, ny1, nz1)
//vertex color:		(r1, g1, b1)
bool swTriangleG(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3);


////Phong Shading
bool swTriangleP(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3);

void printZBuffer();

void toggleInterNormal();
void toggleRealGouraud();

#endif                  /* __swgl_h__ */
