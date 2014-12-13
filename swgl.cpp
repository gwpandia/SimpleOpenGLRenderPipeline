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

#include "swgl.h"
#include <iostream>
#include <cmath>
#include <stack>
#include <climits>
#include <fstream>
#include "math3d.h"
#define PI 3.14159265358979323846
#define DOUBLE_MIN -99999999.0
#define DOUBLE_MAX 99999999.0

using namespace std;

GLdouble CTM_MV[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};	//Current Transformation Matrix: ModelView
GLdouble CTM_P[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};		//Current Transformation Matrix: Projection
GLenum CURRENT_MODE = GL_MODELVIEW;
GLint Xw = 0, Yw = 0, gWidth = 0, gHeight = 0;
GLdouble **zbuffer = NULL;
int Z_WIDTH = -1, Z_HEIGHT = -1;
int *ScannedX = NULL, *ScannedY = NULL;
GLdouble *ScannedR = NULL, *ScannedG = NULL, *ScannedB = NULL, *ScannedZ = NULL,
    *ScannedLX = NULL, *ScannedVX = NULL, *ScannedNX = NULL, 
    *ScannedLY = NULL, *ScannedVY = NULL, *ScannedNY = NULL, 
    *ScannedLZ = NULL, *ScannedVZ = NULL, *ScannedNZ = NULL, 
    *ScannedLW = NULL, *ScannedVW = NULL, *ScannedNW = NULL;
int ScannedLength = 0;
bool isZBUFFER = false;
GLdouble COP[4];

typedef struct matrix{
    GLdouble matrix[16];
} Matrix;

stack<Matrix> STACK_MV;
stack<Matrix> STACK_P;

bool debug = false;
bool isInterNormal = false;
bool isRealGouraud = false;

bool swTransformation(const GLdouble h[4], GLdouble w[4])
{
    
	//p = CTM_P*CTM_MV*h
    GLdouble temp[16];
    GLdouble p[4];
    swMatrixMulplication(CTM_P, CTM_MV, temp, 4, 4, 4);
    swMatrixMulplication(temp, h, p, 4, 4, 1);
    
	//prespective division
    //cout << "p3:" << p[3] << endl;
    w[2] = p[2];
    if(p[3] != 0){
        for(int i = 0; i < 4; i++){
		  p[i] = p[i] / p[3];
        }
    }

	//viewport transformation
    
    w[0] = (p[0]+1)*(gWidth/2)+Xw;
	w[1] = (p[1]+1)*(gHeight/2)+Yw;
	
	w[3] = p[3];
    //w[2] = (p[2]+1.0)/2.0;
    //std::cout << "w2" << w[2] << std::endl;
	return true;
}

void swTranslated(GLdouble x, GLdouble y, GLdouble z){
    GLdouble translate[16] = {1,0,0,0,0,1,0,0,0,0,1,0,x,y,z,1};
    swMultMatrixd(translate);
}

void swScaled(GLdouble x, GLdouble y, GLdouble z){
    GLdouble scale[16] = {x,0,0,0,0,y,0,0,0,0,z,0,0,0,0,1};
    swMultMatrixd(scale);
}

void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z){
    GLdouble radian = angle * PI / 180.0;
    GLdouble dist = sqrt(x*x + y*y + z*z);
    GLdouble C = cos(radian);
    GLdouble S = sin(radian);
    
    // Normalizing
    if(dist == 0){
        std::cerr << "Rotated Error" << std::endl;
        return;    
    }
    x /= dist;
    y /= dist;
    z /= dist;

    GLdouble rotate[16] = {
        x*x*(1-C)+C ,y*x*(1-C)+z*S, x*z*(1-C)-y*S, 0,
        x*y*(1-C)-z*S, y*y*(1-C)+C, y*z*(1-C)+x*S, 0,
        x*z*(1-C)+y*S, y*z*(1-C)-x*S, z*z*(1-C)+C, 0,
        0, 0, 0, 1};
    
    swMultMatrixd(rotate);
}

void swMatrixMode(GLenum mode){
    CURRENT_MODE = mode;
}

void swLoadIdentity(void){
    GLdouble identity[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    swLoadMatrixd(identity);
}

void swLoadMatrixd(const GLdouble * m){
    if(CURRENT_MODE == GL_PROJECTION){
        copyMatrix(m, CTM_P, 16);
    }
    else if(CURRENT_MODE == GL_MODELVIEW){
        copyMatrix(m, CTM_MV, 16);
    }
    else{}
}

void copyMatrix(const GLdouble *source, GLdouble *dest, int size){
    for(int i = 0; i < size; i++){
        dest[i] = source[i];
    }
}

GLdouble innerProduct(GLdouble *a, int mid){
    GLdouble result = 0;
    for(int i = 0; i < mid; i++){
        result += (a[i] * a[i+mid]);
    }    
    return result;
}

void swMatrixMulplication(const GLdouble *first, const GLdouble *m, GLdouble *result, int fR, int mid, int lC){
    
    if((fR==0) || (mid==0) || (lC==0)){
        std::cerr << "Divided by zero" << std::endl;
        return;    
    }
    
    for(int i = 0; i < fR*lC; i++){
        GLdouble temp[2*mid];
        for(int j = 0; j < mid; j++){
            temp[j] = first[(i%fR) + fR*j];
            temp[j+mid] = m[(i/fR)*mid + j];
        }
        
        result[i] = innerProduct(temp, mid);
    }
    
}

void swMatrixMulplication(const GLdouble *first, const GLfloat *m, GLfloat *result, int fR, int mid, int lC){
    
    if((fR==0) || (mid==0) || (lC==0)){
        std::cerr << "Divided by zero" << std::endl;
        return;    
    }
    
    for(int i = 0; i < fR*lC; i++){
        GLdouble temp[2*mid];
        for(int j = 0; j < mid; j++){
            temp[j] = first[(i%fR) + fR*j];
            temp[j+mid] = (double)m[(i/fR)*mid + j];
        }
        
        result[i] = (float)innerProduct(temp, mid);
    }
    
}

void swMultMatrixd(const GLdouble * m){
    GLdouble* first;
    GLdouble result[16];
    if(CURRENT_MODE == GL_PROJECTION){
        first = CTM_P;
    }
    else if(CURRENT_MODE == GL_MODELVIEW){
        first = CTM_MV;
    }
    else{}
    
    swMatrixMulplication(first, m, result, 4, 4, 4);
    copyMatrix(result, first, 16);
}

void swPushMatrix(void){
    if(CURRENT_MODE == GL_PROJECTION){
        Matrix temp;
        copyMatrix(CTM_P, temp.matrix, 16);
        STACK_P.push(temp);
    }
    else if(CURRENT_MODE == GL_MODELVIEW){
        Matrix temp;
        copyMatrix(CTM_MV, temp.matrix, 16);
        STACK_MV.push(temp);
    }
    else{}
}

void swPopMatrix(void){
    if(CURRENT_MODE == GL_PROJECTION && !STACK_P.empty()){
        copyMatrix(STACK_P.top().matrix, CTM_P, 16);
        STACK_P.pop();
    }
    else if(CURRENT_MODE == GL_MODELVIEW && !STACK_MV.empty()){
        copyMatrix(STACK_MV.top().matrix, CTM_MV, 16);
        STACK_MV.pop();
    }
    else{}
}

void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ){
    
    GLdouble f[3] = {centerX-eyeX, centerY-eyeY, centerZ-eyeZ};
    GLdouble up[3] = {upX, upY, upZ};
    GLdouble dist_f = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    GLdouble dist_up = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    COP[0] = eyeX;
    COP[1] = eyeY;
    COP[2] = eyeZ;
    COP[3] = 0;
    if((dist_f==0) || (dist_up==0)){
        std::cerr << "Divided by zero" << std::endl;    
        return;
    }
    
    // Normalization
    for(int i = 0; i < 3; i++){
        f[i] /= dist_f;
        up[i] /= dist_up;    
    }
    
    GLdouble s[] = { f[1]*up[2] - f[2]*up[1], f[2]*up[0] - f[0]*up[2] , f[0]*up[1]-f[1]*up[0] };
    GLdouble u[] = { s[1]*f[2] - s[2]*f[1], s[2]*f[0] - s[0]*f[2] , s[0]*f[1]-s[1]*f[0] };
    
    GLdouble M[16] = {s[0], u[0], -f[0], 0, 
                      s[1], u[1], -f[1], 0, 
                      s[2], u[2], -f[2], 0, 
                      0, 0, 0, 1};
    
    swMultMatrixd(M);
    swTranslated(-eyeX, -eyeY, -eyeZ);
    
}

void swFrustum(	GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble nearVal, GLdouble farVal){
    
    if((right==left) || (top==bottom) || (farVal==nearVal)){
        std::cerr << "Divided by zero" << std::endl;    
        return;    
    }
    GLdouble A = (right+left)/(right-left);
    GLdouble B = (top+bottom)/(top-bottom);
    GLdouble C = -(farVal+nearVal)/(farVal-nearVal);
    GLdouble D = -(2*farVal*nearVal)/(farVal-nearVal);
    GLdouble temp1 = (2*nearVal)/(right-left);
    GLdouble temp2 = (2*nearVal)/(top-bottom);
    GLdouble M[16] = {temp1, 0, 0, 0, 0, temp2, 0, 0, A, B, C, -1, 0, 0, D, 0};
    swMultMatrixd(M);
}

void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar){
    
    if(zFar==zNear){
        std::cerr << "Divided by zero" << std::endl;    
        return;    
    }
    GLdouble radian = fovy * PI / 180.0;
    GLdouble f = 1/tan(radian/2.0);
    GLdouble temp1 = (zFar+zNear)/(zNear-zFar);
    GLdouble temp2 = (2*zFar*zNear)/(zNear-zFar);
    GLdouble M[16] = {f/aspect, 0, 0, 0,
                        0, f, 0, 0,
                        0, 0, temp1, -1, 
                        0, 0, temp2, 0};
                        
    swMultMatrixd(M);
}

void swViewport(GLint x, GLint y, GLsizei width, GLsizei height){    
    //GLint coor[4];
    //glGetIntegerv(GL_VIEWPORT, coor);
    
    //Xw = (coor[0]+1)*(width/2) + x;
    //Yw = (coor[1]+1)*(height/2) + y;
    //swClearZbuffer();
    //swInitZbuffer(width*2, height);
    Xw = x;
    Yw = y;
	gWidth = width;
	gHeight = height;
}

//---------------------------------------------------------------------------
//cghw2
//---------------------------------------------------------------------------

void writepixel(int x, int y, GLdouble r, GLdouble g, GLdouble b)
{
	GLubyte map[1]={255};

	glColor3d(r, g, b);
	glWindowPos2i(x, y);
	glBitmap(1, 1, 0, 0, 0, 0, map);
}

bool BresenhamLine(int x1, int y1, int x2, int y2, GLdouble r, GLdouble g, GLdouble b)
{
    ScanLine(x1, y1, x2, y2);
    for(int i = 0; i < x2-x1+1; i++){
        writepixel(ScannedX[i], ScannedY[i], r, g, b);
        //std::cout << ScannedX[i] << " " << ScannedY[i] << " " << r << g << b << std::endl;
    }
	return true;
}

void ScanLine(int &x1, int &y1, int &x2, int &y2){
    bool steep = abs(y2-y1) > abs(x2-x1);
    if(steep){
        swap(x1, y1);
        swap(x2, y2);
    }
    if(x1 > x2){
        swap(x1, x2);
        swap(y1, y2);
    }
    int deltax = x2 - x1;
    int deltay = abs(y2 - y1);
    int error = deltax / 2;
    int ystep = (y1 < y2)? 1 : -1;
    int y = y1;
    doScanLine(x1, x2, y, ystep, deltax, deltay, error, steep);
}

void doScanLine(int x1, int x2, int y, int ystep, int deltax, int deltay, int error, bool steep){
    delete [] ScannedX;
    delete [] ScannedY;
    ScannedLength = x2 - x1 + 1;
    ScannedX = new int[ScannedLength];
    ScannedY = new int[ScannedLength];
    
    for(int x = x1; x <= x2; x++){
        if(steep){
            ScannedX[x-x1] = y;
            ScannedY[x-x1] = x;
        }
        else{
            ScannedX[x-x1] = x;
            ScannedY[x-x1] = y;
        }
        error -= deltay;
        if(error < 0){
            y += ystep;
            error += deltax;
        }
    }
}

bool swTriangle(GLdouble x1, GLdouble y1, GLdouble z1,
			 GLdouble x2, GLdouble y2, GLdouble z2,
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble r, GLdouble g, GLdouble b)
{
    GLdouble points[][3] = {{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}};
    int Max_x = round(max(x1, max(x2, x3)));
    int Min_x = round(min(x1, min(x2, x3)));
    int YBoundLength = Max_x-Min_x+1;
    int *Max_YBound = new int[YBoundLength];
    int *Min_YBound = new int[YBoundLength];
    GLdouble *Max_ZBound = new GLdouble[YBoundLength];
    GLdouble *Min_ZBound = new GLdouble[YBoundLength];
    int _x1, _x2, _y1, _y2;
    GLdouble _z1, _z2;
    
    for(int i = 0; i < YBoundLength; i++){
        Max_YBound[i] = INT_MIN;    
        Min_YBound[i] = INT_MAX;
        Max_ZBound[i] = DOUBLE_MAX; //MIN
        Min_ZBound[i] = DOUBLE_MAX; //MIN
    }
    

    for(int i = 0; i < 3; i++){
        _x1 = round(points[i%3][0]);
        _x2 = round(points[(i+1)%3][0]);
        _y1 = round(points[i%3][1]);
        _y2 = round(points[(i+1)%3][1]);
        _z1 = points[i%3][2];
        _z2 = points[(i+1)%3][2];
        ScanLine(_x1, _y1, _x2, _y2);
        
        GLdouble zstep = (_z2-_z1) / ScannedLength;
        
        for(int j = 0; j < ScannedLength; j++){
            GLdouble cur_z = _z1+ zstep*j;
            if(ScannedY[j] > Max_YBound[ScannedX[j]-Min_x]){
                Max_YBound[ScannedX[j]-Min_x] = ScannedY[j];
                Max_ZBound[ScannedX[j]-Min_x] = cur_z;
            }
            if(ScannedY[j] < Min_YBound[ScannedX[j]-Min_x]){
                Min_YBound[ScannedX[j]-Min_x] = ScannedY[j];
                Min_ZBound[ScannedX[j]-Min_x] = cur_z;
            }
        }
    }
    
    for(int i = 0; i < YBoundLength; i++){
        int cur_x = Min_x + i;
        GLdouble zstep = (Max_ZBound[i] - Min_ZBound[i]) / (Max_YBound[i] - Min_YBound[i] + 1);
        for(int j = Min_YBound[i]; j <= Max_YBound[i]; j++){
            GLdouble cur_z = Min_ZBound[i] + zstep * (j-Min_YBound[i]);
            if(isZBUFFER && isInBound(j, 0, Z_HEIGHT) && isInBound(cur_x, 0, Z_WIDTH) && (cur_z <= zbuffer[j][cur_x])){ // >=
                writepixel(cur_x, j, r, g, b);
                zbuffer[j][cur_x] = cur_z;
            }
            else if(!isZBUFFER){
                writepixel(cur_x, j, r, g, b);                
            }
        }
    }
    
    delete [] Max_YBound;
    delete [] Min_YBound;
    delete [] Max_ZBound;
    delete [] Min_ZBound;
    
	return true;
}

bool isInBound(int test, int min, int max){
    return (test >= min) && (test < max);
}

bool swInitZbuffer(int width, int height){
    if(!isZBUFFER){
        Z_WIDTH = width;
        Z_HEIGHT = height;
        zbuffer = new GLdouble*[height];
        for(int i = 0; i < Z_HEIGHT; i++){
            zbuffer[i] = new GLdouble[Z_WIDTH];
            for(int j = 0; j < Z_WIDTH; j++){
                zbuffer[i][j] = DOUBLE_MAX; //MIN
            }
        }
        
        isZBUFFER = true;
    }
    return isZBUFFER;
}

bool swClearZbuffer(){
    if(isZBUFFER){
        for(int i = 0; i < Z_HEIGHT; i++ ){
            delete [] zbuffer[i];
        }
        delete [] zbuffer;
        zbuffer = NULL;
        Z_WIDTH = Z_HEIGHT = -1;
        isZBUFFER = false;
    }
    return !isZBUFFER;
}

//---------------------------------------------------------------------------
//cghw2 - software z buffer
//---------------------------------------------------------------------------
//GLdouble *ZBUFFER;

bool swZbufferResize(int w, int h)
{
	if( w!=Z_WIDTH || h!=Z_HEIGHT ) {
		swClearZbuffer();
		swInitZbuffer(w, h);
	}
	return true;
}

GLdouble swZbufferValue(int x, int y)
{
	if(x>=0 && x<Z_WIDTH && y>=0 && y<Z_HEIGHT)
		return zbuffer[y][x];
	else
		return 1.0;
}

void swClear(GLbitfield mask)
{
	for(int i = 0; i<Z_HEIGHT; ++i){
        for(int j = 0; j < Z_WIDTH; ++j){
            zbuffer[i][j] = DOUBLE_MAX;
        }    
    }
		
}


//---------------------------------------------------------------------------
//cghw3
//---------------------------------------------------------------------------


GLfloat  _ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
GLfloat  _diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };
GLfloat  _specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat  _specref[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat  _shininess = 128.0f;

GLfloat _lightPos[] = { -50.f, 50.0f, 100.0f, 1.0f };
GLfloat cur_lightPos[] = {0.f, 0.f, 0.f, 0.f};

void writepixelfast(int x, int y, GLdouble r, GLdouble g, GLdouble b)
{
	glBegin(GL_POINTS);
		glColor3d(r, g, b);
		glVertex2i(x, y);
	glEnd();
}



bool swNormalTransformation(const GLdouble h[4], GLdouble w[4])
{
	//apply transformation
	GLdouble mv[16], tmp[16];
	m3dInvertMatrix44(mv, CTM_MV);
	//m3dTransposeMatrix44(mv, tmp);
	swMatrixMulplication(h, mv, w, 1, 4, 4);

	//Normalize
    normalize(w);
	return true;
}


void swLightfv(GLenum light, GLenum pname, const GLfloat *params)
{
	switch(pname) {
		case GL_AMBIENT:
			_ambientLight[0]=params[0];
			_ambientLight[1]=params[1];
			_ambientLight[2]=params[2];
			_ambientLight[3]=params[3];
			break;

		case GL_DIFFUSE:
			_diffuseLight[0]=params[0];
			_diffuseLight[1]=params[1];
			_diffuseLight[2]=params[2];
			_diffuseLight[3]=params[3];
			break;

		case GL_SPECULAR:
			_specular[0]=params[0];
			_specular[1]=params[1];
			_specular[2]=params[2];
			_specular[3]=params[3];
			break;

		case GL_POSITION:
			_lightPos[0]=params[0];
			_lightPos[1]=params[1];
			_lightPos[2]=params[2];
			_lightPos[3]=params[3];

			//------------------------------------------------------
			//ADD transformation to eye coordinate
            swMatrixMulplication(CTM_MV, _lightPos, cur_lightPos, 4, 4, 1);

			break;
	}
}

void swMaterialfv (GLenum face, GLenum pname, const GLfloat *params)
{
	switch(pname) {
		case GL_SPECULAR:
			_specref[0]=params[0];
			_specref[1]=params[1];
			_specref[2]=params[2];
			_specref[3]=params[3];
			break;
	}
}

void swMateriali (GLenum face, GLenum pname, GLint param)
{
	switch(pname) {
		case GL_SHININESS:
			_shininess=param;
			break;
	}
}

void normalize(GLdouble *w){
    GLdouble length = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    for(int i = 0; i < 3; i++){
        w[i] /= length;
    }
}

GLdouble innerProduct(GLdouble *a, GLdouble *b, int size){
    GLdouble ret = 0;
    for(int i = 0; i < size; i++){
        ret += a[i] * b[i];    
    }
    return ret;
}

void gouraudShad(GLdouble *L, GLdouble *V, GLdouble *N, GLdouble *C){
    normalize(L);
    normalize(V);
    GLdouble H[4];
    for(int i = 0; i < 3; i++){
        H[i] = L[i] + V[i];    
    }
    H[3] = 0;
    normalize(H);
    GLdouble nh = innerProduct(N, H, 4);
	if(nh < 0){
		nh = 0;
	}
    GLdouble diff = innerProduct(L, N, 4);
    GLdouble spec = pow(nh, (GLdouble)_shininess);
    if(diff < 0){
        diff = 0;
    }
    
    if(spec < 0){
        spec = 0;
    }

    for(int i = 0; i < 3; i++){
        C[i] = ( C[i] * ( _diffuseLight[i] * diff + _ambientLight[i]) + _specular[i] * _specref[i] * spec);   
        //C[i] *= _diffuseLight[i] * diff + _specular[i] * _specref[i] * spec + _ambientLight[i];
    }
}

void ScanLineG(int &x1, int &y1, int &x2, int &y2, GLdouble &z1, GLdouble &z2, 
	GLdouble &r1, GLdouble &g1, GLdouble &b1, GLdouble &r2, GLdouble &g2, GLdouble &b2, 
    GLdouble &nx1, GLdouble &ny1, GLdouble &nz1, GLdouble &nw1,
    GLdouble &nx2, GLdouble &ny2, GLdouble &nz2, GLdouble &nw2,
    GLdouble &vx1, GLdouble &vy1, GLdouble &vz1, GLdouble &vw1,
    GLdouble &vx2, GLdouble &vy2, GLdouble &vz2, GLdouble &vw2,
    GLdouble &lx1, GLdouble &ly1, GLdouble &lz1, GLdouble &lw1,
    GLdouble &lx2, GLdouble &ly2, GLdouble &lz2, GLdouble &lw2){
    bool steep = abs(y2-y1) > abs(x2-x1);
    if(steep){
        swap(x1, y1);
        swap(x2, y2);
    }
    if(x1 > x2){
        swap(x1, x2);
        swap(y1, y2);
        swap(z1, z2);
        swap(r1, r2);
        swap(g1, g2);
        swap(b1, b2);
        if(isInterNormal){
            swap(nx1, nx2);
            swap(ny1, ny2);
            swap(nz1, nz2);
            swap(nw1, nw2);
            swap(vx1, vx2);
            swap(vy1, vy2);
            swap(vz1, vz2);
            swap(vw1, vw2);
            swap(lx1, lx2);
            swap(ly1, ly2);
            swap(lz1, lz2);
            swap(lw1, lw2);
        }
    }
    int deltax = x2 - x1;
    int deltay = abs(y2 - y1);
    int error = deltax / 2;
    int ystep = (y1 < y2)? 1 : -1;
    int y = y1;
    doScanLineG(x1, x2, y, ystep, deltax, deltay, error, steep, z1, z2, r1, g1, b1, r2, g2, b2,
        nx1, ny1, nz1, nw1, nx2, ny2, nz2, nw2, vx1, vy1, vz1, vw1, vx2, vy2, vz2, vw2,
        lx1, ly1, lz1, lw1, lx2, ly2, lz2, lw2);
}

void doScanLineG(int x1, int x2, int y, int ystep, int deltax, int deltay, int error, bool steep, 
	GLdouble z1, GLdouble z2, GLdouble r1, GLdouble g1, GLdouble b1, GLdouble r2, GLdouble g2, GLdouble b2,
    GLdouble nx1, GLdouble ny1, GLdouble nz1, GLdouble nw1,
    GLdouble nx2, GLdouble ny2, GLdouble nz2, GLdouble nw2,
    GLdouble vx1, GLdouble vy1, GLdouble vz1, GLdouble vw1,
    GLdouble vx2, GLdouble vy2, GLdouble vz2, GLdouble vw2,
    GLdouble lx1, GLdouble ly1, GLdouble lz1, GLdouble lw1,
    GLdouble lx2, GLdouble ly2, GLdouble lz2, GLdouble lw2){
    delete [] ScannedX;
    delete [] ScannedY;
    delete [] ScannedZ;
	delete [] ScannedR;
    delete [] ScannedG;
    delete [] ScannedB;
    
    if(isInterNormal){
        delete [] ScannedNX;
        delete [] ScannedNY;
        delete [] ScannedNZ;
        delete [] ScannedNW;
        delete [] ScannedVX;
        delete [] ScannedVY;
        delete [] ScannedVZ;
        delete [] ScannedVW;
        delete [] ScannedLX;
        delete [] ScannedLY;
        delete [] ScannedLZ;
        delete [] ScannedLW;
    }
    
    ScannedLength = x2 - x1 + 1;
    ScannedX = new int[ScannedLength];
    ScannedY = new int[ScannedLength];
    ScannedZ = new double[ScannedLength];
	ScannedR = new double[ScannedLength];
    ScannedG = new double[ScannedLength];
    ScannedB = new double[ScannedLength];
    
    if(isInterNormal){
        ScannedNX = new double[ScannedLength];
        ScannedNY = new double[ScannedLength];
        ScannedNZ = new double[ScannedLength];
        ScannedNW = new double[ScannedLength];
        ScannedVX = new double[ScannedLength];
        ScannedVY = new double[ScannedLength];
        ScannedVZ = new double[ScannedLength];
        ScannedVW = new double[ScannedLength];
        ScannedLX = new double[ScannedLength];
        ScannedLY = new double[ScannedLength];
        ScannedLZ = new double[ScannedLength];
        ScannedLW = new double[ScannedLength];
    }
    
    GLdouble rstep = (r2-r1) / ScannedLength;
    GLdouble gstep = (g2-g1) / ScannedLength;
    GLdouble bstep = (b2-b1) / ScannedLength;
    GLdouble zstep = (z2-z1) / ScannedLength;
    GLdouble nxstep = (nx2-nx1) / ScannedLength;
    GLdouble nystep = (ny2-ny1) / ScannedLength;
    GLdouble nzstep = (nz2-nz1) / ScannedLength;
    GLdouble nwstep = (nw2-nw1) / ScannedLength;
    GLdouble vxstep = (vx2-vx1) / ScannedLength;
    GLdouble vystep = (vy2-vy1) / ScannedLength;
    GLdouble vzstep = (vz2-vz1) / ScannedLength;
    GLdouble vwstep = (vw2-vw1) / ScannedLength;
    GLdouble lxstep = (lx2-lx1) / ScannedLength;
    GLdouble lystep = (ly2-ly1) / ScannedLength;
    GLdouble lzstep = (lz2-lz1) / ScannedLength;
    GLdouble lwstep = (lw2-lw1) / ScannedLength;
    //cout << "ZZZ: " << zstep << endl;
    //cout << z1 << " " << z2 << endl;
    for(int x = x1; x <= x2; x++){
        if(steep){
            ScannedX[x-x1] = y;
            ScannedY[x-x1] = x;
        }
        else{
            ScannedX[x-x1] = x;
            ScannedY[x-x1] = y;
        }
        ScannedZ[x-x1] = z1 + zstep * (x-x1);
		ScannedR[x-x1] = r1 + rstep * (x-x1);
        ScannedG[x-x1] = g1 + gstep * (x-x1);
        ScannedB[x-x1] = b1 + bstep * (x-x1);
        
        if(isInterNormal){
            ScannedNX[x-x1] = nx1 + nxstep * (x-x1);
            ScannedNY[x-x1] = ny1 + nystep * (x-x1);
            ScannedNZ[x-x1] = nz1 + nzstep * (x-x1);
            ScannedNW[x-x1] = nw1 + nwstep * (x-x1);
            ScannedVX[x-x1] = vx1 + vxstep * (x-x1);
            ScannedVY[x-x1] = vy1 + vystep * (x-x1);
            ScannedVZ[x-x1] = vz1 + vzstep * (x-x1);
            ScannedVW[x-x1] = vw1 + vwstep * (x-x1);
            ScannedLX[x-x1] = lx1 + lxstep * (x-x1);
            ScannedLY[x-x1] = ly1 + lystep * (x-x1);
            ScannedLZ[x-x1] = lz1 + lzstep * (x-x1);
            ScannedLW[x-x1] = lw1 + lwstep * (x-x1);
        }
        error -= deltay;
        if(error < 0){
            y += ystep;
            error += deltax;
        }
    }
}


//Gouraud shading
bool swTriangleG(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3)
{
    GLdouble points[][4] = {{x1, y1, z1, 1.0}, {x2, y2, z2, 1.0}, {x3, y3, z3, 1.0}};
    GLdouble n1[4] = {nx1, ny1, nz1, 0.0};
    GLdouble n2[4] = {nx2, ny2, nz2, 0.0};
    GLdouble n3[4] = {nx3, ny3, nz3, 0.0};
    GLdouble tranCOP[4] = {0.0, 0.0, 0.0, 0.0};
    GLdouble colors[][3] = {{r1, g1, b1}, {r2, g2, b2}, {r3, g3, b3}};
    GLdouble trP[][4] = {{x1, y1, z1, 1.0}, {x2, y2, z2, 1.0}, {x3, y3, z3, 1.0}};
    GLdouble light[3][4];
    GLdouble view[3][4];
    GLdouble trN[3][4];
    
	//transformation all data(vertex, normal, light vector, COP) to eye coordiante
    swMatrixMulplication(CTM_MV, points[0], trP[0], 4, 4, 1);
    swMatrixMulplication(CTM_MV, points[1], trP[1], 4, 4, 1);
    swMatrixMulplication(CTM_MV, points[2], trP[2], 4, 4, 1);
    swMatrixMulplication(CTM_MV, COP, tranCOP, 4, 4, 1);
    swNormalTransformation(n1, trN[0]);
    swNormalTransformation(n2, trN[1]);
    swNormalTransformation(n3, trN[2]);
    
    //Light and V
    for(int i = 0; i < 4; i++){
        light[0][i] = cur_lightPos[i] - trP[0][i];
        light[1][i] = cur_lightPos[i] - trP[1][i];
        light[2][i] = cur_lightPos[i] - trP[2][i];
        view[0][i] = tranCOP[i] - trP[0][i];
        view[1][i] = tranCOP[i] - trP[1][i];
        view[2][i] = tranCOP[i] - trP[2][i];
    }
    
	//modified Pong shading equation in each vertex
    GLdouble trGN[4];
    for(int i = 0; i < 4; i++){
        trGN[i] = (trN[0][i] + trN[1][i] + trN[2][i]) / 3.0;
    }

    if(!isRealGouraud){
        gouraudShad(light[0], view[0], trN[0], colors[0]);
        gouraudShad(light[1], view[1], trN[1], colors[1]);
        gouraudShad(light[2], view[2], trN[2], colors[2]);
    }
    else{
        gouraudShad(light[0], view[0], trGN, colors[0]);
        gouraudShad(light[1], view[1], trGN, colors[1]);
        gouraudShad(light[2], view[2], trGN, colors[2]);
    }
    swTransformation(points[0], trP[0]);
    swTransformation(points[1], trP[1]);
	swTransformation(points[2], trP[2]);

	//Raterization: 
    int Max_x = round(max(trP[0][0], max(trP[1][0], trP[2][0])));
    int Min_x = round(min(trP[0][0], min(trP[1][0], trP[2][0])));
    int YBoundLength = Max_x-Min_x+1;
    int *Max_YBound = new int[YBoundLength];
    int *Min_YBound = new int[YBoundLength];
    GLdouble *Max_ZBound = new GLdouble[YBoundLength];
    GLdouble *Min_ZBound = new GLdouble[YBoundLength];
    GLdouble (*Max_CBound)[3] = new GLdouble[YBoundLength][3];
    GLdouble (*Min_CBound)[3] = new GLdouble[YBoundLength][3];
    GLdouble (*Max_VBound)[4] = new GLdouble[YBoundLength][4];
    GLdouble (*Min_VBound)[4] = new GLdouble[YBoundLength][4];
    GLdouble (*Max_LBound)[4] = new GLdouble[YBoundLength][4];
    GLdouble (*Min_LBound)[4] = new GLdouble[YBoundLength][4];
    GLdouble (*Max_NBound)[4] = new GLdouble[YBoundLength][4];
    GLdouble (*Min_NBound)[4] = new GLdouble[YBoundLength][4];
    
    int _x1, _x2, _y1, _y2;
    GLdouble _z1, _z2, _r1, _r2, _g1, _g2, _b1, _b2;
    GLdouble _nx1, _ny1, _nz1, _nw1, _nx2, _ny2, _nz2, _nw2;
    GLdouble _vx1, _vy1, _vz1, _vw1, _vx2, _vy2, _vz2, _vw2;
    GLdouble _lx1, _ly1, _lz1, _lw1, _lx2, _ly2, _lz2, _lw2;
    
    for(int i = 0; i < YBoundLength; i++){
        Max_YBound[i] = INT_MIN;    
        Min_YBound[i] = INT_MAX;
        Max_ZBound[i] = DOUBLE_MAX; //MIN
        Min_ZBound[i] = DOUBLE_MAX; //MIN
    }
    
    for(int i = 0; i < 3; i++){
        _x1 = round(trP[i%3][0]);
        _x2 = round(trP[(i+1)%3][0]);
        _y1 = round(trP[i%3][1]);
        _y2 = round(trP[(i+1)%3][1]);
        _z1 = trP[i%3][2];
        _z2 = trP[(i+1)%3][2];
        _r1 = colors[i%3][0];
        _r2 = colors[(i+1)%3][0];
        _g1 = colors[i%3][1];
        _g2 = colors[(i+1)%3][1];
        _b1 = colors[i%3][2];
        _b2 = colors[(i+1)%3][2];
        if(isInterNormal){
            _nx1 = trN[i][0];
            _ny1 = trN[i][1];
            _nz1 = trN[i][2];
            _nw1 = trN[i][3];
            _nx2 = trN[(i+1)%3][0];
            _ny2 = trN[(i+1)%3][1];
            _nz2 = trN[(i+1)%3][2];
            _nw2 = trN[(i+1)%3][3];
            _vx1 = view[i][0];
            _vy1 = view[i][1];
            _vz1 = view[i][2];
            _vw1 = view[i][3];
            _vx2 = view[(i+1)%3][0];
            _vy2 = view[(i+1)%3][1];
            _vz2 = view[(i+1)%3][2];
            _vw2 = view[(i+1)%3][3];
            _lx1 = light[i][0];
            _ly1 = light[i][1];
            _lz1 = light[i][2];
            _lw1 = light[i][3];
            _lx2 = light[(i+1)%3][0];
            _ly2 = light[(i+1)%3][1];
            _lz2 = light[(i+1)%3][2];
            _lw2 = light[(i+1)%3][3];
        }
        
        ScanLineG(_x1, _y1, _x2, _y2, _z1, _z2, _r1, _g1, _b1, _r2, _g2, _b2,
                _nx1, _ny1, _nz1, _nw1, _nx2, _ny2, _nz2, _nw2,
                _vx1, _vy1, _vz1, _vw1, _vx2, _vy2, _vz2, _vw2,
                _lx1, _ly1, _lz1, _lw1, _lx2, _ly2, _lz2, _lw2);

        for(int j = 0; j < ScannedLength; j++){
            if(ScannedY[j] > Max_YBound[ScannedX[j]-Min_x]){
                Max_YBound[ScannedX[j]-Min_x] = ScannedY[j];
                Max_ZBound[ScannedX[j]-Min_x] = ScannedZ[j];
                Max_CBound[ScannedX[j]-Min_x][0] = ScannedR[j];//cur_r;
                Max_CBound[ScannedX[j]-Min_x][1] = ScannedG[j];//cur_g;
                Max_CBound[ScannedX[j]-Min_x][2] = ScannedB[j];//cur_b;
                if(isInterNormal){
                    Max_NBound[ScannedX[j]-Min_x][0] = ScannedNX[j];
                    Max_NBound[ScannedX[j]-Min_x][1] = ScannedNY[j];
                    Max_NBound[ScannedX[j]-Min_x][2] = ScannedNZ[j];
                    Max_NBound[ScannedX[j]-Min_x][3] = ScannedNW[j];
                    Max_LBound[ScannedX[j]-Min_x][0] = ScannedLX[j];
                    Max_LBound[ScannedX[j]-Min_x][1] = ScannedLY[j];
                    Max_LBound[ScannedX[j]-Min_x][2] = ScannedLZ[j];
                    Max_LBound[ScannedX[j]-Min_x][3] = ScannedLW[j];
                    Max_VBound[ScannedX[j]-Min_x][0] = ScannedVX[j];
                    Max_VBound[ScannedX[j]-Min_x][1] = ScannedVY[j];
                    Max_VBound[ScannedX[j]-Min_x][2] = ScannedVZ[j];
                    Max_VBound[ScannedX[j]-Min_x][3] = ScannedVW[j];
                }
            }
            if(ScannedY[j] < Min_YBound[ScannedX[j]-Min_x]){
                Min_YBound[ScannedX[j]-Min_x] = ScannedY[j];
                Min_ZBound[ScannedX[j]-Min_x] = ScannedZ[j];
                Min_CBound[ScannedX[j]-Min_x][0] = ScannedR[j];//cur_r;
                Min_CBound[ScannedX[j]-Min_x][1] = ScannedG[j];//cur_g;
                Min_CBound[ScannedX[j]-Min_x][2] = ScannedB[j];//cur_b;
                if(isInterNormal){
                    Min_NBound[ScannedX[j]-Min_x][0] = ScannedNX[j];
                    Min_NBound[ScannedX[j]-Min_x][1] = ScannedNY[j];
                    Min_NBound[ScannedX[j]-Min_x][2] = ScannedNZ[j];
                    Min_NBound[ScannedX[j]-Min_x][3] = ScannedNW[j];
                    Min_LBound[ScannedX[j]-Min_x][0] = ScannedLX[j];
                    Min_LBound[ScannedX[j]-Min_x][1] = ScannedLY[j];
                    Min_LBound[ScannedX[j]-Min_x][2] = ScannedLZ[j];
                    Min_LBound[ScannedX[j]-Min_x][3] = ScannedLW[j];
                    Min_VBound[ScannedX[j]-Min_x][0] = ScannedVX[j];
                    Min_VBound[ScannedX[j]-Min_x][1] = ScannedVY[j];
                    Min_VBound[ScannedX[j]-Min_x][2] = ScannedVZ[j];
                    Min_VBound[ScannedX[j]-Min_x][3] = ScannedVW[j];
                }
            }
        }
    }
    
    for(int i = 0; i < YBoundLength; i++){
        int cur_x = Min_x + i;
        GLdouble zstep = (Max_ZBound[i] - Min_ZBound[i]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble rstep = (Max_CBound[i][0] - Min_CBound[i][0]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble gstep = (Max_CBound[i][1] - Min_CBound[i][1]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble bstep = (Max_CBound[i][2] - Min_CBound[i][2]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble nxstep = (Max_NBound[i][0] - Min_NBound[i][0]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble nystep = (Max_NBound[i][1] - Min_NBound[i][1]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble nzstep = (Max_NBound[i][2] - Min_NBound[i][2]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble nwstep = (Max_NBound[i][3] - Min_NBound[i][3]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble vxstep = (Max_VBound[i][0] - Min_VBound[i][0]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble vystep = (Max_VBound[i][1] - Min_VBound[i][1]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble vzstep = (Max_VBound[i][2] - Min_VBound[i][2]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble vwstep = (Max_VBound[i][3] - Min_VBound[i][3]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble lxstep = (Max_LBound[i][0] - Min_LBound[i][0]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble lystep = (Max_LBound[i][1] - Min_LBound[i][1]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble lzstep = (Max_LBound[i][2] - Min_LBound[i][2]) / (Max_YBound[i] - Min_YBound[i] + 1);
        GLdouble lwstep = (Max_LBound[i][3] - Min_LBound[i][3]) / (Max_YBound[i] - Min_YBound[i] + 1);
        
        for(int j = Min_YBound[i]; j <= Max_YBound[i]; j++){
            GLdouble cur_z = Min_ZBound[i] + zstep * (j-Min_YBound[i]);
            GLdouble cur_r = Min_CBound[i][0] + rstep * (j-Min_YBound[i]);
            GLdouble cur_g = Min_CBound[i][1] + gstep * (j-Min_YBound[i]);
            GLdouble cur_b = Min_CBound[i][2] + bstep * (j-Min_YBound[i]);
            GLdouble cur_nx = Min_NBound[i][0] + nxstep * (j-Min_YBound[i]);
            GLdouble cur_ny = Min_NBound[i][1] + nystep * (j-Min_YBound[i]);
            GLdouble cur_nz = Min_NBound[i][2] + nzstep * (j-Min_YBound[i]);
            GLdouble cur_nw = Min_NBound[i][3] + nwstep * (j-Min_YBound[i]);
            GLdouble cur_vx = Min_VBound[i][0] + vxstep * (j-Min_YBound[i]);
            GLdouble cur_vy = Min_VBound[i][1] + vystep * (j-Min_YBound[i]);
            GLdouble cur_vz = Min_VBound[i][2] + vzstep * (j-Min_YBound[i]);
            GLdouble cur_vw = Min_VBound[i][3] + vwstep * (j-Min_YBound[i]);
            GLdouble cur_lx = Min_LBound[i][0] + lxstep * (j-Min_YBound[i]);
            GLdouble cur_ly = Min_LBound[i][1] + lystep * (j-Min_YBound[i]);
            GLdouble cur_lz = Min_LBound[i][2] + lzstep * (j-Min_YBound[i]);
            GLdouble cur_lw = Min_LBound[i][3] + lwstep * (j-Min_YBound[i]);
            GLdouble L[4] = {cur_lx, cur_ly, cur_lz, cur_lw};
            GLdouble V[4] = {cur_vx, cur_vy, cur_vz, cur_vw};
            GLdouble N[4] = {cur_nx, cur_ny, cur_nz, cur_nw};
            GLdouble C[3] = {cur_r, cur_g, cur_b};
            gouraudShad(L, V, N, C);
            if(isZBUFFER && isInBound(j, 0, Z_HEIGHT) && isInBound(cur_x, 0, Z_WIDTH) && (cur_z <= zbuffer[j][cur_x])){ // >=   
                if(isInterNormal){
                    writepixelfast(cur_x, j, C[0], C[1], C[2]);
                }
                else{
                    writepixelfast(cur_x, j, cur_r, cur_g, cur_b);
                }
                zbuffer[j][cur_x] = cur_z;
            }
            else if(!isZBUFFER){
                if(isInterNormal){
                    writepixelfast(cur_x, j, C[0], C[1], C[2]);
                }
                else{
                    writepixelfast(cur_x, j, cur_r, cur_g, cur_b);
                }
            }
        }
    }
    
    delete [] Max_YBound;
    delete [] Min_YBound;
    delete [] Max_ZBound;
    delete [] Min_ZBound;
    /*for(int i = 0; i < YBoundLength; i++){
        //delete Max_CBound[i];
        //delete Min_CBound[i];
        //delete Max_VBound[i];
        //delete Min_VBound[i];
        //delete Max_NBound[i];
        //delete Min_NBound[i];
        //delete Max_LBound[i];
        //delete Min_LBound[i];
    }*/
    delete [] Max_CBound;
    delete [] Min_CBound;
    delete [] Max_NBound;
    delete [] Min_NBound;
    delete [] Max_VBound;
    delete [] Min_VBound;
    delete [] Max_LBound;
    delete [] Min_LBound;
	return true;
}

void printZBuffer(){
    FILE *fp = fopen("z.txt", "w");
    for(int i = 0; i < Z_HEIGHT; i++){
        for(int j = 400; j < Z_WIDTH; j++){
            fprintf(fp, "%10.2f ", zbuffer[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void toggleInterNormal(){
    isInterNormal = !isInterNormal;    
}

void toggleRealGouraud(){
    isRealGouraud = !isRealGouraud;
}
