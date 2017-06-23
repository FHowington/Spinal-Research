/*
  mathtools.c - useful functions for manipulating vectors and matrices.
  vectors are in the form of VECTOR structures.  Matrices are 3X3 double
  arrays (e.g. m[row][col]).  All parameters are passed as pointers.

  The file "mathtools.h" must be #include'd in any program calling
  these routines.

  Routines:
  
  vangle       - returns the angle between two vectors
  vmag         - returns the magnitude of a vector
  dot          - returns dot product of 2 vectors
  distance     - returns distance between two vector endpoints
  norm         - returns unit vector in the direction of the input vector
  vproj        - calculates the projection of one vector onto another vector
  cross        - vector cross product
  copy_vector  - copies one vector to another
  add_vector   - adds two vectors
  sub_vector   - subtracts one vector from another
  sxvector     - multiplies a vector by a scalar
  matrix_det   - returns determinant of a 3x3 matrix
  transpose    - determines transpose of a 3x3 matrix
  copy_matrix  - copies one matrix to another
  mult_matrix  - multiplies two 3x3 matrices
  matxvec      - premultiplies a vector structure by a 3x3 matrix
  vecxmat      - postmultiplies a vector structure by a 3x3 matrix
  build_vector - creates a vector strucure from a 3-el. array
  strip_vector - creates a 3-el. array from a vector structure
  build_matrix - creates a 3x3 matrix from three vectors
  strip_matrix - extracts 3 column vectors from a matrix
  convert_vector- converts the location of a vector from 1 coordinate system to another
  vproj_value  - project one vector onto another, return + or - mag. of projection depending on direction
  
**********************************************************************/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mathtools.h"


/**********************************************************************/

double vangle(v1, v2)

VECTOR *v1, *v2;

{
  double angle, cosine,a,b,c;
  if(vmag(v1)== 0.0000 || vmag(v2) == 0.0000) {
    angle=0.00;
//    printf("zero angle");
  }
  else{
    a=dot(v1,v2);
    b=vmag(v1);
    c=vmag(v2);
    cosine=dot(v1,v2)/vmag(v1)/vmag(v2);
    if(cosine>0.9999999) cosine=1.;
    if(cosine<-0.9999999) cosine=-1.;

    angle=acos(cosine);
   if(angle<0.0)angle=-angle;
  }
  return (angle);
}

/**********************************************************************/

double vmag(v)

VECTOR *v;

{
  double vm;
  vm = sqrt(dot(v,v));
  if (vm < 0.0001)
    return (0);
  return (vm);
}

/**********************************************************************/

double dot(v1,v2)

VECTOR *v1;
VECTOR *v2;

{
  return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z);
}

/**********************************************************************/

double distance(v1,v2)

VECTOR *v1;
VECTOR *v2;

{
  return (sqrt((v2->x-v1->x)*(v2->x-v1->x) + (v2->y-v1->y)*(v2->y-v1->y) +
	       (v2->z-v1->z)*(v2->z-v1->z)));
}

/**********************************************************************/

void norm(v,normv)

VECTOR *v;
VECTOR *normv;

{
  double mag;

  mag = sqrt(dot(v,v));
  if (mag < .0001)
    return;
  normv->x = v->x/mag;
  normv->y = v->y/mag;
  normv->z = v->z/mag;
}

/**********************************************************************/
void vproj(v1,v2,v3)

VECTOR *v1;
VECTOR *v2;
VECTOR *v3;

{
  double theta, b;

  theta=vangle(v1,v2);
  norm(v2,v3);
  b=vmag(v1)*cos(theta);
  v3->x = v3->x*b;
  v3->y = v3->y*b;
  v3->z = v3->z*b;
}

/**********************************************************************/
double vproj_value(v1,v2)

VECTOR *v1;
VECTOR *v2;

{
  VECTOR v4;

  vproj(v1,v2,&v4);
  if(dot(v1,v2)>0.0)
    return(vmag(&v4));
  else
    return(-vmag(&v4));
}


/**********************************************************************/

void cross(v1,v2,vc)

VECTOR *v1;
VECTOR *v2;
VECTOR *vc;

{
  vc->x = v1->y*v2->z - v1->z*v2->y;
  vc->y = v1->z*v2->x - v1->x*v2->z;
  vc->z = v1->x*v2->y - v1->y*v2->x;
}


/***********************************************************************

  copy_vector - copies one vector structure to another: vout = v1

*/

void copy_vector(v1,vout)

VECTOR *v1;
VECTOR *vout;

{
  vout->x = v1->x;
  vout->y = v1->y;
  vout->z = v1->z;
}

/***********************************************************************

  add_vector - adds two vector structures: vout = v1 + v2

*/

void add_vector(v1,v2,vout)

VECTOR *v1,*v2;
VECTOR *vout;

{
  vout->x = v1->x + v2->x;
  vout->y = v1->y + v2->y;
  vout->z = v1->z + v2->z;
}


/***********************************************************************

  sub_vector - subtracts one vector structure from another: vout = v1 - v2

*/

void sub_vector(v1,v2,vout)

VECTOR *v1,*v2;
VECTOR *vout;

{
  vout->x = v1->x - v2->x;
  vout->y = v1->y - v2->y;
  vout->z = v1->z - v2->z;
}


/***********************************************************************

  sxvector - multiplies a vector structure by a scalar: vout = s*v1

*/

void sxvector(s,v,vout)

double s;
VECTOR *v;
VECTOR *vout;

{
  vout->x = s*v->x;
  vout->y = s*v->y;
  vout->z = s*v->z;
}


/***********************************************************************

  transpose - determines transpose of a 3x3 matrix

*/

void transpose(m,mout)

double m[3][3];
double mout[3][3];

{
  int i, j;

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      mout[i][j] = m[j][i];
}





/***********************************************************************

  matrix_det - returns the determinant of a 3x3 matrix

*/

double matrix_det(a)

double a[3][3];

{
  return(a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] +
	 a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] -
	 a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0]);
}

/***************************************************************************

  mult_matrix - multiplies two 3x3 matrices: mprod = m1xm2

*/

void mult_matrix(m1,m2,mprod)

double m1[3][3];
double m2[3][3];
double mprod[3][3];

{
  int i;
  int nrow, ncol;

  for (nrow=0; nrow<3; nrow++)
    for (ncol=0; ncol<3; ncol++)
      mprod[nrow][ncol] = m1[nrow][0]*m2[0][ncol] +
	                  m1[nrow][1]*m2[1][ncol] +
	                  m1[nrow][2]*m2[2][ncol];
}

/***************************************************************************

  copy_matrix - copies contents of one 3x3 matrix into another: mout = m1

*/

void copy_matrix(m1,mout)

double m1[3][3];
double mout[3][3];

{
  int nrow, ncol;

  for (nrow=0; nrow<3; nrow++)
    for (ncol=0; ncol<3; ncol++)
      mout[nrow][ncol] = m1[nrow][ncol];
}


/***************************************************************************

  matxvec - premultiplies a 3-el. vector structure by a 3x3 matrix
            vout = m x v  (v, vout are assumed to be column vectors)

*/

void matxvec(m,v,vout)

VECTOR *v;
VECTOR *vout;
double m[3][3];

{
  vout->x = v->x*m[0][0] + v->y*m[0][1] + v->z*m[0][2];
  vout->y = v->x*m[1][0] + v->y*m[1][1] + v->z*m[1][2];
  vout->z = v->x*m[2][0] + v->y*m[2][1] + v->z*m[2][2];
}


/***************************************************************************

  vecxmat - postmultiplies a 3-el. vector structure by a 3x3 matrix
            vout = v x m  (v, vout are assumed to be row vectors)

	    note: this is the same as vout = transpose(m) x transpose(v)

*/

void vecxmat(m,v,vout)

VECTOR *v;
VECTOR *vout;
double m[3][3];

{
  vout->x = v->x*m[0][0] + v->y*m[1][0] + v->z*m[2][0];
  vout->y = v->x*m[0][1] + v->y*m[1][1] + v->z*m[2][1];
  vout->z = v->x*m[0][2] + v->y*m[1][2] + v->z*m[2][2];
}

/***************************************************************************/


/***************************************************************************

  vec2xmat - use this function when you know the location of a 3D point (or vector) in
             a coordinate system and the coordinate system changes.
             Input the point location relative to the coordinate system (v)
             and the orthogonal coordinate system.  The location of the point 
             in the new coordinate system is output.

             -- --     -- --     --       --
             | x |     | x |     | 00 01 02 |
       vout  | y | = v | y | x m | 10 11 12 |
             | z |     | z |     | 20 21 22 |
             -- --     -- --     --        --
            
             where v[x,y,z] is the loction of the point in the original coordinate system
             and m[0][0], m[1][0], m[2][0] is the unit vector in the x direction
*/



/*****************************************************/

void vecxmat2(m,v,vout)

VECTOR *v;
VECTOR *vout;
double m[3][3];

{
  vout->x = v->x*m[0][0] + v->y*m[0][1] + v->z*m[0][2];
  vout->y = v->x*m[1][0] + v->y*m[1][1] + v->z*m[1][2];
  vout->z = v->x*m[2][0] + v->y*m[2][1] + v->z*m[2][2];
}


/***************************************************************************

  build_vector - loads values from v[3] into a vector structure

*/

void build_vector(v,vout)

double v[3];
VECTOR *vout;

{
  vout->x = v[0];
  vout->y = v[1];
  vout->z = v[2];
}

/***************************************************************************

  strip_vector - loads values from  a vector structure into vout[3]

*/

void strip_vector(v,vout)

double vout[3];
VECTOR *v;

{
  vout[0] = v->x;
  vout[1] = v->y;
  vout[2] = v->z;
}


/***************************************************************************

  build_matrix - creates a matrix (double 3x3 array) whose columns consist
  of three input vectors in the order supplied.

*/

void build_matrix(v1,v2,v3,m)

VECTOR *v1, *v2, *v3;
double m[3][3];

{
  m[0][0] = v1->x;
  m[1][0] = v1->y;
  m[2][0] = v1->z;

  m[0][1] = v2->x;
  m[1][1] = v2->y;
  m[2][1] = v2->z;

  m[0][2] = v3->x;
  m[1][2] = v3->y;
  m[2][2] = v3->z;
}

/***************************************************************************

  build_matrix2 - creates a matrix (double 3x3 array) whose rows consist
  of three input vectors in the order supplied.

*/

void build_matrix2(v1,v2,v3,m)

VECTOR *v1, *v2, *v3;
double m[3][3];

{
  m[0][0] = v1->x;
  m[1][0] = v2->x;
  m[2][0] = v3->x;

  m[0][1] = v1->y;
  m[1][1] = v2->y;
  m[2][1] = v3->y;

  m[0][2] = v1->z;
  m[1][2] = v2->z;
  m[2][2] = v3->z;
}



/***************************************************************************

  strip_matrix - creates 3 vectors from the columns of a matrix (double 3x3
  array)

*/

void strip_matrix(v1,v2,v3,m)

VECTOR *v1, *v2, *v3;
double m[3][3];

{
  v1->x =  m[0][0];
  v1->y =  m[0][1];
  v1->z =  m[0][2];
  
  v2->x = m[1][0];
  v2->y = m[1][1];
  v2->z = m[1][2];
  
  v3->x = m[2][0];
  v3->y = m[2][1];
  v3->z = m[2][2];
}

void strip_matrix2(v1,v2,v3,m)

VECTOR *v1, *v2, *v3;
double m[3][3];

{
  v1->x =  m[0][0];
  v1->y =  m[1][0];
  v1->z =  m[2][0];
  
  v2->x = m[0][1];
  v2->y = m[1][1];
  v2->z = m[2][1];
  
  v3->x = m[0][2];
  v3->y = m[1][2];
  v3->z = m[2][2];
}

/**********************************************************************

  convert_vector - converts a vector from 1 coordinate system 
                   to another coordiante system

  v1: origin of final CS
  v2: vector to be converted to final CS
  m:  xyz directions of final CS, x-axis column 0, y-axis column 1, z-axis coulmn 2
  vout: v2 in final CS 

*/

void convert_vector(v1,v2,m,vout)

double v1[3], v2[3], vout[3];
double m[3][3];

{
  int xyz;
  double p[3];
  VECTOR r1, r2;

  for(xyz=0;xyz<3;xyz++)
    p[xyz]=v2[xyz]-v1[xyz];
  build_vector(p,&r1);
  vecxmat(m,&r1,&r2); 
  strip_vector(&r2,vout);
}

/**********************************************************************

  find_centroid - given the 3 vertices of a triangle, 
  finds the centroid of the triangle

  v1: vertex 1
  v2: vertex 2
  v3: vertex 3

*/

void find_centroid(v1,v2,v3,center)

double v1[3], v2[3], v3[3];
double center[3];

{
  int xyz;
  
  for(xyz=0;xyz<3;xyz++)
    center[xyz]=(v1[xyz]+v2[xyz]+v3[xyz])/3.;

//  printf("%6.4f,%6.4f,%6.4f\n",center[0],center[1],center[2]);
  

}

/**********************************************************************

  find_triangle_area - given the 3 vertices of a triangle, 
  finds the area of the triangle

  v1: vertex 1
  v2: vertex 2
  v3: vertex 3
  area: area

  algorithm: Sides of triangle are lengths a,b,c.  s is (a+b+c)/2
  Area is sqrt[s*(s-a)(s-b)*(s-c)]

*/

float find_triangle_area(v1,v2,v3)

double v1[3], v2[3], v3[3];


{
  VECTOR r1, r2, r3, vout1, vout2, vout3;
  double a, b, c, s;
  float area;
	
  build_vector(v1,&r1);
  build_vector(v2,&r2);
  build_vector(v3,&r3);
  sub_vector(&r1,&r2,&vout1);
  sub_vector(&r2,&r3,&vout2);
  sub_vector(&r1,&r3,&vout3);
  a=vmag(&vout1);
  b=vmag(&vout2);
  c=vmag(&vout3);
  s=(a+b+c)/2.;

  area=sqrt(s*(s-a)*(s-b)*(s-c));
/*  printf("%6.2f,%6.2f%6.2f\n",v1[0],v1[1],v1[2]);
  printf("%6.2f,%6.2f%6.2f\n",v2[0],v2[1],v2[2]);
  printf("%6.2f,%6.2f%6.2f\n",v3[0],v3[1],v3[2]);
  printf("%6.4f\n",area);*/
  return (area);


}


/**********************************************************************

  print_matrix - prints out a matrix (double 3x3 array)

*/

void print_matrix(m)

double m[3][3];

{
  int  nrow;

  for(nrow=0;nrow<3;nrow++)
    printf("\n %f, %f, %f",m[nrow][0] ,m[nrow][1] ,m[nrow][2] );
  printf("\n");
}


/**********************************************************************

  vector_triangle_intersection - find the intersection of a ray and triangle

  v1: start point for ray
  v2: direction of ray
  v3: 1st vertex of triangle
  v4: 2nd vertex of triangle
  v5: 3rd vertex of triangle
  vout: intersection of ray and triangle

  algorithm: adopted from http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

*/

int ray_triangle_intersection(v1, v2, v3, v4, v5, vout )

VECTOR *v1, *v2, *v3, *v4, *v5, *vout;


{
  VECTOR u, v, n;
  VECTOR dir, w0, w, ray_mag, ray;
  double r, a, b;
  double uu, uv, vv, wu, wv, D;
  float s, t;
	
//	printf("\t%6.2f %6.2f %6.2f\n",v4->x,v4->y,v4->z);
  sub_vector(v4, v3, &u); // get triangle edges
  sub_vector(v5, v3, &v);
  cross(&u,&v,&n); // normal to the triangle
  if(vmag(&n) < 0.000001)
	return (-1);
	
   copy_vector(v2,&dir); // ray direction vector
   sub_vector(v1,v3,&w0);
   a = -dot(&n,&w0);
   b = dot(&n,&dir);
   
//   printf("\t%6.2f %6.2f\t",a,b);
   
   if(fabs(b) < 0.0000001){ // ray is parallel to triangle plane
	  if(a == 0) // ray lies in triangle plane
		return (2);
	  else
		return (0);
	}
	// get intersect point of ray with triangle plane
	r = a/b;
	if(r<0.0) // ray goes away from triangle
		return (0); // no intersect
	
	sxvector(r,&dir,&ray_mag);
	add_vector(v1,&ray_mag,vout); // intersect point of ray and plane
//	printf("\t%6.2f %6.2f %6.2f\n",vout->x,vout->y,vout->z);		
	// test to see if vout is inside the triangle
	
	uu = dot(&u,&u);
	uv = dot(&u,&v);
	vv = dot(&v,&v);
	sub_vector(vout,v3,&w);
	wu = dot(&w,&u);
	wv = dot(&w,&v);
	D = uv*uv - uu*vv;
	
	// get and test parametric coords
	s = (uv*wv - vv*wu)/D;
	if (s < 0.0 || s > 1.0) // vout is outside triangle
		return (0);
	t = (uv*wu - uu*wv)/D;
	if(t<0.0 || (s+t)>1.0) // vout is outside triangle
		return (0);

	return (1);   // vout is in triangle
		
	

}

/**********************************************************************

  vpoint_to_line - find the perpendicular distance from a point to a line

  p1: one point for line
  p2: second point on line
  p3: a point

*/

double point_to_line(p1, p2, p3)

VECTOR *p1, *p2, *p3;

{
	VECTOR v1, v2, v3, v4, v5;
	double d;
	
	sub_vector(p2,p1,&v1); // vector in direction of the line
	sub_vector(p3,p1,&v2); // vector from point 1 on line to the point
	vproj(&v2,&v1,&v3);
    sub_vector(&v2,&v3,&v4);
    return(vmag(&v4));
    /*
	norm(&v1,&v1);
	sxvector(vmag(&v3),&v1,&v4);
	add_vector(p1,&v4,&v5);
//	printf("%6.3f\t%6.3f\t%6.3f\n",p3->x,p3->y,p3->z);
//	printf("from subroutine\t");
	d=distance(&v5,p3);

//	printf("%6.3f\t%6.3f\n",d,distance(&v5,p3));
	return(d);
     */
}

/*********************************************
  plane_triangle_intersection - find the intersection points between a plane and the edges of a triangle
 
  inputs:
    v1: first point making up triangle
    v2: second point making up triangle
    v3: third point making up triangle
    v4: a point on the plane
    vnorm: vector normal to the plane
  outputs:
    start_point: the start of the intersection line segment 
    end_point: the end of the intersection line segment
 **********************************************/

void plane_triangle_intersection(v1, v2, v3, v4, vnorm, start_point, end_point)


VECTOR *v1, *v2, *v3, *v4, *vnorm, *start_point, *end_point;
#define ABSVAL(a)  (a>0. ? (a):(-a))

{
    double dp1, dp2, dp3, dd, n, si;
    VECTOR start1, start2, end1, u_vector, w_vector, u_vector_norm, w_vector_norm, u_vector_out, segment_start, segment_end;
    
    sub_vector(v1,v4,v1); // build a vector from the point on the plane to each vertex of the triangle
    sub_vector(v2,v4,v2);
    sub_vector(v3,v4,v3);
    dp1=dot(v1,vnorm); // find the dot product between each vector and the normal to the plane
    dp2=dot(v2,vnorm);
    dp3=dot(v3,vnorm);
    // first identify which markers lie on opposite sides of the plane using "if" loop below
    
    if(dp1>0 && dp2>0 && dp3<0){
        copy_vector(v1,&start1);
        copy_vector(v2,&start2);
        copy_vector(v3,&end1);
    }
    else if (dp1>0 && dp3>0 && dp2<0){
        copy_vector(v1,&start1);
        copy_vector(v3,&start2);
        copy_vector(v2,&end1);
    }
    else if (dp2>0 && dp3>0 && dp1<0){
        copy_vector(v2,&start1);
        copy_vector(v3,&start2);
        copy_vector(v1,&end1);
    }
    else if(dp1<0 && dp2<0 && dp3>0){
        copy_vector(v1,&start1);
        copy_vector(v2,&start2);
        copy_vector(v3,&end1);
    }
    else if (dp1<0 && dp3<0 && dp2>0){
        copy_vector(v1,&start1);
        copy_vector(v3,&start2);
        copy_vector(v2,&end1);
    }
    else if (dp2<0 && dp3<0 && dp1>0){
        copy_vector(v2,&start1);
        copy_vector(v3,&start2);
        copy_vector(v1,&end1);
    }
    else
        printf("cannot find intersection of triangle and plane\n");
        //find the first endpoint of the line created by the triangle/plane intersection
        sub_vector(&end1,&start1,&u_vector);
        copy_vector(&start1,&w_vector);
        norm(&u_vector,&u_vector_norm);
        norm(&w_vector,&w_vector_norm);
        dd = dot(vnorm,&u_vector_norm);//
        n = -dot(vnorm,&w_vector);//
        if(ABSVAL(dd)<.000001)
            printf("segment is parallel to the plane or lies in the plane\n");
        else{
            si = n/dd;
            sxvector(si,&u_vector_norm,&u_vector_out);
            add_vector(&start1,&u_vector_out,&segment_start);
            add_vector(&segment_start,v4,&segment_start);
            copy_vector(&segment_start,start_point);
                //				printf("Bone: %i triangle %i start_count %i start [%6.2f, %6.2f, %6.2f] \n",bone,left_foramen_plane_loop_triangle[bone][i],marker_count,segment_start.x,segment_start.y,segment_start.z);
        }
    
    
    //find the second endpoint of the line created by the triangle/plane intersection
    sub_vector(&end1,&start2,&u_vector);
    copy_vector(&start2,&w_vector);
    norm(&u_vector,&u_vector_norm);
    norm(&w_vector,&w_vector_norm);
    dd = dot(vnorm,&u_vector_norm);
    n = -dot(vnorm,&w_vector);
    if(ABSVAL(dd)<.000001)
        printf("segment is parallel to the plane or lies in the plane\n");
    else{
        si = n/dd;
        sxvector(si,&u_vector_norm,&u_vector_out);
        add_vector(&start2,&u_vector_out,&segment_end);
        add_vector(&segment_end,v4,&segment_end);
        copy_vector(&segment_end,end_point);
            //				printf("Bone: %i triangle %i end_count %i  end [%6.2f, %6.2f, %6.2f]\n",bone,left_foramen_plane_loop_triangle[bone][i],marker_count,segment_end.x,segment_end.y,segment_end.z);
    }

}
/*
 icr.c
 
 This program calculates the instant center of roation parameters following Spoor and Veldpaus (J. Biomechanics, 1980)
 
 INPUT
 +the desired plane of intersection (type=1 = xy, type=2 = yz, type=3 = xz), the XYZ CS of the "moving" bone relative to the "fixed" bone at the "start" frame, the XYZ CS of the "moving" bone relative to the "fixed" bone at the "end" frame, the origin of the "moving" bone at the "start" frame, the origin of the "moving" bone at the "end" frame
 
 OUTPUT
 +the amount of rotation between bodies (phi), the direction of the helical axis (screw_vec), the translation along the helical axis (t), and the intersection of the helical axis and the given plane in the coordinate system of the "fixed" bone (intersection_point)
 */

void icr(int type, double BONE1A_R_BONE2[3][3], double BONE1B_R_BONE2[3][3], double locationA[3], double locationB[3], double *phi, VECTOR *screw_vec, double *t, double intersection_point[3])

{
    int xyz;
    
    double translation[3], RT[3][3], rot[3][3], sin_phi, normal[3], scale;
    
    VECTOR v1, v2, v3, v4, v6, v7, trans, pt, point_rotated, point_transformed, s_star;
    
       // print_matrix(BONE1A_R_BONE2);
       // print_matrix(BONE1B_R_BONE2);
       // printf("locationA: [%6.2f,%6.2f,%6.2f]",locationA[0], locationA[1], locationA[2]);
       // printf("locationB: [%6.2f,%6.2f,%6.2f]",locationB[0], locationB[1], locationB[2]);
    
    transpose(BONE1A_R_BONE2, RT);
    mult_matrix(BONE1B_R_BONE2, RT, rot);
    
    //    print_matrix(rot);
    
    for(xyz=0;xyz<3;xyz++)
        translation[xyz]=locationB[xyz]-locationA[xyz];
    build_vector(translation,&trans);
    
    
    sin_phi = 0.5* sqrt(pow((rot[2][1]-rot[1][2]),2) + pow((rot[0][2]-rot[2][0]),2) + pow((rot[1][0]-rot[0][1]),2));
    *phi = asin(sin_phi);
    
    screw_vec->x = 0.5*(rot[2][1]-rot[1][2])/sin_phi;
    screw_vec->y = 0.5*(rot[0][2]-rot[2][0])/sin_phi;
    screw_vec->z = 0.5*(rot[1][0]-rot[0][1])/sin_phi;
    
    //    printf("screw vec: [%6.2f,%6.2f,%6.2f]",screw_vec->x, screw_vec->y, screw_vec->z);
    //    printf("translation: [%6.2f,%6.2f,%6.2f]",translation[0], translation[1], translation[2]);
    *t = dot(screw_vec, &trans);
    //    printf("t: %6.2f",*t);
    cross(screw_vec,&trans,&v1);
    cross(screw_vec,&v1,&v2);
    sxvector(-0.5,&v2,&v3);
    sxvector(sin(*phi)/(2.*(1.-cos(*phi))),&v1,&v4);
    add_vector(&v3,&v4,&pt);
//    printf("s: [%6.2f,%6.2f,%6.2f]\t",pt.x,pt.y,pt.z);
    norm(screw_vec,screw_vec);
//    printf("screw_vec1: [%6.2f,%6.2f,%6.2f]\t",screw_vec->x,screw_vec->y,screw_vec->z);
    
    // convert point on helical axis to inferior bone CS
    vecxmat(BONE1A_R_BONE2,&pt,&point_rotated);
//    printf("point rotated: [%6.2f,%6.2f,%6.2f]\t",point_rotated.x,point_rotated.y,point_rotated.z);

    build_vector(locationA,&v7);
//    printf("locationA: [%6.2f,%6.2f,%6.2f]\t",v7.x,v7.y,v7.z);

    add_vector(&point_rotated,&v7,&point_transformed);
//    printf("point transformed: [%6.2f,%6.2f,%6.2f]\t",point_transformed.x,point_transformed.y,point_transformed.z);
    //
    // I think screw_vec also needs to be converted into inferior bone CS
    //
    vecxmat(BONE1A_R_BONE2,screw_vec,screw_vec);
//    printf("screw_vec2: [%6.2f,%6.2f,%6.2f]\n",screw_vec->x,screw_vec->y,screw_vec->z);
   
    // calculate the intersection of a vector and a plane

    if (type==3){
        normal[0]=0.0;
        normal[1]=0.0; 
        normal[2]=1.0; 
    }
    else if (type==1){
        normal[0]=1.0;
        normal[1]=0.0;
        normal[2]=0.0;
    }
    else{
        normal[0]=0.0;
        normal[1]=1.0; // should be 1.0- modified for c-spine to 0.0
        normal[2]=0.0;// should be 0.0- modified for c-spine to 1.0
    }
    
    build_vector(normal,&v6);
    scale=-dot(&v6,&point_transformed)/dot(&v6,screw_vec); // this equation from http://geomalgorithms.com/
    sxvector(scale,screw_vec,&v7);
    add_vector(&point_transformed,&v7,&s_star);  // s_star is the intersection of the helical axis and plane of interest
    //    printf("intersection point: [%6.2f,%6.2f,%6.2f]\n",s_star.x,s_star.y,s_star.z);
    strip_vector(&s_star,intersection_point);
    
    
}


/**********************************************************************
 
 screw_angles - find the angles between the screw axis (or any vector) and the XY, YZ and XZ planes
 
 v1: the input vector
 xy_ang: the angle between v1 and the xy-plane
 yz_ang: the angle between v1 and the yz-plane
 xz_ang: the angle between v1 and the xz-plane
 
 */

void screw_angles(v1, xy_ang, yz_ang, xz_ang)

VECTOR *v1;
double *xy_ang, *yz_ang, *xz_ang;

{
    double s[3], sign, pi=3.1415926535;
    VECTOR vnew;
    
	strip_vector(v1,s);
//    printf("%6.3f\t%6.3f\t%6.3f\t",s[0],s[1],s[2]);
    sign = s[1] > 0.00 ? 1.0:-1.0;
    s[1]=0.000;
    build_vector(s,&vnew);
    *xz_ang=sign*vangle(v1,&vnew)*180./pi;

    strip_vector(v1,s);
    sign = s[0] > 0.00 ? 1.0:-1.0;
    s[0]=0.000;
    build_vector(s,&vnew);
    *yz_ang=sign*vangle(v1,&vnew)*180./pi;

    strip_vector(v1,s);
    sign = s[2] > 0.00 ? 1.0:-1.0;
    s[2]=0.000;
    build_vector(s,&vnew);
    *xy_ang=sign*vangle(v1,&vnew)*180./pi;
//    printf("%6.3f\t%6.3f\t%6.3f\n",*xy_ang,*yz_ang,*xz_ang);

}























