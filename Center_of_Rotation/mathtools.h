#ifndef MATHTOOLS_H
#define MATHTOOLS_H

/* defs for mathtools.c */

typedef struct {		         /* x,y,z vector coords */
  double x;
  double y;
  double z;
} VECTOR;

void icr();
void norm();
void cross();
void copy_vector();
void add_vector();
void sub_vector();
void build_matrix();
void build_matrix2();
void strip_matrix();
void mult_matrix();
void build_vector();
void strip_vector();
void matxvec();
void vecxmat();
void vecxmat2();
void transpose();
void copy_matrix();
void vproj();
double vproj_value();
double dot();
double distance();
double matrix_det();
double vmag();
double vangle();
void sxvector();
void build_cs();
void build_vector_cs();
void find_centroid();
void print_matrix();
float find_triangle_area(double[3], double[3], double[3]);
void find_centroid(double[3], double[3], double[3], double[3]);
int ray_triangle_intersection();
void plane_triangle_intersection();
void rot_to_dcm();
void dcm_to_rot();
double point_to_line();
VECTOR ArbitraryRotate2();
void screw_angles();
#endif


































































