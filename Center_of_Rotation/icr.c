 /* 
icr.c
  
  This program calculates the instant center of roation parameters following Spoor and Veldpaus (J. Biomechanics, 1980)
  
  INPUT 
+the desired plane of intersection (type=1 = xy, type=2 = yz, type=3 = xz), the XYZ CS of the "moving" bone, the XYZ CS of the "fixed" bone, the center of the "moving" bone, the center of the "fixed" bone
  
  OUTPUT
+the amount of rotation between bodies (phi), the direction of the helical axis (screw_vec), the translation along the helical axis (t), and the intersection of the helical axis and the given plane in the coordinate system of the "fixed" bone (intersection_point) 
*/
#include <math.h>
#include "mathtools.h"


  void icr(int type, double BONE1A_R_BONE2[3][3], double BONE1B_R_BONE2[3][3], double locationA[3], double locationB[3], double *phi, VECTOR *screw_vec, double *t, double intersection_point[3])

{
    int xyz;
    
    double translation[3], RT[3][3], rot[3][3], sin_phi, normal[3], scale;
    
    VECTOR v1, v2, v3, v4, v6, v7, trans, pt, point_rotated, point_transformed, s_star;
    
//    print_matrix(BONE1A_R_BONE2);
//    print_matrix(BONE1B_R_BONE2);
//    printf("locationA: [%6.2f,%6.2f,%6.2f]",locationA[0], locationA[1], locationA[2]);
//    printf("locationB: [%6.2f,%6.2f,%6.2f]",locationB[0], locationB[1], locationB[2]);

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
//    printf("s: [%6.2f,%6.2f,%6.2f]\n",s->x,s->y,s->z);
    norm(screw_vec,screw_vec);
    
    // convert point on helical axis to inferior bone CS
    vecxmat(RT,&pt,&point_rotated);
    build_vector(locationA,&v7);
    add_vector(&point_rotated,&v7,&point_transformed);
//    printf("point transformed: [%6.2f,%6.2f,%6.2f]\n",point_transformed.x,point_transformed.y,point_transformed.z);
    
    if (type==1){
        normal[0]=0.0;
        normal[1]=0.0;
        normal[2]=1.0;
    }
    else if (type==2){
        normal[0]=1.0;
        normal[1]=0.0;
        normal[2]=0.0;
    }
    else{
        normal[0]=0.0;
        normal[1]=1.0;
        normal[2]=0.0;
    }
    
    build_vector(normal,&v6);
    scale=-dot(&v6,&point_transformed)/dot(&v6,screw_vec); // this equation from http://geomalgorithms.com/
    sxvector(scale,screw_vec,&v7);
    add_vector(&point_transformed,&v7,&s_star);  // s_star is the intersection of the helical axis and plane of interest
//    printf("intersection point: [%6.2f,%6.2f,%6.2f]\n",s_star.x,s_star.y,s_star.z);
    strip_vector(&s_star,intersection_point);
    
    
}
