#include "include.h"
#include "csvparser.h"

int main() {
    int i,j,k,frame;

    char C3_file[]="HomoTransMatrices_C123_-to-Lab.csv"; 
    char C4_file[]="HomoTransMatrices_C4_-to-Lab.csv";
    char C5_file[]="HomoTransMatrices_C5_-to-Lab.csv";
    char C6_file[]="HomoTransMatrices_C6_-to-Lab.csv";
    char C7_file[]="HomoTransMatrices_C7T1_-to-Lab.csv";
    char C3r_file[]="HomoTransMatrices_Lab-to-C123_.csv"; 
    char C4r_file[]="HomoTransMatrices_Lab-to-C4_.csv";
    char C5r_file[]="HomoTransMatrices_Lab-to-C5_.csv";
    char C6r_file[]="HomoTransMatrices_Lab-to-C6_.csv";
    char C7r_file[]="HomoTransMatrices_Lab-to-C7T1_.csv";

    double** cor(char *, char *, int *,int);
    double* flexion_all = malloc(4*500*sizeof(double*));
    FILE *dest;
    CsvRow *row;
    int *angles;
    double **results;

    dest = fopen("results.csv","wb");

    CsvParser *csvparser_angles = CsvParser_new("KinematicsMeasurementReport.csv", ",", 0);
    frame=0;
    while ((row = CsvParser_getRow(csvparser_angles)) ) {
        const char **rowFields = CsvParser_getFields(row);
        if(frame>1){
            flexion_all[4*frame]=atof(rowFields[2]);
            flexion_all[4*frame+1]=atof(rowFields[5]);
            flexion_all[4*frame+2]=atof(rowFields[8]);
            flexion_all[4*frame+3]=atof(rowFields[11]);
        }
        CsvParser_destroy_row(row);
        frame++;
    }
    CsvParser_destroy(csvparser_angles);


    // This is C3 relative to C4 
    fprintf(dest,"C3 to C4, X,Y,Z\n");

    // This finds the next frame with an absolute difference > 2 degrees
    angles = malloc(2*500*sizeof(int));
    k=0;
    for(i=2;i<frame;i++){
        for(j=i;j<frame;j++){
            if(flexion_all[4*i]+2 < flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-2 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C3_file,C4r_file,angles,k);
    for(i=0;i<k;i++){
        fprintf(dest, "%f,%f,%f,%f\n",flexion_all[4*(i+2)],results[i][0],results[i][1],results[i][2]); 
    }
    free(results);
    free(angles);


    // This is C4 relative to C5
    fprintf(dest,"\n\nC4 to C5, X,Y,Z\n");

    // This finds the next frame with an absolute difference > 2 degrees
    angles = malloc(2*500*sizeof(int));
    k=0;
    for(i=2;i<frame;i++){
        for(j=i;j<frame;j++){
            if(flexion_all[4*i+1]+2 < flexion_all[4*j+1]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-2 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C4_file,C5r_file,angles,k);
    for(i=0;i<k;i++){
        fprintf(dest, "%f,%f,%f,%f\n",flexion_all[4*(i+2)+1],results[i][0],results[i][1],results[i][2]); 
    }
    free(results);
    free(angles);



    // This is C5 relative to C6
    fprintf(dest,"\n\nC5 to C6, X,Y,Z\n");

    // This finds the next frame with an absolute difference > 2 degrees
    angles = malloc(2*500*sizeof(int));
    k=0;
    for(i=2;i<frame;i++){
        for(j=i;j<frame;j++){
            if(flexion_all[4*i+2]+2 < flexion_all[4*j+2]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-2 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C5_file,C6r_file,angles,k);
    for(i=0;i<k;i++){
        fprintf(dest, "%f,%f,%f,%f\n",flexion_all[4*(i+2)+2],results[i][0],results[i][1],results[i][2]); 
    }
    free(results);
    free(angles);


    // This is C6 relative to C7
    fprintf(dest,"\n\nC6 to C7, X,Y,Z\n");

    // This finds the next frame with an absolute difference > 2 degrees
    angles = malloc(2*500*sizeof(int));
    k=0;
    for(i=2;i<frame;i++){
        for(j=i;j<frame;j++){
            if(flexion_all[4*i+3]+2 < flexion_all[4*j+3]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-2 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C6_file,C7r_file,angles,k);
    for(i=0;i<k;i++){
        fprintf(dest, "%f,%f,%f,%f\n",flexion_all[4*(i+2)+3],results[i][0],results[i][1],results[i][2]); 
    }
    free(results);
    free(angles);

    return 0;
}


double** cor(char *file_1, char *file_2, int *angles, int count){
    int i = 0;
    int j = 0;
    int k = 0;
    int frame = 0;
    int type = 3;
    double total;
    int framecount;
    double bone1_CS[500][4][4];
    double bone2_CS[500][4][4];
    double bone12_CS[500][4][4];
    double bone12_rot[500][3][3];
    double position12[500][3];

    double phi;
    double t;
    double intersection[3];
    double** results = malloc(500*sizeof(double*));
    FILE *dest;
    VECTOR screw_vec;
    CsvRow *row;

    for (i=0; i<500; i++)
        results[i] = (double *)calloc(3,sizeof(double));

    CsvParser *csvparser1 = CsvParser_new(file_1, ",", 0); 

    while ((row = CsvParser_getRow(csvparser1)) ) {
        const char **rowFields = CsvParser_getFields(row);
        j=k=0;
        for (i = 2 ; i < CsvParser_getNumFields(row) ; i++) {
            //Need to filter out every fourth column beginning with first.
            if(k==4){
                k=0;
                j++;
            }
            bone1_CS[frame][j][k]=atof(rowFields[i]);
            k++;
        }
        frame++;
        CsvParser_destroy_row(row);
    }
    CsvParser_destroy(csvparser1);


    CsvParser *csvparser2 = CsvParser_new(file_2, ",", 0);
    frame=0;
    while ((row = CsvParser_getRow(csvparser2)) ) {
        const char **rowFields = CsvParser_getFields(row);
        j=k=0;
        for (i = 2 ; i < CsvParser_getNumFields(row) ; i++) {
            //Need to filter out every fourth column beginning with first.
            if(k==4){
                k=0;
                j++;
            }
            bone2_CS[frame][j][k]=atof(rowFields[i]);
            k++;
        }
        frame++;
        CsvParser_destroy_row(row);
    }

    CsvParser_destroy(csvparser2);

    // Performing matrix multiplication to transform into C4-C3 reference
    // This yields B*A. B should be Lab to C_lower
    // A should be C_higher to Lab
    for(framecount=0;framecount<frame;framecount++){
        for(i=0;i<4;i++){
            for(j=0;j<4;j++){
                total=0;
                for(k=0;k<4;k++){
                    total+=bone1_CS[framecount][k][j]*bone2_CS[framecount][i][k];
                }
                bone12_CS[framecount][i][j]=total;
            }
        }
    }


    // Split the matrix into its respective rotation matrix and translation vector
    for(framecount=0;framecount<frame;framecount++){
        for(i=0;i<4;i++)
            for(j=0;j<4;j++){
                if(i<3 && j<3){
                    bone12_rot[framecount][i][j]=bone12_CS[framecount][i][j];
                }
                else if(j<3 && i==3){
                    position12[framecount][j]=bone12_CS[framecount][i][j];
                }
            }
    }

    // At this point all data is loaded into memory in the appropriate format.
    // Let's output the rotation angles of bone1 relative to bone2
    
    for(i=0;i<count;i++){
        icr(type,bone12_rot[(angles[2*i]-1)],bone12_rot[(angles[2*i+1]-1)],position12[(angles[2*i]-1)], position12[(angles[2*i+1]-1)],&phi,&screw_vec,&t,intersection);

        for(j=0;j<3;j++){
            results[i][j]=intersection[j];
        }
    }
    return results;
}
