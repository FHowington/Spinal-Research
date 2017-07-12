#include "include.h"
#include "csvparser.h"

double Determinant(double **,int);
void CoFactor();

int main() {
    int i,j,k,frame;

    char C3_file[]="HomoTransMatrices_C123_-to-Lab.csv"; 
    char C4_file[]="HomoTransMatrices_C4_-to-Lab.csv";
    char C5_file[]="HomoTransMatrices_C5_-to-Lab.csv";
    char C6_file[]="HomoTransMatrices_C6_-to-Lab.csv";
    char C7_file[]="HomoTransMatrices_C7T1_-to-Lab.csv";

    double** cor(char *, char *,int *,int);
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
            if(flexion_all[4*i]+4 < flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                //printf("\nAngles are %d\t%d",i,j);
                k++;
                break;
            }
            if(flexion_all[4*i]-4 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C3_file,C4_file,angles,k);
    for(i=0;i<k;i++){
        fprintf(dest, "%f,%f,%f,%f\n",flexion_all[4*(i+2)],results[i][0],results[i][1],results[i][2]); 
    }
    free(results);
    free(angles);


    // This is C4 relative to C5
    fprintf(dest,"\n\nC4 to C5, X,Y,Z\n");

    // This finds the next frame with an absolute difference > 4 degrees
    angles = malloc(2*500*sizeof(int));
    k=0;
    for(i=2;i<frame;i++){
        for(j=i;j<frame;j++){
            if(flexion_all[4*i+1]+4 < flexion_all[4*j+1]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-4 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C4_file,C5_file,angles,k);
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
            if(flexion_all[4*i+2]+4 < flexion_all[4*j+2]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-4 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C5_file,C6_file,angles,k);
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
            if(flexion_all[4*i+3]+4 < flexion_all[4*j+3]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
            if(flexion_all[4*i]-4 > flexion_all[4*j]){
                angles[2*k]=i;
                angles[2*k+1]=j;
                k++;
                break;
            }
        }
    }

    results = cor(C6_file,C7_file,angles,k);
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
    int framecount;
    double bone1_rot[500][3][3];
    double bone2_rot[500][3][3];
    double bone2_invrot[500][3][3];
    double position1[500][3];
    double position2[500][3];
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
            if(k==4){
                k=0;
                j++;
            }
            if(j==3 && k!=3){
                position1[frame][k]=atof(rowFields[i]);
            }
            else if(k!=3){
                bone1_rot[frame][j][k]=atof(rowFields[i]);
            }
            k++;
        }
        frame++;
        CsvParser_destroy_row(row);
    }
    CsvParser_destroy(csvparser1);

    CsvParser *csvparser3 = CsvParser_new(file_2, ",", 0);
    frame=0;
    while ((row = CsvParser_getRow(csvparser3)) ) {
        const char **rowFields = CsvParser_getFields(row);
        j=k=0;
        for (i = 2 ; i < CsvParser_getNumFields(row) ; i++) {
            //Need to filter out every fourth column beginning with first.
            if(k==4){
                k=0;
                j++;
            }
            if(j==3 && k!=3){
                position2[frame][k]=atof(rowFields[i]);
            }
            else if(j<3 && k!=3){
                bone2_invrot[frame][j][k]=atof(rowFields[i]);
            }
            k++;
        }
        frame++;
        CsvParser_destroy_row(row);
    }
    CsvParser_destroy(csvparser3);

    for(framecount=1;framecount<frame;framecount++){
        double temp[3][3];
        CoFactor(bone2_invrot[framecount],3,temp);
        double det=matrix_det(bone2_invrot[framecount]);
        transpose(temp,bone2_rot[framecount]);
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                bone2_rot[framecount][i][j]/=det;
            }
        }
    }


    for(framecount=0;framecount<frame;framecount++){
        mult_matrix(bone2_rot[framecount],bone1_rot[framecount],bone12_rot[framecount]);
    }


    for(framecount=0;framecount<frame;framecount++){
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                printf("%f\t",bone12_rot[framecount][i][j]);
            }
            printf("\n");
        }
        printf("\n%d\n",framecount);
    }

    for(framecount=0;framecount<frame;framecount++){
        for(j=0;j<3;j++){
            position12[framecount][j]=position1[framecount][j]-position2[framecount][j];
        }
    }


    // At this point all data is loaded into memory in the appropriate format.
    for(i=0;i<count;i++){
        //printf("\n\n%d\t%d",angles[2*i+1]-1,angles[2*i]-1);
        icr(type,bone12_rot[(angles[2*i+1]-1)],bone12_rot[(angles[2*i]-1)],position12[(angles[2*i+1]-1)], position12[(angles[2*i]-1)],&phi,&screw_vec,&t,intersection);

        for(j=0;j<3;j++){
            results[i][j]=intersection[j];
        }
    }
    return results;
}


/*
 *    Find the cofactor matrix of a square matrix
 */

void CoFactor(a,n,b)
    double a[3][3];
    int n;
    double b[3][3];
{
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;

    c = malloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
        c[i] = malloc((n-1)*sizeof(double));

    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {

            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }

            /* Calculate the determinate */
            det = Determinant(c,n-1);

            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        free(c[i]);
    free(c);
}



/*
 *    Recursive definition of determinate using expansion by minors.
 */
double Determinant(double **a,int n)
{
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */

    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = malloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
            for (i=0;i<n-1;i++)
                free(m[i]);
            free(m);
        }
    }
    return(det);
}
