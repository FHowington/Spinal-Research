#include "include.h"
#include "csvparser.h"

int main() {
///REDO
//FOR C4 relative to C3 do: Inverse(C4) cross C3




    char C3_file[]="HomoTransMatrices_C123_-to-Lab.csv"; 
    char C4_file[]="HomoTransMatrices_C4_-to-Lab.csv";
    char C5_file[]="HomoTransMatrices_C5_-to-Lab.csv";
    char C6_file[]="HomoTransMatrices_C6_-to-Lab.csv";
    char C7_file[]="HomoTransMatrices_C7T1_-to-Lab.csv";
    double** cor(char *, char *);
    double* flexion_all = malloc(4*500*sizeof(double*));
    double *** results_all = malloc(4*sizeof(double**));
    FILE *dest;
    CsvRow *row;

    int rownum=5;
    
    dest = fopen("results.csv","wb");
    

    CsvParser *csvparser_angles = CsvParser_new("KinematicsMeasurementReport.csv", ",", 0);

    int i=0;
    while ((row = CsvParser_getRow(csvparser_angles)) ) {
        const char **rowFields = CsvParser_getFields(row);
        flexion_all[4*ii]=atof(rowFields[2]);
        flexion_all[4*ii+1]=atof(rowFields[5]);
        flexion_all[4*ii+2]=atof(rowFields[8]);
        flexion_all[4*ii+3]=atof(rowFields[11]);
        CsvParser_destroy_row(row);
    }
    CsvParser_destroy(csvparser_angles);

    double ** results = cor(C3_file,C4_file);
    results_all[0]=results;
   
    results = cor(C4_file,C5_file);
    results_all[1]=results;
    
    results = cor(C5_file,C6_file);
    results_all[2]=results;

    results = cor(C6_file,C7_file);
    results_all[3]=results;

    for(int i=0;i<500;i++){
        if(results[i][1]==0.0 && results[i][2]==0.0)
            break;

        fprintf(dest, "%f,%f,%f, ,%f,%f,%f, ,%f,%f,%f, ,%f,%f,%f\n", results_all[0][i][0],results_all[0][i][1],results_all[0][i][2]) 
    } 

    return 0;
}


double** cor(char *file_1, char *file_2){
    int i = 0;
    int j = 0;
    int k = 0;
    int frame = 0;
    int type = 1;
    double bone1_CS[500][3][3];
    double position1[500][3];
    double bone2_CS[500][3][3];
    double position2[500][3];
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
            if(k!=3 && j<3){
                bone1_CS[frame][j][k]=atof(rowFields[i]);
            }
            if(j==3 && k<3){
                position1[frame][k]=atof(rowFields[i]);
            }
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

            if(k!=3 && j<3){
                bone2_CS[frame][j][k]=atof(rowFields[i]);
            }
            if(j==3 && k<3){
                position2[frame][k]=atof(rowFields[i]);
            }
            k++;
        }
        frame++;
        CsvParser_destroy_row(row);
    }

    CsvParser_destroy(csvparser2);

    //At this point all data is loaded into memory in the appropriate format.
    //Let's output the rotation angles of bone1 relative to bone2
    for(i=2;i<frame;i++){
        icr(type,bone1_CS[i],bone2_CS[i],position1[i],position2[i],&phi,&screw_vec,&t,intersection);
        
        for(j=0;j<3;j++)
            results[i-2][j]=intersection[j];
    }

    return results;
}
