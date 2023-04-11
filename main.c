#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TMAX 240
#define MASSE1 264000000.0
#define RAIDEUR_K1 225000000.0
#define MASSE2 660000.0
#define RAIDEUR_K2 510000.0
#define CONSTANTE_C 52000.0

void def_euleur(double pasT, double position1[], double vitesse1[], double position2[], double vitesse2[],int NMAX){

    FILE* fichier_vitesse = fopen("TAPEI_200000euleur_vitesseTx025.txt", "a"); //Fichier ecriture des positions
    FILE* fichier_position = fopen("TAPEI_200000euleur_positionTx025.txt", "a"); //Fichier ecriture des positions

    double K1_x1 = 0.0;//UN SEUL K PAR POSITION
    double K2_x2 = 0.0;
    double K1_v1 = 0.0;
    double K2_v2 = 0.0;

    position1[0] = 0.25;
    position2[0] = 0.0;
    vitesse1[0] = 0.0; 
    vitesse2[0] = 0.0;

    double vec_U1x = position1[0];
    double vec_U1v = 0.0;
    double vec_U2x = 0.0;
    double vec_U2v = 0.0;

    for(int i = 0; i<NMAX; i++){

        K1_x1 = vitesse1[i];
        K2_x2 = vitesse2[i];
        K1_v1 = -((RAIDEUR_K1+RAIDEUR_K2)/MASSE1)*position1[i]+(RAIDEUR_K2/MASSE1)*position2[i];
        K2_v2 = (RAIDEUR_K2/MASSE2)*position1[i]-(RAIDEUR_K2/MASSE2)*position2[i]+(CONSTANTE_C/MASSE2)*vitesse1[i]-(CONSTANTE_C/MASSE2)*vitesse2[i];
        vec_U1x += pasT*K1_x1;
        vec_U2x += pasT*K2_x2;
        vec_U1v += pasT*K1_v1;
        vec_U2v += pasT*K2_v2;
        position1[i+1] = vec_U1x;
        position2[i+1] = vec_U2x;
        vitesse1[i+1] = vec_U1v;
        vitesse2[i+1] = vec_U2v;
        //printf("K1 : %f\n %f\n %f\n %f\n", K1c1, K1c2, K1c3, K1c4);
        //printf("vec U : %f\n %f\n %f\n %f\n", vec_U1, vec_U2, vec_U3, vec_U4);
        fprintf(fichier_vitesse, "%f %f %f\n", i*pasT, vitesse1[i], vitesse2[i]);
        fprintf(fichier_position, "%f %f %f\n", i*pasT, position1[i], position2[i]);
    }

}

void def_RK4(double pasT, double position1[], double vitesse1[], double position2[], double vitesse2[], int NMAX){

    FILE* fichier_vitesse = fopen("TAPEI_2000RK4_vitesseT.txt", "a"); //Fichier ecriture des positions
    FILE* fichier_position = fopen("TAPEI_2000RK4_positionT.txt", "a"); //Fichier ecriture des positions

    double K1_x1, K1_v1, K2_x1, K2_v1, K3_x1, K3_v1, K4_x1, K4_v1 = 0.0;//4 K PAR POSITION
    double K1_x2, K1_v2, K2_x2, K2_v2, K3_x2, K3_v2, K4_x2, K4_v2 = 0.0;

    position1[0] = 3;
    position2[0] = 0.0;
    vitesse1[0] = 0.0;
    vitesse2[0] = 0.0;

    double vec1_Uy, vec2_Ux, vec2_Uy = 0.0;
    double vec1_Ux = position1[0];

    for(int i = 0; i<NMAX; i++){

        K1_x1 = vitesse1[i];
        K1_x2 = vitesse2[i];

        K1_v1 = -((RAIDEUR_K1+RAIDEUR_K2)/MASSE1)*position1[i]+(RAIDEUR_K2/MASSE1)*position2[i];
        K1_v2 = (RAIDEUR_K2/MASSE2)*position1[i]-(RAIDEUR_K2/MASSE2)*position2[i]+(CONSTANTE_C/MASSE2)*vitesse1[i]-(CONSTANTE_C/MASSE2)*vitesse2[i];

        K2_x1 = vitesse1[i] + (pasT/2.0)*K1_x1;
        K2_x2 = vitesse2[i] + (pasT/2.0)*K1_x2;

        K2_v1 = -((RAIDEUR_K1+RAIDEUR_K2)/MASSE1)*(position1[i]+(pasT/2.0)*K1_x1)+(RAIDEUR_K2/MASSE1)*(position2[i]+(pasT/2.0)*K1_x1);
        K2_v2 = (RAIDEUR_K2/MASSE2)*(position1[i]+(pasT/2.0)*K1_x1)-(RAIDEUR_K2/MASSE2)*(position2[i]+(pasT/2.0)*K1_x2)+(CONSTANTE_C/MASSE2)*(vitesse1[i]+(pasT/2.0)*K1_v1)-(CONSTANTE_C/MASSE2)*(vitesse2[i]+(pasT/2.0)*K1_v2);

        K3_x1 = vitesse1[i] + (pasT/2.0)*K2_x1;
        K3_x2 = vitesse2[i] + (pasT/2.0)*K2_x2;

        K3_v1 = -((RAIDEUR_K1+RAIDEUR_K2)/MASSE1)*(position1[i]+(pasT/2.0)*K2_x1)+(RAIDEUR_K2/MASSE1)*(position2[i]+(pasT/2.0)*K2_x1);
        K3_v2 = (RAIDEUR_K2/MASSE2)*(position1[i]+(pasT/2.0)*K2_x1)-(RAIDEUR_K2/MASSE2)*(position2[i]+(pasT/2.0)*K2_x2)+(CONSTANTE_C/MASSE2)*(vitesse1[i]+(pasT/2.0)*K2_v1)-(CONSTANTE_C/MASSE2)*(vitesse2[i]+(pasT/2.0)*K2_v2);

        K4_x1 = vitesse1[i] + pasT*K3_x1;
        K4_x2 = vitesse2[i] + pasT*K3_x2;

        K4_v1 = -((RAIDEUR_K1+RAIDEUR_K2)/MASSE1)*(position1[i]+pasT*K3_x1)+(RAIDEUR_K2/MASSE1)*(position2[i]+pasT*K3_x1);;
        K4_v2 = (RAIDEUR_K2/MASSE2)*(position1[i]+pasT*K3_x1)-(RAIDEUR_K2/MASSE2)*(position2[i]+pasT*K3_x2)+(CONSTANTE_C/MASSE2)*(vitesse1[i]+pasT*K3_v1)-(CONSTANTE_C/MASSE2)*(vitesse2[i]+pasT*K3_v2);

        vec1_Ux = vec1_Ux + (pasT/6)*(K1_x1 + 2*K2_x1 + 2*K3_x1 + K4_x1);
        vec1_Uy = vec1_Uy + (pasT/6)*(K1_v1 + 2*K2_v1 + 2*K3_v1 + K4_v1);

        vec2_Ux = vec2_Ux + (pasT/6)*(K1_x2 + 2*K2_x2 + 2*K3_x2 + K4_x2);
        vec2_Uy = vec2_Uy + (pasT/6)*(K1_v2 + 2*K2_v2 + 2*K3_v2 + K4_v2);

        position1[i+1] = vec1_Ux;
        vitesse1[i+1] = vec1_Uy;
        position2[i+1] = vec2_Ux;
        vitesse2[i+1] = vec2_Uy;

        fprintf(fichier_vitesse, "%f %f %f\n", i*pasT, vitesse1[i], vitesse2[i]);
        fprintf(fichier_position, "%f %f %f\n", i*pasT, position1[i], position2[i]);

    }
}

int main(){ 

    FILE* fichier_vitesse = fopen("TAPEI_200000euleur_vitesseTx025.txt", "w+"); //Fichier ecriture des positions
    FILE* fichier_position = fopen("TAPEI_200000euleur_positionTx025.txt", "w+"); //Fichier ecriture des positions
    fclose(fichier_position);
    fclose(fichier_vitesse);
    FILE* fichier_vitesseRK4 = fopen("TAPEI_2000RK4_vitesseT.txt", "w+"); //Fichier ecriture des positions
    FILE* fichier_positionRK4 = fopen("TAPEI_2000RK4_positionT.txt", "w+"); //Fichier ecriture des positions
    fclose(fichier_positionRK4);
    fclose(fichier_vitesseRK4);

    //int NMAX[3] = {2000, 20000, 200000};
    int NMAX[1] = {2000};
    /*DEFINITION du pas*/
    for(int w = 0; w<1; w++){
        double pasT = (float)TMAX/(float)NMAX[w];
        printf("pasT%d : %f\n", NMAX[w], pasT);

    /*RESOLUTION EULEUR POSITION ET VITESSE DEBUT*/
            /*DEBUT creation du tableau x(t)*/
            // double *euleur_position1 = NULL;
            // euleur_position1 = malloc((NMAX[w]+1)*sizeof(double));

            // double *euleur_position2 = NULL;
            // euleur_position2 = malloc((NMAX[w]+1)*sizeof(double));
            // /*FIN creation du tableau x(t)*/

            // /*DEBUT creation du tableau v(t)*/
            // double *euleur_vitesse1 = NULL;
            // euleur_vitesse1 = malloc((NMAX[w]+1)*sizeof(double));
            // double *euleur_vitesse2 = NULL;
            // euleur_vitesse2 = malloc((NMAX[w]+1)*sizeof(double)); 
            // /*FIN creation du tableau v(t)*/

            // def_euleur(pasT, euleur_position1, euleur_vitesse1, euleur_position2, euleur_vitesse2, NMAX[w]);
    /*RESOLUTION EULEUR POSITION ET VITESSE DEBUT*/

    /*RESOLUTION RK4 POSITION ET VITESSE DEBUT*/
            double *RK4_position1 = NULL;
            RK4_position1 = malloc((NMAX[w]+1)*sizeof(double));

            double *RK4_position2 = NULL;
            RK4_position2 = malloc((NMAX[w]+1)*sizeof(double));
            /*FIN creation du tableau x(t)*/

            /*DEBUT creation du tableau v(t)*/
            double *RK4_vitesse1 = NULL;
            RK4_vitesse1 = malloc((NMAX[w]+1)*sizeof(double));
            double *RK4_vitesse2 = NULL;
            RK4_vitesse2 = malloc((NMAX[w]+1)*sizeof(double)); 
            /*FIN creation du tableau v(t)*/

            def_RK4(pasT, RK4_position1, RK4_vitesse1, RK4_position2, RK4_vitesse2, NMAX[w]);
    /*RESOLUTION RK4 POSITION ET VITESSE FIN*/

    }
    return 0;
}
