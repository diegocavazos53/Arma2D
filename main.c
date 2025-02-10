#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Programa que resuelve armaduras planas.

//#define nj 6
//#define nm 10
//#define ns 3
#define ncjt 2
//#define nmp 2
//#define ncp 3
//#define njl 3

#define tam 4
double gk[tam][tam];
double s[100][100];
double p[100];
double v[tam];
double t[tam][tam];
double u[tam];
double bk[tam][tam];
double q[tam];
double f[tam];
double r[100];
double maf[100];

void mstiffg(double modelast, double area, double bl, double cx, double cy){
    double z = modelast * area / bl;
    double z1 = z * (pow(cx,2));
    double z2 = z * (pow(cy,2));
    double z3 = z * cx * cy;
    gk[0][0] = z1;
    gk[1][0] = z3;
    gk[2][0] = -z1;
    gk[3][0] = -z3;
    gk[0][1] = z3;
    gk[1][1] = z2;
    gk[2][1] = -z3;
    gk[3][1] = -z2;
    gk[0][2] = -z1;
    gk[1][2] = -z3;
    gk[2][2] = z1;
    gk[3][2] = z3;
    gk[0][3] = -z3;
    gk[1][3] = -z2;
    gk[2][3] = z3;
    gk[3][3] = z2;
}

void stores(int jb, int je, int ndof, int nsc[]){
    int i1 = 0;
    int n1 = 0;
    int n2 = 0;
    for (int i = 0; i < 2 * ncjt; i++){
        if (i < ncjt){
            i1 = (jb - 1) * ncjt + i;
        }
        else{
            i1 = (je - 1) * ncjt + (i - ncjt);
        }
        n1 = nsc[i1];
        if(n1 <= ndof){
            for(int j = 0; j < 2 * ncjt; j++){
                if (j < ncjt){
                    i1 = (jb - 1) * ncjt + j;
                }
                else{
                    i1 = (je - 1) * ncjt + (j - ncjt);
                }
                n2 = nsc[i1];
                if(n2 <= ndof){
                    s[n1 - 1][n2 - 1] = s[n1-1][n2-1]+gk[i][j];
                }
            }
        }
    }
}

void mdispg(int jb, int je, int ndof, int nsc[]){
    int j = 0;
    int n = 0;
    for(int i = 0; i < 2 * ncjt; i++){
        v[i] = 0;
    }
    j = (jb - 1) * ncjt;
    for(int i = 0; i < ncjt; i++){
        j = j + 1;
        n = nsc[j - 1];
        if(n <= ndof){
            v[i] = p[n-1];
        }
    }
    j = (je - 1) * ncjt;
    for(int i = ncjt + 1; i <= 2 * ncjt; i++){
        j = j + 1;
        n = nsc[j - 1];
        if(n <= ndof){
            v[i - 1] = p[n - 1];
        }
    }
}

void mtrans(double cx, double cy){
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            t[i][j] = 0;
        }
    }
    t[0][0] = cx; t[1][0] = -cy;
    t[0][1] = cy; t[1][1] = cx;
    t[2][2] = cx; t[3][2] = -cy;
    t[2][3] = cy; t[3][3] = cx;
}

void mdispl(){
    for (int i = 0; i < 2 * ncjt; i++){
        u[i] = 0;
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            u[i] = u[i] + t[i][j] * v[j];
        }
    }
}

void mstiffl(double modelast, double area, double bl){
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            bk[i][j] = 0;
        }
    }
    double z = modelast * area / bl;
    bk[0][0] = z; bk[2][0] = -z;
    bk[0][2] = -z; bk[2][2] = z;
}

void mforcel(){
    for (int i = 0; i < 2 * ncjt; i++){
        q[i] = 0;
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            q[i] = q[i] + bk[i][j] * u[j];
        }
    }
}

void mforceg(){
    for (int i = 0; i < 2 * ncjt; i++){
        f[i] = 0;
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            f[i] = f[i] + t[j][i] * q[j];
        }
    }
}

void storer(int jb, int je, int ndof, int nsc[]){
    int i1 = 0;
    int n = 0;
    for (int i = 1; i <= 2 * ncjt; i++){
        if(i <= ncjt){
            i1 = (jb - 1) * ncjt + i;
        }
        else{
            i1 = (je - 1) * ncjt + (i - ncjt);
        }
        n = nsc[i1 - 1];
        if(n > ndof){
            r[n - ndof - 1] = r[n - ndof - 1] + f[i - 1];
        }
    }
}

int main(int argc, char** argv)
{
    printf("Programa para resolver armaduras planas\n\n");
    printf("Argumentos: %d\n",argc); // Número de argumentos en la consola.
    // El argumento [0] es la dirección del ejecutable.
    // El argumento [1] lo utilizaremos para indicar un archivo a leer.
    printf("Argumentos: %s\n",argv[1]);
    char nomArchivo[15];
    if (argc > 1){
        strcpy(nomArchivo, argv[1]);
    }
    else{
        printf("¿Qué archivo desea leer?\n");
        scanf("%s", nomArchivo);
    }

    //char nomArchivo[15] = "ar1.a2d";
    printf("Leyendo: %s\n\n",nomArchivo);

    FILE *pFile;
    pFile = fopen(nomArchivo,"r");
    if(!pFile){
        printf("ERROR: No se pudo leer el archivo.\n");
        return 1;
    }

    int nj = 0;
    fscanf(pFile,"%d", &nj);

    double coord[nj][2];
    for (int i = 0; i < nj; i++){
        fscanf(pFile,"%lf",&coord[i][0]);
        fscanf(pFile,"%lf",&coord[i][1]);
    }

    int ns = 0;
    fscanf(pFile,"%d", &ns);

    int msup[ns][ncjt + 1];
    for (int i = 0; i < ns; i++){
        for(int j = 0; j < ncjt + 1; j++){
            fscanf(pFile,"%d",&msup[i][j]);
        }
    }

    int nmp = 0;
    fscanf(pFile,"%d", &nmp);

    double em[nmp];
    for (int i = 0; i < nmp; i++){
        fscanf(pFile,"%lf",&em[i]);
    }

    int ncp = 0;
    fscanf(pFile,"%d", &ncp);

    double cp[ncp];
    for (int i = 0; i < ncp; i++){
        fscanf(pFile,"%lf",&cp[i]);
    }

    int nm = 0;
    fscanf(pFile,"%d", &nm);

    int mprp[nm][4];
    for (int i = 0; i < nm; i++){
        for(int j = 0; j < 4; j++){
            fscanf(pFile,"%d",&mprp[i][j]);
        }
    }

    int njl = 0;
    fscanf(pFile,"%d", &njl);

    int jp[njl];
    double pj[njl][ncjt];
    for(int i = 0; i < njl; i++){
        fscanf(pFile,"%d", &jp[i]);
        for(int j = 0; j < ncjt; j++){
            fscanf(pFile,"%lf",&pj[i][j]);
        }
    }

    printf("El archivo se leyó correctamente.\n");

    if(fclose(pFile) != 0){
        printf("Error: No se pudo cerrar el archivo.\n");
    }
    printf("\n");

    //double coord[6][2] = {{0,0},{288,0},{576,0},{864,0},{288,216},{576,216}};
    //int msup[ns][ncjt + 1] = {{1,1,1},{3,0,1},{4,0,1}};
    //double em[nmp] = {29000,10000};
    //double cp[ncp] = {8,12,16};
    //int mprp[nm][4] = {{1,2,1,1},{2,3,1,1},{3,4,2,3},{5,6,1,1},{2,5,1,1},{3,6,1,1},{1,5,1,2},{2,6,1,2},{3,5,1,2},{4,6,2,3}};
    //int jp[njl] = {2,5,6};
    //double pj[njl][ncjt] = {{0,-75},{25,0},{0,-60}};

    printf("     **********************************\n");
    printf("     *      Programa de Armaduras     *\n");
    printf("     *              por               *\n");
    printf("     *      Diego Cavazos de Lira     *\n");
    printf("     **********************************\n\n");
    printf("        =============================\n");
    printf("        DATOS ESTRUCTURALES GENERALES\n");
    printf("        =============================\n\n");

    printf("Tipo de estructura: Armadura Plana\n");
    printf("Número de nodos: %d\n", nj);
    printf("Número de miembros: %d\n", nm);
    printf("Número de conjuntos de propiedades de materiales: %d\n", nmp);
    printf("Número de conjuntos de propiedades de sección transversal: %d\n\n", ncp);

    printf("        =============================\n");
    printf("          Coordenadas de los nodos\n");
    printf("        =============================\n");
    printf("No. Nodo     Coordenada X         Coordenada Y\n");
    printf("--------     ------------         ------------\n");
    for(int i = 0; i < nj; i++){
        printf("%4d %20.4lf %20.4lf\n", i+1, coord[i][0], coord[i][1]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("                   Apoyos\n");
    printf("        =============================\n");
    printf("No. Nodo      Restricción X    Restricción Y\n");
    printf("--------      -------------    -------------\n");
    for(int i = 0; i < ns; i++){
        printf("%4d %16d %16d\n", msup[i][0], msup[i][1], msup[i][2]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("        Propiedades de los materiales\n");
    printf("        =============================\n");
    printf("Número        Módulo de\n");
    printf("material     Elasticidad\n");
    printf("--------     -----------\n");
    for(int i = 0; i < nmp; i++){
        printf("%4d %19.4lf\n", i+1, em[i]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("        Propiedades de secc. transv.\n");
    printf("        =============================\n");
    printf("No. Prop.      Area (A)\n");
    printf("---------     ---------\n");
    for(int i = 0; i < ncp; i++){
        printf("%4d %18.4lf\n", i+1, cp[i]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("         Información de los miembros\n");
    printf("        =============================\n");
    printf("Miembro  Nodo de  Nodo de  Número    Número Prop.\n");
    printf("número   Inicio      Fin   material  secc. transv.\n");
    printf("-------  -------  -------  --------  -------------\n");
    for(int i = 0; i < nm; i++){
        printf("%4d %8d %8d %9d %11d\n", i+1, mprp[i][0], mprp[i][1], mprp[i][2], mprp[i][3]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("                Cargas nodales\n");
    printf("        =============================\n");
    printf("No. Nodo       Fuerza X           Fuerza Y\n");
    printf("--------       --------           --------\n");
    for(int i = 0; i <njl; i++){
        printf("%4d %18.4lf %18.4lf\n", jp[i], pj[i][0], pj[i][1]);
    }
    printf("\n");
    printf("****** Fin de la información de entrada ******\n\n");


    // Parte 7
    int nr = 0;
    int ndof = 0;
    for (int i = 0; i < ns; i++){
        for (int j = 1; j < ncjt + 1; j++){
            if (msup[i][j] == 1){
                nr = nr + 1;
            }
        }
    }
    ndof = ncjt * nj - nr;
    printf("El número de coordenadas restringidas NR es: %d\n", nr);
    printf("El número de grados de libertad NDOF es: %d\n", ndof);

    // Parte 8
    int nsc[ncjt*nj];
    int icount = 0;
    int i3 = 0;
    int j = 0;
    int k = ndof;
    for(int i = 1; i <= nj; i++){
        icount = 0;
        for(int i1 = 0; i1 < ns; i1++){
            if(msup[i1][0] == i){
                icount = 1;
                for(int i2 = 0; i2 < ncjt; i2++){
                    i3 = (i - 1) * ncjt + i2;
                    if(msup[i1][i2+1] == 1){
                        k = k + 1;
                        nsc[i3] = k;
                    }
                    else{
                        j = j + 1;
                        nsc[i3] = j;
                    }
                }
            }
        }
        if (icount == 0){
            for(int i2 = 0; i2 < ncjt; i2++){
                i3 = (i-1) * ncjt + i2;
                j = j + 1;
                nsc[i3] = j;
            }
        }
    }
    printf("El vector de números de coordenadas de la estructura NSC, es:\n");
    for(int i = 0; i < ncjt * nj; i++){
        printf("%4d\n", nsc[i]);
    }

    // Parte 9
    int jb = 0; int je = 0;
    int imatnum = 0; int isecctransvnum = 0;
    double modelast = 0;
    double area = 0;
    double xb = 0; double xe = 0;
    double yb = 0; double ye = 0;
    double bl = 0;
    double cx = 0; double cy = 0;

    for (int im = 0; im < nm; im++){
        jb = mprp[im][0];
        je = mprp[im][1];
        imatnum = mprp[im][2];
        isecctransvnum = mprp[im][3];
        modelast = em[imatnum - 1];
        area = cp[isecctransvnum - 1];
        xb = coord[jb - 1][0];
        yb = coord[jb - 1][1];
        xe = coord[je - 1][0];
        ye = coord[je - 1][1];
        bl = sqrt(pow(xe-xb,2)+pow(ye-yb,2));
        cx = (xe-xb)/bl;
        cy = (ye-yb)/bl;

        mstiffg(modelast, area, bl, cx, cy);
        printf("La matriz de rigidez en el SCG de la barra %d es:\n", im + 1);
        for(int i = 0; i < 2*ncjt; i++){
            for(int j = 0; j < 2* ncjt; j++){
                printf("%15.3f ", gk[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        stores(jb, je, ndof, nsc);
    }
    printf("La matriz de rigidez del sistema es:\n");
    for(int i = 0; i < ndof; i++){
        for(int j = 0; j < ndof; j++){
            printf("%15.3f ", s[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // Parte 10
    for (int i = 0; i < njl; i++){
        int i1 = jp[i];
        int i2 = (i1-1) * ncjt;
        int n = 0;
        for (int j = 0; j < ncjt; j++){
            i2 = i2 + 1;
            n = nsc[i2 - 1];
            if (n <= ndof){ // Se agregó un =
                p[n - 1] = p[n - 1] + pj[i][j];
            }
        }
    }

    printf("El vector de fuerzas P, es:\n");
    for (int i = 0; i < ndof; i++){
        printf("%15.3lf\n", p[i]);
    }
    printf("\n");

    // Parte 11
    double z1 = 0;
    double z = 0;
    for(int i = 0; i < ndof; i++){
        z1 = s[i][i];
        for(int j = 0; j < ndof; j++){
            s[i][j] = s[i][j] / z1;
        }
        p[i] = p[i] / z1;
        for(int k = 0; k < ndof; k++){
            if(k == i){
            }
            else{
                z = s[k][i];
                for(int m = i; m < ndof; m++){
                    s[k][m] = s[k][m] - s[i][m] * z;
                }
                p[k] = p[k] - p[i] * z;
            }
        }
    }

    printf("El vector de desplazamientos d, es:\n");
    for (int i = 0; i < ndof; i++){
        printf("%15.3lf\n", p[i]);
    }
    printf("\n");

    // Parte 12
    for (int im = 0; im < nm; im++){
        jb = mprp[im][0];
        je = mprp[im][1];
        imatnum = mprp[im][2];
        isecctransvnum = mprp[im][3];
        modelast = em[imatnum - 1];
        area = cp[isecctransvnum - 1];
        xb = coord[jb - 1][0];
        yb = coord[jb - 1][1];
        xe = coord[je - 1][0];
        ye = coord[je - 1][1];
        bl = sqrt(pow(xe-xb,2)+pow(ye-yb,2));
        cx = (xe-xb)/bl;
        cy = (ye-yb)/bl;

        mdispg(jb, je, ndof, nsc);
        printf("El vector de desplazamientos globales V de la barra no. %d, es:\n",im + 1);
        for (int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", v[i]);
        }
        printf("\n");

        mtrans(cx, cy);
        printf("La matriz de transformación de la barra %d es:\n",im+1);
        for (int i = 0; i < 2 * ncjt; i++){
            for (int j = 0; j < 2 * ncjt; j++){
                printf("%15.3lf", t[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        mdispl();
        printf("El vector de desplazamientos en SCL, u de la barra %d es:\n", im+1);
        for (int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n",u[i]);
        }
        printf("\n");

        mstiffl(modelast, area, bl);
        printf("La matriz de rigidez bk de la barra no. %d es:\n",im+1);
        for (int i = 0; i < 2 * ncjt; i++){
            for(int j = 0; j < 2 * ncjt; j++){
                printf("%15.3lf", bk[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        mforcel();
        printf("El vector de fuerzas en el SCL Q de la barra no. %d es:\n", im+1);
        for(int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", q[i]);
        }
        printf("\n");

        mforceg();
        printf("El vector de fuerzas en el SCG F de la barra no. %d es:\n", im+1);
        for(int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", f[i]);
        }
        printf("\n");

        storer(jb, je, ndof, nsc);
        maf[im] = q[0];
    }

    printf("El vector de reacciones en los apoyos R es:\n");
    for(int i = 0; i < nr; i++){
        printf("%15.3lf\n", r[i]);
    }
    printf("\n");

// Impresión de salida
    printf("***************************************************\n");
    printf("*             Resultados del Análisis             *\n");
    printf("***************************************************\n\n");

    printf("============================\n");
    printf("Desplazamientos de los nodos\n");
    printf("============================\n\n");

    printf("No. Nodo     Traslación X    Traslación Y\n");
    printf("========     ============    ============\n");
    double vecdespl[nj*2];
    int contador = 0;
    for(int i = 0; i < nj * 2 + 1; i++){
        if(nsc[i] > ndof){
            vecdespl[i] = 0;
        }
        if(nsc[i] <= ndof){
            vecdespl[i] = p[contador];
            contador = contador + 1;
        }
    }

    int nnodo = 1;
    for(int i = 0; i < nj * 2; i++){
        printf("%4d %20.4lf %15.4lf\n",nnodo,vecdespl[i],vecdespl[i+1]);
        i = i + 1;
        nnodo = nnodo + 1;
    }
    printf("\n");

    printf("===============================\n");
    printf("Fuerzas Axiales en los Miembros\n");
    printf("===============================\n\n");

    printf("Miembro       Fuerza Axial (Qa)\n");
    printf("=======       =================\n");

    for (int i = 0; i < nm; i++){
        char tenscompr;
        if(maf[i] > 0){
            tenscompr = 'C';
        }
        else if(maf[i] < 0){
            tenscompr = 'T';
        }
        else{
            tenscompr = '0';
        }
        printf("%4d %19.3lf (%c)\n", i+1, fabs(maf[i]),tenscompr);
    }
    printf("\n");

    printf("==============================\n");
    printf("   Reacciones en los apoyos   \n");
    printf("==============================\n\n");

    printf("No. Nodo     Fuerza X       Fuerza Y\n");
    printf("========     ========       ========\n");
    double vecreacc[ns*2];
    int contadorreacc = 0;
    int contadorceros = 0;
    for (int i = 0; i < ns; i++){
        if(msup[i][1] == 1){
            vecreacc[contadorreacc + contadorceros] = r[contadorreacc];
            contadorreacc = contadorreacc + 1;
        }
        else{
            vecreacc[contadorreacc + contadorceros] = 0;
            contadorceros = contadorceros + 1;
        }
        if(msup[i][2] == 1){
            vecreacc[contadorreacc + contadorceros] = r[contadorreacc];
            contadorreacc = contadorreacc + 1;
        }
        else{
            vecreacc[contadorreacc + contadorceros] = 0;
            contadorceros = contadorceros + 1;
        }
    }

    int indnodo = 0;
    for(int i = 0; i < ns * 2; i++){
        printf("%4d %15.3lf %15.3lf\n",msup[indnodo][0],vecreacc[i],vecreacc[i+1]);
        i = i + 1;
        indnodo = indnodo + 1;
    }
    printf("\n");

    printf("********** Fin del Análisis **********\n");
    getchar();

    return 0;
}
