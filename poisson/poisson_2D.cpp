//
// Método de diferencias finitas para resolver
// Ecuaciones en del tipo
// aUxx + bUyy + cUzz + dUx + eUy + fUz + hU = g
// Todos los coeficientes son funciones de (x,y,z) y g(x,y,z)
//

#include <stdio.h>
#include <math.h>
#include <gmm/gmm.h>

#define pi 3.14159
#define Interno 0
#define Dirichlet 1
#define Neumann 2

// La función tiene solución analítica
// 0 no tiene solución analítica
// 1 tiene solución analítica
#define SOL 0

// La ecuación es simétrica si
// los coeficientes d, e, f y sus respectivas derivadas
// no aparecen en la ecuación
// SIMETRICA 0 si la matriz resultante no será simétrica
// SIMETRICA 1 si la matriz resultante será simétrica
#define SIMETRICA 1

// Dimensión del problema
#define NDIM 2


// Ejemplo 1 
// Dos dimensiones
// Para n = 4
// -Uxx - Uyy = 2*4^2 pi*pi sen(4 pi x)sen(4 pi y)
// [-1, -1] x [1, 1]
// u(0,y) = 0
// u(x,0) = 0
// u(x,0.5) = 0
// u(0.5,y) = 0

/**
 * Coeficientes de la ecuación diferencial parcial
 * aUxx + bUyy + cUzz + dUx + eUy + fUz +hU = g
 * Todos los coeficientes son funciones de (x,y,z) y g(x,y,z)
 */
double a (double x, double y, double z) {
    return -1;
}

double b (double x, double y, double z) {

    return -1;
}
double c (double x, double y, double z) {

    return 0;
}
double d (double x, double y, double z) {

    return 0;
}
double e (double x, double y, double z) {

    return 0;
}
double f (double x, double y, double z) {

    return 0;
}
double h (double x, double y, double z) {

    return 0;
}

double g(double x, double y, double z) {
    double val = 32*pi*pi*sin(4*pi*x)*sin(4*pi*y);
    return val;
}

/**
 * Solución analítica de la EDP
 * en caso de tener, para poder comparar resultados
 */
double sol_analitica(double x, double y, double z) {
    double val = 0; 

    return val;
}

/**
 * Valores de las fronteras
 */
double valor_frontera(double x, double y, double z, // Coordenadas
                      double a_, double b_,
                      double c_, double d_,
                      double e_, double f_) {

    double val = 0;

    switch (NDIM) {

      case 1:

        if (x == a_) { //alpha_1
            printf("en alpha 1\n");
            val = 0;
        }

        if (x == b_) { //alpha_2
            printf("en alpha 2\n");
            val = 0;
        }

        printf("Valor de la frontera %f\n", val);

        return val;
      break;

      case 2:
        if (x == a_) { //alpha_4
            printf("en alpha 4\n");
            val = 0;
        }

        if (x == b_) { //alpha_2
            printf("en alpha 2\n");
            val = 0;
        }

        if (y == c_) { //alpha_1
            printf("en alpha 1\n");
            val = 0;
        }

        if (y == d_) { //alpha_3
            printf("en alpha 3\n");
            val = 0;
        }

        printf("Valor de la frontera %f\n", val);

        return val;
      break;

      case 3:
        if (y == c_) { //alpha_1
            printf("en alpha 1\n");
            val = 0;
        }

        if (x == b_) { //alpha_2
            printf("en alpha 2\n");
            val = 0;
        }

        if (y == d_) { //alpha_3
            printf("en alpha 3\n");
            val = 0;
        }

        if (x == a_) { //alpha_4
            printf("en alpha 4\n");
            val = 0;
        }
        
        if (z == e_) { //alpha_5
            printf("en alpha 5\n");
            val = 0;
        }

        if (z == f_) {//alpha_6
            printf("en alpha 6\n");
            val = 0;
        }
        printf("Valor de la frontera %f\n", val);

        return val;

      break;
    }

  return val;

}

/**
 * Imprime los resultados a archivo
 */
void print(gmm::dense_matrix<double> &xx, int nNodos, std::vector<double> &w, int nNodosInternos, int a_, int b_,int c_, int d_, int e_, int f_) {

    FILE *fp;
    fp = fopen("MDF.dat", "w");

    int j = 0;
    double val_front; // valor de la frontera
    for (int i = 0; i < nNodos; i++) {

        if (xx(i,3) != Interno)  {

            val_front = valor_frontera(xx(i,0), xx(i,1), xx(i,2), a_,b_,c_,d_,e_,f_);

            switch (NDIM) {
                case 1:
                    fprintf(fp, "%.2f %.3f\n", xx(i,0), val_front);
                    break;
                case 2:
                    fprintf(fp, "%.2f %.2f %.3f\n", xx(i,0), xx(i,1), val_front);
                    break;
                case 3:
                    fprintf(fp, "%.2f %.2f %.2f %.6f\n", xx(i,0), xx(i,1), xx(i,2), val_front);
                    break;

            } // fin switch

        } else {
            switch (NDIM) {
                case 1:
                    fprintf(fp, "%.2f %.3f\n", xx(i,0), w[j]);
                    break;
                case 2:
                    fprintf(fp, "%.2f %.2f %.3f\n", xx(i,0), xx(i,1), w[j]);
                    break;
                case 3:
                    fprintf(fp, "%.2f %.2f %.2f %.6f\n", xx(i,0), xx(i,1), xx(i,2), w[j]);
                    break;

            } // fin switch
            j++;
            } // else

    }

    fclose(fp);
}


/**
 * Imprime los resultados a archivo para ser graficados en matlab
 */
void print_matlab(gmm::dense_matrix<double> &xx, int nNodos, std::vector<double> &w, int nNodosInternos, int a_, int b_,int c_, int d_, int e_, int f_, int n1, int n2, std::vector<double> &x, std::vector<double> &y) {

    int i = 0;
    // Imprimimos el espacio
    FILE *fpx;
    fpx = fopen("EspacioMatlabX.dat", "w");
    for (i = 0; i < n1; i++) 
        fprintf(fpx, "%.3f\n", x[i]);
    fclose(fpx);
    FILE *fpy;
    fpy = fopen("EspacioMatlabY.dat", "w");
    for (i = 0; i < n1; i++) 
        fprintf(fpy, "%.3f\n", y[i]);
    fclose(fpy);

    FILE *fp;
    fp = fopen("Matlab.dat", "w");

    int j = 0;
    int aux = 0;
    double val_front; // valor de la frontera
    for (i = 0; i < nNodos; i++) {

        if (aux == n1) {
            fprintf(fp, "\n");
            aux = 0;
        }
        aux++;

        if (xx(i,3) != Interno)  {

            val_front = valor_frontera(xx(i,0), xx(i,1), xx(i,2), a_,b_,c_,d_,e_,f_);

                    fprintf(fp, "%.4f ", val_front);

        } else {
            fprintf(fp, "%.4f ", w[j]);
            j++;
        } // else

    }

    fclose(fp);
}


void creacionNodos(gmm::dense_matrix<double> &xx, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                   int nNodos, int n1, int n2, int n3, int &nNodosInternos,
                   double a_, double b_, double c_, double d_, double f_, double e_) {
    int k = 0;

    switch (NDIM) {

    case 1:
        for (int ii = 0; ii < n1; ii++) {
            xx(k,0) = x[ii];
            xx(k,1) = 0;
            xx(k,2) = 0;
            printf("xx(%d, x,y,z) = %f, %f, %f\n", k, xx(k,0), xx(k,1), xx(k,2));

            if (xx(k,0) == a_) { // alfa1
                printf("En la frontera alfa1, nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }
            if (xx(k,0) == b_) { // alfa2
                printf("En la frontera alfa2 , nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }

            if (xx(k,3) != 1) {
                printf("Nodo interno %d\n", k);
                nNodosInternos++;
            }
            k++;
        }
        break;

    case 2:
        for (int jj = 0; jj < n2; jj++) {
            for (int ii = 0; ii < n1; ii++) {
                xx(k,0) = x[ii];
                xx(k,1) = y[jj];
                xx(k,2) = 0;
                printf("xx(%d, x,y,z) = %f, %f, %f\n", k, xx(k,0), xx(k,1), xx(k,2));

            if (xx(k,0) == a_) { //alpha_1
                printf("En la frontera alfa1, nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }

            if (xx(k,0) == b_) { //alpha_2
                printf("En la frontera alfa1, nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }

            if (xx(k,1) == c_) { //alpha_3
                printf("En la frontera alfa1, nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }

            if (xx(k,1) == d_) { //alpha_4
                printf("En la frontera alfa1, nodo %d\n", k);
                xx(k,3) = 1;
                xx(k,4) = Dirichlet;
            }
            if (xx(k,3) != 1) {
                printf("Nodo interno %d\n", k);
                nNodosInternos++;
            }
            k++;

            }
        }
        break;

    case 3:
        for (int kk = 0; kk < n3; kk++) {
            for (int jj = 0; jj < n2; jj++) {
                for (int ii = 0; ii < n1; ii++) {
                    xx(k,0) = x[ii];
                    xx(k,1) = y[jj];
                    xx(k,2) = z[kk];
                    printf("xx(%d, x,y,z) = %f, %f, %f\n", k, xx(k,0), xx(k,1), xx(k,2));

                if (xx(k,1) == c_) { // alfa1
                    printf("frontera alfa1, nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,0) == b_) { // alfa2
                    printf("En la frontera alfa2 , nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,1) == d_) { // alfa3
                    printf("En la frontera alfa3, nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,0) == a_) { // alfa4
                    printf("En la frontera alfa4, nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,2) == e_) { // alfa5
                    printf("En la frontera alfa5, nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,2) == f_) { // alfa6
                    printf("En la frontera alfa6, nodo %d\n", k);
                    xx(k,3) = 1;
                    xx(k,4) = Dirichlet;
                }

                if (xx(k,3) != 1) {
                    printf("Nodo interno %d\n", k);
                    nNodosInternos++;
                }
                k++;
                }
            }
        }

        break;
    }
}


/**
 */
int main(int argc, const char *argv[]) {

    // Creación de una malla de 3D
    // [a,b] x [c,d] x [e,f]
    double a_ = -1;
    double b_ = 1;
    double c_ = -1;
    double d_ = 1;
    double e_ = 0;
    double f_ = 0;
    switch (NDIM) {
    case 1:
        printf("Dominio:\n[%.2f, %.2f]\n", a_, b_);
        break;
    case 2:
        printf("Dominio:\n[%.2f, %.2f] x [%.2f, %.2f]\n", a_, b_, c_, d_);
        break;
    case 3:
        printf("Dominio:\n[%.2f, %.2f] x [%.2f, %.2f] x [%.2f , %.2f]\n", a_, b_, c_, d_, e_, f_);
        break;
    default:
        fprintf(stderr, "NDIM debe ser 1, 2 o 3\n");
        exit(0);
    }

    // Numero de nodos para cada eje
    // n1 --> x
    // n2 --> y
    // n3 --> z
    int n1, n2, n3;
    switch(NDIM) {
    case 1:
        n1 = 0;
        n2 = 0;
        n3 = 0;
        break;
    case 2:
        n1 = 100;
        n2 = 100;
        n3 = 0;
        break;
    case 3:
        n1 = 10;
        n2 = 10;
        n3 = 10;
        break;
    }

    // Número total de nodos
    int nNodos;
    switch(NDIM) {
    case 1:
        if(n1 >= 4) nNodos = n1;
        else {
            fprintf(stderr, "n1 debe ser mayor o igual a 4\n");
            exit(0);
        }
        break;

    case 2:
        if (n1 >= 4 && n2 >= 4) nNodos = n1*n2;
        else {
            fprintf(stderr, "n1 y n2 deben ser mayores o iguales a 4\n");
            exit(0);
        }
        break;

    case 3:
        if (n1 >= 4 && n2 >= 4 && n3 >= 4) nNodos = n1*n2*n3;
        else {
            fprintf(stderr, "n1, n2 y n3 deben ser mayores o iguales a 4\n");
            exit(0);
        }
        break;
    }
    printf("nNodos %d\n", nNodos);

    // Tamaños de los intervalos
    double h1, h2, h3;
    switch (NDIM) {
    case 1:
        h1 =  (b_ - a_)/(n1-1);
        h2 = 1;
        h3 = 1;
        printf("h1 %f\n", h1);
        break;
    case 2:
        h1 =  (b_ - a_)/(n1-1);
        h2 =  (d_ - c_)/(n2-1);
        h3 = 1;
        printf("h1 %f\th2 %f\n", h1, h2);
        break;
    case 3:
        h1 =  (b_ - a_)/(n1-1);
        h2 =  (d_ - c_)/(n2-1);
        h3 =  (f_ - e_)/(n3-1);
        printf("h1 %f\th2 %f\th3 %f\n", h1, h2, h3);
        break;
    }


    // Vectores
    // donde se guardan la coordenada de la partición
    // por cada eje x,y,z
    std::vector<double> x(n1), y(n2), z(n3);

    switch(NDIM) {
    case 1:
        for (int i = 0; i < n1; i++) x[i] = a_ + i*h1;
        break;
    case 2:
        for (int i = 0; i < n1; i++) x[i] = a_ + i*h1;
        for (int j = 0; j < n2; j++) y[j] = c_ + j*h2;
        break;
    case 3:
        for (int i = 0; i < n1; i++) x[i] = a_ + i*h1;
        for (int j = 0; j < n2; j++) y[j] = c_ + j*h2;
        for (int k = 0; k < n3; k++) z[k] = e_ + k*h3;
        break;
    }

    // Matriz con la información de los nodos:
    // Entradas 0-2 son las coordenadas
    // Definimos las coordenadas de los nodos y si son forntera
    // Creación de los puntos de la malla
    // Cada nodo tiene tres coordenadas
    // {0,1} no es frontera o sí es respectivamente
    // {1,2} si la frontera es Dirichlet o Neumann
    gmm::dense_matrix<double> xx(nNodos, 5);
    int nNodosInternos = 0;
    creacionNodos(xx, x, y, z, nNodos, n1, n2, n3, nNodosInternos, a_, b_, c_, d_, f_, e_);

    printf("\nNodos Internos %d\n", nNodosInternos);


    // Empezamos la construcción de la matriz A y del vector zz
    // A*w = zz

    // Matriz dispersa
    gmm::row_matrix< gmm::rsvector<double> > A(nNodosInternos, nNodosInternos);
    std::vector<double> w(nNodosInternos), zz(nNodosInternos);

    // Vector para calcular la solución analítica
    // Contiene las coordenadas de los nodos internos
    gmm::dense_matrix<double> espacio(nNodosInternos, 3);

    int i = 0; // Contador de los nodos internos
    for (int ii = 0; ii < nNodos ; ii++) {

        if (xx(ii,3) == 0) {
            printf("\n\nNodo %d\n", ii);
            printf("Valor de i %d\n", i);

            // Nodo i
            A(i,i) = -2*(a(xx(ii,0), xx(ii,1), xx(ii,2))/(h1*h1) + b(xx(ii,0), xx(ii,1), xx(ii,2))/(h2*h2) + c(xx(ii,0), xx(ii,1), xx(ii,2))/(h3*h3) - h(xx(ii,0), xx(ii,1), xx(ii,2))/2);

            zz[i] = g(xx(ii,0), xx(ii,1), xx(ii,2));

            espacio(i,0) = xx(ii,0);
            espacio(i,1) = xx(ii,1);
            espacio(i,2) = xx(ii,2);

            // Nodo i-1
            switch((int)xx(ii-1,4)) {
            case Interno:
                printf("Nodo i-1 (%d) Nodo incognita\n", ii-1);
                A(i,i-1) = a(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(h1*h1) - d(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(2*h1);
                break;

            case Dirichlet:
                printf("Nodo i-1 (%d) Dirichlet\n", ii-1);
                zz[i] = zz[i] - valor_frontera(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2), a_,b_,c_,d_,e_,f_)*(a(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(h1*h1) - d(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(2*h1));
                break;

            case Neumann:
                break;

            default:
                fprintf(stderr, "Nodo desconocido %d\n", ii-1);
                exit(0);
            }

            // Nodo i+1
            switch((int)xx(ii+1,4)) {
            case Interno:
                printf("Nodo i+1 (%d) Nodo incognita\n", ii+1);
                A(i,i+1) = a(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(h1*h1) + d(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(2*h1);
                break;

            case Dirichlet:
                printf("Nodo i+1 (%d) Dirichlet\n", ii+1);
                zz[i] = zz[i] - valor_frontera(xx(ii+1,0), xx(ii+1,1), xx(ii+1,2), a_,b_,c_,d_,e_,f_)*(a(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(h1*h1) + d(xx(ii-1,0), xx(ii-1,1), xx(ii-1,2))/(2*h1));
                break;

            case Neumann:
                break;

            default:
                fprintf(stderr, "Nodo desconocido %d\n", ii+1);
                exit(0);
            }


            if (NDIM >= 2) {
                // Nodo i-n1
                switch((int) xx(ii-n1,4)) {
                case Interno:
                    printf("Nodo i-(n1) (%d) Nodo incognita\n", ii-(n1));
                    A(i,i-(n1-2)) = b(xx(ii-n1,0), xx(ii-n1,1), xx(ii-n1,2))/(h2*h2) - e(xx(ii-n1,0), xx(ii-n1,1), xx(ii-n1,2))/(2*h2);
                    break;

                case Dirichlet:
                    printf("Nodo i-(n1) (%d) Dirichlet\n", ii-(n1));
                    zz[i] = zz[i] - valor_frontera(xx(ii-n1,0), xx(ii-n1,1), xx(ii-n1,2), a_,b_,c_,d_,e_,f_)*(b(xx(ii-n1,0), xx(ii-n1,1), xx(ii-n1,2))/(h2*h2) - e(xx(ii-n1,0), xx(ii-n1,1), xx(ii-n1,2))/(2*h2));
                    break;

                case Neumann:
                    break;

                default:
                    fprintf(stderr, "Nodo desconocido %d\n", ii-n1);
                    exit(0);
                }


                // Nodo i+n1
                switch((int)xx(ii+n1,4)) {
                case Interno:
                    printf("Nodo i+n1 (%d) Nodo incognita\n", ii+n1);
                    A(i,i+(n1-2)) = b(xx(ii+n1,0), xx(ii+n1,1), xx(ii+n1,2))/(h2*h2) + e(xx(ii+n1,0), xx(ii+n1,1), xx(ii+n1,2))/(2*h2);
                    break;

                case Dirichlet:
                    printf("Nodo i+n1 (%d) Dirichlet\n", ii+n1);
                    zz[i] = zz[i] - valor_frontera(xx(ii+n1,0), xx(ii+n1,1), xx(ii+n1,2), a_,b_,c_,d_,e_,f_)*(b(xx(ii+n1,0), xx(ii+n1,1), xx(ii+n1,2))/(h2*h2) + e(xx(ii+n1,0), xx(ii+n1,1), xx(ii+n1,2))/(2*h2));
                    break;

                case Neumann:
                    break;

                default:
                    fprintf(stderr, "Nodo desconocido %d\n", ii+n1);
                    exit(0);
                }


                if (NDIM == 3) {
                    // Nodo i-n1n2
                    switch((int)xx(ii-(n1*n2),4)) {
                    case Interno:
                        printf("Nodo i-n1n2 (%d) Nodo incognita\n", ii-n1*n2);
                        A(i,i-(n1-2)*(n2-2)) = c(xx(ii-n1*n2,0), xx(ii-n1*n2,1), xx(ii-n1*n2,2))/(h3*h3) - f(xx(ii-n1*n2,0), xx(ii-n1*n2,1), xx(ii-n1*n2,2))/(2*h3);
                        break;

                    case Dirichlet:
                        printf("Nodo i-n1n2 (%d) Dirichlet\n", ii-(n1*n2));
                        zz[i] = zz[i] - valor_frontera(xx(ii-n1*n2,0),xx(ii-n1*n2,1),xx(ii-n1*n2,2), a_,b_,c_,d_,e_,f_)*(c(xx(ii-n1*n2,0), xx(ii-n1*n2,1), xx(ii-n1*n2,2))/(h3*h3) - f(xx(ii-n1*n2,0), xx(ii-n1*n2,1), xx(ii-n1*n2,2))/(2*h3));
                        break;

                    case Neumann:
                        break;

                    default:
                        fprintf(stderr, "Nodo desconocido %d\n", ii-n1*n2);
                        exit(0);
                    }


                    // Nodo i-n1n2
                    switch((int)xx(ii+n1*n2,4)) {
                    case Interno:
                        printf("Nodo i+n1n2 (%d) Nodo incognita\n", ii+n1*n2);
                        A(i,i+(n1-2)*(n2-2)) = c(xx(ii+n1*n2,0), xx(ii+n1*n2,1), xx(ii+n1*n2,2))/(h3*h3) + f(xx(ii+n1*n2,0), xx(ii+n1*n2,1), xx(ii+n1*n2,2))/(2*h3);
                        break;

                    case Dirichlet:
                        printf("Nodo i+n1n2 (%d) Dirichlet\n", ii+n1*n2);
                        zz[i] = zz[i] - valor_frontera(xx(ii+n1*n2,0),xx(ii+n1*n2,1),xx(ii+n1*n2,2), a_,b_,c_,d_,e_,f_)*(c(xx(ii+n1*n2,0), xx(ii+n1*n2,1), xx(ii+n1*n2,2))/(h3*h3) + f(xx(ii+n1*n2,0), xx(ii+n1*n2,1), xx(ii+n1*n2,2))/(2*h3));
                        break;

                    case Neumann:
                        break;

                    default:
                        fprintf(stderr, "Nodo desconocido %d\n", ii-n1*n2);
                        exit(0);
                    }
                }
            }
            i++;

        }

    }
    std::cout << "espacio "<< espacio << std::endl;
    std::cout << "Matriz"<< w << std::endl;
    std::cout << A << std::endl;
    std::cout << "zz "<< zz << std::endl;

    // Si la función tiene solución analítica
    // para verificar el resultado dado por
    // el método de diferencias finitas
    if (SOL != 0) {
        std::vector<double> sol(nNodosInternos);
        for (i = 0; i < nNodosInternos; i++)
            sol[i] = sol_analitica(espacio(i, 0), espacio(i, 1), espacio(i, 2));
        std::cout << "sol "<< sol << std::endl;
    }

    if (SIMETRICA) {
      gmm::identity_matrix PS;   // Optional scalar product for cg
      gmm::identity_matrix PR;   // Optional preconditioner
      gmm::iteration iter(10E-6);// Iteration object with the max residu
      //gmm::copy(A,AA);

      gmm::cg(A, w, zz, PS, PR, iter); // Conjugate gradient
      std::cout << "CGM "<< w << std::endl;

    } else {
      // computation of a preconditioner (ILUT)
      size_t restart = 50;       // restart parameter for GMRES
      gmm::iteration iter(10E-6);// Iteration object with the max residu
      gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > Pre(A, 10, 1e-4);
      gmm::gmres(A, w, zz, Pre, restart, iter);  // execute the GMRES algorithm
      std::cout << "GMRES preconditiones ILUT"<< w << std::endl;
    }

    print(xx, nNodos, w, nNodosInternos, a_,b_,c_,d_,e_,f_);
    //if (NDIM == 2)
      //  print_matlab(xx, nNodos, w, nNodosInternos, a_,b_,c_,d_,e_,f_, n1, n2, x, y);

    return 0;
}
