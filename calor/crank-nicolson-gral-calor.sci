// Crank-Nicolson
// Para ecuaciones del tipo
// u_t - a(x)u_xx + b(x)u_x +c(x)u = f

// Ejemplo
// Ecuación de calor 1D
// u_t = u_xx

// 0<x<l 
// 0<t<T
// u(0,t) = u(l,t) = constante  0 < t < T cond de frontera
// u(x,0) = g(x) 0 <=x <= l

// Datos de entrada
// extremo l
// tiempo máximo T
// la constante alpha
// entero m>=3 
// entero N>=1

// Salida
// aproximaciones w_ij a u(x_i,t_j) para toda i=1,...,m-1 y j=1,2,...,N


//
// Funciones de los coeficientes
function y = a(x)
    y = 1;
endfunction

function y = b(x)
    y = 0;
endfunction

function y = c(x)
    y = 0;
endfunction

function y = f(x)
    y = 0;
endfunction

//
// función de la condicón inicial
function y = g(x)
    y = sin(x * %pi);
endfunction

function y = solucionAnalitica(x,t)
    y = exp(-2*t*%pi)*sin(%pi*x)
endfunction

// Datos de entrada
l = 1; // intervalo [0,1]
T = 0.5; // intervalo del tiempo [0,0.5]
m = 200; 
N = 50;

// Abrimos el archivo para imprimir
// en caso de querer guardar los valores de la matriz
// que se genera en la discretización, así como 
// los valores de la solución
//[fd,err]=mopen("CN_Gral_Calor.txt", "w");
// Función norma infinito

function [m, indice] = normaInfinito(v,w, dim)
    indice = 1
    m = abs(v(1) - w(1));
    
    for i = 2:dim
        m2 = abs(v(i) - w(i))
        if m2 > m then
            m = m2
            indice = i
        end
    end
    
endfunction

// División del espacio y del tiempo
h = l/m;
k = T/N;

// Aw^(j+1) = Bw^j + f^j
// creo el vector w donde se guradará la solución para cada tiempo
// A matriz del lado izquierdo
// B matriz del lado derecho
// B_prima para guardar el resultado de Bw^j
// ff que es el vector de los valores de f en cada nodo
w = zeros(m-1,1);
ff = zeros(m-1,1)
A = zeros(m-1,m-1);
B = zeros(m-1,m-1);
B_prima = zeros(m,1);

espacio = zeros(m-1,1);
sol_analitica = zeros(m-1,1);
// Paso 2
for i=1:m-1
    xx = i*h;
    espacio(i) = xx;
    ff(i) = f(i);
    w(i) = g(xx) + ff(i);
end


// primer renglon de cada matriz
A(1,1) = 1/k + a(h)/(h*h); // falta condición de frontera
A(1,2) = b(h)/(4*h) - a(h)/(2*h*h);

B(1,1) = 1/k - a(h)/(h*h) - c(h);
B(1,2) = a(h)/(2*h*h) - b(h)/(4*h);

// se completa las matrices desde el renglon 2 hasta el m-1
for i = 2:m-2
    A(i, i-1) = -b(i*h)/(4*h) - a(i*h)/(2*h*h);
    A(i,i) = 1/k + a(i*h)/(h*h);
    A(i,i+1) = b(i*h)/(4*h) - a(i*h)/(2*h*h);

    B(i, i-1) = a(i*h)/(2*h*h) + b(i*h)/(4*h);
    B(i,i) = 1/k - a(i*h)/(h*h) - c(i*h);
    B(i,i+1) = a(i*h)/(2*h*h) - b(i*h)/(4*h);
end

// último renglon de cada matriz
A(m-1,m-2) = -b(l-h)/(4*h) - a(l-h)/(2*h*h); 
A(m-1,m-1) = 1/k + a(l-h)/(h*h);

B(m-1,m-2) = a(l-h)/(2*h*h) + b(l-h)/(4*h);
B(m-1,m-1) = 1/k - a(l-h)/(h*h) - c(l-h);

// muestro las matrices
//printf("Matriz A\n");
//disp(A);
//printf("Matriz B\n");
//disp(B);

// tiempo inicial 
//mfprintf(fd, "t = 0.0\n");
//for kk = 1:m-1
//    mfprintf(fd, "w(%d) %f\n", kk, w(kk));
//end

// Resolvemos el sistema iterativamente
for j=1:N
    t = j*k;
    B_prima = B * w + ff;
    w = inv(A) * B_prima;

    // Cálculo del vector de la solución analítica para calcular  
    // la norma infinito
    for i = 1:m-1
        sol_analitica(i) = solucionAnalitica(espacio(i),t)
    end
    //disp(w)
    plot(espacio, w);
    //    plot(espacio, sol_analitica)

    if j == 13 | j == 39 then
        [norma, ind] = normaInfinito(w,solucionAnalitica,m-1);
        printf("tiempo t =  %f\n", t);
        printf("h %f\nk %f\n", h, k)
        printf("norma infinito %.8f\nEn el índice %d\n", norma, ind);
    end
    // mfprintf(fd,"t = %f\n", t);    
    // for kk = 1:m-1
    //  mfprintf(fd, "w(%d) %f\n", kk, w(kk));
    // end
end

// Cerramos el archivo
//mclose(fd);
