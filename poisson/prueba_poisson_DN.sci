// Diferencias Finitas

// Ecuación
// Uxx = -pi*pi*cos(pi*x)
// x in [0,0.5]
// u(0) = 1
// u(0.5) = -pi

//
function y=g(x)
    y = -%pi*%pi*cos(%pi*x);
endfunction
// Solución analitica
function y = sol_analitica(punto)
    y = cos(%pi*punto)
endfunction

// dominio
a_ = 0;
b_ = 0.5;

N = 11; // Número de nodos

// tamaño del intervalo
h = (b_ - a_) / (N-1);

// Para resolver Ax=y
A = zeros(N-1, N-1);
y = zeros(N-1,1);


// Primer renglon de la matriz A y 
// vector y
A(1,1)= -2/(h*h);
A(1,2)= 1/(h*h);
y(1)= g(a_ + h) - 1/(h*h); // Frontera Dirichlet

// Renglones intermedios de la matriz A y 
// vector b
for i=2:N-2
 A(i,i-1)=1/(h*h);
 A(i,i)=-2/(h*h);
 A(i,i+1)=1/(h*h)
 y(i)= g(i*h);
end

// Renglon final de la matriz A y 
// vector y
// debidos a la frontera Neumann
A(N-1,N-2) = -1/h;
A(N-1,N-1)= 1/h;
y(N-1)= -%pi;

// Resuelve el sistema lineal Ax=y
printf("Matriz A\n")
disp(A)
printf("Vector y\n")
disp(y)
x=inv(A)*y;
printf("Vector Solucion x\n")
disp(x)


// Prepara la graficacion
// con la solución analítica
xx=zeros(N-1,1);
zz=zeros(N-1,1);
for i=1:N-1
 xx(i)=i*h;
 zz(i)=sol_analitica(xx(i));
end
printf("vector solucion analitica\n")
disp(zz)


// Grafica la solucion de la Ecuacion Diferencial Parcial en 1D
// Solucion de Dif Finitas
plot(xx,x)    // azul
// Solucion Analitica
plot(xx,zz, 'r')    // roja