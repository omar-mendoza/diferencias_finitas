//Programa general
//Método de Diferencias Finitas
//Con Condiciones Neumman o Dirichlet


// Ecuaciones tipo
// a(x)Uxx + b(x)Ux + c(x)U = g(x)
// En un intervalo [alp, bet]

/////////////////////////////
// Def de los coeficientes //
/////////////////////////////

// Coeficientes de la ecuació diferencial
function y=a(x)
    y = 0;
endfunction

function y=b(x)
    y = 0.5;
endfunction

function y=c(x)
    y = 0;
endfunction

// 
function y=g(x)
    y = 0
endfunction


// Solución analítica para 
// comparar la solución numérica
function y=SolucionAnalitica(x)
    t = 2.0
 y= exp(-(x-t)^2 / 10*10);
endfunction

// Definición de la frontera
// Frontera
// Definir la frontera
// Dirichlet = 1
// Neumman = 2
function y = tipoFrontera(valNodo)
    
    //y = 1
    
    if valNodo == alp then
        y = 1;
    end
    
    if valNodo == bet  then
        y = 1;
    end
endfunction


// Valores de las fronteras
function y = valor_frontera(nodo, x)
    
    //y = 0;
    
    if nodo == 1 then
        y = 0;
    elseif nodo == nNodos
        y = 0;
    end
    
endfunction


// Dominio 
// x in [alpha, beta]
alp = 0;
bet = 10;


// Número total de nodos en el dominio
nNodos = 101;

// Valor del incremento
h = (bet-alp)/(nNodos-1);
printf("h = %f", h);

// Creamos la malla
// La primer entrada será la coordenada
// La segunda si es frontera, cero significa que es nodo interior
nodos = zeros(nNodos, 2)
nNodosIncognita = nNodos;

for i=1:nNodos
    nodos(i, 1) = alp + (i-1)*h
    if i==1 | i==nNodos then
        tipo = tipoFrontera(nodos(i,1))
        if tipo <> 0 & tipo == 1 then
            nodos(i,2) = tipo  
            nNodosIncognita = nNodosIncognita - 1;        
        end
        if tipo <> 0 & tipo == 2 then
            nodos(i,2) = tipo
                        
        end
    end
end

printf("Nodos\n");
disp(nodos)
printf("nNodosIncognita %d", nNodosIncognita)


// Creamos la matriz A 
// y el vector y para guardar los valores
//y resolver el sistema
// A*x = y
A = zeros(nNodosIncognita, nNodosIncognita);
y = zeros(nNodosIncognita);

// Vector para la solución analítica
zz = zeros(nNodosIncognita, 1);
espacio = zeros(nNodosIncognita, 1);

//
// Armamos la matriz y el vector
//

i = 1;// Contador de los renglones de la matriz 'A' y del vector 'y'

for j=1:nNodos
      
    if nodos(j,2) == 0 then // el nodo j no es frontera
         printf("\nValor de nodos j %d %d\n", j, nodos(j,2));
        
        espacio(i) = nodos(j,1);
        zz(i)=SolucionAnalitica(nodos(j,1));
        
        A(i, i) = (-2*a(nodos(j,1))/(h*h)) + c(nodos(j,1));
        y(i) = g(nodos(j,1));
        
        if (j-1) >= 1 then
            select nodos(j-1,2)
            case 0 then 
                 A(i,i-1) = a(nodos(j,1))/(h*h) - b(nodos(j,1))/(2*h);               
            case 1 then // Frontera Dirichlet por la izquierda
                y(i) = y(i) - valor_frontera(1, nodos(j,1))*(a(nodos(j,1))/(h*h) - b(nodos(j,1))/(2*h));
            end
        end
    
        if j < nNodos then
            select nodos(j+1,2)
            case 0 then
                A(i, i+1) = a(nodos(j,1))/(h*h) + b(nodos(j,1))/(2*h);
            case 1 then   // Frontera Dirichlet por la derecha 
                y(i) = y(i) - valor_frontera(1, nodos(j,1))*(a(nodos(j,1))/(h*h) + b(nodos(j,1))/(2*h));
            end
        end
       
        i = i + 1;
            
    elseif nodos(j,2) == 2 then
        espacio(i) = nodos(j,1);
        if j == 1 then //Frontera Neumann por la izquierda
            A(i,i) = -1/h;
            A(i,i+1) = 1/h;
            y(i) = valor_frontera(1, nodos(1,1))
            i = i + 1
        elseif j == nNodos then  // Frontera Neumann por la derecha
            A(i, i-1) = -1/h;
            A(i,i) = 1/h;
            y(i) = valor_frontera(nNodos, nodos(nNodos,1));
            i = i + 1;
        end
    end
end

// Resuelve el sistema lineal Ax=y
printf("Matriz A")
disp(A)
printf("Vector y")
disp(y)
x=inv(A)*y;
printf("Vector Soucion x")
disp(x)

// Mostramos la solución analítica
printf("Solución analítica");
disp(zz);
printf("Espacio\n");
disp(espacio);
plot(espacio,x);
plot(espacio,zz, 'r');
