// Ejemplo para la ecuación de onda en 1D
// u_tt - 4 u_xx = 0
// x in [0,1]
// u(0,t) = u(1,t) = 0
// u(x,0) = sen(%pi * x)
// u_t(x,0) = 0

// Solución analítica
// u(x,t) = sen(%pi * x) * cos(2* %pi * t)

// Dominio
a_ = -1
b_ = 1
// Particion en x
m = 1001; // número de nodos
h = (b_ - a_)/(m-1)
dt = 0.001 // salto del tiempo

// Pasos de tiempo
N = 500

// Para que sea estable se debe cumplir que
// sqrt(a) <= h/dt

// Coeficiente
function y = a(x)
    y = -4;
endfunction

// Condición inicial
function y = inicial(x)
    y = sin(%pi * x)
endfunction

function y = u_t(x)
    y = 0;
endfunction

// Solución analítica
function y = analitica(x,t)
    y = sin(%pi * x) * cos(2* %pi * t)
endfunction


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
////////////////////////////////////
// Aw^(j+1) = Bw^j
// creo el vector w_sol donde se guradará la solución para cada tiempo
// B matriz del lado derecho
// B_prima para guardar el resultado de Bw^j
w_sol = zeros(m,1); //
w = zeros(m,1); //
w_ = zeros(m,1); //
B = zeros(m,m); //
espacio = zeros(m,1) //
sol = zeros(m,1); //


// primer renglon de cada matriz
B(1,1) = 2*a(a_)*dt*dt/(h*h)
B(1,2) = -a(a_)*dt*dt/(h*h)

// se completa las matrices desde el renglon 2 hasta el m-1
for i = 2:m-1    
    B(i, i-1) = -a(i*h)*dt*dt/(h*h)
    B(i,i) = 2*a(i*h)*dt*dt/(h*h)
    B(i,i+1) = -a(i*h)*dt*dt/(h*h)
end

// último renglon de cada matriz
B(m,m-1) = -a(b_)*dt*dt/(h*h)
B(m,m) = 2*a(b_)*dt*dt/(h*h)

// muestro la matriz
//printf("Matriz B\n");
//disp(B);

for i=1:m
    xx = (i-1)*h;
    espacio(i) = a_ + xx;

    w(i) = inicial(espacio(i)); // Condiciones iniciales
    w_(i) = inicial(espacio(i)) + u_t(espacio(i)) * dt
end
//
//disp(espacio)
//disp(w)

////////////
//Para t = 0
B_prima = B * w;
for i = 1:m
    w_sol(i) = w_(i) + B_prima(i);
end
////////////
//printf("w para t = 0\n");
//disp(w_sol);

for i = 1:m
    sol(i) = analitica(espacio(i), 0)
end

//printf("Solucion analitica para t = 0\n")
//disp(sol)

//plot(espacio,w_sol)
//plot(espacio,sol,'r')



w_sol_temp = w_sol;
w_temp = w

for i=1:N
    t = i*dt
    B_prima = B * w_sol_temp;
    w_ = 2 * w_sol_temp - w_temp

    w_sol = w_ + B_prima; 
    //
    //    if i == 130 | i == 390 then
    //        [norma, ind] = normaInfinito(w_sol,sol,m);
    //        printf("tiempo t =  %f\n", t);
    //        printf("norma infinito %.8f\nEn el índice %d\n", norma, ind);
    //    end   

    w_temp = w_sol_temp
    w_sol_temp = w_sol

    if i == 5 | i == 50 | i == 100 | i == 150 | i == 200 | i == 250 | i == 300 | i== 350 | i == 400 | i == 450 | i == 500 then
        plot(espacio,w_sol)
        for j = 1:m      
            sol(j) = analitica(espacio(j), t)  
        end 
        plot(espacio,sol,'r')
    end


end





