clc
clear
mode(-1)


//---------------------------------------------------------------------
// Nomes:
//Nicole Maia Argondizzi
//
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//           Metodo_Newton_Raphson_para_sistema_de_equações            

// Dados Utilizados
Fe = 100 // L/min
Fs = 100 // L/min
k = 0.01 // L^(1/2)/mol^(1/2)*min
T = 500 // K
P = 1485 // kPa
V = 5000 // L

// CONCENTRAÇÕES INICIAIS
Cae = 0.1 // mol/L
Cbe = 0.054 // mol/L
Cce = 0 // mol/L
Cie = 0.2032 // mol/L



// FUNÇÕES 1
function f1 = func1(Ca, Cb, Cc, Ci)
    f1 = ((Cae*Fe - Ca*Fs)/(V)) - (k*Ca*Cb^(1/2))
endfunction

function f2 = func2(Ca, Cb, Cc, Ci)
    f2 = ((Cbe*Fe - Cb*Fs)/(V)) - (k*Ca*Cb^(1/2)*1/2)
endfunction

function f3 = func3(Ca, Cb, Cc, Ci)
    f3 = ((- Cc*Fs)/(V)) + (k*Ca*Cb^(1/2))
endfunction

function f4 = func4(Ca, Cb, Cc, Ci)
    f4 = ((Cie*Fe-Ci*Fs)/(V))
endfunction

// FUNÇÕES 2
function f1 = func1_(Xi)
    Ca = Xi(1)
    Cb = Xi(2)
    Cc = Xi(3)
    Ci = Xi(4)
    f1 = ((Cae*Fe - Ca*Fs)/(V)) - (k*Ca*Cb^(1/2))
endfunction

function f2 = func2_(Xi)
    Ca = Xi(1)
    Cb = Xi(2)
    Cc = Xi(3)
    Ci = Xi(4)
    f2 = ((Cbe*Fe - Cb*Fs)/(V)) - (k*Ca*Cb^(1/2)* (1/2))
endfunction

function f3 = func3_(Xi)
    Ca = Xi(1)
    Cb = Xi(2)
    Cc = Xi(3)
    Ci = Xi(4)
    f3 = ((- Cc*Fs)/(V)) + (k*Ca*Cb^(0.5))
endfunction

function f4 = func4_(Xi)
    Ca = Xi(1)
    Cb = Xi(2)
    Cc = Xi(3)
    Ci = Xi(4)
    f4 = ((Cie*Fe-Ci*Fs)/(V))
endfunction

function J = jacobiano(Xi, erro)
    Ca = Xi(1)
    Cb = Xi(2)
    Cc = Xi(3)
    Ci = Xi(4)
    J(1,:) = [(func1(Ca+erro, Cb, Cc, Ci)-func1(Ca, Cb, Cc, Ci))/erro, ...
              (func1(Ca, Cb+erro, Cc, Ci)-func1(Ca, Cb, Cc, Ci))/erro, ...
              (func1(Ca, Cb, Cc+erro, Ci)-func1(Ca, Cb, Cc, Ci))/erro, ...
              (func1(Ca, Cb, Cc, Ci+erro)-func1(Ca, Cb, Cc, Ci))/erro]
              
     J(2,:) = [(func2(Ca+erro, Cb, Cc, Ci)-func2(Ca, Cb, Cc, Ci))/erro, ...
               (func2(Ca, Cb+erro, Cc, Ci)-func2(Ca, Cb, Cc, Ci))/erro, ...
               (func2(Ca, Cb, Cc+erro, Ci)-func2(Ca, Cb, Cc, Ci))/erro, ...
               (func2(Ca, Cb, Cc, Ci+erro)-func2(Ca, Cb, Cc, Ci))/erro]
               
     J(3,:) = [(func3(Ca+erro, Cb, Cc, Ci)-func3(Ca, Cb, Cc, Ci))/erro, ...
               (func3(Ca, Cb+erro, Cc, Ci)-func3(Ca, Cb, Cc, Ci))/erro, ...
               (func3(Ca, Cb, Cc+erro, Ci)-func3(Ca, Cb, Cc, Ci))/erro, ...
               (func3(Ca, Cb, Cc, Ci+erro)-func3(Ca, Cb, Cc, Ci))/erro]
               
     J(4,:) = [(func4(Ca+erro, Cb, Cc, Ci)-func4(Ca, Cb, Cc, Ci))/erro, ...
               (func4(Ca, Cb+erro, Cc, Ci)-func4(Ca, Cb, Cc, Ci))/erro, ...
               (func4(Ca, Cb, Cc+erro, Ci)-func4(Ca, Cb, Cc, Ci))/erro, ...
               (func4(Ca, Cb, Cc, Ci+erro)-func4(Ca, Cb, Cc, Ci))/erro]
endfunction

erro = 1e-5
Xi = [Cae; Cbe; Cce; Cie]
Y = [func1_(Xi); func2_(Xi); func3_(Xi); func4_(Xi)]
J = jacobiano(Xi, erro)
Xf = Xi - inv(J)*Y
contador = 0

while abs(max(Xf - Xi)) > erro do
    contador = contador + 1
    Xi = Xf
    Y = [func1_(Xi); func2_(Xi); func3_(Xi); func4_(Xi)]
    J = jacobiano(Xi, erro)
    Xf = Xi - inv(J)*Y
end

mprintf(" Concentrações em estado estacionário:\n")
mprintf(" SO2: %f mol/L\n O2: %f mol/L\n SO3: %f mol/L\n Inerte: Nitrogênio[N2]: %f mol/L\n", Xf(1), ...
        Xf(2), Xf(3), Xf(4))
mprintf("\n")
mprintf(" Prova real de que as funções zeraram")
disp([func1_(Xf); func2_(Xf); func3_(Xf); func4_(Xf)])
mprintf("\n")
mprintf(" O número de iterações foram: %f\n", contador)
//---------------------------------------------------------------------
//           Estudo da relação dos volumes e concentrações             

V_vector = [5000:100:100000]
for i = 1:1:length(V_vector)
    V = V_vector(i)
    erro = 1e-5
    Xi = [Cae; Cbe; Cce; Cie]
    Y = [func1_(Xi); func2_(Xi); func3_(Xi); func4_(Xi)]
    J = jacobiano(Xi, erro)
    Xf = Xi - inv(J)*Y
    contador = 0
    
    while abs(max(Xf - Xi)) > erro do
        contador = contador + 1
        Xi = Xf
        Y = [func1_(Xi); func2_(Xi); func3_(Xi); func4_(Xi)]
        J = jacobiano(Xi, erro)
        Xf = Xi - inv(J)*Y
    end
    
    result(:,i) = Xf
end

scf(0)
clf(0)
xtitle("Concentração de cada componente", "$Volume [L]$", "$Concentração [mol/L]$" )
plot(V_vector, result(1,:), "r")
plot(V_vector, result(2,:), "b")
plot(V_vector, result(3,:), "y")
plot(V_vector, result(4,:), "g")
legend(["Dióxido de enxofre [SO2]", "Oxigênio [O2]", "Trióxido de Enxofre [SO3]", "Inerte:Nitrogênio[N2]"],5)


























 






