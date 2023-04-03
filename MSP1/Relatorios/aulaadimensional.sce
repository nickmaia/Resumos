//------------------------------------------------------------------------------
//                                DESCRIÇÃO
//------------------------------------------------------------------------------
//
//Autor: Nicole Maia Argondizzi
//Data: 09/03/2023
//------------------------------------------------------------------------------
//                             PRÉ-PROCESSAMENTO
//------------------------------------------------------------------------------
clc
clear
mode(-1)
//------------------------------------------------------------------------------
//                             FUNÇÕES
//------------------------------------------------------------------------------
//Modelo dimensional
function f=ee(x,F,A1,A2,B1,B2)
    h1=x(1)
    h2=x(2)
    f(1)=F/A1-(B1*sqrt(h1-h2))/A1
    f(2)=(B1*sqrt(h1-h2))/A2-(B2*sqrt(h2))/A2
endfunction

function dxdt=din(t,x,F,A1,A2,B1,B2)
    h1=x(1)
    h2=x(2)
    dxdt(1)=F/A1-(B1*sqrt(h1-h2))/A1
    dxdt(2)=(B1*sqrt(h1-h2))/A2-(B2*sqrt(h2))/A2
endfunction

//modelo adimensional 
function f = ee_ad(x,B1,B2,sol,F,A1,A2)
    h1star=x(1)
    h2star=x(2)
    f(1) = 1-B1 * sqrt(sol(1)) * sqrt(h1star-h2star)/F
    f(2) = B1* A1 *sqrt(sol(1)) * sqrt(h1star-h2star)/(F*A2) - B2*A1*sqrt(sol(1))*sqrt(h2star)/(F*A2)
endfunction

function f = din_ad(t,x,B1,B2,sol,F,A1,A2)
    h1star=x(1)
    h2star=x(2)
    f(1) = 1-B1 * sqrt(sol(1)) * sqrt(h1star-h2star)/F
    f(2) = B1* A1 *sqrt(sol(1)) * sqrt(h1star-h2star)/(F*A2) - B2*A1*sqrt(sol(1))*sqrt(h2star)/(F*A2)
endfunction
//------------------------------------------------------------------------------
//                              PROGRAMA PRINCIPAL
//------------------------------------------------------------------------------

//----------->Entrada de dados
F=5;
A1=5;
A2=10;
B1=2.5;
B2=5/sqrt(6);
//----------->Processamento
//Modelo dimensional
lista=list(ee,F,A1,A2,B1,B2)
xchute=[8;5];
[sol,v,info]=fsolve(xchute,lista)
if info==1 then
    disp('Solução:')
    disp(sol)
else
    disp('Tente novamente!')
end

lista2=list(din,F,A1,A2,B1,B2)
xinicial=[0;0]
t=[0:0.1:300]
x=ode(xinicial,t(1),t,lista2)

//Adimensional
//EE
xchute_ad=xchute/sol(1)
lista3=list(ee_ad,B1,B2,sol,F,A1,A2)
[sol_ad,v,info]=fsolve(xchute_ad,lista3)
if info==1 then
    disp('Solução adimensional:')
    disp(sol_ad)
else
    disp('Tente novamente')
end

//dinâmico


lista4=list(din_ad,B1,B2,sol,F,A1,A2)
xinicial_ad=xinicial/sol(1)
tau=t*F/(A1*sol(1))
x_ad=ode(xinicial_ad,tau(1),tau,lista4)

//----------->Saída de dados
scf(0)
clf()
subplot(2,2,1),plot(t,x(1,:),'r')
xtitle('Dimensional','$t[min]$','$h_1[ft]$')
subplot(2,2,2),plot(t,x(2,:),'b')
xtitle('Dimensional','$t[min]$','$h_2[ft]$')
subplot(2,2,3),plot(tau,x_ad(1,:),'r')
xtitle('Adimensional','$tau[-]$','$h_1^*[-]$')
subplot(2,2,4),plot(tau,x_ad(2,:),'b')
xtitle('Adimensional','$tau[-]$','$h_2^*[-]$')
//----------->Fim do programa
disp('***FIM***')



