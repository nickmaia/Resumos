//------------------------------------------------------------------------------
//                                DESCRIÇÃO
//------------------------------------------------------------------------------
//Relatório 01
//Autor: Nicole Maia Argondizzi
//Autor: Carina Estela Aquino
//Autor: Iara Solimar da Silva

//------------------------------------------------------------------------------
//                             PRÉ-PROCESSAMENTO
//------------------------------------------------------------------------------
clc
clear
mode(-1)

//------------------------------------------------------------------------------
//                             COLETANDO DADOS
//------------------------------------------------------------------------------
Cao = 1 //mol/L
k = 1.3 //L^n-1 mol^n-1 /s

// Ordem da reação

n1 = 1
n2 = 1.5
n3 = 2
n4 = 0.5

//------------------------------------------------------------------------------
//                             CONDIÇÕES INICIAIS
//------------------------------------------------------------------------------
Cai = Cao
Cbi = 0

//------------------------------------------------------------------------------
//                             FUNÇÕES
//------------------------------------------------------------------------------

//-------------------------Estado estacionário
function f=ee(xi)
     
    Ca1 = xi(1)
    Ca2 = xi(2)
    Ca3 = xi(3)
    Ca4 = xi(4)
    
    Cb1 = xi(5)
    Cb2 = xi(6)
    Cb3 = xi(7)
    Cb4 = xi(8)
    
    //Ca
    f(1)= -k * Ca1^n1
    f(2)= -k * Ca2^n2
    f(3)= -k * Ca3^n3
    f(4)= -k * Ca4^n4
    
    //Cb
    f(5)= k * Ca1^n1
    f(6)= k * Ca2^n2
    f(7)= k * Ca3^n3
    f(8)= k * Ca4^n4
    
    f=real(f)
    
endfunction

//Determinando a solução para o estado Estacionário

xchute = [Cai;Cai;Cai;Cai;Cbi;Cbi;Cbi;Cbi] 
[sol,v,info]=fsolve(xchute,ee)


if info==1 then
    mprintf(' Solução estado estacionário: \n\n')
    // n=1
    mprintf(' Ca1 : %f\n', sol(1))
    mprintf(' Cb1 : %f\n', sol(5))
    // n=1.5
    mprintf(' Ca2 : %f\n', sol(2))
    mprintf(' Cb2 : %f\n', sol(6))
    // n=2
    mprintf(' Ca3 : %f\n', sol(3))
    mprintf(' Cb3 : %f\n', sol(7))
    // n=0.5
    mprintf(' Ca4 : %f\n', sol(4))
    mprintf(' Cb4 : %f\n', sol(8))
    
else 
    disp("Tente novamente!")
end



//-----------------------------Estado Dinâmico --------------------------------
function f=din(t,xi)
    
    
    Ca1 = xi(1)
    Ca2 = xi(2)
    Ca3 = xi(3)
    Ca4 = xi(4)
    
    Cb1 = xi(5)
    Cb2 = xi(6)
    Cb3 = xi(7)
    Cb4 = xi(8)
    
    //Ca
    f(1)= -k * Ca1^n1
    f(2)= -k * Ca2^n2
    f(3)= -k * Ca3^n3
    f(4)= -k * Ca4^n4
    
    //Cb
    f(5)= k * Ca1^n1
    f(6)= k * Ca2^n2
    f(7)= k * Ca3^n3
    f(8)= k * Ca4^n4
    
    f=real(f)
    
endfunction

// Determinando resultados para diferentes n's

t=0:0.1:30
xinicial=[Cai;Cai;Cai;Cai;Cbi;Cbi;Cbi;Cbi] 
x=ode(xinicial,t(1),t,din)

//-----------------------------------------------------------------------------
//                          Conversão 𝑋𝐴's
//

xa1 = (Cao-x(1,:))./Cao
xa2 = (Cao-x(2,:))./Cao
xa3 = (Cao-x(3,:))./Cao
xa4 = (Cao-x(4,:))./Cao

//---------------------------------------------------------------------------
//           Gráfico da concentração com seu respectivo n
//

scf()
clf()

// n=1
subplot(1,1,1),plot(t,x(1,:),"m")
subplot(1,1,1),plot(t,x(5,:),"b")
xtitle('Ci [mol/L] vs Tempo [s] n=1','t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"])

scf()
clf()

// n=1.5
subplot(1,1,1),plot(t,x(2,:),"m")
subplot(1,1,1),plot(t,x(6,:),"b")
xtitle('Ci [mol/L] vs Tempo [s] n=1.5','t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"])

scf()
clf()
//n=2
subplot(1,1,1),plot(t,x(3,:),"m")
subplot(1,1,1),plot(t,x(7,:),"b")
xtitle('Ci [mol/L] vs Tempo [s] n=2','t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"],5)

scf()
clf()
// n=0.5
subplot(1,1,1),plot(t,x(4,:),"m")
subplot(1,1,1),plot(t,x(8,:),"b")
xtitle('Ci [mol/L] vs Tempo [s] n=0.5','t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"])

//----------------------------------------------------------------------------
//                Gráfico da Conversão para cada n determinado
//

scf()
clf()
xgrid(0)
subplot(1,1,1),plot(t,xa1,"m") //n=1
subplot(1,1,1),plot(t,xa2,"r") //n=1.5
subplot(1,1,1),plot(t,xa3,"b") //n=2
subplot(1,1,1),plot(t,xa4,"g") //n=0.5
xtitle('Conversão [Xa] vs Tempo [s]','t[s]','Xai')
legend(['Xa1','Xa2','Xa3','Xa4'],4)
//-----------------------------------------------------------------------------
// Questão 05 
//-----------------------------------------------------------------------------
// Determinando a resposta para 3 tipos de Cao's
Cao1=0.5
Cao2=1
Cao3=1.5

// Condições iniciais
Cai1 = Cao1
Cai2 = Cao2
Cai3 = Cao3

n=2
//-----------------------------------------------------------------------------
//                            Estado Dinâmico 
//
function f=y(t,xi)
    
    Ca1 = xi(1)
    Ca2 = xi(2)
    Ca3 = xi(3)
    
    //Cao's
    f(1)= -k * Ca1^n
    f(2)= -k * Ca2^n
    f(3)= -k * Ca3^n
    
    f=real(f)
endfunction

t=0:0.1:10
xinicial=[Cai1;Cai2;Cai3] 
xi=ode(xinicial,t(1),t,y)

//-----------------------------Determinando Xa 

xa1i = (Cao1-xi(1,:))./Cao1
xa2i = (Cao2-xi(2,:))./Cao2
xa3i = (Cao3-xi(3,:))./Cao3

//----------------------------------Gráficos-----------------------------------
scf()
clf()
xgrid(0)
subplot(1,1,1),plot(t,xi(1,:),"m")
subplot(1,1,1),plot(t,xi(2,:),"b")
subplot(1,1,1),plot(t,xi(3,:),"g")

xtitle('Ci [mol/L] vs Tempo [s]','t[s]','Cai[mol/L]')
legend(["Cao:0.5 [mol/L]","Cao:1 [mol/L]","Cao:1.5 [mol/L]"])

scf()
clf()
xgrid(0)
subplot(1,1,1),plot(t,xa1i,"m")
subplot(1,1,1),plot(t,xa2i,"b")
subplot(1,1,1),plot(t,xa3i,"g")

xtitle('Xai [mol/L] vs Tempo [s]','t[s]','Xai[mol/L]')
legend(["Cao:0.5 [mol/L]","Cao:1 [mol/L]","Cao:1.5 [mol/L]"],4)
