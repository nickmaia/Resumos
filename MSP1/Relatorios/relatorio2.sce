//------------------------------------------------------------------------------
//                                DESCRIÇÃO
//------------------------------------------------------------------------------
//Relatório 02
//Autor: Nicole Maia Argondizzi


//------------------------------------------------------------------------------
//                             PRÉ-PROCESSAMENTO
//------------------------------------------------------------------------------
clc
clear
mode(-1)

//------------------------------------------------------------------------------
//                             COLETANDO DADOS
//------------------------------------------------------------------------------
Cao = 10 //kmol m^-3
Ea = 11843 //kcal kmol^-1
tau = 1 * 3600 //s
ko = 14825 //s^-1
R = 1.987 //kcal kmol^-1 K^-1
To = 298 //K
Tj = 298 //K
Ua_V = 250/3600 //kcal m^-3 s^-1 K^-1 
deltaH = -5215 //kcal kmol^-1
pcP = 500 //kcal m^-3 K^-1

//------------------------------------------------------------------------------
//                             DADOS INICIAIS
//------------------------------------------------------------------------------

Cai = Cao //kmol m^-3
Tri = To //K

//------------------------------------------------------------------------------
//                             FUNÇÕES
//----------------------------------------------------------------------------


function f=ee(xi)
     
    Ca = xi(1)   
    Tr = xi(2)
    
    //Ca
    f(1)= (Cao/tau) - (Ca/tau) - ((ko * exp(-Ea/(R*Tr))) * Ca)
    
    //Tr
    f(2) = ((To-Tr)/tau) - (((ko * exp(-Ea/(R*Tr)))* Ca * deltaH)/pcP) - (Ua_V * (Tr -Tj)/pcP)
    
    
    
endfunction

//---------------------------------------------------------------------------
//1. Determinando a solução para o estado Estacionário
//---------------------------------------------------------------------------
xchute = [Cai;Tri] 
[sol,v,info]=fsolve(xchute,ee)

//[1].Calcule o estado estacionário.

if info==1 then
    mprintf(' Solução estado estacionário: \n\n')
    // Ca
    mprintf(' Ca : %f\n', sol(1))
    //Tr
    mprintf(' Tr : %f\n', sol(2))

    
else 
    disp("Tente novamente!")
end


function f=din1(t,xi)
     
    Ca = xi(1)   
    Tr = xi(2)
    
    //Ca
    f(1)= (Cao/tau) - (Ca/tau) - ((ko * exp(-Ea/(R*Tr))) * Ca)
    
    //Tr
    f(2) = ((To-Tr)/tau) - (((ko * exp(-Ea/(R*Tr)))* Ca * deltaH)/pcP) - (Ua_V * (Tr -Tj)/pcP)
    
    f=real(f)
    
endfunction

t=0:0.1:10*3600
xinicial=[Cai;Tri] 
xi=ode(xinicial,t(1),t,din1)







//---------------------------------------------------------------------------
//2. O comportamento da concentração de A e da temperatura do reator com o tempo.
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//           Gráfico da concentração pelo tempo
//
//Ao longo do tempo, a concentração de A no reator diminui exponencialmente até alcançar o estado estacionário, o qual é atingido após 10 horas de reação.

scf()
clf()

subplot(1,1,1),plot(t,xi(1,:),"m")
xtitle('Ca [kmol/m^-3] vs Tempo [s]','t[s]','Ci[kmol/m^-3]')
legend(["Ca [kmol/m^-3]"])
xgrid(0)


//---------------------------------------------------------------------------
//           Gráfico da temperatura pelo tempo
//
//De acordo com o gráfico gerado, a temperatura começa em torno de 298 K e aumenta até alcançar um estado estacionário próximo de 315 K, como esperado para essa reação química.


scf()
clf()

subplot(1,1,1),plot(t,xi(2,:),"m")
xtitle('Tr [K] vs Tempo [s]','t[s]','Tr [K]')
legend(["Tr [K]"],4)
xgrid(0)







//---------------------------------------------------------------------------
//3. Avalie a influência da concentração da corrente de alimentação, CA0, nasvariáveis CA e TR. Discuta os resultados.
//---------------------------------------------------------------------------

// Concentração

Cao1 = 8
Cao2 = 10
Cao3 = 12

//Condições iniciais

Cai1 = Cao1
Cai2 = Cao2
Cai3 = Cao3

Tri = To

function f=din2(t,x)
     
    Ca1 = x(1)   
    Tr1 = x(2)
    
     
    Ca2 = x(3)   
    Tr2 = x(4)
    
     
    Ca3 = x(5)   
    Tr3 = x(6)
    
    
    //Cao1
    f(1)= (Cao1/tau) - (Ca1/tau) - ((ko * exp(-Ea/(R*Tr1))) * Ca1)
    
    //Tr1
    f(2) = ((To-Tr1)/tau) - (((ko * exp(-Ea/(R*Tr1)))* Ca1 * deltaH)/pcP) - (Ua_V * (Tr1 -Tj)/pcP)
    
    //Cao2
    f(3)= (Cao2/tau) - (Ca2/tau) - ((ko * exp(-Ea/(R*Tr2))) * Ca2)
    
    //Tr2
    f(4) = ((To-Tr2)/tau) - (((ko * exp(-Ea/(R*Tr2)))* Ca2 * deltaH)/pcP) - (Ua_V * (Tr2 -Tj)/pcP)
    
    //Cao3
    f(5)= (Cao3/tau) - (Ca3/tau) - ((ko * exp(-Ea/(R*Tr3))) * Ca3)
    
    //Tr3
    f(6) = ((To-Tr3)/tau) - (((ko * exp(-Ea/(R*Tr3)))* Ca3 * deltaH)/pcP) - (Ua_V * (Tr3 -Tj)/pcP)
    
    f=real(f)
    
endfunction

//10 horas
t=0:0.1:10*3600
xinicial=[Cai1;Tri;Cai2;Tri;Cai3;Tri] 
x=ode(xinicial,t(1),t,din2)


scf()
clf()

//Concentração
subplot(2,1,1),plot(t,x(1,:),"m")
subplot(2,1,1),plot(t,x(3,:),"g")
subplot(2,1,1),plot(t,x(5,:),"b")
xtitle('Ca [kmol/m^-3] vs Tempo [s]','t[s]','Ca[kmol/m^-3]')
legend(["Cao: 5 [kmol/m^-3]","Cao: 10 [kmol/m^-3]", "Cao: 15 [kmol/m^-3]"],3)
xgrid(0)

//Temperatura
subplot(2,1,2),plot(t,x(2,:),"m")
subplot(2,1,2),plot(t,x(4,:),"g")
subplot(2,1,2),plot(t,x(6,:),"b")
xtitle('Temperatura [K] vs Tempo [s]','t[s]','Toi [K]')
legend(["Cao: 5 [kmol/m^-3]","Cao: 10 [kmol/m^-3]", "Cao: 15 [kmol/m^-3]"],2)
xgrid(0)






//---------------------------------------------------------------------------
//4. Avalie a influência da temperatura da corrente de alimentação, T0, nas variáveis CA e TR. Discuta os resultados.
//---------------------------------------------------------------------------

// Temperaturas

To1 = 290 //K
To2 = 298 //K
To3 = 310 //K

//Condições iniciais

Cai = Cao

Tri1 = To1
Tri2 = To2
Tri3 = To3

function f=din3(ti,y)
     
    Ca1 = y(1)   
    Tr1 = y(2)
    
     
    Ca2 = y(3)   
    Tr2 = y(4)
    
     
    Ca3 = y(5)   
    Tr3 = y(6)
    
    
    //T01
    f(1)= (Cao/tau) - (Ca1/tau) - ((ko * exp(-Ea/(R*Tr1))) * Ca1)
    f(2) = ((To1-Tr1)/tau) - (((ko * exp(-Ea/(R*Tr1)))* Ca1 * deltaH)/pcP) - (Ua_V * (Tr1 -Tj)/pcP)
    
    //T02
    f(3)= (Cao/tau) - (Ca2/tau) - ((ko * exp(-Ea/(R*Tr2))) * Ca2)
    f(4) = ((To2-Tr2)/tau) - (((ko * exp(-Ea/(R*Tr2)))* Ca2 * deltaH)/pcP) - (Ua_V * (Tr2 -Tj)/pcP)
    
    //T03
    f(5)= (Cao/tau) - (Ca3/tau) - ((ko * exp(-Ea/(R*Tr3))) * Ca3)
    f(6) = ((To3-Tr3)/tau) - (((ko * exp(-Ea/(R*Tr3)))* Ca3 * deltaH)/pcP) - (Ua_V * (Tr3 -Tj)/pcP)
    
    f=real(f)
    
endfunction

//10 Horas

ti=0:0.1:12*3600
xinicial=[Cai;Tri1;Cai;Tri2;Cai;Tri3] 
y=ode(xinicial,t(1),ti,din3)


scf()
clf()

//Concentração
subplot(2,1,2),plot(ti,y(1,:),"m")
subplot(2,1,2),plot(ti,y(3,:),"g")
subplot(2,1,2),plot(ti,y(5,:),"b")
xtitle('Ca [kmol/m^-3] vs Tempo [s]','t[s]','Ca[kmol/m^-3]')
legend(["To: 290 [K]","To: 298 [K]", "To: 310 [K]"],4)
xgrid(0)

//Temperatura
subplot(2,1,1),plot(ti,y(2,:),"m")
subplot(2,1,1),plot(ti,y(4,:),"g")
subplot(2,1,1),plot(ti,y(6,:),"b")
xtitle('Temperatura [K] vs Tempo [s]','t[s]','Toi [K]')
legend(["To: 290 [K]","To: 298 [K]", "To: 310 [K]"],1)
xgrid(0)







//---------------------------------------------------------------------------
//5. Avalie a influência da temperatura do fluido refrigerante, Tj, nas variáveis CA e TR. Discuta os resultados.
//---------------------------------------------------------------------------

//Condições iniciais

Cai = Cao
Tri = To

//Temperatura do fluido refrigerante
Tj1 = 290 //K
Tj2 = 298 //K
Tj3 = 310 //K

function f=din4(t,yi)
     
    Ca1 = yi(1)   
    Tr1 = yi(2)
    
     
    Ca2 = yi(3)   
    Tr2 = yi(4)
    
     
    Ca3 = yi(5)   
    Tr3 = yi(6)
    
    
    //Tj1
    f(1)= (Cao/tau) - (Ca1/tau) - ((ko * exp(-Ea/(R*Tr1))) * Ca1)
    f(2) = ((To-Tr1)/tau) - (((ko * exp(-Ea/(R*Tr1)))* Ca1 * deltaH)/pcP) - (Ua_V * (Tr1 -Tj1)/pcP)
    
    //Tj2
    f(3)= (Cao/tau) - (Ca2/tau) - ((ko * exp(-Ea/(R*Tr2))) * Ca2)
    f(4) = ((To-Tr2)/tau) - (((ko * exp(-Ea/(R*Tr2)))* Ca2 * deltaH)/pcP) - (Ua_V * (Tr2 -Tj2)/pcP)
    
    //Tj3
    f(5)= (Cao/tau) - (Ca3/tau) - ((ko * exp(-Ea/(R*Tr3))) * Ca3)
    f(6) = ((To-Tr3)/tau) - (((ko * exp(-Ea/(R*Tr3)))* Ca3 * deltaH)/pcP) - (Ua_V * (Tr3 -Tj3)/pcP)
    
    f=real(f)
    
endfunction

//10 Horas

t=0:0.1:10*3600
xinicial=[Cai;Tri;Cai;Tri;Cai;Tri] 
yi=ode(xinicial,t(1),t,din4)

// ----------------------------------------------------------------------------
// --------------------------------Gráficos
// ----------------------------------------------------------------------------
scf()
clf()

//Concentração
subplot(2,1,2),plot(t,yi(1,:),"m")
subplot(2,1,2),plot(t,yi(3,:),"g")
subplot(2,1,2),plot(t,yi(5,:),"b")
xtitle('Ca [kmol/m^-3] vs Tempo [s]','t[s]','Ca[kmol/m^-3]')
legend(["Tj: 290 [K]","Tj: 298 [K]", "Tj: 310 [K]"],3)
xgrid(0)

//Temperatura
subplot(2,1,1),plot(t,yi(2,:),"m")
subplot(2,1,1),plot(t,yi(4,:),"g")
subplot(2,1,1),plot(t,yi(6,:),"b")
xtitle('Temperatura [K] vs Tempo [s]','t[s]','Toi [K]')
legend(["Tj: 290 [K]","Tj: 298 [K]", "Tj: 310 [K]"],2)
xgrid(0)








//---------------------------------------------------------------------------
//6. Avalie o comportamento do sistema caso após 4h ocorrer a interrupção na
//vazão de água de resfriamento na jaqueta, ou seja, Qrem=0. Discuta os resultados.
//---------------------------------------------------------------------------

//Condições iniciais
Cai = Cao
Tri = To

function f=din5(t,j)
    
    Ca = j(1)   
    Tr = j(2)
    
    //4 horas
    if t<=14400 then
        f(1)= (Cao/tau) - (Ca/tau) - ((ko * exp(-Ea/(R*Tr))) * Ca)
        f(2) = ((To-Tr)/tau) - (((ko * exp(-Ea/(R*Tr)))* Ca * deltaH)/pcP) - (Ua_V * (Tr -Tj)/pcP)
    else
        f(1)= (Cao/tau) - (Ca/tau) - ((ko * exp(-Ea/(R*Tr))) * Ca)
        f(2) = ((To-Tr)/tau) - (((ko * exp(-Ea/(R*Tr)))* Ca * deltaH)/pcP)
    end
    f=real(f)

endfunction

//10 Horas
t=0:0.1:10*3600
xinicial=[Cai;Tri] 
j=ode(xinicial,t(1),t,din5)

// ----------------------------------------------------------------------------
// --------------------------------Gráficos
// ----------------------------------------------------------------------------
scf()
clf()

//Concentração

subplot(2,1,1),plot(t,j(1,:),"m")
xtitle('Ca [kmol/m^-3] vs Tempo [s]','t[s]','Ci[kmol/m^-3]')
legend(["Ca [kmol/m^-3]"])
xgrid(0)

//Temperatura

subplot(2,1,2),plot(t,j(2,:),"b")
xtitle('Tr [K] vs Tempo [s]','t[s]','Tr [K]')
legend(["Tr [K]"],4)
xgrid(0)









