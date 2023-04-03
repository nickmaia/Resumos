//------------------------------------------------------------------------------
//                                DESCRIÇÃO
//------------------------------------------------------------------------------
//Adimensional
//------------------------------------------------------------------------------
//                             PRÉ-PROCESSAMENTO
//------------------------------------------------------------------------------
clc
clear
mode(-1)

//------------------------------------------------------------------------------
//                             DADOS
//------------------------------------------------------------------------------

Cao = 1 //mol/L
k = 1.3 //1/s
Cai = Cao
Cbi = 0

//------------------------------------------------------------------------------
//                             FUNÇÕES
//------------------------------------------------------------------------------
//Modelo dimensional
function f=ee(x)
    Ca=x(1)
    Cb=x(2)
   
    f(1)= -k * Ca
    f(2)= k * Ca
endfunction

function dxdt=din(t,x)
    Ca=x(1)
    Cb=x(2)
    dxdt(1)= -k * Ca
    dxdt(2)= k * Ca
endfunction


//Adimensional
function f = ee_adimensional(x)
    Ca_ad = x(1)
    Cb_ad = x(2)
    Ca_star = Ca_ad/Cao
    Cb_star = Cb_ad/Cao
    f(1) = -Ca_star
    f(2) = Ca_star
endfunction


function dxdt = din_adimensional(t, x)
    Ca_ad = x(1)
    Cb_ad = x(2)
    Ca_star = Ca_ad/Cao
    Cb_star = Cb_ad/Cao
    dxdt(1) = -Ca_star
    dxdt(2) = Ca_star
    
endfunction

//------------------------------------------------------------------------------
//                              PROGRAMA PRINCIPAL
//------------------------------------------------------------------------------



//----------->Processamento

//Modelo dimensional
xi=[Cai; Cbi]
[sol,v,info]=fsolve(xi, ee)

if info==1 then
     mprintf(' Solução estado estacionário dimensionalizado: \n\n')
    mprintf(' Ca : %f\n', sol(1))
    mprintf(' Cb : %f\n', sol(2))
else
    disp('Tente novamente!')
end


xi=[Cai; Cbi]
t=[0:0.1:30]
xf = ode(xi, t(1), t, din)

xa = (Cao-xf(1,:))./Cao

//Modelo adimensional

xi=[Cai; Cbi]/Cao
[sol2,v,info]=fsolve(xi, ee_adimensional)

if info==1 then
     mprintf(' Solução estado estacionário dimensionalizado: \n\n')
    mprintf(' Ca : %f\n', sol2(1))
    mprintf(' Cb : %f\n', sol2(2))
else
    disp('Tente novamente!')
end


xi=[Cai; Cbi]/Cao
t=[0:0.1:30]
tau = t*k
xf2 = ode(xi, tau(1), tau, din_adimensional)

xa_ad = (Cao-xf2(1,:))./Cao
//----------->Saída de dados
scf(0)
clf()
// Dimensionais
subplot(2,2,1),plot(t,xf(1,:),'r')
subplot(2,2,1),plot(t,xf(2,:),'b')
xtitle("Dimensional",'t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"])
xgrid(0)

subplot(2,2,3),plot(t,xa,"m")
xtitle('Dimensional','t[s]','Xa')
xgrid(0)

// Adimensionais
subplot(2,2,2),plot(tau,xf2(1,:),'r')
subplot(2,2,2),plot(tau,xf2(2,:),'b')
xtitle('Adimensional','tau[-]','Ci*[-]')
legend(["Ca* [-]","Cb* [-]"])
xgrid(0)

subplot(2,2,4),plot(tau,xa_ad,"m")
xtitle('Adimensional','tau[-]','Xa')
xgrid(0)


//----------->Fim do programa
disp('***FIM***')



