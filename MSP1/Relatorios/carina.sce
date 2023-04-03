//------------------------------------------------------------------------------
//                                DESCRIÇÃO
//------------------------------------------------------------------------------
//Relatório 03
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
k = 1.3 //1/s
//------------------------------------------------------------------------------
//                             CONDIÇÕES INICIAIS
//------------------------------------------------------------------------------
Cai = Cao
Cbi = 0

//------------------------------------------------------------------
//                             Dimensional
//------------------------------------------------------------------

//-------------------------Estado estacionário
function f=ee(xi)
     
    Ca = xi(1)    
    Cb = xi(2)

    
    //Ca
    f(1)= -k * Ca
    //Cb
    f(2)= k * Ca

    f=real(f)
    
endfunction

//Determinando a solução para o estado Estacionário
xchute = [Cai;Cbi] 
[sol,v,info]=fsolve(xchute,ee)


if info==1 then
    mprintf(' Solução estado estacionário dimensionalizado: \n\n')
    mprintf(' Ca : %f\n', sol(1))
    mprintf(' Cb : %f\n', sol(2))
    
else 
    disp("Tente novamente!")
end


//-----------------------------Estado Dinâmico
function f=din(t,xi)   
      
    Ca = xi(1)    
    Cb = xi(2)

    
    //Ca
    f(1)= -k * Ca
    //Cb
    f(2)= k * Ca

    f=real(f)
    
    
endfunction

// Determinando resultados

t=0:0.01:15
xinicial=[Cai;Cbi] 
x=ode(xinicial,t(1),t,din)


//----------------- Conversão 𝑋𝐴
xa = (Cao-x(1,:))./Cao

//-------------------------------------------------------
//           Gráfico da concentração
//

scf()
clf()

subplot(1,1,1),plot(t,x(1,:),"m")
subplot(1,1,1),plot(t,x(2,:),"b")
xtitle("Dimensional",'t[s]','Ci[mol/L]')
legend(["Ca [mol/L]","Cb [mol/L]"])
xgrid(0)
//----------------------------------------------------------------------------
//                Gráfico da Conversão
//
scf()
clf()

subplot(1,1,1),plot(t,xa,"m")
xtitle('Dimensional','t[s]','Xa')
xgrid(0)

//------------------------------------------------------------------
//                          Adimensional
//------------------------------------------------------------------

//-------------------------Estado estacionário
function f=ee_ad(ji)
     
    Ca = ji(1)    
    Cb = ji(2)
    Ca_star = Ca/Cao
    Cb_star = Cb/Cao

    f(1)= -Ca_star
    f(2)= Ca_star
    //f=real(f)
    
endfunction

//Determinando a solução para o estado Estacionário

xchute_ad=xchute/Cao
[sol_ad,v_ad,info_ad]=fsolve(xchute_ad,ee_ad)


if info_ad==1 then
    mprintf('\n Solução estado estacionário adimensionalizado: \n\n')
    mprintf(' Ca Adimensional: %f\n', sol_ad(1))
    mprintf(' Cb Adimensional: %f\n', sol_ad(2))
    
else 
    disp("Tente novamente!")
end


//--------------Estado Dinâmico
function f=din_ad(t,ji)   
    Ca = ji(1)    
    Cb = ji(2)
    Ca_star = Ca/Cao
    Cb_star = Cb/Cao

    f(1)= -Ca_star
    f(2)= Ca_star
    
    //f=real(f)
    
endfunction

// Determinando resultados

tau=t*k
xinicial_ad=xinicial/Cao
x_ad=ode(xinicial_ad,tau(1),tau,din_ad)

//------------------------------------------------------
//            Conversão 𝑋𝐴
xa_ad = (Cao-x_ad(1,:))./Cao
//-------------------------------------------------------
//           Gráfico da concentração
//

scf()
clf()

subplot(1,1,1),plot(tau,x_ad(1,:),"m")
subplot(1,1,1),plot(tau,x_ad(2,:),"b")
xtitle('Adimensional','tau[-]','Ci*[-]')
legend(["Ca* [-]","Cb* [-]"])
xgrid(0)
//----------------------------------------------------------------------------
//                Gráfico da Conversão
//
scf()
clf()
subplot(1,1,1),plot(tau,xa_ad,"m")
xtitle('Adimensional','tau[-]','Xa')
xgrid(0)
