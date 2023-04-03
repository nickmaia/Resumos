clc
clear 
mode(-1)

//            Funções

function f=din(t,T,deltax,q,k,alfa,TN_1,N)
    To=T(2)+2*deltax*q/k
    f(1)=alfa*(T(2)-2*T(1)+To)/deltax^2
    
    for i=2:N-1
        f(i)=alfa*(T(i+1)-2*T(i)+T(i-1))/deltax^2
    end
    f(N)=alfa*(TN_1-2*T(N)+T(N-1))/deltax^2
    
endfunction

//                Entrada de dados
Ti=70
q=300
k=1
alfa=0.04
TN_1=Ti
N=50
xo=0
xf=5
deltax=(xf-xo)/N
x=0:deltax:xf

//                Processamento

xinicial=Ti*ones(N,1)
t=0:0.1:20
lista=list(din,deltax,q,k,alfa,TN_1,N)
T=ode(xinicial,t(1),t,lista)


//                Saída de dados

scf(0)
clf()
plot(t,T(1,:),'r')
plot(t,T(25,:),'b')
plot(t,T($,:),'m')
xtitle('','$t[h]$','$T[^oF]$')
legend(['N=1','N=25','N=50'])
xgrid(0)

//  Vetor Espacial
TT = [T;Ti*ones(1,length(t))]

scf(1)
clf()
plot(x,TT(:,11),'r')
plot(x,TT(:,111),'b')
plot(x,TT(:,$),'m')
xtitle('','$x[ft]$','$T[^oF]$')
legend(['t=1 min','t=10 min','t=20 min'])
xgrid(0)



