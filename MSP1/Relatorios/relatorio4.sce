clc
clear 
mode(-1)
//Nomes: Nicole Maia, Carina Aquino, Iara Solimar



// Funções

function f=difusao(t, CA, deltay, CA0, DAB, N)
    f(1) = DAB*(CA(2)-2*CA(1)+CA0)/deltay^2
    
    for i=2:N-1
        f(i) = DAB*(CA(i+1)-2*CA(i)+CA(i-1))/deltay^2
    end
    
    f(N) = DAB*(2*CA(N-1)-2*CA(N))/deltay^2
endfunction

// Entrada de dados
CA0 = 0.01
DAB = 2e-6
N = 50
yo = 0
yf = 0.1
deltay = (yf-yo)/N
y = yo:deltay:yf

// Processamento
CAinicial = zeros(N,1)
t = 0:1:5000
lista = list(difusao, deltay, CA0, DAB, N)
CA = ode(CAinicial, t(1), t, lista)

// Saída de dados

scf(2)
clf()
plot(t,CA(1,:),'r')
plot(t,CA(25,:),'b')
plot(t,CA($,:),'m')
xtitle('','$t[s]$','$C_A[mol/m^3]$')
legend(['y=0[m]','y=0.05[m]','y=0.1[m]'],4)
xgrid(0)

// Vetor Espacial
CAA = [CA;CA0*ones(1,length(t))]
scf(3)
clf()
plot(y,CAA(:,801),'r')
plot(y,CAA(:,3001),'b')
plot(y,CAA(:,$),'m')
xtitle('','$y[m]$','$C_A[mol/m^3]$')
legend(['t=800 s','t=3000 s','t=5000 s'],3)
xgrid(0)
