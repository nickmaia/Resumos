from math import exp

# Definindo os valores das variáveis
delta_cp = -7
tau = 0.1228
tr = 68
cpi = 403.445
to = 75
delta_hr = -36400
k1 = 10
ea = 32400
R = 1.987
t1 = 68

# Definindo a função f(t) a ser resolvida
def f(t):
    return (delta_cp*tau*tr - cpi*tau*to - delta_hr*tau) / (-delta_cp*tau + cpi*(k1*exp((ea/R) * (1/t1 - 1/t))) + cpi*tau) - t

# Definindo os limites inferior e superior do intervalo onde a raiz está
a = 0.1
b = 1000.0

# Definindo a tolerância de erro
tol = 1e-6

# Iniciando o processo de bissecção
while (b - a) / 2 > tol:
    c = (a + b) / 2  # selecionando o ponto médio do intervalo
    if f(c) == 0:
        break  # encontrou a raiz exatamente
    elif f(c) * f(a) < 0:  # selecionando o subintervalo apropriado
        b = c
    else:
        a = c

# Imprimindo o valor de "t"
t = (a + b) / 2
print(t)
