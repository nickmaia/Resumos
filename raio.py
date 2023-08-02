import time
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


class Cilindro:
    def __init__(self, raio, altura, discretizacao) -> None:
        self.raio = raio
        self.altura = altura
        self.delta = discretizacao
        self.pontos = []
        self.erro = 10**-9
        for r in np.arange(0, self.raio + self.erro, self.delta):
            for z in np.arange(0, self.altura + self.erro, self.delta):
                self.pontos.append([round(r, 3), round(z, 3)])

    def get_points(self):
        return self.pontos

    def edge_degrees(self, lateral=None):
        for ponto in self.pontos:
            if ponto[0] == self.raio:
                ponto.append(lateral)
        return self.pontos

    def distancia(self, p1, p2):
        r1 = p1[0]
        r2 = p2[0]
        z1 = p1[1]
        z2 = p2[1]
        dr = round(abs(r1 - r2), 3)
        dz = round(abs(z1 - z2), 3)
        if dr == self.delta and dz == 0:
            return True
        elif dr == 0 and dz == self.delta:
            return True
        else:
            return False

    def descovering_degrees(self):
        # DEFINING POINTS
        self.valide_points = []
        self.result_points = []
        self.edge_points = []
        for p in self.pontos:
            if len(p) == 2:
                self.valide_points.append(p.copy())
                self.result_points.append(p.copy())
            else:
                self.edge_points.append(p)

        num = 0
        for p in self.valide_points:
            p.append(num)
            num += 1

        # DEFINING MATRIX AND VECTOR
        matrix_A = []
        for ponto in self.valide_points:
            a = [0 for i in range(len(self.valide_points))]
            a[ponto[2]] = 1
            for p in self.valide_points:
                if self.distancia(ponto, p):
                    a[p[2]] += -1 / 4
            matrix_A.append(a)

        vector_b = []
        for ponto in self.valide_points:
            b = 0
            for p in self.edge_points:
                if self.distancia(ponto, p):
                    b += p[2] / 4
            vector_b.append(b)

        # SAVING RESULTS
        self.matrix_A = matrix_A
        self.vector_b = vector_b

        # CALCULATING THE ANSWERS
        matrix_A = np.array(matrix_A)
        vector_b = np.array(vector_b)
        inversa_A = np.linalg.inv(matrix_A)
        temperaturas = np.dot(inversa_A, vector_b)
        for i, T in enumerate(temperaturas):
            self.result_points[i].append(T)
        self.result = self.result_points + self.edge_points


# Create a cylinder and calculate the temperature distribution
cilindro = Cilindro(0.5, 0.5, 0.05)
cilindro.edge_degrees(lateral=500)
cilindro.descovering_degrees()

# Recreate the temperature distribution in cylindrical coordinates
r = np.array([r[0] for r in cilindro.result])
z = np.array([z[1] for z in cilindro.result])
T = np.array([T[2] for T in cilindro.result])

# Create a grid in the theta direction (from 0 to 2pi)
theta = np.linspace(0, 2*np.pi, 100)

# Use meshgrid to create 2D arrays of r, theta, z
r_2D, theta_2D = np.meshgrid(r, theta)
T_2D, _ = np.meshgrid(T, theta)

# Create the figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, polar=True)

# Plot the surface
c = ax.pcolormesh(theta_2D, r_2D, T_2D, cmap='hot')
fig.colorbar(c, ax=ax)

plt.show()
