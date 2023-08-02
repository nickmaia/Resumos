import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

class CilindroComprimento:
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

    def edge_degrees(self, lateral=None, topo=None, base=None):
        for ponto in self.pontos:
            if ponto[0] == self.raio or ponto[1] == 0 or ponto[1] == self.altura:
                if ponto[0] == self.raio:
                    ponto.append(lateral)
                elif ponto[1] == 0:
                    ponto.append(base)
                elif ponto[1] == self.altura:
                    ponto.append(topo)
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

    def plot_points(self, c):
        cores = [
            plt.cm.jet,
            plt.cm.hot,
            plt.cm.plasma,
            plt.cm.magma,
            plt.cm.OrRd,
        ]
        self.descovering_degrees()
        r, z, T = (
            [r[0] for r in self.result],
            [z[1] for z in self.result],
            [T[2] for T in self.result],
        )
        if c >= len(cores) or c < 0:
            print("Escolha uma cor entre 0 e 4")
            return None
        print(f"A cor escolhida foi {cores[c].name}")

        color = cores[c]
        fig = plt.figure(figsize=(16, 12))
        for i, a in enumerate(range(30, 121, 30), start=1):
            ax = plt.subplot(2, 2, i, projection="3d")
            surf = ax.plot_trisurf(r, z, T, cmap=color, linewidth=0.2)
            fig.colorbar(surf, shrink=0.5, aspect=5)
            ax.view_init(45, a)
            ax.set_title(f"Cilindro, angulo {a}")
            ax.set_xlabel("Raio do cilindro (m)")
            ax.set_ylabel("Comprimento do cilindro (m)")
            ax.set_zlabel("Temperatura do cilindro")

        plt.show()


# Implementing the class and calculating the temperature distribution
cilindro_comprimento = CilindroComprimento(0.5, 1, 0.05)
cilindro_comprimento.edge_degrees(lateral=100, topo=50, base=50)
cilindro_comprimento.descovering_degrees()

# Plotting the temperature distribution
cilindro_comprimento.plot_points(0)
