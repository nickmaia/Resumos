import numpy as np
import matplotlib.pyplot as plt

class Esfera:
    def __init__(self, raio, discretizacao):
        self.raio = raio
        self.delta = discretizacao
        self.pontos = []
        self.erro = 10**-9
        for r in np.arange(0, self.raio + self.erro, self.delta):
            self.pontos.append([round(r, 3)])

    def get_points(self):
        return self.pontos

    def edge_degrees(self, superficie=None):
        for ponto in self.pontos:
            if ponto[0] == self.raio:
                ponto.append(superficie)
        return self.pontos

    def distancia(self, p1, p2):
        r1 = p1[0]
        r2 = p2[0]
        dr = round(abs(r1 - r2), 3)
        if dr == self.delta:
            return True
        else:
            return False

    def descovering_degrees(self):
        # DEFINING POINTS
        self.valide_points = []
        self.result_points = []
        self.edge_points = []
        for p in self.pontos:
            if len(p) == 1:
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
            a[ponto[1]] = -2
            for p in self.valide_points:
                if self.distancia(ponto, p):
                    a[p[1]] = 1
            matrix_A.append(a)

        vector_b = []
        for ponto in self.valide_points:
            b = 0
            for p in self.edge_points:
                if self.distancia(ponto, p):
                    b -= p[1]
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

    def plot(self):
        r = [p[0] for p in self.result]
        T = [p[1] for p in self.result]
        plt.figure(figsize=(10, 5))
        plt.plot(r, T)
        plt.xlabel('Raio (m)')
        plt.ylabel('Temperatura (°C)')
        plt.title('Distribuição de Temperatura na Esfera')
        plt.grid(True)
        plt.show()


# Creating an instance of the Sphere class
esfera = Esfera(1, 0.05)
esfera.edge_degrees(superficie=100)
esfera.descovering_degrees()

# Plotting the temperature distribution
esfera.plot()
