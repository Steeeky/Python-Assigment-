# Python-Assigment-
Final Assignment
# Python-Assigment-
Final Assignment
#Background
# 1.Briefly describe what Computational Physics is (A Brief Bakground)
# 2.Link it to the task you have been given
# 3.Outline your objectives for this study

          #SOLUTIONS
1.Computational physics is the study of scientific problems using computational methods.it combines computer science, physics and applied mathematics to develop scientific solutions to complex problems.Computational physics has become a very essential method in dealing with complex analysics or calculation in some areas of physics such as quantum physics and astrophysics. From last two decades it has also been very essencial tools in engineering(materials, mechanical, automobile).

2. The fundamental equation that describes how particles behave in a specific potential is known as the Schrodinger equation in the quantum realm. By employing a graphical representation, the Schrodinger equation can be solved in Python to provide a clear picture of the particle's behavior. One instance of a computational technique in physics is the Schrodinger equation. 

3. Get a solution for the wave function eigenvalue in a one-dimensional well using numerical methods.
   Examine the behavior of the wave function throughout the potential well and how it relates to the energy of the particle.


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

class ParticleInABox:
    def __init__(self, L=1e-9, N=1000):
        self.L = L  # length of the box
        self.N = N  # number of points to discretize the box
        self.dx = L / N  # step size
        self.x = np.linspace(0, L, N)  # position array

        self.hbar = 1.055e-34  # Planck's constant divided by 2π
        self.m = 9.11e-31  # mass of electron

        # Define the potential energy function
        self.V = lambda x: 0

        # Define the Hamiltonian matrix
        self.H = np.zeros((N, N))
        for i in range(1, N-1):
            self.H[i, i] = -2 / (self.dx ** 2) + self.V(self.x[i])
            self.H[i, i-1] = 1 / (self.dx ** 2)
            self.H[i, i+1] = 1 / (self.dx ** 2)
        self.H[0, 0] = -2 / (self.dx ** 2) + self.V(self.x[0])
        self.H[0, 1] = 1 / (self.dx ** 2)
        self.H[-1, -1] = -2 / (self.dx ** 2) + self.V(self.x[-1])
        self.H[-1, -2] = 1 / (self.dx ** 2)

    def solve(self):
        # Diagonalize the Hamiltonian matrix
        self.E, self.psi = eig(self.H)
        idx = self.E.argsort()  # sort the eigenvalues and eigenvectors
        self.E = self.E[idx]
        self.psi = self.psi[:, idx]

    def plot_wave_functions(self, num_functions=5):
        # Plot the wave functions
        fig, ax = plt.subplots()
        for i in range(num_functions):
            ax.plot(self.x, self.psi[:, i], label=f"n={i+1}")
        ax.set_xlabel("x")
        ax.set_ylabel("ψ(x)")
        ax.legend()

    def plot_energy_levels(self):
        # Plot the energy levels
        fig, ax = plt.subplots()
        ax.plot(range(1, len(self.E)+1), self.E, 'bo')
        ax.set_xlabel("n")
        ax.set_ylabel("Energy (J)")
        ax.set_title("Particle in a one-dimensional box")

        plt.show()

# Example usage
p = ParticleInABox()
p.solve()
p.plot_wave_functions()
p.plot_energy_levels()

#Conclusion
# 1. Give a brief conclusion.
# 2. Also, discuss your learning journey. Highlight your challenges and great moments.
 """"            SOLUTIONS
1. In order to determine the schrodinger equation of a particle. One needs to be familiar with both Python programming, including functions, loops, and lists, as well as quantum mechanics theory in order to utilize Python to solve this equation for a one-dimensional potential well. 

2. While learning Python's syntax and creating a few loops was fun, dealing with the data structure on its own was challenging. particularly if you haven't experimented with the kind of problem. However, looking back, it appears to be very simple and straightforward once you've been able to offer the right answers.   """"

Recomendations
The course should provide basic project in order for students to known how to handle a problem the way thye it should be handle.
There are data structure books at the college library which are easy to understand.
    
