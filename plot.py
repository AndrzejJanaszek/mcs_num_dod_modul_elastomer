import numpy as np
import matplotlib.pyplot as plt

# Współczynniki wielomianu
coeffs = [-0.000000000106297, 0.735922667334247, 0.122236351715270, -0.006190810898055, 0.000125815278125, -0.000001195081766, 0.000000004219488]

# Tworzenie funkcji wielomianowej
def polynomial(x, coeffs):
    return sum(c * x**i for i, c in enumerate(coeffs))

# Zakres wartości x
x_values = np.linspace(0, 60, 200)

# Obliczanie wartości wielomianu
y_values = polynomial(x_values, coeffs)

# Tworzenie wykresu
plt.plot(x_values, y_values, label="Wielomian")
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.xlabel("odkształcenie [%]")
plt.ylabel("naprężenie [MPa]")
plt.title("Wykres wielomianu")
plt.show()