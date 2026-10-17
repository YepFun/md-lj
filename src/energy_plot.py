import numpy as np
import matplotlib.pyplot as plt

# Читаем данные
try:
    data = np.loadtxt("data/diag.txt", skiprows=1)
except FileNotFoundError:
    print("Error: File data/diag.txt was not found.")
    exit()

# Разделяем данные на столбцы
steps = data[:, 0]
K = data[:, 1]
U = data[:, 2]
E = data[:, 3]

plt.figure(figsize=(8, 6))

# Рисуем три линии
plt.plot(steps, K, label='Kinetic Energy (K)', color='red')
plt.plot(steps, U, label='Potential Energy (U)', color='blue')
plt.plot(steps, E, label='Total Energy (E = K + U)', color='black')

plt.title('Energy Conservation')
plt.xlabel('Step')
plt.ylabel('Energy (LJ units)')
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()

plt.savefig("data/plots/energy_plot.png", dpi=300)