import numpy as np
import matplotlib.pyplot as plt

# Читаем данные
try:
    data = np.loadtxt("data/maxwell.txt", skiprows=1)
except FileNotFoundError:
    print("Error: File data/maxwell.txt not found.")
    exit()

# Разделяем данные на столбцы
velocity = data[:, 0]
sim_density = data[:, 1]
theory_density = data[:, 2]

# Вычисляем ширину столбца
bin_width = velocity[1] - velocity[0]

plt.figure(figsize=(8, 6))

# Рисуем столбики симуляции
plt.bar(velocity, sim_density, width=bin_width, alpha=0.6, color='skyblue', edgecolor='black', label='Simulation')

# Рисуем теоретическую кривую
plt.plot(velocity, theory_density, 'r-o', linewidth=2, markersize=4, label='Maxwell Theory')

plt.title('Maxwell-Boltzmann Velocity Distribution')
plt.xlabel('Velocity vx (LJ units)')
plt.ylabel('Probability Density')
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()

plt.savefig("data/plots/maxwell_plot.png", dpi=300)