import numpy as np
import matplotlib.pyplot as plt

# Definición de parámetros para la ecuación de advección 1D
L = 1.0             # Longitud del dominio
Nx = 30             # Número de puntos en la malla
dx = L / Nx         # Tamaño de la malla
x = np.linspace(0, L, Nx)  # Puntos de la malla
T = 0.3             # Tiempo total de simulación
c = 1.0             # Velocidad de advección constante

# Solución analítica inicial (onda gaussiana)
def initial_condition(x):
    return np.exp(-100 * (x - 0.5)**2)

# Solución analítica
def analytical_solution(x, t):
    return np.exp(-100 * (x - 0.5 - c * t)**2)

# Configuraciones de los números de Courant
courant_numbers = [ 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1]
dt_values = [courant * dx / c for courant in courant_numbers]  # Calculamos dt para cada número de Courant
Nt_values = [int(T / dt) for dt in dt_values]  # Calculamos el número de pasos de tiempo

# Preparación para almacenar resultados
results = {}

# Esquema Upwind
for i, courant in enumerate(courant_numbers):
    dt = dt_values[i]
    Nt = Nt_values[i]
    
    # Inicializar la condición inicial
    u = initial_condition(x).copy()
    u_new = np.zeros_like(u)
    
    # Avanzar en el tiempo con el esquema Upwind
    for n in range(Nt):
        for j in range(1, Nx):  # j empieza en 1 para evitar acceso fuera de rango
            u_new[j] = u[j] - courant * (u[j] - u[j-1])
        u = u_new.copy()
    
    # Almacenar resultados
    results[courant] = u

# Solución analítica para comparar al tiempo final
u_analytical = analytical_solution(x, T)

# Calcular errores L2 para cada Courant
errors = {courant: np.sqrt(np.sum((results[courant] - u_analytical) ** 2) * dx) for courant in courant_numbers}

# Mostrar los resultados de error
print("Errores L2 para diferentes números de Courant:")
for courant, error in errors.items():
    print(f"Courant = {courant}: Error L2 = {error}")

# Graficar las soluciones
plt.figure(figsize=(10, 6))
plt.plot(x, u_analytical, label='Solución Analítica', linestyle='--', color='black')
for courant in courant_numbers:
    plt.plot(x, results[courant], label=f'Courant = {courant}')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Comparación de Soluciones Numéricas y Analítica')
plt.legend()
plt.grid(True)
plt.show()
