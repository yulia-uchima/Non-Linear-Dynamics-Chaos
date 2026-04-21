import numpy as np
import matplotlib.pyplot as plt
from plotUtils import plot_runge_kutta,calculate_error,f
from Classes import RKG


if __name__ == '__main__':
    '''
    This part of code is made to solve a differential equation using the classes.py with RKG.

    df_dx= x-y  with x_0=0, y_0=0 X_f=1 and h=0.00001

    
    '''

    
    x_0 = 0
    y_0 = np.array([0])
    x_f = 1
    h = 0.0001
    order = 2
    plot_runge_kutta(f, x_0, y_0, x_f, h,"./SolutionDifferentialEquation.png")


    # List of values of h for which we are going to calculate error
    h_vals = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005]
    errores = {order: [calculate_error(h, order) for h in h_vals] for order in [1, 2, 4]}

    # Convergence graph
    plt.figure(figsize=(10, 6))
    for order, error_vals in errores.items():
        plt.plot(h_vals, error_vals, marker='o', label=f"RKG Order {order}")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Step h (log)')
    plt.ylabel('Error (log)')
    plt.title('Convergence of Method RKG for $dy/dx = x - y$ with $y(0) = 0$')
    plt.legend()
    plt.grid(True)
    plt.savefig("./ConvergenceGraph.png")
    
