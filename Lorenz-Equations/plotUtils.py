
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


from Classes import RKG,RKGA

"""
I define the functions used in main1 and main to visualize the solutions of the differential equation 
systems, such as the exact solutions, the graphs,
and the functions that define the problems. 
"""

# Exact solution using sympy-----------------------------------------------------
x_sym = sp.Symbol('x')
y_sym = sp.Function('y')(x_sym)
ode = sp.Eq(y_sym.diff(x_sym), x_sym - y_sym)
sol_particular = sp.dsolve(ode, y_sym, ics={y_sym.subs(x_sym, 0): 0})
sol_particular_simplificada = sp.simplify(sol_particular.rhs)
sol_func_particular = sp.lambdify(x_sym, sol_particular_simplificada, "numpy")


# Funtion for Graph------------------------------------------------------------------

def plot_runge_kutta(f, x_0, y_0, x_f, h,ruta:str):
    orders = [1, 2, 4]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    # Exact Solution using sympy
    x_sym = sp.Symbol('x')
    y_sym = sp.Function('y')(x_sym)
    ode = sp.Eq(y_sym.diff(x_sym), f(x_sym, y_sym))
    sol_particular = sp.dsolve(ode, y_sym, ics={y_sym.subs(x_sym, 0): y_0[0]})
    sol_func_particular = sp.lambdify(x_sym, sol_particular.rhs, "numpy")
    
    # range for independet Variable
    x_vals = np.linspace(x_0, x_f, 100)
    y_vals_particular = sol_func_particular(x_vals)

    for i, order in enumerate(orders):
        rk_solver = RKG(f, x_0, x_f, y_0, h, order)
        x_rk, y_rk = rk_solver.soluciones()
        
        axs[i].plot(x_vals, y_vals_particular, label="Exact Solution (sympy)", color='olive')
        axs[i].plot(x_rk, y_rk[0], label=f"Runge-Kutta Order {order}", color='darkolivegreen', linestyle='--')
        axs[i].set_xlabel('x')
        axs[i].set_ylabel('y(x)')
        axs[i].set_title(f'Order {order}')
        axs[i].legend()
        axs[i].grid()

    plt.tight_layout()
    plt.savefig(ruta)

# Differential Equation------------------------------------------------------
def f(x, y):
    return x - y


def calculate_error(h, order=2):
    x_0 = 0
    y_0 = np.array([0])
    x_f = 1

    # Create the solver and obtain the numerical solution using RKG
    rk_solver = RKG(f, x_0, x_f, y_0, h, order)
    x_, y_ = rk_solver.soluciones()

    # Solución exacta usando sympy
    sol_particular_simplificada = sp.simplify(sp.dsolve(sp.Eq(y_sym.diff(x_sym), x_sym - y_sym), y_sym, ics={y_sym.subs(x_sym, 0): 0}).rhs)
    sol_func_particular = sp.lambdify(x_sym, sol_particular_simplificada, "numpy")

    # Evaluate the exact solution at the points x_
    y_exacta = sol_func_particular(x_)

    # Calculate the error (standard of the difference between the exact and numerical solution).
    error = np.linalg.norm(y_exacta - y_[0], ord=2)  # Norma L2 del error
    return error


#System of differential equations----------------------------------------------------------
def lorenz_x(t, y, sigma=10):
    return sigma * (y[1] - y[0])

def lorenz_y(t, y, rho=28):
    return y[0] * (rho - y[2]) - y[1]

def lorenz_z(t, y, beta=8/3):
    return y[0] * y[1] - beta * y[2]

# Graph for solutions-------------------------------------------------------
def graph_lorenz_solutions(t, y, condiciones_iniciales, parameters,ruta:str):
    fig = plt.figure(figsize=(10,10))
    fig.suptitle("Lorenz Equations", size=19)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4, projection='3d')

    # Gráficas
    ax1.plot(y[0], y[1], color='goldenrod')
    ax1.set_xlabel('x(t)')
    ax1.set_ylabel('y(t)')

    ax2.plot(y[0], y[2], color='olivedrab')
    ax2.set_xlabel('x(t)')
    ax2.set_ylabel('z(t)')

    ax3.plot(y[1], y[2], color='olive')
    ax3.set_xlabel('y(t)')
    ax3.set_ylabel('z(t)')

    ax4.plot(y[0], y[1], y[2], color='darkolivegreen')
    ax4.set_xlabel('x(t)')
    ax4.set_ylabel('y(t)')
    ax4.set_zlabel('z(t)')

    # Add initial conditions and parameters to the figure caption
    fig.text(0.5, 0.06, condiciones_iniciales, ha='center', fontsize=12)
    fig.text(0.5, 0.04, parameters, ha='center', fontsize=12)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(ruta)

# Function to solve the system with different initial conditions and parameters-------------------
def solve_lorenz_system(y0, sigma, rho, beta,ruta:str, t0=0, tf=100, h=0.01, orden=4):
    
    def lorenz_x(t, y):
        return sigma * (y[1] - y[0])

    def lorenz_y(t, y):
        return y[0] * (rho - y[2]) - y[1]

    def lorenz_z(t, y):
        return y[0] * y[1] - beta * y[2]

    # Initializing the RKGA solver for the Lorenz system                  <<<<<================== RKGA
    rk_solver = RKGA([lorenz_x, lorenz_y, lorenz_z], t0, tf, y0, h, orden)

    #Solutions
    t, y = rk_solver.soluciones()

    # Captions
    condiciones_iniciales = f"Initial Conditions: x0={y0[0]}, y0={y0[1]}, z0={y0[2]}"
    parameters = f"Parameters: σ={sigma}, ρ={rho}, β={beta}"
    
    # Graph results
    graph_lorenz_solutions(t, y, condiciones_iniciales, parameters,ruta)






