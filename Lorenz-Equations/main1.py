
import numpy as np

from plotUtils import solve_lorenz_system


if __name__ == '__main__':
    '''
    Solve the lorenz problem using the classes.py with RKGA.
    '''

    # Solution for Case 1-------------------------------------------------------------------------------
    y0 = np.array([0.0, 0.0, 0.0])  # Initial conditions
    sigma = 10
    rho = 28
    beta = 8/3
    solve_lorenz_system(y0, sigma, rho, beta,ruta="./SolutionCase1.png")


    # Solution for Case 2
    y0 = np.array([0.0, 0.01, 0.0])  # Initial conditions
    sigma = 10
    rho = 28
    beta = 8/3
    solve_lorenz_system(y0, sigma, rho, beta,ruta="./SolutionCase2.png")

    # Solution for Case 3------------------------------------------------------------
    y0 = np.array([0.0, 0.0, 0.0])  # Initial conditions
    sigma = 16
    rho = 45
    beta = 4
    solve_lorenz_system(y0, sigma, rho, beta,ruta="./SolutionCase3.png")

    # Solution for Case 4-------------------------------------------------------------------------
    y0 = np.array([0.0, 0.001, 0.0])  # Initial conditions
    sigma = 16
    rho = 45
    beta = 4
    solve_lorenz_system(y0, sigma, rho, beta,ruta="./SolutionCase4.png")