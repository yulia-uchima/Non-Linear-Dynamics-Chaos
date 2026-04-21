
import numpy as np
import matplotlib.pyplot as plt

"""
I create two classes to solve differential equations using the generalized Runge-Kutta method.
The first class solves the problem for a single differential equation, while the second finds
the solutions to a system of differential equations.
"""

#=====================================================================================================
class RKG:
    def __init__(self,func: callable,t_0:float, t_f:float,y_0:float,h:float,order) -> None:
        '''

        Params:
            func: Is the funtion that define the ODE. In explicit is f(x,y)= dy/dx

                      def func(x:float ,y:float ,params:list):
                        return dy/dx
            
            y_0: Initial value of the function that is sln for thr EDO.
            t_0: Initial value for the independent variable.                    

            t_f: Is the final point for the independent variable.

            h: Step size 

            order: Is 1 or 2 or 4 depending on the order of approach

        '''

        ''' Methods
              -------
            ks : array
                  returns an array with the ks at a point (t,y)
            soluciones : array
                  returns the solutions at each point tn
        '''


        self.func = func
        self.t0 =t_0
        self.y0 =y_0
        self.h = h
        self.t = np.arange(self.t0,t_f+h,h)# array for store data for the variable t
        self.y = np.zeros((len(y_0), len(self.t))) #Matrix for store solutions (n_equations, n_steps)

       
        
        #-------------------------------------------------------------------------------
        self.t[0] = self.t0   # set initial value for independent varible
        self.y[:,0] = self.y0 # set initial value for dependent varible
        self.solved = False   # keep information about the system has been solved or not
        self.order = order
        
 # Selection of the coefficients according to the order of the coefficients.

        '''a,b,c: a is a matrix, b and c are vectos, and together form the Butcher Tableau.
            
            Depending on the method these have different inputs and values.
            https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge%E2%80%93Kutta_methods

            
            '''
        if order == 4:
            # Coefficients for Runge-Kutta 4
            self.c = np.array([0, 0.5, 0.5, 1])
            self.a = np.array([
                [0,   0,   0,   0],
                [0.5, 0,   0,   0],
                [0,   0.5, 0,   0],
                [0,   0,   1,   0]
            ])
            self.b = np.array([1/6, 1/3, 1/3, 1/6])

        elif order == 2:
            # Coefficients for Runge-Kutta of order 2 (Midpoint method)
            self.c = np.array([0, 0.5])
            self.a = np.array([
                [0,   0],
                [0.5, 0]
            ])
            self.b = np.array([0, 1])

        elif order == 1:
            # Coefficients for Runge-Kutta of order 1 (Euler's method)
            self.c = np.array([0])
            self.a = np.array([[0]])
            self.b = np.array([1])

        else:
            raise ValueError("The order of the method: 1, 2 o 4.")



    def ks(self,i):
        
        '''
        Calculate the slopes k for each step, in the generalized Runge-Kutta method we must
        find K1,...,Ks for each step in the iteration.

        k1=f( tn, yn)
        k2=f( tn + c_2 h  , yn + (a_12 k_1) h )
        k3=f( tn + c_3 h  , yn + ( a_31 K_1  + a_32 K_2) h ) 
        .
        .
        .
        ks=f(tn + c_s h  , yn + ( a_s1 K_1  + a_s2 K_2+ ...+ a_s(s-1) K_(s-1) ) h ) 

        In this case, they are stoes in  nxs matrix, where:
        Each row represents all the k_i values for one dependent varible at i step.                                
        '''
        k = np.zeros((len(self.y0),len(self.b)))  #Matrix for store ks values
        
        for j in range(len(self.b)):
            k[:,j] =self.func(
                self.t[i] + self.h * self.c[j], 
                self.y[:, i] + self.h * np.dot(k[:, :j], self.a[j, :j])
            )
        
        return k
    
    
    def soluciones(self):
        '''
        We store all the solution of the n dependent variables in a nxN matrix, where n in the amount of dependent variables
        and N is the amount of steps from x0 to xf.
        
        Which means, each row is the solution for one dependent variable step by step.
        
        This methos calls the previous method to calculate the Ki values for each step and then use it here to calculate the solutions.
        '''
        for i in range(0,len(self.t)-1):
            k = self.ks(i)
            self.y[:,i+1] = self.y[:,i] + self.h*np.dot(k,self.b)

        
        self.solved = True  
        return self.t, self.y
    

#===================================================================================================


# Inheriting class to solve systems of coupled equations---------------------------------------------
class RKGA(RKG):
    def __init__(self, funcs: list, t_0: float, t_f: float, y_0: np.ndarray, h: float, orden: int):
        '''
        Params:
            funcs: List of functions defining the system of coupled ODEs.
            t_0: Initial value of the independent variable common to all equations.
            t_f: Final value of the independent variable, comon to all equations.
            y_0: Initial values of the system (a vector with the initial conditions for all the functions)
            h: Step size
            order: Runge-Kutta method order (1, 2 or 4).
        '''
        super().__init__(funcs, t_0, t_f, y_0, h, orden)

    def ks(self, i):
        '''
        Calculate the values of k1, i for each function of the system.
        '''
            

        num_eqs = len(self.func)  # Number of equations of the sytem
        k = np.zeros((num_eqs, len(self.b)))  # Matrix for store K values
        
        '''Matrix for store K values (n_equations, n_steps)
        For each equation in the system you have a set of data that depends on the number of steps.
        On the number of steps, each row of the matrix will store the data set associated to each function
        of the system .
        '''
        
        for j in range(len(self.b)):
            # Calculation of k_j for each function in the system
            for eq_index, func in enumerate(self.func):
                # Each function is evaluated with the current state of the system.
                k[eq_index, j] = func(
                    self.t[i] + self.h * self.c[j], 
                    self.y[:, i] + self.h * np.dot(k[:, :j], self.a[j, :j])
                )
        return k


