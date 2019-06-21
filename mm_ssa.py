from random import uniform
from math import log
import matplotlib.pyplot as plt
from time import sleep
import numpy as np

class MM_SSA():
    '''
    Gillespie SSA implementation for Michaelis Menten enzyme kinetics
    
    E + S -> C (k1)
    C -> E + S (k2)
    C -> E + P (k3)
    
    Input parameters:
        E0
        S0
        C0
        P0
        k1
        k2
        k3
        max_time (default set to 100)
        initial_time (default set to 0)
        
    Default usage:
        model=MM_SSA(1, 10, 1, max_time =20, E0 = 100, S0 = 10, C0 = 0, P0 = 0)
    '''
    __author__ = 'Stefan Alexandru Obada'
    __date__ = '09/06/2019'
    
    def __init__(self, k1, k2, k3, max_time = 100, E0 = 100, S0 = 100, C0 = 0, P0 = 0, initial_time = 0):
        self.E0 = E0
        self.S0 = S0
        self.C0 = C0
        self.P0 = P0
        self.initial_time = initial_time
        self.E = [E0]
        self.S = [S0]
        self.C = [C0]
        self.P = [P0]
        self.time = [initial_time]
        self.max_time = max_time
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.equilibrium = 0
    
    def SSA(self):
        while(self.time[-1] < self.max_time):
            self.step()
#        if self.equilibrium == 1:
#            print('Equilibrium state reached')
    
    def step(self):
        #Generate r1, r2
        r1 = uniform(0, 1)
        r2 = uniform(0, 1)
        
        #Calculate propensities and and their sum. ie alfa[]
        alfa = [0, 0, 0, 0]
        alfa[1] = self.k1 * self.E[-1] * self.S[-1]
        alfa[2] = self.k2 * self.C[-1]
        alfa[3] = self.k3 * self.C[-1]
        alfa[0] = alfa[1] + alfa[2] + alfa[3]
        
        #Set 'tau'
        if(alfa[0] > 0):
            tau = 1/alfa[0] * log(1 / r1)
            #Computing
            if ( r2 < (alfa[1]/alfa[0]) ):
                self.E.append(self.E[-1] - 1)
                self.S.append(self.S[-1] - 1)
                self.C.append(self.C[-1] + 1)
                self.P.append(self.P[-1])
                self.time.append(self.time[-1] + tau)
            elif ( r2 >= alfa[1]/alfa[0] and r2 < (alfa[1] + alfa[2]) / alfa[0] ):
                self.E.append(self.E[-1] + 1)
                self.S.append(self.S[-1] + 1)
                self.C.append(self.C[-1] - 1)
                self.P.append(self.P[-1])
                self.time.append(self.time[-1] + tau)
            elif ( r2 >= (alfa[1] + alfa[2])/alfa[0] and r2 < 1 ):
                self.E.append(self.E[-1] + 1)
                self.S.append(self.S[-1])
                self.C.append(self.C[-1] - 1)
                self.P.append(self.P[-1] + 1)
                self.time.append(self.time[-1] + tau)
        else:
            self.E.append(self.E[-1])
            self.S.append(self.S[-1])
            self.C.append(self.C[-1])
            self.P.append(self.P[-1])
            self.time.append(self.time[-1] + 1)
            self.equilibrium = 1
            
    def clear_SSA(self):
        self.E = [self.E0]
        self.S = [self.S0]
        self.C = [self.C0]
        self.P = [self.P0]
        self.time = [self.initial_time]
        
    def get_trajectory(self):
        'Returns a 4-dimensional numpy array'
        return np.array(list(zip(self.E, self.S, self.C, self.P)))
    
    def multi_plot(self, reactant, N_plots = 10, scaler = 2):
        '''Provide N_plots = number of plots (10 by default)
        scaler = scale of Ox axis(2 by default)
        reactant = {Choose one of E,S,C,P}. input as 'E' '''
        
        for i in range(N_plots):
            self.clear_SSA()
            self.SSA()
            _reactant = self._get_reactant(reactant)
            color = {'E' : 'red', 'S' : 'blue', 'C' : 'orange', 'P': 'magenta'}
            plt.plot(self.time[::scaler], _reactant[::scaler], linewidth = 0.7, c = color[reactant])
        
        plt.xlim((self.initial_time, self.max_time))
        plt.xlabel('TIME')
        plt.ylabel('Reactant')
        
        plt.grid()
        
    def _get_reactant(self, reactant):
        if reactant == 'E':
            return self.E
        elif reactant == 'S':
            return self.S
        elif reactant == 'C':
            return self.C
        elif reactant == 'P':
            return self.P
    
    #The following is just a test function for the form of distribution. INCOMPLETE
    '''def fi_E(self, n, N_max = 50):
        values = []
        for i in range(N_max):
            self.clear_SSA()
            self.SSA()
            values.append(len([x for x in self.E if x == n])/len(self.E))
        return(sum(values)/len(values))'''

                
    