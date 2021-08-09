import numpy as np
import matplotlib.pyplot as plt
from functions import *

#Funktion
def f(x, t):
    w = 2*np.pi
    k = 2*np.pi
    # np.exp(-k*x**2)
    return np.exp(-k*x**2)


def waveequation(x, t, h_x, h_t, bc):

    if bc == "periodic":
        phi_0 = np.zeros(len(x) + 2)    # Phi(x, t = 0)
        phi_0[1:-1] = f(x, 0)
        phi_0[0] = phi_0[-3]            # add Ghostpoints
        phi_0[-1] = phi_0[2]

        pi_0 = np.zeros(len(x)) 	   # Pi(x, t = 0)
        pi_0 = -1*first_derivative(phi_0, h_x)      # -1 für richtung der welle


    elif bc == "open_i":
        phi_0 = np.zeros(len(x) + 2)    # Phi(x, t = 0)
        phi_0[1:-1] = f(x, 0)
        phi_0[0] = 2*phi_0[1] - phi_0[2]            # add Ghostpoints
        phi_0[-1] = 2*phi_0[-2] - phi_0[-3]

        pi_0 = np.zeros(len(x)) 	   # Pi(x, t = 0)
        pi_0 = -1*first_derivative(phi_0, h_x)      # -1 für richtung der welle

    elif bc == "open_ii":
        phi_0 = np.zeros(len(x) + 2)    # Phi(x, t = 0)
        phi_0[1:-1] = f(x, 0)
        phi_0[0] = phi_0[1] - phi_0[2]            # add Ghostpoints
        phi_0[-1] = phi_0[-2] - phi_0[-3]

        pi_0 = np.zeros(len(x)) 	   # Pi(x, t = 0)
        pi_0 = -1*first_derivative(phi_0, h_x)      # -1 für richtung der welle

    solution = [phi_0[1:-1], pi_0]
#-------------------------------------------------------------------------------

    def dgl(t, solution):                  # PDE of waveequation + boundary conditions
        if bc == "periodic":
            phi = np.zeros(len(solution[0]) + 2)
            phi[1:-1] = solution[0]
            phi[0] = solution[0][-3]
            phi[-1] = solution[0][2]

        if bc == "open_i":
            phi = np.zeros(len(solution[0]) + 2)
            phi[1:-1] = solution[0]
            phi[0] = phi[1] - phi[2]
            phi[-1] = phi[-2] - phi[-3]

        if bc == "open_ii":
            phi = np.zeros(len(solution[0]) + 2)
            phi[1:-1] = solution[0]
            pi = solution[1]

            phi[0] = phi[2] - 2*h_x*pi[0]
            phi[-1] = phi[-3] - 2*h_x*pi[-1]

        dphidt = solution[1]
        dpidt =  second_derivative(phi, h_x)
        u_punkt = np.array([dphidt, dpidt])
        return u_punkt

    yout, t = rungekutta(solution, t, dgl)

    return yout, t

#------------------------------------------------------------------------------
#convergence maximumsnorm Z = phi(x) - phi(x/2)  und N = phi(x/2) - phi(x/4)

def time_convergence(h_x, h_t, x_start, x_end, t_start, t_end, bc):
    h_xc = [h_x, h_x/2, h_x/4]
    h_tc = [h_t, h_t/2, h_t/4]

    N_xc = [int((x_end - x_start)/h_xc[0] + 1),int((x_end - x_start)/h_xc[1] + 1),int((x_end - x_start)/h_xc[2] + 1)]
    N_tc = [int((t_end - t_start)/h_tc[0] + 1),int((t_end - t_start)/h_tc[1] + 1),int((t_end - t_start)/h_tc[2] + 1)]

    xc = [np.linspace(x_start, x_end, N_xc[0]),np.linspace(x_start, x_end, N_xc[1]),np.linspace(x_start, x_end, N_xc[2])]
    tc = [np.linspace(t_start, t_end, N_tc[0]),np.linspace(t_start, t_end, N_tc[1]),np.linspace(t_start, t_end, N_tc[2])]

    yout1, t1 = waveequation(xc[0], tc[0], h_xc[0], h_tc[0], bc)
    yout2, t2 = waveequation(xc[1], tc[1], h_xc[1], h_tc[1], bc)

    #here the convergence for every period T is calculated
    #if you change the time in Global parameters you have to change the points here
    Q2 = convergence(yout1[1,0,:], yout1[2000,0,:], yout2[4000,0,::2])
    Q3 = convergence(yout1[1,0,:], yout1[4000,0,:], yout2[8000,0,::2])
    Q4 = convergence(yout1[1,0,:], yout1[6000,0,:], yout2[12000,0,::2])
    Q5 = convergence(yout1[1,0,:], yout1[8000,0,:], yout2[16000,0,::2])
    Q6 = convergence(yout1[1,0,:], yout1[10000,0,:], yout2[20000,0,::2])

    return  Q2, Q3, Q4, Q5, Q6
#-----------------------------------------------------------------------------
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#---------------------------------------------------------------------
#Globale parameter
#start and end point of x
x_start = -1
x_end = 1

#start and end of t
#if you change timeperiods you have to calculate by hand the right timesteps for full periods and write them in the time_convergence function
t_start = 0
t_end = 10

#step scize Delta x and Delta t
h_x = 0.01
h_t = 0.001

#velocity
c = 1

alpha = c*h_t/h_x
print("alpha = ", alpha)

#define spacial grdid and time grid
N_x = int((x_end - x_start)/h_x + 1)
N_t = int((t_end - t_start)/h_t + 1)

x = np.linspace(x_start, x_end, N_x)
t = np.linspace(t_start, t_end, N_t)

#---------------------------------------------------------------------
#print the wave: choose the boundary conditions between 'periodic', 'open_i' and 'open_ii'
yout, t = waveequation(x, t, h_x, h_t, 'periodic')

#Q2, Q3, Q4, Q5, Q6 = convergence at T= 2, 4, 6, 8, 10
#if you change timeperiods you have to calculate by hand the right timesteps for full periods and write them in the time_convergence function
Q2, Q3, Q4, Q5, Q6 = time_convergence(h_x, h_t, x_start, x_end, t_start, t_end, "periodic")


#to plot it
x_achse = [1,2,3,4,5]
Q = [ Q2, Q3, Q4, Q5, Q6]

plt.figure()
plt.pcolor(x, t, yout[:,0], label = '$\\alpha$ = %1.3f' %alpha)
plt.ylabel('t')
plt.xlabel('x')
plt.colorbar()
plt.title('$\\phi(x, t)$')
plt.legend()


plt.figure()
plt.plot(x_achse, Q, 'bo', label = 'Q(T)')
plt.ylim(0,3)
plt.xlabel('period T')
plt.ylabel('Q')
plt.legend()
plt.grid()
plt.ylim(1.75,2.4)

plt.show()
