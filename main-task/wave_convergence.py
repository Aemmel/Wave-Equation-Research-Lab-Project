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
def v_convergence(h_x, h_t, x_start, x_end, t_start, t_end, bc, stept):
    h_xc = [h_x, h_x/2, h_x/4]
    h_tc = [h_t, h_t/2, h_t/4]

    N_xc = [int((x_end - x_start)/h_xc[0] + 1),int((x_end - x_start)/h_xc[1] + 1),int((x_end - x_start)/h_xc[2] + 1)]
    N_tc = [int((t_end - t_start)/h_tc[0] + 1),int((t_end - t_start)/h_tc[0] + 1),int((t_end - t_start)/h_tc[0] + 1)]

    xc = [np.linspace(x_start, x_end, N_xc[0]),np.linspace(x_start, x_end, N_xc[1]),np.linspace(x_start, x_end, N_xc[2])]
    tc = [np.linspace(t_start, t_end, N_tc[0]),np.linspace(t_start, t_end, N_tc[1]),np.linspace(t_start, t_end, N_tc[2])]

    yout1, t1 = waveequation(xc[0], tc[0], h_xc[0], h_tc[0], bc)
    yout2, t2 = waveequation(xc[1], tc[1], h_xc[1], h_tc[1], bc)
    yout4, t4 = waveequation(xc[2], tc[2], h_xc[2], h_tc[2], bc)

    Z1, N1 = difference(yout1[stept,0,:], yout2[stept,0,::2], yout4[stept,0,::4])
    #Z2, N2 = maxnorm(yout1[2,0,:], yout2[2,0,::2], yout4[2,0,::4])

    return Z1, N1, yout1, stept
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Globale parameter
#start and end point of x
x_start = -1
x_end = 1

#start and end of t
t_start = 0
t_end = 2

#step scize Delta x and Delta t
h_x = 0.1
h_t = 0.01

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

#calculate time convergence
#stept = insert t at which you want to view the time convergence
stept = 200
Z1, N1, yout1, stept = v_convergence(h_x, h_t, x_start, x_end, t_start, t_end, 'periodic', stept)




plt.figure()
plt.plot(x, yout1[stept,0,:], color = 'black', label = 't = ' )
plt.xlabel('x')
plt.ylabel('$\\phi(x) $')
plt.legend(loc='upper right')
plt.grid()

plt.figure()
plt.plot(x,Z1, color = 'darkgreen', linewidth=2, label = '$\\Delta \\phi_{1} = \\phi(\\Delta(x)) - \\phi(\\Delta(x/2)) $')
plt.plot(x,N1, color = '#1f77b4', linewidth=2, label = '$\\Delta \\phi_{2} = \\phi(\\Delta(x/2)) - \\phi(\\Delta(x/4)) $')
plt.plot(x,2*N1, ':', color = 'darkorange', linewidth=3, label = 'Q = 2 $\\cdot \\Delta \\phi_2$ ')
plt.legend(loc='upper right')
plt.xlabel('x')
#plt.ylim(0, 0.01)
plt.yscale("log")
plt.ylabel('$\\Delta\\phi_{1,2}$')
plt.grid()

plt.show()
