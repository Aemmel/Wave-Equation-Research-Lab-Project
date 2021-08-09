import numpy as np
import matplotlib.pyplot as plt


def first_derivative(u, h):
    return (u[2:] - u[:-2])/(2*h)


def second_derivative(u, h):
    return (u[:-2] - 2*u[1:-1] + u[2:])/(h**2)


def convergence(ana, num1, num2):
    Z = (ana - num1)**2
    N = (ana - num2)**2
    Q = np.sqrt(sum(Z))/np.sqrt(sum(N))
    return Q


def selfconvergence(num1, num2, num3):
    Z = (num1 - num2)**2
    N = (num2 - num3)**2
    Q_self = np.sqrt(sum(Z))/np.sqrt(sum(N))
    return Q_self

def difference(num1, num2, num3):
    Z = abs(num1[:] - num2[:])
    N = abs(num2[:] - num3[:])
    return Z, N

#selfconvergence maxnorm
def maxnorm(num1, num2, num3):
    Z = np.max(num1 - num2)
    N = np.max(num2 - num3)
    Q_max = Z/N
    return Q_max

#convergence maxnorm
def maxnorm2(num1, num2, num3):
    Z = np.max(num1 - num2)
    N = np.max(num1 - num3)
    Q_max = Z/N
    return Q_max



def rungekutta(y, x, f):
    h = x[1] - x[0]
    y_out = np.zeros(((len(x),) + np.shape(y)))
    for i in range(len(x)):
        y_out[i] = y
        k1 = h*f(x[i], y)
        k2 = h*f(x[i] + 0.5*h, y + 0.5*k1)
        k3 = h*f(x[i] + 0.5*h, y + 0.5*k2)
        k4 = h*f(x[i] + h, y + k3)

        y = y + 1/6*(k1 + 2*k2 + 2*k3 + k4)

    return y_out, x


def error(numerical_values, analytical_values, x):
    error = abs(numerical_values - np.array(analytical_values))
    return error
