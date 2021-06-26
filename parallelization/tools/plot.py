import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

def plot_Q6_standard():
    data = np.genfromtxt("out/HO_2dsolver.csv", delimiter=",", skip_header=1, comments="#").transpose()

    start = data[0][0]
    stop = data[0][-1]
    x = np.arange(start, stop, 0.1)
    y_a = np.cos(x)
    plt.plot(data[0], data[1], label="num")
    plt.plot(x, y_a, label="ana", linestyle="dashed")

    plt.xlabel("t")
    plt.ylabel("q(t)")

    plt.ylim((-1.5, 1.5))

    plt.title("dt=0.01, m=1")

    plt.legend(loc="best")
    plt.show()

def plot_1d_wave():
    file_names = ["0.5", "1.0", "1.5", "2.0", "2.5", "3.0"]
    file_names = ["out/t=" + t + ".dat" for t in file_names]

    x = np.linspace(-5, 5, 1000)
    opacity = np.linspace(0, 1, len(file_names) + 1)
    for i, file in enumerate(file_names):
        y = np.genfromtxt(file, delimiter=",")
        plt.plot(x, y, color=(0, 0, 0, opacity[i + 1]))
    
    plt.show()

def plot_3d(file):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    X = np.linspace(-5, 5, 100)
    Y = np.linspace(-5 ,5, 100)

    Z = 0 # read in data

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

    # # Customize the z axis.
    # ax.set_zlim(-1.01, 1.01)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # # A StrMethodFormatter is used automatically
    # ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

if __name__ == "__main__": 
    plot_3d("test")