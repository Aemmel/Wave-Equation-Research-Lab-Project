import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    data = pd.read_csv("out/test.csv")
    data = np.array(data)

    plt.plot(data[0], data[1])
    plt.show()