import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import numpy.matlib
import keyboard
import math

def main():
    axes = [5, 5]
    data = np.ones(axes)
    data[0][0] = 0
    rate = 1
    steps = 1

    plot(axes, data, rate, steps)

def plot(axes, data, rate, steps):
    if rate == 0:
        return

    i = 0
    while True:
        # plt.title("Time step: " + str(i))
        plot_update(axes, data)
        event = keyboard.read_event()
        if event.event_type == "down":
            match event.name:
                case "1":
                    i = (i + 1) % steps
                case "2":
                    i = (i - 1) % steps
                case other:
                    print("Kiitos ohjelman käytöstä.")
                    return
 
def plot_update(axes, data):

    dim = np.size(axes)
    if dim == 1:
        axes[1] = 1
        dim = 2

    x = range(axes[0])
    y = range(axes[1])
    if dim == 3:
        z = range(axes[2])

    # ax.clear()

    if dim == 2:
        # x = np.matlib.repmat(range(0,np.size(data)-1) % axes[0], 5)
        # x = x[(data != 0)][:]
        # y = np.matlib.repmat(math.floor(range(0,np.size(data)-1)/axes[0]) % axes[1], 5)
        # y = y[(data != 0)][:]

        cmap = colors.ListedColormap(['blue', 'red'])
        plt.figure(figsize=(6,6))
        plt.pcolor(data[::-1],cmap=cmap,edgecolors='k', linewidths=3)
        plt.show()

    if dim == 3:
        # x = np.matlib.repmat(range(0,np.size(data)-1) % axes[0], 5)
        # x = x[(data != 0)][:]
        # y = np.matlib.repmat(math.floor(range(0,np.size(data)-1)/axes[0]) % axes[1], 5)
        # y = y[(data != 0)][:]
        # z = np.matlib.repmat(math.floor(range(0,np.size(data)-1)/(axes[0]*axes[1])) % axes[2], 5)
        # z = z[(data != 0)][:]

        Colors = np.empty(axes + [4], dtype=np.float32)
        # Control Transparency
        alpha = .9
        Colors[:] = [0, 0, 1, alpha]
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.voxels(data, facecolors=Colors, edgecolors='black')
        plt.show()

main()