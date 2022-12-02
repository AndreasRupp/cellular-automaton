import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler

 
def plot_update(axes, data, ax):
  dim = np.size(axes)
  if dim == 1:
    axes[1] = 1
    dim = 2

  if dim == 2:
    data = (data != 0)
    cmap = colors.ListedColormap(['blue', 'red'])
    # plt.figure(figsize=(6,6))
    ax.pcolor(data[::-1],cmap=cmap,edgecolors='k', linewidths=3)
    # ax.show()
  elif dim == 3:
    data = (data != 0)
    Colors = np.empty(axes + [4], dtype=np.float32)
    # Control Transparency
    alpha = .9
    Colors[:] = [0, 0, 1, alpha]
    # fig = plt.figure(figsize=(6,6))
    # ax = fig.add_subplot(111, projection='3d')
    ax.voxels(data, facecolors=Colors, edgecolors='black')
    # plt.show()
  return ax



def plot(axes, save_data, time_step):
  if time_step < 0 or time_step >= len(save_data):
    print("Time step for plotting does not exist!")
    return

  root = Tk()
  animated_cam(root, axes, save_data, time_step)
  root.mainloop()


def config_plot():
    fig, ax = plt.subplots()
    ax.set(xlabel='time (s)', ylabel='voltage (mV)',
           title='Graph One')
    return (fig, ax)

class animated_cam:
  def __init__(self, master, axes, save_data, time_step):
    self.master = master
    self.frame = Frame(self.master)
    self.fig, self.ax = config_plot()
    self.canvas = FigureCanvasTkAgg(self.fig, self.master)  
    self.config_window()
    self.axes      = axes
    self.save_data = save_data
    self.time_step = time_step
    self.n_steps   = len(save_data)
    self.draw_cam()
    self.frame.pack(expand=YES, fill=BOTH)

  def config_window(self):
    self.canvas.mpl_connect("key_press_event", self.on_key_press)
    toolbar = NavigationToolbar2Tk(self.canvas, self.master)
    toolbar.update()
    self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    self.button = Button(self.master, text="Quit", command=self._quit)
    self.button.pack(side=BOTTOM)
    self.button_next = Button(self.master, text="Next", command=self.plot_next)
    self.button_next.pack(side=BOTTOM)
    self.button_prev = Button(self.master, text="Previous", command=self.plot_previous)
    self.button_prev.pack(side=BOTTOM)

  def draw_cam(self):
    self.ax.clear() # clear current axes
    data = self.save_data[self.time_step]
    self.ax = plot_update(self.axes, np.reshape(data, self.axes), self.ax)
    self.ax.set(title='Graph One')
    self.canvas.draw()

  def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, self.canvas, toolbar)

  def _quit(self):
    self.master.quit()  # stops mainloop

  def plot_next(self):
    self.time_step = (self.time_step + 1 ) % self.n_steps
    self.draw_cam()

  def plot_previous(self):
    self.time_step = (self.time_step - 1 + self.n_steps ) % self.n_steps
    self.draw_cam()
