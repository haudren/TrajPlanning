#! /usr/bin/env python
#

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

if __name__ == '__main__':
  file = sys.argv[1]
  obsCost = []
  smCost = []
  speedCost = []
  bodyPos = []
  execfile(file)

  fig = plt.figure()
  ax = fig.add_subplot(111)

  line, = ax.plot([], [], 'o-', lw=2)
  time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

  def animate(i):
    X = np.array(bodyPos[i])[:,0]
    Y = np.array(bodyPos[i])[:,1]
    line.set_data(X, Y)
    return line, time_text

  ani = animation.FuncAnimation(fig, animate, np.arange(1, len(bodyPos)), interval=200)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(obsCost)
  ax.plot(smCost)
  ax.plot(speedCost)
  plt.show()
