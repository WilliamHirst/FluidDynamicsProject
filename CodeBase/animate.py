import time, matplotlib.animation
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera

data = open("particleDensity.txt", "r")



values = []
for val in data.readline().split():
    values.append(float(val))
width = int(values[0])
height = int(values[1])
timeSteps = int(values[2])
rho = np.zeros((timeSteps, height, width))
xi = 0
yi = 0
t = 0
counter = 0
for line in data:
    t = counter // timeSteps
    xi = 0
    if t<timeSteps:
        for element in line.split():
            rho[t,xi,yi] = element
            xi += 1
            counter += 1
        yi += 1
        if yi == height:
            yi = 0

fig = plt.figure()
camera = Camera(fig)

for i in range(timeSteps):
    plt.imshow(rho[i])
    camera.snap()

animation = camera.animate()
animation.save('celluloid_subplots.gif', writer = 'imagemagick')
