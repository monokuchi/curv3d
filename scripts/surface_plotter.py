
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import re



fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection="3d")


points = []
x_coords = []
y_coords = []
z_coords = []


file = open("points.txt", "r")

for line in file:
    point = [item.strip() for item in re.findall("(?:\d{1,10}(?:,\d{10})*|\d+)(?:\.\d{1,10})?", line.strip())]
    # print(point)
    points.append(point)
    
for point in points:
    x_coords.append(float(point[0]))
    y_coords.append(float(point[1]))
    z_coords.append(float(point[2]))

file.close()


ax.scatter3D(x_coords, y_coords, z_coords, color="green", marker='o')
plt.title("Bezier Surface")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.savefig("surface_plot.jpg")
