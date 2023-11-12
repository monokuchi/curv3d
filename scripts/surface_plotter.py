
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection="3d")


points = []
x_coords = []
y_coords = []
z_coords = []

control_points = []
control_x_coords = []
control_y_coords = []
control_z_coords = []


with open("points.txt", "r") as file:
    control_points_flag = False

    for line in file:
        line = line.strip()
        if line != "ControlPoints" and not control_points_flag:
            point = [item.strip() for item in line.split(",")]
            points.append(point)
        elif control_points_flag:
            control_point = [item.strip() for item in line.split(",")]
            control_points.append(control_point)
        elif not control_points_flag:
            control_points_flag = True
            continue
            
    
for point in points:
    x_coords.append(float(point[0]))
    y_coords.append(float(point[1]))
    z_coords.append(float(point[2]))

for control_point in control_points:
    control_x_coords.append(float(control_point[0]))
    control_y_coords.append(float(control_point[1]))
    control_z_coords.append(float(control_point[2]))


# Plot Bezier Surface
ax.scatter(x_coords, y_coords, z_coords, c=z_coords, cmap="jet", alpha=.3)

# Plot Control Points
ax.scatter(control_x_coords, control_y_coords, control_z_coords, c="black", marker="*", s=200)
# ax.plot(control_x_coords, control_y_coords, control_z_coords, c="red")

plt.title("Bezier Surface with Control Points")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# plt.show()
plt.savefig("surface_plot.jpg")
