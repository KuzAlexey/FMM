import math
import random

def generate_points_in_circle(n, radius):
    points = []
    for i in range(n):
        angle = random.uniform(0, 2*math.pi)
        bangle = random.uniform(0, 2*math.pi)
        d = random.uniform(0, radius)
        x = (radius) * math.cos(angle)
        y = (radius) * math.sin(angle)
        z = (radius) * math.sin(angle)
        m = random.uniform(1, 1)
        v_x = 0
        v_y = 0
        v_z = 0
        points.append((x, y, m))
    return points

def gggenerate_points_in_square(n):
    points = []
    used = []
    for i in range(n):
        x = random.uniform(-100, 100)
        y = random.uniform(-100, 100)
        m = random.uniform(1, 100)  # generating random value for M
        x = int(x) + 0.5
        y = int(y) + 0.5
        if (x, y) in used:
            continue
        used.append((x, y))
        points.append((x, y, m))
    return points

def generate_points_in_square(n):
    points = []
    for i in range(n):
        x = random.uniform(-10000, 10000)
        y = random.uniform(-10000, 10000)
        m = random.uniform(1, 100)  # generating random value for M
        points.append((x, y, m))
    return points

n = 5000
def write_points_to_file(points, filename):
    global n
    with open(filename, 'w') as file:
        file.write(str(len(points))+"\n")
        for point in points:
            file.write(' '.join(str(value) for value in point) + '\n')

def my_write_points_to_file(points, filename):
    global n
    points = [(30, 50, 0, 1, 0, 0, 0),
              (80, 10, 0, 1, 0, 0, 0),
              (90, 10, 0, 1, 0, 0, 0),
              (90, 20, 0, 1, 0, 0, 0),
              (90, 50, 0, 1, 0, 0, 0),
              (60, 90, 0, 1, 0, 0, 0),]
    with open(filename, 'w') as file:
        file.write(str(len(points))+"\n")
        for point in points:
            file.write(' '.join(str(value) for value in point) + '\n')

# Generate 10 points within a circle of radius 5
radius = 500
points = generate_points_in_square(n)
# points = generate_points_in_square(n)

# Write the points to a file
filename = 'input.txt'
write_points_to_file(points, filename)
print(f"Points data written to {filename} successfully.")