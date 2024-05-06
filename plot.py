import matplotlib.pyplot as plt


file_path ="Computational Methods/iris.txt"
data = []

with open(file_path, 'r') as file:
    for line in file:
        values = list(map(float, line.split()))
        data.append(values)

x = [row[0] for row in data]
y = [row[1] for row in data]
colors = [row[2] for row in data]

color_map = {
    0: 'red',
    1: 'green',
    2: 'blue',
}

mapped_colors = [color_map[color] for color in colors]

plt.figure(figsize=(10, 6))
plt.scatter(x, y, c=mapped_colors, label='Data Points', alpha=0.7)
plt.title('Scatter Plot with Color-Coded Points')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.grid(True)
plt.show()
