import numpy as np
from matplotlib import pyplot as plt


# моделирование равномерного распределения в круге радиуса 10
def get_vertices_in_circle(amount):
    radii = np.sqrt(np.random.uniform(0, 100, amount))
    phi = np.random.uniform(0, 2 * np.pi, amount)

    xvert = radii * np.cos(phi)
    yvert = radii * np.sin(phi)

    dots = np.concatenate((xvert, yvert), axis=0)
    dots = np.reshape(dots, newshape=(2, amount))

    return dots


# вычисляет синус половинного угла между векторами v1 и v2
def sine(v1, v2):
    return 1 - np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))


# строит выпуклую оболочку для точек dots
def build_coverance(dots):
    dots_sorted = dots[:, np.argsort(dots[1, :])]
    coverance = [dots_sorted[:, 0]]
    pivot = np.array([-1, 0])
    i = 0
    inext = 0
    while inext != 0 or i == 0:
        i += 1
        sindistances = np.apply_along_axis(lambda x: sine(x - coverance[i - 1], pivot), axis=0, arr=dots_sorted)
        inext = np.argmax(np.nan_to_num(sindistances))
        coverance.append(dots_sorted[:, inext])
        pivot = coverance[i - 1] - coverance[i]
    return np.array(coverance)


# триангулирует выпуклый многоугольник
def triangularize(coverance):
    return np.array([[coverance[0], coverance[i], coverance[i + 1]] for i in range(1, coverance.shape[0] - 2)])


# вычисляет площадь треугольника triangle, заданного своими точками
def get_area(triangle):
    v1 = triangle[1] - triangle[0]
    v2 = triangle[2] - triangle[0]
    return np.abs(v1[0] * v2[1] - v1[1] * v2[0]) / 2


# моделирование равномерного распределения в треугольнике triangle, заданного своими точками
def triangle_uniform(triangle):
    v1 = triangle[1] - triangle[0]
    v2 = triangle[2] - triangle[0]
    alpha1, alpha2 = np.random.uniform(), np.random.uniform()
    alpha1, alpha2 = (alpha1, alpha2) if alpha1 + alpha2 < 1 else (1 - alpha1, 1 - alpha2)
    beta1 = v1[0] * alpha1 + v2[0] * alpha2 + triangle[0, 0]
    beta2 = v1[1] * alpha1 + v2[1] * alpha2 + triangle[0, 1]
    return [beta1, beta2]


# моделирование равномерного распределения в многоугольнике
def polygon_uniform(triangles, amount):
    areas = np.array([get_area(triangle) for triangle in triangles])
    triangles_probs = areas / np.sum(areas)

    itriangles = np.random.choice(triangles.shape[0], size=amount, p=triangles_probs)

    generated = np.array([triangle_uniform(triangles[i]) for i in itriangles])
    return generated


def plot_dots(dots):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.scatter(dots[0, :], dots[1, :], color='g')
    ax.set_aspect(1)
    ax.grid()


def plot_dots_and_triangles(dots, triangles, coverance):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(coverance[:, 0], coverance[:, 1])
    for triangle in triangles:
        plt.plot(triangle[:, 0], triangle[:, 1])
    plt.scatter(dots[0, :], dots[1, :], color='g')
    ax.set_aspect(1)
    ax.grid()


def plot_coverance(dots, coverance):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.scatter(dots[0, :], dots[1, :], color='g')
    plt.plot(coverance[:, 0], coverance[:, 1])
    ax.set_aspect(1)
    ax.grid()


def plot_triangles(coverance, triangles):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(coverance[:, 0], coverance[:, 1])
    for triangle in triangles:
        plt.plot(triangle[:, 0], triangle[:, 1])
    ax.set_aspect(1)
    ax.grid()
