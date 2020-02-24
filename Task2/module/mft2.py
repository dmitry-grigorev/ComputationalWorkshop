tune_b = 10 ** (-3)  # Настройка метода бисекции
tune_n = 10 ** (-3)  # Настройка метода Ньютона


class PolLagrange:
    def __init__(self, nodes, values):  # Инициализация класса, описывающего интерполяционный
        self.nodes = nodes  # многочлен Лагранжа
        val = values
        for k in range(len(nodes)):
            for l in range(len(nodes)):
                if k != l:
                    val[k] /= (nodes[k] - nodes[l])
        self.f_xk_div_omega_k = val

    def get_value(self, pivot):              # Вычисляем значение в конкретной точке
        s = 0
        val = self.f_xk_div_omega_k
        for k in range(len(self.nodes)):
            top = 1
            for l in range(len(self.nodes)):
                if k != l:
                    top *= (pivot - self.nodes[l])
            s += val[k] * top

        return s


def bisection(polynom, val, a, b, eps):      # Функция, выполняющая метод бисекции
    while b - a >= eps:
        c = (a + b) / 2
        if (polynom.get_value(c) - val) * (polynom.get_value(a) - val) < 0:
            a, b = (a, c)
        else:
            a, b = (c, b)

    return Newtons_method(polynom, a, b, eps * tune_b)

def Newtons_method(polynom, a, b, eps):      # Функция, выпоняющая метод Ньютона
    imax = 2000  # число итераций, максимальное
    i = 1
    x_n = b
    h = eps * tune_n
    f_xn = polynom.get_value(x_n)
    derf_xn = (polynom.get_value(x_n + h) - f_xn) / h
    x_n1 = x_n - f_xn / derf_xn

    while True:
        i += 1
        if x_n1 < a or x_n1 > b:
            while x_n1 < a or x_n1 > b:
                x_n1 = (x_n + x_n1) / 2
        if abs(x_n1 - x_n) < eps or i >= imax:
            return x_n1, i
        x_n = x_n1
        f_xn = polynom.get_value(x_n)
        derf_xn = (f_xn - polynom.get_value(x_n + h)) / h
        x_n1 = x_n - f_xn / derf_xn
