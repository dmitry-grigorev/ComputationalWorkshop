tune = 10 ** (-3)  # Настройка eps


class PolLagrange:
    def __init__(self, nodes, values):  # Инициализация класса, описывающего интерполяционный
        self.nodes = nodes  # многочлен Лагранжа
        val = values
        for k in range(len(nodes)):
            for l in range(len(nodes)):
                if k != l:
                    val[k] /= (nodes[k] - nodes[l])
        self.f_xk_div_omega_k = val

    def get_value(self, pivot):  # Вычисляем значение в конкретной точке
        s = 0
        val = self.f_xk_div_omega_k
        for k in range(len(self.nodes)):
            top = 1
            for l in range(len(self.nodes)):
                if k != l:
                    top *= (pivot - self.nodes[l])
            s += val[k] * top

        return s


def bisection(polynom, val, a, b, eps):  # Функция, выполняющая метод бисекции
    answer = []
    if abs(polynom.get_value(a) - val) < eps * tune:
        answer += [[a, 0]]
        a = a + eps
    if abs(polynom.get_value(b) - val) < eps * tune:
        answer += [[b, 0]]
        b = b - eps
    while b - a >= eps:
        c = (a + b) / 2
        if (polynom.get_value(c) - val) * (polynom.get_value(a) - val) < 0:
            if (polynom.get_value(c) - val) * (polynom.get_value(b) - val) < 0:
                answer += bisection(polynom, val, c, b, eps)
            else:
                if (polynom.get_value(c) - val) * (polynom.get_value((b + c) / 2) - val) < 0:
                    answer += bisection(polynom, val, c, (b + c) / 2, eps)
                if (polynom.get_value(b) - val) * (polynom.get_value((b + c) / 2) - val) < 0:
                    answer += bisection(polynom, val, (b + c) / 2, b, eps)
            a, b = (a, c)
        else:
            if (polynom.get_value(c) - val) * (polynom.get_value(a) - val) < 0:
                answer += bisection(polynom, val, a, c, eps)
            else:
                if (polynom.get_value(a) - val) * (polynom.get_value((a + c) / 2) - val) < 0:
                    answer += bisection(polynom, val, a, (a + c) / 2, eps)
                if (polynom.get_value(c) - val) * (polynom.get_value((a + c) / 2) - val) < 0:
                    answer += bisection(polynom, val, (a + c) / 2, c, eps)
            a, b = (c, b)
    answer += Newtons_method(polynom, a, b, val, eps * tune)
    return answer


def Newtons_method(polynom, a, b, val, eps):  # Функция, выпоняющая метод Ньютона
    imax = 2000  # число итераций, максимальное
    i = 1
    x_n = (a + b) / 2
    h = eps * tune
    while True:
        i += 1
        f_xn = polynom.get_value(x_n)
        derf_xn = (polynom.get_value(x_n + h) - f_xn) / h
        x_n1 = x_n - f_xn / derf_xn
        if x_n1 < a or x_n1 > b:
            while x_n1 < a or x_n1 > b:
                x_n1 = (x_n + x_n1) / 2
        if abs(x_n1 - x_n) < eps or i >= imax:
            return [[x_n1, i]]
        x_n = x_n1
