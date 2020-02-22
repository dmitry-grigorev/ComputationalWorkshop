class PolLagrange:
    def __init__(self, nodes, values):
        self.nodes = nodes
        self.values = values
        val = values
        for k in range(len(nodes)):
            for l in range(len(nodes)):
                if k != l:
                    val[k] /= (nodes[k] - nodes[l])
        self.f_xk_div_omega_k = val

    def polynom_Lagrange(self, nodes, values,pivot):
        s = 0
        val = self.f_xk_div_omega_k
        for k in range(len(nodes)):
            top = 1
            for l in range(len(nodes)):
                if k != l:
                    top *= (pivot - nodes[l])
            s += val[k] * top

        return s  # Собственно,вычисляем значением