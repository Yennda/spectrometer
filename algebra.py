class la:
    @staticmethod
    def plus(a, b):
        return [a[i] + b[i] for i in range(3)]

    @staticmethod
    def minus(a, b):
        return [a[i] - b[i] for i in range(3)]

    @staticmethod
    def x(a, x):
        return [a[i] * x for i in range(len(a))]

    @staticmethod
    def i(a):
        return [-a[i] for i in range(len(a))]

    @staticmethod
    def dot(a, b):
        s = 0
        for i in range(3):
            s += a[i] * b[i]
        return s

    @staticmethod
    def dotmv(A, b):
        out = list()
        for r in A:
            out.append(la.dot(r, b))
        return out

    @staticmethod
    def norm(a):
        s = 0
        for i in a:
            s += i ** 2
        return s ** 0.5

    @staticmethod
    def normalize(a):
        return la.x(a, 1 / la.norm(a))
