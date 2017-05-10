class la:
    @staticmethod
    def plus(a, b):
        return [a[i] + b[i] for i in range(3)]

    @staticmethod
    def minus(a, b):
        return [a[i] - b[i] for i in range(3)]

    @staticmethod
    def minustwo(a, b):
        return [a[i] - b[i] for i in range(2)]

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
    def vec(a, b):
        return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

    @staticmethod
    def area(a, b):
        return la.norm(la.vec(a, b))

    @staticmethod
    def cos(a, b):
        return la.dot(la.normalize(a), la.normalize(b))

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
