import timeit


def test():
    a = [1.2, 2.04, 3.7]
    b = [4.5, 3.2, 7.8]
    return [a[i] + b[i] for i in range(len(a))]


print(timeit.timeit('np.array([1.2,2.04,3.7])+np.array([4.5,3.2,7.8])', setup='import numpy as np', number=10000))

print(timeit.timeit('np.matrix([1.2,2.04,3.7])+np.matrix([4.5,3.2,7.8])', setup='import numpy as np', number=10000))

print(timeit.timeit('test()', setup='from __main__ import test', number=10000))
