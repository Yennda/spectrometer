import timeit
import numpy as np
from algebra import la
from tools import Tools

def test_add():
    a = [1.2, 2.04, 3.7]
    b = [4.5, 3.2, 7.8]
    if len(a) != len(b):
        raise ValueError
    return [a[i] + b[i] for i in range(len(a))]


def test_multiply():
    a = [1.2, 2.04, 3.7]
    return [a[i] * 3.14 for i in range(len(a))]


def test_sp():
    a = [1.2, 2.04, 3.7]
    b = [4.5, 3.2, 7.8]
    if len(a) != len(b):
        raise ValueError
    s = 0
    for i in range(len(a)):
        s += a[i] * b[i]
    return s

def rock_curve(x):
    return tl.gauss(x, mi=self.bragg, s=0.0014544410433286077 / 3)

print('add')
print('array')
print(timeit.timeit('np.array([1.2,2.04,3.7])+np.array([4.5,3.2,7.8])', setup='import numpy as np', number=10000))
print('matrix')
print(timeit.timeit('np.matrix([1.2,2.04,3.7])+np.matrix([4.5,3.2,7.8])', setup='import numpy as np', number=10000))
print('list')
print(timeit.timeit('test_add()', setup='from __main__ import test_add', number=10000))
print('LA')
print(timeit.timeit('la.plus([1.2,2.04,3.7],[4.5,3.2,7.8])', setup='from algebra import la', number=10000))
print('------------------------')
print('multiply')
print('array')
print(timeit.timeit('np.array([1.2,2.04,3.7])*3.14', setup='import numpy as np', number=10000))
print('matrix')
print(timeit.timeit('np.matrix([1.2,2.04,3.7])*3.14', setup='import numpy as np', number=10000))
print('list')
print(timeit.timeit('test_multiply()', setup='from __main__ import test_multiply', number=10000))
print('LA')
print(timeit.timeit('la.x([1.2,2.04,3.7],3.14)', setup='from algebra import la', number=10000))
print('------------------------')
print('scalar product')
print('array')
print(
    timeit.timeit('np.array([1.2,2.04,3.7]).dot(np.array([4.5,3.2,7.8]))', setup='import numpy as np', number=10000))
print('matrix')
print(timeit.timeit('np.matrix([1.2,2.04,3.7])*np.matrix([4.5,3.2,7.8]).T', setup='import numpy as np', number=10000))
print('list')
print(timeit.timeit('test_sp()', setup='from __main__ import test_sp', number=10000))
print('LA')
print(timeit.timeit('la.dot([1.2,2.04,3.7],[4.5,3.2,7.8])', setup='from algebra import la', number=10000))
print('------------------------')
print('normalize')
print('Tools')
print(timeit.timeit('Tools.normalize(np.array([1.2,2.04,3.7]))', setup='from tools import Tools;import numpy as np',
                    number=10000))
print('LA')
print(timeit.timeit('la.normalize([1.2,2.04,3.7])', setup='from algebra import la', number=10000))
print('------------------------')
print('matrix x vector')
print('Matrix')
print(timeit.timeit('np.matrix([[1.4,2.7,3],[4.3,5.8,6.4],[7.3,8.3,9.7]])*np.matrix([2.6,3.2,6.5]).T',
                    setup='import numpy as np', number=10000))
print('LA')
print(timeit.timeit('la.dotmv([[1.4,2.7,3],[4.3,5.8,6.4],[7.3,8.3,9.7]],[2.6,3.2,6.5])', setup='from algebra import la',
                    number=10000))

print('------------------------')
print('object x dict')
print('object')
print(timeit.timeit('dp.loc',
                    setup='from tracing import DetectorPoint;dp=DetectorPoint([2.2,6.7,7.8])', number=10000))
print('dict')
print(timeit.timeit('dp["loc"]', setup='dp={"loc":[2.2,6.7,7.8]}',
                    number=10000))

print('rocking_curves')
print('gauss')
print(timeit.timeit('tl.gauss(1.383, mi=1.382543, s=0.0014544410433286077 / 3)',
                    setup='from tools import Tools as tl', number=10000))
print('FILE')
print(timeit.timeit('curveSi.curve(1.383)', setup='from curves import Curves;curveSi = Curves("reflectivity.csv")',
                    number=10000))

print('------------------------')
print('solve x matrix')
print('matrix')
print(timeit.timeit('np.matrix([[2,4,1],[2,5,1],[0,4,5]]).I*np.matrix([1,2,3]).T',
                    setup='import numpy as np; ', number=10000))
print('solve np.array')
print(timeit.timeit('np.linalg.solve(np.array([[2,4,1],[2,5,1],[0,4,5]]),np.array([1,2,3]))', setup='import numpy as np',
                    number=10000))
print('solve list')
print(timeit.timeit('np.linalg.solve([[2,4,1],[2,5,1],[0,4,5]],[1,2,3])', setup='import numpy as np',
                    number=10000))
print('########################')
