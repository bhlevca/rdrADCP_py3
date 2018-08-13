
from matplotlib.patches import Patch
from pylab import *

def func3(x, y):
    return (1 - x / 2 + x ** 5 + y ** 3) * exp(-x ** 2 - y ** 2)


# make these smaller to increase the resolution
dx, dy = 0.05, 0.05

x = arange(-3.0, 3.0001, dx)
y = arange(-3.0, 3.0001, dy)
X, Y = meshgrid(x, y)

Z = func3(X, Y)
pcolor(X, Y, Z)
colorbar()
axis([-3, 3, -3, 3])

show()

Z = rand(6, 10)

subplot(2, 1, 1)
c = pcolor(Z)
title('default: no edges')

subplot(2, 1, 2)
c = pcolor(Z, edgecolors = 'k', linewidths = 4)
title('thick edges')

show()


import matplotlib.pyplot as plt
t = np.arange(256)
sp = np.fft.fft(np.sin(t))
tt = t.shape[-1]
freq = np.fft.fftfreq(tt)
plt.plot(freq, sp.real, freq, sp.imag)

plt.show()
