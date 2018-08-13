"""
FFT using Cooley-Tukey algorithm where N = 2^n.  The choice of
e^{-j2\pi/N} or e^{j2\pi/N} is made by 'sign=-1' or 'sign=1'
respectively.  Since I prefer Engineering convention, I chose
'sign=-1' as the default.

FFT is performed as follows:
1. bit-reverse the array.
2. partition the data into group of m = 2, 4, 8, ..., N data points.
3. for each group with m data points,
    1. divide into upper half (section A) and lower half (section B),
	each containing m/2 data points.
    2. divide unit circle by m.
    3. apply "butterfly" operation 
	    |a| = |1  w||a|	or	a, b = a+w*b, a-w*b
	    |b|   |1 -w||b|
	where a and b are data points of section A and B starting from
	the top of each section, and w is data points along the unit
	circle starting from z = 1+0j.
FFT ends after applying "butterfly" operation on the entire data array
as whole, when m = N.
"""
def fft(x, sign=-1):
    from cmath import pi, exp
    N, W = len(x), []
    for i in range(N):		# exp(-j...) is default
	W.append(exp(sign * 2j * pi * i / N))
    x = bitrev(x)
    m = 2
    while m <= N:
	for s in range(0, N, m):
	    for i in range(m/2):
		n = i * N / m
	    	a, b = s + i, s + i + m/2
	        x[a], x[b] = x[a] + W[n % N] * x[b], x[a] - W[n % N] * x[b]
        m = m * 2
    return x
    
"""
Inverse FFT with normalization by N, so that x == ifft(fft(x)) within
round-off errors.
"""
def ifft(X):
    N, x = len(X), fft(X, sign=1)	# e^{j2\pi/N}
    for i in range(N):
	x[i] = x[i] / float(N)
    return x