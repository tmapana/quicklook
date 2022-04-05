#!/usr/bin/env python

"""
Prog to interpolate DGPS data
code: J.M. Horrell 1999-07-28
"""

import sys, Polynomial, Gnuplot

# Misc
order = 3
pts = [1.,2.5,3,4,6,8]
vals = [2,2.4,5.,5,2,6]

# fit polynomial
fit = Polynomial.fitPolynomial(order,pts,vals)
print fit.coeff

fitX, fitY = [], [];
steps = 20
incr = (max(pts)-min(pts))/(steps-1)
x = min(pts)
for b in range(steps):
  fitX.append(x)
  fitY.append(fit(x))
  x = x + incr

# plot using gnuplot
g1 = Gnuplot.Gnuplot(debug=1)
d1a = Gnuplot.Data(pts,vals,title='Orig Data',with='linespoints')
d1b = Gnuplot.Data(fitX,fitY,title='Polynomial Fit',with='linespoints')
g1.plot(d1a,d1b)
print 'Press return to continue...'
sys.stdin.readline()
del g1,d1a,d1b
