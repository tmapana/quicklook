#!/usr/bin/env python

"""
=========================================================
COPYRIGHT: UCT Radar Remote Sensing Group (1999)
FILE NAME: grnn.py
CODE CONTROLLER: Jasper Horrell (jasper@eng.uct.ac.za)
DESCRIPTION:
General regression neural network class. Expects initialization
with two lists of data of equal length (x and y values). The 
variance of the Gaussians used determines the fit to the data -
small variance implies fit approaches nearest neighbour, large 
variance implies more smoothing. Unit amplitude Gaussians are used.

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-08-09
COMMENTS:
Initial version. Updated to use stddev (1999-08-11).

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================================
"""

import sys, math, time, Gnuplot

class GenRegNN:

    def __init__(self,xvals=[],yvals=[]):
        if len(xvals) != len(yvals):
            print 'ERROR (in GenRegNN init) xy lists of different length!'
        self.datx = xvals
        self.daty = yvals
        self.stddev = 1.0     # std deviation for gaussians

    # Evaluate the regression curve at specified x value using
    # full data set (for faster operation use fastEval()). 
    def eval(self,x):
        numSum = 0.0
        denomSum = 0.0
        var = float(self.stddev) * float(self.stddev)
        for i in range(len(self.datx)):
            tmpG = x - self.datx[i]
            tmp = math.exp(-tmpG*tmpG/(2.0*var))
            numSum = numSum + self.daty[i]*tmp
            denomSum = denomSum + tmp   
        if denomSum == 0.0:
            return (0.0)
        else:
            return (numSum/denomSum)

    # Evaluate with speedup (potentially less accurate, as does not
    # use full data set for each evaluation). The 'ptsToUse' will depend
    # on the variance used (larger variance => more pts, I guess).
    # Note that this method automatically weights the end points
    # by using these values multiple times when approaching the
    # bounds of the array. This is not necessarily undesirable as it
    # reduces the GRNN assymptotic end effect. 

    def fastEval(self,x,ptsToUse):

        # find index close to the supplied x value
        n = 0
        maxindx = len(self.datx)-1
        while (self.datx[n] < x) and (n < maxindx):
            n = n+1
        indx = int(n - ptsToUse/2.0)  # start index    
        if indx < 0:   # ensure still in bounds
            indx = 0

        # Do actual calc
        numSum = 0.0
        denomSum = 0.0
        var = float(self.stddev) * float(self.stddev)
        for i in range(ptsToUse):
            tmpG = x - self.datx[indx]
            tmp = math.exp(-tmpG*tmpG/(2.0*var))
            numSum = numSum + self.daty[indx]*tmp
            denomSum = denomSum + tmp   
            indx = indx+1
            if indx > maxindx:  # ensure still in bounds
                indx = maxindx
        if denomSum == 0.0:
            return (0.0)
        else:
            return (numSum/denomSum)
    

    # Plot the orig data and the regression curve with specified no. steps
    # flag - a string which is either sleep, file, wait
    def display(self,numsteps,flag):
        stepsize = (max(self.datx)-min(self.datx))/(numsteps-1)
        indx = min(self.datx)
        x,y = [],[]
        for i in range(numsteps):
            x.append(indx)
            y.append(self.eval(indx))
            indx = indx + stepsize
        pl = Gnuplot.Gnuplot()
        dat0 = Gnuplot.Data(self.datx,self.daty,with='lines',\
                            title='Orig data')
        grnnTitle='GRNN (std dev '+str(self.stddev)+')'                     
        dat1 = Gnuplot.Data(x,y,with='linespoints',title=grnnTitle)
        pl.plot(dat0,dat1)
        if flag=='sleep':
            time.sleep(5.0)
        elif flag=='file':
            time.sleep(5.0)
            pl.hardcopy(filename='grnn_display.ps',color=1,fontsize=20)
        elif flag=='wait':
            print 'Press return to continue...'
            sys.stdin.readline()
        return 0      
            
    def displayFast(self,numsteps,ptsToUse,flag):  # to plot the speed up approx
        stepsize = (max(self.datx)-min(self.datx))/(numsteps-1)
        indx = min(self.datx)
        x,y,fasty = [],[],[]
        for i in range(numsteps):
            x.append(indx)
            y.append(self.eval(indx))
            fasty.append(self.fastEval(indx,ptsToUse))
            indx = indx + stepsize
        pl = Gnuplot.Gnuplot()
        dat0 = Gnuplot.Data(self.datx,self.daty,with='linespoints',\
                            title='Orig data')
        grnnTitle='GRNN (std dev '+str(self.stddev)+')'                     
        dat1 = Gnuplot.Data(x,y,with='linespoints',title=grnnTitle)
        dat2 = Gnuplot.Data(x,fasty,with='linespoints',title='GRNN with speedup')
        pl.plot(dat0,dat1,dat2)
        time.sleep(5.0)
        pl.hardcopy(filename='grnn_display.ps',color=1,fontsize=20)
        if flag=='sleep':
            time.sleep(5.0)
        elif flag=='file':
            time.sleep(5.0)
            pl.hardcopy(filename='grnn_display.ps',color=1)
        elif flag=='wait':
            print 'Press return to continue...'
            sys.stdin.readline()
        return 0
          
  
# demo code
if __name__ == '__main__':          
    # create some data  
    x = [1.0,2.0,3.0,4.2,6.0,7.0]
    y = [2.0,4.0,4.6,6.1,4.0,4.7]

    # set up three GRNN curves with different variances
    nn1,nn2,nn3 = GenRegNN(x,y), GenRegNN(x,y), GenRegNN(x,y)  
    nn1.stddev,nn2.stddev,nn3.stddev = 0.1,0.8,2.0
    nn1x,nn1y,nn2y,nn3y, = [],[],[],[]    # lists used for plots
    steps = 40
    stepsize = (max(x)+1.0)/(steps-1)
    indx = min(x) - 1.0
    for i in range(steps):
        nn1x.append(indx)
        nn1y.append(nn1.eval(indx))
        nn2y.append(nn2.eval(indx))
        nn3y.append(nn3.eval(indx))
        indx = indx + stepsize

    # plot the orig data plus the three GRNN's using Gnuplot
    gp = Gnuplot.Gnuplot()
    d0 = Gnuplot.Data(x,y,with='linespoints',title='Orig data')        
    d1 = Gnuplot.Data(nn1x,nn1y,with='linespoints',title='GRNN (std dev 0.1)')
    d2 = Gnuplot.Data(nn1x,nn2y,with='linespoints',title='GRNN (std dev 0.8)')
    d3 = Gnuplot.Data(nn1x,nn3y,with='linespoints',title='GRNN (std dev 2.0)')
    gp.plot(d0,d1,d2,d3)
    print 'Press return to continue...'
    sys.stdin.readline()
