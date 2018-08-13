import numpy as np
import matplotlib.pyplot as plt


def autolabel(rects, values):
    # attach some text labels
    idx = 0
    for rect in rects:
        #height = rect.get_height()
        height = values[idx]
        if height > 0 :
            placement = height + 25
        else :
            placement = height - 25

        plt.text(rect.get_x() + rect.get_width() / 2., placement, '%d' % int(height),
                ha = 'center', va = 'bottom')
        idx += 1

TrenberthSurface = (333, 161, -17, -80, -396)
TrenberthSurfaceStd = (0, 0, 0, 0, 0)
StephensSurface = (347, 168, -24, -88, -398)
StephensSurfaceStd = (15, 10, 10, 10, 5)

TrenberthTOA = (0, 341, 0, 0, -239, -102)
TrenberthTOAStd = (0, 0, 0, 0, 0, 0)
StephensTOA = (0, 340.2, 0, 0, -238.9, -97.7)
StephensTOAStd = (0, 0.1, 0, 0, 3.3, 2)


N = 5

ind = np.arange(N)  # the x locations for the groups
width = 0.45       # the width of the bars


plt.subplot(121)
rects1 = plt.bar(ind, TrenberthSurface, width,
                    color = 'r',
                    yerr = TrenberthSurfaceStd,
                    error_kw = dict(elinewidth = 6, ecolor = 'pink'))

rects2 = plt.bar(ind + width, StephensSurface, width,
                    color = 'y',
                    yerr = StephensSurfaceStd,
                    error_kw = dict(elinewidth = 6, ecolor = 'yellow'))

# add some
plt.ylabel('heat fluxes W/m^2')
plt.title('Heat Fluxes for Surface')
plt.xticks(ind + width, ('FLW', 'FSW', 'FSH', 'FLH', 'OLR', 'OSW'))

plt.grid()
plt.legend((rects1[0], rects2[0]), ('Trenberth', 'Stephens'))
autolabel(rects1, TrenberthSurface)
autolabel(rects2, StephensSurface)


N = 6

ind = np.arange(N)  # the x locations for the groups
width = 0.45       # the width of the bars

plt.subplot(122)
rects3 = plt.bar(ind, TrenberthTOA, width,
                   color = 'r',
                   yerr = TrenberthTOAStd,
                   error_kw = dict(elinewidth = 6, ecolor = 'pink'))
rects4 = plt.bar(ind + width, StephensTOA, width,
                  color = 'y',
                  yerr = StephensTOAStd,
                  error_kw = dict(elinewidth = 6, ecolor = 'yellow'))

# add some
plt.ylabel('heat fluxes W/m^2')
plt.title('Heat Fluxes for TOA')

plt.xticks(ind + width, ('FLW', 'FSW', 'FSH', 'FLH', 'OLR', 'OSW'))
plt.grid()
plt.legend((rects3[0], rects4[0]), ('Trenberth', 'Stephens'))
autolabel(rects3, TrenberthTOA)
autolabel(rects4, StephensTOA)





plt.show()
