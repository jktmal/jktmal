import matplotlib.pyplot as plt
import numpy as np

def update_line(pct, old, new_data):
    #print([ el[0] for el in new_data ])
    old.set_xdata(new_data[0])
    old.set_ydata(new_data[1])
    #input()
    #hl.set_xdata(numpy.append(hl.get_xdata(), new_data))
    #hl.set_ydata(numpy.append(hl.get_ydata(), new_data))
    pct.canvas.draw()
    pct.canvas.flush_events()
    

'''

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'r-') # Returns a tuple of line objects, thus the comma

for phase in np.linspace(0, 10*np.pi, 500):
    line1.set_ydata(np.sin(x + phase))
    fig.canvas.draw()
    fig.canvas.flush_events()'''
