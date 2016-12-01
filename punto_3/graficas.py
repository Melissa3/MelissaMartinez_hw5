import numpy as np
import matplotlib.pyplot as plt
import corner

data = np.genfromtxt('parametros.csv', delimiter=',', usecols=(1,2,3))

figure = corner.corner(data)
plt.savefig('parametros.pdf')
