import numpy as np
import matplotlib.pyplot as plt
import corner

data = np.genfromtxt('iteraciones.csv', delimiter=',')

figure = corner.corner(data)
plt.savefig('epicentro.pdf')
