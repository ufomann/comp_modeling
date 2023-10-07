import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=4, ncols=4)
fig.set_size_inches(10, 10)
ph1 = 0
ph2 = 0
for i in range(1, 5):
    for j in range(1, 5):
        t = np.linspace(0, 2 * np.pi / np.gcd(i, j), 100)

        x = np.cos(i * t + ph1)
        y = np.sin(j * t + ph2)
        ax[i - 1][j - 1].plot(x, y, 'o-')
        ax[i - 1][j - 1].set_xlabel('x')
        ax[i - 1][j - 1].set_ylabel('y')
        ax[i - 1][j - 1].set_title(f'{i} and {j}')
fig.tight_layout()
fig.savefig('task2/fig.png')
fig.show()
