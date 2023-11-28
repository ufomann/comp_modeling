import numpy as np
import matplotlib.pyplot as plt
import json

class Segment:
    def __init__(self, begin, end):
        self.begin = begin 
        self.end = end

fig, ax = plt.subplots(nrows=2)
fig.set_size_inches(10, 10)

#configuration
config_data = dict()
with open ('week2/config.json', 'r') as outfile:
    config_data = json.load(outfile)
ph1 = config_data['phase1']
ph2 = config_data['phase2']
freq1 = config_data['freq1']
freq2 = config_data['freq2']
count = config_data['count']

t = np.linspace(0, 2 * np.pi / np.gcd(freq1, freq2), count)
x = np.cos(freq1 * t + ph1)
y = np.sin(freq2 * t + ph2)
segs = [Segment(np.array([x[i], y[i]], dtype=float), 
                np.array([x[(i + 1) % count], y[(i + 1) % count]], dtype=float)) for i in range(count)]
ax[0].plot(x, y, 'o-')
ax[0].set_xlabel('x')
ax[0].set_xlim([-1.2, 1.2])
ax[0].set_ylabel('y')
ax[0].set_title(f'{freq1} and {freq2}')
fig.tight_layout()
    
def if_inside(segment, point):
    sgn = np.dot(segment.begin - point, segment.end - point)
    if (sgn < 0):
        return True
    else:
        return False
    
def make_line(seg): #finds the coefficients of a straight line ax + by + c = 0 containing a segment seg
    a = seg.begin[1] - seg.end[1]
    b = seg.end[0] - seg.begin[0]
    c = -(a * seg.begin[0] + b * seg.begin[1])
    return np.array([a, b, c], dtype=float)

NO_INTERSECT = -2

def intersect(seg1, seg2): 
    line1 = make_line(seg1)
    line2 = make_line(seg2)
    #finding intersection of lines
    if (line1[0] * line2[1] == line1[1] * line2[0]): #lines are parallel
        return NO_INTERSECT
    matr = np.array([line1[:2], line2[:2]], dtype=float)
    koef = np.array([-line1[2], -line2[2]], dtype=float)
    point = np.linalg.solve(matr, koef)
    if (if_inside(seg1, point) and if_inside(seg2, point)):
        return point[0]
    else:
        return NO_INTERSECT

points = []
for i in range(count):
    for j in range(i + 2, count):
        pt = intersect(segs[i], segs[j])
        if (pt != NO_INTERSECT):
            points.append(pt)

num_bins = 50
n, bins, patches = ax[1].hist(points, num_bins)

ax[1].set_xlabel('Smarts')
ax[1].set_xlim([-1.2, 1.2])
ax[1].set_ylabel('Probability density')
ax[1].set_title(f'distribution of intersection for {freq1} and {freq2}')

fig.tight_layout()
fig.savefig(f'week2/inters_distr_{freq1}_and_{freq2}.png')