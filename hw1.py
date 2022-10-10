import random

import scipy.stats as st
from scipy.stats import rv_continuous
import numpy as np
theta = {'x1' : 5/11, 'x2': 6/11, 'mu1': 4.5, 'mu2':10.1 }
y1 = np.random.normal(theta['mu1'], 1)
y2 = np.random.normal(theta['mu2'], 1)
