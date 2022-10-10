import numpy as np
import pandas as pd
import json
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
from numpy import asarray
import random
import time
import csv
from numpy import genfromtxt
my_data = genfromtxt('data03.csv', delimiter=';')
x = []
y = []
for arr in my_data:
  arr1 = arr[:2]
  x.append(arr1)
  if arr[-1] == 1:
    y.append(1)
  else:
    y.append(-1)
x = np.array(x)
y = np.array(y)
size_data = len(x)
# for i in range(size_data):
#   if y[i] == 1:
#     plt.scatter(x[i][0], x[i][1], color = 'red' )
#   else:
#     plt.scatter(x[i][0], x[i][1], color = 'blue' )
#
# plt.show()

w = np.array([ 3.655, -0.394, -1.   ])
wT = np.matrix.transpose(w)
new_x = []
for xi in x:
  xi = np.array(xi.tolist() + [1,])
  new_x.append(xi)
new_x = np.array(new_x)
job = 'not_done'
while True:
  if job == 'done':
    break
  m = 0
  for xi, yi in zip(new_x, y):
    if yi * np.matmul(wT, xi) <= 0:
      wT += yi*xi
      m += 1
  if m == 0:
    job = 'done'
    break

print(w)