A = [[0.2, 0.3, 0.5], [0.2, 0.2, 0.6], [0, 0.2, 0.8]]
B = [[0.7, 0.2, 0.1], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]
mu = [0.05, 0.2, 0.75]
y = [0, 2, 1, 0, 2]
Alpha = []
Alpha0 = []
for i in range(3):
    alpha0 = mu[i]*B[i][y[0]]
    Alpha0.append(alpha0)
Alpha.append(Alpha0)
for i in range (1,5):
    Alpha_i = []
    for j in range(3):
        sum = 0
        for k in range(3):
            sum+= Alpha[i-1][k]*A[k][i-1]*B[j][y[i]]
        Alpha_i.append(sum)
    Alpha.append(Alpha_i)
