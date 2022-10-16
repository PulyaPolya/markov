import math
import numpy as np
alphabet  ='абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
alp_dict = {}
for i in range(len(alphabet)):
    alp_dict[alphabet[i]] = i

def alpha_pass_tr(pi, B, observation, T, A, N, c):
    alphas= []
    alpha0 = []
    c_arr = []
    c0 = 0
    for i in range(N):
        alphai = pi[i]*B[i][observation[0]]
        alpha0.append(alphai)
        c0+= alphai
    c0 = 1 / c0
    c_arr.append(c0)
    if c:
        for i in range(N):
            alpha0[i] *= c0
    alphas.append(alpha0)
    for t in range(1, T):
        ct = 0
        alphat = []
        for i in range(N):
            alphati = 0
            for j in range(N):
                alphati +=  alphas[t-1][j]*A[j][i]
            alphati *= B[i][observation[t]]
            ct += alphati
            alphat.append(alphati)

        c_t = 1 / ct
        c_arr.append(c_t)
        if c:
            for i in range(N):
                alphat[i] *= c_t
        alphas.append(alphat.copy())
    return alphas, c_arr


def betta_pass_tr(c_arr, T, N, A, B, observation, c):
    bettas = [[]]*(T)
    bettas_T  = []
    for i in range(N):
        if c:
            betta_T_1 = c_arr[T-1]
        else:
            betta_T_1 =1
        bettas_T.append(betta_T_1)
    bettas[-1] = bettas_T.copy()
    for t in range (T-2, -1, -1):
        bettat = []
        sum = 0
        for i in range(N):
            bettai = 0
            for j in range(N):
                b_temp =A[i][j]*B[j][observation[t+1]]*bettas[t+1][j]
                bettai += b_temp
            sum+= bettai
            if c:
                bettat.append(bettai)
            else:
                bettat.append(bettai)
        for i in range(N):
            bettat[i] /= sum
        bettas[t] = bettat.copy()
    return bettas


def get_gammas(T, A, B, alphas, bettas, v, N):
    gammas = []
    di_gammas = []
    for t in range(T-1):
        #di_gammas_t = [[0] * N] * N
        di_gammas_t = []
        for i in range(N):
            di_gammas_t.append([0]*N)
        denom = 0
        for i in range(N):
            for j in range(N):
                denom +=alphas[t][i]*A[i][j]*B[j][v[t + 1]]*bettas[t+1][j]

        gammas_ti = []
        for i in range(N):
            gamma_ti = 0
            for j in range(N):
                gamma_ij = alphas[t][i]*A[i][j]*B[j][v[t + 1]]*bettas[t+1][j]
                gamma_ij /= denom
                di_gammas_t[i][j] =gamma_ij
                gamma_ti += gamma_ij
            gammas_ti.append(gamma_ti)
        gammas.append(gammas_ti.copy())
        di_gammas.append(di_gammas_t.copy())
    denom = 0
    gammasT_1 = []
    for i in range(N):
        denom += alphas[T-1][i]
    for i in range(N):
        gamma_T_1_i = alphas[T-1][i] / denom
        gammasT_1.append(gamma_T_1_i)
    gammas.append(gammasT_1.copy())
    return gammas, di_gammas

def reestimate(gammas, di_gammas, N, M, A, B, T, v):
    pis = []
    for i in range (N):
        pis.append(gammas[0][i])
    for i in range(N):
        for j in range(N):
            numer = 0
            denom = 0
            for t in range(T-1):
                numer += di_gammas[t][i][j]
                denom += gammas[t][i]
            A[i][j] = numer/ denom
    for i in range(N):
        for j in range(M):
            numer =0
            denom =0
            for t in range(T):
                if v[t] == j:
                    numer += gammas[t][i]
                denom += gammas[t][i]
            B[i][j] =  numer/denom
    return pis, A, B

def compute_probs(T, c_arr):
    logProb = 0
    for i in range(T):
        logProb += math.log2(c_arr[i])
    logProb = - logProb
    return logProb


def viterbi_st(A, pi, B, v):
    I = A.shape[0]    # Number of states
    N = len(v)  # Length of observation sequence

    # Initialize D and E matrices
    D = np.zeros((I, N))
    E = np.zeros((I, N-1)).astype(np.int32)
    D[:, 0] = np.multiply(pi, B[:, v[0]])

    # Compute D and E in a nested loop
    for n in range(1, N):
        for i in range(I):
            temp_product = np.multiply(A[:, i], D[:, n-1])
            D[i, n] = np.max(temp_product) * B[i, v[n]]
            E[i, n-1] = np.argmax(temp_product)

    # Backtracking
    S_opt = np.zeros(N).astype(np.int32)
    S_opt[-1] = np.argmax(D[:, -1])
    for n in range(N-2, -1, -1):
        S_opt[n] = E[int(S_opt[n+1]), n]

    return S_opt, D, E

def learn(N, M,v, A , B, pi, minIters, c):
    T = len(v)
    for i in range(minIters):
        alphas, c_arr = alpha_pass_tr(pi, B, v, T, A, N, c)
        bettas = betta_pass_tr(c_arr, T, N, A, B, v, c)
        gammas, di_gammas = get_gammas(T,A,B, alphas, bettas, v, N)
        pi, A, B = reestimate(gammas, di_gammas, N, M, A, B, T, v)
        #print(A)
    arr_x, _, _ = viterbi_st(np.array(A), np.array(pi), np.array(B), np.array(v))
    return pi, A, B, arr_x


f =open( 'text2', encoding='UTF-8')
observation = list(f)
N = 2
observation = [x for x in observation[0]][:200]
T = len(observation)
v = []
for i in observation:
    v.append(alp_dict[i])
eps= 0.15
A = [[1 / N - eps, 1 / N + eps], [1 / N + eps, 1 / N - eps]]
M = len(alphabet)
eps1 = eps
pi = [1/N - eps1, 1/N +eps1]
B = []
for i in range(N):
    b = []
    for i in range(M):
        b.append(1/M)
    B.append(b)
alphabet = 'абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
alphabet = [x for x in alphabet]

eps1 = 0.01
b = (1,)*M
pi, A, B, arr_x = learn(N, M,v, A , B, pi, 10, c = True)
probs = {}
for i in alphabet:
    probs[i] = [0,0]
for j in range(len(observation)):
    probs[observation[j]][arr_x[j]] +=1
print(probs)
group1 = []
group2 = []
group3 = []
for i in probs:
    if probs[i][0] > probs[i][1]:
        group1.append(i)
    elif probs[i][0] == probs[i][1]:
        group2.append(i)
    else:
        group3.append(i)
print(group1)
print(group2)
print(group3)


print(pi)
print(A)
print(B)

# print(list(arr_x))