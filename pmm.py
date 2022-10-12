import math
import numpy as np
alphabet  ='абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
alp_dict = {}
for i in range(len(alphabet)):
    alp_dict[alphabet[i]] = i
# def alpha_pass(pi, B, observation, T, A):
#     alphas= []
#     alpha0 = []
#     c0 = 0
#     index = alphabet.index(observation[0])
#     for i in range(N):
#         alphai = pi[i]*B[i][index]
#         alpha0.append(alphai)
#         c0+= alphai
#     for i in range(len(alpha0)):
#         alpha0[i] /= c0
#     alphas.append(alpha0)
#     for t in range(1, T):
#         ct = 0
#         alphat = []
#         for i in range(N):
#             alphati = 0
#             for j in range(N):
#                 alphati += alphas[t-1][j]*A[i][j]
#             index = alphabet.index(observation[t])
#             alphati *= B[i][index]
#             ct += alphati
#             alphat.append(alphati)
#         for i in range(N):
#             alphat[i] /= ct
N = 2
def alpha_pass_tr(pi, B, observation, T, A, N):
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
    for i in range(len(alpha0)):
        alpha0[i] *= c0
    alphas.append(alpha0)
    for t in range(1, T):
        ct = 0
        alphat = []
        for i in range(N):
            alphati = 0
            for j in range(N):
                alph =  alphas[t-1][j]*A[j][i]
                alphati += alph
            alphati *= B[i][observation[t]]
            ct += alphati
            alphat.append(alphati)
        alphas.append(alphat)
        c_t = 1 / ct
        c_arr.append(c_t)
        for i in range(N):
            alphat[i] *= c_t
    return alphas, c_arr
def alpha_pass(pi, B, observation, T, A, N):
    alphas= []
    alphas_0 = []
    c_arr = []
    c0 = 0
    index = alp_dict[observation[0]]
    for i in range(N):
        alpha_0i = pi[i] * B[i][index]
        alphas_0.append(alpha_0i)
        c0+= alpha_0i
    c0 = 1/c0
    c_arr.append(c0 )
    for i in range(N):
        alphas_0[i] *= c0
    alphas.append(alphas_0)
    for t in range(1, T):
        c_t = 0
        alphas_t = []
        for i in range(N):
            alpha_ti = 0
            for j in range(N):
                alph =  alphas[t-1][j]*A[j][i]
                alpha_ti += alph
            index = alp_dict[observation[t]]
            alpha_ti *= B[i][index]
            c_t += alpha_ti
            alphas_t.append(alpha_ti)
        alphas.append(alphas_t)
        for i in range(N):
            alphas_t[i] /= c_t
        c_arr.append(1/c_t)
    return alphas, c_arr

def betta_pass_tr(c_arr, T, N, A, B, observation):
    bettas = [[]]*(T)
    bettas_T  = []
    for i in range(N):
        betta_T_1 = c_arr[T-1]
        #betta_T_1 =1
        bettas_T.append(betta_T_1)
    bettas[-1] = bettas_T
    for t in range (T-2, -1, -1):
        bettat = []
        for i in range(N):
            bettai = 0
            for j in range(N):
                b_temp =A[i][j]*B[j][observation[t+1]]*bettas[t+1][j]
                bettai += b_temp
            bettat.append(bettai*c_arr[t])
            bettat.append(bettai)
        bettas[t] = bettat
    return bettas

def betta_pass(c_arr, T, N, A, B, observation):
    bettas = [[]]*(T)
    bettas_T  = []
    for i in range(N):
        betta_T_1 = c_arr[T-1]
        bettas_T.append(betta_T_1)
    bettas[-1] = bettas_T
    for t in range (T-2, -1, -1):
        bettas_t = []
        index = alp_dict[observation[t+1]]
        for i in range(N):
            betta_ti = 0
            for j in range(N):
                b_temp =A[i][j]*B[j][index]*bettas[t+1][j]
                betta_ti += b_temp
            betta_ti *= c_arr[t]
            bettas_t.append(betta_ti)
        bettas[t] = bettas_t
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
                d_temp = alphas[t][i]*A[i][j]*B[j][v[t + 1]]*bettas[t+1][j]
                denom += d_temp
        gammas_ti = []
        for i in range(N):
            gamma_ti = 0
            for j in range(N):
                gamma_ij = alphas[t][i]*A[i][j]*B[j][v[t + 1]]*bettas[t+1][j]
                gamma_ij /= denom
                di_gammas_t[i][j] =gamma_ij
                gamma_ti += gamma_ij
            gammas_ti.append(gamma_ti)
        gammas.append(gammas_ti)
        di_gammas.append(di_gammas_t)
    denom = 0
    gammasT_1 = []
    for i in range(N):
        denom += alphas[T-1][i]
    for i in range(N):
        gamma_T_1_i = alphas[T-1][i] / denom
        gammasT_1.append(gamma_T_1_i)
    gammas.append(gammasT_1)
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

def viterbi1(A, B, pi, observation, N, T):
    deltas = []
    phis = []
    deltas_0 = []
    for i in range(N):
        index = alp_dict[observation[0]]
        delta_0i = math.log(pi[i]*B[i][index])
        deltas_0.append(delta_0i)
    deltas.append(deltas_0)
    phis.append([0,0])
    for t in range(1, T):
        deltas_t = []
        phis_t = []
        index = alp_dict[observation[t]]
        for i in range(N):
            delta_t_i_arr = []
            for j in range(N):
                #delta_t_i_j = deltas[t-1][j] + math.log(A[i][j]) + math.log(B[i][index])
                delta_t_i_j = deltas[t - 1][j] *A[i][j]**B[i][index]
                delta_t_i_arr.append(delta_t_i_j)
            delta_t_i = max(delta_t_i_arr)
            phi_t_i = delta_t_i_arr.index(delta_t_i)
            deltas_t.append(delta_t_i)
            phis_t.append(phi_t_i)
        deltas.append(deltas_t)
        phis.append(phis_t)
    arr_x = [0]*T
    delta_star = max(deltas[T-1])
    phi_star = deltas[T-1].index(delta_star)
    arr_x[T-1] = phi_star
    for t in range(T-2, -1, -1):
        x_t_star = phis[t+1][arr_x[t+1]]
        arr_x[t] = x_t_star
    return arr_x

def viterbi_edit(A, B, pi, v, N, T):
    deltas = []
    phis = []
    deltas_0 = []
    for i in range(N):
        #index = alp_dict[v[0]]
        delta_0i = pi[i]*B[i][v[0]]
        deltas_0.append(delta_0i)
    deltas.append(deltas_0)
    phis.append([110,1110,200])
    for t in range(1, T):
        deltas_t = []
        phis_t = []
        #index = alp_dict[v[t]]
        for i in range(N):
            delta_t_i_arr = []
            for j in range(N):
                #delta_t_i_j = deltas[t-1][j] + math.log(A[i][j]) + math.log(B[i][index])
                delta_t_i_j = deltas[t - 1][j] *A[j][i]
                delta_t_i_arr.append(delta_t_i_j)
            delta_t_i = max(delta_t_i_arr)
            phi_t_i = delta_t_i_arr.index(delta_t_i)
            deltas_t.append(delta_t_i*B[i][v[t]])
            phis_t.append(phi_t_i)
        deltas.append(deltas_t)
        phis.append(phis_t)
    arr_x = [0]*T
    delta_star = max(deltas[T-1])
    phi_star = deltas[T-1].index(delta_star)
    arr_x[T-1] = phi_star
    for t in range(T-2, -1, -1):
        x_t_star = phis[t+1][arr_x[t+1]]
        arr_x[t] = x_t_star
    return arr_x


def learn():
    eps = 0.01
    N = 2
    M = 34
    f =open( 'text.txt', encoding='UTF')
    observation = list(f)
    observation = [x for x in observation[0]][:5000]
    T = len(observation)
    v = []
    for i in observation:
        v.append(alp_dict[i])
    A = [[1 / N + eps, 1 / N - eps], [1 / N - eps, 1 / N + eps]]
    pi = [1/N + eps, 1/N -eps]
    B = []
    alphabet = 'абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
    alphabet = [x for x in alphabet]
    for i in range(N):
        b = []
        for j in range(M-1):
            b.append(1/M - eps)
        b.append(1/M + (M-1)*eps)
        B.append(b)
    minIters = 100
    e = 0.01
    iters= 0
    oldLogProb = -100000
    for i in range(minIters):
        alphas, c_arr = alpha_pass_tr(pi, B, v, T, A, N)
        bettas = betta_pass_tr(c_arr, T, N, A, B, v)
        gammas, di_gammas = get_gammas(T,A,B, alphas, bettas, v, N)
        pi, A, B = reestimate(gammas, di_gammas, N, M, A, B, T, v)
    # v = []
    # print(B)
    # for i in observation:
    #     v.append(alp_dict[i])
    # v = np.array(v)
    arr_x= viterbi_edit(A, B, pi, v, N, T)
    return pi, A, B, arr_x


# observation = [1,3,2,1,3]
# classes = [1,2,3]
# v = []
# for i in observation:
#     v.append(classes.index(i))
# A = [[0.2, 0.3, 0.5], [0.2, 0.2, 0.6], [0, 0.2, 0.8]]
# B = [[0.7, 0.2, 0.1], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]
# mu = [0.05, 0.2, 0.75]
# alphas, _= alpha_pass_tr(mu, B, v.copy(), 5, A,3)
# # alphas = [[0.035, 0.06, 0], [0.0019, 0.04815, 0.00075], [0.000346, 0.00462, 0.004352],
# #           [0.00069, 0.00032, 0], [0.00002,0.000081, 0.0006459]]
# bettas = betta_pass_tr(5, 3, A, B, v.copy())
# gammas, di_gammas = get_gammas(5, A, B, alphas, bettas, v, 3)
# pis, A, B = reestimate(gammas, di_gammas, 3, 3, A, B, 5, v)
# print(pis)
# print(A)
# print(B)
# print(alphas)
# viterbi_edit(A, B,mu, v.copy(), 3, 5)

# print(bettas)
# betta_pass
pi, A, B, arr_x = learn()
print(A)
print(B)
f = open('text.txt', encoding='UTF')
observation = list(f)
observation = [x for x in observation[0]][:5000]
# print(arr_x)
group1 = []
group2 = []
for i in range(len(arr_x)):
    if arr_x[i] == 0 and observation[i] not in group1 and observation[i] not in group2:
        group1.append(observation[i])
    elif  arr_x[i] == 1 and observation[i] not in group2 and observation[i] not in group1:
        group2.append(observation[i])
print(group1)
print(group2)

