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
        if c:
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

def betta_pass_tr(c_arr, T, N, A, B, observation, c):
    bettas = [[]]*(T)
    bettas_T  = []
    for i in range(N):
        if c:
            betta_T_1 = c_arr[T-1]
        else:
            betta_T_1 =1
        bettas_T.append(betta_T_1)
    bettas[-1] = bettas_T
    for t in range (T-2, -1, -1):
        bettat = []
        for i in range(N):
            bettai = 0
            for j in range(N):
                b_temp =A[i][j]*B[j][observation[t+1]]*bettas[t+1][j]
                bettai += b_temp
            if c:
                bettat.append(bettai*c_arr[t])
            else:
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


def viterbi_edit(A, B, pi, v, N, T):
    deltas = []
    phis = []
    deltas_0 = []
    e = 0.000000001
    for i in range(N):
        delta_0i = pi[i]*B[i][v[0]]
        #delta_0i = math.log(pi[i] * B[i][v[0]]+e)
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
                #delta_t_i_j = deltas[t-1][j] + math.log(A[i][j] + e) + math.log(B[i][v[t]] +e)
                delta_t_i_j = deltas[t - 1][j] *A[j][i]
                delta_t_i_arr.append(delta_t_i_j)
            delta_t_i = max(delta_t_i_arr)
            phi_t_i = delta_t_i_arr.index(delta_t_i)
            deltas_t.append(delta_t_i*B[i][v[t]])
            #deltas_t.append(delta_t_i)
            phis_t.append(phi_t_i)
        deltas.append(deltas_t)
        phis.append(phis_t)
    arr_x = [-9]*T
    delta_star = max(deltas[T-1])
    phi_star = deltas[T-1].index(delta_star)
    arr_x[T-1] = phi_star
    for t in range(T-2, -1, -1):
        x_t_star = phis[t+1][arr_x[t+1]]
        arr_x[t] = x_t_star
    return arr_x

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
    oldLogProb = -100000
    for i in range(minIters):
        alphas, c_arr = alpha_pass_tr(pi, B, v, T, A, N, c)
        bettas = betta_pass_tr(c_arr, T, N, A, B, v, c)
        gammas, di_gammas = get_gammas(T,A,B, alphas, bettas, v, N)
        pi, A, B = reestimate(gammas, di_gammas, N, M, A, B, T, v)
        #print(A)
    arr_x, _, _ = viterbi_st(np.array(A), np.array(pi), np.array(B), np.array(v))
    #arr_x= viterbi_edit(A, B, pi, v, N, T)
    return pi, A, B, arr_x


# observation = [1,3,2,1,3]
# classes = [1,2,3]
# v = []
# for i in observation:
#     v.append(classes.index(i))
# A = [[0.2, 0.3, 0.5], [0.2, 0.2, 0.6], [0, 0.2, 0.8]]
# B = [[0.7, 0.2, 0.1], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]
# mu = [0.05, 0.2, 0.75]
# arr_x = viterbi_edit(A, B, mu, v, 3, 5)
# pi, A, B, arr_x = learn(3, 3,v, A , B, mu, 1, c = True)


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
# A = np.random.dirichlet((3,3), 2).tolist()
# pi = np.random.dirichlet((2,2), 1).tolist()

# A = [[0.5, 0.5], [0.5, 0.5]]
# pi = [0.5, 0.5]
M = len(alphabet)
#eps1 = 0.1
pi = [1/N - eps, 1/N +eps]
B = []
for i in range(N):
    b = []
    for i in range(M):
        b.append(1/M)
    B.append(b)
alphabet = 'абвгґдеєжзиіїйклмнопрстуфхцчшщьюя '
alphabet = [x for x in alphabet]

eps1 = 0.01
# for i in range(N):
#     b = []
#     b = np.random.dirichlet((3, 3), 2).tolist()
#     # for j in range(M - 1):
#     #     b.append(1 / M - eps1)
#     # b.append(1 / M + (M - 1) * eps1)
#     B.append(b)
b = (1,)*M
#B = np.random.dirichlet(b,2).tolist()
pi, A, B, arr_x = learn(N, M,v, A , B, pi, 100, c = True)
probs = {}
for i in alphabet:
    probs[i] = [0,0]
for j in range(len(observation)):
    probs[observation[j]][arr_x[j]] +=1
print(probs)
group1 = []
group2 = []
for i in probs:
    if probs[i][0] > probs[i][1]:
        group1.append(i)
    else:
        group2.append(i)
print(group1)
print(group2)


# group1 = []
# group2 = []
# for i in range(len(arr_x)):
#     if arr_x[i] == 0 and observation[i] not in group1:
#         group1.append(observation[i])
#     elif  arr_x[i] == 1 and observation[i] not in group2:
#         group2.append(observation[i])
# print(group1)
# print(group2)

# v = [3, 0, 1]
# A = [[0.5, 0.25, 0.25], [0.3, 0.4, 0.3], [0.25, 0.25, 0.5]]
# B  =[[0.6, 0.2, 0.15, 0.05], [0.25, 0.25, 0.25, 0.25], [0.05, 0.1, 0.35, 0.5]]
# pi = [0.33, 0.33, 0.33]
# #arr_x, _, _ = viterbi_st(np.array(A), np.array(pi), np.array(B), np.array(v))
# pi, A, B, arr_x = learn(3, 3,v, A , B, pi, 1, c = True)

# A = [[0.8, 0.1, 0.1], [0.2, 0.7, 0.1], [0.1, 0.3, 0.6]]
# B= [[0.7, 0, 0.3],[0.1, 0.9, 0],[0.0, 0.2, 0.8]]
# pi = [0.6,0.2, 0.2]
# v = [0, 2, 0, 2, 2, 1]
#arr_x, _, _ = viterbi_st(np.array(A), np.array(pi), np.array(B), np.array(v))
print(pi)
print(A)
print(B)

print(list(arr_x))
