#
#
import math

import math
observation = [1,3,2,1,3]
A = [[0.2, 0.3, 0.5], [0.2, 0.2, 0.6], [0, 0.2, 0.8]]
B = [[0.7, 0.2, 0.1], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]
mu = [0.05, 0.2, 0.75]
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
    c0 = 0
    #index = alphabet.index(observation[0])
    for i in range(N):
        alphai = pi[i]*B[i][observation[0]-1]
        alpha0.append(alphai)
        c0+= alphai
    # for i in range(len(alpha0)):
    #     alpha0[i] /= c0
    alphas.append(alpha0)
    for t in range(1, T):
        ct = 0
        alphat = []
        for i in range(N):
            alphati = 0
            for j in range(N):
                alph =  alphas[t-1][j]*A[j][i]
                alphati += alph
            #index = alphabet.index(observation[t])
            alphati *= B[i][observation[t]-1]
            ct += alphati
            alphat.append(alphati)
        alphas.append(alphat)
        # for i in range(N):
        #     alphat[i] /= ct
    return alphas
def alpha_pass(pi, B, observation, T, A, N):
    alphas= []
    alphas_0 = []
    c0 = 0
    index = alp_dict[observation[0]]
    for i in range(N):
        alpha_0i = pi[i] * B[i][index]
        alphas_0.append(alpha_0i)
        c0+= alpha_0i
    for i in range(N):
        alphas_0[i] /= c0
    alphas.append(alphas_0)
    for t in range(1, T):
        c_t = 0
        alphas_t = []
        for i in range(N):
            alpha_ti = 0
            for j in range(N):
                alph =  alphas[t-1][j]*A[j][i]
                alpha_ti += alph
            index = alp_dict[observation[0]]
            alpha_ti *= B[i][index]
            c_t += alpha_ti
            alphas_t.append(alpha_ti)
        alphas.append(alphas_t)
        for i in range(N):
            alphas_t[i] /= c_t
    return alphas

def betta_pass(cT_1, T, N, A, B, observation):
    bettas = [[]]*(T)
    bettas_T  = []
    for i in range(N):
        betta_T_1 = 1
        bettas_T.append(betta_T_1)
    bettas[-1] = bettas_T
    for t in range (T-2, -1, -1):
        bettat = []
        for i in range(N):
            bettai = 0
            for j in range(N):
                b_temp =A[i][j]*B[j][observation[t+1]-1]*bettas[t+1][j]
                bettai += b_temp
            bettat.append(bettai)
        bettas[t] = bettat
    return bettas
def get_gammas(T, alphas, bettas):
    gammas = []
    gammas_pairs = []
    for t in range(T-1):
        gammas_pairs_t = [[0] * N] * N
        denom = 0
        for i in range(N):
            for j in range(N):
                d_temp = alphas[t][i]*A[i][j]*B[j][observation[t+1]-1]*bettas[t+1][j]
                denom += d_temp
        gammas_ti = []
        for i in range(N):
            gamma_ti = 0
            for j in range(N):
                gamma_ij = (alphas[t][i]*A[i][j]*B[j][observation[t+1]-1]*bettas[t+1][j]) / denom
                gammas_pairs_t[i][j] =gamma_ij
                gamma_ti += gamma_ij
            gammas_ti.append(gamma_ti)
        gammas.append(gammas_ti)
        gammas_pairs.append(gammas_pairs_t)
    denom = 0
    gammasT_1i = []
    for i in range(N):
        denom += alphas[T-1][i]
    for i in range(N):
        gamma_T_1 = alphas[T-1][i] / denom
        gammasT_1i.append(gamma_T_1)
    gammas.append(gammasT_1i)
    return gammas, gammas_pairs

def reestimate(gammas,gammas_pairs, N, M, A, B, T):
    pis = []
    for i in range (N):
        pis.append(gammas[0][i])
    for i in range(N):
        for j in range(N):
            numer = 0
            denom = 0
            for t in range(T-1):
                numer += gammas_pairs[t][i][j]
                denom += gammas[t][i]
            A[i][j] = numer/ denom
    for i in range(N):
        for j in range(M):
            numer =0
            denom =0
            for t in range(T):
                if observation[t]-1 == j:
                    numer += gammas[t][i]
                denom += gammas[t][i]
            B[i][j] =  numer/denom
    return pis, A, B

def compute_probs(T, c):
    logProb = 0
    for i in range(T):
        logProb += math.log2(c[i])
    logProb = - logProb
    return logProb

def viterbi(A, B, pi, observation, N, T):
    deltas = []
    phis = []
    deltas_0 = []
    for i in range(N):
        index = alp_dict[observation[0]]
        delta_0i = math.log2(pi[i]*B[i][index])
        deltas_0.append(delta_0i)
    for t in range(1, T):
        deltas_t = []
        phis_t = []
        index = alp_dict[observation[t]]
        for i in range(N):
            delta_t_i_arr = []
            for j in range(N):
                delta_t_i_j = deltas[t-1][j] + math.log2(A[i][j]) + math.log(B[i][index])
                delta_t_i_arr.append(delta_t_i_j)
            delta_t_i = max(delta_t_i_arr)
            phi_t_i = delta_t_i_arr.index(delta_t_i)
            deltas_t.append(delta_t_i)
            phis_t.append(phi_t_i)
        deltas.append(deltas_t)
        phis.append(phis_t)




def learn():
    eps = 0.01
    N = 2
    M = 34
    f =open( 'text.txt')
    observation = list(f)
    observation = [x for x in observation[0]]
    T = len(observation)
    A = [[1 / N + eps, 1 / N - eps], [1 / N + eps, 1 / N - eps]]
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

#alpha_pass(mu, B, observation.copy(), 5, A,N)
# betta_pass
learn()