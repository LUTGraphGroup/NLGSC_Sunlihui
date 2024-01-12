import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy.stats import kendalltau


def distance(G):
    length = dict(nx.all_pairs_shortest_path_length(G))
    # print(length)
    distance_matrix = np.zeros([len(G), len(G)])
    for key1, value1 in length.items():
        for key2, value2 in value1.items():
            distance_matrix[key1][key2] = length[key1][key2]
    return distance_matrix



def k_shell(G):
    """计算每个节点的KS值"""
    G1 = G.copy()  # 目的是对G1进行删点和边的操作，对G没有影响

    def k_shell_1(G1):
        importance_dict = {}
        level = 1
        while len(G1.degree):
            importance_dict[level] = []
            while True:
                level_node_list = []
                for item in G1.degree:
                    if item[1] <= level:
                        level_node_list.append(item[0])
                G1.remove_nodes_from(level_node_list)  # 从G中移除节点，移除完后为空，导致后续函数调用G报列表索引越界，k_sheel(G)放到最后
                importance_dict[level].extend(level_node_list)
                if not len(G1.degree):
                    return importance_dict
                if min(G1.degree, key=lambda x: x[1])[1] > level:
                    break
            level = min(G1.degree, key=lambda x: x[1])[1]
        # print('importance_dict',importance_dict)
        return importance_dict

    a = k_shell_1(G1)
    # print('a',a)
    H = {}
    for x, y in a.items():
        for z in y:
            H[z] = x
    # print('H',H)
    H_reverse = sorted(H.items(), key=lambda x: x[0])
    # print(dict(H_reverse))
    KS1 = list(dict(H_reverse).values())
    # print('KS1',KS1)
    return KS1




def calculate_Ei_values(G):
    Ei_values = []  # 用于存储计算的Ei值

    for a in np.arange(0, 1.1, 0.1):
        LC_values = []  # 存储每个节点的LC值

        # 计算每个节点的LC值
        for i in G.nodes():
            di = G.degree(i)
            neighbor_nodes = list(G.neighbors(i))
            neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
            neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
            T = nx.triangles(G, i)

            LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                    6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                       neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                           (1 - a) ** 3) * T

            LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的Ai值

        Ai_values = []
        ks = k_shell(G)  # 节点的K壳值
        for i in range(len(G)):
            neighbor_nodes_i = list(G.neighbors(i))
            Ai = sum(ks[j] for j in neighbor_nodes_i)
            Ai_values.append(Ai)
        # print(Ai_values)

        # 计算每个节点的GI值
        dis = distance(G)  #节点之间的最短距离
        GI_values = []
        for i in range(len(G)):
            SS = 0
            for j in range(len(G)):
                if i != j and dis[i][j] != 0:
                    SS += LC_values[j] / dis[i][j]
            GI_values.append(SS)
        # print('GI1',GI)
        # print(GI_values)


        # 计算每个节点的Ei值
        Ei_values_for_a = []
        ks = k_shell(G)
        for i in range(len(G)):
            sum1 = LC_values[i]+sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)) )
            Ei = GI_values[i] * sum1
            Ei_values_for_a.append(Ei)

        Ei_values.append(Ei_values_for_a)
    return Ei_values


def calculate_Ei_values_USAir(G):
    LC_values = []  # 存储每个节点的LC值
    a=0.1
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                    6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                       neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                           (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
        # print(Ai_values)


        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i]+sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)) )
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_Dolphins(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.3
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
        # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)



    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_Jazz(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.3
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
        # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values




def calculate_Ei_values_Email(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.6
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_Stelzl(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.4
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_Hamster(G):
    LC_values = []  # 存储每个节点的LC值
    a = 1.0
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_Facebook(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.1
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values

def calculate_Ei_values_PGP(G):
    LC_values = []  # 存储每个节点的LC值
    a = 1.0
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values


def calculate_Ei_values_Yeast(G):
    LC_values = []  # 存储每个节点的LC值
    a = 0.4
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                 neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                     (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)

        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i] + sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)))
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values




def calculate_Ei_values_Power(G):
    LC_values = []  # 存储每个节点的LC值
    a=0.4
    # 计算每个节点的LC值
    for i in G.nodes():
        di = G.degree(i)
        neighbor_nodes = list(G.neighbors(i))
        neighbor_degrees_sum = sum(G.degree(j) for j in neighbor_nodes)
        neighbor_degree_sum_square = sum([G.degree(j) ** 2 for j in neighbor_nodes])
        T = nx.triangles(G, i)

        LC = (a ** 3) * (di ** 3) + 3 * a * ((1 - a) ** 2) * (di ** 2) + (
                    6 * (a ** 2) - 3 * a - 2 * (a ** 3)) * di + 3 * (a ** 3) * (
                       neighbor_degree_sum_square) + 3 * a * (a ** 2 - 4 * a + 2) * neighbor_degrees_sum + 6 * (
                           (1 - a) ** 3) * T

        LC_values.append(LC)
        # print(LC_values)

        # 计算每个节点的GI值
    dis = distance(G)
    GI_values = []
    for i in range(len(G)):
        SS = 0
        for j in range(len(G)):
            if i != j and dis[i][j] != 0:
                SS += LC_values[j] / dis[i][j]
        GI_values.append(SS)
    # print(GI_values)

    # 计算每个节点的Ai值
    Ai_values = []
    ks = k_shell(G)
    for i in range(len(G)):
        neighbor_nodes_i = list(G.neighbors(i))
        Ai = sum(ks[j] for j in neighbor_nodes_i)
        Ai_values.append(Ai)
    #     # print(Ai_values)


        # 计算每个节点的Ei值
    Ei_values = []  # 用于存储计算的Ei值
    ks = k_shell(G)
    for i in range(len(G)):
        sum1 = LC_values[i]+sum(LC_values[j] * (ks[i] / Ai_values[j]) for j in list(G.neighbors(i)) )
        Ei = GI_values[i] * sum1
        Ei_values.append(Ei)

    # print(Ei_values)
    return Ei_values




def sckendall_1(a, b):
    kendall_tau_2, p_value = kendalltau(a, b)  # kendalltau：系统自带肯德尔系数
    # kendall_tau = sckendall(a, b)#sckendall：自己定义的肯德尔系数，好像
    # print(kendall_tau_2)
    # print(kendall_tau)
    return kendall_tau_2
# generalized_centrality_parameter(G, A)