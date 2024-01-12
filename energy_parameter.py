import pandas as pd
import random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import kendalltau
import xlrd as xd
from get_eneralized import *

# '''读取MATLAB中的邻接矩阵Excel文件'''
G = nx.read_gml("data/Dolphins.gml", label='id')
A = nx.to_scipy_sparse_array(G).todense()#构造邻接矩阵
G = nx.to_networkx_graph(A)


# MM = pd.read_excel("data/Dolphins.xlsx",header=None)
# N = nx.from_numpy_matrix(np.array(MM))
# A = nx.to_scipy_sparse_array(N).todense()#构造邻接矩阵
# G = nx.to_networkx_graph(A)


H_1 =calculate_Ei_values(G)

data = xd.open_workbook('data/data_set_1.xls')  # 打开excel表所在路径
sheet = data.sheet_by_name('Dolphins')  # 读取数据，以excel表名来打开
SIR = []
for r in range(sheet.ncols):  # 将表中数据按列逐步添加到列表中，最后转换为list结构
    data1 = []
    for c in range(sheet.nrows):
        data1.append(sheet.cell_value(c, r))
    SIR.append(list(data1))


SIR_matrix = np.zeros([len(H_1), len(SIR)])
for i in range(len(H_1)):
    for j in range(len(SIR)):
        SIR_matrix[i][j]=sckendall_1(H_1[i],SIR[j])
# print(SIR_matrix)
result = pd.DataFrame(SIR_matrix)


# 声明一个读写对象
writer = pd.ExcelWriter("out/Parameter experiment1.xlsx",engine='openpyxl',mode='a')
result.to_excel(writer,sheet_name='Dolphins',index=False)
# # writer.save()#保存读写的内容
writer.close()