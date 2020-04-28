import networkx as nx
import matplotlib.pyplot as plt
import random

# WS network
# 生成20个节点的小世界网络
# 每个节点有四个邻居
# 每条连边以0.3的概率随机重置链接
WS = nx.random_graphs.watts_strogatz_graph(30000, 10, 0.5)
print(len(WS.edges))
# pos = nx.circular_layout(WS)
# nx.draw(WS, pos, with_labels = False, node_size = 30)
# plt.show()

file = "./gen/test_data.txt"
f = open(file, "w")
for it in WS.edges:
    x = random.randint(5, 300)
    f.write("{},{},{}\n".format(it[0], it[1], x))
    f.write("{},{},{}\n".format(it[1], it[0], x))
