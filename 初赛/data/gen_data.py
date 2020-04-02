import random

n = 100
file = "./gen1/test_data.txt"

edges = []
for u in range(1, n):
    x = random.randint(5, 15)
    for i in range(x):
        v = random.randint(1, n)
        w = random.randint(5, 300)
        edges.append((u, v, w))
random.shuffle(edges)
f = open(file, "w")
for it in edges:
    f.write("{},{},{}\n".format(it[0], it[1], it[2]))
