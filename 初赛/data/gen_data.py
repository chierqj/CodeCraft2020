import random

n, m, scc = 1000000, 280000, 8000
file = "./gen/test_data.txt"


def get_mid_point():
    block = int(n / scc)
    start = 0
    ans = []
    for i in range(scc):
        rd = random.randint(start, start+block-1)
        ans.append((rd, start, start+block-1))
        start += block
    return ans


def generate_data(mid_points):
    count = int(m / scc)
    print(count)
    ans = []
    st = set()
    st.add((-1, -1))
    for it in mid_points:
        rd = random.randint(0, min(10, count))
        st.add((it[0], -1))
        st.add((-1, it[0]))

        # [mid, x] 任意rd条边
        for _ in range(rd):
            x = -1
            while (it[0], x) in st:
                x = random.randint(it[1], it[2])
            st.add((it[0], x))
            ans.append((it[0], x, random.randint(5, 300)))

        # [l, r]区间内随机count-rd条边
        for _ in range(count-rd):
            u, v = -1, -1
            while (u, v) in st:
                u = random.randint(it[1], it[2])
                v = random.randint(it[1], it[2])
            st.add((u, v))
            ans.append((u, v, random.randint(5, 300)))
        st.add((it[0], -1))
    sz = len(mid_points)
    while len(ans) < m:
        u = mid_points[random.randint(0, sz-1)][0]
        v = mid_points[random.randint(0, sz-1)][0]
        if u == v or (u, v) in st:
            continue
        st.add((u, v))
        ans.append((u, v, random.randint(5, 300)))
    return ans


def save(edges):
    f = open(file, "w")
    for it in edges:
        f.write("{},{},{}\n".format(it[0], it[1], it[2]))


mid = get_mid_point()
edges = generate_data(mid)
save(edges)
