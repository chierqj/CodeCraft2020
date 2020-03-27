f1 = open("../data/result.txt")
f2 = open("../data/answer.txt")
line1 = f1.readlines()
line2 = f2.readlines()

ac, tol = 0, 0

for x1, x2 in zip(line1, line2):
    if x1 == x2:
        ac += 1
    tol += 1

print("@ acc: ", ac / tol)
