import numpy as np

fd = open("OutputSimulation1.txt")

serviceSats = list()
tokens = list()
for line in fd:
    tokens.append(line.strip().split(": "))
fd.close()
tokens.pop(len(tokens)-1)
for token in tokens:
    serviceSats.append(int(token[1]))
count0 = 0
for i in serviceSats:
    if i == 0:
        count0 += 1
#   the outage probability calculation
print("Outage probability is: {}".format(count0/len(serviceSats)))