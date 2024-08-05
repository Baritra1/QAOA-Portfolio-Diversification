import numpy as np

q=2
#bitstring="10001000000100011001"
# bitstring="11111110111111110111"
bitstring="00011001010100111111"
correlation_matrix=np.loadtxt("correlationmatrix_size4.txt")
maxEdgeWeight=0
for i in range(len(correlation_matrix)):
    for j in range(len(correlation_matrix[0])):
        maxEdgeWeight=max(maxEdgeWeight,correlation_matrix[i][j])
n=4

# score = 0  #Score for this bitstring
# for i in range(n):
#     for j in range(n):
#         if bitstring[n*i+j]=='1':
#             score+=correlation_matrix[i][j]

# #penalty for violating constraints
# penalty=3*maxEdgeWeight
# substrings=[bitstring[i:i + n] for i in range(0, len(bitstring)-n, n)]
# clusters=bitstring[n**2:n**2+n]

# #apply penalty for wrong number of clusters

# count=0
# for k in clusters:
#     if(k=='1'):
#         count+=1
# if count!=q:
#     score-=penalty
#     print("grra")

# #apply penalty for multiple similar stock
# for j in substrings:
#     count=0
#     for k in j:
#         if(k=='1'):
#             count+=1
#     if count!=1:
#         score-=penalty
#         print("grrb")

# #apply penalty for x_{ij}>y_j
# for i in range(n):
#     for j in range(n):
#         if bitstring[n*i+j]>bitstring[n**2+j]:
#             score-=penalty
#             print("grrc")


# for j in range(n):
#     if(bitstring[n*j+j]!=bitstring[n**2+j]):
#         score-=penalty
#         print("grrd")
print(score)
