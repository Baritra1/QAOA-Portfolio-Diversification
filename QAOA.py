from qiskit import *
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import DiagonalGate
from qiskit import transpile
from qiskit_aer import Aer
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import math




#Create quantum circuit for QAOA from edges and parmeters
def QAOA(n,edges,p,betas,gammas):  
    #Define quantum and classical registers
    nQubits=n**2+n
    qr = QuantumRegister(nQubits)
    cr = ClassicalRegister(nQubits)
    circuit = QuantumCircuit(qr,cr)
    
    #Initial Hadamards
    for q in range(nQubits):
        circuit.h(q)        
    #For the number of specified iterations
    for P in range(p):

        #cost hamiltonian
        for i in range(n):
            for j in range(n):
                        qubit=n*i+j
                        circuit.z(qubit)
                        circuit.rz(phi=gammas[P]*correlation_matrix[i][j],qubit=qubit)
                        circuit.z(qubit)

        penalty=3*maxEdgeWeight
        #generate penalty unitary matrices 

        #choosing multiple most similar stocks
        entries= [1]*(2**n)
        cur=0
        count=0
        while cur<2**n:
            if(count.bit_count()!=1):
                entries[cur]=0
            cur+=n+1
            count+=1
        for i in range(len(entries)):
            entries[i]=math.e**(entries[i]*gammas[P]*penalty*1j)
        for i in range(n):
            excessSimilarStock = list(range(i*n,(i+1)*n))
            circuit.append(DiagonalGate(entries),excessSimilarStock)

        #choosing the wrong amount of clusters 
        entries= [1]*(2**n)
        cur=0
        count=0
        while cur<2**n:
            if(count.bit_count()!=q):
                entries[cur]=0
            cur+=n+1
            count+=1
        for i in range(len(entries)):
            entries[i]=math.e**(entries[i]*gammas[P]*penalty*1j)
        wrongCluster = list(range(n**2,n**2+n))
        circuit.append(DiagonalGate(entries),wrongCluster)

        #x_{jj}=y_j
        entries=[1,0,0,1]
        for i in range(len(entries)):
            entries[i]=math.e**(entries[i]*gammas[P]*penalty*1j)
        for j in range(n):
            similarCluster=[n*j+j,n**2+j]
            circuit.append(DiagonalGate(entries),similarCluster)

        # #x_{ij}<=y_j
        # entries=[1,1,0,1]
        # for i in range(len(entries)):
        #     entries[i]=math.e**(entries[i]*gammas[P]*penalty*1j)
        # for i in range(n):
        #     for j in range(n):
        #         consistencyConstraint=[n*i+j,n**2+j]
        #         circuit.append(DiagonalGate(entries),consistencyConstraint)

        #mixer hamiltonian
        for q in range(nQubits):
            circuit.rx(theta=2*betas[P],qubit=q)   
    circuit.measure(qr,cr)
    return circuit

#Compute all scores for a set of edges
def computeExpectationValue(counts,edges,n):
    totalScore = 0
    totalSamples = 0
    #For each bitstring measured (keys in counts dictionary)
    for bitstring in counts.keys():
        score = 0  #Score for this bitstring
        for i in range(n):
            for j in range(n):
                if bitstring[n*i+j]=='1':
                    score+=correlation_matrix[i][j]

        #penalty for violating constraints
        penalty=3*maxEdgeWeight
        substrings=[bitstring[i:i + n] for i in range(0, len(bitstring)-n, n)]
        clusters=bitstring[n**2:n**2+n]

        #apply penalty for wrong number of clusters

        count=0
        for k in clusters:
            if(k=='1'):
                count+=1
        if count!=q:
            score-=penalty

        #apply penalty for multiple similar stock
        for j in substrings:
            count=0
            for k in j:
                if(k=='1'):
                    count+=1
            if count!=1:
                score-=penalty

        # #apply penalty for x_{ij}>y_j
        # for i in range(n):
        #     for j in range(n):
        #         if bitstring[n*i+j]>bitstring[n**2+j]:
        #             score-=penalty

        for j in range(n):
            if(bitstring[n*j+j]!=bitstring[n**2+j]):
                score-=penalty

        totalScore += score * counts[bitstring]  #Multiply score times the # of times it was observed
        totalSamples += counts[bitstring]        #Keep track of the number of measurements (samples)
    print("Cost: "+str(totalScore))
    return(totalScore/totalSamples)

#Run the circuit and return counts
def runCKT(params):
    #TODO: Add noise simulation
    simulator = Aer.get_backend('qasm_simulator')
    shots=nSamples
    betas=[]
    gammas=[]
    for i in range(len(params)):
        if(i<len(params)/2):
            betas.append(params[i])
        else:
            gammas.append(params[i])
    circuit=QAOA(n,q,p,betas,gammas)
    transpiled = transpile(circuit,backend=simulator)
    counts=simulator.run(transpiled,shots=shots).result().get_counts(circuit)  
    return counts

#Run a circuit and get expectation value
def ExpectationValue(params):
    #Run circuit and collect counts
    counts = runCKT(params)
    # print(counts)
    #Get the score of the counts
    score = computeExpectationValue(counts,q,n)
    return(-score)

def optimize(n,q,p,nIterations,nSamples,params,rhobeg):
    out = minimize(ExpectationValue,x0=params,method="COBYLA",options={'maxiter':nIterations,'rhobeg':rhobeg})
    print(f'Out: {out}')
    optimal_params = out['x'] 
    counts=runCKT(optimal_params)
    # plt.bar(list(counts.keys()),list(counts.values()),width=0.9)
    # plt.xlabel("bitstrings")
    # plt.ylabel("counts")
    # plt.title("Results")
    # plt.show()
    return max(zip(counts.values(), counts.keys()))[1]


#4 stocks
n = 4
#pick 2 representative stocks
q = 2
#4 copies
p = 4

correlation_matrix=np.loadtxt("correlationmatrix_size4.txt")
#A sufficient number of optimization iterations to solve problem
nIterations = 10000
#Typically need quite a few samples (measurements of quantum circuit) per iteration to 
nSamples = 10000
#initial trust
rhobeg=1000.0
#choose parameters near 0
params=[]
maxEdgeWeight=0
for i in range(len(correlation_matrix)):
    for j in range(len(correlation_matrix[0])):
        maxEdgeWeight=max(maxEdgeWeight,correlation_matrix[i][j])
for i in range(p*2):
    params.append(0.01*np.random.rand())
print("Best bitstring:",optimize(n=n,q=q,p=p,nIterations=nIterations,nSamples=nSamples,params=params,rhobeg=rhobeg))