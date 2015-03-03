''' Searching for optimal 1D partition for 3D corellation modularity.
'''

# Import section
import igraph
import random
import numpy as np
from copy import deepcopy
from pickle import load

def randomG(n):
    out = igraph.Graph(n)
    
    temp = random.getstate()
    random.seed(1234)
    for node in range(n):
        for connection in range(random.randint(1, n/10)):
            delta = int(random.gauss(0.0, n/10))
            e = node + delta
            if e >= 0 and e <= n-1 and delta:
                if not out[node,e] and not out[e, node]:
                    out.add_edge(node, e)
                    
    # Plant module
    for i in range(50,60):
        for j in range(i+1,60):
            out.add_edge(i,j)
            
    random.setstate(temp)
                
    return out

# GA class
class GAsolution:
    def __init__(self, mem, score):
        self.mem = mem
        self.score = score
        
class GAStat:
    def __init__(self):
        self.minscore = 0.0
        self.maxscore = 0.0
        self.medianscore = 0.0
        self.meanscore = 0.0
        self.timestamp = 0
        
class DomainGA:
    def __init__(self, G):
        # igraph instance
        self.G = G
        
        # Infer graph size
        self.length = len(self.G.vs)
        
        # GA stuff [memvect] = fitness
        self.population = []
        
        # GA parameters
        self.popsize = 100
        
        self.rate_mutation = 0.3
        self.rate_insertion = 0.3
        self.rate_deletion = 0.3
        self.rate_xover = 0.1
        
        # Stats
        self.history = []
        self.timestamp = 0
        
        # Best Solution
        self.best = GAsolution([], self.Score([]))
        self.epochbest = []
        
        
    def Score(self, mem):
        return self.G.modularity(self.MembershipVector(mem), weights='wts')
    
    def SeedPopulation(self, mass_extinction=False):
        
        # clear population
        if mass_extinction:
            self.population = self.population[:self.popsize/10]
        else:
            self.population = []
        
        # Iterate
        while len(self.population) < self.popsize:
            mem = self.GetRandomMem()
            self.population.append(GAsolution(mem, self.Score(mem)))
            
        self.Record()
            
    def Select(self, n=2):
        # Tournament
        out = random.randint(0,len(self.population)-1)
        
        while n:
            t = random.randint(0,len(self.population)-1)
            if self.population[t].score > self.population[out].score:
                out = t
            n -= 1
            
        return out
    
    def Roulette(self):
        # Get all scores
        vals = []
        for i in self.population:
            vals.append(i.score)
        
        # random value
        rnd = random.random() * sum(vals)
        
        for i in range(len(vals)):
            if rnd > sum(vals[:i]):
                return i
        return i
    
    def Mutate(self, n):
        # Individual
        parent = self.population[n]
        indiv = deepcopy(parent)
        
        # Locus
        if not indiv.mem:
            return indiv
        
        # multiple mutations
        for mut in range(random.randint(1,len(indiv.mem))):
            locus = random.randint(0,len(indiv.mem)-1)
            
            # Delta
            delta = random.randint(-5,5)
            newbound = max(0,min(indiv.mem[locus] + delta,self.length-1))
            
            # Negative delta
            if locus == 0 or locus == len(indiv.mem)-1:
                pass
            elif newbound > indiv.mem[locus-1] or newbound < indiv.mem[locus+1]:
                indiv.mem[locus] = newbound
            
        return indiv
            
            
    
    def Delete(self, n):
        # Individual
        parent = self.population[n]
        indiv = deepcopy(parent)
        
        # Null case
        if len(indiv.mem) == 0:
            return indiv
        
        # Locus
        locus = random.randint(0,len(indiv.mem)-1)    
        
        # Delete
        indiv.mem.remove(indiv.mem[locus])
        
        return indiv
        
    
    def Insert(self, n):
        # Individual
        parent = self.population[n]
        indiv = deepcopy(parent)
        
        # Locus
        locus = random.randint(0,self.length-1)    
        
        # Delete
        if not locus in indiv.mem:
            indiv.mem.append(locus)
            indiv.mem.sort()
        
        return indiv        
    
    def CrossOver(self, n, x):
        p1 = self.population[n].mem
        p2 = self.population[x].mem
        
        # Locus
        locus = random.randint(0,self.length-1)  
        
        newmem = []
        for i in p1:
            if i <= locus:
                newmem.append(i)
        for i in p2:
            if i > locus:
                newmem.append(i)
        
        # new individual
        return GAsolution(newmem, self.Score(newmem))
            
    def AddToPool(self, indiv):
        # Grow pool
        if len(self.population) < self.popsize:
            self.population.append(indiv)
        else:
            # Randomly replace a genome
            position = random.randint(0,self.popsize-1)
            
            self.population[position] = indiv
        
            
    def IterateOnce(self):
        # Select a genome
        genome_i = self.Roulette()
        
        
        # Select an operator
        x = random.random()
        
        if x <= self.rate_mutation:
            indiv = self.Mutate(genome_i)
        elif x <= self.rate_mutation + self.rate_insertion:
            indiv = self.Insert(genome_i)
        elif x <= self.rate_mutation + self.rate_insertion+self.rate_deletion:
            indiv = self.Delete(genome_i)        
        else:
            indiv = self.CrossOver(genome_i,self.Roulette())
            
        # Make sure that there are no 
        indiv.mem = list(set(indiv.mem))
        indiv.mem.sort()
        indiv.score = self.Score(indiv.mem)
        
        self.AddToPool(indiv)
        
        if indiv.score > self.best.score:
            self.best = indiv
            self.Record()
        
        self.timestamp += 1
        
    def Iterate(self, N):
        # Do the loop
        for i in range(N):
            self.IterateOnce()
            
        self.epochbest.append(self.best)
        self.best = GAsolution([], self.Score([]))

        
    def Record(self):
        # Stats
        stat = GAStat()
        
        # vals
        vals = []
        for i in self.population:
            vals.append(i.score)
            
        a = stat.maxscore = max(vals)
        b = stat.minscore = min(vals)
        c = stat.meanscore = sum(vals)/len(vals)
        vals.sort()
        d = stat.medianscore = vals[int(len(vals)/2)]
        e = stat.timestamp = self.timestamp
        
        self.history.append(stat)
        
        print e, b, c, d, a, self.best.mem, self.best.score
            
            
            
    def GetRandomMem(self):
        mem = []
        for i in range(random.randint(0,10)):
            x = random.randint(0,self.length)
            if not x in mem:
                mem.append(x)
                
        mem.sort()
        
        return mem
    
    def MembershipVector(self, v):
        mem = []
        symbol = 0
        for i in range(self.length):
            if i in v:
                symbol += 1
            mem.append(symbol)
        return mem
    
    def AllTimeBest(self):
        # FInd the best from all epoch
        out = self.epochbest[0]
        
        for i in self.epochbest:
            if i.score > out.score:
                out = i
        
        return self.MembershipVector(out.mem)
    
    def MermaidMagicModularity(self, runlength=10000, epoch=100):
        self.SeedPopulation()
        for e in range(epoch):
            self.Iterate(runlength)
            self.SeedPopulation()
            
        # Last victory lap
        self.population = deepcopy(self.epochbest)
        self.Iterate(100000)
        
        return self.AllTimeBest()        

if __name__ == "__main__":
    # create a random graph for now
    #G = randomG(100)
    
    # Load a graph
    G = load(open('amylase.pickle','rb'))
    
    # Create GA
    ga = DomainGA(G)
    print ga.MermaidMagicModularity()
    
    
    