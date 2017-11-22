import sys

class ReactionPathwayPlugin:
   def input(self, filename):
      self.myfile = filename

   def run(self):
      clusterfile = open(self.myfile+".clusters.noa", 'r')
      pathwayfile = open(self.myfile+".pathways.noa", 'r')
      pathwaynamesfile = open(self.myfile+".pathways.txt", 'r')
      networkfile = open(self.myfile+".csv", 'r')
      
      self.clusters = dict()
      self.pathways = dict()
      self.pathwaynames = dict()
      firstline = clusterfile.readline()
      firstline = pathwayfile.readline()
      firstline = pathwaynamesfile.readline()
      for line in clusterfile:
         chem, cluster = line.split('\t')
         self.clusters[chem] = int(cluster)
      for line in pathwayfile:
         chem, pathway = line.split('\t')
         self.pathways[chem] = int(pathway)
      for line in pathwaynamesfile:
         id, pathway = line.split('\t')
         self.pathwaynames[int(id)] = pathway.strip()

      lastpathway = id
      self.pospath = []
      self.zeropath = []
      self.negpath = []
      self.corrsum = []
      for i in range(int(lastpathway)+1):
         self.pospath.append(0)  
         self.zeropath.append(0)  
         self.negpath.append(0)  
         self.corrsum.append(0)  

      firstline = networkfile.readline()
      self.bacteria = firstline.split(',')
      self.bacteria.remove('\"\"')
      self.n = len(self.bacteria)
      for i in range(self.n):
         self.bacteria[i] = self.bacteria[i].strip()
         self.bacteria[i] = self.bacteria[i][1:len(self.bacteria[i])-1]

      self.ADJ = []
      i = 0
      for line in networkfile:
         contents = line.split(',')
         self.ADJ.append([])
         for j in range(self.n):
            if (contents[j+1] == '-Inf'):
               print "Infinity between ", self.bacteria[i], " and ", self.bacteria[j]
            value = float(contents[j+1])
            if (i != j and value != 0):
               self.ADJ[i].append(value)
            else:
               self.ADJ[i].append(0)
         i += 1






   def output(self, filename):
      # Given two nodes on the same pathway, how often
      # are they in the same cluster
      # How often different clusters
      samecluster = 0
      diffcluster = 0
      positive_edges = 0
      no_edge = 0
      negative_edges = 0
      correlation_sum = 0
      print self.pathways
      print self.clusters
      print self.bacteria
      for chem in self.pathways:
         for chem2 in self.pathways:
            if (chem2 != chem and self.pathways[chem] == self.pathways[chem2]):
               if (self.clusters.has_key(chem) and self.clusters.has_key(chem2) and self.clusters[chem] == self.clusters[chem2]):
                  samecluster += 1
               else:
                  diffcluster += 1
               idx1 = self.bacteria.index(chem)
               idx2 = self.bacteria.index(chem2)
               if (self.ADJ[idx1][idx2] > 0):
                  positive_edges += 1
                  self.pospath[self.pathways[chem]] += 1
               elif (self.ADJ[idx1][idx2] == 0):
                  no_edge += 1
                  self.zeropath[self.pathways[chem]] += 1
               else:
                  negative_edges += 1
                  self.negpath[self.pathways[chem]] += 1
               correlation_sum += self.ADJ[idx1][idx2]
               self.corrsum[self.pathways[chem]] += self.ADJ[idx1][idx2]

      # We counted each pair twice, so divide by two
      samecluster /= 2
      diffcluster /= 2
      positive_edges /= 2
      no_edge /= 2
      negative_edges /= 2
      correlation_sum /= 2

      pos = 0
      neg = 0
      zero = 0
      corr_diff = 0
      corr_diff_count = 0
      for i in range(self.n):
         for j in range(i+1, self.n):
            if (self.ADJ[i][j] > 0):
               pos += 1
            elif (self.ADJ[i][j] == 0):
               zero += 1
            else:
               neg += 1
            if (self.pathways[self.bacteria[i]] != self.pathways[self.bacteria[j]]):
               corr_diff += self.ADJ[i][j]
               corr_diff_count += 1
      

      print ("Metabolites on same pathway: %d in same cluster, %d in different clusters, %f percent" % (samecluster, diffcluster, samecluster / float(samecluster+diffcluster)))
      print ("Positively correlated: %d, Zero correlated: %d, Negatively correlated: %d, %f percent positive" % (positive_edges, no_edge, negative_edges, positive_edges / float(positive_edges+no_edge+negative_edges)))
      print ("Average correlation value (same): %f" % (correlation_sum / (positive_edges+no_edge+negative_edges)))
      print ("Average correlation value (diff): %f" % (corr_diff / corr_diff_count))
      print ("Percent of positive edges consumed by metabolites on same pathway: %f" % (positive_edges / float(pos)))
      print ("Percent of zero edges consumed by metabolites on same pathway: %f" % (no_edge / float(zero)))
      print ("Percent of negative edges consumed by metabolites on same pathway: %f" % (negative_edges / float(neg)))

      filestuff = open(filename, 'w')
      filestuff.write('Pathway\tPositive\tZero\tNegative\tPct Positive\tAverage Correlation\n')
      for key in self.pathwaynames:
          total = float((self.pospath[key]+self.zeropath[key]+self.negpath[key])/2)
          if (total > 0):
            average = self.pospath[key] / 2 / total
            filestuff.write(self.pathwaynames[key]+'\t\t\t\t'+str(self.pospath[key]/2)+'\t'+str(self.zeropath[key]/2)+'\t'+str(self.negpath[key]/2)+'\t'+str(average*100)+'\t'+str(self.corrsum[key]/2/total)+'\n')



