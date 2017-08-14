#This file will contain methods for getting dendrograms of patients/samples

#Current proposed methods:
#   -Hierarchical clustering of omics data
#   -Clustering of individual PCSF runs dynamically 
#   -Pre-determined subtyping
#   



#Including some dummy code for an approximation of method 1 - clustering data
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

#random data for now
a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size=[100,])
b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size=[50,])
X = np.concatenate((a, b),)

Z = linkage(X, 'ward') #returns an array of length n - 1, Z[i] will tell us which clusters were merged in the i-th iteration
# Maybe in this part, we can automatically try different clustering algorithms and select the ones based on goodness of cluster metrics...
# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()

return Z

# computing s and d indexs for each clades(I'm adding this):



