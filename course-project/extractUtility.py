from sklearn.semi_supervised import LabelPropagation
import numpy as np
from random import sample

X = np.genfromtxt('/Users/dgy/Desktop/385project/grouped_GPS.csv',delimiter=',')
X = X[1:(X.shape[0]-1) ,1:3]

rows = X.shape[0]
seedsX = X[sample(range(0, rows), 500),]
seedsY = np.repeat([0,0,0,0,1],100)

lp = LabelPropagation(gamma = 10000)
lp.fit(seedsX,seedsY)

Y = lp.predict_proba(X)
p = Y[0:rows, 1]
(np.where(p > 0.5))[0].shape
np.count_nonzero(np.isnan(p))
np.savetxt("utility.txt", Y, delimiter=" ", fmt="%s")
