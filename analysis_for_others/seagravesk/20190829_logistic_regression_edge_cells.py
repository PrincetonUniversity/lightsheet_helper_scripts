#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 13:50:39 2019

@author: wanglab
"""

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
import pandas as pd, os, matplotlib.pyplot as plt, numpy as np, seaborn as sns
from sklearn import metrics

"""
zscore formula: z = (x(sample mean) - mu(population aka grouth truth mean))/(std. dev. of population)
https://www.datacamp.com/community/tutorials/understanding-logistic-regression-python
"""
src = "/jukebox/wang/zahra/kelly_cell_detection_analysis"

rcells = pd.read_csv(os.path.join(src, "real_cell_stats.csv"))
ecells = pd.read_csv(os.path.join(src, "edge_cell_stats.csv"))

#find variables
rc = rcells.copy().dropna()
params = [xx for xx in rc.columns if xx != "cell_id"]

#find population statistics
mus = {}; sigmas = {}
for param in params:
    mus[param] = rc[param].mean()
    sigmas[param] = rc[param].std()

#calculate z scores for edge cells & make params into 2d array
ec = ecells.copy().dropna()
rc_params = []; ec_params = []

for param in params:
    rc[param] = (rcells[param]-mus[param])/sigmas[param]
    ec[param] = (ecells[param]-mus[param])/sigmas[param]
    rc_params.append(rc[param]); ec_params.append(ec[param])
    
rc_params = np.array(rc_params); ec_params = np.array(ec_params)
#%%
#edge cells
#boxplots of all params
plt.figure()
plt.boxplot(ec_params.T, vert = False, labels = params, sym = "", showcaps = False)
plt.axvline(x=0, color="gray")
plt.xlabel("Z-score")
plt.xlim([-20, 10])
plt.title("Z-scored features of edge cells")

#boxplots excluding chisq
plt.figure()
plt.boxplot(ec_params[6:].T, vert = False, labels = params[6:], sym = "", showcaps = False)
plt.axvline(x=0, color="gray")
plt.xlabel("Z-score")
plt.title("Z-scored features of edge cells w/o Chi stats")

#real cells
#boxplots of all params
plt.figure()
plt.boxplot(rc_params.T, vert = False, labels = params, sym = "", showcaps = False)
plt.axvline(x=0, color="gray")
plt.xlabel("Z-score")
plt.title("Z-scored features of real cells")

#boxplots excluding chisq
plt.figure()
plt.boxplot(rc_params[6:].T, vert = False, labels = params[6:], sym = "", showcaps = False)
plt.axvline(x=0, color="gray")
plt.xlabel("Z-score")
plt.title("Z-scored features of real cells w/o Chi stats")

#%%
#combine both cells into one dataset with a label colum
ecs = ecells.drop(columns=["cell_id"]).dropna()
rcs = rcells.drop(columns=["cell_id"]).dropna()
ecs["label"] = np.zeros(len(ecs))
rcs["label"] = np.ones(len(rcs))
allcells = pd.concat([ecs, rcs])

X = allcells[params] # Features
y = allcells.label # Target variable

X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.25,random_state=0)
# instantiate the model (using the default parameters)
logreg = LogisticRegression()
# fit the model with data
logreg.fit(X_train,y_train)

y_pred = logreg.predict(X_test)
lbls = y_test.values

# import the metrics class
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')

print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
print("Precision:",metrics.precision_score(y_test, y_pred))
print("Recall:",metrics.recall_score(y_test, y_pred))

plt.figure()
y_pred_proba = logreg.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)
plt.plot(fpr,tpr,label="auc="+str(auc))
plt.legend(loc=4)
plt.show()

#%%
prob = logreg.predict_proba(X_test)
true_cell_probs = [xx for i, xx in enumerate(prob[:,1]) if y_test.values[i] == 1]
sns.swarmplot(true_cell_probs, orient = "v")
plt.ylabel("Classifier probability")
plt.xlabel("Real cells (n = 67)")
plt.axhline(y=0.5, color = "gray", linestyle = "--")