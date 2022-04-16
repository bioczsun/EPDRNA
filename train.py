from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
import numpy as np
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
import joblib
from imblearn.combine import SMOTEENN

import warnings
warnings.filterwarnings('ignore')
model_path = r'.\Model\dnamodel.pkl' #model output file
data_path = r'.\Featuredir\train_dna.npy'

data = np.load(data_path)
label = np.loadtxt(r'.\data\traindnalabel.txt')##label txt

print('----------------------------------------Data load----------------------------------------')


print('----------------------------------------Train model----------------------------------------')


cc = SMOTEENN(random_state=12,n_jobs=5)

x, y = cc.fit_resample(data, label)

log_clf = LogisticRegression(n_jobs=5, max_iter=100)
rnd_clf = RandomForestClassifier(n_estimators=100, n_jobs=5)
knn_clf = KNeighborsClassifier(n_jobs=5)
xgbc_clf = XGBClassifier(n_estimators=100, n_jobs=5)

voting_clf = VotingClassifier(
        estimators=[('lc', log_clf), ('rf', rnd_clf,), ('knn', knn_clf), ('xgbc_clf', xgbc_clf)],
        voting='soft', n_jobs=10)

voting_clf.fit(x, y)
joblib.dump(voting_clf,model_path)
print('----------------------------------------Train completed----------------------------------------')


