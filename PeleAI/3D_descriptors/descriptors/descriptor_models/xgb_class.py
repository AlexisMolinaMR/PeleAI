from xgboost import XGBClassifier
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import accuracy_score

train = pd.read_csv("train.csv")
test = pd.read_csv("test.csv")

ic50_train = train['IC50']
ic50_test = test['IC50']

train.drop(['Ligand','IC50'], axis=1, inplace=True)
test.drop(['Ligand','IC50'], axis=1, inplace=True)

train.fillna(value=0.0, inplace=True)
test.fillna(value=0.0, inplace=True)

ic50_train_enc = [1 if ic50 > -7 else 0 for ic50 in ic50_train.values]
ic50_test_enc = [1 if ic50 > -7 else 0 for ic50 in ic50_test.values]


scaler = MinMaxScaler()

train_mmsc = scaler.fit_transform(train)
test_mmsc = scaler.fit_transform(test)

xgb_clf = XGBClassifier()

xgb_clf.fit(train_mmsc, ic50_train_enc)

print(xgb_clf)

y_pred = model.predict(X_test)
predictions = [round(value) for value in y_pred]

accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
