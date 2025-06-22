# baseline.py

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from xgboost import XGBClassifier

# 🔹 Step 1: 读取表达矩阵和 phenotype 数据
expr = pd.read_csv("GSE78220_expr_clean.csv", index_col=0).T
pheno = pd.read_csv("GSE78220_pheno.csv", index_col=0)

# 🔹 Step 2: 对齐样本
X = expr.copy()
y = pheno.loc[X.index, "response_bin"]

# 🔍 可选：检查标签是否为 0 和 1
print("标签分布：")
print(y.value_counts())

# 🔹 Step 3: 划分训练测试集
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42, stratify=y
)

# 🔹 Step 4: 建模（XGBoost）
clf = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
clf.fit(X_train, y_train)

# 🔹 Step 5: 预测与评估
y_pred = clf.predict(X_test)
y_proba = clf.predict_proba(X_test)[:, 1]

print("\n🎯 分类结果：")
print(classification_report(y_test, y_pred, digits=4))

auc = roc_auc_score(y_test, y_proba)
print(f"\n🔍 AUC score: {auc:.4f}")
