from sklearn import tree, model_selection, svm
from sklearn.feature_selection import SelectFromModel
from sklearn import preprocessing
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc
import pandas as pd, pydotplus, goatools, numpy as np, matplotlib.pyplot as plt


df = pd.DataFrame.from_csv("./data/yeast_data.csv")
df.to_pickle("./data/yeast_data.pkl")

df = pd.read_pickle("./data/yeast_data.pkl")
p = goatools.obo_parser.GODag('../go.obo')


print df.shape
print df["Outcome"].value_counts()

feature_names = [(x + ": " + p[x].name + "\n") for x in df.drop("Outcome",1).columns.values]
X = df.drop("Outcome",1).values
y = df["Outcome"].values

skf = model_selection.StratifiedKFold(n_splits=10)
skf.get_n_splits(X, y)

y_pred_dtc = []
y_true = []
y_score_svm = []
y_pred_svm = []

for fold_num, (train, test) in enumerate(skf.split(X, y)):
    clf = tree.DecisionTreeClassifier(random_state=42,
                                          # To calculate a leaf, the feature has to be evident in 5% of the dataset.
                                          min_samples_leaf=0.05,
                                          # For classification purposes, it's usually a good idea to limit the max no
                                          # of features to the sqrt of how many features are being considered
                                          #max_features=int(sqrt(len(feature_names))),
                                          class_weight="balanced")
    clf.fit(X[train],y[train])


    important_features = []
    for index, feature_score in enumerate(clf.feature_importances_):
        important_features.append([feature_names[index], feature_score])

    pd.DataFrame(sorted(important_features, key=lambda x: x[1]),
                 columns=["GO Term", "Gini Score"]).to_excel\
        ("./results/DTC/feature_importances/fold_"+str(fold_num+1)+".xlsx", index=False)

    dot_file = tree.export_graphviz(clf, feature_names=feature_names,
                                    filled=True, rounded=True, class_names=y, out_file=None,
                                    proportion=True, impurity=False, max_depth=10)
    graph = pydotplus.graph_from_dot_data(dot_file)
    graph.write_pdf("./results/DTC/trees/fold_"+str(fold_num+1)+".pdf")
    # Evaluation
    y_pred_dtc.extend(clf.predict(X[test]))
    y_true.extend(y[test])

    fr_X = SelectFromModel(clf, prefit=True).transform(X)
    svm_classifier = svm.SVC(class_weight="balanced", kernel="poly")
    svm_classifier.fit(fr_X[train], y[train])
    y_score_svm.extend(svm_classifier.decision_function(fr_X[test]))
    y_pred_svm.extend(svm_classifier.predict(fr_X[test]))

lb = preprocessing.LabelBinarizer()
false_positive_ratio, true_positive_ratio, _ = roc_curve(lb.fit_transform(y_true), y_score_svm)
area_under_curve = auc(false_positive_ratio, true_positive_ratio)

plt.figure()
plt.plot(false_positive_ratio, true_positive_ratio, color="b", lw=4, label="ROC curve (area = % 0.2f)" % area_under_curve)
plt.xlim([0.0, 1.0])
plt.ylim([0., 1.05])
plt.xlabel("Specificity")
plt.ylabel("Sensitivity")
plt.fill_between(false_positive_ratio, true_positive_ratio, color="k", alpha=0.1)
plt.legend(loc="lower right")
plt.savefig("./results/SVM/ROC.pdf", type="pdf")
###################

with open("./results/DTC/classification_report.txt", "w") as cr_out:
    cr_out.write(classification_report(y_true, y_pred_dtc))
cr_out.close()

with open("./results/SVM/classification_report.txt", "w") as cr_out:
    cr_out.write(classification_report(y_true, y_pred_svm))
cr_out.close()

with open("./results/DTC/confusion_matrix.txt", "w") as cm_out:
    cm_out.write(np.array2string(confusion_matrix(y_true, y_pred_dtc), separator="\t"))
cm_out.close()

with open("./results/SVM/confusion_matrix.txt", "w") as cm_out:
    cm_out.write(np.array2string(confusion_matrix(y_true, y_pred_svm), separator="\t"))
cm_out.close()