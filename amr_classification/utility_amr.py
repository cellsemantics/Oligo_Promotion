
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO

from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif, SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    matthews_corrcoef, accuracy_score, precision_score,
    recall_score, f1_score, cohen_kappa_score,
    roc_auc_score, confusion_matrix,
    auc, roc_curve, classification_report, average_precision_score, precision_recall_curve
)
from sklearn.utils.class_weight import compute_class_weight

from imblearn.pipeline import Pipeline as imbpipeline
from imblearn.over_sampling import SMOTE, SVMSMOTE, BorderlineSMOTE

import seaborn as sns
import pickle
import kaos
import utility

from sklearn.decomposition import PCA
from sklearn.preprocessing import QuantileTransformer, StandardScaler, RobustScaler

# Other imports from your code




def return_model(drug, type_classifier = "rf", type_embedding = "label_embedding"):

    from imblearn.pipeline import Pipeline as imbpipeline
    from imblearn.over_sampling import SMOTE, SVMSMOTE, BorderlineSMOTE
    from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif, SelectFromModel
    from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, GradientBoostingClassifier

    if type_embedding == "oligo_embedding":
        if type_classifier=="rf":
            if (drug.lower() == "cip"):

                pipeline = imbpipeline([
                    ('variance_threshold', VarianceThreshold(threshold=0.1)),
                    ('feature_selection', SelectKBest(score_func=f_classif, k=5000)),
                    ('clf', RandomForestClassifier(random_state=0, n_estimators=500, n_jobs=18, min_samples_leaf=2, bootstrap= False))])

            if (drug.lower() == "ctx"):

                pipeline = imbpipeline([
            ('variance_threshold', VarianceThreshold(threshold=0.1)),
            ('feature_selection', SelectKBest(score_func=f_classif, k=1500)),
            ('clf', RandomForestClassifier(random_state=0, n_estimators=300, n_jobs=18, min_samples_leaf=2, bootstrap= False
                                        , class_weight="balanced"))  ]) 

            if (drug.lower() == "gen"):

                pipeline = imbpipeline([
        ('variance_threshold', VarianceThreshold(threshold=0.1)),
        ('feature_selection', SelectKBest(score_func=f_classif, k=1500)),
        ('clf', RandomForestClassifier(random_state=0, n_estimators=500, n_jobs=18, class_weight="balanced", min_samples_leaf=3, bootstrap= True
                                       ))  ])
                
            if (drug.lower() == "ctz"):

                pipeline =  imbpipeline([
            ('variance_threshold', VarianceThreshold(threshold=0.1)),
            ('feature_selection', SelectKBest(score_func=f_classif, k=2500)),
            ('clf', RandomForestClassifier(random_state=0, n_estimators=500, n_jobs=18, min_samples_leaf=2, bootstrap= False))])
                

    if type_embedding == "label_embedding":
        if type_classifier=="rf":
            pipeline =    imbpipeline([
                        ('smote', SMOTE(random_state=0)),  
                        ('clf', RandomForestClassifier(random_state=0, n_estimators=200, n_jobs=16)) 
                    ])
            
        if type_classifier=="svm":
            pipeline =    imbpipeline([
                        ('smote', SMOTE(random_state=0)),  
                        ('clf', SVC(kernel='linear', probability=True)) ])
            

        if type_classifier=="lr":
            pipeline =    imbpipeline([
                        ('smote', SMOTE(random_state=0)),  
                        ('clf', LogisticRegression(solver = 'lbfgs',max_iter=1000))])



        

    return pipeline


def return_result_for_a_drug(drug, X_log, data_gessin, type_embedding="label_embedding", type_classifier = "rf"):

    y_gen = data_gessin[drug]

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

    lst_model = list()
    lst_precision = list()
    lst_recall = list()
    lst_mcc = list()
    lst_auc = list()
    lst_f1 = list()

    y_cv_true = list()
    y_cv_predicted = list()
    y_pred = list()

    print("\n The result of drug is: ", drug)


    for i, (train_index, val_index) in enumerate(skf.split(X_log, y_gen)):

        X_train_fold, y_train_fold = X_log.iloc[train_index, :], y_gen.iloc[train_index]
        X_val_fold, y_val_fold = X_log.iloc[val_index, :], y_gen.iloc[val_index]

        y_cv_true.extend(y_val_fold)

        if type_embedding=="oligo_embedding":
        
            pipeline = return_model(drug, type_classifier = type_classifier, type_embedding = type_embedding)

        if type_embedding=="label_embedding":

            pipeline = return_model(drug, type_classifier = type_classifier, type_embedding = type_embedding)
        pipeline.fit(np.array(X_train_fold), np.array(y_train_fold))

        lst_model.append(pipeline)

        y_pred.extend(pipeline.predict_proba(np.array(X_val_fold))[:, 1])

        print("Result of val set for fold: ", str(i))
        print("f1_score: ", str(f1_score(y_val_fold, pipeline.predict(np.array(X_val_fold)))), " | mcc_score: ", str(matthews_corrcoef(y_val_fold, pipeline.predict(np.array(X_val_fold)))),
            " |  precision_score: ", str(precision_score(y_val_fold, pipeline.predict(np.array(X_val_fold)))), " | recall_score: ", str(recall_score(y_val_fold, pipeline.predict(np.array(X_val_fold)))),
            )
        
        lst_precision.append(precision_score(np.array(y_val_fold), pipeline.predict(np.array(X_val_fold))))
        lst_recall.append(recall_score(np.array(y_val_fold), pipeline.predict(np.array(X_val_fold))))
        lst_mcc.append(matthews_corrcoef(np.array(y_val_fold), pipeline.predict(np.array(X_val_fold))))
        lst_auc.append(roc_auc_score(y_val_fold, pipeline.predict_proba(np.array(X_val_fold))[:, 1]))
        lst_f1.append(f1_score(y_val_fold, pipeline.predict(np.array(X_val_fold))))

    print("\n The result of mean of cross validation set across the five fold is: ")
    print("The mean and std of precision is: ", str(np.mean(lst_precision)) , str(np.std(lst_precision)))
    print("The mean and std of recall is: ", str(np.mean(lst_recall)) , str(np.std(lst_recall)))
    print("The mean and std of f1 is: ", str(np.mean(lst_f1)) , str(np.std(lst_f1)))
    print("The mean and std of mcc is: ", str(np.mean(lst_mcc)) , str(np.std(lst_mcc)))
    print("The mean and std of roc_auc is: ", str(np.mean(lst_auc)) , str(np.std(lst_auc)))

    return pipeline, y_pred, y_cv_true

    
    
def custom_figure_axis(ax, fontsize=10, show_ticks = True, rotation=90):
    
        
    """
    Customize the appearance of matplotlib axis for a figure.
    
    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis to be customized.
        fontsize (int, optional): Font size for axis labels and ticks. Default is 10.
        show_ticks (bool, optional): Whether to display ticks and labels. Default is True.
    
    Returns:
        matplotlib.axes._subplots.AxesSubplot: The customized axis.
    """
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.2)
    ax.spines['left'].set_linewidth(0.2)
    ax.tick_params(axis='x', labelsize=fontsize, rotation=rotation)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.tick_params(axis='both', which='both', width=0.5)
    ax.xaxis.label.set_fontsize(fontsize)
    ax.yaxis.label.set_fontsize(fontsize)
    
    if show_ticks==False:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
    return ax
