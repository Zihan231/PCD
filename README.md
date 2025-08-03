!wget -O variant_summary.txt.gz
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary
.txt.gz
--2025-08-02 18:22:50--
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary
.txt.gz
Resolving ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)... 130.14.250.7,
130.14.250.10, 130.14.250.13, ...
Connecting to ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)|
130.14.250.7|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 382607393 (365M) [application/x-gzip]
Saving to: ‘variant_summary.txt.gz’
variant_summary.txt 100%[===================>] 364.88M 59.9MB/s in
6.5s
2025-08-02 18:22:56 (56.5 MB/s) - ‘variant_summary.txt.gz’ saved
[382607393/382607393]
import pandas as pd
# Load the compressed tab-delimited file
df = pd.read_csv('variant_summary.txt.gz', sep='\t',
compression='gzip', low_memory=False)
# Filter rows where 'PhenotypeList' contains "Primary ciliary
dyskinesia"
pcd_df = df[df['PhenotypeList'].str.contains('Primary ciliary
dyskinesia', na=False, case=False)]
print(f"Total variants related to Primary Ciliary Dyskinesia:
{len(pcd_df)}")
pcd_df.head()
Total variants related to Primary Ciliary Dyskinesia: 52747
 #AlleleID Type \
507 15302 Duplication
508 15302 Duplication
509 15303 single nucleotide variant
510 15303 single nucleotide variant
511 15304 single nucleotide variant
 Name GeneID GeneSymbol \
507 NM_178452.6(DNAAF1):c.1349dup (p.Pro451fs) 123872 DNAAF1
508 NM_178452.6(DNAAF1):c.1349dup (p.Pro451fs) 123872 DNAAF1
509 NM_178452.6(DNAAF1):c.811C>T (p.Arg271Ter) 123872 DNAAF1
510 NM_178452.6(DNAAF1):c.811C>T (p.Arg271Ter) 123872 DNAAF1
511 NM_178452.6(DNAAF1):c.792C>A (p.Tyr264Ter) 123872 DNAAF1
 HGNC_ID ClinicalSignificance ClinSigSimple LastEvaluated \
507 HGNC:30539 Pathogenic 1 Nov 27, 2023
508 HGNC:30539 Pathogenic 1 Nov 27, 2023
509 HGNC:30539 Pathogenic 1 Jun 27, 2022
510 HGNC:30539 Pathogenic 1 Jun 27, 2022
511 HGNC:30539 Pathogenic 1 Apr 04, 2024
 RS# (dbSNP) ... AlternateAlleleVCF SomaticClinicalImpact \
507 397515339 ... AC -
508 397515339 ... AC -
509 267607225 ... T -
510 267607225 ... T -
511 267607226 ... A -
 SomaticClinicalImpactLastEvaluated ReviewStatusClinicalImpact \
507 - -
508 - -
509 - -
510 - -
511 - -
 Oncogenicity OncogenicityLastEvaluated ReviewStatusOncogenicity \
507 - - -
508 - - -
509 - - -
510 - - -
511 - - -
 SCVsForAggregateGermlineClassification \
507 SCV004297037
508 SCV004297037
509 SCV000624531
510 SCV000624531
511 SCV001228496|SCV002579825|SCV004810007
 SCVsForAggregateSomaticClinicalImpact \
507 -
508 -
509 -
510 -
511 -
 SCVsForAggregateOncogenicityClassification
507 -
508 -
509 -
510 -
511 -
[5 rows x 43 columns]
# Select useful columns
cols = ['GeneSymbol', 'Type', 'ClinicalSignificance', 'Chromosome',
'Start', 'ReferenceAllele', 'AlternateAllele']
pcd_df = pcd_df[cols]
# Quick look at unique ClinicalSignificance values
print("Unique ClinicalSignificance values:")
print(pcd_df['ClinicalSignificance'].unique())
# Simplify ClinicalSignificance to binary labels:
# Pathogenic or Likely pathogenic -> 1
# Benign, Likely benign, Uncertain significance, others -> 0
def label_pathogenicity(cs):
 if isinstance(cs, str):
 cs = cs.lower()
 if 'pathogenic' in cs:
 return 1
 else:
 return 0
 return 0
pcd_df['Pathogenic'] =
pcd_df['ClinicalSignificance'].apply(label_pathogenicity)
# Check class balance
print("\nPathogenic vs Non-Pathogenic counts:")
print(pcd_df['Pathogenic'].value_counts())
pcd_df.head()
Unique ClinicalSignificance values:
['Pathogenic' 'Pathogenic/Likely pathogenic' 'Benign'
 'Uncertain significance' 'Benign/Likely benign' 'Likely pathogenic'
 'Conflicting classifications of pathogenicity' 'Likely benign'
 'not provided']
Pathogenic vs Non-Pathogenic counts:
Pathogenic
0 42221
1 10526
Name: count, dtype: int64
 GeneSymbol Type ClinicalSignificance
Chromosome \
507 DNAAF1 Duplication Pathogenic
16
508 DNAAF1 Duplication Pathogenic
16
509 DNAAF1 single nucleotide variant Pathogenic
16
510 DNAAF1 single nucleotide variant Pathogenic
16
511 DNAAF1 single nucleotide variant Pathogenic
16
 Start ReferenceAllele AlternateAllele Pathogenic
507 84203778 na na 1
508 84170172 na na 1
509 84193349 na na 1
510 84159744 na na 1
511 84193330 na na 1
# Count of Pathogenic (1) vs Non-Pathogenic (0) for Primary Ciliary
Dyskinesia
pcd_counts = df[df['PhenotypeList'].str.contains('Primary ciliary
dyskinesia', na=False)]
pathogenic_counts = pcd_counts['ClinSigSimple'].value_counts()
print(pathogenic_counts)
ClinSigSimple
0 45823
1 6922
Name: count, dtype: int64
from sklearn.preprocessing import LabelEncoder
# Columns to encode
cat_cols = ['GeneSymbol', 'Type', 'Chromosome', 'ReferenceAllele',
'AlternateAllele']
# Create LabelEncoders for each categorical column
encoders = {}
for col in cat_cols:
 le = LabelEncoder()
 pcd_df[col] = le.fit_transform(pcd_df[col].astype(str))
 encoders[col] = le
pcd_df.head()
 GeneSymbol Type ClinicalSignificance Chromosome Start \
507 22 1 Pathogenic 6 84203778
508 22 1 Pathogenic 6 84170172
509 22 7 Pathogenic 6 84193349
510 22 7 Pathogenic 6 84159744
511 22 7 Pathogenic 6 84193330
 ReferenceAllele AlternateAllele Pathogenic
507 1 2 1
508 1 2 1
509 1 2 1
510 1 2 1
511 1 2 1
from sklearn.model_selection import train_test_split
# Features and target
features = pcd_counts[['GeneSymbol', 'Type', 'Chromosome', 'Start',
'ReferenceAllele', 'AlternateAllele']]
target = pcd_counts['ClinSigSimple']
# First, split 70% train, 30% temp (validation + test)
X_train, X_temp, y_train, y_temp = train_test_split(features, target,
test_size=0.3, random_state=42, stratify=target)
# Then split the 30% temp into 15% validation and 15% test (half-half)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp,
test_size=0.5, random_state=42, stratify=y_temp)
print("Train samples:", len(X_train))
print("Validation samples:", len(X_val))
print("Test samples:", len(X_test))
Train samples: 36921
Validation samples: 7912
Test samples: 7912
from sklearn.preprocessing import LabelEncoder
# List of categorical columns
cat_cols = ['GeneSymbol', 'Type', 'Chromosome', 'ReferenceAllele',
'AlternateAllele']
# Initialize LabelEncoders for each categorical column
encoders = {col: LabelEncoder() for col in cat_cols}
# Fit and transform on training set, then transform validation and
test
for col in cat_cols:
 # Factorize train set
 X_train[col], uniques = pd.factorize(X_train[col].astype(str))
 # Map val and test to train categories, assign -1 for unseen
 X_val[col] = X_val[col].astype(str).map({k: i for i, k in
enumerate(uniques)}).fillna(-1).astype(int)
 X_test[col] = X_test[col].astype(str).map({k: i for i, k in
enumerate(uniques)}).fillna(-1).astype(int)
print(X_train.head())
print(X_val.head())
print(X_test.head())
 GeneSymbol Type Chromosome Start ReferenceAllele \
984467 0 0 0 132611349 0
5445033 1 0 1 21818277 0
6596178 1 0 1 21908547 0
3881936 1 0 1 21857948 0
883924 2 0 2 13913886 0
 AlternateAllele
984467 0
5445033 0
6596178 0
3881936 0
883924 0
 GeneSymbol Type Chromosome Start ReferenceAllele \
1824521 10 0 8 34506737 0
1076440 1 0 1 21840850 0
2435814 27 0 0 101252907 0
2930080 5 0 4 72285793 0
897423 5 0 4 72281279 0
 AlternateAllele
1824521 0
1076440 0
2435814 0
2930080 0
897423 0
 GeneSymbol Type Chromosome Start ReferenceAllele \
259268 4 0 4 78061541 0
3240060 1 0 1 21907642 0
1359769 1 0 1 21591188 0
6533044 2 2 2 13900434 0
980285 22 0 3 116627810 0
 AlternateAllele
259268 0
3240060 0
1359769 0
6533044 0
980285 0
from sklearn.preprocessing import LabelEncoder
cat_cols = ['GeneSymbol', 'Type', 'Chromosome', 'ReferenceAllele',
'AlternateAllele']
encoders = {}
for col in cat_cols:
 le = LabelEncoder()
 X_train[col] = le.fit_transform(X_train[col].astype(str))
 encoders[col] = le
 # Convert val/test: set unknowns to a new class index (len of
classes)
 known_classes = set(le.classes_)
 val_mapped = X_val[col].astype(str).map(lambda x:
le.transform([x])[0] if x in known_classes else len(known_classes))
 test_mapped = X_test[col].astype(str).map(lambda x:
le.transform([x])[0] if x in known_classes else len(known_classes))
 X_val[col] = val_mapped.astype(int)
 X_test[col] = test_mapped.astype(int)
import xgboost as xgb
from sklearn.metrics import classification_report, accuracy_score
# Define the classifier with moderate class balancing
xgb_model = xgb.XGBClassifier(
 objective='binary:logistic',
 eval_metric='aucpr',
 use_label_encoder=False,
 scale_pos_weight=2.5, # New weight to help class 1 recall
 learning_rate=0.1,
 n_estimators=100,
 max_depth=6,
 random_state=42
)
# Train the model
xgb_model.fit(
 X_train, y_train,
 eval_set=[(X_train, y_train), (X_val, y_val)],
 verbose=True
)
# Predict on test set
y_pred = xgb_model.predict(X_test)
# Evaluate
accuracy = accuracy_score(y_test, y_pred)
print(f"\nTest Accuracy: {accuracy}")
print(classification_report(y_test, y_pred))
[0] validation_0-aucpr:0.51498 validation_1-aucpr:0.49835
[1] validation_0-aucpr:0.53128 validation_1-aucpr:0.51293
[2] validation_0-aucpr:0.53648 validation_1-aucpr:0.52384
[3] validation_0-aucpr:0.53729 validation_1-aucpr:0.52264
[4] validation_0-aucpr:0.53996 validation_1-aucpr:0.52518
[5] validation_0-aucpr:0.54123 validation_1-aucpr:0.52468
[6] validation_0-aucpr:0.54167 validation_1-aucpr:0.52561
[7] validation_0-aucpr:0.54194 validation_1-aucpr:0.52670
[8] validation_0-aucpr:0.54313 validation_1-aucpr:0.52640
[9] validation_0-aucpr:0.54551 validation_1-aucpr:0.52898
[10] validation_0-aucpr:0.54636 validation_1-aucpr:0.52949
[11] validation_0-aucpr:0.54663 validation_1-aucpr:0.52907
[12] validation_0-aucpr:0.54868 validation_1-aucpr:0.53155
[13] validation_0-aucpr:0.54917 validation_1-aucpr:0.53242
[14] validation_0-aucpr:0.55086 validation_1-aucpr:0.53491
[15] validation_0-aucpr:0.55227 validation_1-aucpr:0.53522
[16] validation_0-aucpr:0.55280 validation_1-aucpr:0.53536
[17] validation_0-aucpr:0.55394 validation_1-aucpr:0.53572
[18] validation_0-aucpr:0.55489 validation_1-aucpr:0.53712
[19] validation_0-aucpr:0.55557 validation_1-aucpr:0.53877
[20] validation_0-aucpr:0.55623 validation_1-aucpr:0.53929
[21] validation_0-aucpr:0.55655 validation_1-aucpr:0.53961
[22] validation_0-aucpr:0.55810 validation_1-aucpr:0.54043
[23] validation_0-aucpr:0.55923 validation_1-aucpr:0.54167
[24] validation_0-aucpr:0.56014 validation_1-aucpr:0.54174
[25] validation_0-aucpr:0.56177 validation_1-aucpr:0.54327
[26] validation_0-aucpr:0.56240 validation_1-aucpr:0.54311
[27] validation_0-aucpr:0.56287 validation_1-aucpr:0.54301
[28] validation_0-aucpr:0.56325 validation_1-aucpr:0.54353
[29] validation_0-aucpr:0.56352 validation_1-aucpr:0.54365
[30] validation_0-aucpr:0.56439 validation_1-aucpr:0.54419
/usr/local/lib/python3.11/dist-packages/xgboost/training.py:183:
UserWarning: [18:25:58] WARNING: /workspace/src/learner.cc:738:
Parameters: { "use_label_encoder" } are not used.
 bst.update(dtrain, iteration=i, fobj=obj)
[31] validation_0-aucpr:0.56476 validation_1-aucpr:0.54370
[32] validation_0-aucpr:0.56532 validation_1-aucpr:0.54525
[33] validation_0-aucpr:0.56584 validation_1-aucpr:0.54492
[34] validation_0-aucpr:0.56618 validation_1-aucpr:0.54401
[35] validation_0-aucpr:0.56679 validation_1-aucpr:0.54406
[36] validation_0-aucpr:0.56702 validation_1-aucpr:0.54458
[37] validation_0-aucpr:0.56765 validation_1-aucpr:0.54491
[38] validation_0-aucpr:0.56812 validation_1-aucpr:0.54472
[39] validation_0-aucpr:0.56872 validation_1-aucpr:0.54524
[40] validation_0-aucpr:0.56892 validation_1-aucpr:0.54574
[41] validation_0-aucpr:0.56929 validation_1-aucpr:0.54592
[42] validation_0-aucpr:0.56943 validation_1-aucpr:0.54606
[43] validation_0-aucpr:0.56980 validation_1-aucpr:0.54597
[44] validation_0-aucpr:0.57013 validation_1-aucpr:0.54593
[45] validation_0-aucpr:0.57025 validation_1-aucpr:0.54634
[46] validation_0-aucpr:0.57043 validation_1-aucpr:0.54618
[47] validation_0-aucpr:0.57068 validation_1-aucpr:0.54668
[48] validation_0-aucpr:0.57112 validation_1-aucpr:0.54664
[49] validation_0-aucpr:0.57151 validation_1-aucpr:0.54658
[50] validation_0-aucpr:0.57187 validation_1-aucpr:0.54637
[51] validation_0-aucpr:0.57209 validation_1-aucpr:0.54656
[52] validation_0-aucpr:0.57244 validation_1-aucpr:0.54657
[53] validation_0-aucpr:0.57303 validation_1-aucpr:0.54676
[54] validation_0-aucpr:0.57342 validation_1-aucpr:0.54690
[55] validation_0-aucpr:0.57365 validation_1-aucpr:0.54686
[56] validation_0-aucpr:0.57413 validation_1-aucpr:0.54778
[57] validation_0-aucpr:0.57440 validation_1-aucpr:0.54724
[58] validation_0-aucpr:0.57475 validation_1-aucpr:0.54728
[59] validation_0-aucpr:0.57503 validation_1-aucpr:0.54728
[60] validation_0-aucpr:0.57538 validation_1-aucpr:0.54772
[61] validation_0-aucpr:0.57583 validation_1-aucpr:0.54776
[62] validation_0-aucpr:0.57613 validation_1-aucpr:0.54768
[63] validation_0-aucpr:0.57623 validation_1-aucpr:0.54763
[64] validation_0-aucpr:0.57667 validation_1-aucpr:0.54787
[65] validation_0-aucpr:0.57692 validation_1-aucpr:0.54753
[66] validation_0-aucpr:0.57715 validation_1-aucpr:0.54807
[67] validation_0-aucpr:0.57740 validation_1-aucpr:0.54795
[68] validation_0-aucpr:0.57762 validation_1-aucpr:0.54775
[69] validation_0-aucpr:0.57771 validation_1-aucpr:0.54766
[70] validation_0-aucpr:0.57802 validation_1-aucpr:0.54799
[71] validation_0-aucpr:0.57832 validation_1-aucpr:0.54781
[72] validation_0-aucpr:0.57848 validation_1-aucpr:0.54775
[73] validation_0-aucpr:0.57864 validation_1-aucpr:0.54771
[74] validation_0-aucpr:0.57882 validation_1-aucpr:0.54799
[75] validation_0-aucpr:0.57886 validation_1-aucpr:0.54800
[76] validation_0-aucpr:0.57913 validation_1-aucpr:0.54807
[77] validation_0-aucpr:0.57944 validation_1-aucpr:0.54816
[78] validation_0-aucpr:0.57971 validation_1-aucpr:0.54788
[79] validation_0-aucpr:0.57992 validation_1-aucpr:0.54787
[80] validation_0-aucpr:0.58000 validation_1-aucpr:0.54764
[81] validation_0-aucpr:0.58015 validation_1-aucpr:0.54770
[82] validation_0-aucpr:0.58060 validation_1-aucpr:0.54769
[83] validation_0-aucpr:0.58081 validation_1-aucpr:0.54795
[84] validation_0-aucpr:0.58111 validation_1-aucpr:0.54806
[85] validation_0-aucpr:0.58129 validation_1-aucpr:0.54839
[86] validation_0-aucpr:0.58155 validation_1-aucpr:0.54863
[87] validation_0-aucpr:0.58177 validation_1-aucpr:0.54847
[88] validation_0-aucpr:0.58203 validation_1-aucpr:0.54850
[89] validation_0-aucpr:0.58226 validation_1-aucpr:0.54878
[90] validation_0-aucpr:0.58237 validation_1-aucpr:0.54870
[91] validation_0-aucpr:0.58246 validation_1-aucpr:0.54869
[92] validation_0-aucpr:0.58257 validation_1-aucpr:0.54876
[93] validation_0-aucpr:0.58272 validation_1-aucpr:0.54874
[94] validation_0-aucpr:0.58277 validation_1-aucpr:0.54866
[95] validation_0-aucpr:0.58299 validation_1-aucpr:0.54879
[96] validation_0-aucpr:0.58342 validation_1-aucpr:0.54877
[97] validation_0-aucpr:0.58352 validation_1-aucpr:0.54866
[98] validation_0-aucpr:0.58359 validation_1-aucpr:0.54855
[99] validation_0-aucpr:0.58369 validation_1-aucpr:0.54852
Test Accuracy: 0.8992669362992922
 precision recall f1-score support
 0 0.92 0.96 0.94 6874
 1 0.66 0.47 0.55 1038
 accuracy 0.90 7912
 macro avg 0.79 0.72 0.75 7912
weighted avg 0.89 0.90 0.89 7912
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
# Initialize Random Forest with some reasonable defaults
rf_model = RandomForestClassifier(
 n_estimators=100,
 max_depth=10, # You can tune this later
 class_weight='balanced', # Helps handle class imbalance
 random_state=42,
 n_jobs=-1 # Use all CPU cores
)
# Train on training set
rf_model.fit(X_train, y_train)
# Predict on test set
y_pred_rf = rf_model.predict(X_test)
# Evaluate
accuracy_rf = accuracy_score(y_test, y_pred_rf)
print(f"\nRandom Forest Test Accuracy: {accuracy_rf:.4f}")
print(classification_report(y_test, y_pred_rf))
Random Forest Test Accuracy: 0.8830
 precision recall f1-score support
 0 0.93 0.93 0.93 6874
 1 0.56 0.54 0.55 1038
 accuracy 0.88 7912
 macro avg 0.74 0.74 0.74 7912
weighted avg 0.88 0.88 0.88 7912
import lightgbm as lgb
from sklearn.metrics import classification_report, accuracy_score
lgb_train = lgb.Dataset(X_train, label=y_train)
lgb_val = lgb.Dataset(X_val, label=y_val)
params = {
 'objective': 'binary',
 'metric': 'auc',
 'is_unbalance': True,
 'learning_rate': 0.1,
 'num_leaves': 31,
 'max_depth': -1,
 'seed': 42
}
model_lgb = lgb.train(
 params,
 lgb_train,
 valid_sets=[lgb_val],
 callbacks=[lgb.early_stopping(stopping_rounds=10)]
)
y_pred_lgb = (model_lgb.predict(X_test) > 0.5).astype(int)
print("LightGBM Test Accuracy:", accuracy_score(y_test, y_pred_lgb))
print(classification_report(y_test, y_pred_lgb))
[LightGBM] [Info] Number of positive: 4845, number of negative: 32076
[LightGBM] [Info] Auto-choosing row-wise multi-threading, the overhead
of testing was 0.000448 seconds.
You can set `force_row_wise=true` to remove the overhead.
And if memory is not enough, you can set `force_col_wise=true`.
[LightGBM] [Info] Total Bins 346
[LightGBM] [Info] Number of data points in the train set: 36921,
number of used features: 4
[LightGBM] [Info] [binary:BoostFromScore]: pavg=0.131226 ->
initscore=-1.890161
[LightGBM] [Info] Start training from score -1.890161
Training until validation scores don't improve for 10 rounds
Early stopping, best iteration is:
[35] valid_0's auc: 0.798434
LightGBM Test Accuracy: 0.8805611729019212
 precision recall f1-score support
 0 0.93 0.93 0.93 6874
 1 0.55 0.54 0.54 1038
 accuracy 0.88 7912
 macro avg 0.74 0.73 0.74 7912
weighted avg 0.88 0.88 0.88 7912
from catboost import CatBoostClassifier
cat_model = CatBoostClassifier(
 iterations=100,
 learning_rate=0.1,
 depth=6,
 eval_metric='AUC',
 random_seed=42,
 early_stopping_rounds=10,
 verbose=False
)
cat_model.fit(X_train, y_train, eval_set=(X_val, y_val))
y_pred_cat = cat_model.predict(X_test)
print("CatBoost Test Accuracy:", accuracy_score(y_test, y_pred_cat))
print(classification_report(y_test, y_pred_cat))
CatBoost Test Accuracy: 0.8945904954499494
 precision recall f1-score support
 0 0.92 0.96 0.94 6874
 1 0.63 0.47 0.54 1038
 accuracy 0.89 7912
 macro avg 0.78 0.71 0.74 7912
weighted avg 0.88 0.89 0.89 7912
pip install --upgrade scikit-learn xgboost catboost joblib
Requirement already satisfied: scikit-learn in
/usr/local/lib/python3.11/dist-packages (1.7.1)
Requirement already satisfied: xgboost in
/usr/local/lib/python3.11/dist-packages (3.0.3)
Requirement already satisfied: catboost in
/usr/local/lib/python3.11/dist-packages (1.2.8)
Requirement already satisfied: joblib in
/usr/local/lib/python3.11/dist-packages (1.5.1)
Requirement already satisfied: numpy>=1.22.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (1.26.4)
Requirement already satisfied: scipy>=1.8.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (1.15.3)
Requirement already satisfied: threadpoolctl>=3.1.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (3.6.0)
Requirement already satisfied: nvidia-nccl-cu12 in
/usr/local/lib/python3.11/dist-packages (from xgboost) (2.21.5)
Requirement already satisfied: graphviz in
/usr/local/lib/python3.11/dist-packages (from catboost) (0.21)
Requirement already satisfied: matplotlib in
/usr/local/lib/python3.11/dist-packages (from catboost) (3.7.2)
Requirement already satisfied: pandas>=0.24 in
/usr/local/lib/python3.11/dist-packages (from catboost) (2.2.3)
Requirement already satisfied: plotly in
/usr/local/lib/python3.11/dist-packages (from catboost) (5.24.1)
Requirement already satisfied: six in /usr/local/lib/python3.11/distpackages (from catboost) (1.17.0)
Requirement already satisfied: mkl_fft in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (1.3.8)
Requirement already satisfied: mkl_random in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (1.2.4)
Requirement already satisfied: mkl_umath in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (0.1.1)
Requirement already satisfied: mkl in /usr/local/lib/python3.11/distpackages (from numpy>=1.22.0->scikit-learn) (2025.2.0)
Requirement already satisfied: tbb4py in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (2022.2.0)
Requirement already satisfied: mkl-service in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (2.4.1)
Requirement already satisfied: python-dateutil>=2.8.2 in
/usr/local/lib/python3.11/dist-packages (from pandas>=0.24->catboost)
(2.9.0.post0)
Requirement already satisfied: pytz>=2020.1 in
/usr/local/lib/python3.11/dist-packages (from pandas>=0.24->catboost)
(2025.2)
Requirement already satisfied: tzdata>=2022.7 in
/usr/local/lib/python3.11/dist-packages (from pandas>=0.24->catboost)
(2025.2)
Requirement already satisfied: contourpy>=1.0.1 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(1.3.2)
Requirement already satisfied: cycler>=0.10 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(0.12.1)
Requirement already satisfied: fonttools>=4.22.0 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(4.58.4)
Requirement already satisfied: kiwisolver>=1.0.1 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(1.4.8)
Requirement already satisfied: packaging>=20.0 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(25.0)
Requirement already satisfied: pillow>=6.2.0 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(11.2.1)
Requirement already satisfied: pyparsing<3.1,>=2.3.1 in
/usr/local/lib/python3.11/dist-packages (from matplotlib->catboost)
(3.0.9)
Requirement already satisfied: tenacity>=6.2.0 in
/usr/local/lib/python3.11/dist-packages (from plotly->catboost)
(8.5.0)
Requirement already satisfied: intel-openmp<2026,>=2024 in
/usr/local/lib/python3.11/dist-packages (from mkl->numpy>=1.22.0-
>scikit-learn) (2024.2.0)
Requirement already satisfied: tbb==2022.* in
/usr/local/lib/python3.11/dist-packages (from mkl->numpy>=1.22.0-
>scikit-learn) (2022.2.0)
Requirement already satisfied: tcmlib==1.* in
/usr/local/lib/python3.11/dist-packages (from tbb==2022.*->mkl-
>numpy>=1.22.0->scikit-learn) (1.4.0)
Requirement already satisfied: intel-cmplr-lib-rt in
/usr/local/lib/python3.11/dist-packages (from mkl_umath-
>numpy>=1.22.0->scikit-learn) (2024.2.0)
Requirement already satisfied: intel-cmplr-lib-ur==2024.2.0 in
/usr/local/lib/python3.11/dist-packages (from intelopenmp<2026,>=2024->mkl->numpy>=1.22.0->scikit-learn) (2024.2.0)
Note: you may need to restart the kernel to use updated packages.
pip install --upgrade scikit-learn
Requirement already satisfied: scikit-learn in
/usr/local/lib/python3.11/dist-packages (1.7.1)
Requirement already satisfied: numpy>=1.22.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (1.26.4)
Requirement already satisfied: scipy>=1.8.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (1.15.3)
Requirement already satisfied: joblib>=1.2.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (1.5.1)
Requirement already satisfied: threadpoolctl>=3.1.0 in
/usr/local/lib/python3.11/dist-packages (from scikit-learn) (3.6.0)
Requirement already satisfied: mkl_fft in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (1.3.8)
Requirement already satisfied: mkl_random in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (1.2.4)
Requirement already satisfied: mkl_umath in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (0.1.1)
Requirement already satisfied: mkl in /usr/local/lib/python3.11/distpackages (from numpy>=1.22.0->scikit-learn) (2025.2.0)
Requirement already satisfied: tbb4py in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (2022.2.0)
Requirement already satisfied: mkl-service in
/usr/local/lib/python3.11/dist-packages (from numpy>=1.22.0->scikitlearn) (2.4.1)
Requirement already satisfied: intel-openmp<2026,>=2024 in
/usr/local/lib/python3.11/dist-packages (from mkl->numpy>=1.22.0-
>scikit-learn) (2024.2.0)
Requirement already satisfied: tbb==2022.* in
/usr/local/lib/python3.11/dist-packages (from mkl->numpy>=1.22.0-
>scikit-learn) (2022.2.0)
Requirement already satisfied: tcmlib==1.* in
/usr/local/lib/python3.11/dist-packages (from tbb==2022.*->mkl-
>numpy>=1.22.0->scikit-learn) (1.4.0)
Requirement already satisfied: intel-cmplr-lib-rt in
/usr/local/lib/python3.11/dist-packages (from mkl_umath-
>numpy>=1.22.0->scikit-learn) (2024.2.0)
Requirement already satisfied: intel-cmplr-lib-ur==2024.2.0 in
/usr/local/lib/python3.11/dist-packages (from intelopenmp<2026,>=2024->mkl->numpy>=1.22.0->scikit-learn) (2024.2.0)
Note: you may need to restart the kernel to use updated packages.
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from catboost import CatBoostClassifier
from sklearn.metrics import classification_report, accuracy_score
# Define base learners
base_learners = [
 ('rf', RandomForestClassifier(n_estimators=100,
class_weight='balanced', random_state=42)),
 ('xgb', XGBClassifier(eval_metric='logloss', scale_pos_weight=2.5,
random_state=42)), # Removed use_label_encoder
 ('cat', CatBoostClassifier(verbose=0, random_state=42))
]
# Define meta learner
meta_learner = LogisticRegression()
# Create the stacking classifier
stacking_model = StackingClassifier(
 estimators=base_learners,
 final_estimator=meta_learner,
 cv=5, # Cross-validation folds
 n_jobs=-1
)
# Train
stacking_model.fit(X_train, y_train)
# Predict
y_pred_stack = stacking_model.predict(X_test)
# Evaluate
print("Stacked Model Accuracy:", accuracy_score(y_test, y_pred_stack))
print(classification_report(y_test, y_pred_stack))
Stacked Model Accuracy: 0.9035642062689585
 precision recall f1-score support
 0 0.92 0.98 0.95 6874
 1 0.73 0.42 0.53 1038
 accuracy 0.90 7912
 macro avg 0.82 0.70 0.74 7912
weighted avg 0.89 0.90 0.89 7912
