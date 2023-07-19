# Generate Datasets
python3 generate_datasets.py --unique_id thesis_dataset1
python3 generate_datasets.py --unique_id thesis_dataset2 --feature_cols filtered_features.txt

# Train Models
python3 train.py --train_path output/train_thesis_dataset1.csv --test_path output/test_thesis_dataset1.csv --unique_id thesis_rf1
python3 train.py --train_path output/train_thesis_dataset1.csv --test_path output/test_thesis_dataset1.csv --unique_id thesis_lr1 --model LogisticRegression

python3 train.py --train_path output/train_thesis_dataset2.csv --test_path output/test_thesis_dataset2.csv --unique_id thesis_rf2
python3 train.py --train_path output/train_thesis_dataset2.csv --test_path output/test_thesis_dataset2.csv --unique_id thesis_lr2 --model LogisticRegression

# Evaluate Models
python3 eval.py output/train_thesis_dataset1.csv --unique_id=paper_model1

# Random Forest
python3 eval.py output/test_thesis_dataset1.csv --unique_id=thesis_test_rf1 --model_path output/RandomForestClassifier2023-07-18_thesis_rf1.joblib --order=DATASET_ORDER
python3 eval.py output/train_thesis_dataset1.csv --unique_id=thesis_train_rf1 --model_path output/RandomForestClassifier2023-07-18_thesis_rf1.joblib --order=DATASET_ORDER

python3 eval.py output/train_thesis_dataset2.csv --unique_id=thesis_train_rf2 --model_path output/RandomForestClassifier2023-07-18_thesis_rf2.joblib --order=DATASET_ORDER
python3 eval.py output/test_thesis_dataset2.csv --unique_id=thesis_test_rf2 --model_path output/RandomForestClassifier2023-07-18_thesis_rf2.joblib --order=DATASET_ORDER

# Logistic Regression
python3 eval.py output/train_thesis_dataset1.csv --unique_id=thesis_train_lr1 --model_path output/LogisticRegression2023-07-18_thesis_lr1.joblib --order=DATASET_ORDER
python3 eval.py output/test_thesis_dataset1.csv --unique_id=thesis_test_lr1 --model_path output/LogisticRegression2023-07-18_thesis_lr1.joblib --order=DATASET_ORDE

python3 eval.py output/train_thesis_dataset2.csv --unique_id=thesis_train_lr2 --model_path output/LogisticRegression2023-07-18_thesis_lr2.joblib --order=DATASET_ORDER
python3 eval.py output/test_thesis_dataset2.csv --unique_id=thesis_test_lr2 --model_path output/LogisticRegression2023-07-18_thesis_lr2.joblib --order=DATASET_ORDER