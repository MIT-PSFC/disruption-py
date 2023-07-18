python3 generate_datasets.py --unique_id thesis_dataset1

python3 train.py --train_path output/train_thesis_dataset1.csv --test_path output/test_thesis_dataset1.csv --unique_id thesis_rf1

python3 train.py --train_path output/train_thesis_dataset1.csv --test_path output/test_thesis_dataset1.csv --unique_id thesis_lr1 --model LogisticRegression

python3 train.py --train_path output/train_thesis_dataset1.csv --test_path output/test_thesis_dataset1.csv --unique_id thesis_svm1 --model SupportVectorMachine

python3 eval.py output/train_thesis_dataset1.csv --unique_id=paper_model1

python3 eval.py output/train_thesis_dataset1.csv --unique_id=thesis_rf1 --model_path output/RandomForestClassifier2023-07-18_thesis_rf1.joblib
