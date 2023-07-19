# Runs
## generate_dataset_83mjrcfu -> train_guil2kwm -> 
Ran this to generate data for presentation. Used the namespace.pkl data from SORI2020 for training, but generated data for eval using generate_dataset.

# RF1 (disruptivity=.85)
## Train
Disruptions warned: 141/158 (89%)
Missed disruptions: 17/158 (10%)
False Alarms: 91/742 (12%)
Confusion Matrix Percentages [[0.72333333 0.10111111]
 [0.01888889 0.15666667]]
F2-Score: 0.8159722222222222
## Test
Disruptions warned: 20/40 (50%)
Missed disruptions: 20/40 (50%)
False Alarms: 28/185 (15%)
Confusion Matrix Percentages [[0.69777778 0.12444444]
 [0.08888889 0.08888889]]
F2-Score: 0.4807692307692308
# LR1 (disruptivity=.85)
## Train
Disruptions warned: 87/158 (55%)
Missed disruptions: 71/158 (44%)
False Alarms: 444/742 (59%)
Confusion Matrix Percentages [[0.33111111 0.49333333]
 [0.07888889 0.09666667]]
F2-Score: 0.37403267411865865
## Test
Disruptions warned: 18/40 (45%)
Missed disruptions: 22/40 (55%)
False Alarms: 113/185 (61%)
Confusion Matrix Percentages [[0.32       0.50222222]
 [0.09777778 0.08      ]]
F2-Score: 0.30927835051546393
# RF2 (disruptivity=.85)
## Train
Disruptions warned: 157/162 (96%)
Missed disruptions: 5/162 (3%)
False Alarms: 2/738 (0%)
Confusion Matrix Percentages [[0.81777778 0.00222222]
 [0.00555556 0.17444444]]
F2-Score: 0.9727385377942999
## Test
Disruptions warned: 16/36 (44%)
Missed disruptions: 20/36 (55%)
False Alarms: 1/189 (0%)
Confusion Matrix Percentages [[0.83555556 0.00444444]
 [0.08888889 0.07111111]]
F2-Score: 0.4968944099378882


