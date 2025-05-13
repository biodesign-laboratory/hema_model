import numpy as np
import csv
from pathlib import Path
import pandas as pd

indexes = [-1*(i - 42) for i in range(1, 42)]

results = np.vstack((np.array(['string_placeholder' for i in range(41)]), np.array(['string_placeholder' for i in range(41)]), np.array(['string_placeholder' for i in range(41)])))
#print(results.shape)

for i in range(1,4):

    sensitivity_loc = Path.cwd() / 'Model_3_SA' / f'final_morris_chronic_{i}' / 'mu_star_vs_sigma' / 'results.csv'

    with sensitivity_loc.open('r', newline='') as file:

        csv_reader = csv.reader(file)
        data = list(csv_reader)

    results[i-1] = [row[1] for row in data[1:]]

top_10_ranked = list(zip(indexes[-10:], results[0][-10:], results[1][-10:], results[2][-10:]))
df = pd.DataFrame(np.array(top_10_ranked)[:, 1:], columns=['Exp_1', 'Exp_2', 'Exp_3'], index=list(map(str, indexes[-10:])))
df.to_csv(Path.cwd() / 'morris_chronic_results_combined.csv')
print("Results saved.")

unique_params = np.unique(np.array(top_10_ranked)[:, 1:])
print(unique_params)

