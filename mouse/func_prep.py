import pandas as pd, csv
from tqdm import tqdm

df = pd.read_pickle("./data/mouse_data.pkl")

func = []

for gene, row in tqdm(df.iterrows()):
    # gene = Gene Name
    if row["Outcome"] == "viable":
        strat = 0
    else:
        strat = 1
    for go_term, go in row.drop("Outcome").iteritems():
        if go == 1:
            func.append([gene, go_term, strat])

with open("./data/mouse_func_data.txt", "wb") as func_output:
    writer = csv.writer(func_output, delimiter="\t", lineterminator="\n")
    for i in tqdm(func):
        writer.writerow(i)
func_output.close()
