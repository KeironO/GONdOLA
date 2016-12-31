import collections, itertools, csv
from orangecontrib.bio import go

phenotype_data = open("./data/phenotype_data.tab", "r")

gene_label_dictonary = collections.defaultdict(set)

for line in phenotype_data.read().splitlines():
    gene_phenotype_info = line.split("\t")
    if gene_phenotype_info[9] == "viable":
        gene_label_dictonary["viable"].update([gene_phenotype_info[3]])
    elif gene_phenotype_info[9] == "inviable":
        gene_label_dictonary["inviable"].update([gene_phenotype_info[3]])

gene_association = open("./data/gene_association.sgd", "r")

gene_ontology = go.Ontology()
ancestor_dictonary = collections.defaultdict(set)
for line in gene_association.readlines():
    if line.startswith("!") == False:
        gene_info = line.split("\t")
        if gene_info[2] in list(set(list(itertools.chain.from_iterable([gene_label_dictonary[x] for x in gene_label_dictonary])))):
            ancestor_dictonary[gene_info[1]].add(gene_info[4])
            try:
                parents = list(gene_ontology.extract_super_graph(gene_info[4]))
                for p in parents:
                    ancestor_dictonary[gene_info[1]].add(p)
            except KeyError:
                pass

all_go_terms = list(set(list(itertools.chain.from_iterable([ancestor_dictonary[x] for x in ancestor_dictonary]))))

with open("./data/yeast_data.csv", "wb") as csv_file:
    writer = csv.writer(csv_file, delimiter=",")
    writer.writerow(["Gene ID", "Outcome"] + all_go_terms)
    for gene in ancestor_dictonary.keys():
        gene_seen = []
        for term in all_go_terms:
            if term in ancestor_dictonary[gene]:
                gene_seen.append(1)
            else:
                gene_seen.append(0)
        if gene in gene_label_dictonary["inviable"]:
            writer.writerow([gene, "inviable"] + gene_seen)
        else:
            writer.writerow([gene, "viable"] + gene_seen)
csv_file.close()