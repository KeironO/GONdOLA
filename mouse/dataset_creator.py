import collections, json, tqdm, itertools, csv
from orangecontrib.bio import go
from tqdm import tqdm

def get_ancestors(gene_association, interesting_genes):
    gene_ontology = go.Ontology()
    ancestor_dictonary = collections.defaultdict(set)
    print "Getting parents of GO terms within the gene association file."
    for line in tqdm(gene_association.readlines()):
        if line.startswith("!") == False:
            gene_info = line.split("\t")
            if gene_info[1] in interesting_genes:
                ancestor_dictonary[gene_info[1]].add(gene_info[4])
                try:
                    parents = list(gene_ontology.extract_super_graph(gene_info[4]))
                    for p in parents:
                        ancestor_dictonary[gene_info[1]].add(p)
                except KeyError:
                    pass

    ancestor_dictonary = collections.defaultdict(list, ((k, list(v)) for k, v in ancestor_dictonary.items()))
    all_go_terms = list(set(list(itertools.chain.from_iterable([ancestor_dictonary[x] for x in ancestor_dictonary]))))
    print "\n\t-> There were a total of %i genes found within this dataset" % len(ancestor_dictonary.keys())
    print "\t-> Took forward a total of %i GO terms \n" % len(all_go_terms)
    return ancestor_dictonary, all_go_terms


def create_final_files(ancestor_dictonary, all_go_terms, gene_label_dict):
    with open("./data/genes_to_go_terms.json", "wb") as out_file:
        json.dump(ancestor_dictonary, out_file)

    with open("./data/mouse_data.csv", "wb") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        writer.writerow(["Gene ID", "Outcome"] + all_go_terms)
        for gene in tqdm(ancestor_dictonary.keys()):
            gene_seen = []
            for term in all_go_terms:
                if term in ancestor_dictonary[gene]:
                    gene_seen.append(1)
                else:
                    gene_seen.append(0)
            if gene in gene_label_dict["lethal"]:
                writer.writerow([gene, "lethal"] + gene_seen)
            else:
                writer.writerow([gene, "viable"] + gene_seen)
    csv_file.close()

def get_labels(phenogenomp, lethal_terms):
    gene_label_dict = collections.defaultdict(set)
    removed_list = []
    for line in phenogenomp.read().splitlines():
        gene_info = line.split("\t")
        # Checking if the third column is a member of the MP ontology, just in case.
        # We are only interested in single-gene mutation
        # -> Not a good idea.
        #   -> George: Both knocked out and became lethal/viable
        #   -> One GO area that gets knocked out.
        #   -> Learning from lethal GO terms that are multipoint.
        #   -> Not enough data.

        if gene_info[3].startswith("MP:") and len(list(gene_info[5].split(","))) == 1:
            if gene_info[3] in lethal_terms:
                gene_label_dict["lethal"].update([x for x in gene_info[5].split(",")])
            else:
                gene_label_dict["viable"].update([x for x in gene_info[5].split(",")])
        else:
            removed_list.append([x for x in gene_info[5].split(",")])
    print "We have removed %i number of genes" %len(removed_list)
    gene_label_dict = collections.defaultdict(list, ((k, list(v)) for k, v in gene_label_dict.items()))
    return gene_label_dict

if __name__ == "__main__":
    phenogenomp = open("./data/MGI_PhenoGenoMP.rpt", "r")

    # Taken from previous work.
    lethal_terms = ["MP:0011400","MP:0011083","MP:0011089","MP:0011087","MP:0011085","MP:0008028","MP:0002083",
        "MP:0011091", "MP:0011092", "MP:0011098", "MP:0011093", "MP:0011094", "MP:0011095", "MP:0011096", "MP:0011097",
        "MP:0011099", "MP:0011100", "MP:0011111", "MP:0013139"]

    gene_label_dict = get_labels(phenogenomp, lethal_terms)

    gene_association = open("./data/gene_association.mgi", "r")
    ancestor_dictonary, all_go_terms = get_ancestors(gene_association,
                  # HORRIBLE LAZY HACK, PLEASE FORGIVE ME.
                  list(set(list(itertools.chain.from_iterable([gene_label_dict[x] for x in gene_label_dict]))))
                  )
    create_final_files(ancestor_dictonary, all_go_terms, gene_label_dict)