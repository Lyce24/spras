egfr = "phosphosite-irefindex13.0-uniprot.txt"
prizes = "tps-egfr-prizes.txt"

network_tsv = 'network_1.tsv'
index_gene = 'network_1_index_gene.tsv'
network_edge_list = 'network_1_edge_list.tsv'
score_1_tsv = 'scores_1.tsv'

index_gene_dict = {}
gene_set = set()

with open(egfr, 'r') as f1:
    for line in f1:
        temp = line.strip().split('\t')
        gene_set.add(temp[0])
        gene_set.add(temp[1])
        
with open(index_gene, 'w') as f:
    for i, gene in enumerate(gene_set):
        f.write(str(i) + "\t" + gene + "\n")
        index_gene_dict[gene] = i

with open(prizes, 'r') as f1:
        with open(score_1_tsv, "w") as f3:
            for i, line in enumerate(f1):
                if i == 0:
                    continue
                temp = line.strip().split('\t')
                
                f3.write(temp[0] + "\t" + temp[1] + "\n")

with open(egfr, 'r') as f1:
    with open(network_tsv, 'w') as f2:
        with open(network_edge_list, 'w') as f3:
            for line in f1:
                temp = line.strip().split('\t')
                f2.write(temp[0] + "\t" + temp[1] + "\n")
                f3.write(str(index_gene_dict[temp[0]]) + "\t" + str(index_gene_dict[temp[1]]) + "\n")
                    