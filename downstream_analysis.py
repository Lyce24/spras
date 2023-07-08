from pathlib import Path

source_gene_file = Path('./input/source_gene.txt')
receptor_gene_file = Path('./input/receptor.txt')
parent_dir = Path("./output/amigo2")
flybase_parent_dir = Path("./FlyBase")
# Check whether the directory FlyBase exists

cell_target_gene_file = Path('./input/cell_target_gene.txt')
muscle_target_gene_file = Path('./input/muscle_target_gene.txt')

source_gene = set()
receptor_gene = set()
cell_target_gene = set()
muscle_target_gene = set()

oi1_output_file = Path('./downstream_analysis/omicsintegrator1.txt')
oi1_source_file = Path('./downstream_analysis/source_omicsintegrator1.txt')

oi1_dict = {}
oi2_dict = {}
cell_pathlinker_dict = {}
cell_rwr_dict = {}
muscle_pathlinker_dict = {}
muscle_rwr_dict = {}

reseult_oi1_dict = {}
reseult_source_oi1_dict = {}

oi1_runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator1-params')]

def generate_nodes(nodes_file: Path, gene_set) -> None:
    """
    This function is for generating the nodes from the path to the source/target file
    """
    with open(nodes_file, 'r') as f:
        for i, line in enumerate(f):
            # if the first line is the title, skip it
            if i == 0:
                continue
            line = line.strip()
            endpoints = line.split("\t")
            gene_set.add(endpoints[0])

generate_nodes(source_gene_file, source_gene)
generate_nodes(receptor_gene_file, receptor_gene)


for run in oi1_runs:
    temp_set = set()
    source_temp_set = set()
    with open(run / 'pathway.txt') as f:
        for line in f:
            line = line.strip()
            endpoints = line.split("\t")
            if endpoints[0] not in source_gene:
                temp_set.add(endpoints[0])
            else:
                source_temp_set.add(endpoints[0])
            if endpoints[1] not in source_gene:
                temp_set.add(endpoints[1])
            else:
                source_temp_set.add(endpoints[1])
    for gene in temp_set:
        if gene not in reseult_oi1_dict:
            reseult_oi1_dict[gene] = 1
        else:
            reseult_oi1_dict[gene] += 1
    for gene in source_temp_set:
        if gene not in reseult_source_oi1_dict:
            reseult_source_oi1_dict[gene] = 1
        else:
            reseult_source_oi1_dict[gene] += 1
        

sorted_oi1_list = sorted(reseult_oi1_dict.items(), key=lambda x: x[1], reverse=True)
sorted_source_oi1_list = sorted(reseult_source_oi1_dict.items(), key=lambda x: x[1], reverse=True)


with open(oi1_output_file, 'w') as f:
    f.write('gene\tcount\tpercentage\treceptor\n')
    for item in sorted_oi1_list:
        if item[0] in receptor_gene:
            f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(oi1_runs)), 4) * 100) + '\t' + 'True' + '\n')
        else:
            f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(oi1_runs)), 4) * 100) + '\t' + 'False' + '\n')
            
with open(oi1_source_file, 'w') as f:
    f.write('gene\tcount\tpercentage\n')
    for item in sorted_source_oi1_list:
        f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(oi1_runs)), 4) * 100)  + '\n')