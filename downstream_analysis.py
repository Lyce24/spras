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
oi1_source_output_file = Path('./downstream_analysis/source_omicsintegrator1.txt')
oi2_output_file = Path('./downstream_analysis/omicsintegrator2.txt')
oi2_source_output_file = Path('./downstream_analysis/source_omicsintegrator2.txt')

result_oi1_dict = {}
result_oi2_dict = {}
result_source_oi1_dict = {}
result_source_oi2_dict = {}

oi1_runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator1-params')]
oi2_runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator2-params')]

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

def generate_diffusion_analysis(analysis_file: Path, gene_analysis_file, runs, result_dict, result_gene_dict) -> None:
    for run in runs:
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
            if gene not in result_dict:
                result_dict[gene] = 1
            else:
                result_dict[gene] += 1
        for gene in source_temp_set:
            if gene not in result_gene_dict:
                result_gene_dict[gene] = 1
            else:
                result_gene_dict[gene] += 1

    sorted_list = sorted(result_dict.items(), key=lambda x: x[1], reverse=True)
    sorted_source_list = sorted(result_gene_dict.items(), key=lambda x: x[1], reverse=True)

    with open(analysis_file, 'w') as f:
        f.write('gene\tcount\tpercentage\treceptor\n')
        for item in sorted_list:
            if item[0] in receptor_gene:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'True' + '\n')
            else:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'False' + '\n')
                
    with open(gene_analysis_file, 'w') as f:
        f.write('gene\tcount\tpercentage\n')
        for item in sorted_source_list:
            f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4))  + '\n')
            
def generate_fusion_analysis(process, algo) -> None:
    analysis_file = Path(f'./downstream_analysis/{process}-{algo}.txt')
    gene_analysis_file = Path(f'./downstream_analysis/{process}-{algo}-source.txt')
    runs = [x for x in parent_dir.iterdir() if x.is_dir() and x.name.startswith(f'{process}-{algo}-params')]
    if process == 'cell-cell-fusion':
        target_gene = cell_target_gene
    elif process == 'muscle-development':
        target_gene = muscle_target_gene
        
    result_dict = {}
    result_gene_dict = {}
    
    for run in runs:
        temp_set = set()
        gene_temp_set = set() 
        with open(run / 'pathway.txt') as f:
            for line in f:
                line = line.strip()
                endpoints = line.split("\t")
                if endpoints[0] in source_gene or endpoints[0] in target_gene:
                    gene_temp_set.add(endpoints[0])
                else:
                    temp_set.add(endpoints[0])
                if endpoints[1] in source_gene or endpoints[1] in target_gene:
                    gene_temp_set.add(endpoints[1])
                else:
                    temp_set.add(endpoints[1])

        for gene in temp_set:
            if gene not in result_dict:
                result_dict[gene] = 1
            else:
                result_dict[gene] += 1
        for gene in gene_temp_set:
            if gene not in result_gene_dict:
                result_gene_dict[gene] = 1
            else:
                result_gene_dict[gene] += 1

    sorted_list = sorted(result_dict.items(), key=lambda x: x[1], reverse=True)
    sorted_gene_list = sorted(result_gene_dict.items(), key=lambda x: x[1], reverse=True)
    
    with open(analysis_file, 'w') as f:
        f.write('gene\tcount\tpercentage\treceptor\n')
        for item in sorted_list:
            if item[0] in receptor_gene:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'True' + '\n')
            else:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'False' + '\n')
                
    with open(gene_analysis_file, 'w') as f:
        f.write('gene\tcount\tpercentage\tgene\n')
        for item in sorted_gene_list:
            if item[0] in source_gene:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'Source' + '\n')
            else:
                f.write(item[0] + '\t' + str(item[1]) + '\t' + str(round((item[1] / len(runs)), 4)) + '\t' + 'Target' + '\n')

def generate_runs_info(process, algo):
    analysis_file = Path(f'./downstream_analysis/{process}-{algo}-runs.txt')
    if process == 'cell-cell-fusion' or process == 'muscle-development':
        runs = [x for x in parent_dir.iterdir() if x.is_dir() and x.name.startswith(f'{process}-{algo}-params')]
    elif process == 'flybase':
        runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir() and x.name.startswith(f'{process}-{algo}-params')]
        
    if process == 'cell-cell-fusion':
        target_gene = cell_target_gene
    elif process == 'muscle-development':
        target_gene = muscle_target_gene
    else:
        target_gene = set()
    
    run_dict = {}
    
    for run in runs:
        nodes = set()
        gene_temp_set = set()
        receptor_temp_set = set() 
        with open(run / 'pathway.txt') as f:
            for line in f:
                line = line.strip()
                endpoints = line.split("\t")
                nodes.add(endpoints[0])
                nodes.add(endpoints[1])
                if endpoints[0] in source_gene or endpoints[0] in target_gene:
                    gene_temp_set.add(endpoints[0])
                elif endpoints[0] in receptor_gene:
                    receptor_temp_set.add(endpoints[0])
                if endpoints[1] in source_gene or endpoints[1] in target_gene:
                    gene_temp_set.add(endpoints[1])
                elif endpoints[1] in receptor_gene:
                    receptor_temp_set.add(endpoints[1])
        run_dict[str(run).split('-')[-1]] = (len(nodes), len(receptor_temp_set), len(gene_temp_set))

    # sorted the dictionary by the number of receptors
    sorted_runs = sorted(run_dict.items(), key=lambda x: x[1][1], reverse=True)

    with open(analysis_file, 'w') as f:
        f.write('RunID\t# of Nodes\tReceptor Nodes\tPrize Nodes\n')
        for item in sorted_runs:
            f.write(item[0] + '\t' + str(item[1][0]) + '\t' + str(item[1][1]) + '\t' + str(item[1][2]) + '\n')
            
generate_nodes(source_gene_file, source_gene)
generate_nodes(receptor_gene_file, receptor_gene)
generate_nodes(cell_target_gene_file, cell_target_gene)
generate_nodes(muscle_target_gene_file, muscle_target_gene)
generate_diffusion_analysis(oi1_output_file, oi1_source_output_file, oi1_runs, result_oi1_dict, result_source_oi1_dict)
generate_diffusion_analysis(oi2_output_file, oi2_source_output_file, oi2_runs, result_oi2_dict, result_source_oi2_dict)

generate_runs_info('flybase', 'omicsintegrator1')
generate_runs_info('flybase', 'omicsintegrator2')

for i in ['cell-cell-fusion', 'muscle-development']:
    for j in ['pathlinker', 'rwr']:
        generate_fusion_analysis(i, j)
        generate_runs_info(i, j)