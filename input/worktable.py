interactome = 'interactome-flybase-collapsed-weighted.txt'
flybase_interactome = 'flybase_interactome.txt'
source_gene_file = 'source_gene.txt'
flybase_dict = 'flybase_dict.txt'
myoblast_fusion_components = 'myoblast_fusion_components.txt'
cell_cell_fusion_file = 'cell-cell-fusion.txt'
muscle_development_file = 'muscle-development.txt'
cell_target_gene_file = 'cell_target_gene.txt'
muscle_target_gene_file = 'muscle_target_gene.txt'

source_gene = set()
gene_dict = {}

cell_cell_fusion_genes = set()
muscle_gene = set()

cell_valid_genes = set()
muscle_valid_genes = set()

cell_target_gene = set()
muscle_target_gene = set()

with open(source_gene_file, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        temp = line.strip().split('\t')
        source_gene.add(temp[0])
        
with open(flybase_dict, 'r') as f:
    for line in f:
        temp = line.strip().split('\t')
        gene_dict[temp[0]] = temp[1]

with open(cell_cell_fusion_file, 'r') as f:
    for line in f:
        temp = line.strip().split('\t')
        flybase_id = temp[0].split(':')[1]
        cell_cell_fusion_genes.add(flybase_id)
        
with open(muscle_development_file, 'r') as f:
    for line in f:
        temp = line.strip().split('\t')
        flybase_id = temp[0].split(':')[1]
        muscle_gene.add(flybase_id)
        
for i in cell_cell_fusion_genes:
    if i in gene_dict:
        cell_valid_genes.add(gene_dict[i])
    
for i in muscle_gene:
    if i in gene_dict:
        muscle_valid_genes.add(gene_dict[i])
    
print(len(cell_valid_genes))
print(len(muscle_valid_genes))
    
for i in cell_valid_genes:
    if i not in source_gene:
        cell_target_gene.add(i)
        
for i in muscle_valid_genes:
    if i not in source_gene:
        muscle_target_gene.add(i)
    
print(len(cell_target_gene))
print(len(muscle_target_gene))

with open(cell_target_gene_file, 'w') as f:
    with open(muscle_target_gene_file, 'w') as g:
        f.write("NODEID\tprize\n")
        g.write("NODEID\tprize\n")
        
        for i in cell_target_gene:
            f.write(f"{i}\t10\n")
        
        for i in muscle_target_gene:
            g.write(f"{i}\t10\n")