interactome = 'interactome-flybase-collapsed-weighted.txt'
flybase_interactome = 'flybase_interactome.txt'
source_gene_file = 'source_gene.txt'
flybase_dict = 'flybase_dict.txt'
myoblast_fusion_components = 'myoblast_fusion_components.txt'
gene_dict = {}

with open(interactome, 'r') as f:
    with open(flybase_interactome, 'w') as f2:
        with open(flybase_dict, 'w') as f3:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                # split the line by tabs
                temp = line.strip().split('\t')
                if temp[0] == 'nan':
                    temp[0] = 'TRPV'
                if temp[1] == 'nan':
                    temp[1] = 'TRPV'
                if temp[4] not in gene_dict:
                    gene_dict[temp[4]] = temp[0]
                    f3.write(temp[4] + '\t' + temp[0] + '\n')
                if temp[5] not in gene_dict:
                    gene_dict[temp[5]] = temp[1]
                    f3.write(temp[5] + '\t' + temp[1] + '\n')
                f2.write(temp[0] + '\t' + temp[1] + '\t' + temp[2] + '\n')


with open(myoblast_fusion_components, 'r') as f:
    with open(source_gene_file, 'w') as f2:
        f2.write('NODEID\tprize\n')
        for line in f:
            # split the line by tabs
            temp = line.strip().split('\t')
            f2.write(gene_dict[temp[1]] + '\t' + '10' + '\n')


