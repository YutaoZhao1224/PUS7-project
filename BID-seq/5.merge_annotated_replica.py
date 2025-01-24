
### usage: python 5.merge_annotated_replica.py mRNA-shctr-BID-seq.pU_site.annotated.txt mRNA-shctr-BID-seq.pU_site.txt mRNA-shN-myc-BID-seq.pU_site.txt

import sys

site_info = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        gene_info = line.split("\t")
        info = gene_info[0] + "\t" + gene_info[1] + "\t" + gene_info[2] + "\t" + gene_info[3] + "\t" + \
            gene_info[4] + "\t" + gene_info[5] + "\t" + gene_info[6] + "\t" + gene_info[7] + "\t" + \
            gene_info[8] + "\t" + gene_info[9]
        site_info[gene_info[0]+gene_info[1]+gene_info[2]] = info

argv_num = len(sys.argv)
i = argv_num
site_mutation = []
while (i > 2):
    mutation_dir = {}
    with open(sys.argv[argv_num - i + 2]) as f:
        for line in f:
            mutation_info = line.split("\t")
            mutation_dir[mutation_info[0]+mutation_info[2]+mutation_info[5].strip()] = mutation_info[3] + "\t" + mutation_info[4]
    site_mutation.append(mutation_dir)
    i = i - 1

a = site_mutation[0].keys()
b = []
j = 0
while(j < argv_num - 3):
    c = []
    b = site_mutation[j+1].keys()
    for k in a:
        if k in b:
            c.append(k)
    a = c
    j = j + 1

keys = site_info.keys()
final_overlap_key = []
for i in a:
    if i in keys:
        final_overlap_key.append(i)

for i in final_overlap_key:
    a = ''
    for j in site_mutation:
        a = a + j[i] + "\t"
    print(site_info[i] + "\t" + a)








