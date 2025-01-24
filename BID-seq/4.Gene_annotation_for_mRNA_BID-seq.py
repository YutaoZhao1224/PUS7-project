
### usage: cat mRNA-shctr-BID-seq.pU_site.txt | python ~/Genome/tools_custom/BID-seq/Gene_annotation_for_mRNA_BID-seq.py 1.txt > mRNA-shctr-BID-seq.pU_site.annotated.txt
from Bio import SeqIO
import sys

# make a seqence dictionary
genome_seq_dict = SeqIO.to_dict(SeqIO.parse('/home/yutaozhao/Genome/hg38_UCSC.fa', "fasta"))


# define a class for gene annotation file
class gene_annotation():
    def __init__(self, refseq_id, Chr, gene_direction, txStart_site, txEnd_site, cdsStart_site, cdsEnd_site, exon_list,
                 gene_name,gene_class):
        self.refseq_id = refseq_id
        self.Chr = Chr
        self.gene_direction = gene_direction
        self.txStart_site = int(txStart_site)
        self.txEnd_site = int(txEnd_site)
        self.cdsStart_site = int(cdsStart_site)
        self.cdsEnd_site = int(cdsEnd_site)
        self.exon_list = exon_list
        self.gene_name = gene_name
        self.gene_class = gene_class

### to get relative coordinates of CDS
    def mRNA_transcript(self, pos):
        pos_relative = -1 ### using this to judge if pos is in exons
        if self.gene_direction == '+':
            exon_len = []
            for i in range(0, int(len(self.exon_list)/2)):
                exon_len.append(self.exon_list[2*i+1] - self.exon_list[2*i])
                if self.exon_list[2*i] <= self.cdsStart_site <= self.exon_list[2*i+1]:
                    cds_start_relative = self.cdsStart_site - self.exon_list[2*i] + sum(exon_len[0:i])
                if self.exon_list[2*i] <= self.cdsEnd_site <= self.exon_list[2*i+1]:
                    cds_end_relative = self.cdsEnd_site - self.exon_list[2*i] + sum(exon_len[0:i])
                if self.exon_list[2*i] <= pos <= self.exon_list[2*i+1]:
                    pos_relative = pos - self.exon_list[2*i] + sum(exon_len[0:i])
            return [0, cds_start_relative, cds_end_relative, sum(exon_len)], pos_relative
        if self.gene_direction == '-':
            exon_len = []
            for i in range(0, int(len(self.exon_list)/2)):
                exon_len.append(self.exon_list[2*i] - self.exon_list[2*i+1])
                if self.exon_list[2*i+1] <= self.cdsStart_site <= self.exon_list[2*i]:
                    cds_start_relative = self.exon_list[2*i] - self.cdsStart_site + sum(exon_len[0:i])
                if self.exon_list[2*i+1] <= self.cdsEnd_site <= self.exon_list[2*i]:
                    cds_end_relative = self.exon_list[2*i] - self.cdsEnd_site + sum(exon_len[0:i])
                if self.exon_list[2*i+1] <= pos <= self.exon_list[2*i]:
                    pos_relative = self.exon_list[2*i] - pos + sum(exon_len[0:i])
            return [0, cds_start_relative, cds_end_relative, sum(exon_len)], pos_relative
### to make metagene, by calculating metagene score(0---100---200---300)
    def mRNA_metagene_score(self, pos):
        mRNA_pos_list, judge = self.mRNA_transcript(pos)
        if mRNA_pos_list[1] == mRNA_pos_list[2]:
            return 'NA'
        else:
            if judge != -1:
                if 0 < judge <= mRNA_pos_list[1]:
                    return round(judge/mRNA_pos_list[1]*100)
                if mRNA_pos_list[1] <= judge <= mRNA_pos_list[2]:
                    return round((judge - mRNA_pos_list[1])/(mRNA_pos_list[2]-mRNA_pos_list[1])*100 + 100)
                if mRNA_pos_list[2] <= judge <= mRNA_pos_list[3]:
                    return round((judge - mRNA_pos_list[2])/(mRNA_pos_list[3]-mRNA_pos_list[2])*100 + 200)
            else:
                return 'NA'

    def call_relation_pos_in_mRNA(self, pos):
        mRNA_pos_list, judge = self.mRNA_transcript(pos)
        if mRNA_pos_list[1] == mRNA_pos_list[2]:
            return 'NA'
        else:
            a = '0'
            b = format(mRNA_pos_list[1]/mRNA_pos_list[3], '.4f')
            c = format(mRNA_pos_list[2] / mRNA_pos_list[3], '.4f')
            d = '1'
            judge_ratio = format(judge / mRNA_pos_list[3], '.4f')
            if judge != -1:
                if 0 <= judge <= mRNA_pos_list[1]:
                    return a + '----{' + judge_ratio + '}----' + b + '----' + c + '----' + d
                if mRNA_pos_list[1] <= judge <= mRNA_pos_list[2]:
                    return a + '----' + b + '----{' + judge_ratio + '}----' + c + '----' + d
                if mRNA_pos_list[2] <= judge <= mRNA_pos_list[3]:
                    return a + '----' + b + '----' + c + '----{' + judge_ratio + '}----' + d
            else:
                return 'NA'


    def call_pos_info(self, pos):  # determine which exon or intron or utr the given site is located in
        if self.gene_direction == '+':
            if self.cdsEnd_site == self.cdsStart_site:
                for i in range(0, len(self.exon_list)-1):
                    if self.exon_list[i + 1] >= pos >= self.exon_list[i]:
                        if i % 2 == 0:
                            return format(i / 2 + 1, '.0f') + "_exon" + "_of_total_" + format(
                                len(self.exon_list) / 2, '.0f') + "_exon"
                        elif i % 2 == 1:
                            return format((i + 1) / 2, '.0f') + "_intron" + "_of_total_" + format(
                                len(self.exon_list) / 2 - 1, '.0f') + "_intron"
                        else:
                            return "NA"
            elif self.exon_list[0] <= pos <= self.exon_list[-1]:
                for i in range(0, len(self.exon_list)):
                    if self.exon_list[i + 1] >= pos >= self.exon_list[i]:
                        if i % 2 == 0:
                            if self.txStart_site <= pos <= self.cdsStart_site:
                                return "5'-UTR"
                            elif self.txEnd_site >= pos >= self.cdsEnd_site:
                                return "3'-UTR"
                            return format(i / 2 + 1, '.0f') + "_exon" + "_of_total_" + format(len(self.exon_list) / 2, '.0f') + "_exon"
                        elif i % 2 == 1:
                            return format((i + 1) / 2, '.0f') + "_intron" + "_of_total_" + format(
                                len(self.exon_list) / 2 - 1, '.0f') + "_intron"
            else:
                return "NA"
        elif self.gene_direction == '-':
            if self.cdsStart_site == self.cdsEnd_site:
                for i in range(0, len(self.exon_list)):
                    if self.exon_list[i + 1] <= pos <= self.exon_list[i]:
                        if i % 2 == 0:
                            return format(i / 2 + 1, '.0f') + "_exon" + "_of_total_" + format(
                                len(self.exon_list) / 2, '.0f') + "_exon"
                        elif i % 2 == 1:
                            return format((i + 1) / 2, '.0f') + "_intron" + "_of_total_" + format(
                                len(self.exon_list) / 2 - 1, '.0f') + "_intron"
                        else:
                            return "NA"

            elif self.exon_list[-1] <= pos <= self.exon_list[0]:
                for i in range(0, len(self.exon_list)):
                    if self.exon_list[i + 1] <= pos <= self.exon_list[i]:
                        if i % 2 == 0:
                            if self.txStart_site >= pos >= self.cdsStart_site:
                                return "5'-UTR"
                            elif self.txEnd_site <= pos <= self.cdsEnd_site:
                                return "3'-UTR"
                            return format(i / 2 + 1, '.0f') + "_exon" + "_of_total_" + format(len(self.exon_list) / 2,
                                                                                              '.0f') + "_exon"
                        elif i % 2 == 1:
                            return format((i + 1) / 2, '.0f') + "_intron" + "_of_total_" + format(
                                len(self.exon_list) / 2 - 1, '.0f') + "_intron"
            else:
                return "NA"

    def call_relative_pos(self, pos):
        if self.gene_direction == '+':
            if self.txEnd_site == self.cdsStart_site:
                a = '0'
                b = '1'
                c = format((pos - self.txStart_site) / (self.txEnd_site - self.txStart_site), '.4f')
                return a + "---{" + c + "}---" + b
            gene_struc = [self.txStart_site, self.cdsStart_site, self.cdsEnd_site, self.txEnd_site]
            a = '0'
            b = format((gene_struc[1] - gene_struc[0]) / (gene_struc[3] - gene_struc[0]), '.4f')
            c = format((gene_struc[2] - gene_struc[0]) / (gene_struc[3] - gene_struc[0]), '.4f')
            d = '1'
            e = format((pos - gene_struc[0]) / (gene_struc[3] - gene_struc[0]), '.4f')
            if gene_struc[0] <= pos <= gene_struc[1]:
                return a + "---{" + e + "}---" + b + "---" + c + "---" + d
            elif gene_struc[1] <= pos <= gene_struc[2]:
                return a + "---" + b + "---{" + e + "}---" + c + "---" + d
            elif gene_struc[2] <= pos <= gene_struc[3]:
                return a + "---" + b + "---" + c + "---{" + e + "}---" + d
            else:
                return "NA"
        if self.gene_direction == '-':
            if self.txStart_site == self.cdsEnd_site:
                a = '0'
                b = '1'
                c = format((self.txStart_site - pos) / (self.txStart_site - self.txEnd_site), '.4f')
                return a + "---{" + c + "}---" + b
            gene_struc = [self.txStart_site, self.cdsStart_site, self.cdsEnd_site, self.txEnd_site]
            a = '0'
            b = format((gene_struc[0] - gene_struc[1]) / (gene_struc[0] - gene_struc[3]), '.4f')
            c = format((gene_struc[0] - gene_struc[2]) / (gene_struc[0] - gene_struc[3]), '.4f')
            d = '1'
            e = format((gene_struc[0] - pos) / (gene_struc[0] - gene_struc[3]), '.4f')
            if gene_struc[1] <= pos <= gene_struc[0]:
                return a + "---{" + e + "}---" + b + "---" + c + "---" + d
            elif gene_struc[2] <= pos <= gene_struc[1]:
                return a + "---" + b + "---{" + e + "}---" + c + "---" + d
            elif gene_struc[3] <= pos <= gene_struc[2]:
                return a + "---" + b + "---" + c + "---{" + e + "}---" + d
            else:
                return "NA"

gene_annotation_dict_f = {}
gene_annotation_dict_r = {}

def transeq(a):
    trans = ''
    for i in a:
        if i == 'A':
            trans = trans + 'T'
        if i == 'T':
            trans = trans + 'A'
        if i == 'C':
            trans = trans + 'G'
        if i == 'G':
            trans = trans + 'C'
    return trans


with open('/home/yutaozhao/Genome/hg38-UCSC-ensembl', 'r') as f:
    for line in f:
        gene_info = list(item for item in line.split('\t'))
        refseq_id = gene_info[0]
        Chr = gene_info[1]
        gene_direction = gene_info[2]
        if gene_direction == "+":
            txStart_site = gene_info[3]
            txEnd_site = gene_info[4]
            cdsStart_site = gene_info[5]
            cdsEnd_site = gene_info[6]
            exon_num = int(gene_info[7])
            exStart = [int(i) if i != '' else 0 for i in gene_info[8].split(',')][0:-1]
            exEnd = [int(i) if i != '' else 0 for i in gene_info[9].split(',')][0:-1]
            exon_list = []

            i = 0
            while (exon_num > 0):
                exon_list.append(exStart[i])
                exon_list.append(exEnd[i])
                i = i + 1
                exon_num = exon_num - 1
            gene_name = gene_info[11].strip()
            gene_class = gene_info[12].strip()
            refseq_gene = gene_annotation(refseq_id, Chr, gene_direction, txStart_site, txEnd_site, cdsStart_site,
                                          cdsEnd_site, exon_list, gene_name, gene_class)

            if refseq_gene.Chr not in gene_annotation_dict_f.keys():
                gene_annotation_dict_f[refseq_gene.Chr] = [
                    [refseq_gene, refseq_gene.txStart_site, refseq_gene.txEnd_site]]
            else:
                gene_annotation_dict_f[refseq_gene.Chr].append(
                    [refseq_gene, refseq_gene.txStart_site, refseq_gene.txEnd_site])

        elif gene_direction == "-":
            txStart_site = gene_info[4]
            txEnd_site = gene_info[3]
            cdsStart_site = gene_info[6]
            cdsEnd_site = gene_info[5]
            exon_num = int(gene_info[7])
            exStart = [int(i) if i != '' else 0 for i in gene_info[8].split(',')][0:-1]
            exEnd = [int(i) if i != '' else 0 for i in gene_info[9].split(',')][0:-1]
            exon_list = []
            while (exon_num > 0):
                exon_list.append(exEnd[exon_num - 1])
                exon_list.append(exStart[exon_num - 1])
                exon_num = exon_num - 1
            gene_name = gene_info[11].strip()
            gene_class = gene_info[12].strip()
            refseq_gene = gene_annotation(refseq_id, Chr, gene_direction, txStart_site, txEnd_site, cdsStart_site,
                                          cdsEnd_site, exon_list, gene_name,gene_class)

            if refseq_gene.Chr not in gene_annotation_dict_r.keys():
                gene_annotation_dict_r[refseq_gene.Chr] = [
                    [refseq_gene, refseq_gene.txStart_site, refseq_gene.txEnd_site]]
            else:
                gene_annotation_dict_r[refseq_gene.Chr].append(
                    [refseq_gene, refseq_gene.txStart_site, refseq_gene.txEnd_site])

with open(sys.argv[1], 'w') as f:
    0
for line in sys.stdin:
    sites_info = line.split('\t')
    chr_info = sites_info[0]
    pos = int(sites_info[2])
    reads_distri = sites_info[3]
    del_ratio = float(sites_info[4])
    direction = sites_info[5].strip()
    save_unidentified = True
    if direction == "+" and chr_info in gene_annotation_dict_f.keys():
        for gene in gene_annotation_dict_f[chr_info]:
            if gene[2] >= pos >= gene[1]:
                motif = genome_seq_dict[chr_info].seq[pos - 3:pos + 2].upper().replace('T', 'U')
                if motif[1] == 'U':
                    motif1 = genome_seq_dict[chr_info].seq[pos - 4:pos + 1].upper().replace('T', 'U')
                    motif = motif + '|' + motif1
                    if motif[0] == 'U':
                        motif2 = genome_seq_dict[chr_info].seq[pos - 5:pos + 0].upper().replace('T', 'U')
                        motif = motif + '|' + motif2
                if motif[3] == 'U':
                    motif3 = genome_seq_dict[chr_info].seq[pos - 2:pos + 3].upper().replace('T', 'U')
                    motif = motif + '|' + motif3
                    if motif[4] == 'U':
                        motif4 = genome_seq_dict[chr_info].seq[pos - 1:pos + 4].upper().replace('T', 'U')
                        motif = motif + '|' + motif4
                print(
                    f'{chr_info}\t{pos}\t{direction}\t{gene[0].refseq_id}\t{gene[0].gene_name}\t{gene[0].gene_class}\t{gene[0].call_pos_info(pos)}\t{gene[0].call_relation_pos_in_mRNA(pos)}\t{gene[0].mRNA_metagene_score(pos)}\t{motif}\t{del_ratio}\t{reads_distri}')

                save_unidentified = False
                break
        if save_unidentified == False:
            continue
    if direction == "-" and chr_info in gene_annotation_dict_r.keys():
        for gene in gene_annotation_dict_r[chr_info]:
            if gene[1] >= pos >= gene[2]:
                motif = genome_seq_dict[chr_info].seq[pos - 3:pos + 2].upper()
                motif = transeq(motif)[::-1].replace('T', 'U')
                if motif[1] == 'U':
                    motif1 = genome_seq_dict[chr_info].seq[pos - 2:pos + 3].upper()
                    motif1 = transeq(motif1)[::-1].replace('T', 'U')
                    motif = motif + '|' + motif1
                    if motif[0] == 'U':
                        motif2 = genome_seq_dict[chr_info].seq[pos - 1:pos + 4].upper()
                        motif2 = transeq(motif2)[::-1].replace('T', 'U')
                        motif = motif + '|' + motif2
                if motif[3] == 'U':
                    motif3 = genome_seq_dict[chr_info].seq[pos - 4:pos + 1].upper()
                    motif3 = transeq(motif3)[::-1].replace('T', 'U')
                    motif = motif + '|' + motif3
                    if motif[4] == 'U':
                        motif4 = genome_seq_dict[chr_info].seq[pos - 5:pos + 0].upper()
                        motif4 = transeq(motif4)[::-1].replace('T', 'U')
                        motif = motif + '|' + motif4
                    #a, b = gene[0].mRNA_transcript(pos)
                print(
                    f'{chr_info}\t{pos}\t{direction}\t{gene[0].refseq_id}\t{gene[0].gene_name}\t{gene[0].gene_class}\t{gene[0].call_pos_info(pos)}\t{gene[0].call_relation_pos_in_mRNA(pos)}\t{gene[0].mRNA_metagene_score(pos)}\t{motif}\t{del_ratio}\t{reads_distri}')
                save_unidentified = False
                break
        if save_unidentified == False:
            continue
    if save_unidentified:
        with open(sys.argv[1], 'a') as f:
            f.write(line)












