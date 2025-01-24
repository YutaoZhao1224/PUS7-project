import sys

for line in sys.stdin:
    info = line.split()
    counts = int(info[3])
    Forward_seq_mutation = [0, 0, 0, 0, 0, 0, 0] # match, skip, mismatch_A, mismatch_G, mismatch_C, mismatch_T, INDEL
    Reverse_seq_mutation = [0, 0, 0, 0, 0, 0, 0]
    start = 0
    stop = 0
    number_str = 0
    jump_Del = 0
    jump_Ins = 0
    if counts >= 20:
        for i in info[4]:
            number_str += 1

            if jump_Ins != 0:
                jump_Ins -= 1
                continue
            if jump_Del != 0:
                jump_Del -= 1
                continue

            if i == '.':
                Forward_seq_mutation[0] += 1
                continue
            if i == ',':
                Reverse_seq_mutation[0] += 1
                continue
            if i == '>':
                Forward_seq_mutation[1] += 1
                continue
            if i == '<':
                Reverse_seq_mutation[1] += 1
                continue
            if i == 'A':
                Forward_seq_mutation[2] += 1
                continue
            if i == 'a':
                Reverse_seq_mutation[2] += 1
                continue
            if i == 'G':
                Forward_seq_mutation[3] += 1
                continue
            if i == 'g':
                Reverse_seq_mutation[3] += 1
                continue
            if i == 'C':
                Forward_seq_mutation[4] += 1
                continue
            if i == 'c':
                Reverse_seq_mutation[4] += 1
                continue
            if i == 'T':
                Forward_seq_mutation[5] += 1
                continue
            if i == 't':
                Reverse_seq_mutation[5] += 1
                continue
            if i == '-':
                jump_Del = int(info[4][number_str]) + 1
                if info[4][number_str + 1] == 'T':
                    Forward_seq_mutation[6] += 1
                if info[4][number_str + 1] == 'a':
                    Reverse_seq_mutation[6] += 1
                continue
            if i == '+':
                jump_Ins = int(info[4][number_str]) + 1
                continue
            if i == '^':
                jump_Ins = 1
                start += 1
                continue
            if i == '$':
                stop += 1
                continue

        total_counts_forward = Forward_seq_mutation[0] + Forward_seq_mutation[2] + Forward_seq_mutation[3] + \
                               Forward_seq_mutation[4] + Forward_seq_mutation[5]
        total_counts_Reverse = Reverse_seq_mutation[0] + Reverse_seq_mutation[2] + Reverse_seq_mutation[3] + \
                               Reverse_seq_mutation[4] + Reverse_seq_mutation[5]

        try:
            Deletion_counts_forward = Forward_seq_mutation[6]
            Deletion_ratio_forward = round(Deletion_counts_forward / total_counts_forward, 4)
            if total_counts_forward >= 20 and Deletion_counts_forward / total_counts_forward <= 0.02:
                print(f'{info[0]}\t{int(info[1]) + 1}\t+\t{Deletion_ratio_forward}\t{total_counts_forward}\t{Deletion_counts_forward}\tup1_mutation_A/G/C={Forward_seq_mutation[2]}:{Forward_seq_mutation[3]}:{Forward_seq_mutation[4]}\t{stop}')
        except:
            0
        try:
            Deletion_counts_reverse = Reverse_seq_mutation[6]
            Deletion_ratio_reverse = round(Deletion_counts_reverse / total_counts_Reverse, 4)
            if total_counts_Reverse >= 20 and Deletion_counts_reverse / total_counts_Reverse <= 0.02:
                print(f'{info[0]}\t{int(info[1]) + 1}\t-\t{Deletion_ratio_reverse}\t{total_counts_Reverse}\t{Deletion_counts_reverse}\tdown1_mutation_A/G/C={Reverse_seq_mutation[5]}:{Reverse_seq_mutation[4]}:{Reverse_seq_mutation[3]}\t{stop}')
        except:
            0