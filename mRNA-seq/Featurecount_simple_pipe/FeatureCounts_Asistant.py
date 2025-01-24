import sys
def print_list(a):
    b = ''
    for i in a:
        b = b + str(i) + "\t"
    b = b.strip() + "\n"
    return(b)
header_1 = ''
for i in range(int(sys.argv[4])):
    header_1 = header_1 + sys.argv[3] + "_" + str(i+1) + "\t"
for i in range(int(sys.argv[6])):
    header_1 = header_1 + sys.argv[5] + "_" + str(i+1) + "\t"
header_1 = header_1.strip()
header_1 = "Gene_name" + "\t" + header_1 + "\n"
#print(header_1)
### read read_number in
read_num = []
with open(sys.argv[2],'r') as f:
    for line in f:
        read_num.append(int(line.strip()))
#print(read_num)
### crack Counts file
NR_limit = 0
file_counts = open(f'{sys.argv[3]}_vs_{sys.argv[5]}.Counts.txt', 'w')
file_RPKM = open(f'{sys.argv[3]}_vs_{sys.argv[5]}.RPKM.txt', 'w')
file_counts.write(header_1)
file_RPKM.write(header_1)
with open(sys.argv[1],'r') as f:
    for line in f:
        if NR_limit < 2:
            info = line.split()
            length = len(info)
            NR_limit = NR_limit + 1
            continue
        info = line.split()
        Counts_list = [info[0]]
        RPKM_list = [info[0]]
        for i in range(length - 6):
            Counts_list.append(info[i+6].strip())
            #print(info)

            RPKM_list.append(round(int(info[i+6].strip())/int(info[5])/int(read_num[i])*1000000000,4))
            #print(RPKM_list)
        file_counts.write(print_list(Counts_list))
        file_RPKM.write(print_list(RPKM_list))
file_counts.close()
file_RPKM.close()




