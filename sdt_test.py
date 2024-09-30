import sys
from scipy.stats import binomtest as btest


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 generate_sdt_counts.py input_mutations_file famfile outputfile")
        exit(1)

    infile = sys.argv[1]
    famfile = sys.argv[2]
    outfile = sys.argv[3]

    # 使用with语句读取famfile，确保文件正确关闭
    cases = set()
    controls = set()
    with open(famfile) as f:
        for line in f:
            x = line.strip().split()
            if x[5] == '2':
                cases.add(x[1])
            elif x[5] == '1':
                controls.add(x[1])

    Ncase = len(cases)
    Ncontrol = len(controls)
    Ntotal = Ncase + Ncontrol

    singletonP = Ncase / Ntotal
    doubletonP = 2 * Ncase / Ntotal * Ncontrol / Ntotal
    tripletonP = 1 - (Ncase / Ntotal) ** 3 - (Ncontrol / Ntotal) ** 3

    # 初始化计数器，索引0对应同义变异，索引1对应错义变异
    scase = [0, 0]
    scon = [0, 0]
    dcase = [0, 0]
    dcon = [0, 0]
    dsh = [0, 0]
    tcase = [0, 0]
    tcon = [0, 0]
    tsh = [0, 0]

    with open(infile) as f:
        next(f)  # 跳过表头
        for line in f:
            x = line.strip().split('\t')
            if x[2] == "INDEL":
                continue
            csq = x[6]
            if csq not in {'synonymous_variant', 'missense_variant'}:
                continue
            y = x[23].split(';')
            y_set = set(y)
            ncase = len(y_set & cases)
            ncontrol = len(y_set & controls)
            ll = ncase + ncontrol
            if ll > 3:
                continue
            idx = 0 if csq == 'synonymous_variant' else 1
            if ll == 1:
                if ncase == 1:
                    scase[idx] += 1
                elif ncontrol == 1:
                    scon[idx] += 1
            elif ll == 2:
                if ncase == 2:
                    dcase[idx] += 1
                elif ncontrol == 2:
                    dcon[idx] += 1
                elif ncase == 1 and ncontrol == 1:
                    dsh[idx] += 1
            elif ll == 3:
                if ncase == 3:
                    tcase[idx] += 1
                elif ncontrol == 3:
                    tcon[idx] += 1
                elif (ncase == 2 and ncontrol == 1) or (ncase == 1 and ncontrol == 2):
                    tsh[idx] += 1

    # 计算p值
    syn_singleton_pval = btest(scase[0], scase[0] + scon[0], singletonP, alternative="greater").pvalue
    syn_doubleton_pval = btest(dsh[0], dcase[0] + dcon[0] + dsh[0], doubletonP, alternative="less").pvalue
    syn_tripleton_pval = btest(tsh[0], tcase[0] + tcon[0] + tsh[0], tripletonP, alternative="less").pvalue

    mis_singleton_pval = btest(scase[1], scase[1] + scon[1], singletonP, alternative="greater").pvalue
    mis_doubleton_pval = btest(dsh[1], dcase[1] + dcon[1] + dsh[1], doubletonP, alternative="less").pvalue
    mis_tripleton_pval = btest(tsh[1], tcase[1] + tcon[1] + tsh[1], tripletonP, alternative="less").pvalue

    with open(outfile, "w") as of:
        of.write("Synonymous\n")
        of.write("Singleton\t{}\t{}\t\t{}\t{}\n".format(scase[0], scon[0], singletonP, syn_singleton_pval))
        of.write("Doubleton\t{}\t{}\t{}\t{}\t{}\n".format(dcase[0], dcon[0], dsh[0], doubletonP, syn_doubleton_pval))
        of.write("Tripleton\t{}\t{}\t{}\t{}\t{}\n".format(tcase[0], tcon[0], tsh[0], tripletonP, syn_tripleton_pval))

        of.write("Missense\n")
        of.write("Singleton\t{}\t{}\t\t{}\t{}\n".format(scase[1], scon[1], singletonP, mis_singleton_pval))
        of.write("Doubleton\t{}\t{}\t{}\t{}\t{}\n".format(dcase[1], dcon[1], dsh[1], doubletonP, mis_doubleton_pval))
        of.write("Tripleton\t{}\t{}\t{}\t{}\t{}\n".format(tcase[1], tcon[1], tsh[1], tripletonP, mis_tripleton_pval))


if __name__ == "__main__":
    main()
