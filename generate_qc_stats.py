import sys
from collections import Counter


def is_in_cpg(cpg_locs, mutation_id):
    chromosome, position = mutation_id.split(':')
    position = int(position)
    for start, end in cpg_locs.get(chromosome, []):
        if start <= position <= end:
            return True
    return False


def parse_file(filepath):
    with open(filepath, 'r') as file:
        return [line.strip() for line in file]


def process_lines(names, lines, cpg_locs):
    inmap = {name: [0] * 28 for name in names}
    transition_muts = {"A>G", "G>A", "C>T", "T>C"}

    for line in lines[1:]:
        parts = line.split('\t')
        mutation_id = parts[1]
        mutation_type = parts[2]
        change = parts[3]
        dbsnp_status = parts[4]
        effect = parts[6]
        mutation_indices = [v for v in parts[-1].split(';') if v in names]
        counts = Counter(mutation_indices)

        for name in names:
            current_counts = counts.get(name, 0)
            inmap[name][3] += 1 if current_counts == 0 else 0
            if current_counts:
                base_index = 26 if mutation_type == "INDEL" else 5
                inmap[name][base_index] += 1
                inmap[name][2] += current_counts

                if change in transition_muts:
                    inmap[name][0] += 1 if current_counts == 1 else 0
                    inmap[name][1] += 0 if current_counts == 1 else 0

                if dbsnp_status == '.':
                    inmap[name][6] += 1
                    cpg_influence = 11 if is_in_cpg(cpg_locs, mutation_id) else 18
                    inmap[name][cpg_influence] += 1

                else:
                    inmap[name][7] += 1
                    if is_in_cpg(cpg_locs, mutation_id):
                        inmap[name][10] += 1

        # Process effects
        effect_mapping = {
            'synonymous': 20,
            'missense': 22,
            'stop_gained': 24
        }
        base_index = effect_mapping.get(effect, None)
        if base_index:
            inmap[name][base_index] += 1
            inmap[name][base_index + 1] += 0 if change in transition_muts else 1

    return inmap


def main():
    if len(sys.argv) < 5:
        print("Usage: python stat.py input_mutations_file cpgfile famfile outfile")
        sys.exit(1)

    infile, cpgfile, famfile, outfile = sys.argv[1:5]
    cpg_locs = {}

    for line in parse_file(cpgfile):
        if line[0] != "#":
            chromosome, start, end = line.split('\t')[:3]
            cpg_locs.setdefault(chromosome, []).append((int(start), int(end)))

    names = [line.split(maxsplit=1)[1] for line in parse_file(famfile)]
    mutation_lines = parse_file(infile)
    inmap = process_lines(names, mutation_lines, cpg_locs)

    with open(outfile, "w") as f:
        header = "\t".join(["ID", "R/A_Ts", "R/A_Tv", "Indel", "Ref", "Het", "Hom",
                            "novelSNP", "knownSNP", "A/D_Ts", "A/D_Tv", "knownCpG", "novelCpG",
                            "knownTs", "knownTv", "novelTs", "novelTv", "nCpG-K_Ts", "nCpG-K_Tv",
                            "nCpG-N_Ts", "nCpG-N_Tv", "synonymousTs", "synonymousTv", "missenseTs",
                            "missenseTv", "nonsenseTs", "nonsenseTv", "homINDEL", "hetINDEL"])
        f.write(header + "\n")
        for name, values in inmap.items():
            f.write(f"{name}\t" + "\t".join(map(str, values)) + "\n")


if __name__ == "__main__":
    main()
