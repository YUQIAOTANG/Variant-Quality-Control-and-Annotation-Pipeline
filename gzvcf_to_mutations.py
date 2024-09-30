'''
Author: Yuqiao Tang (tang.yuqi@northeastern.edu
This program prepares the input mutations file expected by the RVTT program from a filtered gzipped vcf file.
Input:
- filtered gzipped vcf file: The vcf file must be annotated with SnpEff, SnpSIFT and dbNSFP v4.0a database.
- tab-separated phenotype file in .fam format
- minimum minor allele frequency cutoff (suggested value = 0)
- maximum minor allele frequency cutoff (suggested value = 0.05)
- prefix of the output file
Output: 
- outprefix_mutations.tsv file that contains the necessary columns for running RVTT
'''

import sys
import gzip
import pandas as pd

# Constants for variant classifications
DAMAGING = set([
	'frameshift_variant', 'start_lost', 'stop_gained', 'stop_lost',
	'splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant',
	'structural_interaction_variant', 'initiator_codon_variant',
	'disruptive_inframe_deletion', 'disruptive_inframe_insertion'
])
EXCLUDE = set([
	'3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant',
	'5_prime_UTR_variant', 'bidirectional_gene_fusion', 'downstream_gene_variant',
	'gene_fusion', 'intergenic_region', 'intron_variant',
	'non_coding_transcript_exon_variant', 'sequence_feature',
	'TF_binding_site_variant', 'upstream_gene_variant'
])
INCLUDE = {'missense_variant', 'synonymous_variant'} | DAMAGING


def parse_allele_frequencies(info_parts):
	af_parts = {}
	keys = ['gnomAD_genomes_POPMAX_AF', 'gnomAD_exomes_POPMAX_AF',
			'gnomAD_genomes_NFE_AF', 'gnomAD_exomes_NFE_AF']
	for part in info_parts:
		key_value = part.split('=')
		if key_value[0] in keys:
			af_parts[key_value[0]] = float(key_value[1].split(',')[0]) if key_value[1] != '.' else 0.0
	return af_parts


def process_variant_line(line, indiv_start, controls, cases):
	parts = line.split('\t')
	chromosome, pos, rsid, ref, alt, info_str = parts[0:6]
	genotypes = parts[indiv_start:]

	info_parts = info_str.split(';')
	allele_frequencies = parse_allele_frequencies(info_parts)
	annotations = next((item for item in info_parts if item.startswith('ANN=')), None)

	if annotations:
		annotations = annotations.split('=')[1].split(',')
		primary_annotation = annotations[0].split('|')
		gene, impact, mutation_type = primary_annotation[3], primary_annotation[2], primary_annotation[1]

		if mutation_type not in INCLUDE:
			return None

		ac = genotypes.count('1/1') * 2 + genotypes.count('0/1')
		an = 2 * (len(genotypes) - genotypes.count('./.'))
		af = ac / an if an > 0 else 0

		ac_case = sum(1 for idx, gt in enumerate(genotypes) if gt in {'0/1', '1/1'} and indivs[idx] in cases)
		ac_control = sum(1 for idx, gt in enumerate(genotypes) if gt in {'0/1', '1/1'} and indivs[idx] in controls)

		mutated_individuals = [indivs[idx] for idx, gt in enumerate(genotypes) if gt in {'0/1', '1/1'}]

		return [gene, chromosome + ':' + pos, ref + '>' + alt, rsid, mutation_type, impact, af, ac_case, ac_control,
				';'.join(mutated_individuals)]


def main():
	if len(sys.argv) < 6:
		print("Usage: python preprocess_gzvcf.py filtered_input.vcf.gz famfile mincutoff maxcutoff outprefix")
		sys.exit(1)

	infile, famfile, mincutoff, maxcutoff, outprefix = sys.argv[1:6]
	mincutoff, maxcutoff = float(mincutoff), float(maxcutoff)

	fam_df = pd.read_csv(famfile, sep='\t', header=None)
	controls = fam_df[fam_df[5] == '1'][1].tolist()
	cases = fam_df[fam_df[5] == '2'][1].tolist()

	header = "#Gene\tMutID\tChange\tRSID\tMutationType\tImpact\tAF\tAC_Case\tAC_Control\tMutatedIndividuals\n"
	with gzip.open(infile, 'rt') as fp, open(outprefix + "_mutations.tsv", 'w') as of:
		of.write(header)
		for line in fp:
			if line.startswith('#'):
				if line.startswith('#CHROM'):
					indivs = line.strip().split('\t')[9:]
				continue
			result = process_variant_line(line, 9, controls, cases)
			if result:
				of.write('\t'.join(map(str, result)) + '\n')


if __name__ == '__main__':
	main()
