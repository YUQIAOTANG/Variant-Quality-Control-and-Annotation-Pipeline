import sys
import gzip

def parse_info_field(info_str):
    """
    Parses the INFO field of a VCF entry and returns a dictionary of key-value pairs.
    """
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def main():
    if len(sys.argv) < 3:
        print("Usage: ./hard_filter.py in.vcf out.vcf")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    with open(infile, 'r') as inf, open(outfile, 'w') as outf:
        for line in inf:
            line = line.strip()
            if line.startswith('#'):
                outf.write(line + '\n')
                continue

            fields = line.split('\t')
            ref = fields[3]
            alt = fields[4]
            info_field = fields[7]

            # Determine variant type: 0 for SNP, 1 for INDEL
            variant_type = 1 if len(ref) > 1 or len(alt) > 1 else 0

            # Parse the INFO field into a dictionary
            info = parse_info_field(info_field)

            # Extract required metrics with default values if not present
            fs = float(info.get('FS', '201'))
            inbrd = float(info.get('InbreedingCoeff', '-0.8'))
            mq = float(info.get('MQ', '40.0'))
            mqrs = float(info.get('MQRankSum', '-12.5'))
            qd = float(info.get('QD', '2.0'))
            rprs = float(info.get('ReadPosRankSum', '-0.8'))
            sor = float(info.get('SOR', '11'))

            print(fs, inbrd, mq, mqrs, qd, rprs, sor)

            # Apply filtering criteria based on variant type
            if variant_type == 0:
                if (fs < 60.0 and inbrd > -0.8 and mq > 40 and mqrs > -12.5 and
                    qd > 2.0 and rprs > -8.0 and sor <= 3.0):
                    outf.write(line + '\n')
                else:
                    print(variant_type, fs < 60, inbrd > -0.8, mq > 40, mqrs > -12.5,
                          qd > 2, rprs > -0.8, sor <= 3.0)
            else:
                if fs < 200 and qd > 2.0 and rprs > -20 and sor <= 10:
                    outf.write(line + '\n')
                else:
                    print(variant_type, fs < 200, qd > 2, rprs > -20, sor <= 10)

if __name__ == "__main__":
    main()
# This script reads a VCF file and applies hard filtering based on variant quality metrics.