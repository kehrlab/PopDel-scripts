from sys import argv

if len(argv) != 2:
    print('Usage: ' + argv[0] + ' <VCF file>')
    exit(1)

vcffilename = argv[1]

def hwe_chi_squared(gt0, gt1, gt2):
    total = gt0 + gt1 + gt2

    pfreq = (gt0 + 0.5 * gt1) / total
    qfreq = (gt2 + 0.5 * gt1) / total

    exp_gt0 = max(pfreq * pfreq * total, 0.0000000000001)
    exp_gt1 = max(2 * pfreq * qfreq * total, 0.0000000000001)
    exp_gt2 = max(qfreq * qfreq * total, 0.0000000000001)

    chi0 = ((gt0 - exp_gt0) * (gt0 - exp_gt0)) / exp_gt0
    chi1 = ((gt1 - exp_gt1) * (gt1 - exp_gt1)) / exp_gt1
    chi2 = ((gt2 - exp_gt2) * (gt2 - exp_gt2)) / exp_gt2

    return chi0 + chi1 + chi2

with open(vcffilename, 'r') as vcffile:
    for line in vcffile:

        if line[0] == '#':
            print(line.strip())
            continue

        fields = line.split()
        info = fields[7].split(';')

        gt0, gt1, gt2 = [0, 0, 0]
        for sample in fields[10:]:
            if sample[0:3] == '0/0':
                gt0 += 1
            if sample[0:3] == '0/1':
                gt1 += 1
            if sample[0:3] == '1/1':
                gt2 += 1

        chi_squared = hwe_chi_squared(gt0, gt1, gt2)
        if chi_squared < 6.635:
            print('\t'.join(fields[0:7] + [fields[7] + ';HWE=' + str(chi_squared)] + fields[8:]))
