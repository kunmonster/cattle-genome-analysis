import cyvcf2

input_vcf = 'input.vcf.gz'
output_vcf = 'output.vcf'

vcf = cyvcf2.VCF(input_vcf)

output = open(output_vcf, 'w')

output.write(str(vcf.raw_header))

for record in vcf:
    if  record.aaf is not None:
        af = record.aaf
        maf = min(af, 1 - af)
        if maf > 0.05:
            output.write(str(record))

output.close()
