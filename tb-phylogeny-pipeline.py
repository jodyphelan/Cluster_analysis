import sys
import subprocess
import argparse
import random
import os
rand_generator = random.SystemRandom()

def get_random_file(prefix = None,extension=None):
	randint = rand_generator.randint(1,999999)
	if prefix:
		if extension:
			return "%s.%s.%s" % (prefix,randint,extension)
		else:
			return "%s.%s.txt" % (prefix,randint)
	else:
		if extension:
			return "%s.tmp.%s" % (randint,extension)
		else:
			return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
	sys.stderr.write("\n"+str(msg)+"\n")
	if ext:
		exit(1)

def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False


class vcf_class:
	def __init__(self,filename,threads=4):
		self.samples = []
		self.filename = filename
		self.threads = threads
		self.prefix = filename[:-7]
		if nofile(filename+".csi"):
			run_cmd("bcftools index  %(filename)s" % vars(self))
		self.temp_file = get_random_file()
		run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
		for l in open(self.temp_file):
			self.samples.append(l.rstrip())
		os.remove(self.temp_file)
	def vcf_to_fasta(self,ref_file,threads=4,chunk_size = 50000):
		self.ref_file = ref_file
		self.chunk_size = chunk_size
		self.cmd_split_chr = "bedtools makewindows -g %(ref_file)s.fai -w %(chunk_size)s -s  %(chunk_size)s | awk '{print $1\":\"$2\"-\"$3}'" % vars(self)
		self.tmp_file = "%s.tmp.txt" % self.prefix
		self.threads = threads
		cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"bcftools view  %(filename)s -r {} -Ou | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' | sed 's/\*[\/|]\*/\.\/\./g' |  datamash transpose > %(prefix)s.{}.tmp.txt\"" % vars(self)
		run_cmd(cmd)
		cmd = "paste `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'` > %(tmp_file)s" % vars(self)
		run_cmd(cmd)
		cmd = "rm `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'`" % vars(self)
		run_cmd(cmd)
		O = open(self.prefix+".snps.fa","w")
		for i,l in enumerate(open(self.tmp_file)):
			row = l.rstrip().split()
			if i==0: continue
			s = self.samples[i-1]
			seq = "".join(row).replace("./.","N")
			O.write(">%s\n%s\n" % ( s,seq))
		O.close()
	def vcf_to_matrix(self,):
		self.matrix_file = self.prefix+".mat"
		self.binary_matrix_file = self.prefix+".mat.bin"
		O = open(self.matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' >> %(matrix_file)s" % vars(self))
		O = open(self.binary_matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%GT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' | sed 's/1\/1/1/g' | sed 's/0\/0/0/g' >> %(binary_matrix_file)s" % vars(self))

def main(args):
	params = {"threads": args.threads, "prefix": args.prefix, "ref": args.ref}
	params["map_file"] = "%s.map" % (args.prefix)
	with open(params["map_file"],"w") as O:
		# Set up list to hold sample names
		samples = []
		# Loop through sample-file and do (1) append samples to list, (2) write sample to map file and (3) check for VCF index
		for line in open(args.sample_file):
			sample = line.rstrip()
			samples.append(sample)
			O.write("%s\t%s/%s%s\n" % (sample, args.vcf_dir, sample, args.vcf_extension))
			if nofile("%s/%s%s.tbi" % (args.vcf_dir, sample, args.vcf_extension)):
				run_cmd("bcftools index --tbi %s/%s%s" % (args.vcf_dir, sample, args.vcf_extension))

	# Create .dict file (GATK fasta index) has been created for the reference
	if nofile("%s.dict" % args.ref.replace(".fasta","").replace(".fa","")):
		run_cmd("gatk CreateSequenceDictionary -R %(ref)s" % params)
	# Create .fai file (SAMtools fasta index) has been created for the reference
	if nofile("%s.fai" % args.ref.replace(".fasta","").replace(".fa","")):
		run_cmd("samtools faidx %(ref)s" % params)
	run_cmd("gatk GenomicsDBImport --genomicsdb-workspace-path %(prefix)s_genomics_db -L Chromosome --sample-name-map %(map_file)s --reader-threads %(threads)s" % params)
	run_cmd("gatk --java-options \"-Xmx40g\" GenotypeGVCFs -R %(ref)s -V gendb://%(prefix)s_genomics_db -O %(prefix)s.raw.vcf.gz" % params)
	run_cmd("bcftools view -V indels %(prefix)s.raw.vcf.gz | bcftools filter -e 'GT=\"het\"' -S . | awk 'length($4)==1 || $0~/^#/' | tr '|' '/' | tr '*' '.' | bcftools view -a | bcftools view -c 1 -Oz -o %(prefix)s.filt.vcf.gz" % params)
	vcf = vcf_class("%s.filt.vcf.gz" % (args.prefix))
	vcf.vcf_to_fasta(args.ref)
	vcf.vcf_to_matrix()
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample-file',help='sample file',required=True)
parser.add_argument('--prefix',help='Prefix for files',required=True)
parser.add_argument('--ref',help='reference file',required=True)
parser.add_argument('--vcf-dir',default="./vcf/", type=str, help='VCF firectory')
parser.add_argument('--vcf-extension',default=".gatk.vcf.gz", type=str, help='VCF extension')
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
