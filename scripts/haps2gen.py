import sys,os,time,argparse,string,gzip,io

baselookup = dict(C='G',G='C',A='T',T='A')

def wopen(filename,buf_in_MB=100):
    if os.path.isfile(filename):   sys.exit(filename+" exists!")
    elif filename[-3:]==".gz": f = gzip.open(filename,"wb")
    else: f = open(filename,"wb")
    return  io.BufferedWriter(f,buf_in_MB*(2**20))


def ropen(filename,buf_in_MB=100):
    if  os.path.isfile(filename+".gz"): filename+=".gz"
    if not os.path.isfile(filename):  sys.exit(filename+" does not exist!")
    elif filename[-3:]==".gz":    f = gzip.open(filename,"rb")
    else: return( open(filename,"rb",buf_in_MB*(2**20)))
    return io.BufferedReader(f,buf_in_MB*(2**20))

parser = argparse.ArgumentParser(description='Converts shapeit2 haps to a .gen file')
parser.add_argument('input', metavar='input', type=str, help='shapeit2 output')
#parser.add_argument('sample', metavar='sample', type=str, help='sample file')
parser.add_argument('chrom', metavar='chromo',type=str, help='chromosome')
parser.add_argument('-output', action ='store',dest='output',metavar='output.vcf.gz', default = '',type=str, help='output file')

args = parser.parse_args()

if args.output=="":
    args.output = args.input #args.haps.split(".")[0]


samplefile = ropen(args.input + ".sample")
outfile = wopen(args.output+".gen.gz")
samplefile.next()
samplefile.next()
sampleids = []
for row in samplefile:
    sampleids.append(row.split()[1])
nsample = len(sampleids)
print nsample,'samples'

genotypes = ["1 0 0","0 1 0","0 0 1"]

infile = ropen(args.input + ".haps")
for linenum,row in enumerate(infile):
    data = row.split()
    gt = [genotypes[int(data[idx])+int(data[idx+1])] for idx in range(5,len(data),2)]
    outfile.write(string.join([data[1],data[1],data[2],data[3],data[4]]+gt," ")+"\n")
