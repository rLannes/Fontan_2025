import argparse
from copy import deepcopy
import logging
import numpy as np
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from multiprocessing import Process
from multiprocessing import Pool
from pathlib import Path
import pickle
import Rust_covpyo3
#from Rust_covpyo3 import get_coverage, get_header

# version april 12 2024


def get_cover_for_a_bam(intervalls, bam, lib, mapq, continous=True):
    """
    intervalls is a list of dict
    list of integer
    """
    coverage_rust  = []
    for e in intervalls:
       #start = time.time()
        cover = Rust_covpyo3.get_coverage(e["start"], e["end"], e["chr"], e["strand"], bam, lib, mapq)
        if continous:
            coverage_rust.extend(cover)
        else:
            coverage_rust.append(cover)
       # print("rust: {} ".format(time.time()-start), len(cover))

    return coverage_rust


def parse_bam_mp(args):
    """ 
        compute the stranded coverage on intervall for a list of bam files. in parrallel 
        [[str], dict] -> dict[str] = [int]
        take a a list of files and a dictionnary with the following fields start -> int; end -> int, chr -> str, strand -> str.
        return a dict[file] -> list of integer
        ### FUTURE add option to use custom binary binary flag for the reads filtering. 
        ### FUTURE support other stranded non stranded
    """

    (bam, lib,  intervalls, mapq, flag_in, flag_out) = args
    #print(bam, intervalls)
    #results =  Rust_covpyo3.get_coverage(intervalls["start"], intervalls["end"], intervalls["chr"], intervalls["strand"], bam, lib, mapq)
    results =  Rust_covpyo3.get_coverage_algo2(intervalls["start"], intervalls["end"], intervalls["chr"], intervalls["strand"],
                                                bam, lib, mapq, flag_in, flag_out)
    return results


class Coverage():
    def __init__(self, intervalls) -> None:
        self.intervalls = intervalls
        self.cover = None
        self.strand = None
        self.cover_intron_compress = None
        self.intervalls_scaled = None
        self.intervalls_scaled_intron_compress = None
        self.max_reads = 0
        self.N = None


    @staticmethod
    def merge_duplicates(list):

        cover = Coverage(deepcopy(list[0].intervalls)) 
        cover.strand = list[0].strand
        cover.sort_intervalls()
        cover.make_intervall_scaled()

        if not all([x for x in list[1:] if x.cover.shape == list[0].cover.shape]):
            raise AssertionError("cannot add because the shape does not match {} {}".format(self.cover.shape, other.shape))
        if not all([x for x in list[1:] if x.intervalls == list[0].intervalls]):
            raise AssertionError("cannot add because the intervalls are not the same {}\n{}".format(self.intervalls, other.intervalls))
        
        for e in list:
            if cover.cover is None:
                cover.cover = np.array(e.cover)
                cover.max_reads = e.max_reads
            else:
                cover.cover = cover.cover + np.array(e.cover)
                cover.max_reads += e.max_reads

        return cover
        



    def __add__(self, other):

        if self.cover.shape != other.cover.shape:
            raise AssertionError("cannot add because the shape does not match {} {}".format(self.cover.shape, other.shape))
        if self.intervalls != other.intervalls:
            raise AssertionError("cannot add because the intervalls are not the same {}\n{}".format(self.intervalls, other.intervalls))

        cover = Coverage(deepcopy(self.intervalls)) 
        cover.strand = self.strand
        cover.sort_intervalls()
        cover.make_intervall_scaled()
        cover.set_cover(self.cover + other.cover)
        return cover

    @property
    def get_normalize(self, factor=1_000_000):
        return self.cover / self.max_reads * factor
    
    def normalize(self, factor=1_000_000):
        self.cover = self.cover / self.max_reads * factor

    def add_intron(self):
        introns = []
        self.sort_intervalls()
        if len (self.intervalls) == 1 or "intron" in self.intervalls[1]:
            return
        
        for i, e in enumerate(self.intervalls[:-1]):
            intron = {}
            intron["start"] = e["end"]
            intron["end"] = self.intervalls[i + 1]["start"]  # e["end"]
            intron["strand"] = e["strand"]
            intron["chr"] = e["chr"]
            intron["type"] = "intron"
            introns.append(intron)
        self.intervalls.extend(introns)
        self.sort_intervalls()

    def sort_intervalls(self):
        self.intervalls = sorted(self.intervalls, key=lambda x: x["start"])

    def set_cover(self, cover):
        self.cover = cover

    def make_intervall_scaled(self):
        new_interval = []
        self.sort_intervalls()
        min_ = self.intervalls[0]["start"]
        for x in self.intervalls:
            item = {"start":  x["start"] - min_, "end": x["end"] - min_, "strand": x["strand"], "chr": x["chr"], "type": x["type"]}
            new_interval.append(item)
        self.intervalls_scaled = new_interval

    def compress_intron(self, max_ratio=0.45):

        sum_exon = sum([x["end"] - x["start"] for x in self.intervalls if x["type"] == "exon"])
        sum_intron = sum([x["end"] - x["start"] for x in self.intervalls if x["type"] == "intron"])
        current_r = sum_intron / (sum_exon + sum_intron)
        if current_r < max_ratio:
            self.cover_intron_compress = self.cover
            self.intervalls_scaled_intron_compress = self.intervalls_scaled
            return

        final_s_intron = sum_exon * max_ratio
        coefficient = final_s_intron / sum_intron

        new_cover = []
        new_interval = []
        start = 0
        for sub_i in self.intervalls_scaled:

            sub_i = deepcopy(sub_i)
            cov = self.cover[sub_i["start"]: sub_i["end"]]
            _len = len(cov)
            if sub_i["type"] == "intron":
                windows_s = max(1, int(_len / (1+(coefficient * _len))))

                cov = [np.mean(cov[i:i + windows_s]) for i in range(0, _len, windows_s)]
                _len = len(cov)

            new_cover.extend(cov)
            sub_i["start"] = start
            sub_i["end"] = start + _len
            new_interval.append(sub_i)
            start += _len

        self.intervalls_scaled_intron_compress = new_interval
        self.cover_intron_compress = new_cover
    
    def get_cover_exon(self):
        r = []
        for i, v in enumerate(self.intervalls_scaled):    
            if v["type"] == "exon":
                r.extend(self.cover[v["start"] : v["end"]])
        return r
    
    

def get_max_read(bam_file):

    #if bam_file.endswith("Aligned.sortedByCoord.out.bam"):
    #    log_file = bam_file.replace("Aligned.sortedByCoord.out.bam", "Log.final.out")
    #else:
    log_file = bam_file.split("Aligned")[0] + "Log.final.out"
        #log_file = bam_file.split(".")[0] + ".Log.final.out"
    #print(log_file)
    flagstat = Path(bam_file).parent / (str(Path(bam_file).stem) + ".flagstat")
    # print(log_file)
    # print(flagstat)
    if Path(log_file).exists():
        # print(log_file)
        with open(log_file) as f_in:
            for l in f_in:
                spt = [x.strip() for x in l.strip().split("|")]
                if spt[0] == "Uniquely mapped reads number":
                    return int(spt[1])
    elif Path(flagstat).exists():
        # print(flagstat)
        with open(flagstat) as f_in:
            for l in f_in:
                if "properly paired" in l:
                    return int(l.split()[0])
    else:
        raise AssertionError("no file found to retirve max reads, I tried: {} and {}".format(log_file, flagstat))


def merge_duplicate(dico_bam_cover):

    results = {}
    for bam, cover_reads in dico_bam_cover.items():

        m = reg_get_sample_lanes.match(Path(bam).stem)
        genotype = m.group(1)
        sample = m.group(2)
        lane = m.group(3)

        id_ = "{}_{}".format(genotype, sample)

        if id_ not in results:
            results[id_] = cover_reads
        else:
            cover = results[id_]
            cover.cover = cover.cover + np.array(cover_reads.cover)
            cover.max_reads += cover_reads.max_reads
            results[id_] = cover
        
    return results


def _get_Coverage_intervall(bam_files, intervalls, lib_scheme, n_thread=6, mapq=13, flag_in=256, flag_out=0):
    """
    """

    strand = intervalls[0]["strand"]
    
    intervalls = sorted(intervalls, key = lambda x: x["start"])
    start = intervalls[0]["start"]
    end = intervalls[-1]["end"]
    range_ = {
                "strand": intervalls[0]["strand"],
                "chr": intervalls[0]["chr"]
             }
    if start < end:
        range_["start"] = start
        range_["end"] = end
    else:
        range_["start"] = end
        range_["end"] = start
    
    n = range_["end"] - range_["start"]
    print("{}".format(n))

    args = [ (x, lib_scheme, range_, mapq, flag_in, flag_out)  for x in bam_files ]
    print("launching pool")
    with Pool(processes=n_thread) as pool: 
        cover_list = pool.map(parse_bam_mp, args)
    

    results = {}
    for i, bam in enumerate(bam_files):
        cover = Coverage(deepcopy(intervalls))
        cover.strand = strand
        cover.add_intron()
        cover.sort_intervalls()
        cover.make_intervall_scaled()
        cover.set_cover(np.array(cover_list[i]))
        results[bam] = cover

    return results


def get_attr(string):
    dico = {}
    spt = string.split(";")
    
    for x in spt:
        if x:
            dico[x.split()[0]] =  x.split()[1].replace('"', "")
    return dico


def gtf_to_dict(gtf_file, main_ = "gene_id"):
    print(main_)
    dico = {}
    with open(gtf_file) as f_in:
        """atm reads only genes"""
        for line in f_in:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            spt = line.strip().split("\t")
            try:
                chr_ = spt[0]
                start = int(spt[3])
                end = int(spt[4])
                strand =  spt[6]
                attr = get_attr(spt[-1])
                gene_id = attr[main_]
                gene_symbol = attr.get("gene_symbol", "None")
                type_ = spt[2]
            except:
                print(line)
                print(spt)
                raise
            
            if gene_id not in dico:
                dico[gene_id] = {
                    "chr": chr_,
                    "symbol" :  gene_symbol,
                    "strand" : strand,
                    "transcript": {}
                }
                
            if type_ == "gene":
    
                dico[gene_id]["start"] = start
                dico[gene_id]["end"] = end
            
            else:
                try:
                    transcript_id  = attr.get("transcript_id", "None")
                    transcript_symbol = attr.get("transcript_symbol", "None")
                    
                    if transcript_id not in dico[gene_id]["transcript"]:
                        dico[gene_id]["transcript"][transcript_id] = {
                        "transcript_symbol" : transcript_symbol,
                        "transcript_id" : transcript_id
                        }
                    
                    if type_ not in dico[gene_id]["transcript"][transcript_id]: 
                        dico[gene_id]["transcript"][transcript_id][type_] = []
                        
                    dico[gene_id]["transcript"][transcript_id][type_].append({
                        "chr": chr_,
                        "start": start, 
                        "end" : end,
                        "strand" : strand,
                    })
                except: 
                    continue
    for k, v in dico.items():
        # get bounds
        if v.get("start"):
            continue
        start, end = None, None
#        print(v)
        for transcript, values in v["transcript"].items():
            if start is None:
                start = values["transcript"][0]["start"]
                end = values["transcript"][0]["end"]
                continue
        
            if start > values["transcript"][0]["start"]:
                start = values["transcript"][0]["start"]
            if end > values["transcript"][0]["end"]:
                end = values["transcript"][0]["end"]
            dico[k]["start"] = start
            dico[k]["end"] = end

    return dico


def plot(dico_cover, N=None, intron="exon", out=None,
         intron_prop=0.45, bg_col="whitesmoke",
         dico_color=None, dico_label=None, title=None, linewidth=1):
    """
    
    """

    fig = plt.figure(figsize=(15, 5))

    # ax = fig.add_axes([0,0,1,1])
    ax = plt.gca()
    dico_legend = {}

    for geno_sample, cover in dico_cover.items():
        strand = cover.strand

        color = dico_color[geno_sample]

        if intron == "exon":
            cov = cover.get_cover_exon()
        elif intron == "intron":
            cov = cover.cover
        elif intron == "intron_partial":
            cover.compress_intron(max_ratio=intron_prop)
            cov = cover.cover_intron_compress
        else:
            raise AssertionError("you need to setup intron")
        
        if strand == "-":
                cov = cov[::-1]
        if N:
            cov = np.convolve([x for x in cov], np.ones(N)/N, mode='valid')
        ax.plot(range(len(cov)), cov, color=color, linewidth=linewidth)
        if dico_label:
            dico_legend["{}".format(dico_label[geno_sample])] = color

    # add bg
    max_=len(cov)
    if intron == "intron_partial":

        for e in cover.intervalls_scaled_intron_compress:
            if e["type"] == "intron":
                if strand == "-":
                    val = (max_ - e["end"],  max_ -  e["start"] )
                else:
                    val = (e["start"], e["end"])
                ax.axvspan(val[0], val[1], alpha=0.5, color="gainsboro")
        
    elif intron == "intron":
        ss = cover.intervalls_scaled[0]['start'] - cover.intervalls_scaled[-1]['end'] 
        print(ss)
        for e in cover.intervalls_scaled:
            if e["type"] == "intron":
                if strand == "-":
                    val = (max_ - e["end"], max_ -  e["start"] )
                else:
                    val = (e["start"], e["end"])
                span_size = val[1] - val[0]
                print(span_size, span_size/ ss)
                ax.axvspan(val[0], val[1], alpha=0.5, color='gainsboro')
 
    elif intron == "exon":
        if strand == "+":
            delta = 0
            cpt = 0
            for i, e in enumerate(cover.intervalls_scaled):
                if e["type"] == "exon":
                    l = e["end"] - e["start"]
                    if strand == "+":
                        val = (delta, l + delta)
                    
                    if cpt % 2 == 1:
                        ax.axvspan(val[0], val[1], alpha=0.5, color='gainsboro')
                    delta += l
                    cpt += 1
        else:
            delta = 0
            cpt = 0
            for i, e in enumerate(cover.intervalls_scaled[::-1]):
                if e["type"] == "exon":
                    l = e["end"] - e["start"]
                    val = (delta, l + delta)
                    if cpt % 2 == 1:
                        ax.axvspan(val[0], val[1], alpha=0.5, color='gainsboro')
                    delta += l
                    cpt += 1
    
    ax = plt.gca()

    plt.margins(x=0)
    plt.margins(y=0)
    if dico_label:
        col = [ (id_, col) for id_, col in sorted(dico_legend.items(), key=lambda x : len(x[0])) ]
        custom_lines = [Line2D([0], [0], color=x[1], lw=4) for x in col]
        plt.gca().legend(handles=custom_lines, labels=[x[0] for x in col],  bbox_to_anchor=(1.05, 1), loc='upper left'   )
    plt.gca().set_facecolor(bg_col)

    plt.title(title)
    plt.tight_layout()
    if out:
        plt.savefig(out)

    pass



# def merge_duplicate(dico_bam_cover):
    
#     reg_get_sample_lanes = re.compile("trimmed_([a-z_]+)_.*_(S\d+)_(L00\d)_VS_\S+(?=Aligned.*)")
#     results = {}
#     for bam, cover_reads in dico_bam_cover.items():

#         m = reg_get_sample_lanes.match(Path(bam).stem)
#         genotype = m.group(1)
#         sample = m.group(2)
#         lane = m.group(3)

#         id_ = "{}_{}".format(genotype, sample)

#         if id_ not in results:
#             results[id_] = cover_reads
#         else:
#             cover = results[id_]
#             cover.cover = cover.cover + np.array(cover_reads.cover)
#             cover.max_reads += cover_reads.max_reads
#             results[id_] = cover
        
#     return results


if __name__ == "__main__":


    parse = argparse.ArgumentParser(description="""""")
    
    group1 = parse.add_argument_group('Input / output options', '')

    group1.add_argument("--bam", "-b", nargs = "+", help="space separated path to bam files", required=True)
    group1.add_argument("--bam_dir",  help="Convenient arg used as base path for all bam, if you bam are all in the same dir,\
                        you can use this argument to gave this dir path only once")
    group1.add_argument("--gene", help="geneid to plot", required=True)
    group1.add_argument("--gtf", "-g", default="/lab/solexa_yamashita/people/Romain/References/MD6/gtf/dmel-all-r6.48.dico.pkl", help="path to a gtf file if you use this object provide a geneid/ it can also be a pickle file\
                        for drosophila you can use it is the defulat value  : /lab/solexa_yamashita/people/Romain/References/MD6/gtf/dmel-all-r6.48.dico.pkl")
    group1.add_argument("--transcript", "-t", help="by default we use the first transctipt we see, use this to provide the transcript id")
    group1.add_argument("--out", "-o", help="the output file")
    group1.add_argument("--lib", help="the sequencing scheme defaut rf", choices=["rf", "fr", "se", "uns"], default="rf")

    group2 = parse.add_argument_group('Plot options', '')
  
    group2.add_argument("--intron", help="behaviour for intron default 'exon' only", choices=["intron_partial", "intron", "exon"], default="exon")
    group2.add_argument("--color", "-c", nargs = "+", help="space separated list of color matching the bam file order, if empty will use tab10, will fail if more than 10 bam and no color specified")
    group2.add_argument("--label", "-l", nargs = "+", help="space separated list of label for legend shoudl match the order to bam files, if empty will use bam file name")
    group2.add_argument("--intron_proportion", help="compress intron to make the plot easy to read value betwene 0 and 1 default 0.45", default=0.45)
    group2.add_argument("--title", help="plot title", default=0.45)

    group3 = parse.add_argument_group('Others options', '')
    group3.add_argument("--thread", help="max cpu usage", type=int, default=6)
    group3.add_argument("--smoothing", type=int)
    group3.add_argument("--max_reads", nargs="+", type=int)
    group3.add_argument("--mapq", type=int, default=13, help="default 13")
    group3.add_argument("--normalize", help="normalize option expect to find log of STAR aligner, or a the result of flagstat.\
                        if those are not found/present you need to give it the uniquely mapped reads number\
                        using the max_reads option. example for sample1.aln_read.out.bam will look for the file sample1.flagstat", action = "store_true")
   # parse.add_argument("")
   # parse.add_argument("--intervalls", "-i", help="coma separated list with in order start, end, chr, strand, id several exons, just add several items")
 # TODO provide JSON input output for batch and rest api.
    args = parse.parse_args()

    try:
        assert (args.gtf and args.gene) or args.intervalls
    except AssertionError:
        logging.error("either provide a path to a gtf file and geneid or provide an intervall")
        raise

    bams = args.bam
    if args.bam_dir:
        bams = [args.bam_dir + "/" + x for x in bams]
    for e in bams:
        try: 
            assert Path(e).exists
        except:
            logging.error("{}  not found".format(e))
            raise
    print(bams[0])
    

    gtf = Path(args.gtf)
    if gtf.suffix == ".pkl":
        with open(args.gtf, "rb") as f_i:
         dico_gtf =  pickle.load(f_i)
    else:
        dico_gtf = gtf_to_dict(args.gtf)
    gene = dico_gtf[args.gene]["symbol"]
    transcript = args.transcript
    if not transcript:
        transcript = list(dico_gtf[args.gene]["transcript"].keys())[0]
    intervalls = dico_gtf[args.gene]["transcript"][transcript]["exon"]

    chr_in_bam = Rust_covpyo3.get_header(bam_path=bams[0])
    # annonation canbe quite inconsistent on chromosome naming
    # either adding chr on not.
    # here we try to rescue 
    chr_ = intervalls[0]["chr"]
    
    if chr_ not in chr_in_bam:
        if chr_.startswith("chr"):
            if chr_[2:] in chr_in_bam:
                for e in intervalls:
                    e["chr"] = e["chr"][2:]
        elif "chr" + chr_ in chr_in_bam:
            for e in intervalls:
                e["chr"] = "chr{}".format(e["chr"])
        else:
            raise AssertionError("the chromosome name do not match he chromosome name in the bam header")

    for e in intervalls:
        e["type"] = "exon"



    thread = min([args.thread, len(bams)])
    print(thread)

    results = _get_Coverage_intervall(bams, intervalls=intervalls, n_thread=thread, lib_scheme=args.lib, mapq=args.mapq)

    dico_colors = None
    dico_labels = None

    if args.color:
        dico_colors = dict((x[1], x[0]) for x in zip(args.color, bams))
    else:
        name = "tab10"
        cmap = mpl.colormaps[name]  # type: matplotlib.colors.ListedColormap
        colors = cmap.colors  # type: list
        dico_colors = dict((x, colors[i]) for i, x in enumerate(bams))

    if args.label:
        dico_labels = dict((x[1], x[0]) for x in zip(args.label, bams))
    
    else:
        dico_labels = dict((x, x.stem.split(".")[0]) for x in  bams)


    mpl.rcParams['pdf.fonttype'] = 42
   # mpl.rcParams['text.usetex'] = True

    if args.normalize:
        if args.max_reads:
            dico_max_read = dict((x[1], x[0]) for x in zip(args.max_reads, bams))
            for bam in results:
                results[bam].max_reads = dico_max_read[bam]
        else:
            for bam in results:
                results[bam].max_reads = get_max_read(bam)
        for bam in results:
            results[bam].normalize()
    title = args.title
    if not title:
        title = args.gene
    plot(results, N=args.smoothing, intron=args.intron, out=args.out, intron_prop=args.intron_proportion,
          dico_color=dico_colors,
           dico_label=dico_labels, title=title)

  