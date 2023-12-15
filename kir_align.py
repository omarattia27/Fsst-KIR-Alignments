import sys
import subprocess as sp
import pickle
import mappy as mp

from mappy_trial import align_with_mappy
from common import timing

class Database:
    def __init__(self, path):
        from common import Region, Allele, Gene
        with open(path, "rb") as fo:
            (self.genes, _, _, self.minor_to_major) = pickle.load(fo)

    def alleles(self):
        for g in self.genes.values():
            yield from g.alleles.values()

def create_file_rep(alleles, db, filename):
    with open(filename, 'w') as file:
        for allele in alleles:
            file.write(f'>{allele.gene.name}.{allele.name}\n')
            file.write(f'{allele.seq}\n')

def create_file(clusters, db, filename):
    with open(filename, 'w') as file:
        for gene in clusters.keys():
            for rep in clusters[gene]:
                file.write(f'>{gene}.{rep}\n')
                file.write(f'{db.genes[gene][rep].seq}\n')

def alignment_score_ops(ops1, ops2, section):
    #No INDEL
    score = 0
    start_index, end_index = section 
    for op in ops1:
        if op not in ops2 and op[0]<=end_index and op[0]>=start_index:
            if "del" not in op[1]:
                score+=1
            else:
                score += len(op[1]) - 3
    for op in ops2:
        if op not in ops1 and op[0]<=end_index and op[0]>=start_index:
            if "del" not in op[1]:
                score+=1
            else:
                score += len(op[1]) - 3
    
    return score

def create_cluster(gene, db):
    sections = [(0,int(len(gene.seq)))]
    reprsentitives = {}
    for section in sections:
        for allele in gene.alleles.values():
            if len(reprsentitives.keys()) != 0:
                reps = reprsentitives.keys()
                added = False
                for rep in reps:
                    rep_allele = db.genes[gene.name].alleles[rep]
                    if alignment_score_ops(allele.ops,rep_allele.ops,section) <= 20:
                        reprsentitives[rep].add(allele)
                        added = True
                if added == False:
                    reprsentitives[allele.name] = set()
                    
            else:
                reprsentitives[allele.name] = set()
                
    return reprsentitives

def align_with_mappy_original(database_fasta, input_query):
    # Load or build the index
    a = mp.Aligner(database_fasta)

    if not a:
        raise Exception("ERROR: failed to load/build index")
    allignments = []
    # Read the input query sequence
    with open(input_query, 'r') as file:
        #query_sequence = file.read().strip()
        index = 0
        current_reading = ""
        for line in file:
            if index%2 == 0:
                current_reading += line
            else:
                current_reading += line
                current_reading = current_reading.strip()
                for hit in a.map(current_reading):
                    #print(f"{hit.ctg}\t{hit.r_st}\t{hit.r_en}\t{hit.cigar_str}")
                    allignments.append(f"{hit.ctg}\t{hit.r_st}\t{hit.r_en}\t{hit.cigar_str}")
                current_reading = ""
            index +=1

    return allignments

def align_with_mappy(database_fasta, a_clusters, input_query):
    # Load or build the index
    a = mp.Aligner(database_fasta)

    if not a:
        raise Exception("ERROR: failed to load/build index")
    allignments = []
    # Read the input query sequence
    with open(input_query, 'r') as file:
        #query_sequence = file.read().strip()
        index = 0
        current_reading = ""
        for line in file:
            if index%2 == 0:
                current_reading += line
            else:
                current_reading += line
                current_reading = current_reading.strip()
                for hit in a.map(current_reading):
                    #print(f"{hit.ctg}\t{hit.r_st}\t{hit.r_en}\t{hit.cigar_str}")
                    for hit_in_cluster in a_clusters[ f'kirdb_{hit.ctg}.fa'].map(current_reading):
                        allignments.append(f"{hit_in_cluster.ctg}\t{hit_in_cluster.r_st}\t{hit_in_cluster.r_en}\t{hit_in_cluster.cigar_str}")
                    allignments.append(f"{hit.ctg}\t{hit.r_st}\t{hit.r_en}\t{hit.cigar_str}")
                current_reading = ""
            index +=1

    return allignments

if __name__ == "__main__":
    with timing("initialization"):
        db = Database("kir.pickle")

        clusters = {}
        for gene in db.genes.keys():
            clusters[gene] = create_cluster(db.genes[gene], db)

        create_file(clusters, db, "kirdb_reps.fa")

        clusters_files_names = []
        for gene in clusters.keys():
            for rep in clusters[gene].keys():
                create_file_rep(clusters[gene][rep], db, f'kirdb_{gene}.{rep}.fa')
                clusters_files_names.append(f'kirdb_{gene}.{rep}.fa')

    total_alignments: list
    a_clusters = {cluster:mp.Aligner(cluster) for cluster in clusters_files_names}

    total_alignments_enhanced = []
    with timing("total (enhanced)"):
        total_alignments_enhanced = align_with_mappy("kirdb_reps.fa", a_clusters, sys.argv[1])
        print(len(total_alignments_enhanced))
        #print(total_alignments_enhanced)
    
    total_alignments = []
    with timing("total"):
        total_alignments = align_with_mappy_original("kirdb.fa", sys.argv[1])
        print(len(total_alignments))
        #print(total_alignments)
    total_alignments_enhanced_set = set(total_alignments_enhanced)
    print((len([allignment for allignment in total_alignments if allignment not in total_alignments_enhanced_set])/len(total_alignments))*100)
    
