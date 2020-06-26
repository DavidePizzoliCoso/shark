import sys

import gffutils

def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtf_path,
                                 dbfn="{}.db".format(gtf_path),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        #print("Done.", file=sys.stderr)
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    return gtf

def build_gt_dict(gtf):
    tr2gene = {}
    for gene in gtf.features_of_type('gene'):
        gene_idx = gene.id
        for tr in gtf.children(gene, featuretype='transcript'):
            tr_idx = tr.id
            tr2gene[tr_idx] = gene_idx
    return tr2gene

def main():
    ssv_path = sys.argv[1]
    bed_path = sys.argv[2]
    gtf_path = sys.argv[3]

    gtf = open_gtf(gtf_path)
    tr2gene = build_gt_dict(gtf)

    truth = set()
    for line in open(bed_path, 'r'):
        truth.add(line.split("\t")[3])

    TP = 0
    FP = 0
    read_idx_tps = set()
    for line in open(ssv_path):
        read_idx, gene_idx = line.strip('\n').split(' ')
        real_tr_idx = read_idx.split(':')[2]
        if real_tr_idx not in tr2gene:
            FP+=1
            continue
        real_gene_idx = tr2gene[real_tr_idx]
        if gene_idx == real_gene_idx and read_idx in truth:
            if read_idx not in read_idx_tps:
                TP+=1
                read_idx_tps.add(read_idx)
        else:
            FP+=1

    FN = len(truth) - TP
    # for idx in truth:
    #     if idx not in read_idx_tps:
    #         print(idx)

    print('TP', 'FP', 'FN', 'P', 'R', sep='\t')
    # print('TP', 'FP', 'FN', 'P', 'R')
    P = round(TP/(TP+FP), 3) if TP+FP != 0 else 0
    R = round(TP/(TP+FN), 3) if TP+FN != 0 else 0
    print(TP, FP, FN, P, R, sep='\t')
    # print(TP, FP, FN, P, R)

if __name__ == '__main__':
    main()
