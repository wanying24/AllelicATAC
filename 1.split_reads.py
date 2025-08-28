
import pysam
import sys
import logging
import random
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def cal_score(read1, read2, indel_penalty=2):
    try:
        nm1 = int(read1.get_tag('NM'))
        nm2 = int(read2.get_tag('NM'))
    except (ValueError, KeyError):
        return None
    
    def analyze_indels(read):
        inserts = deletions = 0
        for operation, length in read.cigartuples:
            if operation == 1:
                inserts += length
            elif operation == 2: 
                deletions += length
        return inserts, deletions
    
    inserts1, deletions1 = analyze_indels(read1)
    inserts2, deletions2 = analyze_indels(read2)
    
    score = (
        (read1.query_alignment_length - nm1 - (inserts1 + deletions1) * indel_penalty) +
        (read2.query_alignment_length - nm2 - (inserts2 + deletions2) * indel_penalty)
    )
    return score

def process_pair(reads):
    best_score = float('-inf')
    best_pair = None
    reads1=[i for i in reads if i.is_read1]
    reads2=[i for i in reads if i.is_read2]
    for read1 in reads1:
        for read2 in reads2:
            if (read1.reference_name == read2.reference_name and
                read1.mapping_quality == read2.mapping_quality and
                read1.template_length == -read2.template_length and abs(read1.tlen) <700 ):
                score = cal_score(read1, read2)
                if score is not None and score > best_score:
                    best_score = score
                    best_pair = (read1, read2)
    
    return best_pair if best_pair else False

def process_list(read_list, cc1, cc2):
    reads_by_chrom = defaultdict(list)
    for read in read_list:
        if read.reference_name.startswith(cc1):
            reads_by_chrom[cc1].append(read)
        elif read.reference_name.startswith(cc2):
            reads_by_chrom[cc2].append(read)
    
    best_scores = {cc1+'_best': float('-inf'), cc2+'_best': float('-inf')}
    
    for chrom in [cc1, cc2]:
        if chrom in reads_by_chrom:
            pair = process_pair(reads_by_chrom[chrom])
            if pair:
                read1, read2 = pair
                score = cal_score(read1, read2)
                if score and score > best_scores[f"{chrom}_best"]:
                    best_scores[f"{chrom}_best"] = score
    
    cc1_score = best_scores[cc1+'_best']
    cc2_score = best_scores[cc2+'_best']
    
    if cc1_score > cc2_score:
        return f"{cc1}_unique"
    elif cc2_score > cc1_score:
        return f"{cc2}_unique"
    elif cc1_score != float('-inf') or cc2_score != float('-inf'):
        return 'common'
    else:
        return None


def process_read_pairs(input_bam, cc1, cc2):
    results = {
        f"{cc1}_unique": set(),
        f"{cc2}_unique": set(),
        "common": set()
    }

    with pysam.AlignmentFile(input_bam, "rb") as bam_file:
        current_query_name = None
        current_reads = []
        
        for read in bam_file.fetch(until_eof=True):
            if not read.is_mapped or 'chrM' in read.reference_name:
                continue
            if current_query_name and read.query_name != current_query_name:
                classification = process_list(current_reads, cc1, cc2)
                if classification:
                    results[classification].add(current_query_name)
                current_reads = []
            
            current_query_name = read.query_name
            current_reads.append(read)
        
        if current_reads:
            classification = process_list(current_reads, cc1, cc2)
            if classification:
                results[classification].add(current_query_name)
    
    return results
def modify_header(samfile, cc):
    original_dict = samfile.header.to_dict()
    new_sqs = []
    for sq in original_dict['SQ']:
        original_name = sq['SN']
        if original_name.startswith(cc):
            new_name = original_name[len(cc):]
            new_sqs.append({'SN': new_name, 'LN': sq['LN']})
        else:
            new_sqs.append(sq)

    original_dict['SQ'] = new_sqs
    return pysam.AlignmentHeader.from_dict(original_dict)
def write_to_file(input_bam, final_results, cc1, cc2, NAME):
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        modified_header_cc1 = modify_header(samfile, cc1)
        modified_header_cc2 = modify_header(samfile, cc2)
        common_part=list(final_results['common'])
        random.shuffle(common_part)
        half_size = len(common_part) // 2
        cc1_common_reads = set(common_part[:half_size])
        cc2_common_reads = set(common_part[half_size:])

        with pysam.AlignmentFile(f"{NAME}_{cc1}_common.bam", "wb", header=modified_header_cc1) as out_common1, \
             pysam.AlignmentFile(f"{NAME}_{cc2}_common.bam", "wb", header=modified_header_cc2) as out_common2, \
             pysam.AlignmentFile(f"{NAME}_{cc1}_unique.bam", "wb", header=modified_header_cc1) as out_unique1, \
             pysam.AlignmentFile(f"{NAME}_{cc2}_unique.bam", "wb", header=modified_header_cc2) as out_unique2:
            
            for read in samfile:
                if read.is_unmapped:
                    continue 

                if read.reference_name.startswith(cc1):
                    read.reference_id = modified_header_cc1.get_tid(read.reference_name[len(cc1):])
                elif read.reference_name.startswith(cc2):
                    read.reference_id = modified_header_cc2.get_tid(read.reference_name[len(cc2):])

                if read.reference_id == -1: 
                    logging.error(f"Reference {read.reference_name} not found in modified header.")
                    continue  
                if read.query_name in cc1_common_reads:
                    if read.reference_name.startswith(cc1):
                        out_common1.write(read)
                elif read.query_name in cc2_common_reads:
                    if read.reference_name.startswith(cc2):
                        out_common2.write(read)
                elif read.query_name in final_results[f"{cc1}_unique"]:
                    if read.reference_name.startswith(cc1):
                        out_unique1.write(read)
                elif read.query_name in final_results[f"{cc2}_unique"]:
                    if read.reference_name.startswith(cc2):
                        out_unique2.write(read)

    logging.info(f"Successfully finished writing BAM files for {NAME}.")


def count_properly_paired_reads(bam_file):
    reads = set()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch(until_eof=True):
            if read.is_mapped:
                reads.add(read.query_name)
    return len(reads)


def main(input_bam, cc1_prefix, cc2_prefix):
    cc1 = cc1_prefix.split('_')[1]
    cc2 = cc2_prefix.split('_')[1]
    NAME = cc1_prefix.split('_')[0]
    output = f"{NAME}_Sstats.txt"

    final_results = process_read_pairs(input_bam, cc1, cc2)


    write_to_file(input_bam, final_results, cc1, cc2, NAME)

    with open(output, 'w') as an:
        an.write(
            f"Mapped reads\t{count_properly_paired_reads(input_bam)}\n"
            f"{cc1}_unique\t{len(final_results[f'{cc1}_unique'])}\n"
            f"{cc2}_unique\t{len(final_results[f'{cc2}_unique'])}\n"
            f"common\t{len(final_results['common'])}\n"
        )

        for suffix in ['unique', 'common']:
            bam_file1 = f"{NAME}_{cc1}_{suffix}.bam"
            bam_file2 = f"{NAME}_{cc2}_{suffix}.bam"

            count1 = count_properly_paired_reads(bam_file1)
            count2 = count_properly_paired_reads(bam_file2)

            an.write(f"{cc1}_{suffix}_final\t{count1}\n")
            an.write(f"{cc2}_{suffix}_final\t{count2}\n")

    logging.info(f"Statistics written to {output}.")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input.bam> <cc1_prefix> <cc2_prefix>")
        sys.exit(1)

    input_bam, cc1_prefix, cc2_prefix = sys.argv[1:]
    main(input_bam, cc1_prefix, cc2_prefix)
