# Handles paired-end FASTQ input, validates R1+R2 agreement, and extracts barcode-context pairs

import os
import regex as re
from barcode_parser import extract_barcode_and_context, extract_barcode_and_context_R2


### Core paired-read processor ###
def process_paired_reads(file_path_R1, file_path_R2, compiled_pattern):
    '''
    Process a pair of fastq file (forward and reverse reads of the same sample):
        - Find the compiled pattern within each read
        - Reverse complement the read if it is a R2 file
        - Compare between the same read in file R1 and file R2

    Returns:
    Validated matched reads

    '''
    if not file_path_R1.endswith(('.fastq', '.fq')) or not file_path_R2.endswith(('.fastq', '.fq')):
        raise ValueError("Expected FastQ files with '.fastq' or '.fq' extension.")

    # Check that files exist before parsing filenames
    if not os.path.isfile(file_path_R1) and not os.path.isfile(file_path_R2):
        raise FileNotFoundError(f"Neither R1 nor R2 files found: {file_path_R1}, {file_path_R2}")
    if not os.path.isfile(file_path_R1):
        raise FileNotFoundError(f"R1 file not found: {file_path_R1}")
    if not os.path.isfile(file_path_R2):
        raise FileNotFoundError(f"R2 file not found: {file_path_R2}")
    
    if not compiled_pattern:
        raise ValueError("Compiled regex pattern is required.")
    if not compiled_pattern.pattern:
        raise ValueError("Compiled regex pattern is empty.")

    matched_reads = []
    total_pairs = 0
    matched_pairs = 0

    with open(file_path_R1, 'r') as f_R1, open(file_path_R2, 'r') as f_R2:
        while True:
            header_R1 = f_R1.readline().strip().split('/')[0]
            header_R2 = f_R2.readline().strip().split('/')[0]
            
            if not header_R1 or not header_R2:
                break   # End of one of the files

            seq_R1 = f_R1.readline().strip()
            plus_R1 = f_R1.readline()     # +
            qual_R1 = f_R1.readline()     # quality

            seq_R2 = f_R2.readline().strip()
            plus_R2 = f_R2.readline()
            qual_R2 = f_R2.readline()

            # Validate FastQ structure
            if not all([seq_R1, plus_R1, qual_R1, seq_R2, plus_R2, qual_R2]):
                raise ValueError(
                    f"Malformed FastQ entry detected in file '{file_path_R1}' or '{file_path_R2}'. Ensure all reads contain 4 lines (header, sequence, +, quality)."
                )

            total_pairs += 1
            
            if header_R1 != header_R2:
                continue    # If headers don’t match- skip

            info_R1 = extract_barcode_and_context(seq_R1, compiled_pattern)
            info_R2 = extract_barcode_and_context_R2(seq_R2, compiled_pattern)

            if not info_R1 or not info_R2:
                continue    # Couldn't extract from both

            if info_R1['barcode'] == info_R2['barcode'] and info_R1['context'] == info_R2['context']:
                # Only keep reads where both R1 and R2 agree
                matched_reads.append({
                    'header': header_R1,
                    'sequence_R1': info_R1['full_sequence'],
                    'sequence_R2': info_R2['full_sequence'],
                    'barcode': info_R1['barcode'],
                    'context': info_R1['context'],
                    'source': 'R1+R2'
                })
                matched_pairs += 1

    print(f"Total paired reads processed: {total_pairs}")
    print(f"Validated matched pairs: {matched_pairs} ({round(matched_pairs / total_pairs * 100, 2)}%)")
    return matched_reads

### Load and process all pairs ###

def load_validated_reads(r1_file, r2_file, compiled_pattern):
    '''
    Load and process paired-end reads from R1 and R2 files.
    Returns:
        - sample name inferred from R1 filename
        - list of dictionaries with validated matched reads, their barcodes and contexts
    '''
    if not r1_file.endswith(('.fastq', '.fq')) or not r2_file.endswith(('.fastq', '.fq')):
        raise ValueError("Expected FastQ files with '.fastq' or '.fq' extension.")

    # Check that files exist before parsing filenames
    if not os.path.isfile(r1_file) and not os.path.isfile(r2_file):
        raise FileNotFoundError(f"Neither R1 nor R2 files found: {r1_file}, {r2_file}")
    if not os.path.isfile(r1_file):
        raise FileNotFoundError(f"R1 file not found: {r1_file}")
    if not os.path.isfile(r2_file):
        raise FileNotFoundError(f"R2 file not found: {r2_file}")
    
    if not compiled_pattern:
        raise ValueError("Compiled regex pattern is required.")
    if not compiled_pattern.pattern:
        raise ValueError("Compiled regex pattern is empty.")


    # Infer sample name from R1 file name and check both file names
    match_r1 = re.search(r'(\w+)_R1', os.path.basename(r1_file))
    if not match_r1:
        raise ValueError("Could not infer sample name from R1 file name.")

    match_r2 = re.search(r'(\w+)_R2', os.path.basename(r2_file))
    if not match_r2:
        raise ValueError("Could not infer sample name from R2 file name.")

    if match_r1.group(1) != match_r2.group(1):
        raise ValueError("R1 and R2 filenames appear to come from different samples.")

    sample = match_r1.group(1)

    print(f"Processing sample: {sample}")
    matched_reads = process_paired_reads(r1_file, r2_file, compiled_pattern)
    print('-' * 50)

    return sample, matched_reads