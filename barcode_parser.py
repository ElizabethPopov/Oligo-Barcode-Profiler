import regex as re

###---- Define anchor sequences for barcode & context extraction ----###
def compile_pattern(anchor1, anchor2, anchor3, barcode_len, context_len, mm1, mm2, mm3):
    '''
    Compile a regex pattern to match the structure of the reads.
    The structure of the relevent part in the reads is as follows:
    [anchor1]-[barcode]-[anchor2]-[context]-[anchor3]
    '''
    for name, anchor in zip(["anchor1", "anchor2", "anchor3"], [anchor1, anchor2, anchor3]):
        if not isinstance(anchor, str) or not anchor:
            raise ValueError(f"{name} must be a non-empty string.")
        if not re.fullmatch(r'[ACGT]+', anchor):
            raise ValueError(f"{name} contains invalid characters. Only A, C, G, T allowed.")

    for name, value in zip(
        ["barcode_len", "context_len", "mm1", "mm2", "mm3"],
        [barcode_len, context_len, mm1, mm2, mm3]
    ):
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a non-negative integer.")

    pattern = (
        f"(?P<anchor1>{anchor1}){{s<={mm1}}}"
        f"(?P<barcode>[ACGT]{{{barcode_len}}})"
        f"(?P<anchor2>{anchor2}){{s<={mm2}}}"
        f"(?P<context>[ACGT]{{{context_len}}})"
        f"(?P<anchor3>{anchor3}){{s<={mm3}}}"
    )
    return re.compile(pattern)


###---- search for anchor sequences in R1 & R2 reads with regex ----###
# This will help to extract the barcode and context sequences from the reads by not losing any reads in the process
# (if the defined anchor sequence is not present in the read because of mismatches, the read will be discarded)

###---- Function to extract barcode and context from R1 read using the defined regex pattern ----###
def extract_barcode_and_context(read, compiled_pattern):
    '''
    Extract barcode and context from a read using the compiled regex pattern.
    Returns a dictionary with the full sequence, barcode, and context if a match is found,
    otherwise returns None.
    '''
    if not isinstance(read, str):
        raise ValueError("The read must be a string.")
    if not compiled_pattern:
        raise ValueError("Compiled regex pattern is required.")
    if not compiled_pattern.pattern:
        raise ValueError("Compiled regex pattern is empty.")
    
    read = read.strip().upper()
    
    match = compiled_pattern.search(read)
    if match:
        return {
            'full_sequence': match.group(0),
            'barcode': match.group(2),
            'context': match.group(4)
        }
    else:
        return None

### Reverse complement function for R2 reads ###
def reverse_complement(seq):
    '''
    Return the reverse complement of a DNA sequence.
    '''
    from Bio.Seq import Seq
    
    if not seq:
        raise ValueError("The sequence cannot be empty.")
    if not isinstance(seq, str):
        raise ValueError("The sequence must be a string.")
    
    seq = seq.strip().upper()

    # Use Bio.Seq to get the reverse complement
    return str(Seq(seq).reverse_complement())


### Extractor for R2 (reverse-complemented) ###
def extract_barcode_and_context_R2(read, pattern):
    '''
    Extract barcode and context from a reverse-complemented read using the compiled regex pattern.
    Returns a dictionary with the full sequence, barcode, and context if a match is found,
    otherwise returns None.
    '''
    if not read:
        raise ValueError("The sequence cannot be empty.")
    if not isinstance(read, str):
        raise ValueError("The read must be a string.")
    if not pattern:
        raise ValueError("Compiled regex pattern is required.")
    if not pattern.pattern:
        raise ValueError("Compiled regex pattern is empty.")
    
    read = read.strip().upper()
    
    read_rc = reverse_complement(read)
    match = pattern.search(read_rc)
    if match:
        return {
            'full_sequence': match.group(0),
            'barcode': match.group(2),
            'context': match.group(4)
        }
    else:
        return None