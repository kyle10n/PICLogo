from Bio import SeqIO
import logging
import re


domain_coord_pattern = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)/(\d+)-(\d+)$')


class Alphabet:
    """Minimal alphabet for sequence type representation."""
    def __init__(self, chars, name):
        self._chars = chars
        self.name = name

    def letters(self):
        return self._chars

    def __len__(self):
        return len(self._chars)

    def __iter__(self):
        return iter(self._chars)

    def __contains__(self, c):
        return c in self._chars

    def __str__(self):
        return self._chars

    def __repr__(self):
        return f'Alphabet({self._chars!r}, name={self.name!r})'


unambiguous_dna_alphabet     = Alphabet('ACGT',                    'dna')
dna_alphabet                 = Alphabet('ACGTRYSWKMBDHVN',         'dna')
reduced_nucleic_alphabet     = Alphabet('ACGTU',                   'dna')
nucleic_alphabet             = Alphabet('ACGTURYSWKMBDHVN',        'dna')
unambiguous_rna_alphabet     = Alphabet('ACGU',                    'rna')
rna_alphabet                 = Alphabet('ACGURYSWKMBDHVN',         'rna')
unambiguous_protein_alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY',    'protein')
reduced_protein_alphabet     = Alphabet('ACDEFGHIKLMNPQRSTVWYBXZ', 'protein')
protein_alphabet             = Alphabet('ACDEFGHIKLMNPQRSTVWYBXZU*', 'protein')

dna_alphabets = [
    unambiguous_dna_alphabet,
    reduced_nucleic_alphabet,
    dna_alphabet,
    nucleic_alphabet,
]
rna_alphabets = [
    unambiguous_rna_alphabet,
    rna_alphabet,
]
protein_alphabets = [
    unambiguous_protein_alphabet,
    reduced_protein_alphabet,
    protein_alphabet,
]

all_alphabets = dna_alphabets + rna_alphabets + protein_alphabets


def match_alphabets(char_set, available_alphabets):
    for alphabet in available_alphabets:
        if char_set.issubset(set(alphabet)):
            return alphabet     # Found it!
    return None


def infer_sequence_type(char_set, args):
    """
    Guess sequence type by matching alphabets.
    If the user gave a command-line instruction, then
    choose an alphabet that best matches the observed
    character set. If the input contains ambiguity
    symbols, we will assume it is a protein sequence
    and the user will need to specify the sequence type.

    We will also suggest a logomaker color scheme string.

    Returns: an Alphabet object and a logomaker color scheme string.
    """
    char_set = char_set - {'-', '.'}   # strip alignment gap characters

    if args.type == 'guess':           # User leaves it to us to choose
        choice = match_alphabets(char_set, all_alphabets)
        if choice is not None and choice.name in ('dna', 'rna'):
            return choice, 'classic'
        choice = match_alphabets(char_set, protein_alphabets)
        if choice is not None:
            return choice, 'hydrophobicity'
        unknown = ', '.join(sorted(char_set))
        raise ValueError(f"Cannot determine sequence type — unrecognised characters: {unknown}")

    elif args.type == 'dna':
        choice = match_alphabets(char_set, dna_alphabets)
        if choice is None:
            unknown = ', '.join(sorted(char_set - set(nucleic_alphabet)))
            raise ValueError(f"Sequences declared as DNA but contain non-DNA characters: {unknown}")
        return choice, 'classic'
    elif args.type == 'rna':
        choice = match_alphabets(char_set, rna_alphabets)
        if choice is None:
            unknown = ', '.join(sorted(char_set - set(nucleic_alphabet)))
            raise ValueError(f"Sequences declared as RNA but contain non-RNA characters: {unknown}")
        return choice, 'classic'
    elif args.type == 'aa':
        choice = match_alphabets(char_set, protein_alphabets)
        if choice is None:
            unknown = ', '.join(sorted(char_set - set(protein_alphabet)))
            raise ValueError(f"Sequences declared as protein but contain unrecognised characters: {unknown}")
        return choice, 'hydrophobicity'
    else:
        raise Exception("Bug in infer_sequence_type. Please report!")


def read_sequences(filename, filetype, args):
    """
    Read the input alignment and perform basic length checks.
    Also, determine what kind of bio sequence we read.

    Returns seq_dict, seq_type, seq_length, and a suggested coloring scheme.
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    seq_dict = {}
    characters = set()

    seq_length = None
    for acc, seq_str in record_dict.items():
        if len(seq_str) != seq_length:
            if seq_length:
                raise Exception(
                    "Sequence lengths from provided fastafile are inconsistent"
                )
            else:
                seq_length = len(seq_str)

        seq_str = seq_str.upper()
        characters.update(set(seq_str))
        if args.ignore_coords:
            m = domain_coord_pattern.match(acc) # Is acc on the form "HUBBA/17-35"?
            if m:
                prot_acc, domain_start, domain_end = m.groups()
                if prot_acc in seq_dict:
                    raise IOError(f"'{prot_acc}' is a protein appearing twice, probably because you have two domains from the same protein in the input. If so, you must submit a tree inferred on the domain sequences, not on the proteins.")
                seq_dict[prot_acc] = seq_str # Map protein accession to seq

        # Always insert for full accession, even if domain start and end are included
        seq_dict[acc] = seq_str      # Map domain accession to seq

    seq_type, coloring_scheme = infer_sequence_type(characters, args)
    return seq_dict, seq_type, seq_length, coloring_scheme


def strip_domain_coords(acc):
    """
    acc: a string containing a sequence accession, possibly with domain 
         coordinates, such as 'ABC_HUMAN/73-122'.
   
    returns: same accession, but without domain coordinates.
    """
    m = domain_coord_pattern.match(acc)
    if m:
        return m[1]
    else:
        return acc


def check_accession_consistency(seq_dict, tree, ignore_domain_coords=False):
    """
    Verify that sequence accessions and tree leaf names agree.

    Raises ValueError if any sequence accession has no corresponding
    leaf in the tree.
    """
    tree_leaves = {clade.name for clade in tree.get_terminals()}
    seq_accessions = set(seq_dict.keys())
    if ignore_domain_coords:
        for acc in seq_dict.keys():
            m = domain_coord_pattern.match(acc)
            if m:
                seq_accessions.remove(acc)

    missing_from_seq = tree_leaves - seq_accessions
    missing_from_tree = seq_accessions - tree_leaves

    messages = []
    if missing_from_seq:
        n = len(missing_from_seq)
        if n > 4:
            name_lst = sorted(list(missing_from_seq))
            names = ", ".join(name_lst[:4])
            messages.append(f"{n} tree leaves not found in sequence file: {names} and {n-3} more.")
        else:
            names = ", ".join(sorted(missing_from_seq))
            messages.append(f"Tree leaves not found in sequence file: {names}")

    if missing_from_tree:
        n = len(missing_from_tree)
        if n > 4:
            name_lst = sorted(list(missing_from_tree))
            names = ", ".join(name_lst[:4])
            messages.append(f"{n} sequence accessions not found in tree: {names} and {n-3} more.")
        else:
            name_lst = sorted(list(missing_from_tree))
            names = ", ".join(name_lst[:4])
            messages.append(f"Sequence accessions not found in tree: {names}")

    if messages:
        raise ValueError("\n".join(messages))
