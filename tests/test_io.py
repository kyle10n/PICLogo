import os
import pytest
from types import SimpleNamespace

from phyinc.io import (
    infer_sequence_type, match_alphabets, read_sequences, check_accession_consistency,
    dna_alphabets, rna_alphabets, protein_alphabets,
    strip_domain_coords,
)

SYNTHETIC_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples', 'synthetic_data')


# --- match_alphabets ---

def test_dna_inference():
    alpha = match_alphabets(set('ACGT'), dna_alphabets)
    assert str(alpha) == 'ACGT'

    for c in 'ACGT':
        alpha = match_alphabets(set(c), dna_alphabets)
        assert str(alpha) == 'ACGT'


def test_match_alphabets_rna():
    alpha = match_alphabets(set('ACGU'), rna_alphabets)
    assert alpha is not None
    assert 'U' in str(alpha)


def test_match_alphabets_protein():
    alpha = match_alphabets(set('ACDEFGHIKLMNPQRSTVWY'), protein_alphabets)
    assert alpha is not None


def test_match_alphabets_subset_of_protein():
    # A small set of amino acids should match a protein alphabet
    alpha = match_alphabets({'A', 'C', 'D'}, protein_alphabets)
    assert alpha is not None


def test_match_alphabets_no_match():
    # Digit characters are not in any biological alphabet
    alpha = match_alphabets({'1', '2'}, dna_alphabets)
    assert alpha is None


def test_match_alphabets_empty_set():
    # An empty set is a subset of any set
    alpha = match_alphabets(set(), dna_alphabets)
    assert alpha is not None


# --- infer_sequence_type ---

def test_infer_sequence_type_dna():
    args = SimpleNamespace(type='dna')
    alphabet, color = infer_sequence_type(set('ACGT'), args)
    assert alphabet is not None
    assert color is not None


def test_infer_sequence_type_rna():
    args = SimpleNamespace(type='rna')
    alphabet, color = infer_sequence_type(set('ACGU'), args)
    assert alphabet is not None


def test_infer_sequence_type_invalid_type():
    args = SimpleNamespace(type='invalid')
    with pytest.raises(Exception, match="Bug in infer_sequence_type"):
        infer_sequence_type(set('ACGT'), args)


# --- read_sequences ---

def test_read_sequences_returns_dict():
    fa_file = os.path.join(SYNTHETIC_DIR, 'ex1.fa')
    args = SimpleNamespace(type='guess', ignore_coords=False)
    seq_dict, seq_type, seq_length, color_scheme = read_sequences(fa_file, 'fasta', args)
    assert isinstance(seq_dict, dict)
    assert len(seq_dict) > 0
    assert seq_type is not None
    assert seq_length == 10  # ex1.fa sequences are 10 residues long


def test_read_sequences_all_accessions_present():
    fa_file = os.path.join(SYNTHETIC_DIR, 'ex1.fa')
    args = SimpleNamespace(type='guess', ignore_coords=False)
    seq_dict, _, _, _ = read_sequences(fa_file, 'fasta', args)
    for name in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']:
        assert name in seq_dict


def test_read_sequences_explicit_type():
    fa_file = os.path.join(SYNTHETIC_DIR, 'ex1.fa')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    seq_dict, seq_type, _, _ = read_sequences(fa_file, 'fasta', args)
    assert seq_type is not None


def test_read_sequences_with_coords_no_domain_accessions():
    # ex1.fa has plain accessions (no /start-end), coords flag should have no extra effect
    fa_file = os.path.join(SYNTHETIC_DIR, 'ex1.fa')
    args_no_coords = SimpleNamespace(type='guess', ignore_coords=False)
    args_coords = SimpleNamespace(type='guess', ignore_coords=True)
    seq_dict_no, _, _, _ = read_sequences(fa_file, 'fasta', args_no_coords)
    seq_dict_yes, _, _, _ = read_sequences(fa_file, 'fasta', args_coords)
    assert set(seq_dict_no.keys()) == set(seq_dict_yes.keys())


def test_read_sequences_with_domain_coords(tmp_path):
    fa = tmp_path / 'domain.fa'
    fa.write_text('>PROT1/10-19\nACDEFGHIKL\n>PROT2/20-29\nACDEFGHIKL\n')
    args = SimpleNamespace(type='aa', ignore_coords=True)
    seq_dict, _, _, _ = read_sequences(str(fa), 'fasta', args)
    # Both full accession and stripped accession should be present
    assert 'PROT1/10-19' in seq_dict
    assert 'PROT1' in seq_dict
    assert 'PROT2/20-29' in seq_dict
    assert 'PROT2' in seq_dict


def test_read_sequences_inconsistent_length(tmp_path):
    fa = tmp_path / 'bad.fa'
    fa.write_text('>seq1\nACGT\n>seq2\nACGTT\n')
    args = SimpleNamespace(type='dna', ignore_coords=False)
    with pytest.raises(Exception, match="inconsistent"):
        read_sequences(str(fa), 'fasta', args)


# --- check_accession_consistency ---

class MockTerminal:
    def __init__(self, name):
        self.name = name

class MockTree:
    def __init__(self, leaf_names):
        self._terminals = [MockTerminal(n) for n in leaf_names]

    def get_terminals(self):
        return self._terminals


def test_check_accession_consistency_ok():
    seq_dict = {'a': 'seq', 'b': 'seq', 'c': 'seq'}
    tree = MockTree(['a', 'b', 'c'])
    check_accession_consistency(seq_dict, tree)  # should not raise


def test_check_accession_consistency_missing_from_seq():
    seq_dict = {'a': 'seq', 'b': 'seq'}
    tree = MockTree(['a', 'b', 'c'])
    with pytest.raises(ValueError, match="Tree leaves not found in sequence file"):
        check_accession_consistency(seq_dict, tree)


def test_check_accession_consistency_missing_from_tree():
    seq_dict = {'a': 'seq', 'b': 'seq', 'extra': 'seq'}
    tree = MockTree(['a', 'b'])
    with pytest.raises(ValueError, match="Sequence accessions not found in tree"):
        check_accession_consistency(seq_dict, tree)


def test_check_accession_consistency_both_directions():
    seq_dict = {'a': 'seq', 'extra_seq': 'seq'}
    tree = MockTree(['a', 'extra_leaf'])
    with pytest.raises(ValueError) as exc_info:
        check_accession_consistency(seq_dict, tree)
    msg = str(exc_info.value)
    assert "Tree leaves not found in sequence file" in msg
    assert "Sequence accessions not found in tree" in msg


def test_check_accession_consistency_empty():
    check_accession_consistency({}, MockTree([]))  # should not raise


def test_check_accession_consistency_names_listed_in_error():
    seq_dict = {'a': 'seq'}
    tree = MockTree(['a', 'missing1', 'missing2'])
    with pytest.raises(ValueError) as exc_info:
        check_accession_consistency(seq_dict, tree)
    msg = str(exc_info.value)
    assert 'missing1' in msg
    assert 'missing2' in msg


def test_strip_domain_coords():
    assert strip_domain_coords('Hubba') == 'Hubba'
    assert strip_domain_coords('Hubba_HUMAN') == 'Hubba_HUMAN'
    assert strip_domain_coords('ABC123_MOUSE') == 'ABC123_MOUSE'

    assert strip_domain_coords('Hubba/1-2') == 'Hubba'
    assert strip_domain_coords('Hubba_HUMAN/100-102') == 'Hubba_HUMAN'
    assert strip_domain_coords('ABC123_MOUSE/73-92') == 'ABC123_MOUSE'
