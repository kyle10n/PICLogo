import os
import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace

from Bio import Phylo
from io import StringIO

from phyinc.phyinc import (
    collapse_bifurcations,
    compute_pic_array,
    compute_ed_array,
    compute_ed_weights,
    compute_yule_array,
    compute_yule_weights,
    convert_matrix_to_array,
)
from phyinc.main import output_info, setup_argparse, validate_path
from phyinc.io import unambiguous_dna_alphabet, unambiguous_protein_alphabet

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples', 'synthetic_data')


class MockClade:
    """Minimal stand-in for a Bio.Phylo clade used in PIC calculations."""
    def __init__(self, branch_length, name=None, clades=None):
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []



# --- collapse_bifurcations ---

def test_collapse_bifurcations_two_children_equal_branches():
    mat_a = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    mat_b = np.array([[0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    seq_dict = {child_a: mat_a, child_b: mat_b}
    branch, seq = collapse_bifurcations([child_a, child_b], seq_dict, np.zeros((4, 2)))
    assert branch == pytest.approx(0.5)
    np.testing.assert_array_almost_equal(seq, 0.5 * (mat_a + mat_b))


def test_collapse_bifurcations_felsenstein_weighting():
    # Child a has a SHORT branch (0.1) and child b has a LONG branch (0.9).
    # Under Felsenstein's model, the short-branch child is a more reliable
    # estimate of the ancestor (higher precision = 1/branch_length), so it
    # should receive weight b_b/(b_a+b_b) = 0.9, NOT its own branch 0.1.
    mat_a = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])  # all A
    mat_b = np.array([[0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])  # all C
    child_a = MockClade(0.1)
    child_b = MockClade(0.9)
    seq_dict = {child_a: mat_a, child_b: mat_b}
    _, seq = collapse_bifurcations([child_a, child_b], seq_dict, np.zeros((4, 2)))
    # Child a (short branch) should dominate: weight = 0.9/1.0 = 0.9
    assert seq[0, 0] == pytest.approx(0.9)   # A fraction at position 0
    assert seq[0, 1] == pytest.approx(0.1)   # C fraction at position 0


def test_collapse_bifurcations_both_zero_branches():
    sentinel = np.zeros((4, 2))
    child_a = MockClade(0.0)
    child_b = MockClade(0.0)
    seq_dict = {child_a: np.ones((4, 2)), child_b: np.ones((4, 2))}
    branch, seq = collapse_bifurcations([child_a, child_b], seq_dict, sentinel)
    assert branch == 0
    assert seq is sentinel


def test_collapse_bifurcations_three_children():
    mat = np.eye(4, 3)
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    child_c = MockClade(1.0)
    seq_dict = {child_a: mat.copy(), child_b: mat.copy(), child_c: mat.copy()}
    branch, seq = collapse_bifurcations([child_a, child_b, child_c], seq_dict, np.zeros((4, 3)))
    assert isinstance(branch, float)
    assert seq.shape == mat.shape


# --- convert_matrix_to_array ---

def test_convert_matrix_to_array_shape():
    matrix = np.zeros((4, 5))
    result = convert_matrix_to_array(matrix, unambiguous_dna_alphabet, 5)
    assert result.shape == (5, 4)


def test_convert_matrix_to_array_values():
    matrix = np.zeros((4, 2))
    matrix[0][0] = 1.0   # A at position 0
    matrix[2][1] = 0.5   # G at position 1
    result = convert_matrix_to_array(matrix, unambiguous_dna_alphabet, 2)
    assert result[0][0] == pytest.approx(1.0)   # position 0, char A (index 0)
    assert result[1][2] == pytest.approx(0.5)   # position 1, char G (index 2)



# --- output_info ---

def test_output_info_default_format():
    args = SimpleNamespace(outfile=None, format=None, seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert str(outfile) == 'seqs.fa_logo.pdf'
    assert formatter is not None


def test_output_info_explicit_format():
    args = SimpleNamespace(outfile=None, format='png', seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert str(outfile) == 'seqs.fa_logo.png'


def test_output_info_explicit_outfile():
    args = SimpleNamespace(outfile='result.eps', format=None, seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert Path(outfile).name == 'result.eps'


def test_output_info_outfile_no_suffix():
    args = SimpleNamespace(outfile='result', format=None, seq_filename='seqs.fa')
    with pytest.raises(ValueError, match="suffix"):
        output_info(args)


def test_output_info_outfile_unsupported_format():
    args = SimpleNamespace(outfile='result.bmp', format=None, seq_filename='seqs.fa')
    with pytest.raises(ValueError, match="not a valid output format"):
        output_info(args)


# --- validate_path ---

def test_validate_path_existing_file(tmp_path):
    f = tmp_path / 'test.txt'
    f.write_text('hello')
    result = validate_path(str(f))
    assert result == Path(str(f))


def test_validate_path_nonexistent():
    with pytest.raises(FileNotFoundError):
        validate_path('/no/such/file.tree')


def test_validate_path_directory(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_path(str(tmp_path))


# --- compute_pic_array ---

def test_compute_pic_array_shape():
    from types import SimpleNamespace
    from Bio import Phylo, SeqIO
    import io as _io
    import phyinc.io as phyinc_io

    fa_file = os.path.join(EXAMPLES_DIR, 'ex1.fa')
    tree_file = os.path.join(EXAMPLES_DIR, 'ex1_t1.tree')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    alignment, seq_type, seq_length, _ = phyinc_io.read_sequences(fa_file, 'fasta', args)
    tree = Phylo.read(tree_file, 'newick')

    array = compute_pic_array(tree, alignment, seq_type, seq_length)

    assert array.shape == (seq_length, len(seq_type))


def test_compute_pic_array_values_sum_to_one():
    from types import SimpleNamespace
    from Bio import Phylo
    import phyinc.io as phyinc_io

    fa_file = os.path.join(EXAMPLES_DIR, 'ex1.fa')
    tree_file = os.path.join(EXAMPLES_DIR, 'ex1_t1.tree')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    alignment, seq_type, seq_length, _ = phyinc_io.read_sequences(fa_file, 'fasta', args)
    tree = Phylo.read(tree_file, 'newick')

    array = compute_pic_array(tree, alignment, seq_type, seq_length)

    # Each row is a weighted average of one-hot vectors, so rows should sum to 1.
    row_sums = array.sum(axis=1)
    np.testing.assert_allclose(row_sums, np.ones(seq_length), atol=1e-6)


# ---------------------------------------------------------------------------
# Helpers shared by ED and Yule tests
# ---------------------------------------------------------------------------

def _parse_newick(s):
    return Phylo.read(StringIO(s), 'newick')


# Balanced binary tree: 4 leaves, symmetric topology.
BALANCED_NWK  = '((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);'

# Caterpillar tree: d is maximally isolated (connected directly to the root).
CATERPILLAR_NWK = '(((a:0.1,b:0.1):0.4,c:0.5):0.5,d:1.0);'


# ---------------------------------------------------------------------------
# compute_yule_weights
# ---------------------------------------------------------------------------

def test_yule_weights_sum_to_one():
    tree = _parse_newick(BALANCED_NWK)
    weights = compute_yule_weights(tree)
    assert sum(weights.values()) == pytest.approx(1.0)


def test_yule_weights_balanced_tree_uniform():
    # A perfectly balanced tree should give every leaf equal weight.
    tree = _parse_newick(BALANCED_NWK)
    weights = compute_yule_weights(tree)
    for w in weights.values():
        assert w == pytest.approx(0.25)


def test_yule_weights_caterpillar_isolated_leaf_highest():
    # Leaf d is attached directly to the root — it represents the most
    # unique evolutionary history and should receive the largest weight.
    tree = _parse_newick(CATERPILLAR_NWK)
    weights = compute_yule_weights(tree)
    assert weights['d'] == pytest.approx(0.5)
    assert weights['d'] > weights['c'] > weights['a']
    assert weights['a'] == pytest.approx(weights['b'])


def test_yule_weights_all_leaves_present():
    tree = _parse_newick(BALANCED_NWK)
    weights = compute_yule_weights(tree)
    assert set(weights.keys()) == {'a', 'b', 'c', 'd'}


def test_yule_weights_single_leaf():
    tree = _parse_newick('lone;')
    weights = compute_yule_weights(tree)
    assert weights == {'lone': pytest.approx(1.0)}


# ---------------------------------------------------------------------------
# compute_ed_weights
# ---------------------------------------------------------------------------

def test_ed_weights_sum_to_one():
    tree = _parse_newick(BALANCED_NWK)
    weights = compute_ed_weights(tree)
    assert sum(weights.values()) == pytest.approx(1.0)


def test_ed_weights_balanced_tree_uniform():
    # Equal branch lengths + symmetric topology → uniform weights.
    tree = _parse_newick(BALANCED_NWK)
    weights = compute_ed_weights(tree)
    for w in weights.values():
        assert w == pytest.approx(0.25)


def test_ed_weights_no_branch_lengths_raises():
    tree = _parse_newick('((a,b),(c,d));')
    with pytest.raises(ValueError, match="no branch lengths"):
        compute_ed_weights(tree)


# ---------------------------------------------------------------------------
# Trees without branch lengths: PIC/ED raise, Yule does not
# ---------------------------------------------------------------------------

NO_BRANCH_NWK = '((a,b),(c,d));'


def test_pic_array_no_branch_lengths_raises():
    tree = _parse_newick(NO_BRANCH_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    with pytest.raises(ValueError, match="no branch lengths"):
        compute_pic_array(tree, alignment, unambiguous_dna_alphabet, 4)


def test_ed_array_no_branch_lengths_raises():
    tree = _parse_newick(NO_BRANCH_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    with pytest.raises(ValueError, match="no branch lengths"):
        compute_ed_array(alignment, tree, unambiguous_dna_alphabet, 4)


def test_yule_array_no_branch_lengths_succeeds():
    tree = _parse_newick(NO_BRANCH_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    array = compute_yule_array(alignment, tree, unambiguous_dna_alphabet, 4)
    assert array.shape == (4, len(unambiguous_dna_alphabet))


# ---------------------------------------------------------------------------
# compute_yule_array
# ---------------------------------------------------------------------------

def _make_dna_alignment(seqs):
    """Return a dict {name: record} from a list of (name, seq) pairs."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    return {name: SeqRecord(Seq(seq), id=name) for name, seq in seqs}


def test_yule_array_shape():
    tree = _parse_newick(BALANCED_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    array = compute_yule_array(alignment, tree, unambiguous_dna_alphabet, 4)
    assert array.shape == (4, len(unambiguous_dna_alphabet))


def test_yule_array_rows_sum_to_one_no_gaps():
    tree = _parse_newick(BALANCED_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    array = compute_yule_array(alignment, tree, unambiguous_dna_alphabet, 4)
    np.testing.assert_allclose(array.sum(axis=1), np.ones(4), atol=1e-6)


def test_yule_array_matches_example_files():
    # Smoke test: array from example data has the right shape and sane values.
    import phyinc.io as phyinc_io
    fa_file   = os.path.join(EXAMPLES_DIR, 'ex1.fa')
    tree_file = os.path.join(EXAMPLES_DIR, 'ex1_t1.tree')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    alignment, seq_type, seq_length, _ = phyinc_io.read_sequences(fa_file, 'fasta', args)
    tree = Phylo.read(tree_file, 'newick')

    array = compute_yule_array(alignment, tree, seq_type, seq_length)

    assert array.shape == (seq_length, len(seq_type))
    assert np.all(array >= 0)
    np.testing.assert_allclose(array.sum(axis=1), np.ones(seq_length), atol=1e-6)


# ---------------------------------------------------------------------------
# compute_ed_array
# ---------------------------------------------------------------------------

def test_ed_array_shape():
    tree = _parse_newick(BALANCED_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    array = compute_ed_array(alignment, tree, unambiguous_dna_alphabet, 4)
    assert array.shape == (4, len(unambiguous_dna_alphabet))


def test_ed_array_rows_sum_to_one_no_gaps():
    tree = _parse_newick(BALANCED_NWK)
    alignment = _make_dna_alignment([('a','ACGT'),('b','ACGT'),('c','TGCA'),('d','TGCA')])
    array = compute_ed_array(alignment, tree, unambiguous_dna_alphabet, 4)
    np.testing.assert_allclose(array.sum(axis=1), np.ones(4), atol=1e-6)


# ---------------------------------------------------------------------------
# -w / --weighting argument
# ---------------------------------------------------------------------------

def test_weighting_default_is_pic():
    ap = setup_argparse()
    args = ap.parse_args(['seqs.fa', 'tree.nwk'])
    assert args.weighting == 'pic'


def test_weighting_short_flag_ed():
    ap = setup_argparse()
    args = ap.parse_args(['seqs.fa', 'tree.nwk', '-w', 'ed'])
    assert args.weighting == 'ed'


def test_weighting_long_flag_yule():
    ap = setup_argparse()
    args = ap.parse_args(['seqs.fa', 'tree.nwk', '--weighting', 'yule'])
    assert args.weighting == 'yule'


def test_weighting_explicit_pic():
    ap = setup_argparse()
    args = ap.parse_args(['seqs.fa', 'tree.nwk', '-w', 'pic'])
    assert args.weighting == 'pic'


def test_weighting_invalid_choice_raises():
    ap = setup_argparse()
    with pytest.raises(SystemExit):
        ap.parse_args(['seqs.fa', 'tree.nwk', '-w', 'gsc'])



