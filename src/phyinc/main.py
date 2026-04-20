import argparse
import logging
import sys

import matplotlib.pyplot as plt

from pathlib import Path
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from importlib.metadata import version

import phyinc.io as io
from phyinc.colorhelper import decide_color_scheme
from phyinc.phyinc import create_logo, compute_pic_array, compute_ed_array, compute_yule_array


output_formats = ["pdf", "eps", "png", "jpg", "svg"]

description_str = """
Make sequence logos using Felsenstein's phylogenetically independent contrast
method to take evolution into account.
"""


def setup_argparse():
    parser = argparse.ArgumentParser(
        description=description_str + f" Version {version('phyinc')}."
    )
    parser.add_argument(
        "seq_filename",
        type=str,
        help="Path to the FASTA file. Sequence names may be in the format ACC/start-end.",
    )
    parser.add_argument(
        "tree_filename", type=str, help="Path to the tree file (Newick format)."
    )
    parser.add_argument(
        "-w",
        "--weighting",
        choices=["pic", "ed", "yule"],
        default="pic",
        dest="weighting",
        help=(
            "Sequence weighting method. "
            "'pic': phylogenetically independent contrasts (Felsenstein 1985, default); "
            "'ed': fair proportion evolutionary distinctiveness (Isaac et al. 2007); "
            "'yule': Yule-model topology weights."
        ),
    )
    parser.add_argument(
        "-c",
        "--color-scheme",
        choices=["monochrome", "nucleotide", "base_pairing", "hydrophobicity",
                 "chemistry", "charge", "taylor", "wesanderson", "skylign_aa", "skylign_dna", "random", "guess"],
        default="guess",
        help="Color scheme. 'guess' picks based on sequence type. Default: %(default)s.",
    )
    parser.add_argument(
        "-ic",
        "--ignore-coords",
        action="store_true",
        help="Ignore domain coordinates from accessions when matching to tree leaves. For example, if an accession is 'ABC_HUMAN/17-33' then '/17-33' is ignored.",
    )
    parser.add_argument(
        "-t",
        "--type",
        choices=["aa", "dna", "rna", "guess"],
        default="guess",
        help="Sequence type. Default: %(default)s.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Output filename. Format is inferred from the file suffix.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=output_formats,
        default="pdf",
        help=f"Output format. Ignored if --outfile has a recognised suffix. Default: %(default)s.",
    )
    parser.add_argument(
        "--export",
        metavar="filename",
        type=str,
        help="Write the weighted frequency matrix to a tab-separated file and exit.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Title to display above the logo.",
    )
    parser.add_argument(
        "-cl",
        "--command-line",
        action="store_true",
        help="Show the phyinc command line beneath the logo, ie how phyinc was called.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Write extra information to stderr.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"phyinc version {version('phyinc')}",
    )
    return parser


def output_info(args):
    """
    Determine output filename and format from args.
    Returns (outfile, graphics_format).
    """
    if args.outfile:
        outfile = Path(args.outfile)
        graphics_format = outfile.suffix.strip(".")
        if not graphics_format:
            raise ValueError(
                f"Give your outfile a suffix that determines output format, one of {output_formats}"
            )
        if graphics_format not in output_formats:
            raise ValueError(
                f'"{graphics_format}" is not a valid output format. Use one of {output_formats}.'
            )
    else:
        graphics_format = args.format if args.format else "pdf"
        outfile = args.seq_filename + "_logo." + graphics_format

    return outfile, graphics_format


def validate_path(filename):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"Cannot open {filename}.")
    return path


def export_weighted_data(tree, alignment, seq_type, seq_length, filename, method='pic'):
    """
    Write the weighted frequency matrix as a tab-separated file.
    The weighting method is selected by the method argument ('pic', 'ed', or 'yule').
    """
    try:
        if method == 'ed':
            array = compute_ed_array(alignment, tree, seq_type, seq_length)
        elif method == 'yule':
            array = compute_yule_array(alignment, tree, seq_type, seq_length)
        else:
            array = compute_pic_array(tree, alignment, seq_type, seq_length)
    except Exception as e:
        print(f"Error computing weighted frequency data: {e}", file=sys.stderr)
        sys.exit(4)

    header = '\t'.join(seq_type.letters())
    try:
        with open(filename, 'w') as f:
            f.write(header + '\n')
            for row in array:
                f.write('\t'.join(f'{v:.6f}' for v in row) + '\n')
    except OSError as e:
        logging.error(f"Error writing to '{filename}': {e}")
        sys.exit(1)


def main():
    ap = setup_argparse()
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="phyinc: %(levelname)s: %(message)s",
    )

    try:
        tree_file = validate_path(args.tree_filename)
        seq_file = validate_path(args.seq_filename)
    except FileNotFoundError as e:
        logging.error(e)
        sys.exit(1)

    try:
        outfilename, graphics_format = output_info(args)
    except ValueError as e:
        logging.error(e)
        sys.exit(1)

    try:
        tree = Phylo.read(tree_file, "newick")
    except IOError as e:
        logging.error(e)
        sys.exit(5)

    try:
        seq_dict, seq_type, seq_length, guessed_color = io.read_sequences(seq_file, "fasta", args)
    except IOError as e:
        logging.error(e)
        sys.exit(3)
    except KeyError as e:
        logging.error(e)
        sys.exit(2)
    except Exception as e:
        logging.error(e)
        sys.exit(4)

    logging.info(f"Alignment width:  {seq_length}")
    logging.info(f"Sequence type:    {seq_type}")
    logging.info(f"Number of leaves: {tree.count_terminals()}")

    try:
        io.check_accession_consistency(seq_dict, tree, args.ignore_coords)
    except ValueError as e:
        logging.error(e)
        sys.exit(4)

    color_scheme = decide_color_scheme(args, guessed_color)

    method = args.weighting

    if args.export:
        export_weighted_data(tree, seq_dict, seq_type, seq_length, args.export, method)
        sys.exit(0)

    try:
        executable = Path(sys.argv[0]).name # Get rid of full path
        if args.command_line:
            command_line = executable + ' ' + ' '.join(sys.argv[1:])
        else:
            command_line = ''
        fig = create_logo(seq_dict, tree, seq_type.name, color_scheme=color_scheme, title=args.title, method=method, footnote=command_line)
        fig.savefig(outfilename, format=graphics_format)
        plt.close(fig)
    except Exception as e:
        logging.error(e)
        sys.exit(2)


if __name__ == "__main__":
    main()
