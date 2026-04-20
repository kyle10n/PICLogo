import pytest
from types import SimpleNamespace

from phyinc.colorhelper import decide_color_scheme, taylor


def test_decide_color_scheme_guess_returns_passed_scheme():
    args = SimpleNamespace(color_scheme='guess')
    result = decide_color_scheme(args, 'classic')
    assert result == 'classic'


def test_decide_color_scheme_guess_preserves_hydrophobicity():
    args = SimpleNamespace(color_scheme='guess')
    result = decide_color_scheme(args, 'hydrophobicity')
    assert result == 'hydrophobicity'


def test_decide_color_scheme_monochrome():
    args = SimpleNamespace(color_scheme='monochrome')
    result = decide_color_scheme(args, 'classic')
    assert result == 'gray'


def test_decide_color_scheme_nucleotide():
    args = SimpleNamespace(color_scheme='nucleotide')
    result = decide_color_scheme(args, 'hydrophobicity')
    assert result == 'classic'


def test_decide_color_scheme_hydrophobicity():
    args = SimpleNamespace(color_scheme='hydrophobicity')
    result = decide_color_scheme(args, 'classic')
    assert result == 'hydrophobicity'


def test_decide_color_scheme_chemistry():
    args = SimpleNamespace(color_scheme='chemistry')
    result = decide_color_scheme(args, 'classic')
    assert result == 'chemistry'


def test_decide_color_scheme_charge():
    args = SimpleNamespace(color_scheme='charge')
    result = decide_color_scheme(args, 'classic')
    assert result == 'charge'


def test_decide_color_scheme_taylor():
    args = SimpleNamespace(color_scheme='taylor')
    result = decide_color_scheme(args, 'classic')
    assert result is taylor


def test_decide_color_scheme_unknown_falls_back_to_gray():
    args = SimpleNamespace(color_scheme='nonexistent')
    result = decide_color_scheme(args, 'classic')
    assert result == 'gray'
