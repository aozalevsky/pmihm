import pytest
from pymol import cmd

from ihm_fixtures import one_group_file
from pmihm.loader import load_ihm


def atoms(sel, *fields):
    out = []
    cmd.iterate(sel, 'out.append((%s,))' % ','.join(fields),
                space={'out': out})
    return out


def test_bead_becomes_pseudo_atom_with_file_radius(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    got = atoms('m and chain A and resi 6',
                'name', 'resn', 'vdw', 'custom', 'text_type')
    assert got == [('CA', 'BEA', 4.0, '6-10', 'BEAD')]
    x, y, z = cmd.get_coords('m and chain A and resi 6')[0]
    assert (x, y, z) == pytest.approx((22.0, 0.0, 0.0))


def test_single_residue_bead_keeps_residue_name(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert atoms('m and chain A and resi 5', 'resn', 'custom') == [
        ('PHE', '5-5')]


def test_bead_rmsf_goes_to_b_factor(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert atoms('m and chain A and resi 5', 'b') == [(1.5,)]


def test_all_beads_present_and_selectable(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert cmd.count_atoms('m and tt. BEAD') == 5
    assert cmd.count_atoms('m and chain B') == 2
