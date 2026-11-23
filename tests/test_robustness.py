"""Inputs that broke the first version (code review, 2026-11-23)."""
import os
import pathlib

import pymol
import pytest
from pymol import cmd

import ihm_fixtures as F
from pmihm.loader import load_ihm


def molecules(enabled_only=0):
    return [n for n in cmd.get_names('objects', enabled_only=enabled_only)
            if cmd.get_type(n) == 'object:molecule']


def values(sel, *fields):
    out = []
    cmd.iterate(sel, 'out.append((%s,))' % ','.join(fields),
                space={'out': out})
    return out


# '?' values (ihm.unknown) -------------------------------------------------

def test_unknown_atom_values_get_pymol_defaults(tmp_path):
    load_ihm(F.unknowns_file(tmp_path / 'u.cif'))
    assert values('u and chain A and name N', 'alt', 'q', 'b') == [
        ('', 1.0, 0.0)]


def test_unknown_rmsf_gives_zero_b(tmp_path):
    load_ihm(F.unknowns_file(tmp_path / 'u.cif'))
    assert values('u and chain A and resi 6 and tt. BEAD', 'b') == [(0.0,)]


def test_unknown_names_are_treated_as_missing(tmp_path):
    # single group: object named after the file, state title not '?'
    load_ihm(F.unknowns_file(tmp_path / 'u.cif'))
    assert molecules() == ['u']
    assert cmd.get_title('u', 1) != '?'


def test_sphere_without_residue_range_is_loaded(tmp_path):
    # non-polymer sphere: seq_id_begin/end are '.'
    load_ihm(F.unknowns_file(tmp_path / 'u.cif'))
    assert values('u and chain C and tt. BEAD', 'resn', 'vdw') == [
        ('HEM', 4.0)]


# object names --------------------------------------------------------------

def test_dots_inside_file_stem_are_kept(tmp_path):
    load_ihm(F.one_group_file(tmp_path / 'm.v2.cif'))
    assert molecules() == [cmd.get_legal_name('m.v2')]


def test_dotfile_gives_visible_names(tmp_path):
    load_ihm(F.two_group_file(tmp_path / '.hidden.cif'))
    names = cmd.get_names('all')
    assert names and all(n and not n.startswith('_') for n in names)


def test_loading_same_multi_group_file_twice_gives_new_objects(tmp_path):
    path = F.two_group_file(tmp_path / 'm.cif')
    load_ihm(path)
    load_ihm(path)
    assert sorted(molecules()) == ['m_2_cluster_1', 'm_2_cluster_2',
                                   'm_cluster_1', 'm_cluster_2']
    assert cmd.count_states('m_2_cluster_2') == 2


def test_existing_object_is_left_untouched(tmp_path):
    load_ihm(F.two_group_file(tmp_path / 'two.cif'), object='obj')
    load_ihm(F.one_group_file(tmp_path / 'one.cif'), object='obj_cluster_1')
    assert cmd.count_states('obj_cluster_1') == 2
    assert cmd.count_states('obj_cluster_1_2') == 1


def test_reserved_name_twice(tmp_path):
    path = F.one_group_file(tmp_path / 'm.cif')
    load_ihm(path, object='x')
    load_ihm(path, object='x')
    assert sorted(molecules()) == ['x_', 'x__2']
    assert cmd.count_atoms('x__2 and tt. BEAD and rep spheres') == 5


# file names ----------------------------------------------------------------

def test_pathlib_filename(tmp_path):
    load_ihm(pathlib.Path(F.one_group_file(tmp_path / 'm.cif')))
    assert cmd.count_atoms('m and tt. BEAD') == 5


def test_tilde_filename(tmp_path, monkeypatch):
    monkeypatch.setenv('HOME', str(tmp_path))
    F.one_group_file(tmp_path / 'm.cif')
    load_ihm('~/m.cif')
    assert cmd.count_atoms('m and tt. BEAD') == 5


# alt ids -------------------------------------------------------------------

def test_multi_character_alt_id_in_multi_model_group(tmp_path):
    # PyMOL keeps one character of alt; matching models must still work
    load_ihm(F.alt_id_file(tmp_path / 'a.cif'))
    assert cmd.count_states('a') == 2
    assert cmd.count_atoms('a and not tt. BEAD', state=2) == 16


@pytest.mark.parametrize('given', ['???', '()', '   ', '___'])
def test_object_name_without_legal_characters(tmp_path, given):
    load_ihm(F.two_group_file(tmp_path / 'm.cif'), object=given)
    names = cmd.get_names('all')
    assert len(names) == 3
    assert all(n and not n.startswith('_') for n in names)


def test_file_name_without_legal_characters(tmp_path):
    load_ihm(F.two_group_file(tmp_path / '???.cif'))
    names = cmd.get_names('all')
    assert all(n and not n.startswith('_') for n in names)


def test_empty_group_names_fall_back_to_numbers(tmp_path):
    load_ihm(F.empty_group_name_file(tmp_path / 'e.cif'))
    assert molecules() == ['e_group1', 'e_group2']


def test_non_utf8_bytes_do_not_stop_loading(tmp_path):
    path = F.one_group_file(tmp_path / 'm.cif')
    data = open(path, 'rb').read()
    open(path, 'wb').write(b'# caf\xe9 (latin-1)\n' + data)
    load_ihm(path)
    assert cmd.count_atoms('m and tt. BEAD') == 5
