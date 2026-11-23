import pymol
import pytest
from pymol import cmd

from ihm_fixtures import one_group_file, two_group_file, two_state_file
from pmihm.loader import load_ihm


def molecules(enabled_only=0):
    return [n for n in cmd.get_names('objects', enabled_only=enabled_only)
            if cmd.get_type(n) == 'object:molecule']


def test_single_group_loads_as_object_named_after_file(tmp_path):
    load_ihm(one_group_file(tmp_path / 'entry.cif'))
    assert molecules() == ['entry']


def test_object_name_can_be_given(tmp_path):
    load_ihm(one_group_file(tmp_path / 'entry.cif'), object='mine')
    assert molecules() == ['mine']


def test_each_model_group_becomes_object_in_file_group(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert molecules() == ['m_cluster_1', 'm_cluster_2']
    assert cmd.get_type('m') == 'object:group'
    assert sorted(cmd.get_object_list('m')) == ['m_cluster_1', 'm_cluster_2']


def test_only_first_group_is_shown(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert molecules(enabled_only=1) == ['m_cluster_1']


def test_models_with_same_layout_become_states(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert cmd.count_states('m_cluster_1') == 2
    bead = 'm_cluster_1 and chain A and resi 6 and tt. BEAD'
    assert cmd.count_atoms(bead) == 1
    assert cmd.get_coords(bead, state=1)[0][2] == pytest.approx(0.0)
    assert cmd.get_coords(bead, state=2)[0][2] == pytest.approx(10.0)


def test_state_titles_are_model_names(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert cmd.get_title('m_cluster_1', 1) == 'best'
    assert cmd.get_title('m_cluster_1', 2) == 'second'


def test_models_with_different_layouts_keep_own_beads(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert cmd.count_states('m_cluster_2') == 2
    radii = {}
    for state in (1, 2):
        got = []
        cmd.iterate_state(state, 'm_cluster_2 and chain B',
                          'got.append((custom, vdw))', space={'got': got})
        radii[state] = sorted(got)
    assert radii[1] == [('1-20', 8.0)]
    assert [r[0] for r in radii[2]] == ['1-10', '12-20']
    assert [r[1] for r in radii[2]] == pytest.approx([4.5, 4.2])


def test_states_of_multi_state_file_get_own_objects(tmp_path):
    load_ihm(two_state_file(tmp_path / 'm.cif'))
    assert molecules() == ['m_open_g', 'm_closed_g']


def test_states_do_not_duplicate_bonds(tmp_path):
    load_ihm(one_group_file(tmp_path / 'one.cif'))
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    assert (len(cmd.get_model('m_cluster_1', state=1).bond)
            == len(cmd.get_model('one').bond))


def test_file_without_coordinates_is_an_error(tmp_path):
    from ihm_fixtures import no_coordinates_file
    with pytest.raises(pymol.CmdException, match='no atoms or spheres'):
        load_ihm(no_coordinates_file(tmp_path / 'm.cif'))
    assert cmd.get_names('all') == []


def test_discrete_keeps_per_model_b_factors(tmp_path):
    from ihm_fixtures import rmsf_per_model_file
    load_ihm(rmsf_per_model_file(tmp_path / 'r.cif'), discrete=1)
    b = {}
    for state in (1, 2):
        got = []
        cmd.iterate_state(state, 'r and chain A and resi 5 and tt. BEAD',
                          'got.append(b)', space={'got': got})
        b[state] = got
    assert b == {1: [1.5], 2: [7.0]}
