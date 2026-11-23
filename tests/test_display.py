from pymol import cmd

import pmihm
from ihm_fixtures import one_group_file, two_group_file
from pmihm.loader import load_ihm


def colors(sel):
    out = set()
    cmd.iterate(sel, 'out.add(color)', space={'out': out})
    return out


def test_beads_shown_as_spheres_only(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert cmd.count_atoms('tt. BEAD and rep spheres') == 5
    assert cmd.count_atoms('tt. BEAD and (rep lines or rep sticks)') == 0
    assert cmd.count_atoms('not tt. BEAD and rep spheres') == 0


def test_atomic_polymer_shown_as_cartoon(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert cmd.count_atoms('not tt. BEAD and rep cartoon') == 16


def test_one_colour_per_chain_distinct_between_chains(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    a = colors('chain A and tt. BEAD')
    b = colors('chain B')
    assert len(a) == 1 and len(b) == 1 and a != b
    assert colors('chain A and name CA and not tt. BEAD') == a


def test_chain_colours_same_in_every_group(tmp_path):
    load_ihm(two_group_file(tmp_path / 'm.cif'))
    for chain in 'AB':
        assert (colors('m_cluster_1 and chain %s and tt. BEAD' % chain)
                == colors('m_cluster_2 and chain %s and tt. BEAD' % chain))


def test_unbonded_ion_is_visible(tmp_path):
    import ihm_fixtures as F
    load_ihm(F.ion_file(tmp_path / 'i.cif'))
    assert cmd.count_atoms('i and elem ZN and rep nb_spheres') == 1


def test_camera_centred_on_shown_group(tmp_path):
    import ihm_fixtures as F
    load_ihm(F.far_groups_file(tmp_path / 'f.cif'))
    assert abs(cmd.get_view()[14]) < 50     # origin z: near group at 0


def test_zoom_0_keeps_camera(tmp_path):
    pmihm.register()
    import ihm_fixtures as F
    before = cmd.get_view()
    cmd.load(F.far_groups_file(tmp_path / 'f.cif'), format='ihm', zoom=0)
    assert cmd.get_view() == before
