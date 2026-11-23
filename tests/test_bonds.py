from pymol import cmd

from ihm_fixtures import one_group_file
from pmihm.loader import load_ihm


def bonded(sel1, sel2):
    """True if the single atoms sel1 and sel2 are directly bonded."""
    return cmd.count_atoms('(%s) and neighbor (%s)' % (sel1, sel2)) == 1


def partners(sel):
    out = []
    cmd.iterate('neighbor (%s)' % sel, 'out.append((chain, resi, name))',
                space={'out': out})
    return sorted(out)


def test_consecutive_beads_in_chain_are_bonded(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert bonded('chain A and resi 5 and tt. BEAD',
                  'chain A and resi 6 and tt. BEAD')
    assert bonded('chain A and resi 6 and tt. BEAD',
                  'chain A and resi 11 and tt. BEAD')


def test_beads_not_bonded_across_sequence_gap(tmp_path):
    # B 1-10 and B 12-20: residue 11 is missing
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert partners('chain B and resi 1') == []


def test_atomic_residues_loaded_with_identity(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    got = []
    cmd.iterate('m and chain A and resi 2 and not tt. BEAD',
                'got.append((resn, name, elem))', space={'got': got})
    assert sorted(got) == [('CYS', 'C', 'C'), ('CYS', 'CA', 'C'),
                           ('CYS', 'N', 'N'), ('CYS', 'O', 'O')]


def test_atomic_part_is_bonded_by_distance(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert partners('chain A and resi 1 and name CA') == [
        ('A', '1', 'C'), ('A', '1', 'N')]


def test_last_atomic_ca_bonded_to_following_bead(tmp_path):
    # Residue 4 is atomic, residue 5 is the first bead; the chain should
    # stay connected across the atomic/coarse boundary.
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert partners('chain A and resi 5 and tt. BEAD') == [
        ('A', '4', 'CA'), ('A', '6', 'CA')]


def test_bead_radius_does_not_create_bonds_to_atoms(tmp_path):
    # Bead A5 (r=2.0) at x=16 sits 3.6 A from residue 4's O; distance-based
    # bonding with bead radii would bond them.
    load_ihm(one_group_file(tmp_path / 'm.cif'))
    assert partners('chain A and resi 4 and name O') == [('A', '4', 'C')]


def test_overlapping_beads_bond_to_the_bead_that_follows(tmp_path):
    # multi-scale: 1-10, 11-20 and 1-20 in one model; 1-10 -> 11-20 only
    import ihm_fixtures as F
    load_ihm(F.spheres_file(tmp_path / 'o.cif', [(1, 10), (11, 20), (1, 20)]))
    assert partners('o and resi 1 and custom 1-10') == [('A', '11', 'CA')]
    assert partners('o and resi 1 and custom 1-20') == []


def test_nucleic_acid_bead_anchored_to_c4prime_not_p(tmp_path):
    import ihm_fixtures as F
    load_ihm(F.rna_file(tmp_path / 'r.cif'))
    assert partners('r and tt. BEAD') == [('C', '1', "C4'")]
