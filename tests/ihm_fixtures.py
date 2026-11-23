"""Build small IHM files with python-ihm for the tests.

The layout of the default system (sequence ACDEFGHIKLMNPQRSTVWY, 20 residues,
two copies, asym A and B):

  A  residues 1-4   atomic (backbone N, CA, C, O)
     residue  5     sphere r=2.0    (single residue, PHE)
     residues 6-10  sphere r=4.0
     residues 11-20 sphere r=5.0
  B  residues 1-10  sphere r=4.5
     residues 12-20 sphere r=4.2    (gap at residue 11)
"""

import gzip

import ihm
import ihm.dumper
import ihm.model
import ihm.protocol
import ihm.representation as rep

SEQ = 'ACDEFGHIKLMNPQRSTVWY'


class Builder:
    def __init__(self):
        s = self.system = ihm.System(title='test')
        self.entity = ihm.Entity(SEQ, description='Prot')
        s.entities.append(self.entity)
        self.A = ihm.AsymUnit(self.entity, details='copy A')
        self.B = ihm.AsymUnit(self.entity, details='copy B')
        s.asym_units.extend([self.A, self.B])
        self.assembly = ihm.Assembly([self.A, self.B], name='all')
        self.representation = rep.Representation([
            rep.AtomicSegment(self.A(1, 4), rigid=True),
            rep.FeatureSegment(self.A(5, 20), rigid=False,
                               primitive='sphere', count=3),
            rep.FeatureSegment(self.B(1, 20), rigid=False,
                               primitive='sphere', count=2)])
        self.protocol = ihm.protocol.Protocol(name='p')
        self.protocol.steps.append(ihm.protocol.Step(
            assembly=self.assembly, dataset_group=None,
            method='Monte Carlo', name='s', num_models_begin=0,
            num_models_end=1))

    def model(self, name=None):
        return ihm.model.Model(assembly=self.assembly,
                               protocol=self.protocol,
                               representation=self.representation,
                               name=name)

    def add_backbone(self, model, dz=0.0):
        """Residues 1-4 of A, backbone only, 3.8 A apart along x."""
        for i in range(1, 5):
            x = 3.8 * (i - 1)
            for name, el, (ax, ay) in (('N', 'N', (0.0, 0.0)),
                                       ('CA', 'C', (1.46, 0.0)),
                                       ('C', 'C', (2.0, 1.42)),
                                       ('O', 'O', (2.0, 2.65))):
                model.add_atom(ihm.model.Atom(
                    asym_unit=self.A, seq_id=i, atom_id=name,
                    type_symbol=el, x=x + ax, y=ay, z=dz))

    def add_default_beads(self, model, dz=0.0):
        A, B = self.A, self.B
        sph = ihm.model.Sphere
        model.add_sphere(sph(A, (5, 5), 16.0, 0.0, dz, 2.0, rmsf=1.5))
        model.add_sphere(sph(A, (6, 10), 22.0, 0.0, dz, 4.0))
        model.add_sphere(sph(A, (11, 20), 31.0, 0.0, dz, 5.0))
        model.add_sphere(sph(B, (1, 10), 0.0, 20.0, dz, 4.5))
        model.add_sphere(sph(B, (12, 20), 9.0, 20.0, dz, 4.2))

    def default_model(self, name=None, dz=0.0):
        m = self.model(name)
        self.add_backbone(m, dz)
        self.add_default_beads(m, dz)
        return m

    def set_groups(self, *groups):
        """groups: (name, [models]) pairs, all in one state."""
        state = ihm.model.State([ihm.model.ModelGroup(models, name=name)
                                 for name, models in groups])
        self.system.state_groups.append(ihm.model.StateGroup([state]))

    def write(self, path, fmt='mmCIF', check=True):
        path = str(path)
        if fmt == 'BCIF':
            with open(path, 'wb') as fh:
                ihm.dumper.write(fh, [self.system], format='BCIF',
                                 check=check)
        elif path.endswith('.gz'):
            with gzip.open(path, 'wt') as fh:
                ihm.dumper.write(fh, [self.system], check=check)
        else:
            with open(path, 'w') as fh:
                ihm.dumper.write(fh, [self.system], check=check)
        return path


def one_group_file(path, fmt='mmCIF'):
    """One group, one model with the default layout."""
    b = Builder()
    b.set_groups(('cluster 1', [b.default_model('best')]))
    return b.write(path, fmt)


def two_group_file(path):
    """cluster 1: two models with the default layout, the second shifted
       +10 A in z.  cluster 2: two models with different bead layouts
       (chain B as one r=8 bead, then the default layout)."""
    b = Builder()
    m1 = b.default_model('best')
    m2 = b.default_model('second', dz=10.0)
    m3 = b.model('coarse')
    b.add_backbone(m3)
    sph = ihm.model.Sphere
    m3.add_sphere(sph(b.A, (5, 5), 16.0, 0.0, 0.0, 2.0))
    m3.add_sphere(sph(b.A, (6, 10), 22.0, 0.0, 0.0, 4.0))
    m3.add_sphere(sph(b.A, (11, 20), 31.0, 0.0, 0.0, 5.0))
    m3.add_sphere(sph(b.B, (1, 20), 5.0, 20.0, 0.0, 8.0))
    m4 = b.default_model('fine')
    b.set_groups(('cluster 1', [m1, m2]), ('cluster 2', [m3, m4]))
    return b.write(path)


def two_state_file(path):
    """Two states ('open', 'closed'), each with one group called 'g'."""
    b = Builder()
    states = []
    for sname, dz in (('open', 0.0), ('closed', 5.0)):
        mg = ihm.model.ModelGroup([b.default_model(dz=dz)], name='g')
        states.append(ihm.model.State([mg], name=sname))
    b.system.state_groups.append(ihm.model.StateGroup(states))
    return b.write(path)


def no_coordinates_file(path):
    """A model group whose only model has no atoms or spheres."""
    b = Builder()
    b.set_groups(('empty', [b.model('nothing')]))
    return b.write(path, check=False)


def unknowns_file(path):
    """'?' (ihm.unknown) wherever python-ihm allows it: model, group and
       state names, atom alt_id/occupancy/B, sphere rmsf; plus a sphere of
       a non-polymer (heme) with no residue range ('.' '.')."""
    b = Builder()
    m = b.model(ihm.unknown)
    b.add_default_beads(m)
    m._spheres[1].rmsf = ihm.unknown          # A 6-10
    m.add_atom(ihm.model.Atom(
        asym_unit=b.A, seq_id=1, atom_id='N', type_symbol='N',
        x=0.0, y=0.0, z=0.0, biso=ihm.unknown, occupancy=ihm.unknown,
        alt_id=ihm.unknown))
    heme = ihm.Entity([ihm.NonPolymerChemComp('HEM')], description='heme')
    b.system.entities.append(heme)
    h = ihm.AsymUnit(heme)
    b.system.asym_units.append(h)
    m.add_sphere(ihm.model.Sphere(h, (None, None), 50.0, 0.0, 0.0, 4.0))
    mg = ihm.model.ModelGroup([m], name=ihm.unknown)
    b.system.state_groups.append(ihm.model.StateGroup(
        [ihm.model.State([mg], name=ihm.unknown)]))
    return b.write(path, check=False)


def ion_file(path):
    """Default model plus one zinc ion (HETATM, no bonds)."""
    b = Builder()
    zn = ihm.Entity([ihm.NonPolymerChemComp('ZN')], description='zinc')
    b.system.entities.append(zn)
    z = ihm.AsymUnit(zn)
    b.system.asym_units.append(z)
    m = b.default_model('m')
    m.add_atom(ihm.model.Atom(asym_unit=z, seq_id=None, atom_id='ZN',
                              type_symbol='ZN', x=5.0, y=5.0, z=5.0,
                              het=True))
    b.set_groups(('g', [m]))
    return b.write(path, check=False)


def spheres_file(path, ranges, asym='A'):
    """One model whose chain `asym` has only spheres covering `ranges`,
       placed 5 A apart along x."""
    b = Builder()
    m = b.model('m')
    for i, r in enumerate(ranges):
        m.add_sphere(ihm.model.Sphere(getattr(b, asym), r, 5.0 * i, 0.0,
                                      0.0, 2.0))
    b.set_groups(('g', [m]))
    return b.write(path, check=False)


def rna_file(path):
    """RNA ACGUA: residue 1 atomic with P listed before C4', residues
       2-5 one bead."""
    b = Builder()
    rna = ihm.Entity('ACGUA', alphabet=ihm.RNAAlphabet, description='rna')
    b.system.entities.append(rna)
    r = ihm.AsymUnit(rna)
    b.system.asym_units.append(r)
    m = b.model('m')
    for name, el, x in (('P', 'P', 0.0), ("C4'", 'C', 4.0)):
        m.add_atom(ihm.model.Atom(asym_unit=r, seq_id=1, atom_id=name,
                                  type_symbol=el, x=x, y=0.0, z=0.0))
    m.add_sphere(ihm.model.Sphere(r, (2, 5), 9.0, 0.0, 0.0, 3.0))
    b.set_groups(('g', [m]))
    return b.write(path, check=False)


def far_groups_file(path):
    """Group 'near' at z=0 and group 'far' at z=500."""
    b = Builder()
    b.set_groups(('near', [b.default_model('a')]),
                 ('far', [b.default_model('b', dz=500.0)]))
    return b.write(path)


def rmsf_per_model_file(path):
    """One group, two models with identical layout; bead A5 has rmsf 1.5
       in the first model and 7.0 in the second."""
    b = Builder()
    m1 = b.default_model('m1')
    m2 = b.default_model('m2')
    m2._spheres[0].rmsf = 7.0
    b.set_groups(('g', [m1, m2]))
    return b.write(path)


def alt_id_file(path, alt='A1'):
    """Two models with identical layout whose atoms all carry alt_id."""
    b = Builder()
    models = [b.default_model('m1'), b.default_model('m2', dz=3.0)]
    for m in models:
        for a in m._atoms:
            a.alt_id = alt
    b.set_groups(('g', models))
    return b.write(path, check=False)


def empty_group_name_file(path):
    """Two groups named '' (written as '')."""
    b = Builder()
    b.set_groups(('', [b.default_model('a')]), ('', [b.default_model('b')]))
    return b.write(path)
