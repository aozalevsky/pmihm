"""Load IHM (integrative/hybrid) mmCIF/BinaryCIF files into PyMOL.

Parsing is done by python-ihm; this module only maps python-ihm objects
onto PyMOL atoms.
"""

import gzip
import os
import re
import shutil
import urllib.request

import ihm
import ihm.reader
import pymol
from chempy import Atom, Bond, models
from pymol import cmd, util

BEAD_TYPE = 'BEAD'      # text_type of every bead; select with "tt. BEAD"
BEAD_RESN = 'BEA'       # residue name of beads covering several residues
BEAD_NAME = 'CA'
# Atoms of an atomic residue that a neighbouring bead is bonded to,
# in order of preference.
ANCHOR_NAMES = ('CA', "C4'", 'P')
# Stripped from file names to make object names.
FILE_EXTENSIONS = ('.gz', '.cif', '.bcif', '.mmcif')

PDB_IHM_URL = 'https://pdb-ihm.org/cif/'
DOWNLOAD_TIMEOUT = 60   # seconds without data before a download fails


def _known(value, default=None):
    """value, or default if it is missing ('.') or unknown ('?')."""
    if value is None or value is ihm.unknown:
        return default
    return value


def _read_system(filename):
    """First data block of an mmCIF or BinaryCIF file, optionally gzipped."""
    name = filename.lower()
    if name.endswith('.gz'):
        opener, name = gzip.open, name[:-3]
    else:
        opener = open
    binary = name.endswith('.bcif')
    try:
        if binary:
            with opener(filename, 'rb') as fh:
                return ihm.reader.read(fh, format='BCIF')[0]
        with opener(filename, 'rt', encoding='utf-8',
                    errors='replace') as fh:
            return ihm.reader.read(fh)[0]
    except OSError as e:
        raise pymol.CmdException('cannot read %s: %s' % (filename, e))


def _residue_name(asym, seq_id):
    seq = asym.entity.sequence if asym.entity else []
    if not seq:
        return 'UNK'
    if seq_id is not None and 1 <= seq_id <= len(seq):
        return seq[seq_id - 1].id
    return seq[0].id


def _sphere_range(sphere):
    """(begin, end) residue range of a sphere; (0, 0) if it has none
       (a sphere of a non-polymer)."""
    begin, end = sphere.seq_id_range
    return _known(begin, 0), _known(end, 0)


def _bead_atom(sphere):
    asym = sphere.asym_unit
    begin, end = _sphere_range(sphere)
    a = Atom()
    a.name = BEAD_NAME
    a.symbol = 'C'
    a.resn = (_residue_name(asym, begin or None) if begin == end
              else BEAD_RESN)
    a.resi = str(begin or 1)
    a.chain = a.segi = asym.id
    a.hetatm = 0
    a.coord = [sphere.x, sphere.y, sphere.z]
    a.vdw = sphere.radius
    a.b = _known(sphere.rmsf, 0.0)
    a.custom = '%d-%d' % (begin, end) if begin else ''
    a.text_type = BEAD_TYPE
    return a


def _atom(atom):
    asym = atom.asym_unit
    seq_id = _known(atom.seq_id)
    a = Atom()
    a.name = atom.atom_id
    a.symbol = _known(atom.type_symbol, '')
    a.resn = _residue_name(asym, seq_id)
    a.resi = str(seq_id or 1)
    a.chain = a.segi = asym.id
    a.alt = _known(atom.alt_id, '')[:1]     # PyMOL keeps one character
    a.hetatm = 1 if atom.het else 0
    a.coord = [atom.x, atom.y, atom.z]
    a.b = _known(atom.biso, 0.0)
    a.q = _known(atom.occupancy, 1.0)
    return a


def _atom_key(atom):
    return (atom.asym_unit.id, _known(atom.seq_id) or 1, atom.atom_id,
            _known(atom.alt_id, '')[:1])


def _chempy_atom_key(a):
    """The _atom_key of the IHM atom that chempy atom a was made from."""
    return (a.chain, a.resi_number, a.name, a.alt)


def _sphere_key(sphere):
    return (sphere.asym_unit.id,) + _sphere_range(sphere)


def _add_bond(m, i, j):
    b = Bond()
    b.index = [i, j]
    b.order = 1
    m.add_bond(b)


def _bonded_atomic_part(ihm_atoms):
    """chempy model of the atomic part, bonded by PyMOL from distances.

       Bonding is done without the beads present: PyMOL's distance-based
       bonding takes atom radii into account and would bond beads to
       everything around them."""
    m = models.Indexed()
    for atom in ihm_atoms:
        m.add_atom(_atom(atom))
    if not ihm_atoms:
        return m
    tmp = cmd.get_unused_name('_ihm_tmp')
    cmd.load_model(m, tmp, zoom=0)
    try:
        cmd.rebond(tmp)
        return cmd.get_model(tmp)
    finally:
        cmd.delete(tmp)


def _anchors(m):
    """(chain, residue number) -> index of the atom a neighbouring bead is
       bonded to, the first of ANCHOR_NAMES present in that residue."""
    best = {}
    for i, a in enumerate(m.atom):
        if a.name in ANCHOR_NAMES:
            key = (a.chain, a.resi_number)
            rank = ANCHOR_NAMES.index(a.name)
            if key not in best or rank < best[key][0]:
                best[key] = (rank, i)
    return {key: i for key, (rank, i) in best.items()}


def _build_model(atoms, spheres):
    """chempy model with atoms, beads and all bonds for one IHM model, and
       the key (_atom_key or _sphere_key) of each of its atoms."""
    m = _bonded_atomic_part(atoms)
    keys = [_chempy_atom_key(a) for a in m.atom]
    anchors = _anchors(m)

    # asym id -> first residue -> indices of the beads starting there
    starts = {}
    beads = []
    for sphere in spheres:
        idx = len(m.atom)
        m.add_atom(_bead_atom(sphere))
        key = _sphere_key(sphere)
        keys.append(key)
        asym_id, begin, end = key
        if begin:
            beads.append((idx, asym_id, begin, end))
            starts.setdefault(asym_id, {}).setdefault(begin, []).append(idx)
    for idx, asym_id, begin, end in beads:
        # the next bead(s) in the chain: those that start right after
        for nxt in starts[asym_id].get(end + 1, ()):
            _add_bond(m, idx, nxt)
        for seq_id in (begin - 1, end + 1):
            if (asym_id, seq_id) in anchors:
                _add_bond(m, anchors[(asym_id, seq_id)], idx)
    return m, keys


def _layout(atoms, spheres):
    """What has to be identical for two models to share one PyMOL object
       as states: the same atoms and the same beads with the same radii."""
    return (sorted(_atom_key(a) for a in atoms),
            sorted(_sphere_key(s) + (s.radius,) for s in spheres))


def _has_unique_keys(atoms, spheres):
    keys = [_atom_key(a) for a in atoms] + [_sphere_key(s) for s in spheres]
    return len(set(keys)) == len(keys)


def _coords(keys, atoms, spheres):
    """Coordinates of a model's atoms and spheres, in the order of keys."""
    xyz = {_atom_key(a): (a.x, a.y, a.z) for a in atoms}
    xyz.update((_sphere_key(s), (s.x, s.y, s.z)) for s in spheres)
    return [xyz[k] for k in keys]


def _load_group(models_, name, discrete):
    """Load one model group as object `name` (which must not exist yet),
       one state per model."""
    data = [(model, list(model.get_atoms()), list(model.get_spheres()))
            for model in models_]
    layout = _layout(*data[0][1:])
    shared = (not discrete
              and all(_has_unique_keys(a, s) for _, a, s in data)
              and all(_layout(a, s) == layout for _, a, s in data[1:]))
    if shared:
        template, keys = _build_model(*data[0][1:])
        for state, (model, atoms, spheres) in enumerate(data, start=1):
            if state > 1:
                for a, xyz in zip(template.atom,
                                  _coords(keys, atoms, spheres)):
                    a.coord = list(xyz)
            cmd.load_model(template, name, state=state, zoom=0)
    else:
        for state, (model, atoms, spheres) in enumerate(data, start=1):
            m, _ = _build_model(atoms, spheres)
            cmd.load_model(m, name, state=state, discrete=1, zoom=0)
    for state, (model, _, _) in enumerate(data, start=1):
        title = _known(model.name)
        if title:
            cmd.set_title(name, state, str(title))


def _has_coordinates(model):
    return (any(True for _ in model.get_spheres())
            or any(True for _ in model.get_atoms()))


def _model_groups(system):
    """(label parts, models) for every model group that has coordinates."""
    states = [s for sg in system.state_groups for s in sg]
    out = []
    for si, state in enumerate(states, start=1):
        for gi, group in enumerate(state, start=1):
            models_ = [m for m in group if _has_coordinates(m)]
            if not models_:
                continue
            parts = []
            if len(states) > 1:
                parts.append(str(_known(state.name) or 'state%d' % si))
            parts.append(str(_known(group.name) or 'group%d' % gi))
            out.append((parts, models_))
    return out


def _chain_colors(system):
    """One colour index per asym, in file order, PyMOL's chain palette."""
    cycle = getattr(util, '_color_cycle', None) or [
        cmd.get_color_index(c) for c in
        ('green', 'cyan', 'magenta', 'yellow', 'salmon', 'slate', 'orange')]
    return {asym.id: cycle[i % len(cycle)]
            for i, asym in enumerate(system.asym_units)}


def _style(name, chain_colors):
    beads = '%s and tt. %s' % (name, BEAD_TYPE)
    atoms = '%s and not tt. %s' % (name, BEAD_TYPE)
    cmd.hide('everything', name)
    cmd.show('spheres', beads)
    cmd.show('cartoon', atoms + ' and not hetatm')
    cmd.show('sticks', atoms + ' and hetatm and not solvent')
    # ions and other atoms without bonds draw nothing as sticks
    cmd.show('nb_spheres',
             atoms + ' and hetatm and not solvent and not (neighbor all)')
    # alter, not "color ... chain X": chain selections ignore case
    cmd.alter(name, 'color = colors.get(chain, color)',
              space={'colors': chain_colors})
    if cmd.count_atoms(atoms):
        util.cnc(atoms)
        cmd.dss(atoms)
    cmd.rebuild(name)


def _object_name(filename):
    """File name without directory and IHM file extensions, never empty
       and not starting with '.' or '_' (PyMOL hides such objects)."""
    base = os.path.basename(filename)
    while base.lower().endswith(FILE_EXTENSIONS):
        base = base[:base.rindex('.')]
    return base.lstrip('._') or 'ihm'


class _Names:
    """Hands out legal PyMOL names that are not in use yet."""

    def __init__(self):
        self.taken = set(cmd.get_names('all'))

    def new(self, name):
        name = cmd.get_legal_name(name).lstrip('._')
        name = cmd.get_legal_name(name or 'ihm')
        out, n = name, 1
        while out in self.taken:
            n += 1
            out = '%s_%d' % (name, n)
        self.taken.add(out)
        return out


def load_ihm(filename, object='', discrete=0, zoom=-1, quiet=1):
    """Load an IHM mmCIF/BinaryCIF file.

       Each model group becomes one object with one state per model; with
       several groups the objects are put in a PyMOL group named after the
       file and only the first is shown. Existing objects are never
       changed: names in use get a _2, _3, ... suffix.

       discrete = 1: one set of atoms per state, so each state keeps its
       own B-factors and occupancies (otherwise models that share a layout
       share atoms, and those values come from the first model).

       Returns the name of the object, or of the group."""
    filename = os.path.expanduser(os.fspath(filename))
    system = _read_system(filename)
    groups = _model_groups(system)
    if not groups:
        raise pymol.CmdException(
            '%s: no atoms or spheres in any model (Gaussian objects are '
            'not read by python-ihm)' % filename)
    discrete = int(discrete) > 0
    colors = _chain_colors(system)
    names = _Names()
    base = names.new(str(object).strip() or _object_name(filename))
    if len(groups) == 1:
        objects = [base]
    else:
        objects = [names.new('_'.join([base] + parts))
                   for parts, _ in groups]
    for obj, (_, models_) in zip(objects, groups):
        _load_group(models_, obj, discrete)
        _style(obj, colors)
    if len(groups) > 1:
        cmd.group(base, ' '.join(objects))
        for obj in objects[1:]:
            cmd.disable(obj)
    zoom = int(zoom)
    if zoom > 0 or (zoom < 0 and cmd.get_setting_int('auto_zoom')):
        cmd.zoom(objects[0])
    if not int(quiet):
        print(' load_ihm: %s loaded as %s' % (filename, ', '.join(objects)))
    return base


def normalize_id(code):
    """Canonical PDB-IHM entry id: PDBDEV_######## or a 4-character PDB id.
       A bare number is taken as a PDBDEV number."""
    code = str(code).strip()
    m = re.fullmatch(r'(?:PDBDEV_)?(\d{1,8})', code, re.IGNORECASE)
    if m:
        return 'PDBDEV_%08d' % int(m.group(1))
    if re.fullmatch(r'[1-9][A-Za-z0-9]{3}', code):
        return code.upper()
    raise ValueError('not a PDB-IHM entry id: %r' % code)


def _download(url, dest, quiet):
    if not quiet:
        print(' fetch_ihm: downloading %s' % url)
    part = dest + '.part'
    try:
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT) as r, \
                open(part, 'wb') as fh:
            shutil.copyfileobj(r, fh)
        os.replace(part, dest)
    except OSError as e:    # includes urllib.error.URLError
        raise pymol.CmdException('cannot download %s: %s' % (url, e))
    finally:
        if os.path.exists(part):
            os.remove(part)


def fetch_ihm(code, name='', path='', discrete=0, zoom=-1, quiet=1,
              base_url=''):
    """Download a PDB-IHM entry into path (PyMOL's fetch_path by default)
       unless it is already there, and load it with load_ihm."""
    try:
        code = normalize_id(code)
    except ValueError as e:
        raise pymol.CmdException(str(e))
    quiet = int(quiet)
    path = os.path.expanduser(path or cmd.get('fetch_path') or '.')
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        raise pymol.CmdException('cannot use %s for downloads: %s'
                                 % (path, e))
    dest = os.path.join(path, code + '.cif')
    if not os.path.exists(dest):
        url = (base_url or PDB_IHM_URL).rstrip('/') + '/' + code + '.cif'
        _download(url, dest, quiet)
    return load_ihm(dest, name or code, discrete=discrete, zoom=zoom,
                    quiet=quiet)
