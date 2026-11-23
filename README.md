# pmihm: integrative (IHM) models in PyMOL

PyMOL reads the ordinary atomic part (`atom_site`) of an IHM/IHMCIF file, but
not the coarse-grained beads (`ihm_sphere_obj_site`) that most integrative
models in [PDB-IHM](https://pdb-ihm.org) are made of. pmihm reads the whole
file with [python-ihm](https://github.com/ihmwg/python-ihm) and loads atoms
and beads together.

## Requirements

- PyMOL 2 or later. Tested with open-source PyMOL 3.1.0 on Python 3.14.
- [python-ihm](https://pypi.org/project/ihm/) 2.11 or later (`ihm` on PyPI), installed into
  the Python that PyMOL runs on. BinaryCIF files also need `msgpack`.

## Install

Install the package into the Python that PyMOL uses, then turn it on in
PyMOL.

1. Install. With a system PyMOL (for example Fedora's), run
   `python3 -m pip install --user pmihm`, or `pip install pmihm[bcif]` to
   read BinaryCIF too. From a checkout of this repository, use
   `pip install -e .` instead.
2. Turn it on by adding this line to `~/.pymolrc`:

   ```
   /import pmihm; pmihm.register()
   ```

Or, without pip, use PyMOL's Plugin Manager. In PyMOL go to Plugin → Plugin
Manager → Install New Plugin and choose `pmihm/__init__.py` from a checkout.
This also adds **IHM: Open file...** and **IHM: Fetch from PDB-IHM...** to the
Plugin menu. python-ihm still has to be installed separately.

## Usage

```
fetch PDBDEV_00000001                 # PDBDEV ids go to PDB-IHM automatically
fetch 9A8R, type=ihm                  # other ids: ask for type=ihm
fetch_ihm 9A8R                        # the same; also takes PDBDEV_1 or 1
load_ihm entry.cif                    # .cif, .cif.gz, .bcif, .bcif.gz
load_ihm entry.cif, myname            # choose the object name
load entry.cif, format=ihm            # through PyMOL's own load command
load_ihm entry.cif, discrete=1        # keep per-model B-factors (see below)
```

Downloads go to PyMOL's `fetch_path`, and a file already there is reused.
Loading never changes existing objects: if a name is taken, the new object
gets `_2`, `_3`, ... appended.

## How IHM data maps onto PyMOL

| IHM | PyMOL |
|---|---|
| model group (`ihm_model_group`) | One object. If the file has several groups, they go into a PyMOL group named after the file, and only the first group is shown. |
| model in a group | One state, titled with the model name. Models with the same atoms and beads share atoms, as in PyMOL's own multi-model loading, so B-factors and occupancies come from the first model. Use `discrete=1` to keep each model's own values. Models with different beads or radii are always loaded discretely. |
| several states (`ihm_multi_state_modeling`) | The state name becomes part of the object name, e.g. `entry_open_cluster_1`. |
| bead (`ihm_sphere_obj_site`) | A pseudo-atom with `name CA`, `vdw` = sphere radius, and `b` = RMSF (0 if not given). |
| residue range of a bead | `resi` = first residue. `custom` holds the full range, e.g. `"6-10"`. |
| residue name of a bead | The real residue name for a 1-residue bead, `BEA` for a multi-residue bead. |
| chain | `chain` and `segi` = `label_asym_id`. All numbering is label numbering (`label_seq_id`), for beads and atoms alike. |
| chain connectivity | Each bead is bonded to the bead that starts on the next residue, so there's no bond across a sequence gap. A bead next to an atomic residue is bonded to that residue's CA (or C4', then P). Atoms are bonded by PyMOL from distances. |

Default display:
- beads as spheres
- protein and nucleic acid as cartoon, with secondary structure from `dss`
- ligands as sticks, ions as non-bonded spheres

Everything is coloured by chain, and a chain has the same colour in every
group's object.

Useful selections and commands:

```
select beads, tt. BEAD                        # every bead
show sticks, tt. BEAD                         # draw the bead chain
label tt. BEAD, custom                        # show each bead's residue range
spectrum b, blue_red, tt. BEAD                # colour beads by RMSF
iterate tt. BEAD and chain A, print(custom, vdw)
```

## Limitations

- Gaussian objects (`ihm_gaussian_obj_site`) are not read by python-ihm, so
  they are not loaded. A file with only Gaussians gives an error.
- Datasets, maps, localization densities, restraints and starting models are
  not loaded.
- `resi` is the first residue of a bead, so `resi 7` does not select a bead
  covering residues 6-10. Use the bead's `custom` range instead.
- The Plugin-menu dialogs have not been tested in a running GUI yet.

## Development

```
pip install -e .[test]
python3 -m pytest tests
```

The tests build small IHM files with python-ihm and run PyMOL headless.
