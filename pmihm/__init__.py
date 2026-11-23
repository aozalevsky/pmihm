"""pmihm: integrative/hybrid (IHM) models from PDB-IHM in PyMOL.

Commands:
    load_ihm filename [, object [, discrete [, zoom [, quiet ]]]]
    fetch_ihm code [, name [, path [, discrete [, zoom [, quiet ]]]]]
    load filename, format=ihm
    fetch code, type=ihm       (PDBDEV_ ids work without type=ihm)
"""

import re

__version__ = '0.1.0'


def _is_pdbdev(code):
    return re.fullmatch(r'PDBDEV_\d+', code, re.IGNORECASE) is not None


def _install_fetch_hook():
    """Make PyMOL's fetch (command and cmd.fetch) handle type=ihm.

       PyMOL's fetch has no hook for new types; it hands every code to
       pymol.importing._multifetch, looked up at call time, so wrapping that
       keeps fetch's own handling of async, fetch_path and flushing."""
    from pymol import importing
    original = importing._multifetch
    if getattr(original, '_pmihm', False):
        return

    def _multifetch(code, name, state, finish, discrete, multiplex, zoom,
                    type, path, file, quiet, _self):
        from .loader import fetch_ihm
        codes = str(code).split()
        if str(type).lower() == 'ihm':
            ihm_codes, others = codes, []
        elif not type:
            ihm_codes = [c for c in codes if _is_pdbdev(c)]
            others = [c for c in codes if not _is_pdbdev(c)]
        else:
            ihm_codes, others = [], codes
        r = None
        for c in ihm_codes:
            r = fetch_ihm(c, name, path, discrete=discrete, zoom=zoom,
                          quiet=quiet)
        if others or not codes:
            r = original(' '.join(others) if ihm_codes else code, name,
                         state, finish, discrete, multiplex, zoom, type,
                         path, file, quiet, _self)
        return r

    _multifetch._pmihm = True
    importing._multifetch = _multifetch


def register():
    """Add the commands, the 'ihm' load format and fetch type to PyMOL."""
    from pymol import cmd, importing
    from .loader import load_ihm, fetch_ihm

    cmd.extend('load_ihm', load_ihm)
    cmd.extend('fetch_ihm', fetch_ihm)
    importing.loadfunctions['ihm'] = load_ihm
    _install_fetch_hook()


def _error_box(error):
    from pymol.Qt import QtWidgets
    QtWidgets.QMessageBox.critical(None, 'pmihm', str(error))


def _open_dialog():
    import pymol
    from pymol.Qt import QtWidgets
    from .loader import load_ihm
    filename, _ = QtWidgets.QFileDialog.getOpenFileName(
        None, 'Open IHM file', '',
        'IHM mmCIF (*.cif *.cif.gz *.bcif *.bcif.gz);;All files (*)')
    if filename:
        try:
            load_ihm(filename, quiet=0)
        except (pymol.CmdException, OSError, ValueError) as e:
            _error_box(e)


def _fetch_dialog():
    from pymol import cmd
    from pymol.Qt import QtWidgets
    from .loader import normalize_id
    code, ok = QtWidgets.QInputDialog.getText(
        None, 'Fetch from PDB-IHM', 'Entry id (e.g. PDBDEV_00000001, 9A8R):')
    if not (ok and code.strip()):
        return
    try:
        normalize_id(code)
    except ValueError as e:
        _error_box(e)
        return
    # in the background, so the GUI stays responsive; errors go to the log
    cmd.fetch(code.strip(), type='ihm', async_=1, quiet=0)


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt, QtNotAvailableError
    register()
    try:
        addmenuitemqt('IHM: Open file...', _open_dialog)
        addmenuitemqt('IHM: Fetch from PDB-IHM...', _fetch_dialog)
    except QtNotAvailableError:
        pass  # command-line PyMOL (-c): commands only
