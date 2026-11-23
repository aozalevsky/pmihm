import pytest
from pymol import cmd

import pymol
from ihm_fixtures import one_group_file, two_group_file
import pmihm
from pmihm.loader import load_ihm, fetch_ihm, normalize_id


def test_loads_binarycif(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.bcif', fmt='BCIF'))
    assert cmd.count_atoms('m and tt. BEAD') == 5


def test_loads_gzipped_cif(tmp_path):
    load_ihm(one_group_file(tmp_path / 'm.cif.gz'))
    assert cmd.count_atoms('m and tt. BEAD') == 5


def test_pymol_load_with_ihm_format(tmp_path):
    pmihm.register()
    cmd.load(one_group_file(tmp_path / 'm.cif'), 'obj1', format='ihm')
    assert cmd.count_atoms('obj1 and tt. BEAD') == 5


def test_load_ihm_command_registered(tmp_path):
    pmihm.register()
    path = one_group_file(tmp_path / 'm.cif')
    cmd.do('load_ihm %s, obj2' % path)
    assert cmd.count_atoms('obj2 and tt. BEAD') == 5


def test_reserved_object_name_is_renamed_and_styled(tmp_path):
    # 'x' is a PyMOL keyword; PyMOL itself renames such objects to 'x_'
    load_ihm(one_group_file(tmp_path / 'm.cif'), object='x')
    assert cmd.count_atoms('x_ and tt. BEAD and rep spheres') == 5


@pytest.mark.parametrize('given, want', [
    ('PDBDEV_00000012', 'PDBDEV_00000012'),
    ('pdbdev_00000012', 'PDBDEV_00000012'),
    ('PDBDEV_12', 'PDBDEV_00000012'),
    ('12', 'PDBDEV_00000012'),
    ('9a8r', '9A8R'),
])
def test_entry_ids_normalized(given, want):
    assert normalize_id(given) == want


def test_bad_entry_id_rejected():
    with pytest.raises(ValueError):
        normalize_id('not an id')


def test_fetch_downloads_once_and_loads(tmp_path):
    server = tmp_path / 'server'
    server.mkdir()
    one_group_file(server / 'PDBDEV_00000007.cif')
    cache = tmp_path / 'cache'
    cache.mkdir()
    base_url = server.as_uri() + '/'

    fetch_ihm('pdbdev_7', path=str(cache), base_url=base_url)
    assert cmd.count_atoms('PDBDEV_00000007 and tt. BEAD') == 5
    assert (cache / 'PDBDEV_00000007.cif').exists()

    # second fetch uses the cached copy
    (server / 'PDBDEV_00000007.cif').unlink()
    fetch_ihm('PDBDEV_00000007', name='again', path=str(cache),
              base_url=base_url)
    assert cmd.count_atoms('again and tt. BEAD') == 5


def test_plugin_init_registers_commands_without_gui():
    pmihm.__init_plugin__()
    assert 'load_ihm' in cmd.keyword and 'fetch_ihm' in cmd.keyword


def test_filename_completion_for_load_ihm(tmp_path):
    pmihm.register()
    one_group_file(tmp_path / 'entry_for_completion.cif')
    from pymol import parser
    got = parser.Parser(cmd).complete('load_ihm %s/entry_for' % tmp_path)
    assert got == 'load_ihm %s/entry_for_completion.cif' % tmp_path


def test_fetch_creates_missing_cache_dir_and_accepts_url_without_slash(
        tmp_path):
    server = tmp_path / 'server'
    server.mkdir()
    one_group_file(server / 'PDBDEV_00000007.cif')
    cache = tmp_path / 'not' / 'there'
    fetch_ihm('7', path=str(cache), base_url=server.as_uri())
    assert cmd.count_atoms('PDBDEV_00000007 and tt. BEAD') == 5


def test_failed_download_is_pymol_error_and_leaves_no_partial_file(
        tmp_path):
    server = tmp_path / 'server'
    server.mkdir()
    cache = tmp_path / 'cache'
    cache.mkdir()
    with pytest.raises(pymol.CmdException):
        fetch_ihm('7', path=str(cache), base_url=server.as_uri())
    assert list(cache.iterdir()) == []


def test_bad_id_is_pymol_error():
    with pytest.raises(pymol.CmdException):
        fetch_ihm('not an id')


def serve(tmp_path, monkeypatch, *codes):
    server = tmp_path / 'server'
    server.mkdir()
    for code in codes:
        one_group_file(server / (code + '.cif'))
    import pmihm.loader
    monkeypatch.setattr(pmihm.loader, 'PDB_IHM_URL', server.as_uri() + '/')
    cache = tmp_path / 'cache'
    cache.mkdir()
    return str(cache)


def test_pymol_fetch_with_type_ihm(tmp_path, monkeypatch):
    cache = serve(tmp_path, monkeypatch, 'PDBDEV_00000007', '9A8R')
    pmihm.register()
    cmd.do('fetch pdbdev_7 9a8r, type=ihm, path=%s' % cache)
    assert cmd.count_atoms('PDBDEV_00000007 and tt. BEAD') == 5
    assert cmd.count_atoms('9A8R and tt. BEAD') == 5


def test_api_fetch_with_type_ihm_and_name(tmp_path, monkeypatch):
    cache = serve(tmp_path, monkeypatch, 'PDBDEV_00000007')
    pmihm.register()
    cmd.fetch('PDBDEV_00000007', 'mine', type='ihm', path=cache)
    assert cmd.count_atoms('mine and tt. BEAD') == 5


def test_pymol_fetch_routes_pdbdev_ids_to_pdb_ihm(tmp_path, monkeypatch):
    cache = serve(tmp_path, monkeypatch, 'PDBDEV_00000007')
    pmihm.register()
    cmd.fetch('PDBDEV_00000007', path=cache)
    assert cmd.count_atoms('PDBDEV_00000007 and tt. BEAD') == 5


def test_other_fetch_types_still_use_pymol(tmp_path):
    # PyMOL's own fetch loads an existing local file without downloading,
    # with PyMOL's own mmCIF reader, which never sets tt. BEAD.
    path = one_group_file(tmp_path / 'plain.cif')
    pmihm.register()
    cmd.fetch('1abc', 'plain', type='cif', file=path)
    assert cmd.count_atoms('plain') > 0
    assert cmd.count_atoms('plain and tt. BEAD') == 0


def test_not_quiet_reports_objects(tmp_path, capsys):
    load_ihm(two_group_file(tmp_path / 'm.cif'), quiet=0)
    out = capsys.readouterr().out
    assert 'm_cluster_1' in out and 'm_cluster_2' in out


def test_missing_file_is_pymol_error(tmp_path):
    with pytest.raises(pymol.CmdException):
        load_ihm(str(tmp_path / 'missing.cif'))


def test_unwritable_fetch_path_is_pymol_error(tmp_path):
    blocker = tmp_path / 'file'
    blocker.write_text('not a directory')
    with pytest.raises(pymol.CmdException):
        fetch_ihm('7', path=str(blocker / 'sub'))


def test_fetch_type_is_case_insensitive(tmp_path, monkeypatch):
    cache = serve(tmp_path, monkeypatch, '9A8R')
    pmihm.register()
    cmd.fetch('9a8r', type='IHM', path=cache)
    assert cmd.count_atoms('9A8R and tt. BEAD') == 5


def test_fetch_mixed_ids_sends_only_pdbdev_to_pdb_ihm(tmp_path, monkeypatch):
    cache = serve(tmp_path, monkeypatch, 'PDBDEV_00000007')
    # a local file for PyMOL's own fetch, so it doesn't download
    import shutil
    shutil.copy(one_group_file(tmp_path / 'x.cif'),
                tmp_path / 'cache' / '1abc.cif')
    pmihm.register()
    cmd.fetch('PDBDEV_00000007 1abc', path=cache, type='')
    assert cmd.count_atoms('PDBDEV_00000007 and tt. BEAD') == 5
    assert cmd.count_atoms('1abc') > 0
    assert cmd.count_atoms('1abc and tt. BEAD') == 0
