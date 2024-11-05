# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['mainwindow.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=['numpy', 'numpy.core._multiarray_umath', 'numpy.core.multiarray', 'numpy.core._dtype'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Relaxyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='BrandIcon.ico',
)
