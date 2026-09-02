# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules


datas = []
binaries = []
hiddenimports = [
    "astroquery.imcce",
    "astropy.io.fits",
    "astropy.time",
    "astropy.wcs",
    "astropy.wcs._wcs",
    "astropy.coordinates",
    "astropy.units",
    "scipy.ndimage",
    "scipy.optimize",
    "exoplanet_core",
    "sep",
    *collect_submodules("astroquery.imcce"),
]

datas += collect_data_files("astropy")
datas += collect_data_files("astropy_iers_data")
datas += collect_data_files("astroquery")
datas += collect_data_files("pyvo")
datas += collect_data_files("scipy")
datas += collect_data_files("exoplanet_core")
datas += [("../voyager_alpha/assets", "voyager_alpha/assets")]
binaries += collect_dynamic_libs("astropy")
binaries += collect_dynamic_libs("scipy")
binaries += collect_dynamic_libs("exoplanet_core")


a = Analysis(
    ["../voyager_alpha/main.py"],
    pathex=[".."],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=["pyinstaller_runtime_hook.py"],
    excludes=[
        "matplotlib",
        "IPython",
        "pandas",
        "tkinter",
        "PyQt6.QtQml",
        "PyQt6.QtQuick",
        "PyQt6.QtWebEngineWidgets",
    ],
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
    name="Voyager-Alpha",
    icon="../voyager_alpha/assets/voyager-alpha.ico",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
