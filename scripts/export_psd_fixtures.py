#!/usr/bin/env python3
"""Export PSD fixtures from TIPTOP .ini configurations for MASTSEL regression tests.

This script runs TIPTOP simulations to collect PSD and geometry inputs used by
`mastsel.mavisPsf.psdSetToPsfSet`, then saves them as compressed NumPy fixtures.
The resulting fixtures can be used in MASTSEL tests without requiring runtime
P3/TIPTOP dependencies.
"""

from __future__ import annotations

import argparse
import json
import os
from configparser import ConfigParser
from pathlib import Path

import numpy as np

from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray


DEFAULT_CASES = [
    "MAVIS",
    "ERIS",
    "SPHERE",
    "SOUL",
    "HARMONI_SCAO",
]


def _ensure_numpy_array(x):
    return np.asarray(cpuArray(x))


def export_case(case_name: str, ini_dir: Path, output_dir: Path) -> Path:
    # Instantiate and run once to let TIPTOP compute PSD and all derived scales.
    simulation = baseSimulation(
        path=str(ini_dir),
        parametersFile=case_name,
        outputDir=str(output_dir),
        outputFile=f"_tmp_{case_name}",
        doConvolve=False,
        doPlot=False,
        addSrAndFwhm=False,
        verbose=False,
        savePSDs=False,
    )
    simulation.doOverallSimulation()

    n_pointings = simulation.pointings.shape[1]
    psd_ho = _ensure_numpy_array(simulation.PSD[0:n_pointings])
    mask = _ensure_numpy_array(simulation.mask.sampling)

    fixture_path = output_dir / f"{case_name}.npz"
    np.savez_compressed(
        fixture_path,
        case_name=np.asarray(case_name),
        input_psds=psd_ho,
        mask=mask,
        wavelength=np.asarray(simulation.wvl, dtype=np.float64),
        N=np.asarray(int(simulation.N)),
        nPixPup=np.asarray(int(simulation.sx)),
        grid_diameter=np.asarray(float(simulation.grid_diameter)),
        freq_range=np.asarray(float(simulation.freq_range)),
        dk=np.asarray(float(simulation.dk)),
        nPixPsf=np.asarray(int(simulation.nPixPSF)),
        wvlRef=np.asarray(float(simulation.wvlRef)),
        oversampling=np.asarray(float(simulation.overSamp)),
        padPSD=np.asarray(bool(simulation.nWvl > 1)),
        has_opd=np.asarray(simulation.opdMap is not None),
        opd_map=(
            _ensure_numpy_array(simulation.opdMap)
            if simulation.opdMap is not None
            else np.asarray([], dtype=np.float64)
        ),
    )
    return fixture_path


def export_multi_wavelength_case(
    case_name: str, wavelengths_m: list[float], ini_dir: Path, output_dir: Path
) -> Path:
    """
    Export a per-wavelength PSD fixture (exactMultiWavelengthPSD=True): P3
    computes one exact PSD grid per science wavelength instead of one shared
    grid, so unlike export_case() the resulting PSD arrays can have different
    shapes per wavelength -- saved as separate input_psds_{i} entries plus the
    matching per-wavelength freq_range/dk/nPixPup, rather than a single
    input_psds array. Used by MASTSEL's per-wavelength psdSetToPsfSet
    regression tests, kept in a separate fixture directory from export_case's
    output so the legacy single-grid fixture-scanning tests
    (report_nyquist_zeroing_impact.py, test_psf_fixture_regression.py) never
    pick it up by accident.
    """
    base_ini_path = ini_dir / f"{case_name}.ini"
    if not base_ini_path.exists():
        raise FileNotFoundError(f"Missing INI file: {base_ini_path}")

    config = ConfigParser()
    config.optionxform = str
    config.read(base_ini_path)
    config.set("sources_science", "Wavelength", str(list(wavelengths_m)))

    tmp_ini_dir = output_dir / "_tmp_ini"
    tmp_ini_dir.mkdir(parents=True, exist_ok=True)
    tmp_basename = f"{case_name}_multiwvl"
    with open(tmp_ini_dir / f"{tmp_basename}.ini", "w") as f:
        config.write(f)

    simulation = baseSimulation(
        path=str(tmp_ini_dir),
        parametersFile=tmp_basename,
        outputDir=str(output_dir),
        outputFile=f"_tmp_{tmp_basename}",
        doConvolve=False,
        doPlot=False,
        addSrAndFwhm=False,
        verbose=False,
        savePSDs=False,
        exactMultiWavelengthPSD=True,
    )
    simulation.doOverallSimulation()

    if not simulation.multiGridPSD:
        raise RuntimeError(
            f"{case_name}: exactMultiWavelengthPSD=True did not produce a "
            f"per-wavelength PSD list (nWvl={simulation.nWvl}) -- need at "
            f"least two distinct wavelengths."
        )

    n_pointings = simulation.pointings.shape[1]
    n_wvl = len(wavelengths_m)

    payload = {
        "case_name": np.asarray(case_name),
        "wavelength": np.asarray(simulation.wvl, dtype=np.float64),
        "n_wvl": np.asarray(n_wvl),
        "mask": _ensure_numpy_array(simulation.mask.sampling),
        "nPixPsf": np.asarray(int(simulation.nPixPSF)),
        "freq_range_per_wvl": np.asarray(simulation.freq_range_per_wvl, dtype=np.float64),
        "dk_per_wvl": np.asarray(simulation.dk_per_wvl, dtype=np.float64),
        "nPixPup_per_wvl": np.asarray(simulation.sx_per_wvl, dtype=np.int64),
        "has_opd": np.asarray(simulation.opdMap is not None),
        "opd_map": (
            _ensure_numpy_array(simulation.opdMap)
            if simulation.opdMap is not None
            else np.asarray([], dtype=np.float64)
        ),
    }
    for i in range(n_wvl):
        payload[f"input_psds_{i}"] = _ensure_numpy_array(simulation.PSD[i][0:n_pointings])

    fixture_path = output_dir / f"{case_name}_multiwvl.npz"
    np.savez_compressed(fixture_path, **payload)
    return fixture_path


def select_psd_subset(psd_ho: np.ndarray, keep_single: bool) -> np.ndarray:
    if keep_single:
        return np.asarray(psd_ho[0:1])
    return psd_ho


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ini-dir",
        default="tiptop/perfTest",
        help="Directory containing TIPTOP .ini files (default: tiptop/perfTest)",
    )
    parser.add_argument(
        "--out-dir",
        default="../MASTSEL/tests/fixtures/psd",
        help="Output directory for .npz fixtures (default: ../MASTSEL/tests/fixtures/psd)",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=DEFAULT_CASES,
        help=f"Case names without .ini (default: {' '.join(DEFAULT_CASES)})",
    )
    parser.add_argument(
        "--single-psd-cases",
        nargs="+",
        default=["MAVIS"],
        help="Cases for which only the first HO PSD is stored (default: MAVIS)",
    )
    parser.add_argument(
        "--multi-wavelength-case",
        default=None,
        help="If set, additionally export ONE exactMultiWavelengthPSD=True "
             "fixture for this case (e.g. MAVIS), using --multi-wavelength-nm.",
    )
    parser.add_argument(
        "--multi-wavelength-nm",
        nargs="+",
        type=float,
        default=[450.0, 550.0, 850.0],
        help="Wavelengths in nm for --multi-wavelength-case (default: 450 550 850)",
    )
    parser.add_argument(
        "--multi-wavelength-out-dir",
        default="../MASTSEL/tests/fixtures/psd_multiwavelength",
        help="Output directory for the per-wavelength fixture (kept separate "
             "from --out-dir so legacy single-grid fixture-scanning tests "
             "never pick it up).",
    )

    args = parser.parse_args()

    ini_dir = Path(args.ini_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    exported = []
    single_psd_cases = set(args.single_psd_cases)
    for case in args.cases:
        ini_path = ini_dir / f"{case}.ini"
        if not ini_path.exists():
            raise FileNotFoundError(f"Missing INI file: {ini_path}")
        fixture_path = export_case(case, ini_dir, out_dir)

        # Optionally replace the fixture content with a single-PSD subset.
        if case in single_psd_cases:
            data = np.load(fixture_path, allow_pickle=False)
            payload = {k: data[k] for k in data.files}
            payload["input_psds"] = select_psd_subset(np.asarray(payload["input_psds"]), keep_single=True)
            np.savez_compressed(fixture_path, **payload)

        exported.append(str(fixture_path))
        print(f"Exported: {fixture_path}")

    metadata_path = out_dir / "manifest.json"
    with metadata_path.open("w", encoding="utf-8") as f:
        json.dump({"fixtures": exported}, f, indent=2)
    print(f"Wrote manifest: {metadata_path}")

    if args.multi_wavelength_case:
        mw_out_dir = Path(args.multi_wavelength_out_dir).resolve()
        mw_out_dir.mkdir(parents=True, exist_ok=True)
        wavelengths_m = [nm * 1e-9 for nm in args.multi_wavelength_nm]
        mw_fixture_path = export_multi_wavelength_case(
            args.multi_wavelength_case, wavelengths_m, ini_dir, mw_out_dir
        )
        print(f"Exported multi-wavelength fixture: {mw_fixture_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
