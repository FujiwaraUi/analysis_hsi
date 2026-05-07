"""
make_mock_only.py
=================
他ファイルに依存せず、単体で動く「モックキューブ生成だけ」のスクリプト。

本ファイルは、元の `claude_make_mock_araki.py` / `plot_araki.py` の
データ読込や検証・可視化を切り離し、最低限のモック生成に限定している。

出力:
  - *.npz : cube と wavelengths のみ
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np


ALIS_WL = np.arange(900.0, 1640.0 + 0.1, 5.0)  # 900–1640nm, 5nm step


@dataclass
class GroupSpec:
    mineral_type: str
    grain_size: Literal["coarse", "fine"]
    grain_size_um: float
    reflectance_1500: float
    gradient_experimental: float


# Table 2 (Araki & Saiki 2025) に対応する値を、モック用途として定数化
GROUP_SPECS: dict[tuple[str, str], GroupSpec] = {
    ("olivine",       "coarse"): GroupSpec("olivine",       "coarse", 215, 0.460, 5.18),
    ("olivine",       "fine"):   GroupSpec("olivine",       "fine",   100, 0.542, 1.83),
    ("plagioclase",   "coarse"): GroupSpec("plagioclase",   "coarse", 215, 0.565, 8.48),
    ("plagioclase",   "fine"):   GroupSpec("plagioclase",   "fine",   100, 0.679, 5.83),
    ("clinopyroxene", "coarse"): GroupSpec("clinopyroxene", "coarse", 215, 0.610, 11.87),
    ("clinopyroxene", "fine"):   GroupSpec("clinopyroxene", "fine",   100, 0.665, 5.91),
    ("mixture",       "coarse"): GroupSpec("mixture",       "coarse", 215, 0.509, 5.71),
    ("mixture",       "fine"):   GroupSpec("mixture",       "fine",   100, 0.634, 3.11),
}


@dataclass
class MockConfig:
    mineral_type: Literal["olivine", "plagioclase", "clinopyroxene", "mixture"] = "olivine"
    grain_size: Literal["coarse", "fine"] = "coarse"
    ny: int = 50
    nx: int = 100
    ice_content_min: float = 0.0
    ice_content_max: float = 2.2
    ice_pattern: Literal["gradient", "patches", "uniform"] = "gradient"
    seed: int = 0
    snr: float = 100.0


def _gauss(wl: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    x = (wl - mu) / sigma
    return np.exp(-0.5 * x * x)


def _dry_spectrum_template(spec: GroupSpec, wl: np.ndarray) -> np.ndarray:
    """
    実測CSVを使わず、鉱物ごとの典型形状を簡易に再現した乾燥スペクトルを作る。
    1500nmでの反射率が spec.reflectance_1500 になるように正規化する。
    """
    base = np.ones_like(wl, dtype=np.float64)

    if spec.mineral_type == "plagioclase":
        # ほぼ平坦
        base *= 1.0 + 0.02 * (wl - 1500.0) / 600.0
    elif spec.mineral_type == "olivine":
        # ~1000nm付近の吸収を簡易に表現
        base *= 1.0 - 0.55 * _gauss(wl, 1050.0, 85.0)
    elif spec.mineral_type == "clinopyroxene":
        # ~900nm と ~1000nm の2つの吸収を簡易に表現
        base *= 1.0 - 0.35 * _gauss(wl, 920.0, 55.0) - 0.28 * _gauss(wl, 1030.0, 70.0)
    elif spec.mineral_type == "mixture":
        # 混合は中間（簡易に plagioclase と olivine と cpx の平均）
        p = _dry_spectrum_template(GROUP_SPECS[("plagioclase", spec.grain_size)], wl)
        o = _dry_spectrum_template(GROUP_SPECS[("olivine", spec.grain_size)], wl)
        c = _dry_spectrum_template(GROUP_SPECS[("clinopyroxene", spec.grain_size)], wl)
        base = (p + o + c) / 3.0
    else:
        raise ValueError(f"Unknown mineral_type: {spec.mineral_type}")

    # 粒径が大きいほど暗くなる傾向を、簡易なスケールで反映
    if spec.grain_size == "coarse":
        base *= 0.92
    else:
        base *= 1.00

    r1500 = float(np.interp(1500.0, wl, base))
    if r1500 <= 0:
        r1500 = 1e-6
    base *= spec.reflectance_1500 / r1500
    return np.clip(base, 0.01, 1.5)


def _unit_absorption_profile(wl: np.ndarray) -> np.ndarray:
    """
    1.5μm水氷吸収帯の形状（単位BDあたりの吸収プロファイル）を簡易に生成。
    peak=1 となるように正規化して返す。
    """
    prof = (
        1.00 * _gauss(wl, 1500.0, 55.0)
        + 0.20 * _gauss(wl, 1450.0, 35.0)
        + 0.10 * _gauss(wl, 1580.0, 60.0)
    )
    m = float(prof.max())
    return prof / (m if m > 0 else 1.0)


def _create_ice_map(cfg: MockConfig, rng: np.random.Generator) -> np.ndarray:
    if cfg.ice_pattern == "uniform":
        return np.full((cfg.ny, cfg.nx), cfg.ice_content_max, dtype=np.float32)

    if cfg.ice_pattern == "gradient":
        yy, xx = np.mgrid[0 : cfg.ny, 0 : cfg.nx]
        cy, cx = (cfg.ny - 1) / 2.0, (cfg.nx - 1) / 2.0
        r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
        r = r / (r.max() if r.max() > 0 else 1.0)
        ice = cfg.ice_content_min + (cfg.ice_content_max - cfg.ice_content_min) * r
        return ice.astype(np.float32)

    if cfg.ice_pattern == "patches":
        ice = np.full((cfg.ny, cfg.nx), cfg.ice_content_min, dtype=np.float64)
        n_patches = 5
        for _ in range(n_patches):
            cy = rng.integers(0, cfg.ny)
            cx = rng.integers(0, cfg.nx)
            rad = rng.uniform(min(cfg.ny, cfg.nx) * 0.08, min(cfg.ny, cfg.nx) * 0.20)
            amp = rng.uniform(0.3, 1.0)
            yy, xx = np.mgrid[0 : cfg.ny, 0 : cfg.nx]
            mask = (yy - cy) ** 2 + (xx - cx) ** 2 <= rad * rad
            ice[mask] += amp
        ice = np.clip(ice, 0.0, 1.0)
        ice = cfg.ice_content_min + (cfg.ice_content_max - cfg.ice_content_min) * ice
        return ice.astype(np.float32)

    raise ValueError(f"Unknown ice_pattern: {cfg.ice_pattern}")


def make_mock_cube(cfg: MockConfig) -> tuple[np.ndarray, np.ndarray]:
    """
    モックキューブだけを生成して返す。
    戻り値:
      cube: (ny, nx, n_bands) float32
      wavelengths: (n_bands,) float64
    """
    rng = np.random.default_rng(cfg.seed)
    key = (cfg.mineral_type, cfg.grain_size)
    if key not in GROUP_SPECS:
        raise KeyError(f"Unknown group key {key}. Valid: {list(GROUP_SPECS)}")
    spec = GROUP_SPECS[key]

    dry = _dry_spectrum_template(spec, ALIS_WL).astype(np.float64)
    unit_abs = _unit_absorption_profile(ALIS_WL).astype(np.float64)
    ice_map = _create_ice_map(cfg, rng).astype(np.float64)

    gradient = float(spec.gradient_experimental)
    target_bd_map = gradient * ice_map / 100.0

    # 乾燥スペクトルに対して、BDに比例した吸収を乗せる（簡易モデル）
    cube = np.empty((cfg.ny, cfg.nx, len(ALIS_WL)), dtype=np.float32)
    noise_sigma = float(spec.reflectance_1500 / cfg.snr) if cfg.snr > 0 else 0.0

    for iy in range(cfg.ny):
        for ix in range(cfg.nx):
            target_bd = float(target_bd_map[iy, ix])
            absorption_factor = np.clip(1.0 - target_bd * unit_abs, 0.0, None)
            spec_px = dry * absorption_factor
            if noise_sigma > 0:
                spec_px = spec_px + rng.normal(0.0, noise_sigma, size=spec_px.shape)
            cube[iy, ix, :] = np.clip(spec_px, 0.01, 1.5).astype(np.float32)

    return cube, ALIS_WL.copy()


def make_mock_only(
    *,
    mineral: Literal["olivine", "plagioclase", "clinopyroxene", "mixture"] = "olivine",
    grain: Literal["coarse", "fine"] = "coarse",
    ny: int = 50,
    nx: int = 100,
    snr: float = 100.0,
    ice_pattern: Literal["gradient", "patches", "uniform"] = "gradient",
    seed: int = 0,
) -> np.ndarray:
    """
    他コードから呼び出す用の薄いラッパー。

    例:
      from make_mock_only import make_mock_only
      cube = make_mock_only(mineral="olivine", grain="fine", ny=50, nx=100, snr=100)
    """
    cfg = MockConfig(
        mineral_type=mineral,
        grain_size=grain,
        ny=ny,
        nx=nx,
        snr=snr,
        ice_pattern=ice_pattern,
        seed=seed,
    )
    cube, _wl = make_mock_cube(cfg)
    return cube


def _save_cube(out_dir: Path, cfg: MockConfig) -> None:
    tag = f"{cfg.mineral_type}_{cfg.grain_size}"
    cube, wl = make_mock_cube(cfg)
    out_npz = out_dir / f"{tag}.npz"
    np.savez_compressed(out_npz, cube=cube, wavelengths=wl)
    print(f"saved: {out_npz} (cube={cube.shape}, bands={len(wl)})")


def main() -> None:
    p = argparse.ArgumentParser(description="Generate mock cubes only (standalone).")
    p.add_argument("--out-dir", type=Path, default=Path(__file__).parent / "alis_mock_output_araki" / "mock_only")
    p.add_argument("--all", action="store_true", help="Generate all 8 groups (default when mineral/grain omitted).")
    p.add_argument("--mineral", choices=["olivine", "plagioclase", "clinopyroxene", "mixture"])
    p.add_argument("--grain", choices=["coarse", "fine"])
    p.add_argument("--ny", type=int, default=50)
    p.add_argument("--nx", type=int, default=100)
    p.add_argument("--snr", type=float, default=100.0)
    p.add_argument("--ice-pattern", choices=["gradient", "patches", "uniform"], default="gradient")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.all or (args.mineral is None and args.grain is None):
        for (mineral, grain) in list(GROUP_SPECS.keys()):
            cfg = MockConfig(
                mineral_type=mineral,  # type: ignore[arg-type]
                grain_size=grain,      # type: ignore[arg-type]
                ny=args.ny,
                nx=args.nx,
                snr=args.snr,
                ice_pattern=args.ice_pattern,  # type: ignore[arg-type]
                seed=args.seed,
            )
            _save_cube(out_dir, cfg)
        return

    if args.mineral is None or args.grain is None:
        raise SystemExit("--mineral and --grain must be specified when not using --all.")

    cfg = MockConfig(
        mineral_type=args.mineral,  # type: ignore[arg-type]
        grain_size=args.grain,      # type: ignore[arg-type]
        ny=args.ny,
        nx=args.nx,
        snr=args.snr,
        ice_pattern=args.ice_pattern,  # type: ignore[arg-type]
        seed=args.seed,
    )
    _save_cube(out_dir, cfg)


if __name__ == "__main__":
    main()
