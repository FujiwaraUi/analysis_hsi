"""
make_3d_cube.py
===============
`claude_make_mock_araki.py` から「モックキューブ生成」だけを抜き出したラッパー。

`main.py` は生成 + 各種プロット保存まで行うが、このスクリプトは
「cube を作って保存する」ことだけに用途を限定する。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from claude_make_mock_araki import ALISMockConfigAraki, DATA_DIR, create_alis_mock_araki


def main() -> None:
    p = argparse.ArgumentParser(description="Generate Araki mock cube only (no plots).")
    p.add_argument("--out-dir", type=Path, default=Path(__file__).parent / "alis_mock_output_araki" / "cube_only")
    p.add_argument("--mineral", choices=["olivine", "plagioclase", "clinopyroxene", "mixture"], default="olivine")
    p.add_argument("--grain", choices=["coarse", "fine"], default="coarse")
    p.add_argument("--ny", type=int, default=50)
    p.add_argument("--nx", type=int, default=100)
    p.add_argument("--snr", type=float, default=100.0)
    p.add_argument("--ice-pattern", choices=["gradient", "patches", "uniform"], default="gradient")
    p.add_argument("--ice-min", type=float, default=0.0)
    p.add_argument("--ice-max", type=float, default=2.2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--with-maps", action="store_true", help="Also store ice_map and bd_map in the npz.")
    args = p.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = ALISMockConfigAraki(
        mineral_type=args.mineral,  # type: ignore[arg-type]
        grain_size=args.grain,      # type: ignore[arg-type]
        ny=args.ny,
        nx=args.nx,
        ice_content_min=args.ice_min,
        ice_content_max=args.ice_max,
        snr=args.snr,
        ice_pattern=args.ice_pattern,  # type: ignore[arg-type]
        seed=args.seed,
    )

    cube, wavelengths, ice_map, bd_map, meta, _gd = create_alis_mock_araki(cfg, DATA_DIR)

    tag = (
        f"{args.mineral}_{args.grain}_ny{args.ny}_nx{args.nx}_snr{args.snr:.0f}_"
        f"{args.ice_pattern}_ice{args.ice_min:g}-{args.ice_max:g}_seed{args.seed}"
    )
    out_npz = out_dir / f"{tag}.npz"

    payload: dict[str, np.ndarray] = {"cube": cube, "wavelengths": wavelengths}
    if args.with_maps:
        payload["ice_map"] = ice_map
        payload["bd_map"] = bd_map
    np.savez_compressed(out_npz, **payload)

    print(f"saved: {out_npz}")
    print(f"cube: {cube.shape} dtype={cube.dtype}")
    print(f"wavelengths: {wavelengths[0]:.0f}–{wavelengths[-1]:.0f} nm ({len(wavelengths)} bands)")
    print(f"gradient (exp/Eq1): {meta['gradient_experimental']:.2f} / {meta['gradient_predicted_eq1']:.2f}")


if __name__ == "__main__":
    main()
