"""
make_3d_cube.py
===============
`make_mock_only.py` で生成したモックキューブを、3Dボリュームとして可視化する。

注意:
  - ここでの Z 軸は空間ではなく「スペクトルバンド（波長方向）」。
  - Plotly を使用する（未インストールの場合はエラー表示して終了）。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from make_mock_only import MockConfig, make_mock_cube


def _maybe_import_plotly():
    try:
        import plotly.graph_objects as go  # type: ignore
    except ModuleNotFoundError as e:
        raise SystemExit(
            "plotly が必要です。`pip install plotly` を実行してから再実行してください。"
        ) from e
    return go


def _select_band_indices_by_step_nm(wavelengths: np.ndarray, step_nm: int) -> np.ndarray:
    wl_min = float(wavelengths[0])
    wl_max = float(wavelengths[-1])
    targets = np.arange(np.ceil(wl_min / step_nm) * step_nm, wl_max + 1e-6, step_nm)
    idx = [int(np.argmin(np.abs(wavelengths - t))) for t in targets]
    # unique & sorted (argmin can collide if step_nm < sampling)
    return np.array(sorted(set(idx)), dtype=int)


def _downsample_cube(cube: np.ndarray, step_xy: int, band_idx: np.ndarray | None) -> np.ndarray:
    cube2 = cube[:: max(1, step_xy), :: max(1, step_xy), :]
    if band_idx is None:
        return cube2
    return cube2[:, :, band_idx]


def visualize_cube_volume(
    cube: np.ndarray,
    wavelengths: np.ndarray,
    *,
    opacity: float = 0.15,
    surface_count: int = 20,
    colorscale: str = "Jet",
):
    go = _maybe_import_plotly()

    h, w, b = cube.shape
    x, y, z_idx = np.mgrid[0:h, 0:w, 0:b]
    z_nm = np.broadcast_to(wavelengths.reshape(1, 1, b), (h, w, b))

    fig = go.Figure(
        data=go.Volume(
            x=x.flatten(),
            y=y.flatten(),
            z=z_nm.astype(np.float32).flatten(),
            value=cube.astype(np.float32).flatten(),
            opacity=opacity,
            surface_count=surface_count,
            colorscale=colorscale,
        )
    )

    fig.update_layout(
        title="Mock Cube 3D Visualization (x,y,band)",
        width=1200,
        height=800,
        scene=dict(
            xaxis_title="Y [spatial]",
            yaxis_title="X [spatial]",
            zaxis_title="Spectral [nm]",
            aspectmode="cube",
        ),
    )
    return fig


def main() -> None:
    p = argparse.ArgumentParser(description="Make mock cube and visualize as 3D volume (Plotly).")
    p.add_argument("--out-dir", type=Path, default=Path(__file__).parent / "alis_mock_output_araki" / "mock_3d")
    p.add_argument("--mineral", choices=["olivine", "plagioclase", "clinopyroxene", "mixture"], default="olivine")
    p.add_argument("--grain", choices=["coarse", "fine"], default="coarse")
    p.add_argument("--ny", type=int, default=64)
    p.add_argument("--nx", type=int, default=64)
    p.add_argument("--snr", type=float, default=100.0)
    p.add_argument("--ice-pattern", choices=["gradient", "patches", "uniform"], default="gradient")
    p.add_argument("--ice-min", type=float, default=0.0, help="Minimum ice content (wt.%) used to build ice_map.")
    p.add_argument("--ice-max", type=float, default=2.2, help="Maximum ice content (wt.%) used to build ice_map.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--step-xy", type=int, default=2, help="Downsample step for spatial dims (speed).")
    p.add_argument(
        "--band-step-nm",
        type=int,
        default=100,
        help="Band sampling step in nm for Z axis (e.g., 100 => 900,1000,1100,...).",
    )
    p.add_argument("--opacity", type=float, default=0.15)
    p.add_argument("--surface-count", type=int, default=20)
    p.add_argument("--colorscale", type=str, default="Jet")
    p.add_argument("--no-show", action="store_true", help="Do not open an interactive window; just save.")
    args = p.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = MockConfig(
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

    cube, wavelengths = make_mock_cube(cfg)
    band_idx = _select_band_indices_by_step_nm(wavelengths, step_nm=args.band_step_nm) if args.band_step_nm > 0 else None
    cube_ds = _downsample_cube(cube, step_xy=args.step_xy, band_idx=band_idx)
    wl_ds = wavelengths if band_idx is None else wavelengths[band_idx]

    fig = visualize_cube_volume(
        cube_ds,
        wl_ds,
        opacity=args.opacity,
        surface_count=args.surface_count,
        colorscale=args.colorscale,
    )

    tag = f"{args.mineral}_{args.grain}_ny{args.ny}_nx{args.nx}_ds{args.step_xy}_nm{args.band_step_nm}"
    out_html = out_dir / f"{tag}.html"
    fig.write_html(out_html)
    print(f"saved: {out_html}")
    print(f"cube: {cube.shape} -> downsampled {cube_ds.shape}")
    print(f"wavelengths: {wl_ds[0]:.0f}–{wl_ds[-1]:.0f} nm ({len(wl_ds)} bands)")

    if not args.no_show:
        fig.show()


if __name__ == "__main__":
    main()
