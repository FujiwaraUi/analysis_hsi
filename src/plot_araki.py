"""
plot_araki.py
==============
claude_make_mock_araki.py が生成したモックキューブの可視化モジュール。
スクリプトとして直接実行すると全グループのキューブを生成してプロットを保存する。

依存パッケージ: numpy, matplotlib
"""

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from claude_make_mock_araki import (
    ALIS_WL,
    GROUP_SPECS,
    build_group_data,
    interp_to_alis,
    _load_csv,
)


def plot_calibration(
    ice_map:  np.ndarray,
    bd_map:   np.ndarray,
    metadata: dict,
    out_path: Path,
) -> None:
    """BD vs ice content の校正線検証プロット (論文 Fig. 8 に対応)"""
    ice_flat      = ice_map.ravel()
    bd_flat       = bd_map.ravel()
    gradient_exp  = metadata["gradient_experimental"]
    gradient_pred = metadata["gradient_predicted_eq1"]
    wl_min, wl_max = metadata["wavelength_range_nm"]

    ice_range    = np.linspace(0, ice_map.max() * 1.05, 100)
    bd_exp_line  = gradient_exp  * ice_range / 100.0
    bd_pred_line = gradient_pred * ice_range / 100.0

    fig, ax = plt.subplots(figsize=(8, 6))
    rng_plot = np.random.default_rng(0)
    idx = rng_plot.choice(len(ice_flat), min(3000, len(ice_flat)), replace=False)
    ax.scatter(ice_flat[idx], bd_flat[idx], s=4, alpha=0.25, color="steelblue",
               label="Mock data")
    ax.plot(ice_range, bd_exp_line,  "r-",  lw=2,
            label=f"Experimental gradient ({gradient_exp:.2f})")
    ax.plot(ice_range, bd_pred_line, "k--", lw=1.5,
            label=f"Eq.(1) predicted ({gradient_pred:.2f})")

    ax.set_xlabel("Water ice content (wt.%)", fontsize=12)
    ax.set_ylabel("1.5-μm Band Depth", fontsize=12)
    ax.set_title(
        f"Calibration: {metadata['mineral_type']} ({metadata['grain_size']})\n"
        f"R(1.5μm)={metadata['reflectance_1500']:.3f}, SNR={metadata['snr']:.0f} "
        f"[{wl_min:.0f}–{wl_max:.0f} nm]",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=-0.01)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_montage(
    cube:     np.ndarray,
    ice_map:  np.ndarray,
    bd_map:   np.ndarray,
    metadata: dict,
    out_path: Path,
) -> None:
    """代表バンド反射率 + 氷量マップ + BD マップのモンタージュ"""
    wl_min, wl_max = metadata["wavelength_range_nm"]

    # 実験データの有効範囲 (ALIS_WL[0]–ALIS_WL[-1]) 内から代表波長を選択
    target_wls   = [950, 1050, 1200, 1350, 1500, 1600]
    band_indices = [int(np.argmin(np.abs(ALIS_WL - wl))) for wl in target_wls]

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.ravel()
    vmin_r, vmax_r = float(np.percentile(cube, 2)), float(np.percentile(cube, 98))

    for i, bidx in enumerate(band_indices):
        wl = ALIS_WL[bidx]
        im = axes[i].imshow(cube[:, :, bidx], origin="lower",
                            cmap="viridis", vmin=vmin_r, vmax=vmax_r)
        axes[i].set_title(f"{wl:.0f} nm", fontsize=11)
        fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04, label="Reflectance")

    im_ice = axes[6].imshow(ice_map, origin="lower", cmap="Blues")
    axes[6].set_title("Ice content (wt.%)", fontsize=11)
    fig.colorbar(im_ice, ax=axes[6], fraction=0.046, pad=0.04, label="wt.%")

    im_bd = axes[7].imshow(bd_map, origin="lower", cmap="magma")
    axes[7].set_title("1.5-μm Band Depth", fontsize=11)
    fig.colorbar(im_bd, ax=axes[7], fraction=0.046, pad=0.04, label="BD")

    fig.suptitle(
        f"ALIS Mock (Araki data): {metadata['mineral_type']} {metadata['grain_size']} | "
        f"gradient={metadata['gradient_experimental']:.2f} (exp) / "
        f"{metadata['gradient_predicted_eq1']:.2f} (Eq.1) | "
        f"{wl_min:.0f}–{wl_max:.0f} nm\n"
        f"{metadata['reference']}",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_spectra_comparison(
    cube:       np.ndarray,
    ice_map:    np.ndarray,
    group_data: dict,
    metadata:   dict,
    out_path:   Path,
) -> None:
    """氷量別スペクトルの比較 (論文 Fig. 4 に対応)"""
    ice_flat   = ice_map.ravel()
    quantiles  = [0.05, 0.25, 0.50, 0.95]
    thresholds = np.quantile(ice_flat, quantiles)
    colors     = ["#2166ac", "#4393c3", "#f4a582", "#b2182b"]

    fig, ax = plt.subplots(figsize=(10, 6))
    for q, thresh, color in zip(quantiles, thresholds, colors):
        idx = int(np.argmin(np.abs(ice_flat - thresh)))
        iy, ix = np.unravel_index(idx, ice_map.shape)
        spec    = cube[iy, ix, :]
        ice_val = ice_map[iy, ix]
        ax.plot(ALIS_WL, spec, color=color, lw=1.2,
                label=f"ice={ice_val:.2f} wt.% (q={q:.0%})")

    ax.plot(ALIS_WL, group_data["dry_mean"], "k--", lw=1.8, label="Dry mean (real data)")

    ax.set_xlabel("Wavelength (nm)", fontsize=12)
    ax.set_ylabel("Reflectance", fontsize=12)
    ax.set_title(
        f"Spectral comparison: {metadata['mineral_type']} {metadata['grain_size']} "
        f"(cf. Araki & Saiki 2025, Fig. 4)\n"
        f"Dry spectra from real {metadata['n_dry_files']} measurements",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(ALIS_WL[0], ALIS_WL[-1])
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_dry_spectra_all_groups(
    data_dir: Path,
    out_path: Path,
) -> None:
    """全グループの平均乾燥スペクトルを 1 枚に重ね描き (Fig. 4 dry lines に対応)"""
    colors_min = {
        "olivine":       "#2ca02c",
        "plagioclase":   "#1f77b4",
        "clinopyroxene": "#d62728",
        "mixture":       "#ff7f0e",
    }
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    titles     = ["180-250 μm (coarse)", "75-125 μm (fine)"]
    grain_keys = ["coarse", "fine"]

    for col, (grain_key, title) in enumerate(zip(grain_keys, titles)):
        ax = axes[col]
        for mineral in ["olivine", "plagioclase", "clinopyroxene", "mixture"]:
            key  = (mineral, grain_key)
            spec = GROUP_SPECS[key]
            try:
                gd = build_group_data(spec, data_dir)
            except FileNotFoundError:
                continue
            color   = colors_min[mineral]
            raw_dir = data_dir / spec.mineral_folder / spec.grain_folder / "raw data"
            for f in sorted(raw_dir.glob("*_dry.csv")):
                wl, refl = _load_csv(f, has_header=False)
                ax.plot(ALIS_WL, interp_to_alis(wl, refl),
                        color=color, alpha=0.15, lw=0.8)
            ax.plot(ALIS_WL, gd["dry_mean"], color=color, lw=2.0, label=mineral)

        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        ax.set_ylabel("Reflectance", fontsize=11)
        ax.set_xlim(ALIS_WL[0], ALIS_WL[-1])
        ax.axvline(1500, color="gray", lw=0.8, ls=":")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    fig.suptitle(
        "Dry mineral spectra (real data, thin=individual, thick=mean)\n"
        "Araki & Saiki (2025) DOI: 10.60574/87068",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_unit_absorption_profiles(
    data_dir: Path,
    out_path: Path,
) -> None:
    """全グループの単位 BD 吸収プロファイルを比較するプロット"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    titles     = ["180-250 μm (coarse)", "75-125 μm (fine)"]
    grain_keys = ["coarse", "fine"]
    colors = {
        "olivine":       "#2ca02c",
        "plagioclase":   "#1f77b4",
        "clinopyroxene": "#d62728",
        "mixture":       "#ff7f0e",
    }

    for col, (grain_key, title) in enumerate(zip(grain_keys, titles)):
        ax = axes[col]
        for mineral in ["olivine", "plagioclase", "clinopyroxene", "mixture"]:
            key  = (mineral, grain_key)
            spec = GROUP_SPECS[key]
            try:
                gd = build_group_data(spec, data_dir)
            except FileNotFoundError:
                continue
            ax.plot(ALIS_WL, gd["unit_abs"], color=colors[mineral], lw=2.0, label=mineral)

        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        ax.set_ylabel("Unit absorption (per BD)", fontsize=11)
        ax.set_xlim(1200, ALIS_WL[-1])
        ax.set_ylim(bottom=-0.1)
        ax.axvline(1500, color="gray", lw=0.8, ls=":", label="1500 nm")
        ax.axhline(1.0,  color="gray", lw=0.8, ls="--")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    fig.suptitle(
        "Unit absorption profiles extracted from real data\n"
        "(= fractional absorption per unit BD, averaged per group)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
