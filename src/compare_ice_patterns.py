"""
compare_ice_patterns.py
=======================
ice_pattern = "gradient" / "patches" / "uniform" の 3 パターンを
同一条件で生成し、空間分布・BD 分布・スペクトルを並べて比較する。

出力 (alis_mock_output_araki/):
  pattern_comparison_spatial.png  — 氷量マップ・BD マップの空間比較
  pattern_comparison_bd_vs_ice.png — BD–ice 散布図と BD ヒストグラム
  pattern_comparison_spectra.png  — 代表ピクセルのスペクトル比較
"""

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from alis_mock_generator import (
    ALIS_WL,
    ALISMockConfig,
    DATA_DIR,
    generate_alis_mock_cube,
)

PATTERNS: list[str] = ["gradient", "patches", "uniform"]
PATTERN_COLORS: dict[str, str] = {
    "gradient": "#1f77b4",
    "patches":  "#d62728",
    "uniform":  "#2ca02c",
}
PATTERN_LABELS: dict[str, str] = {
    "gradient": "Gradient (center-high)",
    "patches":  "Patches (random clusters)",
    "uniform":  "Uniform (constant)",
}

MINERAL_TYPE = "mixture"
GRAIN_SIZE   = "coarse"
SNR          = 100
SEED         = 42


def generate_all_patterns(data_dir: Path) -> dict[str, dict]:
    results: dict[str, dict] = {}
    for pattern in PATTERNS:
        print(f"\n[{pattern}] Generating...")
        config = ALISMockConfig(
            mineral_type=MINERAL_TYPE,
            grain_size=GRAIN_SIZE,
            ice_pattern=pattern,  # type: ignore[arg-type]
            snr=SNR,
            seed=SEED,
        )
        cube, _, ice_map, bd_map, meta, gd = generate_alis_mock_cube(config, data_dir)
        results[pattern] = {
            "cube": cube,
            "ice_map": ice_map,
            "bd_map": bd_map,
            "meta": meta,
            "group_data": gd,
        }
    return results


def plot_spatial_comparison(results: dict[str, dict], out_path: Path) -> None:
    """Row 1: 氷量マップ、Row 2: BD マップ を 3 パターン並べて比較"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    ice_vmax = max(r["ice_map"].max() for r in results.values())
    bd_vmax  = max(r["bd_map"].max()  for r in results.values())

    for col, pattern in enumerate(PATTERNS):
        r = results[pattern]

        im_ice = axes[0, col].imshow(
            r["ice_map"], origin="lower", cmap="Blues", vmin=0, vmax=ice_vmax,
        )
        axes[0, col].set_title(PATTERN_LABELS[pattern], fontsize=11, fontweight="bold")
        axes[0, col].set_xlabel("x (pixel)")
        if col == 0:
            axes[0, col].set_ylabel("Ice content (wt.%)\ny (pixel)", fontsize=10)
        fig.colorbar(im_ice, ax=axes[0, col], fraction=0.046, pad=0.04, label="wt.%")

        im_bd = axes[1, col].imshow(
            r["bd_map"], origin="lower", cmap="magma", vmin=0, vmax=bd_vmax,
        )
        axes[1, col].set_xlabel("x (pixel)")
        if col == 0:
            axes[1, col].set_ylabel("1.5-μm Band Depth\ny (pixel)", fontsize=10)
        fig.colorbar(im_bd, ax=axes[1, col], fraction=0.046, pad=0.04, label="BD")

    meta0 = next(iter(results.values()))["meta"]
    fig.suptitle(
        f"Ice Pattern Comparison — {meta0['mineral_type']} ({meta0['grain_size']})\n"
        "Top: ice content maps  |  Bottom: 1.5-μm band depth maps",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_bd_vs_ice(results: dict[str, dict], out_path: Path) -> None:
    """BD vs ice 散布図 (左) と BD ヒストグラム (右) を並べる"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    rng_plot  = np.random.default_rng(0)

    # --- 散布図 ---
    ax = axes[0]
    for pattern in PATTERNS:
        r        = results[pattern]
        ice_flat = r["ice_map"].ravel()
        bd_flat  = r["bd_map"].ravel()
        idx = rng_plot.choice(len(ice_flat), min(2000, len(ice_flat)), replace=False)
        ax.scatter(
            ice_flat[idx], bd_flat[idx],
            s=5, alpha=0.35,
            color=PATTERN_COLORS[pattern],
            label=PATTERN_LABELS[pattern],
        )

    meta0     = next(iter(results.values()))["meta"]
    grad_exp  = meta0["gradient_experimental"]
    grad_pred = meta0["gradient_predicted_eq1"]
    ice_range = np.linspace(0, max(r["ice_map"].max() for r in results.values()) * 1.05, 100)
    ax.plot(ice_range, grad_exp  * ice_range / 100, "k-",  lw=2.0,
            label=f"Exp. gradient ({grad_exp:.2f})")
    ax.plot(ice_range, grad_pred * ice_range / 100, "k--", lw=1.5,
            label=f"Eq.(1) ({grad_pred:.2f})")
    ax.set_xlabel("Water ice content (wt.%)", fontsize=12)
    ax.set_ylabel("1.5-μm Band Depth", fontsize=12)
    ax.set_title("BD vs Ice content", fontsize=12)
    ax.legend(fontsize=9, markerscale=3)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=-0.005)

    # --- BD ヒストグラム ---
    ax2 = axes[1]
    for pattern in PATTERNS:
        bd_flat = results[pattern]["bd_map"].ravel()
        ax2.hist(
            bd_flat, bins=60, alpha=0.55,
            color=PATTERN_COLORS[pattern],
            label=PATTERN_LABELS[pattern],
            density=True,
        )
    ax2.set_xlabel("1.5-μm Band Depth", fontsize=12)
    ax2.set_ylabel("Density", fontsize=12)
    ax2.set_title("BD distribution", fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        f"BD–Ice Relationship and BD Distribution\n"
        f"{meta0['mineral_type']} ({meta0['grain_size']}), SNR={meta0['snr']:.0f}",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_spectra_comparison(results: dict[str, dict], out_path: Path) -> None:
    """高・中・低 氷量の代表ピクセルを 3 パターン重ね描き"""
    quantiles = [0.95, 0.50, 0.05]
    q_labels  = ["High ice (95th pctile)", "Mid ice (50th pctile)", "Low ice (5th pctile)"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for ax, q, qlabel in zip(axes, quantiles, q_labels):
        for pattern in PATTERNS:
            r        = results[pattern]
            ice_flat = r["ice_map"].ravel()
            thresh   = float(np.quantile(ice_flat, q))
            idx      = int(np.argmin(np.abs(ice_flat - thresh)))
            iy, ix   = np.unravel_index(idx, r["ice_map"].shape)
            ice_val  = r["ice_map"][iy, ix]
            ax.plot(
                ALIS_WL, r["cube"][iy, ix, :],
                color=PATTERN_COLORS[pattern], lw=1.5,
                label=f"{PATTERN_LABELS[pattern]}\n  ({ice_val:.2f} wt.%)",
            )

        gd = next(iter(results.values()))["group_data"]
        ax.plot(ALIS_WL, gd["dry_mean"], "k--", lw=1.5, alpha=0.8, label="Dry mean")

        ax.axvline(1500, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel("Reflectance", fontsize=11)
        ax.set_title(qlabel, fontsize=11, fontweight="bold")
        ax.set_xlim(ALIS_WL[0], ALIS_WL[-1])
        ax.legend(fontsize=7.5, loc="lower right")
        ax.grid(True, alpha=0.2)

    meta0 = next(iter(results.values()))["meta"]
    fig.suptitle(
        f"Spectra comparison by ice pattern — {meta0['mineral_type']} ({meta0['grain_size']})\n"
        "Representative pixels at 95th / 50th / 5th ice-content percentile",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    out_dir = Path(__file__).parent / "alis_mock_output_araki_comp"
    out_dir.mkdir(exist_ok=True)

    results = generate_all_patterns(DATA_DIR)

    plot_spatial_comparison(results, out_dir / "pattern_comparison_spatial.png")
    plot_bd_vs_ice(results,         out_dir / "pattern_comparison_bd_vs_ice.png")
    plot_spectra_comparison(results, out_dir / "pattern_comparison_spectra.png")

    print(f"\nDone. All outputs saved to: {out_dir}/")
