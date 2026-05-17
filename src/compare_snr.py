"""
compare_snr.py
==============
ice_pattern="gradient" 固定で SNR を 0, 50, 100, 150, 200 と変えて
ノイズ量がキューブ・BD・スペクトルに与える影響を比較する。

SNR=0 はノイズなし (noise_sigma=0) の理想ケース。
SNR が大きいほど noise_sigma = R(1.5μm)/SNR が小さくなる（高品質）。

出力 (alis_mock_output_araki_comp/):
  snr_comparison_spatial.png   — 氷量マップ(共通) + 各 SNR の BD マップ
  snr_comparison_bd_vs_ice.png — BD–ice 散布図と BD ヒストグラム
  snr_comparison_spectra.png   — 代表ピクセルのスペクトル比較
"""

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import numpy as np
from pathlib import Path

from alis_mock_generator import (
    ALIS_WL,
    ALISMockConfig,
    DATA_DIR,
    generate_alis_mock_cube,
)

SNR_VALUES: list[int] = [0, 50, 100, 150, 200]
MINERAL_TYPE = "mixture"
GRAIN_SIZE   = "coarse"
SEED         = 42

_cmap = mcm.get_cmap("plasma_r", len(SNR_VALUES))
SNR_COLORS: dict[int, tuple] = {snr: _cmap(i) for i, snr in enumerate(SNR_VALUES)}


def snr_label(snr: int) -> str:
    return "No noise (SNR=∞)" if snr == 0 else f"SNR={snr}"


def generate_all_snr(data_dir: Path) -> dict[int, dict]:
    results: dict[int, dict] = {}
    for snr in SNR_VALUES:
        print(f"\n[SNR={snr}] Generating...")
        config = ALISMockConfig(
            mineral_type=MINERAL_TYPE,
            grain_size=GRAIN_SIZE,
            ice_pattern="gradient",
            snr=snr,
            seed=SEED,
        )
        cube, _, ice_map, bd_map, meta, gd = generate_alis_mock_cube(config, data_dir)
        results[snr] = {
            "cube": cube,
            "ice_map": ice_map,
            "bd_map": bd_map,
            "meta": meta,
            "group_data": gd,
        }
    return results


def plot_spatial_comparison(results: dict[int, dict], out_path: Path) -> None:
    """氷量マップ(1枚共通) + 各 SNR の BD マップを 2×3 グリッドで表示"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes_flat = axes.ravel()

    # ice_map は全 SNR で同一 (同 seed / 同パターン)
    ice_map  = results[SNR_VALUES[0]]["ice_map"]
    bd_vmax  = max(r["bd_map"].max() for r in results.values())
    bd_vmin  = 0.0

    im_ice = axes_flat[0].imshow(ice_map, origin="lower", cmap="Blues",
                                  vmin=0, vmax=ice_map.max())
    axes_flat[0].set_title("Ice content (shared)", fontsize=11, fontweight="bold")
    axes_flat[0].set_xlabel("x (pixel)")
    axes_flat[0].set_ylabel("y (pixel)")
    fig.colorbar(im_ice, ax=axes_flat[0], fraction=0.046, pad=0.04, label="wt.%")

    for ax_idx, snr in enumerate(SNR_VALUES, start=1):
        bd_map = results[snr]["bd_map"]
        im_bd  = axes_flat[ax_idx].imshow(
            bd_map, origin="lower", cmap="magma", vmin=bd_vmin, vmax=bd_vmax,
        )
        axes_flat[ax_idx].set_title(snr_label(snr), fontsize=11, fontweight="bold")
        axes_flat[ax_idx].set_xlabel("x (pixel)")
        axes_flat[ax_idx].set_ylabel("y (pixel)")
        fig.colorbar(im_bd, ax=axes_flat[ax_idx], fraction=0.046, pad=0.04, label="BD")

    meta0 = next(iter(results.values()))["meta"]
    fig.suptitle(
        f"SNR Comparison — {meta0['mineral_type']} ({meta0['grain_size']}), "
        "ice_pattern=gradient\n"
        "Top-left: ice content (common)  |  Others: 1.5-μm band depth",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_bd_vs_ice(results: dict[int, dict], out_path: Path) -> None:
    """BD vs ice 散布図 (左) と BD ヒストグラム (右) を各 SNR で重ね合わせ"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    rng_plot  = np.random.default_rng(0)

    meta0     = next(iter(results.values()))["meta"]
    grad_exp  = meta0["gradient_experimental"]
    grad_pred = meta0["gradient_predicted_eq1"]
    ice_range = np.linspace(
        0,
        max(r["ice_map"].max() for r in results.values()) * 1.05,
        100,
    )

    # --- 散布図 ---
    ax = axes[0]
    for snr in SNR_VALUES:
        r        = results[snr]
        ice_flat = r["ice_map"].ravel()
        bd_flat  = r["bd_map"].ravel()
        idx      = rng_plot.choice(len(ice_flat), min(1500, len(ice_flat)), replace=False)
        ax.scatter(
            ice_flat[idx], bd_flat[idx],
            s=5, alpha=0.35,
            color=SNR_COLORS[snr],
            label=snr_label(snr),
        )

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

    # --- BD ヒストグラム ---
    ax2 = axes[1]
    for snr in SNR_VALUES:
        bd_flat = results[snr]["bd_map"].ravel()
        ax2.hist(
            bd_flat, bins=60, alpha=0.55,
            color=SNR_COLORS[snr],
            label=snr_label(snr),
            density=True,
        )
    ax2.set_xlabel("1.5-μm Band Depth", fontsize=12)
    ax2.set_ylabel("Density", fontsize=12)
    ax2.set_title("BD distribution", fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        f"BD–Ice Relationship and BD Distribution\n"
        f"{meta0['mineral_type']} ({meta0['grain_size']}), ice_pattern=gradient",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_spectra_comparison(results: dict[int, dict], out_path: Path) -> None:
    """高・中・低氷量の代表ピクセルに対して各 SNR のスペクトルを重ね描き"""
    quantiles = [0.95, 0.50, 0.05]
    q_labels  = ["High ice (95th pctile)", "Mid ice (50th pctile)", "Low ice (5th pctile)"]

    # 代表ピクセル座標は SNR=0 の ice_map (全 SNR 共通) から決める
    ice_map_ref = results[SNR_VALUES[0]]["ice_map"]
    ice_flat    = ice_map_ref.ravel()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for ax, q, qlabel in zip(axes, quantiles, q_labels):
        thresh = float(np.quantile(ice_flat, q))
        idx    = int(np.argmin(np.abs(ice_flat - thresh)))
        iy, ix = np.unravel_index(idx, ice_map_ref.shape)
        ice_val = ice_map_ref[iy, ix]

        for snr in SNR_VALUES:
            spec = results[snr]["cube"][iy, ix, :]
            ax.plot(
                ALIS_WL, spec,
                color=SNR_COLORS[snr], lw=1.4,
                label=snr_label(snr),
            )

        gd = next(iter(results.values()))["group_data"]
        ax.plot(ALIS_WL, gd["dry_mean"], "k--", lw=1.5, alpha=0.8, label="Dry mean")

        ax.axvline(1500, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel("Reflectance", fontsize=11)
        ax.set_title(f"{qlabel}\nice = {ice_val:.2f} wt.%", fontsize=11, fontweight="bold")
        ax.set_xlim(ALIS_WL[0], ALIS_WL[-1])
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(True, alpha=0.2)

    meta0 = next(iter(results.values()))["meta"]
    fig.suptitle(
        f"Spectra comparison by SNR — {meta0['mineral_type']} ({meta0['grain_size']}), "
        "ice_pattern=gradient\n"
        "Same pixel, different noise levels",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    out_dir = Path(__file__).parent / "alis_mock_output_araki_comp"
    out_dir.mkdir(exist_ok=True)

    results = generate_all_snr(DATA_DIR)

    plot_spatial_comparison(results, out_dir / "snr_comparison_spatial.png")
    plot_bd_vs_ice(results,          out_dir / "snr_comparison_bd_vs_ice.png")
    plot_spectra_comparison(results, out_dir / "snr_comparison_spectra.png")

    print(f"\nDone. All outputs saved to: {out_dir}/")
