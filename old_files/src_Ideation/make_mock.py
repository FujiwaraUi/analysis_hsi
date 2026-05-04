import numpy as np
from pathlib import Path
import os


def save_band_montage(
    cube: np.ndarray,
    out_dir: str | Path = "mock_dta_output",
    step: int = 100,
):
    """
    キューブの特定バンドを並べて可視化して保存する。

    出力ファイル名は `mock_plot_0_to_{nlam}.png`。
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Matplotlib のキャッシュ/フォントキャッシュ先を確実に書き込み可能な場所へ
    mpl_config_dir = out_dir / ".mplconfig"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    xdg_cache_home = out_dir / ".cache"
    xdg_cache_home.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache_home))

    try:
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "matplotlib が必要です。`uv sync` または `pip install matplotlib` を実行してください。"
        ) from e

    if cube.ndim != 3:
        raise ValueError("cube must be (ny, nx, nlam)")

    ny, nx, nlam = cube.shape
    indices = list(range(0, nlam, step))
    if indices[-1] != nlam - 1:
        indices.append(nlam - 1)

    nplots = len(indices)
    ncols = min(3, nplots)
    nrows = int(np.ceil(nplots / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)
    axes_flat = axes.ravel()

    for ax, band_idx in zip(axes_flat, indices, strict=False):
        img = cube[:, :, band_idx]
        vmin, vmax = np.percentile(img, [1, 99])
        im = ax.imshow(img, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set_title(f"band {band_idx}")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes_flat[nplots:]:
        ax.axis("off")

    out_path = out_dir / f"mock_plot_0_to_{nlam}.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    return out_path

def generate_hsi_cube(
    nx=60,
    ny=60,
    nlam=500,
    mode="peak",  # "peak" or "mixture"
    add_noise=True,
    seed=0,
):
    """
    HSIキューブ生成

    Parameters
    ----------
    nx, ny : int
        空間サイズ
    nlam : int
        波長チャンネル数
    mode : str
        "peak"（ランダムピーク）または "mixture"（2成分混合）
    add_noise : bool
        ノイズ付加の有無
    seed : int
        乱数シード

    Returns
    -------
    cube : (ny, nx, nlam)
    noise_cube : 同形状
    """
    rng = np.random.default_rng(seed)

    # =========================
    # 1. スペクトル生成
    # =========================
    lam = np.arange(nlam)

    if mode == "peak":
        spectrum = 1.0 + rng.normal(0, 0.001, nlam)

        # ランダムピーク
        for _ in range(5):
            center = rng.integers(50, 150)
            amp = rng.uniform(10, 20)
            width = rng.uniform(2, 5)
            spectrum += amp * np.exp(-0.5 * ((lam - center) / width) ** 2)

    elif mode == "mixture":
        # 簡易テンプレート（fallback相当）
        x = (lam - lam.mean()) / (np.ptp(lam) + 1e-6)

        spec_old = (1.2 - 0.8 * x)
        spec_young = (1.1 - 0.4 * x)

        # 吸収線（簡易モデル）
        spec_old *= (1 - 0.08 * np.exp(-0.5 * ((lam - 200) / 5) ** 2))
        spec_young *= (1 - 0.06 * np.exp(-0.5 * ((lam - 250) / 8) ** 2))

        spectrum = (spec_old, spec_young)

    else:
        raise ValueError("mode must be 'peak' or 'mixture'")

    # =========================
    # 2. 空間構造
    # =========================
    xg, yg = np.meshgrid(np.arange(nx) - nx // 2,
                         np.arange(ny) - ny // 2)

    q = 0.6
    Reff = 12.0

    r = np.sqrt(xg**2 + (yg / q)**2)

    brightness = 10.0 * np.exp(-7.67 * ((r / Reff)**0.25 - 1))
    brightness = np.clip(brightness, 0.1, None)

    # =========================
    # 3. キューブ生成
    # =========================
    cube = np.zeros((ny, nx, nlam), dtype=np.float32)
    noise_cube = np.zeros_like(cube)

    for iy in range(ny):
        for ix in range(nx):

            if mode == "peak":
                spec = brightness[iy, ix] * spectrum

            else:  # mixture
                frac_old = np.clip(1.0 - 0.5 * (r[iy, ix] / Reff), 0.2, 1.0)
                frac_young = 1.0 - frac_old

                spec_old, spec_young = spectrum
                spec = brightness[iy, ix] * (
                    frac_old * spec_old + frac_young * spec_young
                )

            # ノイズ
            noise_std = np.sqrt(np.clip(spec, 0, None) + 1.0)

            if add_noise:
                noise = rng.normal(0, noise_std)
                cube[iy, ix, :] = spec + noise
            else:
                cube[iy, ix, :] = spec

            noise_cube[iy, ix, :] = noise_std

    return cube, noise_cube


# =========================
# 実行例
# =========================
if __name__ == "__main__":
    cube, noise = generate_hsi_cube(
        mode="mixture",
        add_noise=True
    )

    print("cube shape:", cube.shape)
    print("noise shape:", noise.shape)
#    out_path = save_band_montage(noise, out_dir="mock_dta_output", step=100)
    out_path = save_band_montage(cube, out_dir="mock_dta_output", step=100)
    print("saved plot:", out_path)
