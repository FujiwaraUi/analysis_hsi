import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

def save_lunar_montage_fixed_scale(cube, wavelengths, vmin=0.35, vmax=0.75, out_dir="lunar_ice_mock_fixed"):
    """全パネルで同じcolorbar範囲に固定"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ny, nx, nlam = cube.shape
    step = 40
    indices = list(range(0, nlam, step))
    if indices[-1] != nlam - 1:
        indices.append(nlam - 1)
    
    nplots = len(indices)
    ncols = min(4, nplots)
    nrows = int(np.ceil(nplots / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.2 * nrows), squeeze=False)
    axes = axes.ravel()
    
    for i, (ax, bidx) in enumerate(zip(axes, indices)):
        wl = wavelengths[bidx]
        img = cube[:, :, bidx]
        
        im = ax.imshow(img, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_title(f'{wl:.0f} nm')
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Reflectance')
    
    for ax in axes[nplots:]:
        ax.axis('off')
    
    fig.suptitle('Lunar PSR Water Ice Mock\n(Fixed Color Scale: 0.35 - 0.75)', fontsize=16)
    fig.tight_layout()
    out_path = out_dir / "lunar_ice_fixed_scale.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    return out_path


def create_lunar_ice_mock(
    nx: int = 120,
    ny: int = 120,
    nlam: int = 350,
    base_reflectance: float = 0.55,
    seed: int = 42,
    ice_content_map: np.ndarray | None = None,
):
    """
    論文「Estimation of calibration lines for water ice content...」に基づいた
    月極域水氷ハイパースペクトルイメージキューブを生成する。
    
    主な特徴:
    - 900-1650nm範囲
    - 1.5μm (1500nm) 付近に水氷の吸収バンド
    - バンド深さが水氷含有量にほぼ比例（論文の校正線を模擬）
    """
    rng = np.random.default_rng(seed)
    lam = np.linspace(900, 1650, nlam)  # nm
    
    # 乾燥基材の連続スペクトル（やや減少傾向）
    continuum = base_reflectance * (1.05 - 0.0003 * (lam - 900))
    
    # 水氷吸収パラメータ（1.5μm中心）
    ice_center = 1500.0
    ice_width = 80.0   # nm
    
    # 水氷含有量マップ（空間的に変化）
    if ice_content_map is None:
        x, y = np.meshgrid(np.arange(nx), np.arange(ny))
        # 模擬的な分布（波状 + ランダム）
        ice_content_map = 1.2 + 1.1 * np.sin(2*np.pi*x/45) * np.cos(2*np.pi*y/38)
        ice_content_map += rng.normal(0, 0.4, (ny, nx))
        ice_content_map = np.clip(ice_content_map, 0.3, 2.5)
    
    cube = np.zeros((ny, nx, nlam), dtype=np.float32)
    
    for iy in range(ny):
        for ix in range(nx):
            ice_wt = ice_content_map[iy, ix]
            
            # 論文に基づく近似：バンド深さ ≈ k × 水氷量
            band_depth = 0.085 * (ice_wt / 2.2)   # 最大で約0.085程度
            
            spec = continuum.copy()
            
            # ガウス型吸収バンド
            absorption = band_depth * np.exp(-0.5 * ((lam - ice_center) / ice_width)**2)
            spec = spec * (1.0 - absorption)
            
            # 鉱物種・粒径による微小変動
            mineral_noise = 1.0 + rng.normal(0, 0.018, nlam)
            spec *= mineral_noise
            
            # 全体の明るさ変動（影・照射角効果模擬）
            brightness = 0.85 + 0.3 * rng.random()
            spec *= brightness
            
            # クリップ（物理的に妥当な反射率範囲）
            cube[iy, ix, :] = np.clip(spec, 0.05, 1.25)
    
    return cube, lam, ice_content_map


def save_lunar_montage(cube, wavelengths, out_dir="lunar_ice_mock", step=40):
    """論文向けモンタージュ保存"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ny, nx, nlam = cube.shape
    indices = list(range(0, nlam, step))
    if indices[-1] != nlam - 1:
        indices.append(nlam - 1)
    
    nplots = len(indices)
    ncols = min(4, nplots)
    nrows = int(np.ceil(nplots / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)
    axes = axes.ravel()
    
    for i, (ax, bidx) in enumerate(zip(axes, indices)):
        wl = wavelengths[bidx]
        img = cube[:, :, bidx]
        vmin, vmax = np.percentile(img, [2, 98])
        im = ax.imshow(img, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_title(f'{wl:.0f} nm')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    for ax in axes[nplots:]:
        ax.axis('off')
    
    fig.suptitle('Lunar PSR Water Ice Mock (1.5μm Absorption)', fontsize=16)
    fig.tight_layout()
    out_path = out_dir / "lunar_ice_montage.png"
    fig.savefig(out_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    
    return out_path


# =========================
# 実行例
# =========================
if __name__ == "__main__":
    cube, wavelengths, ice_map = create_lunar_ice_mock(
        nx=100,
        ny=100,
        seed=123
    )
    
    print("Cube shape :", cube.shape)
    print("Wavelength range:", wavelengths[0], "〜", wavelengths[-1], "nm")
    print("Ice content range:", f"{ice_map.min():.2f} 〜 {ice_map.max():.2f} wt%")
    
#    out_path = save_lunar_montage(cube, wavelengths)
    out_path = save_lunar_montage_fixed_scale(cube, wavelengths)
    print("Montage saved:", out_path)