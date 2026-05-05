"""
ALIS モック ハイパースペクトルデータ生成器
=============================================

LUPEX 搭載 ALIS (Advanced Lunar Imaging Spectrometer) の仕様と
Araki & Saiki (2025, Geochemical Journal) の校正線モデルに基づき、
月極域 PSR の水氷混合レゴリスのモックハイパースペクトルキューブを生成する。

ALIS 仕様（佐伯ほか, 2025 宇宙科学技術連合講演会 3G08）:
  - 観測波長: 750–1650 nm
  - 波長分解能: < 5 nm
  - 検出器: SONY IMX990
  - iFOV: 2 m × 1 cm @ 5 m（水平視野 22°）

Araki & Saiki (2025) 校正線モデル:
  - BD = gradient × ice_content  (原点通過の線形フィット)
  - gradient = 0.0517 * grain_size_um + 26.0 * reflectance_1500 - 17.4  (Eq. 1)
  - BD = (Rc - Rb) / Rc  (Clark & Roush, 1984)

参照データ:
  - Table 2: 鉱物種ごとの粒径・反射率・校正線勾配
  - 校正線勾配の実験値範囲: 1.83–11.87
  - 水氷量範囲: 0.3–2.2 wt.%
"""

import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from typing import Literal


# ======================================================================
# 鉱物パラメータ定義 (Araki & Saiki 2025, Table 2)
# ======================================================================

@dataclass
class MineralParams:
    """鉱物種ごとのスペクトル・校正パラメータ"""
    name: str
    reflectance_1500_coarse: float   # 180–250 μm での 1.5 μm 反射率
    reflectance_1500_fine: float     # 75–125 μm での 1.5 μm 反射率
    gradient_coarse: float           # 校正線勾配（実測, 180–250 μm）
    gradient_fine: float             # 校正線勾配（実測, 75–125 μm）
    # 1 μm 付近の吸収帯パラメータ（簡易モデル）
    abs_1um_center: float            # nm
    abs_1um_width: float             # nm（半値幅相当）
    abs_1um_depth: float             # 吸収深さ（0–1）


# Table 2 の値を使用
MINERALS = {
    "olivine": MineralParams(
        name="Olivine (Fo90)",
        reflectance_1500_coarse=0.460,
        reflectance_1500_fine=0.542,
        gradient_coarse=5.18,
        gradient_fine=1.83,
        abs_1um_center=1050,    # olivine の ~1 μm 吸収帯
        abs_1um_width=200,
        abs_1um_depth=0.25,
    ),
    "plagioclase": MineralParams(
        name="Plagioclase (An60)",
        reflectance_1500_coarse=0.565,
        reflectance_1500_fine=0.679,
        gradient_coarse=8.48,
        gradient_fine=5.83,
        abs_1um_center=1200,    # plagioclase の浅く広い吸収帯
        abs_1um_width=300,
        abs_1um_depth=0.08,
    ),
    "clinopyroxene": MineralParams(
        name="Clinopyroxene (Diopside)",
        reflectance_1500_coarse=0.610,
        reflectance_1500_fine=0.665,
        gradient_coarse=11.87,
        gradient_fine=5.91,
        abs_1um_center=1000,    # CPX の ~1 μm 吸収帯
        abs_1um_width=150,
        abs_1um_depth=0.30,
    ),
    "mixture": MineralParams(
        name="Three-mineral mixture",
        reflectance_1500_coarse=0.509,
        reflectance_1500_fine=0.634,
        gradient_coarse=5.71,
        gradient_fine=3.11,
        abs_1um_center=1050,    # 混合のため olivine 寄り
        abs_1um_width=250,
        abs_1um_depth=0.18,
    ),
}


# ======================================================================
# 定数クラス
# ======================================================================

class ALISConstants:
    """ALIS 仕様関連定数（Araki & Saiki 2025, 佐伯ほか 2025）
    
    観測波長: 750–1650 nm、波長分解能: < 5 nm に準拠。
    """
    WAVELENGTH_MIN = 750.0    # nm
    WAVELENGTH_MAX = 1650.0   # nm
    SPECTRAL_RESOLUTION = 5.0 # nm


class IceAbsorptionConstants:
    """水氷吸収関連定数
    
    1.5 μm 吸収帯（Fig. 4: 約 1350–1600 nm）をガウス型でモデル化。
    """
    ABSORPTION_CENTER = 1500.0  # nm
    ABSORPTION_SIGMA = 60.0     # nm（半値幅 ~140 nm に対応）


def alis_wavelengths() -> np.ndarray:
    """ALIS の波長軸を生成する (5 nm 刻み, 750–1650 nm)"""
    return np.arange(ALISConstants.WAVELENGTH_MIN, ALISConstants.WAVELENGTH_MAX + 0.1,
                     ALISConstants.SPECTRAL_RESOLUTION)


# ======================================================================
# Eq. (1): 校正線勾配の予測
# ======================================================================

def predict_gradient(grain_size_um: float, reflectance_1500: float) -> float:
    """
    Araki & Saiki (2025) Eq. (1) に基づく校正線勾配の予測。

    Parameters
    ----------
    grain_size_um : float
        レゴリス粒径の中央値 (μm)
    reflectance_1500 : float
        乾燥鉱物粉末の 1.5 μm における反射率

    Returns
    -------
    float
        校正線勾配 (BD / ice_content)
    """
    return 0.0517 * grain_size_um + 26.0 * reflectance_1500 - 17.4


# ======================================================================
# 乾燥鉱物のベーススペクトル生成（簡易モデル）
# ======================================================================

def dry_mineral_spectrum(
    wavelengths: np.ndarray,
    mineral: MineralParams,
    grain_size_um: float,
) -> np.ndarray:
    """
    乾燥鉱物粉末の反射率スペクトルを簡易生成する。

    論文 Fig. 4 の乾燥試料スペクトルの定性的特徴を再現:
    - 1.5 μm 付近の反射率レベルを Table 2 の値に合わせる
    - 鉱物種固有の 1 μm 付近吸収帯を含む
    - 750–900 nm は反射率がやや低い傾向

    注意: 実データの代用として使用する簡易モデルであり、
    精密な再現には Araki & Saiki の公開データ (DOI: 10.60574/87068)
    を読み込むべきである。

    Parameters
    ----------
    wavelengths : array
        波長配列 (nm)
    mineral : MineralParams
        鉱物パラメータ
    grain_size_um : float
        粒径中央値 (μm)。粗粒(215)か細粒(100)かで反射率を内挿。

    Returns
    -------
    array
        反射率スペクトル
    """
    # 粒径に応じた 1.5 μm 反射率の線形内挿
    t = (grain_size_um - 100.0) / (215.0 - 100.0)
    t = np.clip(t, 0.0, 1.0)
    ref_1500 = (1 - t) * mineral.reflectance_1500_fine + \
               t * mineral.reflectance_1500_coarse

    # コンティニュアム: 750 nm 付近で ref_1500 よりやや低く、
    # 1500 nm 付近で ref_1500 に達し、1650 nm でほぼ同レベル
    # → 短波長側が低い傾斜を持つ
    continuum = ref_1500 * (0.85 + 0.15 * (wavelengths - 750) / (1500 - 750))
    continuum = np.clip(continuum, 0.0, None)

    # 鉱物固有の 1 μm 付近吸収帯（ガウス型）
    abs_1um = mineral.abs_1um_depth * np.exp(
        -0.5 * ((wavelengths - mineral.abs_1um_center) / mineral.abs_1um_width) ** 2
    )
    spectrum = continuum * (1.0 - abs_1um)

    return spectrum


# ======================================================================
# 水氷吸収の適用と BD 計算
# ======================================================================

def apply_ice_absorption(
    dry_spectrum: np.ndarray,
    wavelengths: np.ndarray,
    ice_content_wt: float,
    gradient: float,
) -> np.ndarray:
    """
    乾燥スペクトルに水氷吸収を適用する。

    校正線モデル BD = gradient × ice_content から目標 BD を算出し、
    IceAbsorptionConstants で定義されたガウス型吸収プロファイルを適用する。

    BD の定義（Clark & Roush, 1984）: BD = (Rc - Rb) / Rc
    吸収バンド中心での減衰を表現する。

    Parameters
    ----------
    dry_spectrum : array
        乾燥鉱物の反射率スペクトル
    wavelengths : array
        波長配列 (nm)
    ice_content_wt : float
        水氷含有量 (wt.%)
    gradient : float
        校正線勾配 (Eq. 1 による予測値 or 実測値)

    Returns
    -------
    array
        水氷吸収適用後の反射率スペクトル
    """
    target_bd = gradient * ice_content_wt / 100.0
    # BD は通常 0–1 の範囲。論文の実験では最大 ~0.14 程度
    target_bd = np.clip(target_bd, 0.0, 0.5)

    # ガウス型吸収プロファイル（IceAbsorptionConstants から中心・幅を取得）
    gauss = np.exp(
        -0.5 * ((wavelengths - IceAbsorptionConstants.ABSORPTION_CENTER) / IceAbsorptionConstants.ABSORPTION_SIGMA) ** 2
    )
    absorption_profile = target_bd * gauss

    # スペクトルへの適用: R_ice = R_dry * (1 - absorption_profile)
    wet_spectrum = dry_spectrum * (1.0 - absorption_profile)

    return wet_spectrum


def compute_band_depth(
    spectrum: np.ndarray,
    wavelengths: np.ndarray,
    shoulder_left_nm: float = 1350.0,
    shoulder_right_nm: float = 1650.0,
    band_center_nm: float = 1500.0,
) -> float:
    """
    コンティニュアム除去による 1.5 μm バンド深さを計算する。

    簡易実装: 左右の肩の反射率から線形コンティニュアムを定義し、
    吸収帯中心の反射率との比から BD を算出。

    BD = (Rc - Rb) / Rc

    Parameters
    ----------
    spectrum : array
        反射率スペクトル
    wavelengths : array
        波長配列 (nm)
    shoulder_left_nm, shoulder_right_nm : float
        コンティニュアムの左右肩波長 (nm)
    band_center_nm : float
        吸収帯中心波長 (nm)

    Returns
    -------
    float
        バンド深さ BD
    """
    idx_left = np.argmin(np.abs(wavelengths - shoulder_left_nm))
    idx_right = np.argmin(np.abs(wavelengths - shoulder_right_nm))
    idx_center = np.argmin(np.abs(wavelengths - band_center_nm))

    r_left = spectrum[idx_left]
    r_right = spectrum[idx_right]

    # 線形コンティニュアム: 左肩→右肩を直線で結ぶ
    frac = (wavelengths[idx_center] - wavelengths[idx_left]) / \
           (wavelengths[idx_right] - wavelengths[idx_left])
    rc = r_left + frac * (r_right - r_left)
    rb = spectrum[idx_center]

    if rc <= 0:
        return 0.0
    return (rc - rb) / rc


# ======================================================================
# モックキューブ生成
# ======================================================================

@dataclass
class ALISMockConfig:
    """ALIS モックキューブの生成設定
    
    Araki & Saiki (2025) の実験条件に基づくパラメータ範囲を採用。
    """
    # 空間グリッド
    nx: int = 100                # 空間方向ピクセル数（スキャン方向）
    ny: int = 50                 # 空間方向ピクセル数（スリット方向）

    # 鉱物種マップの設定
    mineral_type: str = "mixture"  # "olivine", "plagioclase", "clinopyroxene", "mixture"
    grain_size_um: float = 150.0   # 粒径中央値 (μm)、75–250 の範囲

    # 水氷分布の設定
    ice_content_min: float = 0.0   # wt.%
    ice_content_max: float = 2.2   # wt.%
    ice_pattern: Literal["gradient", "patches", "uniform"] = "gradient"

    # ノイズ
    snr: float = 100.0             # 信号対雑音比（1.5 μm 付近での目標値）
    seed: int = 42


def create_ice_map(config: ALISMockConfig, rng: np.random.Generator) -> np.ndarray:
    """水氷含有量の空間分布マップを生成する (wt.%)"""
    ny, nx = config.ny, config.nx

    if config.ice_pattern == "uniform":
        ice_map = np.full((ny, nx),
                          (config.ice_content_min + config.ice_content_max) / 2)

    elif config.ice_pattern == "gradient":
        # PSR の縁から中心に向かって氷量が増加するモデル
        x_norm = np.linspace(0, 1, nx)
        y_norm = np.linspace(0, 1, ny)
        xx, yy = np.meshgrid(x_norm, y_norm)
        # 中心からの距離で氷量を決定
        r = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)
        r_norm = r / r.max()
        ice_map = config.ice_content_max * (1.0 - r_norm)
        ice_map += rng.normal(0, 0.1, (ny, nx))

    elif config.ice_pattern == "patches":
        # パッチ状の氷分布
        ice_map = np.full((ny, nx), config.ice_content_min)
        n_patches = 5
        for _ in range(n_patches):
            cx = rng.integers(0, nx)
            cy = rng.integers(0, ny)
            radius = rng.integers(5, min(nx, ny) // 4)
            xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
            mask = ((xx - cx) ** 2 + (yy - cy) ** 2) < radius ** 2
            patch_ice = rng.uniform(0.5, config.ice_content_max)
            ice_map[mask] = patch_ice
    else:
        raise ValueError(f"Unknown ice_pattern: {config.ice_pattern}")

    return np.clip(ice_map, config.ice_content_min, config.ice_content_max)


def create_alis_mock(config: ALISMockConfig | None = None):
    """
    ALIS 仕様に基づくモックハイパースペクトルキューブを生成する。

    Returns
    -------
    cube : ndarray, shape (ny, nx, nlam)
        反射率キューブ
    wavelengths : ndarray, shape (nlam,)
        波長配列 (nm)
    ice_map : ndarray, shape (ny, nx)
        入力した水氷含有量マップ (wt.%)
    bd_map : ndarray, shape (ny, nx)
        計算した 1.5 μm バンド深さマップ
    metadata : dict
        生成パラメータ
    """
    if config is None:
        config = ALISMockConfig()

    rng = np.random.default_rng(config.seed)
    wavelengths = alis_wavelengths()
    nlam = len(wavelengths)
    ny, nx = config.ny, config.nx

    mineral = MINERALS[config.mineral_type]

    # Eq. (1) による校正線勾配の予測
    # 粒径に応じた反射率を内挿
    t = (config.grain_size_um - 100.0) / (215.0 - 100.0)
    t = np.clip(t, 0.0, 1.0)
    ref_1500 = (1 - t) * mineral.reflectance_1500_fine + \
               t * mineral.reflectance_1500_coarse
    gradient_predicted = predict_gradient(config.grain_size_um, ref_1500)

    # 乾燥スペクトル（全ピクセル共通のベース）
    dry_spec = dry_mineral_spectrum(wavelengths, mineral, config.grain_size_um)

    # 水氷分布マップ
    ice_map = create_ice_map(config, rng)

    # キューブ・BD マップの初期化
    cube = np.zeros((ny, nx, nlam), dtype=np.float32)
    bd_map = np.zeros((ny, nx), dtype=np.float32)

    for iy in range(ny):
        for ix in range(nx):
            ice_wt = ice_map[iy, ix]

            # 水氷吸収の適用
            spec = apply_ice_absorption(
                dry_spec, wavelengths, ice_wt, gradient_predicted
            )

            # 検出器ノイズの付加
            if config.snr > 0:
                noise_sigma = ref_1500 / config.snr
                spec = spec + rng.normal(0, noise_sigma, nlam)

            spec = np.clip(spec, 0.01, 1.5).astype(np.float32)
            cube[iy, ix, :] = spec

            # BD の計算（ノイズ込みスペクトルから復元）
            bd_map[iy, ix] = compute_band_depth(spec, wavelengths)

    metadata = {
        "instrument": "ALIS (Advanced Lunar Imaging Spectrometer)",
        "wavelength_range_nm": (ALISConstants.WAVELENGTH_MIN, ALISConstants.WAVELENGTH_MAX),
        "spectral_resolution_nm": ALISConstants.SPECTRAL_RESOLUTION,
        "n_bands": nlam,
        "mineral_type": mineral.name,
        "grain_size_um": config.grain_size_um,
        "reflectance_1500": ref_1500,
        "gradient_predicted_eq1": gradient_predicted,
        "ice_content_range_wt": (config.ice_content_min, config.ice_content_max),
        "snr": config.snr,
        "reference": "Araki & Saiki (2025) Geochem. J., 59, 174-191",
    }

    return cube, wavelengths, ice_map, bd_map, metadata


# ======================================================================
# 検証: BD vs ice content の校正線プロット
# ======================================================================

def validate_calibration(
    ice_map: np.ndarray,
    bd_map: np.ndarray,
    metadata: dict,
    out_path: str | Path = "calibration_validation.png",
) -> Path:
    """
    生成データの BD vs ice content をプロットし、
    論文の校正線（Fig. 8）と比較する。
    """
    out_path = Path(out_path)

    ice_flat = ice_map.ravel()
    bd_flat = bd_map.ravel()

    # 論文の校正線（予測勾配）
    gradient = metadata["gradient_predicted_eq1"]
    ice_range = np.linspace(0, ice_map.max(), 100)
    bd_theory = gradient * ice_range / 100.0

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # 散布図（サンプリング）
    n_samples = min(2000, len(ice_flat))
    idx = np.random.default_rng(0).choice(len(ice_flat), n_samples, replace=False)
    ax.scatter(ice_flat[idx], bd_flat[idx], s=4, alpha=0.3, color="steelblue",
               label="Mock data (from cube)")

    # 理論校正線
    ax.plot(ice_range, bd_theory, "r-", linewidth=2,
            label=f"Eq. (1) prediction (gradient={gradient:.2f})")

    ax.set_xlabel("Water ice content (wt.%)", fontsize=12)
    ax.set_ylabel("1.5-μm Band Depth", fontsize=12)
    ax.set_title(
        f"Calibration validation: {metadata['mineral_type']}, "
        f"grain size = {metadata['grain_size_um']:.0f} μm\n"
        f"R(1.5 μm) = {metadata['reflectance_1500']:.3f}, "
        f"SNR = {metadata['snr']:.0f}",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=-0.01)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ======================================================================
# モンタージュ可視化
# ======================================================================

def save_montage(
    cube: np.ndarray,
    wavelengths: np.ndarray,
    ice_map: np.ndarray,
    bd_map: np.ndarray,
    metadata: dict,
    out_path: str | Path = "alis_mock_montage.png",
) -> Path:
    """モンタージュ図（代表バンド画像 + 氷量マップ + BD マップ）を保存"""
    out_path = Path(out_path)

    # 代表波長の選択
    target_wls = [800, 1000, 1200, 1400, 1500, 1600]
    band_indices = [np.argmin(np.abs(wavelengths - wl)) for wl in target_wls]

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.ravel()

    # 反射率画像
    vmin_r, vmax_r = np.percentile(cube, [2, 98])
    for i, bidx in enumerate(band_indices):
        wl = wavelengths[bidx]
        im = axes[i].imshow(cube[:, :, bidx], origin="lower", cmap="viridis",
                            vmin=vmin_r, vmax=vmax_r)
        axes[i].set_title(f"{wl:.0f} nm", fontsize=11)
        fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04,
                     label="Reflectance")

    # 氷量マップ
    im_ice = axes[6].imshow(ice_map, origin="lower", cmap="Blues")
    axes[6].set_title("Ice content (wt.%)", fontsize=11)
    fig.colorbar(im_ice, ax=axes[6], fraction=0.046, pad=0.04,
                 label="wt.%")

    # BD マップ
    im_bd = axes[7].imshow(bd_map, origin="lower", cmap="magma")
    axes[7].set_title("1.5-μm Band Depth", fontsize=11)
    fig.colorbar(im_bd, ax=axes[7], fraction=0.046, pad=0.04,
                 label="BD")

    fig.suptitle(
        f"ALIS Mock: {metadata['mineral_type']} | "
        f"grain={metadata['grain_size_um']:.0f} μm | "
        f"gradient={metadata['gradient_predicted_eq1']:.2f}\n"
        f"Ref: {metadata['reference']}",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ======================================================================
# スペクトル比較プロット
# ======================================================================

def plot_spectra_comparison(
    cube: np.ndarray,
    wavelengths: np.ndarray,
    ice_map: np.ndarray,
    metadata: dict,
    out_path: str | Path = "spectra_comparison.png",
) -> Path:
    """氷量の異なるピクセルのスペクトルを比較する（論文 Fig. 4 に対応）"""
    out_path = Path(out_path)

    # 氷量で4分位のピクセルを選択
    ice_flat = ice_map.ravel()
    quantiles = [0.05, 0.25, 0.50, 0.95]
    thresholds = np.quantile(ice_flat, quantiles)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    colors = ["#2166ac", "#4393c3", "#f4a582", "#b2182b"]

    for q, thresh, color in zip(quantiles, thresholds, colors):
        # 該当する氷量に最も近いピクセルを探す
        idx = np.argmin(np.abs(ice_flat - thresh))
        iy, ix = np.unravel_index(idx, ice_map.shape)
        spec = cube[iy, ix, :]
        ice_val = ice_map[iy, ix]
        ax.plot(wavelengths, spec, color=color, linewidth=1.2,
                label=f"ice = {ice_val:.2f} wt.% (q={q:.0%})")

    # 乾燥スペクトル（参考）
    mineral = MINERALS[metadata["mineral_type"].split("(")[0].strip().lower()
                       if "(" in metadata["mineral_type"]
                       else list(MINERALS.keys())[
                           [m.name for m in MINERALS.values()].index(
                               metadata["mineral_type"])]]
    dry_spec = dry_mineral_spectrum(wavelengths, mineral, metadata["grain_size_um"])
    ax.plot(wavelengths, dry_spec, "k--", linewidth=1.5, label="Dry (no ice)")

    ax.set_xlabel("Wavelength (nm)", fontsize=12)
    ax.set_ylabel("Reflectance", fontsize=12)
    ax.set_title(
        f"Spectral comparison: {metadata['mineral_type']} "
        f"(cf. Araki & Saiki 2025, Fig. 4)",
        fontsize=12,
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(750, 1650)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ======================================================================
# メイン実行
# ======================================================================

if __name__ == "__main__":
    # --- 全鉱物種での生成・検証 ---
    out_dir = Path("alis_mock_output")
    out_dir.mkdir(exist_ok=True)

    configs = [
        ALISMockConfig(mineral_type="olivine", grain_size_um=215,
                       ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="plagioclase", grain_size_um=215,
                       ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="clinopyroxene", grain_size_um=215,
                       ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="mixture", grain_size_um=150,
                       ice_pattern="patches", snr=80),
    ]

    for cfg in configs:
        print(f"\n{'='*60}")
        print(f"Generating: {cfg.mineral_type}, grain={cfg.grain_size_um} μm")
        print(f"{'='*60}")

        cube, wl, ice_map, bd_map, meta = create_alis_mock(cfg)

        print(f"  Cube shape       : {cube.shape}")
        print(f"  Wavelength range : {wl[0]:.0f}–{wl[-1]:.0f} nm "
              f"({len(wl)} bands)")
        print(f"  Gradient (Eq. 1) : {meta['gradient_predicted_eq1']:.2f}")
        print(f"  R(1.5 μm)        : {meta['reflectance_1500']:.3f}")
        print(f"  Ice range        : {ice_map.min():.2f}–{ice_map.max():.2f} wt.%")
        print(f"  BD range         : {bd_map.min():.4f}–{bd_map.max():.4f}")

        prefix = f"{cfg.mineral_type}_{cfg.grain_size_um:.0f}um"

        p1 = validate_calibration(
            ice_map, bd_map, meta,
            out_dir / f"{prefix}_calibration.png")
        print(f"  Calibration plot : {p1}")

        p2 = save_montage(
            cube, wl, ice_map, bd_map, meta,
            out_dir / f"{prefix}_montage.png")
        print(f"  Montage          : {p2}")

    # mixture の詳細スペクトル比較
    cfg_mix = ALISMockConfig(mineral_type="mixture", grain_size_um=150,
                             ice_pattern="gradient", snr=200)
    cube, wl, ice_map, bd_map, meta = create_alis_mock(cfg_mix)
    p3 = plot_spectra_comparison(
        cube, wl, ice_map, meta,
        out_dir / "mixture_spectra_comparison.png")
    print(f"\n  Spectra comparison: {p3}")

    print(f"\nAll outputs saved to: {out_dir}/")