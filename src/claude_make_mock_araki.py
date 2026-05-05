"""
claude_make_mock_araki.py
==========================
Araki & Saiki (2025) 実測データを使った ALIS モックハイパースペクトルデータ生成器

既存の claude_make_mock.py (簡易ガウスモデル) との違い:
  - 乾燥スペクトル : _dry.csv の実測値を読み込み・平均化して使用
  - 水氷吸収形状   : fitted ice CSV と dry CSV の差分から実測プロファイルを抽出
  - 校正線        : Table 2 の実験勾配値を使用 (Eq.1 予測値も参照として出力)

参照:
  Araki, R. and Saiki, K. (2025) Geochem. J., 59, 174-191.
  DOI: 10.2343/geochemj.GJ25010
  データ DOI: 10.60574/87068  (Osaka University OUKA, 134 CSV ファイル)

依存パッケージ: numpy, matplotlib (標準 venv)
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np


# ======================================================================
# パス設定
# ======================================================================

# このファイルの 2 階層上が analysis_hsi/ プロジェクトルート
_SRC_DIR = Path(__file__).parent
DATA_DIR = _SRC_DIR.parent / "Mock_Data" / "Experimental_Data"


# ======================================================================
# ALIS 波長グリッド (750-1650nm, 5nm)
# ======================================================================

ALIS_WL: np.ndarray = np.arange(750.0, 1651.0, 5.0)  # 181 bands


# ======================================================================
# Table S1: 実験ファイルごとの水氷含有量 (wt.%)
# Araki & Saiki (2025) Data paper Table S1
# ======================================================================

EXPERIMENT_ICE_CONTENT: dict[str, float] = {
    # Olivine 180-250um (#1-#5)
    "OC210315": 0.65,  "OC210317": 0.63,  "OC210318_1": 0.91,
    "OC210318_2": 0.86, "OC210405_2": 2.01,
    # Olivine 75-125um (#6-#11)
    "OF210319": 0.31,  "OF210323_1": 0.45, "OF210323_2": 0.65,
    "OF210324": 1.82,  "OF210407": 1.14,   "OF230125": 0.41,
    # Plagioclase 180-250um (#12-#16)
    "PC210329": 1.56,  "PC210330_1": 1.45, "PC210331": 1.22,
    "PC210402": 0.69,  "PC210405_1": 0.33,
    # Plagioclase 75-125um (#17-#22)
    "PF210325": 0.49,  "PF210326_1": 0.59, "PF210326_2": 0.49,
    "PF210327_2": 0.73, "PF210330_2": 2.02, "PF210401": 0.73,
    # CPX 180-250um (#23-#27)
    "CC211018_1": 0.86, "CC211019": 0.36,  "CC211116": 0.50,
    "CC211130": 0.38,   "CC220119": 1.04,
    # CPX 75-125um (#28-#35)
    "CF211016": 0.87,  "CF211018_2": 0.59, "CF211025": 0.35,
    "CF211026": 1.77,  "CF211101": 0.65,   "CF211102": 0.85,
    "CF211115": 0.42,  "CF211119": 1.49,
    # Mixture 180-250um (#36-#40)
    "MC220427": 1.76,  "MC220430": 0.74,   "MC220514_2": 1.23,
    "MC220516": 1.40,  "MC220518": 0.73,
    # Mixture 75-125um (#41-#46)
    "MF220428": 0.68,  "MF220509": 0.52,   "MF220514_1": 2.15,
    "MF220517": 0.83,  "MF221104": 0.81,   "MF221111": 1.67,
}

# fitted ファイル名ステム → raw ファイル名ステムへのマッピング
# (fitted と raw のファイル名が一致しないケースのみ記載)
# None = raw データが存在しない (fitted のみ)
FITTED_TO_RAW_STEM: dict[str, str | None] = {
    "OC210405_2": "OC210405",   # raw に _2 なし
    "PC210330_1": "PC210330",   # raw に _1 なし
    "PC210405_1": "PC210405",   # raw に _1 なし
    "PF210327_2": "PF210327",   # raw に _2 なし
    "PF210330_2": "PF210330",   # raw に _2 なし
    "MF221104":   None,          # raw データなし (fitted のみ)
    "MF221111":   None,          # raw データなし (fitted のみ)
}


# ======================================================================
# グループ定義 (Table 2: Araki & Saiki 2025)
# ======================================================================

@dataclass
class GroupSpec:
    """鉱物種×粒径グループの実験パラメータ"""
    mineral_folder: str   # データディレクトリ内のフォルダ名
    grain_folder: str     # 粒径フォルダ名
    grain_size_um: float  # 粒径中央値 (μm)
    reflectance_1500: float  # 乾燥鉱物の 1.5μm 反射率 (Table 2)
    gradient: float       # 校正線勾配の実験値 (Table 2)


GROUP_SPECS: dict[tuple[str, str], GroupSpec] = {
    ("olivine",       "coarse"): GroupSpec("olivine",      "180-250um", 215, 0.460, 5.18),
    ("olivine",       "fine"):   GroupSpec("olivine",      "75-125um",  100, 0.542, 1.83),
    ("plagioclase",   "coarse"): GroupSpec("plagioclase",  "180-250um", 215, 0.565, 8.48),
    ("plagioclase",   "fine"):   GroupSpec("plagioclase",  "75-125um",  100, 0.679, 5.83),
    ("clinopyroxene", "coarse"): GroupSpec("cpx",          "180-250um", 215, 0.610, 11.87),
    ("clinopyroxene", "fine"):   GroupSpec("cpx",          "75-125um",  100, 0.665, 5.91),
    ("mixture",       "coarse"): GroupSpec("mixture",      "180-250um", 215, 0.509, 5.71),
    ("mixture",       "fine"):   GroupSpec("mixture",      "75-125um",  100, 0.634, 3.11),
}


# ======================================================================
# Eq. (1): 校正線勾配の予測 (Araki & Saiki 2025)
# ======================================================================

def predict_gradient(grain_size_um: float, reflectance_1500: float) -> float:
    """
    gradient = 0.0517 × grain_size_um + 26.0 × reflectance_1500 − 17.4
    """
    return 0.0517 * grain_size_um + 26.0 * reflectance_1500 - 17.4


# ======================================================================
# CSV 読み込み・補間ユーティリティ
# ======================================================================

def _load_csv(path: Path, has_header: bool) -> tuple[np.ndarray, np.ndarray]:
    """CSV から (wavelengths, reflectances) を読み込む"""
    data = np.loadtxt(path, delimiter=",", skiprows=1 if has_header else 0)
    return data[:, 0], np.clip(data[:, 1], 0.0, None)


def interp_to_alis(wl: np.ndarray, refl: np.ndarray) -> np.ndarray:
    """
    任意の波長グリッドの反射率を ALIS グリッド (750-1650nm, 5nm) に変換する。
    - 750-800nm: データ端の傾きから線形外挿
    - 800-1640nm: 線形補間
    - 1640-1650nm: 最右端値を使用
    """
    result = np.interp(ALIS_WL, wl, refl, left=np.nan, right=np.nan)

    # 750-800nm 外挿 (wl[0]=800nm から左へ線形延長)
    mask_low = ALIS_WL < wl[0]
    if mask_low.any():
        n = min(5, len(wl))
        slope = (refl[n - 1] - refl[0]) / (wl[n - 1] - wl[0])
        result[mask_low] = refl[0] + slope * (ALIS_WL[mask_low] - wl[0])

    # 1640-1650nm 外挿 (右端を定数延長)
    mask_high = ALIS_WL > wl[-1]
    if mask_high.any():
        result[mask_high] = refl[-1]

    return np.clip(result, 0.0, None)


# ======================================================================
# グループデータの構築
# ======================================================================

def build_group_data(spec: GroupSpec, data_dir: Path) -> dict:
    """
    指定グループの代表的な乾燥スペクトルと単位吸収プロファイルを構築する。

    乾燥スペクトル : グループ内の全 _dry.csv を ALIS グリッドに補間し平均
    単位吸収プロファイル:
        各 (dry_i, fitted_ice_i) ペアで
            frac_abs_i(λ) = (R_dry_i(λ) − R_ice_i(λ)) / R_dry_i(λ)
            unit_abs_i(λ)  = frac_abs_i(λ) / expected_BD_i
            expected_BD_i  = gradient_exp × ice_content_i / 100
        の平均を取る。これにより「単位 BD あたりの吸収形状」が得られ、
        任意の水氷量に対するスペクトルを:
            R_mock(λ) = R_dry_mean(λ) × (1 − target_BD × unit_abs(λ))
        で生成できる。

    Returns
    -------
    dict with:
        "dry_mean"  : ndarray (n_alis,) 乾燥スペクトル平均
        "dry_std"   : ndarray (n_alis,) 乾燥スペクトル標準偏差
        "unit_abs"  : ndarray (n_alis,) 単位 BD あたりの吸収プロファイル
        "n_dry"     : int    使用した dry ファイル数
        "n_pairs"   : int    使用した (dry, fitted) ペア数
    """
    raw_dir    = data_dir / spec.mineral_folder / spec.grain_folder / "raw data"
    fitted_dir = data_dir / spec.mineral_folder / spec.grain_folder / "fitted data"

    # --- 乾燥スペクトルの読み込みと平均化 ---
    dry_spectra: list[np.ndarray] = []
    for f in sorted(raw_dir.glob("*_dry.csv")):
        wl, refl = _load_csv(f, has_header=False)
        dry_spectra.append(interp_to_alis(wl, refl))

    if not dry_spectra:
        raise FileNotFoundError(f"No _dry.csv found in {raw_dir}")

    dry_arr  = np.array(dry_spectra)  # (n_dry, n_alis)
    dry_mean = dry_arr.mean(axis=0)
    dry_std  = dry_arr.std(axis=0)

    # --- 単位吸収プロファイルの抽出 ---
    unit_abs_list: list[np.ndarray] = []

    for fitted_file in sorted(fitted_dir.glob("*.csv")):
        stem = fitted_file.stem

        # 対応する dry ファイルを解決
        if stem in FITTED_TO_RAW_STEM:
            raw_stem = FITTED_TO_RAW_STEM[stem]
        else:
            raw_stem = stem

        if raw_stem is None:
            dry_this = dry_mean.copy()
        else:
            dry_path = raw_dir / f"{raw_stem}_dry.csv"
            if dry_path.exists():
                wl, refl = _load_csv(dry_path, has_header=False)
                dry_this = interp_to_alis(wl, refl)
            else:
                dry_this = dry_mean.copy()

        # fitted ice スペクトルの読み込み
        wl_ice, refl_ice = _load_csv(fitted_file, has_header=True)
        ice_alis = interp_to_alis(wl_ice, refl_ice)

        # 水氷含有量と期待 BD
        ice_wt = EXPERIMENT_ICE_CONTENT.get(stem)
        if ice_wt is None or ice_wt <= 0:
            continue
        expected_bd = spec.gradient * ice_wt / 100.0
        if expected_bd <= 0:
            continue

        # 分率吸収の計算と正規化
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_abs = np.where(dry_this > 1e-4,
                                (dry_this - ice_alis) / dry_this, 0.0)
        frac_abs = np.clip(frac_abs, -0.5, 1.0)
        unit_abs_list.append(frac_abs / expected_bd)

    if unit_abs_list:
        unit_abs_mean = np.clip(np.mean(unit_abs_list, axis=0), 0.0, None)
    else:
        # フォールバック: ガウス型 (既存コードと同じ)
        unit_abs_mean = np.exp(-0.5 * ((ALIS_WL - 1500.0) / 60.0) ** 2)

    return {
        "dry_mean": dry_mean,
        "dry_std": dry_std,
        "unit_abs": unit_abs_mean,
        "n_dry": len(dry_spectra),
        "n_pairs": len(unit_abs_list),
    }


# ======================================================================
# モックキューブ設定
# ======================================================================

@dataclass
class ALISMockConfigAraki:
    """ALIS モックキューブの生成設定 (Araki 実データ版)"""
    nx:  int   = 100
    ny:  int   = 50
    mineral_type: str  = "mixture"   # "olivine","plagioclase","clinopyroxene","mixture"
    grain_size:   str  = "coarse"    # "coarse" (180-250μm) / "fine" (75-125μm)
    ice_content_min: float = 0.0     # wt.%
    ice_content_max: float = 2.2     # wt.%
    ice_pattern: Literal["gradient", "patches", "uniform"] = "gradient"
    snr:  float = 100.0
    seed: int   = 42
    add_dry_variation: bool = True   # 乾燥スペクトルにピクセル間ばらつきを加える


# ======================================================================
# 氷量マップと BD 計算
# ======================================================================

def create_ice_map(config: ALISMockConfigAraki, rng: np.random.Generator) -> np.ndarray:
    ny, nx = config.ny, config.nx

    if config.ice_pattern == "uniform":
        return np.full((ny, nx), (config.ice_content_min + config.ice_content_max) / 2.0)

    elif config.ice_pattern == "gradient":
        x_n = np.linspace(0, 1, nx)
        y_n = np.linspace(0, 1, ny)
        xx, yy = np.meshgrid(x_n, y_n)
        r = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)
        ice_map = config.ice_content_max * (1.0 - r / r.max())
        ice_map += rng.normal(0, 0.1, (ny, nx))

    elif config.ice_pattern == "patches":
        ice_map = np.full((ny, nx), config.ice_content_min)
        for _ in range(5):
            cx = rng.integers(0, nx)
            cy = rng.integers(0, ny)
            radius = rng.integers(5, max(6, min(nx, ny) // 4))
            xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
            mask = (xx - cx) ** 2 + (yy - cy) ** 2 < radius ** 2
            ice_map[mask] = rng.uniform(0.5, config.ice_content_max)
    else:
        raise ValueError(f"Unknown ice_pattern: {config.ice_pattern}")

    return np.clip(ice_map, config.ice_content_min, config.ice_content_max)


def compute_band_depth(
    spectrum:         np.ndarray,
    wavelengths:      np.ndarray,
    shoulder_left_nm: float = 1350.0,
    shoulder_right_nm: float = 1640.0,
    band_center_nm:   float = 1500.0,
) -> float:
    """BD = (Rc − Rb) / Rc 、2 点線形コンティニュアムによる計算"""
    idx_l = int(np.argmin(np.abs(wavelengths - shoulder_left_nm)))
    idx_r = int(np.argmin(np.abs(wavelengths - shoulder_right_nm)))
    idx_c = int(np.argmin(np.abs(wavelengths - band_center_nm)))
    r_l, r_r = spectrum[idx_l], spectrum[idx_r]
    frac = (wavelengths[idx_c] - wavelengths[idx_l]) / (wavelengths[idx_r] - wavelengths[idx_l])
    rc = r_l + frac * (r_r - r_l)
    if rc <= 0:
        return 0.0
    return float((rc - spectrum[idx_c]) / rc)


# ======================================================================
# モックキューブ生成メイン関数
# ======================================================================

def create_alis_mock_araki(
    config:   ALISMockConfigAraki | None = None,
    data_dir: Path = DATA_DIR,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict, dict]:
    """
    Araki & Saiki (2025) の実測データを基にした ALIS モックキューブを生成する。

    Parameters
    ----------
    config   : 生成設定 (省略時はデフォルト値)
    data_dir : 実験データのルートディレクトリ

    Returns
    -------
    cube       : ndarray (ny, nx, n_bands)  反射率キューブ
    wavelengths: ndarray (n_bands,)         波長配列 (nm)
    ice_map    : ndarray (ny, nx)           入力水氷量 (wt.%)
    bd_map     : ndarray (ny, nx)           1.5μm バンド深さマップ
    metadata   : dict                       生成パラメータと参照情報
    group_data : dict                       乾燥スペクトル等 (検証用)
    """
    if config is None:
        config = ALISMockConfigAraki()

    rng  = np.random.default_rng(config.seed)
    key  = (config.mineral_type, config.grain_size)
    if key not in GROUP_SPECS:
        raise KeyError(f"Unknown group key {key}. Valid: {list(GROUP_SPECS)}")
    spec = GROUP_SPECS[key]

    print(f"  Loading group data: {spec.mineral_folder} / {spec.grain_folder} ...", end="", flush=True)
    group_data = build_group_data(spec, data_dir)
    print(f" done ({group_data['n_dry']} dry, {group_data['n_pairs']} pairs)")

    dry_mean = group_data["dry_mean"]
    dry_std  = group_data["dry_std"]
    unit_abs = group_data["unit_abs"]
    gradient = spec.gradient
    gradient_pred = predict_gradient(spec.grain_size_um, spec.reflectance_1500)

    ice_map = create_ice_map(config, rng)
    ny, nx  = config.ny, config.nx
    nlam    = len(ALIS_WL)
    noise_sigma = spec.reflectance_1500 / config.snr if config.snr > 0 else 0.0

    cube   = np.zeros((ny, nx, nlam), dtype=np.float32)
    bd_map = np.zeros((ny, nx),       dtype=np.float32)

    for iy in range(ny):
        for ix in range(nx):
            ice_wt = ice_map[iy, ix]

            # 乾燥スペクトル (ピクセルごとのばらつきを加える場合)
            if config.add_dry_variation:
                perturb = rng.normal(0.0, 0.5, nlam) * dry_std
                dry_px  = np.clip(dry_mean + perturb, 0.01, None)
            else:
                dry_px = dry_mean.copy()

            # 水氷吸収の適用
            target_bd = gradient * ice_wt / 100.0
            spec_px   = dry_px * (1.0 - target_bd * unit_abs)
            spec_px   = np.clip(spec_px, 0.01, None)

            # 検出器ノイズ
            if noise_sigma > 0:
                spec_px += rng.normal(0.0, noise_sigma, nlam)
            spec_px = np.clip(spec_px, 0.01, 1.5).astype(np.float32)

            cube[iy, ix, :]  = spec_px
            bd_map[iy, ix]   = compute_band_depth(spec_px, ALIS_WL)

    metadata = {
        "instrument":           "ALIS (Advanced Lunar Imaging Spectrometer)",
        "data_source":          "Araki & Saiki (2025), DOI: 10.60574/87068",
        "mineral_type":         config.mineral_type,
        "grain_size":           config.grain_size,
        "grain_size_um":        spec.grain_size_um,
        "reflectance_1500":     spec.reflectance_1500,
        "gradient_experimental": gradient,
        "gradient_predicted_eq1": gradient_pred,
        "ice_content_range_wt": (config.ice_content_min, config.ice_content_max),
        "snr":                  config.snr,
        "n_bands":              nlam,
        "n_dry_files":          group_data["n_dry"],
        "n_pairs":              group_data["n_pairs"],
        "reference":            "Araki & Saiki (2025) Geochem. J., 59, 174-191",
    }
    return cube, ALIS_WL, ice_map, bd_map, metadata, group_data


# ======================================================================
# 可視化関数群
# ======================================================================

def plot_calibration(
    ice_map:  np.ndarray,
    bd_map:   np.ndarray,
    metadata: dict,
    out_path: Path,
) -> None:
    """BD vs ice content の校正線検証プロット (論文 Fig. 8 に対応)"""
    ice_flat = ice_map.ravel()
    bd_flat  = bd_map.ravel()
    gradient_exp  = metadata["gradient_experimental"]
    gradient_pred = metadata["gradient_predicted_eq1"]

    ice_range   = np.linspace(0, ice_map.max() * 1.05, 100)
    bd_exp_line = gradient_exp  * ice_range / 100.0
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
        f"R(1.5μm)={metadata['reflectance_1500']:.3f}, SNR={metadata['snr']:.0f}",
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
    cube:        np.ndarray,
    ice_map:     np.ndarray,
    bd_map:      np.ndarray,
    metadata:    dict,
    out_path:    Path,
) -> None:
    """代表バンド反射率 + 氷量マップ + BD マップのモンタージュ"""
    target_wls  = [800, 1000, 1200, 1400, 1500, 1600]
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
        f"{metadata['gradient_predicted_eq1']:.2f} (Eq.1)\n"
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
        spec = cube[iy, ix, :]
        ice_val = ice_map[iy, ix]
        ax.plot(ALIS_WL, spec, color=color, lw=1.2,
                label=f"ice={ice_val:.2f} wt.% (q={q:.0%})")

    # 実測乾燥スペクトル (平均)
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
    ax.set_xlim(750, 1650)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_dry_spectra_all_groups(
    data_dir: Path,
    out_path: Path,
) -> None:
    """全グループの平均乾燥スペクトルを 1 枚に重ね描き (Fig. 4 dry lines に対応)"""
    linestyles = ["-", "--", "-.", ":"]
    colors_min = {
        "olivine": "#2ca02c",
        "plagioclase": "#1f77b4",
        "clinopyroxene": "#d62728",
        "mixture": "#ff7f0e",
    }
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    titles = ["180-250 μm (coarse)", "75-125 μm (fine)"]
    grain_keys = ["coarse", "fine"]

    for col, (grain_key, title) in enumerate(zip(grain_keys, titles)):
        ax = axes[col]
        for mineral in ["olivine", "plagioclase", "clinopyroxene", "mixture"]:
            key = (mineral, grain_key)
            spec = GROUP_SPECS[key]
            try:
                gd = build_group_data(spec, data_dir)
            except FileNotFoundError:
                continue
            color = colors_min[mineral]
            # 個別スペクトルを薄く描画
            raw_dir = data_dir / spec.mineral_folder / spec.grain_folder / "raw data"
            for f in sorted(raw_dir.glob("*_dry.csv")):
                wl, refl = _load_csv(f, has_header=False)
                ax.plot(ALIS_WL, interp_to_alis(wl, refl),
                        color=color, alpha=0.15, lw=0.8)
            # 平均を太線で描画
            ax.plot(ALIS_WL, gd["dry_mean"],
                    color=color, lw=2.0, label=mineral)

        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        ax.set_ylabel("Reflectance", fontsize=11)
        ax.set_xlim(750, 1650)
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
    titles = ["180-250 μm (coarse)", "75-125 μm (fine)"]
    grain_keys = ["coarse", "fine"]
    colors = {
        "olivine": "#2ca02c",
        "plagioclase": "#1f77b4",
        "clinopyroxene": "#d62728",
        "mixture": "#ff7f0e",
    }

    for col, (grain_key, title) in enumerate(zip(grain_keys, titles)):
        ax = axes[col]
        for mineral in ["olivine", "plagioclase", "clinopyroxene", "mixture"]:
            key = (mineral, grain_key)
            spec = GROUP_SPECS[key]
            try:
                gd = build_group_data(spec, data_dir)
            except FileNotFoundError:
                continue
            ax.plot(ALIS_WL, gd["unit_abs"],
                    color=colors[mineral], lw=2.0, label=mineral)

        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Wavelength (nm)", fontsize=11)
        ax.set_ylabel("Unit absorption (per BD)", fontsize=11)
        ax.set_xlim(1200, 1650)
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


# ======================================================================
# メイン実行
# ======================================================================

if __name__ == "__main__":
    out_dir = Path(__file__).parent / "alis_mock_output_araki"
    out_dir.mkdir(exist_ok=True)

    # --- 全グループの乾燥スペクトルと吸収プロファイルの可視化 ---
    print("Plotting dry spectra (all groups)...")
    plot_dry_spectra_all_groups(DATA_DIR, out_dir / "dry_spectra_all_groups.png")

    print("Plotting unit absorption profiles (all groups)...")
    plot_unit_absorption_profiles(DATA_DIR, out_dir / "unit_absorption_profiles.png")

    # --- 各設定でモックキューブ生成 ---
    configs = [
        ALISMockConfigAraki(mineral_type="olivine",       grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="olivine",       grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="plagioclase",   grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="plagioclase",   grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="clinopyroxene", grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="clinopyroxene", grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfigAraki(mineral_type="mixture",       grain_size="coarse",
                            ice_pattern="patches",  snr=80),
        ALISMockConfigAraki(mineral_type="mixture",       grain_size="fine",
                            ice_pattern="patches",  snr=80),
    ]

    for cfg in configs:
        tag = f"{cfg.mineral_type}_{cfg.grain_size}"
        print(f"\n{'='*60}")
        print(f"Generating: {cfg.mineral_type} ({cfg.grain_size})")
        print(f"{'='*60}")

        cube, wl, ice_map, bd_map, meta, gd = create_alis_mock_araki(cfg, DATA_DIR)

        print(f"  Cube shape           : {cube.shape}")
        print(f"  Gradient (exp / Eq1) : {meta['gradient_experimental']:.2f} / "
              f"{meta['gradient_predicted_eq1']:.2f}")
        print(f"  R(1.5μm)             : {meta['reflectance_1500']:.3f}")
        print(f"  Ice range            : {ice_map.min():.2f}–{ice_map.max():.2f} wt.%")
        print(f"  BD range             : {bd_map.min():.4f}–{bd_map.max():.4f}")
        print(f"  Dry files used       : {meta['n_dry_files']} dry, {meta['n_pairs']} pairs")

        plot_calibration(ice_map, bd_map, meta, out_dir / f"{tag}_calibration.png")
        print(f"  calibration.png -> {out_dir / f'{tag}_calibration.png'}")

        plot_montage(cube, ice_map, bd_map, meta, out_dir / f"{tag}_montage.png")
        print(f"  montage.png     -> {out_dir / f'{tag}_montage.png'}")

        plot_spectra_comparison(cube, ice_map, gd, meta,
                                out_dir / f"{tag}_spectra.png")
        print(f"  spectra.png     -> {out_dir / f'{tag}_spectra.png'}")

    print(f"\nAll outputs saved to: {out_dir}/")
