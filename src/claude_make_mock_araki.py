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

依存パッケージ: numpy
可視化は plot_araki.py を参照。

Docstring 規約:
  主要関数 (build_group_data, create_alis_mock_araki) は NumPy style で記述。
  内部ユーティリティ (_load_csv 等) は一行 docstring ＋インラインコメント。
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np


# ======================================================================
# パス設定
# ======================================================================

_SRC_DIR = Path(__file__).parent
DATA_DIR = _SRC_DIR.parent / "Mock_Data" / "Experimental_Data"


# ======================================================================
# 波長グリッド
# Araki & Saiki (2025) 実験データのカバー範囲 900–1640nm に合わせた 5nm グリッド (149 bands)
# ALIS 機器仕様 (750–1650nm) のうち実測根拠のある領域のみを使用する
# ======================================================================

ALIS_WL: np.ndarray = np.arange(900.0, 1641.0, 5.0)  # 149 bands


# ======================================================================
# 物理定数・ノイズパラメータ
# ======================================================================

ICE_MAP_NOISE_STD: float = 0.1       # 氷量マップに加算するガウスノイズの標準偏差 (wt.%)
DRY_SPEC_PERTURB_SCALE: float = 0.5  # 乾燥スペクトルのピクセル間ばらつき倍率
REFL_MIN: float = 0.01               # 反射率の物理的下限 (負値・ゼロ近傍を除外)
REFL_MAX: float = 1.5                # 反射率の物理的上限 (理想 Lambertian を超える値を除外)


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
    mineral_folder: str      # データディレクトリ内のフォルダ名
    grain_folder: str        # 粒径フォルダ名
    grain_size_um: float     # 粒径中央値 (μm)
    reflectance_1500: float  # 乾燥鉱物の 1.5μm 反射率 (Table 2)
    gradient: float          # 校正線勾配の実験値 (Table 2)


# TODO: エントリ増加時は mineral_type / grain_size を Literal または Enum で制約する
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
    """gradient = 0.0517 × grain_size_um + 26.0 × reflectance_1500 − 17.4"""
    return 0.0517 * grain_size_um + 26.0 * reflectance_1500 - 17.4


# ======================================================================
# CSV 読み込み・補間ユーティリティ
# ======================================================================

def _load_csv(path: Path, has_header: bool) -> tuple[np.ndarray, np.ndarray]:
    """CSV から (wavelengths, reflectances) を読み込む。
    has_header=True: fitted CSV（1行目にカラム名あり）、False: raw CSV（ヘッダなし）。
    反射率を np.clip(..., 0.0, None) で非負化するのは、負値が物理的に無意味なため。
    """
    data = np.loadtxt(path, delimiter=",", skiprows=1 if has_header else 0)
    return data[:, 0], np.clip(data[:, 1], 0.0, None)


def interp_to_alis(wl: np.ndarray, refl: np.ndarray) -> np.ndarray:
    """
    任意の波長グリッドの反射率を ALIS_WL (900–1640nm, 5nm) に補間する。
    データが ALIS_WL より狭い場合:
      - 左端より短波長側: データ端の傾きから線形外挿
      - 右端より長波長側: 最右端値で定数延長

    Notes
    -----
    短波長側に線形外挿を、長波長側に定数延長を採用するのは、
    実験データの短波長端は単調な傾きを持ち外挿精度が高いのに対し、
    長波長端は ALIS 帯域端 (1640nm) と一致するため傾きの推定根拠が乏しいためである。
    """
    result = np.interp(ALIS_WL, wl, refl, left=np.nan, right=np.nan)

    mask_low = ALIS_WL < wl[0]
    if mask_low.any():
        n = min(5, len(wl))
        slope = (refl[n - 1] - refl[0]) / (wl[n - 1] - wl[0])
        result[mask_low] = refl[0] + slope * (ALIS_WL[mask_low] - wl[0])

    mask_high = ALIS_WL > wl[-1]
    if mask_high.any():
        result[mask_high] = refl[-1]

    return np.clip(result, 0.0, None)


# ======================================================================
# グループデータの構築
# ======================================================================

def build_group_data(spec: GroupSpec, data_dir: Path) -> dict[str, Any]:
    """
    指定グループの代表的な乾燥スペクトルと単位吸収プロファイルを構築する。

    乾燥スペクトル : グループ内の全 _dry.csv を ALIS_WL に補間し平均
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

        wl_ice, refl_ice = _load_csv(fitted_file, has_header=True)
        ice_alis = interp_to_alis(wl_ice, refl_ice)

        ice_wt = EXPERIMENT_ICE_CONTENT.get(stem)
        if ice_wt is None or ice_wt <= 0:
            continue
        # 論文の検量線は BD = gradient × f（f = ice_wt/100, 質量分率）で定義されている。
        expected_bd = spec.gradient * ice_wt / 100.0
        if expected_bd <= 0:
            continue

        with np.errstate(divide="ignore", invalid="ignore"):
            frac_abs = np.where(dry_this > 1e-4,
                                (dry_this - ice_alis) / dry_this, 0.0)
        frac_abs = np.clip(frac_abs, -0.5, 1.0)
        unit_abs_list.append(frac_abs / expected_bd)

    if unit_abs_list:
        unit_abs_mean = np.clip(np.mean(unit_abs_list, axis=0), 0.0, None)
    else:
        unit_abs_mean = np.exp(-0.5 * ((ALIS_WL - 1500.0) / 60.0) ** 2)

    # unit_abs の再正規化:
    # frac_abs は乾燥スペクトル基準の分率吸収だが BD はコンティニュアム基準のため、
    # 小さな吸収での線形近似から peak_bd を見積もって正規化する。
    # compute_band_depth は rc <= 0 のとき 0.0 を返し、alpha は非ゼロ定数のため
    # この計算で例外は発生しない。peak_bd == 0 のケースは直下の if 文で処理する。
    alpha = 0.035
    bd_alpha = compute_band_depth(1.0 - alpha * unit_abs_mean, ALIS_WL)
    peak_bd = bd_alpha / alpha

    if peak_bd > 0:
        unit_abs_mean = unit_abs_mean / peak_bd
    else:
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
    """
    氷含有量の 2 次元空間分布マップを生成する。

    Parameters
    ----------
    config : ALISMockConfigAraki
        生成設定。ice_content_min/max とパターン種別を参照する。
    rng : np.random.Generator
        乱数生成器。

    Returns
    -------
    np.ndarray, shape (ny, nx)
        各ピクセルの水氷含有量 (wt.%)。値は [ice_content_min, ice_content_max] にクリップされる。

    Notes
    -----
    パターン:
      - "uniform"  : 全ピクセル一定値（min と max の中間値）。空間一様な実験室環境を模擬。
      - "gradient" : 中心からの放射距離に比例した減衰＋ガウスノイズ。PSR 周辺の
                     氷濃度が中心ほど高く周辺に向かって漸減する地形を模擬。
      - "patches"  : 背景値（ice_content_min）に 5 個のランダム円形高濃度領域を重畳。
                     局所的な氷露出パッチ（クレーター底や露頭）を模擬。
    """
    ny, nx = config.ny, config.nx

    if config.ice_pattern == "uniform":
        return np.full((ny, nx), (config.ice_content_min + config.ice_content_max) / 2.0)

    elif config.ice_pattern == "gradient":
        x_n = np.linspace(0, 1, nx)
        y_n = np.linspace(0, 1, ny)
        xx, yy = np.meshgrid(x_n, y_n)
        r = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)
        ice_map = config.ice_content_max * (1.0 - r / r.max())
        ice_map += rng.normal(0, ICE_MAP_NOISE_STD, (ny, nx))

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
    spectrum:          np.ndarray,
    wavelengths:       np.ndarray,
    shoulder_left_nm:  float = 1350.0,
    shoulder_right_nm: float = 1640.0,
    band_center_nm:    float = 1500.0,
) -> float:
    """BD = (Rc − Rb) / Rc、2点線形コンティニュアムによる計算。
    デフォルト値 (1350nm, 1640nm, 1500nm) は Araki & Saiki (2025) の 1.5μm 水氷吸収帯定義に基づく。
    """
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
    cube       : ndarray (ny, nx, n_bands)  反射率キューブ [900–1640nm]
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

    # TODO: バッチ実行・外部呼び出しを想定する場合は logging.info() への移行が望ましい
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
    # SNR は 1.5μm 反射率を基準とした波長一様なガウスノイズ (InGaAs ノイズの簡略化)
    noise_sigma = spec.reflectance_1500 / config.snr if config.snr > 0 else 0.0

    cube         = np.zeros((ny, nx, nlam), dtype=np.float32)
    bd_map       = np.zeros((ny, nx),       dtype=np.float32)
    bd_map_ideal = np.zeros((ny, nx),       dtype=np.float32)

    for iy in range(ny):
        for ix in range(nx):
            ice_wt = ice_map[iy, ix]

            if config.add_dry_variation:
                perturb = rng.normal(0.0, DRY_SPEC_PERTURB_SCALE, nlam) * dry_std
                dry_px  = np.clip(dry_mean + perturb, REFL_MIN, None)
            else:
                dry_px = dry_mean.copy()

            target_bd         = gradient * ice_wt / 100.0
            absorption_factor = np.clip(1.0 - target_bd * unit_abs, 0.0, None)
            spec_px           = dry_px * absorption_factor

            if noise_sigma > 0:
                spec_px += rng.normal(0.0, noise_sigma, nlam)
            spec_px = np.clip(spec_px, REFL_MIN, REFL_MAX).astype(np.float32)

            cube[iy, ix, :]      = spec_px
            # BD は乾燥スペクトル正規化後の擬似反射率から計算し unit_abs スケールと整合させる
            bd_map_ideal[iy, ix] = compute_band_depth(absorption_factor, ALIS_WL)
            bd_map[iy, ix]       = compute_band_depth(spec_px / (dry_px + 1e-12), ALIS_WL)

    # target_bd と compute_band_depth の整合性チェック
    target_bd_map = (gradient * ice_map / 100.0).astype(np.float64)
    bd_map64      = bd_map_ideal.astype(np.float64)
    mask = target_bd_map > 1e-9
    if mask.any():
        rel_err = (bd_map64[mask] - target_bd_map[mask]) / target_bd_map[mask]
        rms_rel = float(np.sqrt(np.mean(rel_err ** 2)))
        verdict = "OK" if rms_rel <= 0.05 else "NG"
        # TODO: バッチ実行・外部呼び出しを想定する場合は logging.info() への移行が望ましい
        print(f"  BD consistency (RMS relative error): {rms_rel * 100.0:.2f}% ({verdict}, threshold=5%)")
    else:
        print("  BD consistency: skipped (all target_bd are ~0)")

    metadata = {
        "instrument":              "ALIS (Advanced Lunar Imaging Spectrometer)",
        "data_source":             "Araki & Saiki (2025), DOI: 10.60574/87068",
        "mineral_type":            config.mineral_type,
        "grain_size":              config.grain_size,
        "grain_size_um":           spec.grain_size_um,
        "reflectance_1500":        spec.reflectance_1500,
        "gradient_experimental":   gradient,
        "gradient_predicted_eq1":  gradient_pred,
        "ice_content_range_wt":    (config.ice_content_min, config.ice_content_max),
        "snr":                     config.snr,
        "n_bands":                 nlam,
        "n_dry_files":             group_data["n_dry"],
        "n_pairs":                 group_data["n_pairs"],
        "wavelength_range_nm":     (float(ALIS_WL[0]), float(ALIS_WL[-1])),
        "reference":               "Araki & Saiki (2025) Geochem. J., 59, 174-191",
    }
    return cube, ALIS_WL, ice_map, bd_map, metadata, group_data
