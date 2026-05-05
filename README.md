# ALIS Mock Hyperspectral Data Generator

LUPEX 搭載近赤外画像分光装置 ALIS の仕様と、Araki & Saiki (2025) の実測スペクトルデータに基づく、月極域 PSR 水氷混合レゴリスのモックハイパースペクトルキューブ生成器。

## 目的

空間ビニングツール（修士論文テーマ α）の入力テストデータとして、ALIS が取得するであろうハイパースペクトルキューブを模擬する。水氷含有量と 1.5 μm バンド深さ（BD）の関係を論文の実測データに基づいて生成し、ビニング前後の BD 統計量（BD > 3σ 停止基準）の評価に用いることを想定している。

## 依拠する文献

| 略称 | 書誌情報 | 本コードでの使用箇所 |
|------|---------|-------------------|
| Araki & Saiki (2025) | Araki, R. and Saiki, K. (2025) Estimation of calibration lines for water ice content using 1.5-μm absorption band depth for lunar polar exploration. *Geochem. J.*, **59**, 174–191. DOI: 10.2343/geochemj.GJ25010 | 校正線モデル Eq. (1)、Table 2 の鉱物パラメータ、BD 定義、実測スペクトル CSV |
| 佐伯ほか (2025) | 佐伯和人ほか, 月極域探査 LUPEX 搭載 ALIS の開発状況, 宇宙科学技術連合講演会 3G08, 2025. | ALIS 仕様（波長範囲・分解能・検出器・iFOV） |
| Clark & Roush (1984) | Clark, R. N. and Roush, T. L. (1984) *J. Geophys. Res.*, **89**, 6329–6340. | BD の定義式 BD = (Rc − Rb) / Rc |

## ALIS 仕様の反映

本コードは以下の ALIS 仕様を反映している（佐伯ほか, 2025, 表1）。

| パラメータ | 値 | コードでの実装 |
|-----------|-----|-------------|
| 観測波長 | 750–1650 nm | `ALIS_WL`: `np.arange(750, 1651, 5)` → 181 バンド |
| 波長分解能 | < 5 nm | 5 nm 刻みで離散化 |
| 検出器 | SONY IMX990 | 検出器ノイズモデルは SNR パラメータで簡易近似 |
| iFOV | 2 m × 1 cm @ 5 m | 空間グリッドの物理スケールとして意識、ピクセル数は可変 |

## Araki & Saiki (2025) 実測データの使用

### 実測スペクトルデータ（DOI: 10.60574/87068）

134 本の CSV ファイルを `Mock_Data/Experimental_Data/` に格納し、以下の用途で使用している。

| ファイル種別 | 形式 | 使用方法 |
|------------|------|---------|
| `*_dry.csv` | ヘッダなし、800–1699 nm | 乾燥スペクトルとして読み込み、グループ内で平均化 |
| `fitted/*.csv` | ヘッダあり、スムージングスプライン済み | 単位 BD 吸収プロファイルの抽出に使用 |

### 単位吸収プロファイルの抽出

各 (dry, fitted_ice) ペアから以下の手順で「単位 BD あたりの吸収形状」を求め、グループ内で平均する。

```
frac_abs(λ) = (R_dry(λ) − R_ice(λ)) / R_dry(λ)
unit_abs(λ)  = frac_abs(λ) / (gradient_exp × ice_content / 100)
```

モックスペクトルの生成式:

```
R_mock(λ) = R_dry_mean(λ) × (1 − target_BD × unit_abs(λ))
```

### Eq. (1): 校正線勾配の予測（参照用）

```
gradient = 0.0517 × grain_size_μm + 26.0 × reflectance_1500nm − 17.4
```

スペクトル生成には Table 2 の実験勾配値を使用する。Eq. (1) 予測値は出力プロットに参照線として併記される。

### Table 2 の鉱物パラメータ

`GROUP_SPECS` に 4 鉱物種 × 2 粒径 = 8 グループの実験値を格納している。

| 鉱物 | 粒径 | R(1.5 μm) | 勾配（実測） |
|------|------|-----------|------------|
| Olivine (Fo90) | coarse (180–250 μm) / fine (75–125 μm) | 0.460 / 0.542 | 5.18 / 1.83 |
| Plagioclase (An60) | coarse / fine | 0.565 / 0.679 | 8.48 / 5.83 |
| Clinopyroxene (Diopside) | coarse / fine | 0.610 / 0.665 | 11.87 / 5.91 |
| Three-mineral mixture | coarse / fine | 0.509 / 0.634 | 5.71 / 3.11 |

## ファイル構成

```
src/
  claude_make_mock_araki.py         # 本体（実測データ版）
  alis_mock_output_araki/           # 実行時に自動生成
    ├── dry_spectra_all_groups.png  # 全グループの乾燥スペクトル比較
    ├── unit_absorption_profiles.png # 単位 BD 吸収プロファイル比較
    ├── olivine_coarse_calibration.png
    ├── olivine_coarse_montage.png
    ├── olivine_coarse_spectra.png
    ├── olivine_fine_calibration.png
    ├── olivine_fine_montage.png
    ├── olivine_fine_spectra.png
    ├── plagioclase_coarse_calibration.png
    ├── ...（4 鉱物 × 2 粒径 × 3 種類 = 24 ファイル）
    └── mixture_fine_spectra.png
Mock_Data/
  Experimental_Data/                # Araki & Saiki (2025) DOI: 10.60574/87068
    olivine/
      180-250um/
        raw data/   *_dry.csv, *_ice.csv
        fitted data/ *.csv
      75-125um/
        ...
    plagioclase/ cpx/ mixture/
      ...
```

## コード構造（クラス・関数の役割）

- `src/claude_make_mock_araki.py`
  - `GroupSpec`
    - 鉱物種×粒径グループの実験パラメータ（フォルダ名・粒径・R(1.5μm)・実験勾配）を保持するデータコンテナ。
  - `GROUP_SPECS`
    - 8 グループ分の `GroupSpec` を `(mineral_type, grain_size)` キーで参照できる辞書。
  - `EXPERIMENT_ICE_CONTENT`
    - 46 実験ファイルの水氷含有量 (wt.%) を格納した辞書。
  - `FITTED_TO_RAW_STEM`
    - fitted ファイル名と raw ファイル名が不一致な 7 ケースのマッピング。`None` は raw データが存在しないことを示す。
  - `predict_gradient(grain_size_um, reflectance_1500)`
    - Araki & Saiki (2025) の Eq. (1) に基づき、校正線勾配の予測値を返す（参照用）。
  - `_load_csv(path, has_header)`
    - CSV から波長・反射率を読み込み、負の反射率を 0 にクリップして返す。
  - `interp_to_alis(wl, refl)`
    - 任意の波長グリッドを ALIS グリッド（750–1650 nm, 5 nm）に変換する。750–800 nm は線形外挿。
  - `build_group_data(spec, data_dir)`
    - グループ内の全 `_dry.csv` を読み込んで平均乾燥スペクトルを構築し、fitted ice CSV から単位吸収プロファイルを抽出する。
  - `ALISMockConfigAraki`
    - 空間サイズ・鉱物種・粒径カテゴリ・氷分布パターン・SNR・seed などを保持する設定オブジェクト。
  - `create_ice_map(config, rng)`
    - `ice_pattern` に応じて、ユニフォーム、勾配、またはパッチ状の水氷含有量マップを生成する。
  - `compute_band_depth(spectrum, wavelengths, ...)`
    - 2 点線形コンティニュアムから BD = (Rc − Rb) / Rc を計算する。
  - `create_alis_mock_araki(config, data_dir)`
    - 実測乾燥スペクトル・単位吸収プロファイル・氷分布・ノイズを組み合わせてハイパースペクトルキューブを生成し、キューブ・波長・氷量マップ・BD マップ・メタデータ・グループデータを返す。
  - `plot_calibration(ice_map, bd_map, metadata, out_path)`
    - BD vs 水氷含有量の散布図に実験勾配線と Eq. (1) 予測線を重ねて描画する（論文 Fig. 8 対応）。
  - `plot_montage(cube, ice_map, bd_map, metadata, out_path)`
    - 代表波長 6 枚・氷量マップ・BD マップの 2×4 モンタージュを保存する。
  - `plot_spectra_comparison(cube, ice_map, group_data, metadata, out_path)`
    - 氷量分位点のピクセルスペクトルと実測乾燥スペクトル平均を重ね描きする（論文 Fig. 4 対応）。
  - `plot_dry_spectra_all_groups(data_dir, out_path)`
    - 全 8 グループの個別・平均乾燥スペクトルを粗粒/細粒の 2 パネルで比較する。
  - `plot_unit_absorption_profiles(data_dir, out_path)`
    - 全 8 グループの単位 BD 吸収プロファイルを粗粒/細粒の 2 パネルで比較する。
  - `__main__` ブロック
    - 8 グループ全てのモック生成・可視化を実行し、26 ファイルを `src/alis_mock_output_araki/` に出力する。

## 使用方法

### 依存パッケージ

```
numpy
matplotlib
```

### 基本実行

```bash
python src/claude_make_mock_araki.py
```

8 グループ分の calibration・montage・spectra プロットと全体比較プロット 2 枚が `src/alis_mock_output_araki/` に出力される。

### モジュールとしての使用

```python
from src.claude_make_mock_araki import ALISMockConfigAraki, create_alis_mock_araki

config = ALISMockConfigAraki(
    nx=100, ny=50,
    mineral_type="mixture",   # "olivine" / "plagioclase" / "clinopyroxene" / "mixture"
    grain_size="coarse",      # "coarse" (180-250 μm) / "fine" (75-125 μm)
    ice_content_min=0.0,
    ice_content_max=2.2,      # wt.%
    ice_pattern="gradient",   # "gradient" / "patches" / "uniform"
    snr=100,
    seed=42,
)

cube, wavelengths, ice_map, bd_map, metadata, group_data = create_alis_mock_araki(config)

# cube:       shape (50, 100, 181), 反射率キューブ
# wavelengths: shape (181,), 波長配列 (nm)
# ice_map:    shape (50, 100), 入力水氷量 (wt.%)
# bd_map:     shape (50, 100), 復元 BD
# metadata:   dict, 生成パラメータと参照情報
# group_data: dict, 乾燥スペクトル・単位吸収プロファイル（検証用）
```

## 出力の解釈

### 乾燥スペクトル比較（`dry_spectra_all_groups.png`）

粗粒/細粒の 2 パネルに、4 鉱物の個別スペクトル（薄い線）と平均（太い線）を重ね描きしたもの。実データ由来の乾燥スペクトルの形状が確認できる。

### 単位吸収プロファイル（`unit_absorption_profiles.png`）

各グループで抽出した「単位 BD あたりの吸収形状」。1.5 μm 付近でピークを持ち、鉱物・粒径によって形状が異なることが確認できる。

### 校正線検証プロット（`*_calibration.png`）

横軸が水氷含有量 (wt.%)、縦軸が復元 BD。実験勾配による校正線（実線）と Eq. (1) 予測線（破線）を併記。散布点が実験勾配線に沿って分布していれば、スペクトル生成→BD 復元のパイプラインが整合していることを示す。

### モンタージュ（`*_montage.png`）

6 波長（800, 1000, 1200, 1400, 1500, 1600 nm）の反射率画像・入力氷量マップ・復元 BD マップの 8 パネル構成。1500 nm 付近で氷量の多い領域の反射率が低下する様子が確認できる。

### スペクトル比較（`*_spectra.png`）

氷量の 5/25/50/95 パーセンタイルに対応するピクセルのスペクトルと実測乾燥スペクトル平均（破線）を重ね描きしたもの。論文 Fig. 4 に対応する。

## 既知の制約

1. **粒径は離散値のみ**: 実験データが `"coarse"` (180–250 μm) と `"fine"` (75–125 μm) の 2 グループしかないため、中間粒径の指定はできない。
2. **コンティニュアム決定法**: 2 点線形コンティニュアム（肩波長 1350 nm, 1640 nm）を使用。論文ではスムージングスプライン＋反復接線法を用いており、手法の差が BD の復元精度に影響しうる。
3. **750–800 nm の外挿**: 実測データの下限は 800 nm であり、ALIS の 750–800 nm 域は線形外挿で補完している。
4. **低氷量域の非線形性**: 論文では 0.85 wt.% 以下で BD が校正線を下回る傾向が報告されている（Hapke モデルによる再現あり）。本コードはこの非線形性を実装していない。
5. **ALIS-L（能動照明）未対応**: 影領域観測時の光源特性は未実装。
6. **Eq. (1) の適用範囲**: 粒径 75–250 μm、乾燥反射率 0.4–0.7 の範囲で検証されたモデルであり、範囲外への外挿は信頼性が低い。

## ライセンスと引用

本コードは研究・教育目的で作成された。Araki & Saiki (2025) は CC BY 4.0 ライセンスで公開されている。本コードを使用する場合は、上記文献への適切な引用を行うこと。

## 用語集

| 略語 | 意味 |
|------|------|
| BD | Band Depth（バンド深さ）: BD = (Rc − Rb) / Rc |
| PSR | Permanently Shadowed Region（永久影領域） |
| ALIS | Advanced Lunar Imaging Spectrometer（LUPEX 搭載分光器） |
| LUPEX | Lunar Polar Exploration（JAXA-ISRO 月極域探査ミッション） |
| wt.% | weight percent（質量パーセント） |
