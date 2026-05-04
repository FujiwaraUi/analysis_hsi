# ALIS Mock Hyperspectral Data Generator

LUPEX 搭載近赤外画像分光装置 ALIS の仕様と、Araki & Saiki (2025) の校正線モデルに基づく、月極域 PSR 水氷混合レゴリスのモックハイパースペクトルキューブ生成器。

## 目的

空間ビニングツール（修士論文テーマ α）の入力テストデータとして、ALIS が取得するであろうハイパースペクトルキューブを模擬する。水氷含有量と 1.5 μm バンド深さ（BD）の関係を論文の校正線モデルに基づいて生成し、ビニング前後の BD 統計量（BD > 3σ 停止基準）の評価に用いることを想定している。

## 依拠する文献

| 略称 | 書誌情報 | 本コードでの使用箇所 |
|------|---------|-------------------|
| Araki & Saiki (2025) | Araki, R. and Saiki, K. (2025) Estimation of calibration lines for water ice content using 1.5-μm absorption band depth for lunar polar exploration. *Geochem. J.*, **59**, 174–191. DOI: 10.2343/geochemj.GJ25010 | 校正線モデル Eq. (1)、Table 2 の鉱物パラメータ、BD 定義 |
| 佐伯ほか (2025) | 佐伯和人ほか, 月極域探査 LUPEX 搭載 ALIS の開発状況, 宇宙科学技術連合講演会 3G08, 2025. | ALIS 仕様（波長範囲・分解能・検出器・iFOV） |
| Clark & Roush (1984) | Clark, R. N. and Roush, T. L. (1984) *J. Geophys. Res.*, **89**, 6329–6340. | BD の定義式 BD = (Rc − Rb) / Rc |

## ALIS 仕様の反映

本コードは以下の ALIS 仕様を反映している（佐伯ほか, 2025, 表1）。

| パラメータ | 値 | コードでの実装 |
|-----------|-----|-------------|
| 観測波長 | 750–1650 nm | `alis_wavelengths()`: `np.arange(750, 1650.1, 5)` → 181 バンド |
| 波長分解能 | < 5 nm | 5 nm 刻みで離散化 |
| 検出器 | SONY IMX990 | （検出器ノイズモデルは SNR パラメータで簡易近似） |
| iFOV | 2 m × 1 cm @ 5 m | （空間グリッドの物理スケールとして意識、ピクセル数は可変） |

## Araki & Saiki (2025) 校正線モデルの実装

### Eq. (1): 校正線勾配の予測

```
gradient = 0.0517 × grain_size_μm + 26.0 × reflectance_1500nm − 17.4
```

`predict_gradient(grain_size_um, reflectance_1500)` として実装。入力する粒径と乾燥反射率から校正線勾配を予測し、BD = gradient × ice_content / 100 の関係でバンド深さを決定する。

### Table 2 の鉱物パラメータ

`MINERALS` 辞書に以下の 4 鉱物種 × 2 粒径の実験値を格納している。

| 鉱物 | 粒径 (μm) | R(1.5 μm) | 勾配（実測） |
|------|-----------|-----------|------------|
| Olivine (Fo90) | 215 / 100 | 0.460 / 0.542 | 5.18 / 1.83 |
| Plagioclase (An60) | 215 / 100 | 0.565 / 0.679 | 8.48 / 5.83 |
| Clinopyroxene (Diopside) | 215 / 100 | 0.610 / 0.665 | 11.87 / 5.91 |
| Three-mineral mixture | 215 / 100 | 0.509 / 0.634 | 5.71 / 3.11 |

中間粒径を指定した場合は、粗粒（215 μm）と細粒（100 μm）の値を線形内挿する。

## ファイル構成

```
alis_mock.py                  # 本体（単一ファイル、外部依存は numpy + matplotlib のみ）
alis_mock_output/             # 実行時に自動生成
  ├── olivine_215um_calibration.png      # BD vs ice content 検証プロット
  ├── olivine_215um_montage.png          # 代表バンド画像 + 氷量マップ + BD マップ
  ├── plagioclase_215um_calibration.png
  ├── plagioclase_215um_montage.png
  ├── clinopyroxene_215um_calibration.png
  ├── clinopyroxene_215um_montage.png
  ├── mixture_150um_calibration.png
  ├── mixture_150um_montage.png
  └── mixture_spectra_comparison.png     # 氷量別スペクトル比較（Fig. 4 対応）
```

## 使用方法

### 依存パッケージ

```
numpy
matplotlib
```

### 基本実行

```bash
python alis_mock.py
```

4 鉱物種の校正線検証プロットとモンタージュが `alis_mock_output/` に出力される。

### モジュールとしての使用

```python
from alis_mock import ALISMockConfig, create_alis_mock, compute_band_depth

config = ALISMockConfig(
    nx=100, ny=50,
    mineral_type="mixture",    # "olivine" / "plagioclase" / "clinopyroxene" / "mixture"
    grain_size_um=150,         # 75–250 の範囲
    ice_content_min=0.0,
    ice_content_max=2.2,       # wt.%
    ice_pattern="gradient",    # "gradient" / "patches" / "uniform"
    snr=100,
    seed=42,
)

cube, wavelengths, ice_map, bd_map, metadata = create_alis_mock(config)

# cube:        shape (50, 100, 181), 反射率キューブ
# wavelengths: shape (181,), 波長配列 (nm)
# ice_map:     shape (50, 100), 入力水氷量 (wt.%)
# bd_map:      shape (50, 100), 復元 BD
# metadata:    dict, 生成パラメータと参照情報
```

## 出力の解釈

### 校正線検証プロット（`*_calibration.png`）

横軸が水氷含有量 (wt.%)、縦軸が復元 BD。赤線が Eq. (1) による理論校正線。散布点が理論線周辺に分布していれば、スペクトル生成→BD 復元のパイプラインが校正線モデルと整合していることを示す。散布の幅は設定した SNR に対応する。

### モンタージュ（`*_montage.png`）

6 波長の反射率画像（800, 1000, 1200, 1400, 1500, 1600 nm）、入力氷量マップ、復元 BD マップの 8 パネル構成。1500 nm 付近で氷量の多い領域の反射率が低下する（吸収帯が深い）ことが視覚的に確認できる。

### スペクトル比較（`mixture_spectra_comparison.png`）

氷量の異なるピクセルのスペクトルを 1 つのプロットに重ねたもの。論文 Fig. 4 に対応する。乾燥スペクトル（破線）に対して、氷量が増えるほど 1500 nm 付近で反射率が低下する傾向が確認できる。

## 既知の制約

1. **乾燥スペクトルは簡易モデル**: 鉱物種ごとのコンティニュアム形状と 1 μm 吸収帯をガウス関数で近似している。Araki & Saiki の公開データ（DOI: 10.60574/87068, 134 CSV ファイル）を読み込んで置き換えれば精度が向上する。
2. **吸収帯形状**: 1.5 μm 水氷吸収をガウス型（σ = 60 nm）で近似している。論文の実測スペクトルでは非対称性がありうる。
3. **コンティニュアム決定法**: 本コードは 2 点線形コンティニュアム（肩波長 1350 nm, 1650 nm）を使用。論文ではスムージングスプライン＋反復接線法を用いており、手法の差が BD の復元精度に影響しうる。
4. **ALIS-L（能動照明）未対応**: 影領域観測時の光源特性は未実装。能動照明下の SNR モデルは ALIS 側の詳細仕様の公開待ちである。
5. **低氷量域の非線形性**: 論文では 0.85 wt.% 以下で BD が校正線を下回る傾向が報告されている（Hapke モデルによる再現あり）。本コードはこの非線形性を実装していない。
6. **Eq. (1) の適用範囲**: 粒径 75–250 μm、乾燥反射率 0.4–0.7 の範囲で検証されたモデルであり、範囲外への外挿は信頼性が低い。

## ライセンスと引用

本コードは研究・教育目的で作成された。Araki & Saiki (2025) は CC BY 4.0 ライセンスで公開されている。本コードを使用する場合は、上記文献への適切な引用を行うこと。
