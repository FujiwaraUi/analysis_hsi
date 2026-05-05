# ALIS Mock Hyperspectral Data Generator

LUPEX 搭載近赤外撮像分光計 ALIS の仕様と、Araki & Saiki (2025) の実測スペクトルデータに基づく、月極域永久影領域（PSR）の水氷混合レゴリスを模擬したハイパースペクトルキューブ生成器。

## 目的

空間ビニングツール（修士論文テーマ）の入力テストデータとして、ALIS が月面 PSR で取得するであろうハイパースペクトルキューブを模擬する。水氷含有量と 1.5 μm バンド深さ（BD）の関係を論文の実測データに基づいて生成し、ビニング前後の BD 統計量（BD > 3σ 停止基準）の評価に用いることを想定している。

## 依拠する文献

| 略称 | 書誌情報 | 本コードでの使用箇所 |
|------|---------|-------------------|
| Araki & Saiki (2025) | Araki, R. and Saiki, K. (2025) Estimation of calibration lines for water ice content using 1.5-μm absorption band depth for lunar polar exploration. *Geochem. J.*, **59**, 174–191. DOI: 10.2343/geochemj.GJ25010 | 検量線モデル Eq. (1)、Table 2 の鉱物パラメータ、BD 定義、実測スペクトル CSV |
| 佐伯ほか (2025) | 佐伯和人ほか, 月極域探査 LUPEX 搭載 ALIS の開発状況, 宇宙科学技術連合講演会 3G08, 2025. | ALIS 仕様（波長範囲・分解能・検出器・iFOV） |
| Clark & Roush (1984) | Clark, R. N. and Roush, T. L. (1984) *J. Geophys. Res.*, **89**, 6329–6340. | BD の定義式 BD = (Rc − Rb) / Rc |

---

## 背景: Araki & Saiki (2025) の研究概要

### 研究目的

月の PSR に捕捉されている水氷の含有量を、近赤外分光観測からその場で推定する手法を確立すること。本研究は、月 PSR で想定される低水氷含有量（0.3–2.2 wt.%）の条件で検量線を初めて体系的に構築した。

### 実験概要

月レゴリスのアナログとして 4 種の鉱物粉末（オリビン Fo90、斜長石 An60、単斜輝石ディオプサイド、3 種等質量混合物）を 75–125 μm（細粒）と 180–250 μm（粗粒）の 2 粒径範囲に篩分けた計 8 試料を準備した。各試料表面に直径約 10 μm の水氷粒子を付着させ、900–1640 nm の近赤外反射スペクトルを測定した。

### バンド深さ（BD）の定義

水氷の 1.5 μm 吸収バンドの強度は Clark & Roush (1984) に基づく BD で定量される。

```
BD = (Rc − Rb) / Rc
```

Rb は吸収バンド底（約 1500 nm）の実測反射率、Rc は同波長における**コンティニュアム**（吸収バンドの両肩を結ぶ直線）上の値。BD は 0（吸収なし）から 1（完全吸収）の範囲をとる無次元量。

```
反射率
  ↑
  │  ●左肩(1350nm)─────────────●右肩(1640nm)  ← コンティニュアム
  │        ╲                  ╱
  │         ╲    Rc ←────────┤
  │          ╲      │        ╱
  │           ╲   BD│(Rc-Rb)╱
  │            ╲    │      ╱
  │             ●──Rb(バンド底)           ← 1500nm の実測値
  └──────────────────────────→ 波長 (nm)
```

本コードの `compute_band_depth` では、肩を 1350 nm と 1640 nm に固定した 2 点線形コンティニュアムを使用する。論文ではスムージングスプライン処理後に反復接線法で肩を自動決定しており、手法の差が BD の精度に影響しうる。

### 検量線（Calibration Line）

各鉱物種×粒径の組み合わせについて、水氷含有量 *f*（質量分率 = wt.%/100）と BD の間に原点を通る線形関係が得られる。

```
BD = gradient × f
```

**gradient**（検量線勾配）は鉱物の粒径と 1.5 μm 乾燥反射率で予測できる（Eq. (1)）。

```
gradient = 0.0517 × grain_size [μm] + 26.0 × reflectance_1500 − 17.4
```

---

## ALIS 仕様の反映

本コードは以下の ALIS 仕様を反映している（佐伯ほか, 2025, 表1）。

| パラメータ | 値 | コードでの実装 |
|-----------|-----|-------------|
| 観測波長 | 750–1650 nm | `ALIS_WL`: `np.arange(900, 1641, 5)` → 149 バンド（実験データカバー範囲）|
| 波長分解能 | < 5 nm | 5 nm 刻みで離散化 |
| 検出器 | SONY IMX990 | 検出器ノイズは SNR パラメータで簡易近似 |
| iFOV | 2 m × 1 cm @ 5 m | 空間グリッドの物理スケールとして意識、ピクセル数は可変 |

ALIS の機器仕様は 750–1650 nm だが、Araki & Saiki (2025) の実験データがカバーするのは 900–1640 nm のため、本コードでは実測根拠のある 900–1640 nm（149 バンド）のみを出力する。

---

## 実測スペクトルデータの種類と役割

Araki & Saiki (2025) の公開データ（DOI: 10.60574/87068）は、各実験回について 3 種類の CSV ファイルで構成される。`CF211016`（CPX 細粒、2021年10月16日の実験）を例にとる。

| ファイル | 格納場所 | 内容 |
|---------|---------|------|
| `CF211016_dry.csv` | `raw data/` | 乾燥試料（水氷なし）の実測反射率スペクトル |
| `CF211016_ice.csv` | `raw data/` | 着氷試料の実測反射率スペクトル（生データ、本コードでは不使用） |
| `CF211016.csv` | `fitted data/` | 着氷スペクトルにスムージングスプライン（λ = 0.4）を適用した平滑化スペクトル |

本コードでは `_dry.csv`（乾燥スペクトル）と `fitted data/*.csv`（平滑化済み着氷スペクトル）のペアを使用する。生の着氷スペクトルではなく平滑化済みデータを使う理由は、吸収形状の抽出がノイズに影響されにくいためである。

ファイル名の命名規則：
```
[鉱物略号][日付]_[番号]
  CF211016   → CPX (C), Fine (F), 2021年10月16日
  OC210318_1 → Olivine (O), Coarse (C), 2021年3月18日, 1回目
```

---

## モックキューブ生成の処理フロー

### 第 1 段階: 実測データから材料を準備する（`build_group_data`）

指定された鉱物種×粒径の組み合わせについて、2 種類の材料を構築する。

**材料 1: 乾燥スペクトルの平均 `dry_mean(λ)`**

グループ内の全 `*_dry.csv` を 900–1640 nm グリッドに補間し、全ファイルで平均した代表乾燥スペクトル。標準偏差 `dry_std(λ)` もピクセル間ばらつきの付加に使用する。

**材料 2: 単位吸収プロファイル `unit_abs(λ)`**

実験データから抽出した**水氷吸収の波長形状テンプレート**。スペクトル合成式：

```
R_mock(λ) = dry(λ) × (1 − target_bd × unit_abs(λ))
```

において `unit_abs(λ)` が吸収の波長分布を担い、`target_bd` を変えることで任意の水氷量のスペクトルを生成できる。この式が target_bd に対して正しく線形に動くには、`unit_abs(λ)` が適切なスケールに正規化されている必要がある（後述）。

抽出手順：
1. 各実験ペア (dry_i, fitted_ice_i) について波長ごとの分率吸収を計算する：
   ```
   frac_abs_i(λ) = (dry_i(λ) − ice_i(λ)) / dry_i(λ)
   ```
   これは乾燥スペクトルを基準とした「各波長での反射率の低下割合」。
2. `EXPERIMENT_ICE_CONTENT` の水氷含有量から期待 BD を計算し正規化する：
   ```
   expected_bd_i = gradient × ice_wt / 100
   unit_abs_i(λ) = frac_abs_i(λ) / expected_bd_i
   ```
3. 全ペアで平均する → `unit_abs(λ)`（中間状態）
4. **スケール補正（再正規化）**：ここで得られた `unit_abs(λ)` をそのまま合成式に使うと、`compute_band_depth` で計算した BD が `target_bd` からずれる。原因は基準面の違いで、`frac_abs` は乾燥スペクトルを分母とするが、`compute_band_depth` はコンティニュアム直線を分母とする。そのため、小さな係数 α で `1 − α × unit_abs(λ)` を作って `compute_band_depth` を呼び出し、「このスケールのとき実際に出る BD」を逆算して割り戻す：
   ```python
   peak_bd = compute_band_depth(1 − α × unit_abs) / α
   unit_abs ← unit_abs / peak_bd
   ```
   この補正後、合成式の `target_bd` と `compute_band_depth` の返り値がほぼ一致するようになる。

> `unit_abs(λ)` という量は本コード独自の中間変数であり、Araki & Saiki (2025) には直接の記述はない。

### 第 2 段階: ピクセルごとにスペクトルを合成する（`create_alis_mock_araki`）

```
(a) 氷量マップ生成（create_ice_map）
    gradient: 中心が最大、外縁に向けて減少する同心円状分布
    patches:  ランダムな円形パッチ
    uniform:  一様

(b) 各ピクセルでスペクトルを合成:
    dry_px(λ) = dry_mean(λ) + rng × dry_std(λ)   ← ピクセル間ばらつき
    target_bd = gradient × ice_wt / 100
    R_mock(λ) = dry_px(λ) × (1 − target_bd × unit_abs(λ))
    R_mock(λ) += rng × noise_sigma                ← 検出器ノイズ

(c) BD を再計算して bd_map に格納
```

### データフロー図

```
実測 CSV 群
  ├─ *_dry.csv  →  ALIS グリッドに補間 → 平均化 → dry_mean(λ), dry_std(λ)
  └─ fitted/*.csv → dry との差分 → frac_abs(λ)
                         │
                         └─ ice_wt で正規化 → 平均 → 再正規化 → unit_abs(λ)
                                                        ↓
各ピクセル:
  dry_mean + ばらつき → dry_px(λ)
  ice_map(iy, ix)    → target_bd = gradient × ice_wt / 100
                                   ↓
  R_mock(λ) = dry_px(λ) × (1 − target_bd × unit_abs(λ)) + noise
                                   ↓
                             cube(iy, ix, λ)
```

---

## Table 2 の鉱物パラメータ

`GROUP_SPECS` に 4 鉱物種 × 2 粒径 = 8 グループの実験値を格納している。

| 鉱物 | 粒径 | R(1.5 μm) | gradient（実測） |
|------|------|-----------|----------------|
| Olivine (Fo90) | coarse (180–250 μm) / fine (75–125 μm) | 0.460 / 0.542 | 5.18 / 1.83 |
| Plagioclase (An60) | coarse / fine | 0.565 / 0.679 | 8.48 / 5.83 |
| Clinopyroxene (Diopside) | coarse / fine | 0.610 / 0.665 | 11.87 / 5.91 |
| Three-mineral mixture | coarse / fine | 0.509 / 0.634 | 5.71 / 3.11 |

gradient は単斜輝石 > 斜長石 > 混合物 > オリビンの順に大きく、1.5 μm 乾燥反射率の大小関係と一致する。スペクトル生成には実験 gradient 値を使用し、Eq. (1) 予測値は出力プロットに参照線として併記する。

---

## ファイル構成

```
src/
  claude_make_mock_araki.py   # コア: データ読み込み・スペクトル合成（matplotlib 非依存）
  plot_araki.py               # 可視化: plot_* 関数群
  main.py                     # エントリーポイント: 全グループ生成・プロット保存
  alis_mock_output_araki/     # 実行時に自動生成
    ├── dry_spectra_all_groups.png
    ├── unit_absorption_profiles.png
    ├── olivine_coarse_calibration.png
    ├── olivine_coarse_montage.png
    ├── olivine_coarse_spectra.png
    └── ...（4 鉱物 × 2 粒径 × 3 種類 = 24 ファイル）
Mock_Data/
  Experimental_Data/          # Araki & Saiki (2025) DOI: 10.60574/87068
    olivine/ plagioclase/ cpx/ mixture/
      180-250um/ 75-125um/
        raw data/    *_dry.csv, *_ice.csv
        fitted data/ *.csv
```

---

## 使用方法

### 依存パッケージ

```
numpy
matplotlib
```

### 基本実行

```bash
python src/main.py
```

全グループ分の calibration・montage・spectra プロットと全体比較プロット 2 枚が `src/alis_mock_output_araki/` に出力される。

### モジュールとしての使用

```python
from claude_make_mock_araki import ALISMockConfigAraki, create_alis_mock_araki

config = ALISMockConfigAraki(
    nx=100, ny=50,
    mineral_type="mixture",   # "olivine" / "plagioclase" / "clinopyroxene" / "mixture"
    grain_size="coarse",      # "coarse" (180–250 μm) / "fine" (75–125 μm)
    ice_content_min=0.0,
    ice_content_max=2.2,      # wt.%
    ice_pattern="gradient",   # "gradient" / "patches" / "uniform"
    snr=100,
    seed=42,
)

cube, wavelengths, ice_map, bd_map, metadata, group_data = create_alis_mock_araki(config)

# cube:        shape (50, 100, 149)  反射率キューブ [900–1640 nm]
# wavelengths: shape (149,)          波長配列 (nm)
# ice_map:     shape (50, 100)       入力水氷量 (wt.%)
# bd_map:      shape (50, 100)       復元 BD
# metadata:    dict                  生成パラメータと参照情報
# group_data:  dict                  乾燥スペクトル・単位吸収プロファイル（検証用）
```

---

## 出力の解釈

### 乾燥スペクトル比較（`dry_spectra_all_groups.png`）

粗粒/細粒の 2 パネルに、4 鉱物の個別スペクトル（薄い線）と平均（太い線）を重ね描き。`unit_abs(λ)` 抽出の元となる実測スペクトルの形状を確認するためのプロット。

### 単位吸収プロファイル（`unit_absorption_profiles.png`）

各グループで抽出した `unit_abs(λ)` の比較。1.5 μm 付近でピークを持ち、鉱物・粒径によって形状が異なる。

### 検量線検証プロット（`*_calibration.png`）

横軸が水氷含有量 (wt.%)、縦軸が復元 BD。散布点が実験 gradient 線（実線）に沿って分布していれば、スペクトル生成 → BD 復元のパイプラインが整合していることを示す。Eq. (1) 予測線（破線）は参照として併記。

### モンタージュ（`*_montage.png`）

950、1050、1200、1350、1500、1600 nm の反射率画像・入力氷量マップ・復元 BD マップの 8 パネル構成。1500 nm 付近で氷量の多い領域の反射率が低下する様子が確認できる。

### スペクトル比較（`*_spectra.png`）

氷量の 5/25/50/95 パーセンタイルに対応するピクセルのスペクトルと実測乾燥スペクトル平均（破線）を重ね描き。論文 Fig. 4 に対応する。

---

## 既知の制約

1. **粒径は離散値のみ**: 実験データが coarse（180–250 μm）と fine（75–125 μm）の 2 グループしかないため、中間粒径の指定はできない。
2. **コンティニュアム決定法の簡略化**: 本コードでは 2 点線形コンティニュアム（肩: 1350 nm, 1640 nm 固定）を使用する。論文ではスムージングスプライン＋反復接線法により肩の位置を自動決定しており、手法の差が BD 精度に影響しうる。
3. **低氷量域の非線形性**: 論文では 0.85 wt.% 以下で BD が検量線を下回る傾向が報告されている（Hapke モデルによる再現あり）。本コードはこの非線形性を実装していない。
4. **ノイズモデルの簡略化**: 1.5 μm 反射率基準の波長一様ガウスノイズとして適用。実際の InGaAs 検出器の波長依存ノイズ特性は未実装。
5. **Eq. (1) の適用範囲**: 粒径 75–250 μm、乾燥反射率 0.4–0.7 の範囲で検証されたモデルであり、範囲外への外挿は信頼性が低い。

---

## ライセンスと引用

本コードは研究・教育目的で作成された。Araki & Saiki (2025) は CC BY 4.0 ライセンスで公開されている。本コードを使用する場合は、上記文献への適切な引用を行うこと。

---

## 用語集

| 略語 | 意味 |
|------|------|
| ALIS | Advanced Lunar Imaging Spectrometer（LUPEX 搭載近赤外撮像分光計） |
| BD | Band Depth（バンド深さ）: BD = (Rc − Rb) / Rc |
| gradient | Calibration-line gradient（検量線勾配）: BD = gradient × f における比例定数 |
| LUPEX | Lunar Polar Exploration（JAXA-ISRO 月極域探査ミッション） |
| PSR | Permanently Shadowed Region（永久影領域） |
| wt.% | weight percent（質量パーセント）。質量分率 f = wt.%/100 |
