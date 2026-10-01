# 実験担当者向け計算コード共有ガイド

## 共有すべき最小ファイル

| 優先度 | ファイル | 役割 |
|---|---|---|
| 必須 | `code/data/ho_fig1e_absorption.csv` | Ho et al. 公開吸収曲線の再構成データ |
| 必須 | `code/ho_spectrum_model.py` | エネルギー・圧力補間と静水圧吸収カーネル |
| 必須 | `code/ho_odmr_sensitivity.py` | 固定パワー吸収率、感度、約 441 nm 基準 |
| 必須 | `code/nv_model_111.py` | [111] 軸対称応力テンソルと静水圧光学基準の分離 |
| 必須 | `code/report_111_anisotropic.py` | 圧力・応力・有限パワー最適集合のCLIレポート |
| 必須 | `code/fig_alpha_invariance.py` | \(\alpha=0.5\)–0.7 の不変性図の再生成 |
| 必須 | `code/nv_model.py` | 偏差応力補正の差分を与える現象論モデル |
| 必須 | `code/calc_alpha_corrected_optimum.py` | \(\alpha=0.5\)–0.95 の条件付き最適波長計算 |
| 必須 | `code/fig_alpha_optimum_tolerance.py` | \(\alpha\)別最適値、5%許容帯、固定レーザー判定図 |
| 補助 | `code/theory_a1_generalization.py` | 5%許容帯と有限パワーのレベル集合 |
| 補助 | `code/fig_style.py` | 図の共通スタイルとPDF幅検証 |
| 検証 | `code/tests/test_nv_model_111.py` | 応力テンソルと静水圧基準の回帰試験 |
| 検証 | `code/tests/test_fig_alpha_invariance.py` | \(\alpha\)走査と図出力の回帰試験 |
| 検証 | `code/tests/test_calc_alpha_corrected_optimum.py` | \(\alpha=0.6\) の 488.10 nm と静水圧アンカーの回帰試験 |
| 検証 | `code/tests/test_fig_alpha_optimum_tolerance.py` | 許容帯、最悪感度、PNG/PDF出力の回帰試験 |
| 説明 | `docs/experiment/theory_limits_and_required_measurements.md` | 理論限界、必要な実測値、モデル更新条件 |

`nv_model.py` の絶対最適値は採用しない。Ho の 440.65 nm を絶対基準とし、
`calc_alpha_corrected_optimum.py` が同モデルの \(\alpha\) 依存差分だけを用いる。
`nv_model_power.py` と `fig1_*`–`fig3_*` は本計算には不要である。

## 実行方法

`code/` から実行する。

```powershell
python report_111_anisotropic.py --pressure 120 --alpha 0.60 --axial-stress 120
python fig_alpha_invariance.py
python calc_alpha_corrected_optimum.py
python fig_alpha_optimum_tolerance.py
python -m pytest tests/test_nv_model_111.py tests/test_fig_alpha_invariance.py tests/test_calc_alpha_corrected_optimum.py tests/test_fig_alpha_optimum_tolerance.py -q
```

生成図:

- `docs/experiment/figures/alpha_invariance_120gpa.png`
- `docs/experiment/figures/alpha_invariance_120gpa.pdf`
- `docs/experiment/figures/alpha_optimum_tolerance_120gpa.png`
- `docs/experiment/figures/alpha_optimum_tolerance_120gpa.pdf`

## 二つの計算結果

静水圧カーネル固定モデルは全ての \(\alpha\) で約 441 nm を返す。一方、
Ho の静水圧値を \(\alpha=1\) の基準とし、Hilberer の偏差応力係数
\(k_q/k_h=-0.734\) で差分補正すると、120 GPa で次を得る。

| \(\alpha\) | 条件付き最適波長 |
|---:|---:|
| 0.50 | 501.50 nm |
| 0.56 | 493.38 nm |
| 0.60 | 488.10 nm |
| 0.70 | 475.36 nm |
| 0.80 | 463.23 nm |
| 0.90 | 451.67 nm |
| 0.95 | 446.10 nm |
| 1.00 | 440.65 nm |

前者は異方応力結合をゼロと置く帰無仮説、後者は ZPL シフト駆動分だけを含む
対抗仮説である。どちらも異方応力下の最適波長を直接検証した結果ではない。

## 実験担当者へ渡す説明文

\(\sigma_{zz}=120\ \mathrm{GPa}\)、
\(\sigma_{xx}=\sigma_{yy}=\alpha\sigma_{zz}\) として
Ho et al. の静水圧吸収カーネルから、\(\alpha=1\) の絶対基準は約 441 nm となる。
Ho の \(\alpha\approx0.95\) はマイクロピラー検証条件であり、吸収カーネル自体の
\(\alpha\) ではない。Hilberer et al. の通常型 DAC（\(\alpha\approx0.56\)）と
マイクロピラー（\(\alpha\approx0.95\)）の ZPL 応力応答を使った条件付き一般化では、
\(\alpha=0.5\)–0.95 の最適波長を約 502–446 nm、\(\alpha=0.6\) を 488.10 nm と推定する。

文献は ZPL 応力応答との整合性を与えるが、\(\lambda_{\mathrm{opt}}(\alpha)\) を
直接測定していない。特に Huang–Rhys 因子と sideband 形状の異方応力依存性は
未校正であるため、この曲線は確定値ではなく実験計画用の対抗仮説である。

実験では、複数の \(\alpha\) で 441、457、473、488、505、532 nm を低入射パワー・
同一検出条件で比較し、photon rate \(R\)、contrast \(C\)、linewidth
\(\Delta\nu\)、NV⁻/NV⁰ 比、発光スペクトル、局所 \(\alpha\) を同時に記録する。
可能なら 440–505 nm を連続走査し、最適波長の絶対値と
\(d\lambda_{\mathrm{opt}}/d\alpha\) を求める。実測が得られるまで理論探索は停止する。
