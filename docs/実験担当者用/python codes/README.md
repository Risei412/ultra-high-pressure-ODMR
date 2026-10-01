# 120--400 GPa NV 励起波長: 実験担当者用再現パッケージ

このフォルダは、120 GPa の静水圧光学基準と \(\alpha\) 不変という帰無仮説を
再現するパッケージである。追加の偏差応力補正は対抗仮説であり、結果を以下に併記する。

## 前提

- Python 3.10 以上
- `numpy`, `scipy`, `matplotlib`, `pytest`

PowerShell では、このフォルダで次を一度実行する。

```powershell
python -m pip install -r requirements.txt
```

## 再現手順

```powershell
cd codes
python report_111_anisotropic.py --pressure 120 --alpha 0.60 --axial-stress 120
python fig_alpha_invariance.py
python calc_alpha_corrected_optimum.py
python fig_alpha_optimum_tolerance.py
python fig_microNV_wavelength_decision.py
python fig_microNV_wavelength_decision_300_400gpa.py
python -m pytest ../tests -q
```

最後のコマンドで回帰試験を実行する。図は `figures/` に PNG と PDF で生成される。
レーザー選定用の図は `alpha_optimum_tolerance_120gpa.png` / `.pdf` である。
試料室内マイクロNVの圧力別レーザー選定図は
`microNV_wavelength_decision_120_200gpa.png` / `.pdf` および
`microNV_wavelength_decision_300_400gpa.png` / `.pdf` である。

## 期待される結果

- 静水圧光学カーネルの数値最大: 440.65 nm
- 報告値: 約 441 nm
- 静水圧カーネル固定モデルの \(\alpha\) 依存: 0 nm
- 偏差応力補正モデル: \(\alpha=0.5\)–0.95 で約 485–445 nm（等 \(\sigma_{zz}\)）
  または約 502–446 nm（等圧縮）
- \(\alpha=0.60\): 475.18 nm または 488.10 nm
- \(\alpha=0.5\)–0.7 の共通5%帯: 470.4–483.3 nm または487.3–492.6 nm
- 両正規化を含む488 nmの最悪感度ペナルティ: 約8.0%
- \(\alpha=0.6\)–0.7 の条件付き波長域: 120 GPaで466.1–488.1 nm、
  200 GPaで430.0–514.3 nm
- 300 GPaの外挿域: 380.4–510.0 nm、保守的イオン化壁を課すと405.2–510.0 nm
- 400 GPaの外挿域: 346.1–508.4 nm、保守的イオン化壁を課すと405.2–508.4 nm
- 120–400 GPaの共通運用候補: 488 nm

## 重要な解釈上の注意

Ho の吸収カーネルは \(\alpha=1\) の静水圧計算である。Ho の
\(\alpha\approx0.95\) はマイクロピラー検証条件であり、カーネル自体の条件ではない。
図の不変性は、全ての \(\alpha\) に同じ静水圧カーネルを使った帰無仮説の性質である。

Hilberer の434/769の比を等圧縮で比較するか等 \(\sigma_{zz}\) で比較するかにより、
\(k_q/k_h=-0.734\) または \(-0.392\) となる。共有図は両者を系統幅として表示する。
これは ZPL シフト駆動分の条件付き推定であり、異方応力下の最適励起波長の
直接検証ではない。

200 GPa以上では、静水圧吸収カーネルも \(\alpha\) 依存吸収カーネルも未測定である。
従って300–400 GPa図の青帯はモデル外挿の包絡であり、信頼区間や感度最適帯ではない。
緑帯は、現行モデルで最も長波長側となる405.2 nmのイオン化壁を課した保守的部分を
示す。異方応力下の感度は計算しない。

## 感度計算の由来

波長依存吸収は既存の `ho_spectrum_model.py` と `ho_fig1e_absorption.csv` を使う。
固定入射パワーにおける光学限界では

\[
A(\lambda)=\lambda\sigma_{\mathrm{abs}}(\lambda),\qquad
\frac{\eta(\lambda)}{\eta_{\mathrm{opt}}}
=\sqrt{\frac{A_{\max}}{A(\lambda)}}
\]

として感度比を計算する。`fig_alpha_optimum_tolerance.py` の新規仮定は、この既存の
静水圧感度曲線を各 \(\alpha\) の最適波長まで剛体シフトする点である。従って5%帯も
条件付きであり、異方応力下でsideband形状が変われば実測値に置き換える必要がある。

これ以上の理論精密化は行わず、複数の \(\alpha\) で 440–505 nm を走査して
photon rate、ODMR contrast、linewidth、NV⁻/NV⁰ 比、発光スペクトル、局所応力を
同時測定する。詳細は
`theory_limits_and_required_measurements.md` を参照する。
