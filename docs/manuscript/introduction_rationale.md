# Introduction rationale — stress-dependent blue-laser selection

更新: 2026-10-01。主原稿は main 側の設計原稿 `paper/main.tex` とする。
G/M/X を中心とする別原稿は `archive/manuscripts/main_theory_GMX_20260906.tex`
に保存する。

本稿では、青レーザーの選択を局所応力と光学入力に条件づける。

1. 応力比は `alpha = sigma_xx/sigma_zz = sigma_yy/sigma_zz` と定義する。
   同じ軸応力でも alpha が変われば平均応力と偏差応力が変わる。
2. Ho の静水圧120 GPaカーネルは alpha = 1 の基準であり、最適波長は約441 nm。
   このカーネルを全ての alpha にそのまま使う場合の不変性はモデルの仮定である。
3. 異方応力の対抗仮説には現象論モデルの alpha 依存差分だけを加える。
   alpha = 0.60 で475.2–488.1 nm、alpha = 0.95 で444.7–446.1 nmとなる。
   この区間は2つの応力規格化の違いであり、統計的信頼区間ではない。
4. 旧単一モードモデルの475.5 nmと473 nm候補は alpha = 0.95 の旧モデル内の
   設計事前値として残し、静水圧Ho基準や修正モデルの予測と混同しない。
5. alpha 依存図・表を結果の前半に置き、異方応力下の直接測定が未完了である
   ことと、5%帯には曲線の剛体シフトを仮定することを併記する。
6. 実験では励起波長とalphaを変えて photon rate、contrast、linewidth を
   同時測定する。120 GPaの加算基準を別の圧力へ持ち出さない。

計算: `code/calc_alpha_corrected_optimum.py`。
図: `docs/experiment/figures/alpha_optimum_tolerance_120gpa.pdf`。
仮定と出所: `docs/theory/erratum_E5_ho_kernel_geometry.md`。
