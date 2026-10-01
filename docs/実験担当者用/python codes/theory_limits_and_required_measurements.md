# 120 GPa 励起波長理論の限界と必要な実測値

## 現時点の計算結果

Ho et al. の Fig. 1(e) は静水圧（\(\alpha=1\)）の吸収カーネルである。固定光パワー、
低励起、波長に依存しない電荷収率・ODMRコントラスト・線幅を仮定すると、

\[
\lambda_{\mathrm{opt}}
=\arg\max_\lambda[\lambda\sigma_{\mathrm{abs}}(\lambda,120\ \mathrm{GPa})]
=440.65\ \mathrm{nm}
\]

となる。元データの節点間隔を考慮した報告値は約441 nmである。

## \(\alpha\) 一般化の二仮説

\[
\sigma_{xx}=\sigma_{yy}=\alpha\sigma_{zz},\qquad \sigma_{zz}=120\ \mathrm{GPa}
\]

1. **静水圧カーネル固定仮説**：全ての \(\alpha\) に同じカーネルを用いるため、
   最適波長は全域で約441 nmとなる。これはモデル上の恒等式であり、物理的な
   不変性の証明ではない。
2. **偏差応力補正仮説**：Ho の440.65 nmを \(\alpha=1\) の絶対基準とし、
   Hilberer et al. から得た \(k_q/k_h=-0.734\) による現象論モデルの差分を加える。

\[
\lambda_{\mathrm{opt}}(\alpha)=440.65\ \mathrm{nm}
+[\lambda_{\mathrm{opt}}^{\mathrm{NV}}(\alpha)
-\lambda_{\mathrm{opt}}^{\mathrm{NV}}(1)].
\]

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

一般的なDAC領域 \(\alpha=0.5\)–0.7 では約475–502 nm、文献対応を含む
\(\alpha=0.5\)–0.95 では約446–502 nmとなる。

## 文献との対応と限界

- Ho の吸収カーネルは \(\alpha=1\) の絶対基準を与える。
- Ho の \(\alpha\approx0.95\) はマイクロピラー検証条件であり、同論文の
  \(\alpha=1\) と0.95のZPL計算は準静水圧側の応力応答を拘束する。
- Hilberer の通常型DAC（\(\alpha\approx0.56\)）とマイクロピラー
  （\(\alpha\approx0.95\)）の比較は偏差応力係数を拘束する。

これらはZPL応力応答との整合性を与えるが、\(\lambda_{\mathrm{opt}}(\alpha)\) を
直接測定していない。異方応力によるHuang–Rhys因子、フォノンサイドバンド形状、
電荷収率、ODMRコントラストおよび線幅の変化は未校正である。

## 必要な実測

120 GPa付近で、少なくとも \(\alpha\approx0.5,0.6,0.7\)、可能なら
文献対応点の0.56、0.95と静水圧基準1.0を測る。各条件で440–505 nmを走査し、
少なくとも以下を同時に保存する。

- 試料位置での励起波長と入射パワー
- photon rate \(R\)
- ODMR contrast \(C\) と linewidth \(\Delta\nu\)
- NV⁻/NV⁰比と発光スペクトル
- ZPL位置、幅、分裂およびsideband形状
- 局所 \(\sigma_{zz},\sigma_{xx},\sigma_{yy}\) または \(\alpha\)
- 圧力、温度、DAC内の測定位置

\[
\eta(\lambda,\alpha)\propto
\frac{\Delta\nu(\lambda,\alpha)}
{C(\lambda,\alpha)\sqrt{R(\lambda,\alpha)}}.
\]

二仮説を判別する主量は、最適波長の絶対値と
\(d\lambda_{\mathrm{opt}}/d\alpha\) である。

## 理論の停止条件

現状の曲線は実験計画には十分だが確定値ではない。未測定パラメータを追加する
理論精密化はここで停止する。複数 \(\alpha\) の波長走査、または異方応力下の
方向分解吸収スペクトル／第一原理計算が得られた時点でモデルを再開する。

## 参考文献

- K. O. Ho et al., “Optical Stability and Photophysics of NV Centers in
  Diamond up to 120 GPa,” arXiv:2606.02399 (2026).
- A. Hilberer et al., “Enabling quantum sensing under extreme pressure:
  Nitrogen-vacancy magnetometry up to 130 GPa,” Phys. Rev. B 107, L220102
  (2023), doi:10.1103/PhysRevB.107.L220102.
