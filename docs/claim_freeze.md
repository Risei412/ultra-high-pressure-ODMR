# 主張の凍結（8 分発表版）

対象: 8 分・9 枚のスライド発表。数値はすべて `code/` の実測出力であり、
`paper/main.tex` と一致する。以後、ここに書いた文言を発表準備の基準とする。

---

## 1. 中核主張（1 文・動かさない）

> **120 GPa の NV-ODMR は 474 nm で励起すべきであり、市販の 473 nm でそれが買える。
> 緑（532 nm）比で感度 ×2.7、磁気マップ取得時間 ×7。**

価値の源泉は「明日レーザーを発注できる単一の数字」であること。DAC は 1 回のロードに
数週間かかり 1 回の破損で失われるので、感度は指標ではなく**空間分解能・視野・
(P,T) 点数を買う通貨**である（マップ時間 ∝ 画素数 / η²）。

## 2. 支持（3 つ）

| # | 主張 | 根拠 | 出力元 |
|---|---|---|---|
| S1 | 常圧の 2 択はどちらも窓の外。532 nm は吸収端が追い越した（×2.7）、405 nm は基底状態イオン化端そのもの（×3.7） | ZPL +0.40 eV、IP(³A₂) 3.06 eV = 405 nm | `fig6_three_shifts.py`, `fig7_answer_talk.py` |
| S2 | 答えは吸収極大に一致し、**現象論定数に一切依存しない** | ES イオン化と再結合が σ_abs に比例して f₋ から相殺（0.2%）。9 定数の応答が厳密に 0 nm | `fig4_tornado.py`, `tests/test_freeze.py` |
| S3 | 入力は全部他人の測定値。追加調整なしで**既発表 26 件を再現**。較正に使ったのは Dai 2022 の 3 点（定数 2 個）だけで、独立な Hilberer 2023 マイクロピラーとは **×0.6–1.2** で一致 | Ho 2026 / Doherty 2014 / Dai 2022 / Hilberer 2023 / Bhattacharyya 2022 | `repro_literature.py`, `fig8_repro_contrast.py` |

## 3. 境界（2 つ・主張の一部）

| # | 境界 | 内容 | 出力元 |
|---|---|---|---|
| B1 | **反証条件**：71 GPa の符号変化 | 71 GPa 以下では緑が勝つ。50 GPa の「青の優位なし」という既報 null を再現。符号変化は単一圧力の系統誤差では作れない | `fig5_threshold.py` |
| B2 | **適用条件**：準静水圧・低パワー | micropillar（α≈0.95）が前提（平坦キュレットなら 508 nm、だがそこは 40–50 GPa で ODMR が死ぬ）。u ≳ 0.1 で最適は青へ動き、473 nm 固定のコストが ×1.6 | `repro_literature.py`, `fig3_power_sweep.py` |

## 4. 8 分から落とすもの（バックアップ行き）

検出帯 ×1.6（合計 ×4）／圧力追従 0.4–0.6 nm/GPa と 2 波長運用（488 nm below ~110 GPa,
473 nm above）／type Ib アンビルの青光透過問題／校正計画 Part C。
圧力追従だけは `fig7_answer_talk.py` の図中に 1 行として残し、口頭 5 秒で触れる。

## 5. 時間割と図の対応

| # | 時間 | スライド | 図 |
|---|---|---|---|
| 1 | 0:00–0:40 | 結論先出し：474 nm / 473 nm 市販線 / ×2.7 → ×7 | 数字のみ |
| 2 | 0:40–1:50 | **なぜ 120 GPa で局所磁気イメージングなのか**（専門外向けの動機） | T_c–圧力の表 |
| 3 | 1:50–2:40 | 感度が通貨：測定時間 ∝ 画素数/η² | 図なし |
| 4 | 2:40–3:50 | なぜ自明でないか：2 つの端が窓を挟み込む | `talk_three_shifts_st.png` |
| 5 | 3:50–4:50 | 答え | `talk_answer_120GPa_st.png` |
| 6 | 4:50–5:10 | 現象論非依存：フィットできる量は答えを動かさない | `tornado_lambda_opt_120GPa_st.png` |
| 7 | 5:10–6:10 | **先行実験の再現**（追加調整なし） | `repro_contrast_st.png` |
| 8 | 6:10–7:10 | 反証条件：約 70 GPa の符号変化 | `threshold_green_blue_st.png` |
| 9 | 7:10–8:00 | 適用条件とまとめ | 表のみ |

スライド 2 の論理は `slides/intro_deck.html`（¶1–¶3）を圧縮したもの:
狙う圧力帯 → 抵抗では決着しない → 探針をアンビルの内側に作る。

## 6. 凍結した数値と出所

| 量 | 値 | 出所 |
|---|---|---|
| λ_opt @120 GPa | 474.0 nm（MC 16–84%: +3.7 / −4.3 nm） | `repro_literature.py` |
| 5% / 10% 許容窓 | 461.7–486.2 / 456.6–491.1 nm | `repro_literature.py` |
| η/η_opt @120 GPa | 473: 1.00 / 457: 1.10 / 487: 1.06 / 532: 2.66 / 405: 3.67 | `repro_literature.py` |
| 緑→青 交差圧（473 nm） | 70.6 GPa（457 nm: 82、450 nm: 87） | `fig5_threshold.py` |
| 50 GPa での η(532)/η(473) | 0.59（＝緑が有利、既報の null と整合） | `fig5_threshold.py` |
| σ(473)/σ(532) @120 GPa | 10.5 | `fig6_three_shifts.py` |
| 吸収極大の移動 | 586 → 474 nm | `fig6_three_shifts.py` |
| イオン化端の移動 | 463 → 405 nm | `fig6_three_shifts.py` |
| λ_opt 追従 | 484.1 nm @100 GPa → 474.0 nm @120 GPa（0.51 nm/GPa） | `nv_model.lambda_opt` |
| アンカー応答 | ΔE_ZPL(120): ±3.7 nm、S_abs 勾配: ±2.7 nm、他 9 定数: 0.00 nm | `fig4_tornado.py` |
| 文献再現 | 26/26（OPEN 1 件） | `repro_literature.py` |
| 較正に使った実験点 | Dai 2022 の 3 点のみ（C_amb, E_isc の 2 定数） | `fig8_repro_contrast.py` |
| 独立データとの一致 | Hilberer 2023 micropillar 3 点で ×0.58 / ×0.63 / ×1.16 | `fig8_repro_contrast.py` |

## 7. 発表で言わないと危ない留保

- **常圧側の包絡は近似が粗い**。単一有効フォノンモデルなので 0 GPa の吸収極大が
  586 nm に来る（実際の NV⁻ は ~560 nm）。`talk_three_shifts.png` の灰色曲線は
  **シフトの図示**であって常圧の予言ではない。この理由から、図中の比は常圧比ではなく
  **同一圧力（120 GPa）での σ(473)/σ(532) = 10.5** を使っている。
- **交差圧 70.6 GPa には +2.3 GPa の系統差がある**。モデルの励起重みは光子束一定で
  あり、光パワー一定に直すと 72.9 GPa になる（`docs/code_audit_2026-08.md` B-1）。
  発表では「約 70 GPa」と丸めて述べる。
- **120 GPa 超は外挿**。ZPL 以外のアンカーは 120 GPa で clip される（同 B-3）。
- λ_opt は ³E の応力分裂 20–40 meV に対して 4–7 nm の曖昧さを持つ。窓の端ではなく
  中心に近い線を選ぶ理由がこれ。

## 8. 主張を変えるときの手順

このファイルの数値を変える場合は、必ず以下を再実行して差分を反映すること。

```bash
cd code
python -m pytest tests/ -q          # 67 件
python repro_literature.py          # 26/26
python fig4_tornado.py              # アンカー応答
python fig5_threshold.py            # 交差圧
python fig6_three_shifts.py         # 端の移動
python fig7_answer_talk.py          # 窓とペナルティ
python fig8_repro_contrast.py       # 先行実験の再現
python slides/build_pptx.py         # Science Tokyo テンプレートの 9 枚
```
