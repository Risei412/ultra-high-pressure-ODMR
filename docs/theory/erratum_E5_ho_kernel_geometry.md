# 正誤表 E5 — 静水圧最適波長 441 nm を異方応力 DAC へ適用するときの暗黙の仮定

日付: 2026-09-01(**改訂 R1、同日**)
対象: **Addendum A4 の適用範囲**、および `docs/experiment/theory_limits_and_required_measurements.md`
不変: 凍結本体の光学限界表、**A4.1 の数値**、Addendum A1、A2、A3、正誤表 E1–E4、S1
コード: `code/nv_model.py`(`_alpha_factor`, `ZPL_DEV_K`)、`code/ho_odmr_sensitivity.py`
根拠: `docs/references/notes/Ho2026.md` L29–30, L43–45、Hilberer et al., PRB 107, L220102 (2023)
文献確認: 2026-09-01 実施（**改訂 R2** で刷新）。抽出: `code/extract_ho_fig1a_zpl.py`、
`code/data/ho_fig1a_zpl.csv`。計算: `code/calc_alpha_corrected_optimum.py`
参照 PDF（取得済み）: `docs/references/papers/PhysRevB.107.L220102.pdf`（Hilberer）、
`Davies-OpticalStudies1945-1976.pdf`、`2204.05064v1.pdf`（Dai）、`hxtk-vmjq.pdf`（Huang）

---

> **⚠ 改訂 R1 — 初版 §1 は撤回する**
>
> 初版は「Ho カーネルは α ≈ 0.95 のカーネルであり、440.65 nm は静水圧の答えでは
> ない」と主張した。**これは誤りである。**
>
> `Ho2026.md` L29–30 のとおり、Ho らの DFT は有限圧力をダイヤモンドの状態方程式
> (Occelli 2003)による**格子定数のスケーリング**で実現している。したがって
> \(q=0\) は構成上厳密であり、Fig. 1(e) の吸収スペクトル(7 圧力)から再構成した
> \(\sigma_{\rm abs}^{\rm Ho}(\lambda,P)\) は**真の静水圧カーネル**である。
> α ≈ 0.95 は Fig. 5(b) の**検証実験**の DAC の性質であって、カーネルの性質では
> ない。初版はカーネルの出所と検証の出所を取り違えていた。
>
> **したがって「441 nm は間違い」という初版の主張は強すぎる。撤回する。**
> 440.65 nm は静水圧 120 GPa の最適波長として正しく、その位置づけを与えた
> A4.1 の記述も正しい。
>
> 本改訂版が主張するのは、これより**はるかに狭い**一点のみである(§0)。

---

> **⚠ 改訂 R2（2026-09-01）— R1 の「§3 解決」を撤回する**
>
> R1 は、WebFetch 経由で読んだ Hilberer の本文から「横軸は圧縮モル体積だから
> 等圧縮比較、よって k_q/k_h = −0.734 で確定」とした。**確定していない。**
>
> 掲載版 PDF（`PhysRevB.107.L220102.pdf`）を直接読むと、体積は**独立に測られた
> 量ではなく**、「アンビル先端の校正済みダイヤモンド Raman フォノンモード」で測った
> 圧力を静水圧 EOS に入れて得た導出量である。したがって「等体積」は操作的に
> 「等 Raman ゲージ圧」であり、これは等 P_mean より等 sigma_zz に近い。
> つまり**読み A の方が操作的には自然**であり、R1 の結論は逆向きに行きすぎていた。
>
> §3b の三者一致も、再検すると**両方の読みを含む窓**を与えるだけで、
> 選別しない（k の整合窓は −0.89 〜 −0.35）。
>
> **係数は k_q/k_h ∈ [−0.734, −0.391] のままである。**推定帯は R1 初版の
> **466–502 nm** に戻す。ただし minimax 固定線 **473 nm** は、仮説集合を
> 7 つに拡げても変わらない（§7）。結論の実務部分は影響を受けない。
>
> 一方、掲載版の確認で**強化された**点もある（§3）。

---

## 0. 結論（改訂 R2）

> **440.65 nm（報告値 約 441 nm）は静水圧 120 GPa の最適波長として正しい。**
> 誤っているのは数値ではなく、**適用**である。A4 はこの静水圧の値を
> 「sigma_zz = 120 GPa、alpha = 0.5–0.7 の平面キュレット実験の設計事前値」
> として渡した。この変換には k_q/k_h = +2/3 という仮定が必要であり、
> **実測値は負で、符号が逆である。**
>
> 平面キュレットの最適波長の推定は **466–502 nm**（k の不確かさを含む）、
> 実用固定線は **473 nm（minimax 解、最悪 ×1.205）**。
> これは A4 の 441 nm を置き換える確定値ではなく、**同格の対抗仮説**である。

三版の差分。

| | 初版 | 改訂 R1 | **改訂 R2（現行）** |
|---|---|---|---|
| 441 nm 自体 | 誤り | 正しい | **正しい**（静水圧の値として） |
| A4.1 の数値表 | 撤回 | 存続 | **存続** |
| 誤りの所在 | カーネルの幾何 | 変換の欠落 | **変換の欠落** |
| k_q/k_h | 不確定 | −0.734 で確定 | **[−0.734, −0.391]。未確定** |
| 推定帯 | 466–502 nm | 475–502 nm | **466–502 nm** |
| 固定線 | 473 nm（定性） | 473 nm（minimax, 4 仮説） | **473 nm（minimax, 7 仮説）** |

**文献確認で確定したこととしなかったこと。**

| | 状態 |
|---|---|
| k_q の**符号** | **確定。**Davies & Hamer の物理、Hilberer の2幾何比較、Ho の Fig. 1(a) の三者が負を支持 |
| k_q への**帰属** | **確定。**Hilberer 自身が差を「deviatoric stress reduction」に帰している（§3） |
| alpha の圧力一定性 | **確定。**0.56 も 0.95 も圧力に対して一定（§3） |
| Ho の alpha = 0.95 計算の存在 | **確定。**Fig. 1(a) に両曲線。抽出済み（§3b） |
| k_q の**大きさ** | **未確定。**約 2 倍の幅が残る（§3、§3b） |
| 固定線 473 nm | **強化された。**仮説を 7 つに拡げても minimax のまま（§7） |

---

## 1. Ho カーネルは静水圧である

`docs/references/notes/Ho2026.md` §2「手法・理論」:

> - スピン分極 DFT(VASP)。有限圧力はダイヤモンドの状態方程式(Occelli 2003)で
>   格子定数を調整して実現 → **式は静水圧的に構成されている**。

格子定数の等方スケーリングは定義上 \(q=0\) である。よって

\[
\sigma_{\rm abs}^{\rm Ho}(\lambda,P),\qquad P=P_{\rm mean},\quad q\equiv0
\]

であり、第二引数は曖昧さなく静水圧である。**441 nm は静水圧 120 GPa の答えとして
正しい。**

同ノート §4(a) の Fig. 5(b) 検証(532 / 457 nm、0–120 GPa)は、この計算カーネルが
実測と一致することを示す。ただしその実験は α ≈ 0.95 の DAC で行われており、
そこでの偏差応力は

\[
q=(1-0.95)\times120\ {\rm GPa}=6\ {\rm GPa}
\]

にすぎない。**Fig. 5(b) はカーネルを \(q\approx6\) GPa の近傍でのみ検証している。**
平面キュレットの \(q=36\)–60 GPa については何も言っていない。これは E5 が
主張できる範囲を決める重要な制約である。

---

## 2. 「\(\sigma_{zz}=120\) GPa に \(P=120\) GPa を代入」は選択である

ここが本正誤表の実質であり、§1 の訂正では消えない。むしろカーネルが真に静水圧
だと確定したことで、論点は**より鮮明になる**。

カーネルは \(q=0\) の線上でしか定義されていない。\(q\neq0\) の実験に適用するには、
実効圧力を選ばなければならない。光学応答を

\[
\Delta E \;=\; k_h P_{\rm mean} + k_q\,q,
\qquad
P_{\rm eff} \;=\; P_{\rm mean} + \frac{k_q}{k_h}\,q
\]

と書くと、軸対称応力の恒等式

\[
\sigma_{zz}=P_{\rm mean}+\tfrac{2}{3}q,
\qquad
P_{\rm mean}=\frac{1+2\alpha}{3}\sigma_{zz},
\qquad
q=(1-\alpha)\sigma_{zz}
\]

より、A4 が行った「\(\sigma_{zz}=120\) GPa の条件に \(P=120\) GPa のカーネルを
当てる」操作は

\[
\frac{k_q}{k_h}=+\frac{2}{3}
\]

と**数値的に同一**である。すなわち A4 の \(\alpha\) 不変性は「異方応力について
何も仮定しない」結果ではなく、「偏差応力は静水圧成分の 2/3 の重みで**同じ向きに**
効く」という仮定の結果である。

| 選択 | \(k_q/k_h\) | 由来 |
|---|---:|---|
| **A4 の適用(暗黙)** | **+0.667** | \(\sigma_{zz}\) をカーネル座標に代入 |
| A4.4 が退けた naive な代入 | 0 | \(P_{\rm mean}\) をカーネル座標に代入 |
| Hilberer 実測(等 \(\sigma_{zz}\) 読み) | −0.391 | `nv_model.ZPL_DEV_K` |
| Hilberer 実測(等圧縮読み) | −0.734 | §3 |

Davies & Hamer 以来知られるとおり、偏差応力は ³E を**赤方偏移**させ、静水圧成分の
青方偏移を打ち消す。よって \(k_q/k_h<0\) が期待され、Hilberer の2ジオメトリ実測も
そう言っている。**A4 の暗黙値は符号が逆である。**

なお A4.4 が naive な \(P_{\rm mean}\) 代入を「未校正だから採用しない」として
退けたのは、理由としては正しいが不十分だった。退けた選択(0)と採用した選択
(+2/3)のどちらも実測値(負)の側になく、**選択肢を比較する枠組み自体が
欠けていた**。これが A4 の実際の欠陥である。

---

## 3. 正規化は未確定 — k_q/k_h ∈ [−0.734, −0.391]

Hilberer らは同一の ZPL シフトを2つの幾何で測っている。

- \(\alpha=0.56\)(標準平面キュレット): \(-434\pm2\) meV/(cm³ mol⁻¹)
- \(\alpha=0.95\)(FIB 加工マイクロピラー): \(-769\pm4\) meV/(cm³ mol⁻¹)

比 \(434/769=0.564369\)。**2つの幾何を何を揃えて比較したか**で \(k_q/k_h\) が変わる。

**読み A(等 \(\sigma_{zz}\))。**`nv_model.py` の実装規約。
\(f(\alpha)=(1+2\alpha)/3+k(1-\alpha)\) として \(f(0.56)/f(0.95)=0.564369\) を解くと
\(k_q/k_h=-0.391\)。これは `ZPL_DEV_K = -0.39211` を 0.2% で再現する。

**読み B(等圧縮)。**単位が meV/(cm³ mol⁻¹)、すなわち**体積あたり**であることを
素直に取る読み。等方弾性体では体積歪 \(\varepsilon\simeq P_{\rm mean}/K\) なので
比較は等 \(P_{\rm mean}\) で行われる。\(q/P_{\rm mean}=3(1-\alpha)/(1+2\alpha)\) より

\[
\frac{1+0.6226\,k}{1+0.0517\,k}=0.564369
\quad\Longrightarrow\quad
k_q/k_h=-0.734 .
\]

**確定すること。**どちらの読みでも \(k_q/k_h<0\)。符号は物理と測定の両方から
決まる。したがって平面キュレットの最適波長は 441 nm より**赤側**にあると期待される。

**文献確認の結果（2026-09-01、掲載版 PDF）。**

> "While the NV ZPL dependence with pressure is not linear, its evolution
> becomes linear when plotted versus the compressed diamond volume estimated
> using the diamond equation of state [29]. A linear fit gives a slope of
> −769 ± 4 meV/(cm³ mol⁻¹). A similar measurement performed on a nonmodified
> diamond anvil yields a weaker slope of −434 ± 2 meV/(cm³ mol⁻¹).
> **This significant difference in the pressure dependence of the ZPL is
> another indication of the deviatoric stress reduction caused by the
> microstructuration of the anvil tip.**"

### 強化された点

1. **著者自身が差を偏差応力に帰属している。**初版は「原論文は理由を説明して
   いない」と書いたが、**これは誤り**。掲載版は上記のとおり明示している。
   k_q への帰属は本正誤表の独断ではなく、著者の解釈である。
2. **alpha は圧力に対して一定である。**「α = 0.56 that is essentially constant with
   pressure」「α ≃ 0.95 that stays constant within the pressure range tested」。
   これは q/P_mean = 3(1−alpha)/(1+2alpha) が負荷経路上一定であることを意味し、
   読み B の導出が必要とする前提を裏付ける。
3. アンビルは (100)。磁場は「the diamond [100] crystal axis for which all NV centers
   have an equivalent response to stress and magnetic field」に印加されている。

### しかし正規化は確定しない

初版は「横軸が体積なのだから等圧縮比較、よって読み B」とした。これは体積が
独立に測られている場合にのみ成立する。掲載版によれば:

> "Pressure in the DAC was measured using the **calibrated diamond Raman phonon
> mode at the anvil tip** [18,19]."

さらに Fig. 2(b) のキャプションは「as a function of **pressure and diamond
volume**」であり、体積は測定圧力から静水圧 EOS（Shen et al. [29]）を介して
**導出された量**である。したがって「等体積」は操作的には**等 Raman ゲージ圧**で
あり、Akahama–Kawamura の Raman エッジ校正は通常 sigma_zz（法線応力）に対して
与えられる。つまり**等 sigma_zz に近い——読み A 側である**。

さらに同論文自身が、ゲージ自体が応力状態依存であることを述べている:

> "Under hydrostatic conditions, the dependence of the frequency of the Raman
> scattering with diamond volume follows a Grüneisen relation of parameter
> γ = 0.97(1) **whereas the frequency shift is smaller under deviatoric
> stress**."

静水圧校正を平面キュレット側に適用すれば、推定体積に系統誤差が入る。
その大きさを評価するには Supplemental が要る。

**結論。**正規化は**未確定のまま**であり、

\[ k_q/k_h \in [-0.734,\ -0.391] \]

として扱う。確定しているのは**符号が負**であることだけであり、これは
Davies & Hamer（`Davies-OpticalStudies1945-1976.pdf`）の偏差応力赤方偏移と
Hilberer の2幾何比較の両方から得られる。

**`nv_model.py` について。**初版は `ZPL_DEV_K = -0.39211` を「バグ」と断定したが、
上記のとおり等 sigma_zz 読みは操作的に正当でありうる。**バグ判定を取り下げる**。
残る問題は、コメントが「at equal compression」と書きながら実装が等 sigma_zz である
**文言の不整合**であり、これは修正すべきである。

---

## 3b. Ho 自身が alpha = 0.95 を計算していた

**A4 の最大の見落としはここにある。**Ho らは異方応力の場合を自分で計算しており、
Fig. 1(a) にその曲線を載せている。凡例は `expt.` / `thr. alpha = 1` / `thr. alpha = 0.95`。
本文の逐語:

> "We calculated the ZPL position for perfect hydrostatic stress (orange solid
> line) and **the lower branch of the splitted ZPL for alpha = 0.95** (orange
> dashed line). ... modeling with alpha = 0.95 reproduces the lower branch of
> the experimentally observed split ZPL with an excellent agreement."

つまり**吸収カーネルを生んだのと同じ水準の理論で、偏差応力の光学結合が
2点校正されている**。他実験から輸入する必要がない唯一の校正である。

ベクター図から抽出した（`code/extract_ho_fig1a_zpl.py`、15 頂点、曲線適合なし）:

| 量 | 値 |
|---|---:|
| dZPL(alpha=1, 120 GPa) | 0.4651 eV |
| dZPL(alpha=0.95, 120 GPa) | 0.4306 eV |
| P_mean, q | 116, 6 GPa |
| 静水圧曲線の P_mean での値 | 0.4545 eV |
| q に帰属する差 | **−23.9 meV** |
| 静水圧等価圧力 | 107.3 GPa |
| 生の k_effective | **−1.453** |

**抽出の検証。**alpha=1 の 120 GPa 値 0.4651 eV は、正誤表 E4 が
Fig. 1(b),(e) と Pekarian 停留条件という**完全に別経路**で導いた
dE_120 = 0.464 eV と 0.25% で一致する。校正は健全である。

### 生の −1.453 を [111] にそのまま使ってはならない

Ho の破線は **(100) アンビルにおける split ZPL の下枝**である。したがって差
−23.9 meV は、対称（A1）な偏差シフトと、3E 軌道分裂の半分との**和**である。
[111] 軸対称応力下の [111] 配向 NV では、**分裂は対称性により消える**
（`theory_limits_and_required_measurements.md`「期待できること」）。よって [111]
に効くのは A1 成分のみである。

### 三者の独立な一致

| 寄与 | 出所 | q=6 GPa での値 |
|---|---|---:|
| A1 偏差シフト | Hilberer 読み B (k = −0.734) | −11.7 meV |
| 3E 分裂の半分 | 残差 | −12.1 meV → delta = 24.3 meV |
| 3E 分裂の半分 | `analysis_C_lambda.splitting_meV`（独立） | delta = 19.2–36.5 meV |
| **合計** | | **−23.9 meV（Ho の実測差）** |

`analysis_C_lambda` の分裂モデルは NV の軌道歪み感受率とダイヤモンドの弾性定数
だけから作られており、Ho の図も Hilberer の測定も使っていない。その予測
19.2–36.5 meV（弾性率に E=1100 GPa を使うか C44=578 GPa を使うかの幅）は、
残差 24.3 meV を含む。**Ho の DFT、Hilberer の測定、リポジトリ独自の弾性モデルの
三者が、それぞれの不確かさの範囲で整合する。**

### この一致は両読みを選別しない（改訂 R2）

R1 はこの三者一致を読み B の裏付けとして提示したが、分裂モデルの幅をそのまま
伝搬させると A1 残差の幅はこうなる:

| 弾性率の選択 | delta | A1 残差 | 対応する k |
|---|---:|---:|---:|
| E = 1100 GPa | 19.2 meV | −14.3 meV | **−0.893** |
| C44 = 578 GPa | 36.5 meV | −5.6 meV | **−0.351** |

**整合窓は −0.89 〜 −0.35 であり、読み A（−0.391）も読み B（−0.734）も両方含む。**
したがってこの検証が示すのは、**枠組み（符号、桁、各項が足し合わさること）が
正しい**ことであって、係数の値を一意に決めることではない。

逆に、分裂 delta を実測すればこの表は k を一意に定める。**これが§9 の新しい
事前登録項目の根拠である。**

なお (100) 幾何では分裂の寄与が加わるので、同じ alpha でも実効的にさらに赤側へ
寄る。**アンビル配向によって答えが変わる。**

---

## 3c. Huang et al. — [111] は対称性保存であり、コントラストは alpha → 0 で最大

`docs/references/papers/hxtk-vmjq.pdf`。PRL 137, 093801 (2026)。一般応力下の NV 光学
特性を ab initio で扱い、(100)/(110)/(111) の3配向で実験している。

### 規約はリポジトリと一致する

> "σ = α σ_hyd + (1 − α) σ_[111], where α characterizes the degree of
> hydrostaticity ... (i.e., α = 0, 1 respectively)"

alpha = 1 が静水圧、alpha = 0 が一軸 [111]。`nv_model_111.axisymmetric_111_stress` の
規約と一致しており、同関数の docstring はすでにこれを引用している。

### §3b の前提が独立に裏付けられた

Huang は [111] 一軸応力を**対称性保存（symmetry-preserving）**応力に分類し、
対称性を破る応力（一軸 [100] など）と明確に区別している。すなわち
**[111] 軸対称応力下では 3E の分裂が生じない**。これは §3b が
「Ho の生の −1.453 を [111] に使ってはならない、[111] に効くのは A1 成分のみ」と
した根拠を、独立な ab initio 研究から裏付ける。

さらに、(111) カットアンビルのコントラスト向上は

> "intrinsic to the [111]-oriented NV itself, ruling out previous
> interpretations based on the 'darkening' of non-[111] NVs"

であり、向きの揃った NV そのものの性質である。

### 実験の alpha はハンドオフ範囲を挟む

Huang の比較対象は alpha ≈ 1（(100) ナノピラー）、**0.73**（(111) カット、自前の実験）、
**0.57**（(111) カット）。ハンドオフの alpha = 0.5–0.7 はこの実測範囲に入っており、
§5 の補正を適用する領域として妥当である。

### しかし Huang の推奨は alpha → 0 である

> "we motivate the deliberate engineering of stress environments to promote
> metrological sensitivity, i.e., **α = 0**, which may be obtained with a
> uniaxial press"

機構は、横方向 SOC λ⊥ は両方の圧縮で増える一方、振動重なり F(Δ) が
**逆の振る舞いをする**ことである——一軸 [111] 歪みでは増大し、静水圧歪みでは抑制される。
結果として上部 ISC レート、したがってコントラストは alpha → 0 で最大になる。

**これはこの分野の通常の戦略（Hilberer のマイクロピラーによる静水圧の回復）と
逆向きである。**alpha は「キュレットが決めてしまう受動的な雑音」ではなく、
**設計変数**である。

### 本正誤表と合わせると新しい緊張が出る

感度は eta ∝ Δν / (C √R)。

- **Huang（C の項）**: alpha を下げるほどコントラストが上がる
- **E5（R の項）**: alpha を下げるほど lambda_opt が赤側へ逐げる

つまり、**コントラストが最適になる応力状態こそ、励起波長を 441 nm から最も大きく
動かさなければならない状態である**。どちらの論文も単独では述べていない。

### そしてそこでは本モデルが壊れる

モデルの ZPL シフト因子 f(alpha) = (1+2alpha)/3 + k(1−alpha) は、alpha を下げると
符号を変える:

| 読み | f(alpha) = 0 となる alpha | 意味 |
|---|---:|---|
| 等圧縮（k = −0.734） | **0.286** | これ以下で 120 GPa の ZPL が常圧より低エネルギーになる |
| 等 sigma_zz（k = −0.391） | **0.056** | 同上 |

Hilberer の校正は alpha = 0.56 と 0.95 の2点だけであり、その外への線形外推は
ここで物理的に破綻する。**したがって本正誤表の補正は Huang の推奨する
alpha → 0 を評価できない。**`code/calc_alpha_corrected_optimum.py` に
`cancellation_alpha()` を追加し、それ以下の alpha を拒否するようにした。

### §3 を決める最良の経路が入れ替わった

> "for symmetry-preserving stress, σ⊥ = ½(σxx + σyy) and σ∥ = σzz
> **play qualitatively different roles in shifting energy gaps**"

二係数構造そのものである。係数の値は Supplemental Material にあり、本文にはない。
しかしこれは **ab initio で、幾何を分解して** 与えられており、Hilberer の2点実測を
正規化の曖昧さごと解釈するより**直接的**である。
**§8 項目 1b（Hilberer SM）より Huang SM を先に取得すべきである。**

---

## 4. なぜ \(P_{\rm eff}\) をカーネルに直接代入して答えを出せないか

負の結合を入れると \(P_{\rm eff}\) は \(\alpha=0.5\)–0.7 で 36–82 GPa まで落ちる。
そこでは Ho カーネルの大域 argmax が使えない。

| \(P\) [GPa] | Ho カーネル argmax | `nv_model` |
|---:|---:|---:|
| 50 | 578.35 nm | 521.8 nm |
| 55 | 557.70 nm | 517.3 nm |
| 65 | 557.70 nm | 508.9 nm |
| **70** | **480.15 nm** | 505.0 nm |
| **75** | **540.70 nm** | 501.3 nm |
| 80 | 540.70 nm | 497.9 nm |
| 85 | 464.15 nm | 494.6 nm |
| 90 | 461.50 nm | 491.4 nm |
| 100 | 458.65 nm | 485.6 nm |
| 110 | 447.00 nm | 480.3 nm |
| 120 | 440.65 nm | 475.5 nm |

70 GPa で 480 nm、75 GPa で 540 nm。**5 GPa で 60 nm 逆行する。**これは物理では
なく、E2.1 の4つの極大の間を argmax が飛び移っているだけである(A3 定理 X の
分岐交代領域、E3 が \(P^*\) を物理量として撤回した領域)。

**85–90 GPa 以上では両者とも単調**で、約 25–35 nm のほぼ一定オフセットを保って
並走する。したがって α 補正をカーネルの圧力座標経由で行ってはならない。

---

## 5. 対抗仮説の構成 — 差分補正

- **絶対値**はカーネルが有効な 120 GPa の静水圧アンカー 440.65 nm から取る
- **alpha 依存性**は校正済みで滑らかな `nv_model` から**差分**として取る
- カーネルが**理想静水圧**なので、基準は alpha = 1.0（初版の 0.95 ではない）
- 係数は未確定なので**両読みを並走させる**

実装: `code/calc_alpha_corrected_optimum.py`

    lambda_opt(alpha) = 440.65 nm + [ lambda_opt^nv(alpha) - lambda_opt^nv(1.0) ]  at 120 GPa

| alpha | P_mean | q | 読み A（k = −0.391） | 読み B（k = −0.734） |
|---:|---:|---:|---|---|
| 0.50 | 80 GPa | 60 GPa | +44.0 → **484.6 nm** | +60.8 → **501.5 nm** |
| 0.56 | 84.8 GPa | 52.8 GPa | +38.3 → 478.9 nm | +52.7 → 493.4 nm |
| 0.60 | 88 GPa | 48 GPa | +34.5 → 475.2 nm | +47.4 → **488.1 nm** |
| 0.70 | 96 GPa | 36 GPa | +25.4 → **466.1 nm** | +34.7 → 475.4 nm |
| 0.95 | 116 GPa | 6 GPa | +4.1 → 444.7 nm | +5.4 → 446.1 nm |
| 1.00 | 120 GPa | 0 | 0 → 440.6 nm | 0 → 440.6 nm |

**alpha = 0.5–0.7 の推定範囲は両読みを合わせて 466–502 nm**。alpha 依存性は
18.5（読み A）〜 26.1 nm（読み B）。

**整合性の確認。**alpha = 0.95 では補正が +4〜+5 nm にとどまる。Ho の Fig. 5(b)
検証がその幾何（q = 6 GPa）で行われてカーネルと一致したことと矛盾しない。
本構成は**検証済みの領域では小さな補正しか与えず、未検証の領域（q = 36–60 GPa）
でのみ大きくなる**。

**この構成の仮定。**`nv_model` の alpha = 1.0 での絶対値は 471.5 nm で、Ho カーネルの
440.65 nm と 30.9 nm ずれている。そのオフセットが差分で相殺することを仮定している。
§4 の表で 85 GPa 以上の両曲線がほぼ一定オフセットで並走することが根拠だが、
証明ではない。また E4.3 のとおり単一モード Pekarian は極大を１つしか持たないので、
この差分は「ZPL シフトが駆動する分」のみを表す。

温度には依存しない（90 K と 300 K で差分は 0.1 nm 以内）。

---

## 6. A4 のどこが無効になり、どこが残るか(改訂版)

初版は 6 項目を撤回したが、**そのうち 2 項目は撤回を取り消す**。

### 存続(初版の撤回を取り消すもの)

| 箇所 | 内容 | 状態 |
|---|---|---|
| A4.1 | 「静水圧設計事前値 約 441 nm」 | **正しい。**Ho の DFT は静水圧構成 |
| A4.1 | 5% 帯 426.43–457.90 nm、ペナルティ表 | **正しい。**静水圧の帯として有効 |

### 存続(初版から不変)

| 箇所 | 内容 |
|---|---|
| A4.0 | 応力テンソルと不変量の計算 |
| A4.3 | 許容帯 = A2 の level set(\(I/I_c=f^2\))。α と無関係な恒等式 |
| A4.4 | 分岐交代の警告。§4 がその実例で、**強化された** |
| A4.5 | 実測項目表 |
| A4.6 | 未解決事項の列挙 |

### 撤回(縮小後)

| 箇所 | 内容 | 理由 |
|---|---|---|
| A4.0 報告規則8 | 「零シフトはコードが返せる唯一の値」 | 零は \(k_q/k_h=+2/3\) という選択の結果。`nv_model` は校正済みの非零値を返す |
| A4.0 | 「\(\sigma_{zz}=120\) GPa, \(\alpha=0.5\)–0.7 で 441 nm」 | 静水圧の値を無条件に異方応力条件へ移している。§2 |
| A4.2 | 固定線 445 nm、判定余裕 4.60% | 静水圧帯の内部で最適化された選択。対抗仮説の推定帯(466–502 nm)から 21–57 nm 外 |
| A4.7 | 「理論は実験のために凍結された」 | 対抗仮説が未解決のまま残っている |

**撤回されるのは A4 の数値ではなく、A4 の適用範囲の記述である。**

---

## 7. 訂正後の設計事前値

| 量 | A4 | E5 改訂 R2 |
|---|---:|---:|
| 静水圧 120 GPa の最適 | 441 nm | **441 nm（不変）** |
| 静水圧 5% 帯 | 426.4–457.9 nm | **不変** |
| alpha = 0.5–0.7 の最適 | 441 nm（0 nm 変化） | **未確定。**441 nm と 466–502 nm の二択 |
| 実用固定線 | 445 nm | **473 nm** |
| 対照線 | 532 nm | 532 nm（据置） |

### 473 nm は minimax 解であり、仮説集合の拡大に対して頑健

固定線は「正しいから」ではなく、**どの仮説が真でも最悪損失が最小**になることを
基準に選ぶ。帯が剛体的に移動するとして、**7 仮説**（静水圧 + 読み A/B ×
alpha = 0.7/0.6/0.5）に対する光学限界ペナルティの最悪値:

| 線 | 最悪ペナルティ |
|---:|---:|
| 457 nm | 1.510 |
| **473 nm** | **1.205** |
| 476 nm | 1.246 |
| 480 nm | 1.384 |
| 488 nm | 1.538 |
| 505 nm | 3.021 |
| 514.5 nm | 5.108 |

450–510 nm を 0.25 nm 刻みで走査した真の minimax は **472.50 nm（最悪 ×1.198）**であり、
市販レーザー線の **473 nm はそれを 0.6% しか損なわない**。

**この結果は R1（4 仮説）と R2（7 仮説）で完全に一致する。**係数の 2 倍の
不確かさは、固定線の選択を変えない。

**改訂 R2 で取り下げる主張。**R1 は「514.5 nm（ZPL）が次善」としたが、
これは読み B のみの狭い仮説集合によるアーテファクトであった。読み A を含めると
最悪 ×5.108 となり、**最悪の候補になる**。撤回する。

---

## 8. 確認項目の状態（Huang 反映後）

| # | 項目 | 状態 |
|---|---|---|
| 1 | Hilberer の 434/769 の正規化 | **未解決。**体積は Raman ゲージ圧から EOS で導出されており、等 sigma_zz 側である可能性が高い（§3） |
| **1a** | **Huang Supplemental**（未取得） | **新規・最優先。**sigma⊥ と sigma∥ の係数を ab initio で幾何分解して与えている。これが k_q を直接決める（§3c） |
| 1b | Hilberer Supplemental | 降格。体積の導出手順。1a が取れれば不要になる |
| 2 | Ho の異方応力計算の有無 | **解決。**Fig. 1(a)。抽出済み（§3b） |
| 3 | Ho Supplemental S3 | 未取得。優先度中 |
| 4 | 実験のアンビル配向 | **実質解決。**(111) カットが正解。コントラスト向上は [111] 配向 NV 固有であり、[111] 軸対称は対称性保存（§3c） |
| 5 | q が ZPL 以外に与える効果 | **部分解決。**Huang が ISC レートとコントラストの alpha 依存を与えた。sideband 形状は依然未解決 |
| 6 | Huang et al. 本文 | **読了（§3c）** |
| **7** | **alpha をどこに置くか** | **新規。**Huang は alpha → 0 を推奨するが、本モデルは alpha ≤ 0.286 で破綻する（§3c） |

---

## 9. 実験への含意

**判別は容易である。**仮説は 25–61 nm 離れており、簗い波長走査でも分離できる。

| 仮説 | alpha = 0.5–0.7 の lambda_opt | lambda_opt の alpha 依存 |
|---|---|---|
| A4（k = +2/3） | 441 nm | 0 nm |
| E5 読み A（k = −0.391） | 466–485 nm | 約 −93 nm / 単位 alpha |
| E5 読み B（k = −0.734） | 475–502 nm | 約 −131 nm / 単位 alpha |

測るべきは lambda_opt の絶対位置だけでなく、**alpha に対する傾き**である。
傾きは絶対較正を必要としないので判別力が高い。なお傾きの測定は A4 と E5 を
分けるだけでなく、**読み A と B も分ける**（−93 対 −131 nm/alpha）。

線リストは **441 / 457 / 473 / 488 / 505 / 532 nm**、走査域は最低でも
**440–510 nm** を覆う。§4 の分岐交代があるため帯の外側も捨てない。

**事前登録項目（新規、波長走査不要）。**alpha を変えたときの **ZPL 分裂 delta を
同時に記録する**こと。§3b の表のとおり、delta を実測すれば A1 残差が確定し、
**k が一意に決まる**。すなわち分裂の測定が、波長シフトの予測を同じ試料から
独立に与える。A1.2 の T4 と同じ性格を持つ。

---

## 10. 再現

```bash
cd code

# Ho カーネルの argmax が 85 GPa 以下で非単調(§4)
python -c "
from ho_odmr_sensitivity import HoIntegratedODMRModel
m=HoIntegratedODMRModel()
for P in range(50,125,5): print(P, round(m.optimum(float(P)),2))"

# Hilberer 比の2つの読み(§3)
python -c "
from scipy.optimize import brentq
R=434.0/769.0
fA=lambda a,k:(1+2*a)/3.0+k*(1-a)
print('A', brentq(lambda k: fA(0.56,k)/fA(0.95,k)-R,-3,0.9))
qom=lambda a:3*(1-a)/(1+2*a)
print('B', brentq(lambda k:(1+k*qom(0.56))/(1+k*qom(0.95))-R,-1.3,0.9))"

# 差分補正、alpha=1.0 基準(§5)
python -c "
from nv_model import NVModel
ref=NVModel(T=300.,alpha=1.0).lambda_opt(120.)
for a in (0.5,0.56,0.6,0.7,0.95):
    print(a, round(440.65+NVModel(T=300.,alpha=a).lambda_opt(120.)-ref,1))"
```

---

## 11. 状態

**理論探索は二仮説を事前登録した状態で再凍結する。**

- **静水圧層は凍結を維持する。**441 nm、5% 帯、ペナルティ表、A1–A3、A4.3 は
  そのまま有効であり、`code/tests/test_freeze_a4_handoff.py` の静水圧側ロックも維持する。
- **異方応力層は未確定として凍結する。**A4 の 441 nm 不変モデルを帰無仮説、
  E5 の \(\alpha=0.5\)–0.95 に対する 502–446 nm 曲線を対抗仮説とする。
  445 nm、473 nm その他の単一線を、実測前に普遍的な固定線として確定しない。

§8 項目1（\(k_q/k_h=-0.734\)）は達成された。ただし文献が直接拘束するのは
ZPL応力応答であり、\(\lambda_{\mathrm{opt}}(\alpha)\) ではない。異方応力下の
Huang–Rhys因子とsideband形状を仮定で補う追加探索は行わない。再開条件は、
複数 \(\alpha\) における120 GPa波長走査、または異方応力下の方向分解吸収
スペクトル／第一原理計算の取得とする。

反映状況:

1. `docs/experiment/theory_limits_and_required_measurements.md` と実験担当者用コピーに、
   二仮説、\(\alpha=0.5\)–0.95 の表、必要な実測、停止条件を反映済み。
2. `code/calc_alpha_corrected_optimum.py` に、既存モデルを変更しない独立した差分補正を
   実装済み。\(\alpha=0.6\) で488.10 nmを回帰試験する。
3. `theory_freeze_v3_ho_integrated.md` と旧A4回帰試験は静水圧層の履歴として保持する。
   E5を反映した全面更新は、実測後のv4で行う。

---

## 参考文献

取得済みの PDF は `docs/references/papers/` にある。

- K. O. Ho et al., "Optical Stability and Photophysics of NV Centers in
  Diamond up to 120 GPa," arXiv:2606.02399 (2026). — `2606.02399v1.pdf`
- A. Hilberer et al., "Enabling quantum sensing under extreme pressure:
  Nitrogen-vacancy magnetometry up to 130 GPa," Phys. Rev. B 107, L220102
  (2023), doi:10.1103/PhysRevB.107.L220102. — `PhysRevB.107.L220102.pdf`
  （掲載版。§3 の逐語確認に使用）
- G. Davies and M. F. Hamer, "Optical studies of the 1.945 eV vibronic band in
  diamond," Proc. R. Soc. Lond. A 348, 285 (1976). —
  `Davies-OpticalStudies1945-1976.pdf`（§2 の偏差応力赤方偏移）
- B. Huang et al., "Elucidating the Intersystem Crossing of the
  Nitrogen-Vacancy Center up to Megabar Pressures," Phys. Rev. Lett. 137,
  093801 (2026), doi:10.1103/hxtk-vmjq. — `hxtk-vmjq.pdf`
  （§3c。一般応力下の ISC、[111] の対称性保存、alpha → 0 の推奨）
- J.-H. Dai et al., "Quantum sensing with diamond NV centers under megabar
  pressures," arXiv:2204.05064. — `2204.05064v1.pdf`
  （100 GPa での PL のパワー線形性。`test_freeze.py` の根拠）
- F. Occelli, P. Loubeyre, R. LeToullec, "Properties of diamond under
  hydrostatic pressures up to 140 GPa," Nature Materials 2, 151 (2003).
- G. Shen et al., High Press. Res.（Hilberer の ref. [29]、ダイヤモンド EOS）.
- Y. Akahama and H. Kawamura, J. Appl. Phys. 96, 3748 (2004)（Raman 圧力ゲージ）.
