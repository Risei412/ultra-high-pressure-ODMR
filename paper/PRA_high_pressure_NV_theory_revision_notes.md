# 高圧NV-ODMR理論原稿：PRA投稿前の構造・論理修正メモ

対象原稿：`main(2).pdf`  
目的：実験結果を追加する前に、**内容・構成上の欠陥を先に潰し、その後に英語表現・細部の推敲へ移る**ための修正指示書。

---

## 0. 総合判定

現状の原稿は、PRA向けの理論原稿として**骨格は成立している**。

中心構造は

1. **Theorem G**：response structure / identifiability
2. **Theorem M**：power-driven level set / multiplicity ladder
3. **Theorem X**：pressure-driven branch exchange

の三層で整理されており、

- response
- power
- pressure

という三つの軸で理論が積み上がっている。

また、

- Fig. 1：gauge degeneracy
- Fig. 2：level set と multiplicity ladder
- Fig. 3：pressure-driven branch exchange
- prospective / pre-specified tests
- limitations
- kernel validation
- reproducibility code

まで用意されており、**全面的な構成変更は不要**。

ただし、現段階ではまだ「英語表現だけを直す段階」には完全には移らない方がよい。  
まず以下の **A1–A3 を必須修正**、B1–B3をscope整理する。

---

# A. 修正必須

## A1. Eq. (6) の `⇔` が強すぎる

現在の主張：

\[
\lambda_\eta=\lambda_{\mathrm{PL}}
\Longleftrightarrow
\ell_G(\lambda_{\mathrm{PL}})=0
\]

ただし、

\[
G=C/\Delta\nu
\]

である。

### 問題

\[
\ell_G(\lambda_{\mathrm{PL}})=0
\]

から直接言えるのは、基本的には

> \(\lambda_{\mathrm{PL}}\) が \(\Phi=G^2R\) の stationary point にもなる

というところまで。

それだけでは、その点が \(\Phi\) の **global maximum** であることは保証されない。

したがって、`⇔` をそのまま使うなら additional global condition が必要。

### さらに注意

Introduction / Abstract では、

> photoluminescence maximizing wavelength と sensitivity optimum が一致する条件は既に破れている

という強い書き方になっている。

しかし引用している先行実験で直接測られているのは主に **contrast \(C\)** の波長依存であり、

\[
G=C/\Delta\nu
\]

なので、linewidth \(\Delta\nu\) が同じ logarithmic slope を持てば相殺は原理的に可能。

### 推奨修正

主張を次のレベルに弱める：

> The measured wavelength dependence of the ODMR contrast makes coincidence of the photoluminescence and sensitivity optima non-generic. Whether \(G=C/\Delta\nu\) is actually stationary at \(\lambda_{\rm PL}\) must be tested by measuring the linewidth as well.

あるいは Eq. (6) を

- stationary-point condition
- global-optimum condition

に分離する。

### 優先度

**最優先。**

Paper opening の premise に直接関わるため、表現修正ではなく論理修正。

---

## A2. Theorem M の仮定を定理 statement に完全に入れる

現在の Theorem M は主に

- mediation hypothesis (M)
- \(\Phi\) が unique interior maximum \(\Gamma_p^\*\) を持つ

ことから level set と multiplicity ladder を導いている。

一方、proof では途中から

> Regard \(A\) as a Morse function on the analysis interval.

が導入されている。

### 問題1：level set と staircase は別レベルの主張

(M) から

\[
\Lambda_\eta(I)
=
\{\lambda:A(\lambda)=\Gamma_p^\*/(\gamma I)\}
\]

という **level-set structure** が出ることと、

\[
\Delta N=
+2,-2,-1
\]

という **cardinality staircase** が出ることは分けた方がよい。

後者は Morse-type regularity を必要とする。

### 推奨構造

#### Theorem M1 — Level-set theorem

Assume:

1. (M)
2. \(\eta=h(\Gamma_p)\)
3. \(h\) has a unique global minimum at \(\Gamma_p^\*\)

Then

\[
\Lambda_\eta(I)
=
\{\lambda:A(\lambda)=\Gamma_p^\*/(\gamma I)\}.
\]

#### Corollary M2 — Multiplicity staircase

さらに \(A(\lambda)\) が compact interval 上の Morse function で、critical values が nondegenerate なら

- interior maximum crossing: \(\Delta N=+2\)
- interior minimum crossing: \(\Delta N=-2\)
- analysis-window edge crossing: \(\Delta N=-1\)

となる。

### 問題2：\(I<I_c\) の扱い

現在は

> required level exceeds \(A_{\max}\), so the minimum is the single point \(\lambda_{\rm abs}\)

としている。

これは一般の \(h(\Gamma_p)\) に対しては、

- \(\Gamma_p<\Gamma_p^\*\) で \(\Phi\) が単調増加
- あるいは \(\eta\) が単調減少

することが必要。

power-law model では成立するが、一般定理なら仮定に入れるべき。

### 問題3：calibration-free ladder の low-power relation

現在の

\[
R(\lambda)\propto A(\lambda)
\]

を使って low-power PL scan から critical values \(a_k\) を読む部分は、(M) 単独からは出ない。

必要なのは別途、

> linear optical-response / unsaturated low-power regime

という条件。

### 推奨

Theorem M の中心を

> **level set is structural**

に置き、

- Morse staircase
- low-power PL reconstruction
- calibration-free rung prediction

を corollary として積む。

### 優先度

**必須修正。**

---

## A3. Theorem X の excitation bandwidth を effective bandwidth として定義する

現在、

\[
A_{\rm SB}/A_{\rm ZPL}=r(P)W
\]

として、

\[
r(P^\*)W=1
\]

から branch-exchange pressure \(P^\*(W)\) を出している。

### 問題

ZPL がどれだけ励起されるかは単純な laser linewidth のみで決まるわけではない。

少なくとも物理的には

\[
W_{\rm eff}
\sim
\max(W_{\rm laser},\Gamma_{\rm ZPL})
\]

あるいはより正確には

> laser spectral profile と intrinsic ZPL lineshape の convolution

で定義される effective spectral sampling width を使うべき。

原稿自身の Limitations でも、

> If the laser linewidth exceeds the ZPL width ...

と認識している。

### 現状の問題点

\(\Gamma_{\rm ZPL}(P)\) は現在使っている source figure から reliable に読めない。

したがって

> the experimenter selects the exchange pressure through the excitation bandwidth

は現状では少し強い。

### 推奨修正

以下のように regime を明示する：

> In the regime where the excitation bandwidth exceeds the intrinsic ZPL width, the exchange pressure is tunable through the laser bandwidth.

または

\[
W\rightarrow W_{\rm eff}
\]

と全面的に置き換える。

### 実験追加後の強化案

実験で

- \(\Gamma_{\rm ZPL}(P)\)
- laser linewidth
- \(P^\*\)

を同時に測れれば、Theorem X の実験的価値がかなり増す。

### 優先度

**必須修正。**

---

# B. 構成・scope の整理を推奨

## B1. Physical ladder と finite-window ladder を区別する

Fig. 2 / Table III の multiplicity ladder には、

\[
I/I_c\simeq1.86
\]

付近の \(\Delta N=-1\) rung が含まれている。

これは physical absorption edge ではなく、

> reconstruction の source figure が 402 nm で切れている

ことによって生じる **analysis-window edge**。

### 問題

この rung を他の physical critical points と同列に置くと、

> source-figure truncation が物理予言の一部

のように見える。

### 推奨

区別する：

#### Physical multiplicity events

- sideband maxima/minima
- ZPL branch feature

#### Censoring / window event

- 402 nm analysis-window boundary

Fig. 2 / Table III では

- 灰色
- dashed
- `window-censoring event`

などで明確に分ける。

可能なら main-text の「物理 ladder」から除外し、Appendix で finite-window counting として扱う。

---

## B2. ZPL height uncertainty を main claim から切り離す

現在、ZPL peak height uncertainty によって ladder sequence が

\[
2\to4\to6\to4\to3\to5\to3
\]

から別 sequence に変わり得る。

原稿はこの点をかなり正直に扱えているが、main claim を reconstruction-specific ZPL height に背負わせる必要はない。

### 推奨

中心主張を

> critical-value crossings generate a staircase

に固定し、

具体的な 120 GPa sequence は

> worked example conditional on the reconstructed kernel

として一段下げる。

実験では T5 により measured spectrum から \(a_k\) を読むため、この不確実性は本質的でないことを強調。

---

## B3. T1 / T3 / T6 は `(M) holds` ではなく falsification-oriented にする

現在の Table V には、

- collapse ⇒ (M) holds
- equal ⇒ pure (M)
- exactly equal ⇒ (M)

のような表現がある。

### 問題

有限の observable / finite precision での一致は、一般には (M) の厳密な証明ではない。

別の photon-energy-specific mechanism が測定範囲内だけ同じ manifold に乗る可能性は残る。

### 推奨表現

- `consistent with (M)`
- `supports effective mediation by the pump rate`
- `failure falsifies (M)`
- `deviation quantifies the breakdown of (M)`

などに変更。

### より良い論理構造

(M) は

> **strongly falsifiable, but not uniquely provable from a finite set of observables**

と位置付ける。

これにより理論の慎重さが増し、むしろ査読耐性が上がる。

---

## B4. Theorem G の generality を少し絞る

同じ

\[
(E,n)
\]

を持つ response model が同じ \(\eta\) surface を与える、という gauge-degeneracy は非常に良い結果。

ただし composite mechanisms については

- common response scale \(\Gamma\)
- one scale dominant
- scales coincide

などの条件がある。

### 推奨

Abstract / theorem statement では

> within the common-scale power-law response class

など、domain を明示する。

---

## B5. `E=3, n=2` の書き方を Ref. [25] 限定にする

現在の

> Published power-dependence data already fix \(E=3,n=2\).

は高圧 NV response 一般を fix したようにも読める。

実際には、

> Ref. [25] の single-NV CW-ODMR response model / fit

から

\[
E=3,\quad n=2
\]

を読む、という位置付け。

### 推奨

> The power-dependence model fitted to the single-NV measurements of Ref. [25] corresponds to \(E=3\) and \(n=2\).

とする。

実験結果が入った後、自分の高圧条件で T7 により \(E,n\) を再評価する。

---

# C. 構成について

## 現在の大枠は維持してよい

推奨する最終構造：

1. Introduction
2. Framework and the two optima
3. Response exponent / identifiability — Theorem G
4. Level set and multiplicity ladder — Theorem M
5. Pressure-driven branch exchange — Theorem X
6. Experimental protocol / predictions
7. Experimental results
8. Discussion
9. Conclusions

Theory-only draft の現状では

> G → M → X

の順序を維持する。

---

## 実験追加時の原則

Experimental Results も theory と同じ順番にする：

### Result 1 — G

- \(R(I)\)
- \(C(I)\)
- \(\Delta\nu(I)\)
- \(E,n\)
- sensitivity identifiability

### Result 2 — M

- measured low-power \(R(\lambda)\)
- measured \(a_k\)
- predicted \(I_k/I_c\)
- observed multiplicity / degeneracy

### Result 3 — X

- pressure-dependent branch positions
- coexistence
- branch exchange
- bandwidth dependence of \(P^\*\)

これにより theory と experiment が一対一対応する。

---

# D. Prospective tests の圧縮

現在の T1–T11 は理論 draft としては有用だが、実験結果が入った final paper では少し多い。

### 推奨

main text では4本程度に圧縮：

1. **Mediation-collapse / degeneracy test**
2. **Exponent closure test**
3. **Multiplicity-ladder prediction**
4. **Pressure-driven branch-exchange test**

残りは Appendix / Supplemental Material へ。

また formal preregistration を外部で行っていない場合、

`pre-registered tests`

より

- `pre-specified tests`
- `prospective tests`
- `falsification tests`

の方が安全。

---

# E. その他の構造上の注意

## Acquisition time の scaling

Introduction の一部で acquisition time と sensitivity の関係が逆に読める箇所がある。

磁場感度を \(\eta\) とする通常の定義では、

\[
t_{\rm map}\propto \eta^2
\]

である。

Conclusion ではこの形になっているため、Introduction と統一する。

---

# F. Notion PRA level 判定

Notion の内部基準：

- **PRA-1**：complete system-specific physics
- **PRA-2**：mechanism / observable / experimental diagnostic
- **PRA-3**：reusable method / protocol
- **PRA-4**：subfield-general framework / formalism

## 現状

### 判定

**PRA-2.5 ～ PRA-3−**

より言葉で言うと：

> **PRA-3候補だが、現状の evidence closure は PRA-2 相当**

### PRA-2 を超えている理由

中心成果は単なる

> 120 GPa では 440–457 nm が良い

という system-specific optimization ではない。

原稿には、

- level-set theorem
- multiplicity ladder
- calibration-free rung ratios
- mechanism non-identifiability
- experimental null tests
- pressure-driven branch exchange
- bandwidth-dependent exchange signature

があり、**再利用可能な解析・測定原理**を目指している。

### まだ PRA-4 ではない理由

現状では nontrivial application が実質的に

> high-pressure NV-DAC

の一系に集中している。

さらに、

- optical kernel：一つの高圧文献 reconstruction に依存
- response：既存 single-NV model に依存
- 第二の独立な physical platform の実証なし

なので、

> subfield-general formalism demonstrated across multiple nontrivial systems

にはまだ届いていない。

---

# G. 実験結果による PRA level の分岐

## ケース1

単に

> 457 nm excitation が 532 nm より高感度

を実証。

→ **PRA-1～PRA-2**

---

## ケース2

- (M) の null test
- \(C,R,\Delta\nu\) の個別測定
- low-power spectrum から ladder を事前予測
- predicted rung / degeneracy を観測

まで成功。

→ **明確な PRA-3**

---

## ケース3

さらに

- pressure-driven branch exchange
- branch coexistence
- laser/effective bandwidth を変えると \(P^\*\) が系統的に動く

まで実証。

→ **強い PRA-3**

---

## ケース4

別の absorbing quantum sensor / color center /第二モデルでも同じ formalism を非自明に再現。

→ **PRA-4 候補**

---

# H. 修正作業の優先順位

## Phase 1 — Physics / logic freeze

1. Eq. (6) の global / stationary condition 修正
2. Theorem M を level-set theorem と Morse corollary に分離
3. \(I<I_c\) の monotonicity assumption を明示
4. low-power \(R\propto A\) 条件を明示
5. Theorem X の \(W\) を effective bandwidth として整理
6. analysis-window edge と physical critical point を分離
7. (M) tests を falsification-oriented に書き換え
8. Theorem G / \(E=3,n=2\) のscopeを限定

この段階で**理論構造をfreeze**する。

---

## Phase 2 — Manuscript architecture

1. main-text tests を圧縮
2. 実験 section を G → M → X の順で挿入
3. reconstruction-specific numerical detail を Appendix へ整理
4. physical claim と worked-example claim を明確に分離

---

## Phase 3 — 表現修正

Phase 1–2 が終わった後で、

- Abstract
- Introduction
- theorem statements
- Discussion
- Conclusion

の順で英語をPRA向けに研磨する。

---

# 最終判断

現状の原稿は**大改造が必要な状態ではない**。

ただし、

> 「内容・構成に欠陥がないので、あとは表現を直すだけ」

と判断するにはまだ早い。

現在は、

> **かなり良い理論骨格に残っている数本の論理ボルトを締める段階**

である。

A1–A3とscope整理を終えた時点で、physics architecture を凍結してよい。
