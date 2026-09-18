# ultra-high-pressure-ODMR

[![CI](https://github.com/Risei412/ultra-high-pressure-ODMR/actions/workflows/ci.yml/badge.svg)](https://github.com/Risei412/ultra-high-pressure-ODMR/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Author: Risei Abe (Institute of Science Tokyo)

This is the working repository for NV-center ODMR under ultra-high pressure, including the
zero-phonon-line study at high pressure and high temperature. It is public and MIT-licensed:
the analysis code, the tests that pin it, and the manuscripts all live here.

A companion repository,
[NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature](https://github.com/Risei412/NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature),
is a short landing page for the zero-phonon-line study specifically. **This repository is the
one to read for code and results.**

Code-level notes are in [`code/README.md`](code/README.md).

## What is pinned here

The numerical model is frozen and checked against published work rather than against itself:

- [`code/repro_literature.py`](code/repro_literature.py) reproduces results from six published
  sources (Doherty 2014, Ho 2026, Dai 2022, Bhattacharyya 2022, Lyapin 2018, Hilberer 2023).
  The two quantities the model was **fitted** to are marked `[CAL]` in the output, so calibration
  anchors are never presented as independent predictions.
- [`code/repro_yield.py`](code/repro_yield.py) keeps absorption cross-section and detected PL
  separate, which are easy to conflate.
- [`code/dreau_exponent.py`](code/dreau_exponent.py) grounds its assumptions on published data at
  the level of individual equation numbers in Dréau 2011.
- `archive/` is kept as history and is **not** the current implementation. Every figure has the
  code that generates it.

## Reproduce

Install the pinned environment and run the tests:

```bash
python -m pip install -r code/requirements.lock
python -m pytest code/
```

| What | Command | Count | Time |
|---|---|---|---|
| Full suite | `python -m pytest code/` | 243 tests | ~2 min 30 s |
| Frozen values only (fast signal) | `python -m pytest code/tests/test_freeze.py` | 68 tests | ~5 s |

Measured 2026-09-19 on Python 3.14.7 (see `.python-version`) with the versions in
[`code/requirements.lock`](code/requirements.lock). CI runs both jobs on every push.

`code/requirements.txt` keeps loose lower bounds for casual use; **`code/requirements.lock` is
the environment any reproducibility claim about this repository refers to.**

Agreement across environments is reported at a stated tolerance, not bit-for-bit: different CPUs
and BLAS kernels legitimately differ in the last bits. Each CI run therefore logs its CPU, platform
and numpy build configuration, so a disagreement can be attributed rather than guessed at.

## Links

- Researchmap: https://researchmap.jp/risei-abe
- Companion landing page: https://github.com/Risei412/NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature

---

# ultra-high-pressure-ODMR

作成者: 阿部李星（東京科学大学）

このリポジトリは、超高圧下の NV センター ODMR、および高圧・高温下でのゼロフォノン線の研究の作業用本体です。
**公開リポジトリ（MIT ライセンス）** であり、解析コード・それを固定するテスト・原稿はすべてここにあります。

関連リポジトリ
[NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature](https://github.com/Risei412/NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature)
は、ゼロフォノン線の研究に絞った短い紹介ページです。**コードと結果を見るならこのリポジトリのほうです。**

コードの説明は [`code/README.md`](code/README.md) にあります。

## 何を固定しているか

数値モデルは凍結され、自己参照ではなく**公表された研究**に対して照合されています。

- [`code/repro_literature.py`](code/repro_literature.py) は公表された6件（Doherty 2014、Ho 2026、
  Dai 2022、Bhattacharyya 2022、Lyapin 2018、Hilberer 2023）を再現します。モデルの**フィッティングに
  使った**2つの量は出力中に `[CAL]` と明示され、較正アンカーが独立予測として提示されることはありません。
- [`code/repro_yield.py`](code/repro_yield.py) は混同されやすい吸収断面積と検出 PL を分離しています。
- [`code/dreau_exponent.py`](code/dreau_exponent.py) は仮定を Dréau 2011 の式番号の水準で公表データに
  接地しています。
- `archive/` は履歴として保持しているもので、**現行の実装ではありません**。図はすべて生成コード付きです。

## 再現手順

固定環境をインストールしてテストを実行します。

```bash
python -m pip install -r code/requirements.lock
python -m pytest code/
```

| 対象 | コマンド | 件数 | 所要時間 |
|---|---|---|---|
| フルスイート | `python -m pytest code/` | 243 テスト | 約 2 分 30 秒 |
| 凍結値のみ（高速シグナル） | `python -m pytest code/tests/test_freeze.py` | 68 テスト | 約 5 秒 |

2026-09-19 に Python 3.14.7（`.python-version` 参照）と
[`code/requirements.lock`](code/requirements.lock) の各バージョンで実測。CI は push のたびに両方を実行します。

`code/requirements.txt` は簡易利用向けの緩い下限指定です。**再現性の主張が参照するのは
`code/requirements.lock` のほう**です。

環境をまたいだ一致は bit-for-bit ではなく**許容誤差**で報告します。CPU や BLAS 実装が異なれば
末尾の桁は正当に食い違うためです。したがって CI は毎回 CPU・プラットフォーム・numpy のビルド構成を
記録し、食い違いが出たときに推測ではなく切り分けができるようにしています。

## リンク

- Researchmap: https://researchmap.jp/risei-abe
- 関連リポジトリ（紹介ページ）: https://github.com/Risei412/NV-center-Zero-phonon-Line-at-high-pressure-and-high-temperature
