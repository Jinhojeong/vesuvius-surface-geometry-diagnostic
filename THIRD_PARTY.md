# Third-party code in this repository

This repository is MIT licensed, Copyright (c) 2026 Jinho Jeong. One file carries code that
is not mine, listed below with its own licence.

## TAUIL-Abd-Elilah/vesuvius-repro

`scripts/shell_split/geom.py` transcribes the margin and shell definitions from
https://github.com/TAUIL-Abd-Elilah/vesuvius-repro at commit `9afa412`. The point of the
transcription is that the label-side geometry has to be his exactly, otherwise our shell
numbers and his are not comparable. Provenance per function:

| function in `geom.py` | source in his repo |
|---|---|
| `across_sheet_dirs` | `thin_labels.py:43` |
| `relabel_margin` | `margin_relabel.py:44` |
| `distance_profile` | `m7_margin_fp.py:97` |
| `analyse` | `m7_margin_fp.py:136` |

`margin_mask` is a two-line composition of his `relabel_margin` output, matching
`m7_margin_fp.margin_mask`. `descriptives` is mine and has no counterpart in his code.

His licence, reproduced in full:

```
MIT License

Copyright (c) 2026 TAUIL Abd Elillah

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
