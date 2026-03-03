# smart-SXN-screener

A Python-based toolkit for natural product mass spectrometry data analysis, covering theoretical library construction, flavonoid screening, ginkgolide screening, and diagnostic ion-based compound matching.

---

## Module 1: Flavonoid Theoretical Library Builder & Candidate Compound Matcher

A tool for constructing a theoretical flavonoid compound library and matching it against experimental mass spectrometry data using PPM-based filtering.

### Key Features

- **Theoretical Library Generation**: Combines aglycones, sugars, and acyl groups to generate glycosides and acylated glycosides
- **Exact Mass Calculation**: Uses high-precision atomic masses for [M+H]⁺ or [M-H]⁻ ion calculation
- **PPM Matching**: Matches theoretical masses against experimental m/z values within a user-defined PPM threshold
- **Deduplication**: Groups results by molecular formula and merges compound names for clean output
- **Flexible Ion Mode**: Supports both positive ([M+H]⁺) and negative ([M-H]⁻) ionization modes

### Installation
```bash
pip install pandas numpy
```

### Input Format

**aglycones.csv / sugars.csv / acyls.csv** (columns: `Name`, `Formula`):
```csv
Name,Formula
Quercetin,C15H10O7
Kaempferol,C15H10O6
```

**Experimental data CSV** (must contain `NO.` and `m/z` columns):
```csv
NO.,m/z
1,447.093400
2,301.035600
```

### Configuration

Edit the `CONFIG` dictionary at the top of the script:
```python
CONFIG = {
    "aglycones_file":    "path/to/aglycones.csv",
    "sugars_file":       "path/to/sugars.csv",
    "acyls_file":        "path/to/acyls.csv",
    "experimental_file": "path/to/experimental.csv",
    "output_file":       "path/to/output.csv",
    "ion_type":          "M-H",   # "M+H" or "M-H"
    "ppm_threshold":     10
}
```

### Usage
```bash
python module1.py
```

### Output Columns

| Column | Description |
|---|---|
| `NO.` | Experimental peak index |
| `Experimental_mz` | Measured m/z value |
| `Molecular_Formula` | Matched molecular formula |
| `Compound_Name` | All matched compound names (semicolon-separated) |
| `All_Compound_Count` | Number of candidate structures |
| `Theoretical_Mass` | Calculated theoretical m/z |
| `PPM_Difference` | Mass accuracy (ppm) |

---

## Module 2: Flavonoid & Flavonol Analog Screener

A tool for identifying flavonoid and flavonol analogs from LC-MS/MS peak tables by matching aglycone masses and detecting characteristic B-ring fragment ions.

### Key Features

- **Dual Matching Strategy**: Matches aglycone m/z against both precursor ions and MS/MS fragment ions
- **Flavonoid vs. Flavonol Classification**: Distinguishes flavonoids from flavonols using B-ring hydroxylation fragment series
- **Confidence Scoring**: Reports `高度可信 / 可信 / 可能 / 未知` confidence levels based on fragment evidence
- **OH Count Determination**: Estimates the number of B-ring hydroxyl groups from matched fragments

### Installation
```bash
pip install pandas numpy
```

### Input Format

**Aglycone information CSV** (columns: `Name`, `Formula`, `aglycone m/z`):
```csv
Name,Formula,aglycone m/z
Aglycone1,C15H9O5,269.04507
Aglycone2,C15H9O6,285.03999
```

**Peak table CSV** (must contain `m/z`, `RT`, `MS/MS` columns):
```csv
m/z,RT,MS/MS
579.17145,5.23,57.04081:26 153.01824:150 269.04507:800
```

### Configuration

Edit inside the `FlavonoidScreener.__init__` method:
```python
self.mass_tolerance = 0.01
self.aglycone_file  = "path/to/aglycones.csv"
self.peak_file      = "path/to/peaks.csv"
self.output_dir     = "path/to/results/"
```

### Usage
```bash
python module2.py
```

### Fragment Ion Reference

| Ion m/z | Type | OH groups |
|---|---|---|
| 137.02332 | Flavonoid B-ring | 1 |
| 153.01824 | Flavonoid B-ring | 2 |
| 169.01315 | Flavonoid B-ring | 3 |
| 185.00806 | Flavonoid B-ring | 4 |
| 149.02332 | Flavonol B-ring | 1 |
| 165.01824 | Flavonol B-ring | 2 |
| 181.01315 | Flavonol B-ring | 3 |
| 197.00806 | Flavonol B-ring | 4 |

### Output Columns

| Column | Description |
|---|---|
| `Peak_m/z` | Precursor ion m/z |
| `RT` | Retention time |
| `Aglycone_Name` | Matched aglycone identifier |
| `Match_Type` | `母离子匹配` or `碎片匹配` |
| `Flavonoid_Type` | `黄酮` or `黄酮醇` or `未知` |
| `Confidence` | `高度可信 / 可信 / 可能 / 未知` |
| `OH_Count` | Estimated B-ring hydroxyl count |
| `Matched_Fragments_Info` | Detailed fragment match summary |

---

## Module 3: Ginkgolide Analog Screener

A tool for automated screening of ginkgolide-type compounds from LC-MS/MS negative-mode peak tables, using a sequential neutral-loss and fragment-ion decision workflow.

### Key Features

- **Two-stage Screening**: Filters by co-occurrence of characteristic neutral loss (72.99312 Da) and CO/2CO losses, then subclassifies by fragment ions
- **Detailed Classification**: Labels each hit as a known ginkgolide subtype or analog
- **Transparent Output**: Adds diagnostic flag columns for every screening criterion

### Installation
```bash
pip install pandas
```

### Screening Logic
```
Input: precursor m/z + MS/MS fragments
        │
        ▼
[Stage 1] Has neutral loss 72.99312 AND (CO loss OR 2×CO loss)?
        │ NO  → "非银杏内酯类化合物"
        │ YES ↓
[Stage 2] Has fragment m/z 125.02442?
        │ NO  → Has CO₂/H₂O neutral loss? → YES: "银杏内酯类似物"
        │                                    NO:  "未分类银杏内酯"
        │ YES ↓
[Stage 3] Has fragment m/z 141.01933?
             YES → "R₁=OH和R₃=OH型银杏内酯"
             NO  → "可能为R1=OH,C-C16为双键的银杏内酯类似物"
```

### Key Mass Values

| Parameter | Value (Da) |
|---|---|
| CO neutral loss | 27.99546 |
| 2×CO neutral loss | 55.99092 |
| CO₂ neutral loss | 43.99038 |
| H₂O neutral loss | 18.01111 |
| CO₂+H₂O neutral loss | 62.00094 |
| Characteristic neutral loss | 72.99312 |
| Diagnostic fragment 1 | 125.02442 |
| Diagnostic fragment 2 | 141.01933 |
| Mass tolerance | ±0.01 Da |

### Configuration

Edit the file path variables inside `main()`:
```python
input_file  = r"path/to/peak_table.csv"
output_file = r"path/to/ginkgolide_results.csv"
```

Input CSV must contain: `m/z`, `RT`, `MS/MS` columns.

### Usage
```bash
python module3.py
```

### Added Output Columns

| Column | Description |
|---|---|
| `是否为银杏内酯` | `是` / `否` |
| `分类结果` | Subtype classification label |
| `第一级筛选` | Stage 1 pass/fail |
| `是否有特征丢失和CO/2CO丢失` | Stage 1 detail |
| `是否有125.02442` | Stage 2 detail |
| `是否有141.01933` | Stage 3 detail |
| `是否有中性丢失(CO2/H2O/两者)` | Alternative stage 2 detail |

---

## Module 4: Diagnostic Ion-Based Compound Screener

A tool for screening target compound classes from LC-MS/MS peak tables using user-defined diagnostic ions, followed by PPM-based matching against a compound library.

### Key Features

- **Diagnostic Ion Filtering**: Retains only peaks where all specified diagnostic ions co-occur in the MS/MS spectrum
- **Library Matching**: Matches surviving precursor ions to a compound library using a configurable PPM threshold
- **Per-ion Flagging**: Optionally adds a `Has_<ion>` boolean column for each diagnostic ion
- **Automatic RT Detection**: Detects retention time columns from common naming variants
- **Statistical Summary**: Reports match counts, PPM statistics, and RT distribution

### Installation
```bash
pip install pandas numpy
```

### Input Format

**Peak table CSV** (must contain `Precursor m/z` and `MSMS spectrum`):
```csv
Precursor m/z,RT,MSMS spectrum
256.1827,2.34,234.0551:1000 170.0942:800 86.0964:300
```

**Compound library CSV** (must contain `M+H` column):
```csv
Name,Formula,M+H
Compound_A,C10H14N2O3S,243.0747
Compound_B,C12H16N2O4S,285.0853
```

### Configuration

Edit the user configuration section at the bottom of the script:
```python
PEAK_TABLE_PATH        = "path/to/peak_table.csv"
COMPOUND_LIB_PATH      = "path/to/compound_library.csv"
OUTPUT_RESULTS_PATH    = "path/to/matching_results.csv"
OUTPUT_DIAGNOSTIC_PATH = "path/to/diagnostic_peaks.csv"

DIAGNOSTIC_IONS = [234.0551, 170.0942]   # ALL must be present simultaneously
TOLERANCE       = 0.01                    # absolute fragment matching tolerance (Da)
PPM_THRESHOLD   = 10                      # precursor–library matching threshold

MARK_EACH_ION   = True                    # add Has_<ion> columns to output
```

### Usage
```bash
python module4.py
```

### Output Files

**matching_results.csv** — compounds passing both diagnostic ion filter and library PPM match:

| Column | Description |
|---|---|
| `母离子 (Precursor m/z)` | Measured precursor m/z |
| `保留时间 (RT)` | Retention time (if available) |
| `理论M+H` | Theoretical [M+H]⁺ from library |
| `ppm` | Mass accuracy |
| `Has_<ion>` | Per-ion presence flag (if enabled) |
| `MSMS spectrum` | Original MS/MS string |

**diagnostic_peaks.csv** — all peaks passing the diagnostic ion filter.

---

## Requirements

- Python 3.7 or higher
- Windows / macOS / Linux

## Contact

**Version**: v1.0  
**Update Date**: 2025

## License

MIT License
