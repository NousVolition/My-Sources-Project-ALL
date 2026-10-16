# Numerical data downloads

## Hug gap refinement and reciprocal Lorenz tests — October 15, 2026

[Results and reproduction](reports/clay-gap-refinement/README.md) · [Verified raw-array release](https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/hug-gap-refinement-2026-10-15) · [File catalog](reports/clay-gap-refinement/data-release.json). Code, summaries and the saved-state viewer remain in the source tree. Restore this study with `python reports/clay-gap-refinement/restore_data.py`; the helper verifies every byte and refuses conflicting files. This separate addition leaves the dated migration catalog unchanged.


## HUG unloading refinement — October 15, 2026

[Report and verification instructions](reports/hug-molding-refinement/README.md) · [Release](https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/hug-unloading-refinement-2026-10-15) · [Separate file catalog](reports/hug-molding-refinement/data-release.json). The 16.9 MB ZIP preserves all 156 original refinement and parent-comparison files plus a start page and a historical-link context notice. This new catalog leaves the dated migration inventory below unchanged.

**[Open the verified data release](https://github.com/NousVolition/My-Sources-Project-ALL/releases/tag/numerical-data-2026-10-11)**

The 11 October 2026 storage migration moves **1,846 numerical files (1,818,570,776 bytes)** out of the current source tree. Their contents and original paths are preserved in release ZIP archives. No simulation, coefficient, numerical result or evidence level was changed by this move. Summary JSON/CSV, figures, interactive reports and numerical source remain in the repository.

## Read reports and download data

Continue using the normal report pages. Links to moved arrays now download the containing study archive. Each ZIP retains the repository folder structure; unzip into one common root, not separate folders named after the ZIPs. Some studies span several independent ZIPs: download all of that study's ZIPs for its full data.

The [file catalog](data-release/catalog.json) maps every original filename to an archive and records its size, SHA-256 and Git blob identity. [Local archive verification](data-release/archive-verification.json) records the full restore check. The release also contains a SHA256SUMS.txt list.

## Restore data into a current checkout

From the repository folder, with Python 3.10 or later:

```text
python tools/restore_numerical_data.py
python tools/restore_numerical_data.py --verify-only
```

This downloads and checks every numerical archive, restores the exact original paths, then verifies every restored file. It does not rerun simulations or overwrite files with different contents. Matching existing files are preserved. Restored paths are ignored by Git so the large files are not accidentally recommitted.

To restore one study, append its repository path, for example `--study reports/clay-cell-stress`. Omit that option before running the complete repository checks. Automated checks restore their required data before validation; the original package and numerical checks remain in place.

## Keep a complete offline copy of the original publication

The release includes **every one of the 6,627 original files** from commit [`893a778806bf`](https://github.com/NousVolition/My-Sources-Project-ALL/tree/893a778806bff49e400dbab6ce3eea4bd08f4dbc): numerical ZIPs plus `reports-and-source.zip`. That source/report archive preserves all original reports, figures, scripts, checksums and historical records exactly as they were.

Download `restore_numerical_data.py` and `catalog.json` from the release into one folder, then run:

```text
python restore_numerical_data.py --catalog catalog.json --full-snapshot --destination restored-original
```

Use a **new empty destination**. This recreates the original complete publication, including original local report links, without mixing later edited documentation into its historical checksum set. For offline use, download all ZIPs too and add `--assets-dir .`. Manual extraction of all ZIPs into one common folder also restores the original layout; the helper additionally verifies all hashes.

## Archive downloads

| Archive | Files | Bytes |
| --- | ---: | ---: |
| [clay-cell-stress-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/clay-cell-stress-numerical-001.zip) | 65 | 12,706,399 |
| [clay-contact-stress-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/clay-contact-stress-numerical-001.zip) | 104 | 48,428,880 |
| [clay-hug-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/clay-hug-numerical-001.zip) | 53 | 35,461,188 |
| [dynamics-reference-tests-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/dynamics-reference-tests-numerical-001.zip) | 358 | 67,951,104 |
| [fluid-memory-erasure-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/fluid-memory-erasure-numerical-001.zip) | 1 | 103,550 |
| [fluid-organization-pilot-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/fluid-organization-pilot-numerical-001.zip) | 47 | 190,141,771 |
| [hug-dynamics-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/hug-dynamics-numerical-001.zip) | 95 | 8,336,207 |
| [hug-stress-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/hug-stress-numerical-001.zip) | 424 | 76,807,322 |
| [matched-stretch-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/matched-stretch-numerical-001.zip) | 39 | 267,210,442 |
| [matched-stretch-numerical-002.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/matched-stretch-numerical-002.zip) | 63 | 266,401,672 |
| [matched-stretch-numerical-003.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/matched-stretch-numerical-003.zip) | 58 | 265,008,081 |
| [matched-stretch-numerical-004.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/matched-stretch-numerical-004.zip) | 90 | 266,198,177 |
| [matched-stretch-numerical-005.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/matched-stretch-numerical-005.zip) | 1 | 3,390,849 |
| [quantum-water-pilot-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/quantum-water-pilot-numerical-001.zip) | 181 | 116,875,414 |
| [water-biology-study-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/water-biology-study-numerical-001.zip) | 33 | 21,018,395 |
| [water-molecule-influence-numerical-001.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/water-molecule-influence-numerical-001.zip) | 234 | 172,965,183 |
| [reports-and-source.zip](https://github.com/NousVolition/My-Sources-Project-ALL/releases/download/numerical-data-2026-10-11/reports-and-source.zip) | 4,781 | 158,417,534 |

## Preservation and limits

- Original result bytes and historical publication receipts are preserved in the numerical release and the complete original-history backup linked below. Where a current download link changed inside a checked package, its current delivery manifest was updated; the original manifest remains in `reports-and-source.zip`.
- The subsequent [verified history cleanup](HISTORY-CLEANUP.md) also removes archived numerical paths from normal clone history. Original commits and every historical result remain recoverable from the complete Git backup on the release.
- The archive is a dated snapshot. Later studies and completions require their own verified additions and are not implied by this release.
- Reading the kept report pages does not require restoring arrays. Running analyses or package verification that load those arrays does require restoration first.

[Published-download verification](data-release/remote-verification.json) records an independent download of every release asset, exact agreement with the locally restored archives, and validation of every archived member.
