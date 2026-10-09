# Validation record

Validated locally with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8.

- **Study tests: 24 passed.** These include an independent majority-rule oracle for all 512 B starts; rotation, reflection and color symmetry; protected adjacent pairs; initial consensus; failure versus success at exactly 120 s; tied initiators; refusal and initiative boundaries; identity relabeling; copy-credit conservation; balanced orders and crossed positions; known-effect recovery; end-to-end output and corruption detection.
- **Full repository suite: 290 passed, 3 expected failures.** The expected failures belong to existing mathematical-model checks. No new expected failures were added.
- **Static checks:** the repository's fatal-syntax/undefined-name Flake8 selection passed. The staged change passed `git diff --check`; generated SVG whitespace is excluded because it comes from Matplotlib's path formatting.
- **Full reproduction:** 45,852 rounds plus 512 exact B starts were run repeatedly. The first repeat matched all 12 numerical data files byte for byte. After making JSON/CSV line endings explicitly portable, another full run matched all 12 files exactly after normalizing the earlier Windows line endings. The saved final files use LF; the verifier now checks their unmodified hashes.
- **Visual inspection:** all four PNG figures were opened and checked for legible labels, clipping, denominators, and simulation-only labeling. Companion SVG files are included for export. Local report links and packaged source/data hashes are checked by the packaging process.
- **CI:** the repository workflow includes the new tests automatically and now reproduces the complete study without plots, followed by exact numerical-file verification. A local pass is not a claim that remote CI has run.

The first full repository test attempt encountered dependency-access problems in the local Windows sandbox. Once a consistent environment was prepared, six unrelated archive hash checks exposed Git's automatic Windows line-ending conversion. The fresh checkout's five affected archive files were restored to their exact existing Git bytes only after verifying that line endings were their sole difference and their canonical hashes matched the original provenance manifest. No archived source, data, test or expected hash was changed in this branch. The final complete test run passed as stated above.

These checks verify implementation and reproducibility. They do not validate C's behavioral assumptions against people. No volunteer testing, empirical parameter fitting, ethics approval, Freud test, or molecular-dynamics validation has occurred.
