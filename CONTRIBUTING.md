# Contributing

Issues and pull requests are welcome: bug reports, new checkers, and new fault models.

- **Keep the records reproducible.** The campaign must stay bit-identical on x86-64. Run `make test` before opening a pull request. If you change input generation or arithmetic on purpose, regenerate `tests/quick_reference.sha256` and explain why in the pull request.
- **Adding a checker.** If it can be expressed with the stored normalised discrepancies, add it to `fire()` and `clean_fire()` in `analysis/analyze.py`; no new campaign is needed. Otherwise:
  - add the statistic to the records in `src/common.hpp` and `src/campaign.cpp`;
  - implement it independently in `src/verify.cpp`;
  - document it in `docs/methodology.md`.
- **Style.** Use C++17 with no external dependencies. The analysis may use only Python and NumPy.
