## Resubmission

This is a resubmission of the first submission of gdpar 0.1.0. The previous
upload was stopped by the incoming pre-tests with 1 ERROR and 1 NOTE. The
ERROR is fixed; the NOTE items are addressed/explained below.

### Fix for the ERROR (test failures without the Suggests package 'cmdstanr')

The pre-test reported `[ FAIL 29 | WARN 0 | SKIP 168 | PASS 3394 ]`. All 29
failures had the same root cause: several fitting entry points checked for the
optional Suggests package 'cmdstanr' *before* validating their arguments, so on
a machine without 'cmdstanr' the tests that assert input-validation errors
received the "missing dependency" condition first.

The dependency check has been moved so that all pure-R argument validation runs
first and 'cmdstanr' is required only immediately before it is actually used
(i.e. right before `cmdstanr::cmdstan_model()`). No functionality changed when
'cmdstanr' is installed; the package now also reports its own validation errors
when 'cmdstanr' is absent, and the full test suite passes with 'cmdstanr' not
installed. Verified locally by running the suite against a library path from
which only 'cmdstanr' was removed: `[ FAIL 0 | WARN 0 | SKIP 148 | PASS 3874 ]`.

### NOTE

* Possibly misspelled words in DESCRIPTION: "eBird" and "hypernetworks".
  Both are spelled correctly. "eBird" is the proper name of the Cornell Lab of
  Ornithology citizen-science project that supplies the avian abundance data;
  "hypernetworks" is the established machine-learning term for the amortized
  Path-3 estimator.

* Suggested packages not in a mainstream repository: 'cmdstanr' and 'INLA'.
  Both are listed in `Suggests` only and used conditionally
  (`requireNamespace()` guards). Their locations are declared in
  `Additional_repositories` (the Stan r-universe and the INLA download mirror).
  The package builds, checks and runs its tests without either of them.

## Test environments

* Local: Fedora Linux, R 4.6.0
* Local, with 'cmdstanr' removed from the library path (to mirror the CRAN
  check machines): full test suite passes.

## R CMD check results

0 errors | 0 warnings | 1 note

The remaining note is the new-submission / spelling / Additional_repositories
note explained above.
