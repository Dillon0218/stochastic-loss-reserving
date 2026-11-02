#!/usr/bin/env bash
# Download the CAS Loss Reserving Database into data/raw/.
#
# Source: "Loss Reserving Data Pulled From NAIC Schedule P", Casualty Actuarial
# Society, compiled by Glenn Meyers and Peng Shi (the "earlier data sets",
# updated September 1, 2011). Accident years 1988-1997, ten development lags,
# upper AND lower triangles, so a model fitted to the 1997 diagonal can be
# scored against what actually happened.
#
# The CAS asks users to acknowledge that the files were accessed from its web
# page. The page states no redistribution licence, so the raw files are not
# committed; this script rebuilds them.
#
#   ppauto_pos.csv    private passenger auto liability/medical   (suffix _B)
#   comauto_pos.csv   commercial auto/truck liability/medical    (suffix _C)
#   wkcomp_pos.csv    workers' compensation                      (suffix _D)
#   medmal_pos.csv    medical malpractice, claims made           (suffix _F2)
#   othliab_pos.csv   other liability, occurrence                (suffix _h1)
#   prodliab_pos.csv  products liability, occurrence             (suffix _R1)

set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw

BASE="https://www.casact.org/sites/default/files/2021-04"
for f in ppauto_pos comauto_pos wkcomp_pos medmal_pos othliab_pos prodliab_pos; do
  echo "Downloading ${f}.csv"
  curl -fsSL --retry 3 -A "Mozilla/5.0" -o "data/raw/${f}.csv" "${BASE}/${f}.csv"
done

ls -lh data/raw/
echo "Done. Now run: Rscript R/run_all.R"
