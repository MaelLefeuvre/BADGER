#!/usr/bin/env bash

set -euo pipefail

R --slave -e 'devtools::install_github("jonotuke/BREADR@v1.04')'
