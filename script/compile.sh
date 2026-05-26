#!/usr/bin/env sh
set -eu

project_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
build_dir="${project_root}/build"

cmake -S "${project_root}" -B "${build_dir}"
cmake --build "${build_dir}"
