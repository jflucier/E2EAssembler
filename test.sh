#!/bin/bash

set -e

if [[ -z "${E2EAssembler}" ]]; then
    echo "E2EAssembler install path variable must be defined: export E2EAssembler=/path/to/E2EAssembler"
    exit 1
fi

if [ ! -f "${E2EAssembler}/E2EAssembler.config" ]; then
    echo "E2EAssembler config file not found. It should be found in this path: ${E2EAssembler}/E2EAssembler.config"
    exit 1
fi

echo "E2EAssembler=${E2EAssembler}"

source ${E2EAssembler}/E2EAssembler.config

${E2EAssembler}/00_check_environment.sh
