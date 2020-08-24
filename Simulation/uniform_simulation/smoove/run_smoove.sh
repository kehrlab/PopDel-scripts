#!/bin/bash

WD=.

cd $WD

./scripts/smoove/smoove_single.sh ## single sample complete workflow.
./scripts/smoove/smoove_1.sh  ## call
./scripts/smoove/smoove_2.sh  ## merge
./scripts/smoove/smoove_3.sh  ## genotype
./scripts/smoove/smoove_4.sh  ## paste
