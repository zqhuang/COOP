#! /bin/bash
mv configure.in configure_local.in
mv include/user_defined_primordial_power.h include/user_defined_primordial_power_local.h
git stash save --keep-index
git pull origin master
mv configure_local.in configure.in
mv include/user_defined_primordial_power_local.h include/user_defined_primordial_power.h 