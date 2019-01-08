#! /bin/bash
mv configure.in configure_local.in
cp configure_minimum_macOS.in configure.in
git push
cp configure_local.in configure.in
