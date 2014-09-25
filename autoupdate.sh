#! /bin/bash
mv configure.in configure_local.in
git stash save --keep-index
git pull origin master
mv configure_local.in configure.in
