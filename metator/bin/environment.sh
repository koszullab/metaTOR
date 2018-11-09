#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Théo Foutel-Rodier
# Lyam Baudry

#A bunch of functions that are called as needed along the pipeline to ensure the stuff you want to execute is basically there.

set -o pipefail
set -e

#For people who run scripts like this on a cluster - modify as needed
module_script=/local/gensoft2/adm/etc/profile.d/modules.sh

if [ -f $module_script ]; then
  . $module_script
  module purge
fi

#Syntactic sugar to make the other functions below easier to read
function there_is() {
  command -v "$1" >/dev/null 2>&1
}

#More syntactic sugar on top of the syntactic sugar
function check_for() {
  local program=$1
  if ! there_is "$program"; then
    if there_is module; then
      if ! module load "$program"; then
        echo "Loading $program failed. Aborting."
        exit 1
      fi
    else
      echo "$program not found. Aborting."
      exit 1
    fi
  fi
}

#A very dirty function that checks in various folders for tools, then sets the appropriate variables for what to call.
#Most use cases should have been handled there.
function locate_and_set_executable() {

  local _exec_name=$1         #Variable to call whenever executable is needed
  local _exec=$2              #Program to locate
  local _program=${3:-$_exec} #If the executable's name is different from that of the program
  local _path=$4              #If a specific path is specified

  if [ -f "$_path/$_program/$_exec" ]; then
    eval "$_exec_name"=\$_path/\$_program/\$_exec
  elif [ -f "$_path/$_exec" ]; then
    eval "$_exec_name"=\$_path/\$_exec
  elif [ -f "$_path" ] && [ ! -z "${_path}" ]; then
    eval "$_exec_name"=\$_path
  elif [ -f "$tools_dir/$_program/$_exec" ]; then
    eval "$_exec_name"=\$tools_dir/\$_program/\$_exec
  elif there_is "$_exec"; then
    eval "$_exec_name"=\$_exec
  elif there_is module; then
    if ! module load "$_exec"; then
      echo "Loading $program failed. Aborting."
      exit 1
    else
      eval "$_exec_name"=\$_exec
    fi

  else
    echo "Error! I couldn't locate $_exec in any way. Please verify that it is in your \$PATH or specify a folder with the --tools-dir option."
    exit 1
  fi
}

#Makes script die somewhat gracefully if something goes wrong in one of the pipe chains (for bash and zsh),
#in case set -o pipefail becomes problematic in the future
function catch_failures() {
  if [ "$BASH" ]; then
    pipe_status=$PIPESATUS
  elif [ "$ZSH_NAME" ]; then
    pipe_status="$pipestatus"
  fi
  if [ "$pipe_status" ]; then
    for p in "${pipe_status[@]}"; do
      if [ ! "$p" -eq 0 ]; then
        echo >&2 "Something wrong happened. Aborting."
        exit 1
      fi
    done
  fi
}
