#!/bin/bash
inp="$1"
#
if [ ! -r "$inp" ] ; then
  echo "$inp is not readable"
  exit
fi
mv -i "$inp" "${inp}.bak"
sed -e 's|/home/ps/.*/|./|' "${inp}.bak" > "${inp}"
