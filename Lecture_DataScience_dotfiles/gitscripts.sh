#!bin/sh

git add .
git commit -m "fixed"
git push -f origin main

jb build --all .
git checkout main
ghp-import -n -p -f _build/html
rm _build
