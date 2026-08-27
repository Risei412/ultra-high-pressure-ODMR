@echo off
if not exist build mkdir build
latexmk -pdf -interaction=nonstopmode -file-line-error -halt-on-error -outdir=build main_pla.tex
exit /b %errorlevel%