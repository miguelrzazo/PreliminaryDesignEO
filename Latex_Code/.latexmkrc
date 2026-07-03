## latexmk configuration: use pdflatex with -shell-escape and biber
$pdflatex = 'pdflatex -interaction=nonstopmode -shell-escape %O %S';
$biber = 'biber %O %S';
$pdf_mode = 1;
$clean_ext = 'synctex.gz fls aux bbl bcf run.xml';
