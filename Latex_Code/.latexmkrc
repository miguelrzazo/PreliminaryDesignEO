## latexmk configuration: use pdflatex with -shell-escape and biber
$pdflatex = 'pdflatex -interaction=nonstopmode -shell-escape %O %S';
$biber = 'biber %O %S';
$pdf_mode = 1;
$clean_ext = 'synctex.gz fls aux bbl bcf run.xml';

## Los ficheros intermedios (.aux, .toc, .bcf, .bbl, .log...) se escriben fuera
## de la carpeta sincronizada por iCloud, que los corrompe a mitad de escritura
## y genera duplicados del tipo "00main 2.aux". En la carpeta del proyecto solo
## se escribe el PDF final, en una unica operacion al terminar.
$aux_dir = "$ENV{HOME}/Library/Caches/TFG-PreliminaryDesignEO/latex";
$out_dir = '.';

## pdflatex escribe un .aux por cada capitulo incluido, respetando la ruta
## relativa. Esas subcarpetas deben existir de antemano dentro de $aux_dir o la
## compilacion aborta con "I can't write on file".
use File::Find ();
use File::Path ();
File::Path::make_path($aux_dir);
File::Find::find(
    {
        wanted => sub {
            return unless -d $File::Find::name;
            return if $_ eq '.';
            my $rel = $File::Find::name;
            $rel =~ s{^\./}{};
            return if $rel =~ m{^\.};
            File::Path::make_path("$aux_dir/$rel");
        },
        no_chdir => 0,
    },
    '.'
);
