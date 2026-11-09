## latexmk configuration: use pdflatex with -shell-escape and biber
$pdflatex = '/Library/TeX/texbin/pdflatex -interaction=nonstopmode -shell-escape %O %S';
$pdf_mode = 1;
$clean_ext = 'synctex.gz fls aux bbl bcf run.xml';

## Los ficheros intermedios (.aux, .toc, .bcf, .bbl, .log...) se escriben fuera
## de la carpeta sincronizada por iCloud, que los corrompe a mitad de escritura
## y genera duplicados del tipo "00main 2.aux". En la carpeta del proyecto solo
## se escribe el PDF final, en una unica operacion al terminar.
$aux_dir = $ENV{TFG_LATEX_BUILD_DIR} || "$ENV{HOME}/Library/Caches/TFG-PreliminaryDesignEO/latex";
$out_dir = '.';

## pdflatex escribe un .aux por cada capitulo incluido, respetando la ruta
## relativa. Esas subcarpetas deben existir de antemano dentro de $aux_dir o la
## compilacion aborta con "I can't write on file".
use File::Find ();
use File::Path ();
File::Path::make_path($aux_dir);

## MacTeX's universal Biber launcher uses `lipo -extract_family`, which the
## current macOS lipo does not support. Use a cached arm64 slice on Apple Silicon.
my $biber_executable = '/Library/TeX/texbin/biber';
my $machine = `/usr/bin/uname -m`;
chomp $machine;
if ($machine eq 'arm64') {
    my $arm_biber = "$aux_dir/biber-arm64";
    my $source_mtime = (stat($biber_executable))[9] || 0;
    my $cached_mtime = (stat($arm_biber))[9] || 0;
    if (!-x $arm_biber || $source_mtime > $cached_mtime) {
        system('/usr/bin/lipo', '-thin', 'arm64', $biber_executable, '-output', $arm_biber) == 0
            or die "Could not extract arm64 Biber from $biber_executable\n";
        chmod 0755, $arm_biber;
    }
    $biber_executable = $arm_biber;
}
$biber = "'$biber_executable' %O %S";
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
