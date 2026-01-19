#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use Data::Dumper;

my $GZIP = "libdeflate-gzip";

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
my $KMER = shift(@ARGV) // 6;
my $STEP = shift(@ARGV) // 1;
my $MAXSAMP = shift(@ARGV) // 9999991; # for testing
my $MIN_FRAC = 0;
my @DEFER = (); # run at end

say STDERR "Generating: K=$KMER S=$STEP N=$MAXSAMP";
my $name = "kmers.$KMER.$STEP.$MAXSAMP";
say STDERR "Prefix = $name";

open my $OUT, '>', "$name.sh";
select $OUT;

#say "set -uex -o pipefail"; ## won't push through
say "#/usr/bin/env bash";
say "set -eux -o pipefail";
say "CSVTK_THREADS=1";
say "SEQKIT_THREADS=1";
say "LC_ALL=C"; 
say "CPUS=\$(( \$(nproc) / 2))";

my @kfile;
my $PARA = "parallel -j \$CPUS -v";
my $G = 'groups.txt';
say "xargs -a $G mkdir -p";
my @G = sort(path($G)->lines({chomp=>1}));

my @faa;
#my @uniq;
my @gcmd;
for my $g (@G) {
  say STDERR "Prepping: $g";
  my @kmerf;
  my $N=0;
  for my $a ( path("$GTDB/$g/proteomes.txt")->lines({chomp=>1}) ) {
    push @faa, "$g/$a";
    push @kmerf, "$a.$name";
    last if ++$N >= $MAXSAMP;
  }
  my $cmd = species_pepmer($g,$name,$N);
  #push @uniq, "$g/$name";
  push @gcmd, $cmd;
  
  # keep track of subset of files we kmered
  my $kmerf = path("$g/$name.counts.fofn");
  $kmerf->spew(map { "$_\n" } @kmerf);
}
my $faa_fofn = path("$name.faa.fofn");
$faa_fofn->spew(map { "$_\n" } @faa);
say "$PARA -a $faa_fofn ".dq(faa2pep());

# count uniq kmers per isolate
my $gcmd = path("$name.gcmd");
$gcmd->spew(map { "$_\n" } @gcmd);
say "$PARA -a $gcmd";

# compile all the uniqmers
say clean("
  xargs -a $G -I % echo %/$name 
  | tr '\\n' '\\0'
  | LC_ALL=C sort --merge --files0-from=-
  | uniq -u
  > $name.uniq
");

# For each species, find their species-specific ones
say "vmtouch -t $name.uniq";
say "parallel -j 8 -a $G ".dq(clean("
    LC_ALL=C grep -F -f $name.uniq {}/$name.freq
    | tsvtk sort -H -k 3:nr -o {}/$name.freq.uniq
  "));

# collate stats of uniqmers for each species
say "parallel -a $G 'wc -l {}/$name.freq.uniq' > $name.counts";

# finish up
say STDERR "Run this!\n  bash -x $name.sh";
exit(0);

#...........................................
sub species_pepmer {
  my($dir,$name,$num) = @_;
  return clean("
    xargs -a $dir/$name.counts.fofn -I % cat $dir/%
    | tsvtk freq -H --sort-by-key
    | tsvtk mutate2 -H --at 2
            -w 6 -e '\$2 / $num'
    | tee $dir/$name.freq
    | cut -f 1 > $dir/$name
  ");
}
#...........................................
sub dq { 
  $_ = shift;
  s/"/\\"/g;
  return qq/"$_"/;
}
#...........................................
sub faa2pep {
  return clean("
    [[ -f '{}.$name' ]] ||
    $GZIP -dc $GTDB/{}
    | tantan -p
    | seqkit sliding -w 0 -W $KMER -s $STEP
    | grep -v '^>' 
    | grep -v '[a-zXUOW*]'
    | tsvtk uniq -H -o {}.$name
  ");
}
#...........................................
sub clean {
  return join ' ', 
  split ' ', join ' ', @_;
}
#...........................................
