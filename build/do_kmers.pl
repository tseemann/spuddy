#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use Data::Dumper;

my $GZIP = "libdeflate-gzip";
my $CACHE = "vmtouch -q -t";
my $PARA = "parallel";

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
@ARGV or die "$0 KMER STEP MAXSAMP";
my $KMER = shift(@ARGV) // 6;
my $STEP = shift(@ARGV) // 1;
my $MAXSAMP = shift(@ARGV) // 999999; # for testing
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
my $PARAJ = "$PARA -j \$CPUS -v";
my $G = 'groups.txt';
#say "xargs -a $G mkdir -p";
say "$PARAJ -a $G mkdir -p";
my @G = sort(path($G)->lines({chomp=>1}));
#my $GN = scalar(@GN);

my @faa;
#my @uniq;
my @gcmd;
for my $g (@G) {
  say STDERR "Prepping: $g";
  my @kmerf;
  my $N=0;
  for my $a ( path("$GTDB/$g/proteomes.txt")->lines({chomp=>1}) ) {
    push @faa, "$g/$a";
    push @kmerf, "$a.$KMER.$STEP.gz";
    last if ++$N >= $MAXSAMP;
  }
  my $cmd = species_pepmer($g,$name,$N);
  #push @uniq, "$g/$name";
  push @gcmd, $cmd;
  
  # keep track of subset of files we kmered
  my $kmerf = path("$g/$name.list.fofn");
  $kmerf->spew(map { "$_\n" } @kmerf);
}
my $faa_fofn = path("$name.faa.fofn");
$faa_fofn->spew(map { "$_\n" } @faa);
say "$PARAJ -a $faa_fofn ".dq(faa2pep());

# count uniq kmers per isolate
my $gcmd = path("$name.gcmd");
$gcmd->spew(map { "$_\n" } @gcmd);
say "$PARAJ -a $gcmd";

# compile all the uniqmers
say "$PARA -j 1 -a $G $CACHE {}/$name";
say clean("
  xargs -a $G -I % echo %/$name 
  | tr '\\n' '\\0'
  | LC_ALL=C sort --merge --files0-from=-
  | uniq -u
  > $name.uniq
");

# For each species, find their species-specific ones
say "$CACHE -t $name.uniq";

#say "$PARA -v -j 1 -a $G $CACHE {}/$name";
say "$PARAJ -a $G ".dq(clean("
    LC_ALL=C sort --merge $name.uniq {}/$name
    | uniq -d > {}/$name.uniq
    "));

# combine final uniqmers with our porportion data
say "$PARAJ -a $G ".dq(clean("
       tsvtk join -H --left-join
             {}/$name.uniq {}/$name.freq
     | tsvtk sort -H -k 2:nr
              -o {}/$name.uniq.freq
    "));

# collate all the species kmers/socre
say "$PARAJ -a $G ".dq(clean("
       tsvtk mutate2 -H -e \\'{}\\'
             {}/$name.uniq.freq
     | tsvtk filter2 -H -f '\\$2 >= 0.5'
    "))." > $name.database";

# collate stats of uniqmers for each species
say "$PARAJ -a $G 'wc -l {}/$name.uniq' > $name.counts";

# finish up
say "ehco DONE:";
say "echo  less -S $name.counts";
say STDERR "Run this!\n  bash -x $name.sh";
exit(0);

#...........................................
sub species_pepmer {
  my($dir,$name,$num) = @_;
  my $outname = "$dir/$name";
  return clean("
    [[ -f '$outname' ]] ||
    xargs -a $dir/$name.list.fofn -I % $GZIP -dc $dir/%
    | tsvtk freq -H --sort-by-key
    | tsvtk mutate2 -H --at 2
            -w 6 -e '\$2 / $num'
    | tee $dir/$name.freq
    | cut -f 1 
    > $outname
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
  my $outname = "{}.$KMER.$STEP.gz";
  return clean("
    [[ -f '$outname' ]] ||
    $GZIP -dc $GTDB/{}
    | seqkit sliding -W $KMER -s $STEP
    | tr '[:lower:]' '[:upper:]' 
    | grep -v -E '^>|[XUOW*]'
    | tsvtk uniq -H -o $outname
  ");
}
#...........................................
sub clean {
  return join ' ', 
  split ' ', join ' ', @_;
}
#...........................................
