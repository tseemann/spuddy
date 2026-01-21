#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use Data::Dumper;

my $GZIP = "libdeflate-gzip";
my $COMPRESS = "pigz -p \$CPUS";
my $ZCAT = "libdeflate-gzip -d -c -f";
my $CACHE = "vmtouch -q -t";

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
@ARGV or die "$0 KMER STEP MAXSAMP";
my $KMER = shift(@ARGV) // 7;
my $STEP = shift(@ARGV) // 1;
my $MAXSAMP = shift(@ARGV) // 50; # for testing
my $MIN_FRAC = 0;
my @DEFER = (); # run at end

say STDERR "Generating: K=$KMER S=$STEP N=$MAXSAMP";
my $name = "kmers.$KMER.$STEP.$MAXSAMP";
say STDERR "Prefix = $name";

open my $OUT, '>', "$name.sh";
select $OUT;

#say "set -uex -o pipefail"; ## won't push through
say "#/usr/bin/env bash";
say "set -eu -o pipefail";
say "CSVTK_THREADS=1";
say "SEQKIT_THREADS=1";
say "LC_ALL=C"; 
say "CPUS=\$(( \$(nproc) / 2))";

my $LOGFILE = "$name.sh.log";
say "rm -f $LOGFILE";
my $PARA = "parallel --bar --joblog $LOGFILE";
my $PARAJ = "$PARA -j \$CPUS";
my $G = 'groups.txt';

banner("Making folders");
#say "xargs -a $G mkdir -p";
say "$PARAJ -a $G mkdir -p";
my @G = sort(path($G)->lines({chomp=>1}));

my @faa;
my @gcmd;
my @kfile;
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
banner("Generating kmers per isolate");
say "$PARAJ -a $faa_fofn ".dq(faa2pep());

# count uniq kmers per isolate
my $gcmd = path("$name.gcmd");
$gcmd->spew(map { "$_\n" } @gcmd);
banner("Finding unique kmers per isolate");
say "$PARAJ -a $gcmd";

# compile all the uniqmers
banner("Identifying global singleton kmers (slow)");
say "$PARA -j 1 -a $G $CACHE {}/$name";
say clean("
  [[ -r '$name.uniq' ]] ||
  xargs -a $G -I % echo %/$name 
  | tr '\\n' '\\0'
  | LC_ALL=C sort --merge --files0-from=-
  | uniq -u
  > $name.uniq
");

# For each species, find their species-specific ones
say "$CACHE -t $name.uniq";
banner("Identifying kmers unique to each speices");
#say "$PARA -j 1 -a $G $CACHE {}/$name";
say "$PARAJ -a $G ".dq(clean("
    [[ -r '{}/$name.uniq' ]] ||
    LC_ALL=C sort --merge $name.uniq {}/$name
    | uniq -d > {}/$name.uniq
    "));

# combine final uniqmers with our porportion data
banner("Add scores to species-specific kmers");
say "$PARAJ -a $G ".dq(clean("
       tsvtk join -H --left-join
             {}/$name.uniq {}/$name.freq
     | tsvtk sort -H -k 2:nr
              -o {}/$name.uniq.freq
    "));

# collate stats of uniqmers for each species
banner("Collect some species stats");
say "$PARAJ -a $G 'wc -l {}/$name.uniq' | sort -bnr > $name.counts";

# collate all the species kmers/socre
banner("Make final databases");
say "$PARAJ -a $G ".dq(clean("
       tsvtk mutate2 -H -e \\'{}\\'
             {}/$name.uniq.freq
    "))." > $name.db";
#    "))." | $COMPRESS > $name.db.gz";

# generate cut down versions
say "$CACHE $name.db";
say qq{echo -e "50\\n75\\n90\\n95\\n" | $PARAJ -v "tsvtk filter2 -H -f '\\\$2 >= 0.{}' $name.db > $name.db.{}pc"};
say "wc -l $name.db*";
say "\\ls -1sh $name.db*";

# finish up
banner("Pipeline completed!");
say "echo  less -S $name.counts";

system "chmod +x $name.sh";
say STDERR "Run this!\n  ./$name.sh";
exit(0);

#...........................................
sub species_pepmer {
  my($dir,$name,$num) = @_;
  my $outname = "$dir/$name";
  return clean("
    [[ -f '$outname' ]] ||
    xargs -a $dir/$name.list.fofn -I % $ZCAT $dir/%
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
    $ZCAT $GTDB/{}
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
sub banner {
  my($text) = @_;
  my $width = $ENV{'COLUMNS'} || 40; 
  my $line = '='x($width-2);
  say "echo '$line'";
  say "echo '$text'";
  #say "echo '$line'";
}
#...........................................
