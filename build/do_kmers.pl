#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use Data::Dumper;

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
my $KMER = shift(@ARGV) // 7;
my $STEP = shift(@ARGV) // 3;
my $MAXSAMP = shift(@ARGV) // 1E11; # for testing
my $MIN_FRAC = 0;
my @DEFER = (); # run at end

say STDERR "Generating: K=$KMER S=$STEP N=$MAXSAMP";

#say "set -uex -o pipefail"; ## won't push through
say "set -ux";
say "CSVTK_THREAD=3";
say "SEQKIT_THREADS=3";
say "LC_ALL=C"; 

my $name = "kmers.$KMER-$STEP-$MAXSAMP";
my $FOFN = "$name.fofn";
say "rm -f $FOFN";

my @kfile;

for my $g (sort(path("groups.txt")->lines({chomp=>1}))) {
  my @kmers;
  say "mkdir -p $g";
  my $fofn = "$g/$name.fofn";
  say "rm -f $fofn"; 
  my $N=0;
  for my $a ( path("$GTDB/$g/proteomes.txt")->lines({chomp=>1}) ) {
    my $out = "$g/$a.$name.gz";
    say kmers("$GTDB/$g/$a", $out);
    say "echo $out >> $fofn";
    $N++;
    last if $N >= $MAXSAMP;
  }

  my $kmers = "$g/$name";  
  push @kfile, $kmers;
  say "echo $kmers >> $FOFN";
  # --sort-by-key so we can do merge sort later
  defer(clean("
    xargs -a $fofn zcat
    | tsvtk freq -H --sort-by-key
    | tsvtk mutate2 -H --at 2
      -w 6 -e '\$2 / $N'
    | tee $kmers.freq
    | cut -f 1 > $kmers
  "));
}

say $_ for @DEFER;

say clean("
  tr '\\n' '\\0' < $FOFN
  | LC_ALL=C sort --merge --files0-from=-
  | uniq -u
  > $name.uniq
");

say "vmtouch -t $name.uniq";

foreach (@kfile) {
  say clean("
    LC_ALL=C grep -v -F -f $name.uniq $_.freq
    | tsvtk sort -k 3:nr -o $_.uniq
  ");
}

say STDERR "Run: parallel -j 100 -v -k -a ${KMER}-$STEP-$MAXSAMP.sh";

exit(0);

#...........................................
sub defer {
  push @DEFER, @_;
}
#...........................................
sub kmers {
  my($in,$out) = @_;
  return clean("
    zcat $in 
    | tantan -p -x X
    | seqkit sliding -W $KMER  -s $STEP
                     --quiet -t protein
    | seqkit seq -w 0 -u --seq 
                 --quiet -t protein
    | grep -v '[*XUOW]'
    | tsvtk uniq -H -o $out
  ");
}
#...........................................
sub clean {
  return join ' ', 
  split ' ', join ' ', @_;
}
#...........................................
