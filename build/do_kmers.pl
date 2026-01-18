#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use Data::Dumper;

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
my $PREFIX = "ALL";
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

my $name = "kmers.${KMER}_${STEP}";
my $FOFN = "$PREFIX.$name.fofn";
say "rm -f $FOFN";

fsorfor my $g (sort(path("groups.txt")->lines({chomp=>1}))) 
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
  > $PREFIX.$name.uniq
");
#say clean("
#  xargs -a $FOFN cat 
#  | tsvtk cut -H -f 1
#  | tsvtk freq -H -n -r 
#  | tee $PREFIX.$name.freq
#  | tsvtk filter2 -H -f '\$2 == 1'
#  | cut -f 1
#  > $PREFIX.$name.uniq
#");

#  | grep \$'\\t1\$' $PREFIX.$name.freq | cut -f 1 > $PREFIX.$name.uniq";

say STDERR "Run: parallel -j 64 -v -k -a $PREFIX.${KMER}_$STEP.sh";

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
